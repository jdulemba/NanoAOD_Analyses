from numba import njit#, objmode
import numpy as np
from pdb import set_trace
import python.TTBarSolver as ttsolver
import compiled.pynusolver as pynusolver
import awkward as ak

solver = None
def year_to_run(**kwargs):
    global solver
    solver = ttsolver.TTSolver(kwargs["year"])

def run():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="Specify which year to run over")
    args = parser.parse_args()
    year_to_run(**kwargs)

@njit()
def find_best_3j_permutations(njets_array, jets_btag, j1_inds, j2_inds, massdisc, nsdisc, use_merged=False, btag_req=True):
    """
        Inputs:
            1D array with number of jets per event
            2D array of jets (px, py, pz, E, btag pass)
            2D array of leptons (px, py, pz, E)
            2D array of MET (px, py)
        Returns:
            [0]: List of jet assignment ordering
            [1]: List of neutrino solutions
            [2]: List of probabilities (Total, mass discriminant, nu discriminant)
    """
    start = 0
    stop = 0
    evt_idx = 0

    best_perms_ordering = []
    best_perms_probs = []

    #set_trace()
    ## get best perms for 3 jet events
    for njets in njets_array:
        #print("evt idx: ", evt_idx)#, ", njets: ", njets)
        stop += njets

        ## initialize best_perm lists for event
        best_perm_ordering = np.array([-10, -10, -10])
        best_perm_probs = np.array([np.inf, np.inf, np.inf])

        for j0 in range(start, stop):
                # first jet has to be btagged
            if btag_req:
                if jets_btag[j0] < 0.5: continue

            ns_prob = nsdisc[evt_idx][j0-start]

            ## jet info will change for each value of j0
            #set_trace()
            for j1 in range(start, stop):
                if j1 == j0: continue
                    # second jet has to be btagged
                if btag_req:
                    if jets_btag[j1] < 0.5: continue
                for j2 in range(j1+1, stop):
                    if j2 == j0: continue
                    #    # btagged jets can"t be assigned to wjets
                    #if jets_btag[j2, 4] > 0.5: continue
                    #set_trace()
                    for idx in range(j1_inds[evt_idx].size):
                        jet1, jet2 = j1_inds[evt_idx][idx], j2_inds[evt_idx][idx]
                        if ((j0-start == jet1) or (j0-start == jet2)): continue
                        corr_idx = idx
                    mass_prob = massdisc[evt_idx][corr_idx]
                    #print("\tns prob: %f, mass prob: %f, tot prob: %f" % (ns_prob, mass_prob, ns_prob+mass_prob))
                    #print("\tj0: %i, j1: %i, j2: %i" % (j0-start, j1-start, j2-start))
                    if ns_prob+mass_prob < best_perm_probs[0]: ## have to be careful with indexing!!
                        best_perm_probs = np.array([ns_prob+mass_prob, mass_prob, ns_prob])
                        best_perm_ordering = np.array([j0-start, j1-start, j2-start])#, -10]) # set last element to -10 in order to make array have the same shape as 4+ jets

        #print("Ordering:", best_perm_ordering, ", probs:", best_perm_probs)
        best_perms_ordering.append(best_perm_ordering)
        best_perms_probs.append(best_perm_probs)

        start += njets
        evt_idx += 1
    return best_perms_ordering, best_perms_probs


def get_best_3j_permutations(jets, nschi, jets_btagging, btag_req):
        # get combinations of two jets
    j0_comb, j1_comb = ak.unzip(ak.combinations(jets, 2))
    j1_argcombs, j2_argcombs = ak.unzip(ak.argcombinations(jets, 2))
        # invariant mass of two jets
    m2jets = (j0_comb+j1_comb).mass

    mass_discr, ns_discr = solver.ak_solve_3J_lost(m2jets=m2jets, nschi=nschi)

    #set_trace()
    bp_ordering, bp_probs = find_best_3j_permutations(njets_array=ak.num(jets), jets_btag=jets_btagging, j1_inds=ak.to_numpy(j1_argcombs), j2_inds=ak.to_numpy(j2_argcombs),
         massdisc=ak.to_numpy(mass_discr), nsdisc=ns_discr, btag_req=btag_req)

    return bp_ordering, bp_probs


@njit()
def find_best_4pj_permutations(njets_array, jets_btag, mthad_jet_inds, mwhad_jet_inds, mthad_starts_stops, mwhad_starts_stops, massdisc, nsdisc, use_merged=False, btag_req=True):
    """
        Inputs:
            1D array with number of jets per event
            2D array of jets (px, py, pz, E, btag pass)
            2D array of leptons (px, py, pz, E)
            2D array of MET (px, py)
        Returns:
            [0]: List of jet assignment ordering
            [1]: List of neutrino solutions
            [2]: List of probabilities (Total, mass discriminant, nu discriminant)
    """
    start = 0
    stop = 0
    evt_idx = 0
    thad_inds_start, thad_inds_stop = 0, 0
    whad_inds_start, whad_inds_stop = 0, 0
    mdisc_start = 0

    best_perms_ordering = []
    best_perms_probs = []

    #set_trace()
    for njets in njets_array:
        #print("evt idx: ", evt_idx)
        #print("evt idx: ", evt_idx, ", njets: ", njets)
        stop += njets

        #set_trace()
        whad_inds_stop += mwhad_starts_stops[evt_idx]
        thad_inds_stop += mthad_starts_stops[evt_idx]

        ## initialize best_perm lists for event
        best_perm_ordering = np.array([-10, -10, -10, -10])
        best_perm_probs = np.array([np.inf, np.inf, np.inf])

        for j0 in range(start, stop):
                # first jet has to be btagged
            if btag_req:
                if jets_btag[j0] < 0.5: continue

            ns_prob = nsdisc[evt_idx][j0-start]
            if ns_prob == np.inf: continue

            ## jet info will change for each value of j0
            ## get best perms for 4+ jet events
            for j1 in range(start, stop):
                if j1 == j0: continue
                    # second jet has to be btagged
                if btag_req:
                    if jets_btag[j1] < 0.5: continue
                for j2 in range(start, stop):
                    if j2 == j0 or j2 == j1: continue
                    for j3 in range(j2+1, stop):
                        if j3 == j0 or j3 == j1: continue

                        for mw_idx in range(whad_inds_start, whad_inds_stop):
                            mw_j1, mw_j2 = mwhad_jet_inds[mw_idx]
                            if ((j2-start == mw_j1) or (j2-start == mw_j2)) and ((j3-start == mw_j1) or (j3-start == mw_j2)):
                                corr_mw_idx = mw_idx
                                corr_mw_j1_idx = mw_j1
                                corr_mw_j2_idx = mw_j2
                            else: continue

                        for mt_idx in range(thad_inds_start, thad_inds_stop):
                            mt_j1, mt_j2, mt_j3 = mthad_jet_inds[mt_idx]
                            if ((corr_mw_j1_idx == mt_j1) or (corr_mw_j1_idx == mt_j2) or (corr_mw_j1_idx == mt_j3)) and \
                                ((corr_mw_j2_idx == mt_j1) or (corr_mw_j2_idx == mt_j2) or (corr_mw_j2_idx == mt_j3)) and \
                                ((j1-start == mt_j1) or (j1-start == mt_j2) or (j1-start == mt_j3)):
                                corr_mt_idx = mt_idx
                            else: continue

                        mass_prob = massdisc[mdisc_start+(corr_mt_idx-thad_inds_start)*(whad_inds_stop-whad_inds_start)+(corr_mw_idx-whad_inds_start)]
                        #print("\tj0: %i, j1: %i, j2: %i, j3: %i" % (j0-start, j1-start, j2-start, j3-start))
                        #print("\tProb = ", ns_prob+mass_prob, ", MassDiscr = ", mass_prob, ", NuDiscr = ", ns_prob)
                        if ns_prob+mass_prob < best_perm_probs[0]: ## have to be careful with indexing!!
                            best_perm_probs = np.array([ns_prob+mass_prob, mass_prob, ns_prob])
                            best_perm_ordering = np.array([j0-start, j1-start, j2-start, j3-start])#, -10]) # set last element to -10 in order to make array have the same shape as 4+ jets

        #print("evt idx: ", evt_idx, ", Ordering:", best_perm_ordering, ", probs:", best_perm_probs)
        best_perms_ordering.append(best_perm_ordering)
        best_perms_probs.append(best_perm_probs)

        start += njets

        mdisc_start += (thad_inds_stop-thad_inds_start)*(whad_inds_stop-whad_inds_start)
        whad_inds_start += mwhad_starts_stops[evt_idx]
        thad_inds_start += mthad_starts_stops[evt_idx]

        evt_idx += 1

    return best_perms_ordering, best_perms_probs

def get_best_4pj_permutations(jets, nschi, jets_btagging, btag_req):
        # get combinations of two jets
    j0_2comb, j1_2comb = ak.unzip(ak.combinations(jets, 2))
    twoJets_argcombs = ak.argcombinations(jets, 2)
    m2jets = (j0_2comb+j1_2comb).mass # invariant mass of two jets

        # get combinations of three jets
    j0_3comb, j1_3comb, j2_3comb = ak.unzip(ak.combinations(jets, 3))
    threeJets_argcombs = ak.argcombinations(jets, 3)
    m3jets = (j0_3comb+j1_3comb+j2_3comb).mass # invariant mass of two jets

        # repeat values of mThad for each value of mWhad
    mthad_argcombs, mwhad_argcombs = ak.unzip(ak.argcartesian([m3jets, m2jets]))
    mthad_vals, mwhad_vals = ak.unzip(ak.cartesian([m3jets, m2jets]))

    mass_discr, ns_discr = solver.ak_solve_4PJ(mthad=mthad_vals, mwhad=mwhad_vals, nschi=nschi)

    mthad_argcombs_np = np.asarray(ak.to_list(ak.flatten(threeJets_argcombs)))
    mthad_inds_starts_stops = ak.to_numpy(ak.num(threeJets_argcombs, axis=1))
    mwhad_argcombs_np = np.asarray(ak.to_list(ak.flatten(twoJets_argcombs)))
    mwhad_inds_starts_stops = ak.to_numpy(ak.num(twoJets_argcombs, axis=1))

    #set_trace()
    bp_ordering, bp_probs = find_best_4pj_permutations(njets_array=ak.num(jets), jets_btag=jets_btagging, massdisc=ak.to_numpy(ak.flatten(mass_discr)), nsdisc=ns_discr, btag_req=btag_req,
        mthad_jet_inds=mthad_argcombs_np, mwhad_jet_inds=mwhad_argcombs_np, mthad_starts_stops=mthad_inds_starts_stops, mwhad_starts_stops=mwhad_inds_starts_stops)
    #set_trace()

    return bp_ordering, bp_probs


def find_best_permutations(jets, leptons, MET, btagWP, btag_req=True):
    """
    Inputs:
        Jets, leptons, MET, and btag working point for jets
    Returns:
        awkward Table containing
            1: Jet objects (BLeps, BHads, WJas, WJbs)
            2: Leptons, Neutrinos
            3: Calculated probabilities (Total -> Prob, mass discriminant -> MassDiscr, neutrino discrminant -> NuDiscr)
    """

        # make leptons same shape as jets
    leps_rep, jets_rep = ak.unzip(ak.cartesian([leptons, jets]))
    met_rep, jets_rep = ak.unzip(ak.cartesian([MET, jets]))
    lep_inputs = np.stack((ak.to_numpy(ak.flatten(leps_rep.px)), ak.to_numpy(ak.flatten(leps_rep.py)), ak.to_numpy(ak.flatten(leps_rep.pz)), ak.to_numpy(ak.flatten(leps_rep.energy))), axis=1).astype("float64") # one row has (px, py, pyz, E)
    met_inputs = np.stack((ak.to_numpy(ak.flatten(met_rep.px)), ak.to_numpy(ak.flatten(met_rep.py))), axis=1).astype("float64") # one row has (px, py)
    jets_inputs = np.stack((ak.to_numpy(ak.flatten(jets.px)), ak.to_numpy(ak.flatten(jets.py)), ak.to_numpy(ak.flatten(jets.pz)), ak.to_numpy(ak.flatten(jets.energy)), ak.to_numpy(ak.flatten(jets[btagWP]))), axis=1).astype("float64") # one row has (px, py, pyz, E)

    nu_array = np.zeros((jets_inputs.shape[0], 4))
    [pynusolver.run_nu_solver(lep_inputs[idx], jets_inputs[idx], met_inputs[idx], nu_array[idx]) for idx in range(nu_array.shape[0])]

    # Nu
        # convert px, py, pz to pt, eta, phi
    nu_px, nu_py, nu_pz = nu_array[:, 0], nu_array[:, 1], nu_array[:, 2]
    nu_mom, nu_pt = np.sqrt(np.square(nu_px)+np.square(nu_py)+np.square(nu_pz)), np.sqrt(np.square(nu_px)+np.square(nu_py))
    nu_phi = np.arctan2(nu_py, nu_px)
    nu_eta = np.arcsinh(nu_pz/nu_pt)
    Nu = ak.Array({
        "pt" : ak.unflatten(nu_pt, ak.num(jets)),
        "eta" : ak.unflatten(nu_eta, ak.num(jets)),
        "phi" : ak.unflatten(nu_phi, ak.num(jets)),
        "mass" : ak.zeros_like(ak.unflatten(nu_phi, ak.num(jets))),
        "chi2" : ak.unflatten(nu_array[:, 3], ak.num(jets)),
    }, with_name="PtEtaPhiMLorentzVector")

    if ak.num(jets)[0] == 3:
        bp_ordering, bp_probs = get_best_3j_permutations(jets=jets, nschi=np.sqrt(Nu.chi2), jets_btagging=jets_inputs[:,4], btag_req=btag_req)
    else:
        bp_ordering, bp_probs = get_best_4pj_permutations(jets=jets, nschi=np.sqrt(Nu.chi2), jets_btagging=jets_inputs[:,4], btag_req=btag_req)
    #set_trace()


        # convert lists into awkward arrays
    bp_probs = ak.from_iter(bp_probs)
    bp_ordering = ak.from_iter(bp_ordering)

        ## only keep permutations with some sort of solution (prob != infinity)
    valid_evts = ak.to_numpy(bp_probs[:,0] != np.inf)

    # BLep
    blep_inds = ak.unflatten(bp_ordering[valid_evts][:, 0], valid_evts.astype(int))
    best_BLep = ak.Array({
        "pt" : jets[blep_inds].pt,
        "eta": jets[blep_inds].eta,
        "phi": jets[blep_inds].phi,
        "mass": jets[blep_inds].mass,
        "jetIdx": blep_inds,
    }, with_name="PtEtaPhiMLorentzVector")

    # BHad
    bhad_inds = ak.unflatten(bp_ordering[valid_evts][:, 1], valid_evts.astype(int))
    best_BHad = ak.Array({
        "pt" : jets[bhad_inds].pt,
        "eta": jets[bhad_inds].eta,
        "phi": jets[bhad_inds].phi,
        "mass": jets[bhad_inds].mass,
        "jetIdx": bhad_inds,
    }, with_name="PtEtaPhiMLorentzVector")

    # WJa
    wja_inds = ak.unflatten(bp_ordering[valid_evts][:, 2], valid_evts.astype(int))
    best_WJa = ak.Array({
        "pt" : jets[wja_inds].pt,
        "eta": jets[wja_inds].eta,
        "phi": jets[wja_inds].phi,
        "mass": jets[wja_inds].mass,
        "jetIdx": wja_inds,
    }, with_name="PtEtaPhiMLorentzVector")

    # WJb
    if len(bp_ordering[0]) == 4: # WJb exists
        wjb_inds = ak.unflatten(bp_ordering[valid_evts][:, 3], valid_evts.astype(int))
        best_WJb = ak.Array({
            "pt" : jets[wjb_inds].pt,
            "eta": jets[wjb_inds].eta,
            "phi": jets[wjb_inds].phi,
            "mass": jets[wjb_inds].mass,
            "jetIdx": wjb_inds,
        }, with_name="PtEtaPhiMLorentzVector")
    else:
        best_WJb = ak.Array({
            "pt" : ak.zeros_like(best_WJa.pt),
            "eta" : ak.zeros_like(best_WJa.eta),
            "phi" : ak.zeros_like(best_WJa.phi),
            "mass" : ak.zeros_like(best_WJa.mass),
            "jetIdx": -10*ak.ones_like(best_WJa.jetIdx),
        }, with_name="PtEtaPhiMLorentzVector")

    # Nu
    best_Nu = Nu[blep_inds]

    # WHad
    best_WHad = ak.Array({
        "pt" : (best_WJa+best_WJb).pt,
        "eta" : (best_WJa+best_WJb).eta,
        "phi" : (best_WJa+best_WJb).phi,
        "mass" : (best_WJa+best_WJb).mass,
    }, with_name="PtEtaPhiMLorentzVector")

    # THad
    best_THad = ak.Array({
        "pt" : (best_BHad+best_WHad).pt,
        "eta" : (best_BHad+best_WHad).eta,
        "phi" : (best_BHad+best_WHad).phi,
        "mass" : (best_BHad+best_WHad).mass,
    }, with_name="PtEtaPhiMLorentzVector")

    # WLep
    best_WLep = ak.Array({
        "pt" : ak.flatten((leptons+best_Nu).pt, axis=1),
        "eta" : ak.flatten((leptons+best_Nu).eta, axis=1),
        "phi" : ak.flatten((leptons+best_Nu).phi, axis=1),
        "mass" : ak.flatten((leptons+best_Nu).mass, axis=1),
    }, with_name="PtEtaPhiMLorentzVector")

    # TLep
    best_TLep = ak.Array({
        "pt" : (best_BLep+best_WLep).pt,
        "eta" : (best_BLep+best_WLep).eta,
        "phi" : (best_BLep+best_WLep).phi,
        "mass" : (best_BLep+best_WLep).mass,
    }, with_name="PtEtaPhiMLorentzVector")

    # TTbar
    best_TTbar = ak.Array({
        "pt" : (best_TLep+best_THad).pt,
        "eta" : (best_TLep+best_THad).eta,
        "phi" : (best_TLep+best_THad).phi,
        "mass" : (best_TLep+best_THad).mass,
    }, with_name="PtEtaPhiMLorentzVector")

    # Lepton
    best_Lep = ak.Array({key: ak.unflatten(ak.flatten(leptons[key][valid_evts]), valid_evts.astype(int)) for key in leptons.fields}, with_name="PtEtaPhiMLorentzVector")

    # MET
    best_MET = ak.Array({key: ak.unflatten(MET[key][valid_evts], valid_evts.astype(int)) for key in MET.fields}, with_name="PtEtaPhiMLorentzVector")

    best_Top = ak.Array({key: ak.unflatten(ak.flatten(ak.where(np.sign(best_Lep.charge) == 1, best_TLep[key], best_THad[key])), valid_evts.astype(int)) for key in best_TLep.fields}, with_name="PtEtaPhiMLorentzVector")
    best_Tbar = ak.Array({key: ak.unflatten(ak.flatten(ak.where(np.sign(best_Lep.charge) == -1, best_TLep[key], best_THad[key])), valid_evts.astype(int)) for key in best_TLep.fields}, with_name="PtEtaPhiMLorentzVector")

        ## Combine everything into a dictionary
    best_permutations = ak.zip({
        "Top" : best_Top,
        "Tbar" : best_Tbar,
        "BHad" : best_BHad,
        "BLep" : best_BLep,
        "WJa" : best_WJa,
        "WJb" : best_WJb,
        "Lepton" : best_Lep,
        "MET" : best_MET,
        "Nu" : best_Nu,
        "WLep" : best_WLep,
        "TLep" : best_TLep,
        "WHad" : best_WHad,
        "THad" : best_THad,
        "TTbar": best_TTbar,
        "Prob" : ak.unflatten(bp_probs[valid_evts][:, 0], valid_evts.astype(int)),
        "MassDiscr" : ak.unflatten(bp_probs[valid_evts][:, 1], valid_evts.astype(int)),
        "NuDiscr" : ak.unflatten(bp_probs[valid_evts][:, 2], valid_evts.astype(int)),
    }, depth_limit=1)

    return best_permutations
