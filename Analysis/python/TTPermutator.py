from numba import njit, objmode
import numpy as np
from pdb import set_trace
import python.TTBarSolver as ttsolver
import compiled.pynusolver as pynusolver
import awkward as ak

solver = None
def year_to_run(**kwargs):
    global solver
    solver = ttsolver.TTSolver(kwargs['year'])

def run():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument('year', choices=['2016', '2017', '2018'], help='Specify which year to run over')
    args = parser.parse_args()
    year_to_run(**kwargs)

@njit()
def get_permutations(njets_array, jets, leptons, met, use_merged=False, btag_req=True):
    '''
        Inputs:
            1D array with number of jets per event
            2D array of jets (px, py, pz, E, btag pass)
            2D array of leptons (px, py, pz, E)
            2D array of MET (px, py)
        Returns:
            [0]: List of jet assignment ordering
            [1]: List of neutrino solutions
            [2]: List of probabilities (Total, mass discriminant, nu discriminant)
    '''
    #print('Finding best perms')
    start = 0
    stop = 0
    evt_idx = 0

    best_perms_ordering = []
    best_perms_nus = []
    best_perms_probs = []

    #set_trace()
    for njets in njets_array:
        #print('evt idx: ', evt_idx)#, ', njets: ', njets)
        stop += njets

        ## initialize best_perm lists for event
        best_perm_ordering = np.array([-10, -10, -10]) if njets == 3 else np.array([-10, -10, -10, -10])
        best_perm_nu = np.zeros(4)
        best_perm_probs = np.array([np.inf, np.inf, np.inf])

        for j0 in range(start, stop):
                # first jet has to be btagged
            if btag_req:
                if jets[j0, 4] < 0.5: continue
            ## run NS to get nschi2 to be used solver
            nu = np.zeros(4)
            with objmode():
                pynusolver.run_nu_solver(leptons[evt_idx], jets[j0], met[evt_idx], nu) # input order for nusolver is (lepton, jet, met, nu)

            ## lepton and met info won't change inside njets for loop
            ## jet info will change for each value of j0
                ## get best perms for 3 jet events
            if njets == 3:
                #print('3 jets')
                for j1 in range(start, stop):
                    if j1 == j0: continue
                        # second jet has to be btagged
                    if btag_req:
                        if jets[j1, 4] < 0.5: continue
                    for j2 in range(j1+1, stop):
                        if j2 == j0: continue
                        #    # btagged jets can't be assigned to wjets
                        #if jets[j2, 4] > 0.5: continue
                            ## advanced indexing doesn't work with numba at this point, can only add px, py, pz separately
                        mbpjet = np.sqrt( np.sum(jets[:, 3].take([j2, j1]))**2 - (np.sum(jets[:, 0].take([j2, j1]))**2 + np.sum(jets[:, 1].take([j2, j1]))**2 + np.sum(jets[:, 2].take([j2, j1]))**2) ) ## sqrt(E2-p2) of combined j1+j2 4-vector
                        disc_probs = np.zeros(3)
                        if use_merged:
                            maxmjet = max( np.sqrt(jets[j1, 3]**2 - np.sum(jets[j1, 0:3]**2)), np.sqrt(jets[j2, 3]**2 - np.sum(jets[j2, 0:3]**2)) ) ## get max mass between j1 and j2
                            with objmode():
                                probs = solver.solve_3J_merged(maxmjet, mbpjet, nu[3]) ## input order is (maxmjet, mbpjet, nschi2 value)
                                disc_probs[0], disc_probs[1], disc_probs[2] = probs[0], probs[1], probs[2]
                            #print('max(mjet) = ', maxmjet, ', m(b+j) = ', mbpjet, ', nschi = ', nu[3], ', Prob = ', prob, ', MassDiscr = ', massdiscr, ', NuDiscr = ', nudiscr)
                        else:
                            with objmode():
                                probs = solver.solve_3J_lost(mbpjet, nu[3]) ## input order is (mbpjet, nschi2 value)
                                disc_probs[0], disc_probs[1], disc_probs[2] = probs[0], probs[1], probs[2]
                            #print('   end:', j0-start, j1-start, j2-start)
                            #print('   m(b+j) = ', mbpjet, ', nschi = ', nu[3], ', Prob = ', prob, ', MassDiscr = ', massdiscr, ', NuDiscr = ', nudiscr)
                        if disc_probs[0] < best_perm_probs[0]: ## have to be careful with indexing!!
                            best_perm_ordering = np.array([j0-start, j1-start, j2-start])#, -10]) # set last element to -10 in order to make array have the same shape as 4+ jets
                            best_perm_nu = nu
                            best_perm_probs = disc_probs
                            #best_perm_probs = np.array([prob, massdiscr, nudiscr])


                ## get best perms for 4+ jet events
            else:
                #print('4+ jets')
                for j1 in range(start, stop):
                    if j1 == j0: continue
                        # second jet has to be btagged
                    if btag_req:
                        if jets[j1, 4] < 0.5: continue
                    for j2 in range(start, stop):
                        if j2 == j0 or j2 == j1: continue
                        #    # btagged jets can't be assigned to wjets
                        #if jets[j2, 4] > 0.5: continue
                        for j3 in range(j2+1, stop):
                            if j3 == j0 or j3 == j1: continue
                                ## advanced indexing doesn't work with numba at this point, can only add px, py, pz separately
                            mthad = np.sqrt( np.sum(jets[:, 3].take([j1, j2, j3]))**2 - (np.sum(jets[:, 0].take([j1, j2, j3]))**2 + np.sum(jets[:, 1].take([j1, j2, j3]))**2 + np.sum(jets[:, 2].take([j1, j2, j3]))**2) ) ## sqrt(E2-p2) of combined j1+j2+j3 4-vector
                            mwhad = np.sqrt( np.sum(jets[:, 3].take([j2, j3]))**2 - (np.sum(jets[:, 0].take([j2, j3]))**2 + np.sum(jets[:, 1].take([j2, j3]))**2 + np.sum(jets[:, 2].take([j2, j3]))**2) ) ## sqrt(E2-p2) of combined j2+j3 4-vector
                            disc_probs = np.zeros(3)
                            with objmode():
                                probs = solver.solve_4PJ(mthad, mwhad, nu[3]) ## input order is (mthad, mwhad, nschi2 value)
                                disc_probs[0], disc_probs[1], disc_probs[2] = probs[0], probs[1], probs[2]
                            #print('   end:', j0-start, j1-start, j2-start, j3-start)
                            #print('   mthad = ', mthad, ', mwhad = ', mwhad, ', nschi = ', nu[3], ', Prob = ', prob, ', MassDiscr = ', massdiscr, ', NuDiscr = ', nudiscr)
                            if disc_probs[0] < best_perm_probs[0]: ## have to be careful with indexing!!
                                best_perm_ordering = np.array([j0-start, j1-start, j2-start, j3-start])
                                best_perm_nu = nu
                                best_perm_probs = disc_probs
                                #best_perm_probs = np.array([prob, massdiscr, nudiscr])

        #print('    Ordering:', best_perm_ordering)
        #print('Ordering:', best_perm_ordering, ', nu:', best_perm_nu, ', probs:', best_perm_probs)
        best_perms_ordering.append(best_perm_ordering)
        best_perms_nus.append(best_perm_nu)
        best_perms_probs.append(best_perm_probs)

        start += njets
        #if np.mod(evt_idx, 10) == 0:
        #    print('\t', evt_idx+1,'/',len(njets_array), 'events processed')
        evt_idx += 1
    #print('Finished finding best perms')
    return best_perms_ordering, best_perms_nus, best_perms_probs


def find_best_permutations(jets, leptons, MET, btagWP, btag_req=True):
    '''
    Inputs:
        Jets, leptons, MET, and btag working point for jets
    Returns:
        awkward Table containing
            1: Jet objects (BLeps, BHads, WJas, WJbs)
            2: Leptons, Neutrinos
            3: Calculated probabilities (Total -> Prob, mass discriminant -> MassDiscr, neutrino discrminant -> NuDiscr)
    '''

    jets_inputs = np.stack((ak.to_numpy(ak.flatten(jets.px)), ak.to_numpy(ak.flatten(jets.py)), ak.to_numpy(ak.flatten(jets.pz)), ak.to_numpy(ak.flatten(jets.energy)), ak.to_numpy(ak.flatten(jets[btagWP]))), axis=1).astype('float64') # one row has (px, py, pyz, E)
    leptons_inputs = np.stack((ak.to_numpy(ak.flatten(leptons.px)), ak.to_numpy(ak.flatten(leptons.py)), ak.to_numpy(ak.flatten(leptons.pz)), ak.to_numpy(ak.flatten(leptons.energy))), axis=1).astype('float64') # one row has (px, py, pyz, E)
    met_inputs = np.stack((ak.to_numpy(MET.px), ak.to_numpy(MET.py)), axis=1).astype('float64') # one row has (px, py)
    bp_ordering, bp_nus, bp_probs = get_permutations(njets_array=ak.num(jets), jets=jets_inputs, leptons=leptons_inputs, met=met_inputs, btag_req=btag_req)
    ## for testing
    #nj = jets.counts[0:4]
    #bp_ordering, bp_nus, bp_probs = get_permutations(njets_array=nj, jets=jets_inputs[0:sum(nj)], leptons=leptons_inputs[0:len(nj)], met=met_inputs[0:len(nj)])
    ##

        # convert lists into awkward arrays
    bp_probs = ak.from_iter(bp_probs)
    bp_nus = ak.from_iter(bp_nus)
    bp_ordering = ak.from_iter(bp_ordering)

        ## only keep permutations with some sort of solution (prob != infinity)
    valid_evts = ak.to_numpy(bp_probs[:,0] != np.inf)

    # BLep
    blep_inds = ak.unflatten(bp_ordering[valid_evts][:, 0], valid_evts.astype(int))
    best_BLep = ak.Array({
        'pt' : jets[blep_inds].pt,
        'eta': jets[blep_inds].eta,
        'phi': jets[blep_inds].phi,
        'mass': jets[blep_inds].mass,
        'jetIdx': blep_inds,
    }, with_name="PtEtaPhiMLorentzVector")

    # BHad
    bhad_inds = ak.unflatten(bp_ordering[valid_evts][:, 1], valid_evts.astype(int))
    best_BHad = ak.Array({
        'pt' : jets[bhad_inds].pt,
        'eta': jets[bhad_inds].eta,
        'phi': jets[bhad_inds].phi,
        'mass': jets[bhad_inds].mass,
        'jetIdx': bhad_inds,
    }, with_name="PtEtaPhiMLorentzVector")

    # WJa
    wja_inds = ak.unflatten(bp_ordering[valid_evts][:, 2], valid_evts.astype(int))
    best_WJa = ak.Array({
        'pt' : jets[wja_inds].pt,
        'eta': jets[wja_inds].eta,
        'phi': jets[wja_inds].phi,
        'mass': jets[wja_inds].mass,
        'jetIdx': wja_inds,
    }, with_name="PtEtaPhiMLorentzVector")

    # WJb
    if len(bp_ordering[0]) == 4: # WJb exists
        wjb_inds = ak.unflatten(bp_ordering[valid_evts][:, 3], valid_evts.astype(int))
        best_WJb = ak.Array({
            'pt' : jets[wjb_inds].pt,
            'eta': jets[wjb_inds].eta,
            'phi': jets[wjb_inds].phi,
            'mass': jets[wjb_inds].mass,
            'jetIdx': wjb_inds,
        }, with_name="PtEtaPhiMLorentzVector")
    else:
        best_WJb = ak.Array({
            'pt' : ak.zeros_like(best_WJa.pt),
            'eta' : ak.zeros_like(best_WJa.eta),
            'phi' : ak.zeros_like(best_WJa.phi),
            'mass' : ak.zeros_like(best_WJa.mass),
            'jetIdx': -10*ak.ones_like(best_WJa.jetIdx),
        }, with_name="PtEtaPhiMLorentzVector")

    # Nu
        # convert px, py, pz to pt, eta, phi
    nu_px, nu_py, nu_pz = bp_nus[:, 0][valid_evts], bp_nus[:, 1][valid_evts], bp_nus[:, 2][valid_evts]
    nu_mom, nu_pt = np.sqrt(np.square(nu_px)+np.square(nu_py)+np.square(nu_pz)), np.sqrt(np.square(nu_px)+np.square(nu_py))
    nu_phi = np.arctan2(nu_py, nu_px)
    nu_eta = np.arcsinh(nu_pz/nu_pt)
    best_Nu = ak.Array({
        'pt' : ak.unflatten(nu_pt, valid_evts.astype(int)),
        'eta' : ak.unflatten(nu_eta, valid_evts.astype(int)),
        'phi' : ak.unflatten(nu_phi, valid_evts.astype(int)),
        'mass' : ak.zeros_like(ak.unflatten(bp_nus[:, 0][valid_evts], valid_evts.astype(int))),
        'chi2' : ak.unflatten(bp_nus[:, 3][valid_evts], valid_evts.astype(int)),
    }, with_name="PtEtaPhiMLorentzVector")

    # WHad
    best_WHad = ak.Array({
        'pt' : (best_WJa+best_WJb).pt,
        'eta' : (best_WJa+best_WJb).eta,
        'phi' : (best_WJa+best_WJb).phi,
        'mass' : (best_WJa+best_WJb).mass,
    }, with_name="PtEtaPhiMLorentzVector")

    # THad
    best_THad = ak.Array({
        'pt' : (best_BHad+best_WHad).pt,
        'eta' : (best_BHad+best_WHad).eta,
        'phi' : (best_BHad+best_WHad).phi,
        'mass' : (best_BHad+best_WHad).mass,
    }, with_name="PtEtaPhiMLorentzVector")

    # WLep
    best_WLep = ak.Array({
        'pt' : ak.flatten((leptons+best_Nu).pt, axis=1),
        'eta' : ak.flatten((leptons+best_Nu).eta, axis=1),
        'phi' : ak.flatten((leptons+best_Nu).phi, axis=1),
        'mass' : ak.flatten((leptons+best_Nu).mass, axis=1),
    }, with_name="PtEtaPhiMLorentzVector")

    # TLep
    best_TLep = ak.Array({
        'pt' : (best_BLep+best_WLep).pt,
        'eta' : (best_BLep+best_WLep).eta,
        'phi' : (best_BLep+best_WLep).phi,
        'mass' : (best_BLep+best_WLep).mass,
    }, with_name="PtEtaPhiMLorentzVector")

    # TTbar
    best_TTbar = ak.Array({
        'pt' : (best_TLep+best_THad).pt,
        'eta' : (best_TLep+best_THad).eta,
        'phi' : (best_TLep+best_THad).phi,
        'mass' : (best_TLep+best_THad).mass,
    }, with_name="PtEtaPhiMLorentzVector")

    # Lepton
    best_Lep = ak.Array({key: ak.unflatten(ak.flatten(leptons[key][valid_evts]), valid_evts.astype(int)) for key in leptons.fields}, with_name="PtEtaPhiMLorentzVector")

    # MET
    best_MET = ak.Array({key: ak.unflatten(MET[key][valid_evts], valid_evts.astype(int)) for key in MET.fields}, with_name="PtEtaPhiMLorentzVector")

        ## Combine everything into a dictionary
    best_permutations = ak.zip({
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
