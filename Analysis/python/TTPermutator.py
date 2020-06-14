from numba import njit, objmode
import numpy as np
from pdb import set_trace
import python.TTBarSolver as ttsolver
import compiled.pynusolver as pynusolver
import awkward
from coffea.analysis_objects import JaggedCandidateArray
from python.Permutations import make_perm_table

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
def get_permutations(njets_array, jets, leptons, met, use_merged=False):#, btag_req=True)
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


#@njit()
def find_best_permutations(jets, leptons, MET, btagWP):#, btag_req=True):
    '''
    Inputs:
        Jets, leptons, MET, and btag working point for jets
    Returns:
        awkward Table containing
            1: Jet objects (BLeps, BHads, WJas, WJbs)
            2: Leptons, Neutrinos
            3: Calculated probabilities (Total -> Prob, mass discriminant -> MassDiscr, neutrino discrminant -> NuDiscr)
    '''

    jets_inputs = np.stack((jets.p4.x.flatten(), jets.p4.y.flatten(), jets.p4.z.flatten(), jets.p4.energy.flatten(), jets[btagWP].flatten()), axis=1).astype('float64') # one row has (px, py, pyz, E)
    leptons_inputs = np.stack((leptons.p4.x.flatten(), leptons.p4.y.flatten(), leptons.p4.z.flatten(), leptons.p4.energy.flatten()), axis=1).astype('float64') # one row has (px, py, pyz, E)
    met_inputs = np.stack((MET.p4.x.flatten(), MET.p4.y.flatten()), axis=1).astype('float64') # one row has (px, py)
    bp_ordering, bp_nus, bp_probs = get_permutations(njets_array=jets.counts, jets=jets_inputs, leptons=leptons_inputs, met=met_inputs)#, btag_req=btag_req)
    ## for testing
    #nj = jets.counts[0:4]
    #bp_ordering, bp_nus, bp_probs = get_permutations(njets_array=nj, jets=jets_inputs[0:sum(nj)], leptons=leptons_inputs[0:len(nj)], met=met_inputs[0:len(nj)])
    ##

    bp_ordering = np.asarray(bp_ordering)
    bp_nus = np.asarray(bp_nus)
    bp_probs = np.asarray(bp_probs)

    #set_trace()
        ## only keep permutations with some sort of solution (prob != infinity)
    valid_evts = (bp_probs[:,0] != np.inf)

        ## BLep
    blep_inds = bp_ordering[:,0]
    best_BLep = JaggedCandidateArray.candidatesfromcounts(
        counts=valid_evts.astype(int),
        px    = jets[np.arange(len(blep_inds))[valid_evts], blep_inds[valid_evts]].p4.x.flatten(),
        py    = jets[np.arange(len(blep_inds))[valid_evts], blep_inds[valid_evts]].p4.y.flatten(),
        pz    = jets[np.arange(len(blep_inds))[valid_evts], blep_inds[valid_evts]].p4.z.flatten(),
        energy= jets[np.arange(len(blep_inds))[valid_evts], blep_inds[valid_evts]].p4.energy.flatten(),
        jetIdx=blep_inds[valid_evts],
    )

        ## BHad
    bhad_inds = bp_ordering[:,1]
    best_BHad = JaggedCandidateArray.candidatesfromcounts(
        counts=valid_evts.astype(int),
        px    = jets[np.arange(len(bhad_inds))[valid_evts], bhad_inds[valid_evts]].p4.x.flatten(),
        py    = jets[np.arange(len(bhad_inds))[valid_evts], bhad_inds[valid_evts]].p4.y.flatten(),
        pz    = jets[np.arange(len(bhad_inds))[valid_evts], bhad_inds[valid_evts]].p4.z.flatten(),
        energy= jets[np.arange(len(bhad_inds))[valid_evts], bhad_inds[valid_evts]].p4.energy.flatten(),
        jetIdx=bhad_inds[valid_evts],
    )

        ## WJa
    wja_inds = bp_ordering[:,2]
    best_WJa = JaggedCandidateArray.candidatesfromcounts(
        counts=valid_evts.astype(int),
        px    = jets[np.arange(len(wja_inds))[valid_evts], wja_inds[valid_evts]].p4.x.flatten(),
        py    = jets[np.arange(len(wja_inds))[valid_evts], wja_inds[valid_evts]].p4.y.flatten(),
        pz    = jets[np.arange(len(wja_inds))[valid_evts], wja_inds[valid_evts]].p4.z.flatten(),
        energy= jets[np.arange(len(wja_inds))[valid_evts], wja_inds[valid_evts]].p4.energy.flatten(),
        jetIdx=wja_inds[valid_evts],
    )

        ## WJb
    if bp_ordering.shape[1] == 4:
        wjb_inds = bp_ordering[:,3]
        best_WJb = JaggedCandidateArray.candidatesfromcounts(
            counts=valid_evts.astype(int),
            px    = jets[np.arange(len(wjb_inds))[valid_evts], wjb_inds[valid_evts]].p4.x.flatten(),
            py    = jets[np.arange(len(wjb_inds))[valid_evts], wjb_inds[valid_evts]].p4.y.flatten(),
            pz    = jets[np.arange(len(wjb_inds))[valid_evts], wjb_inds[valid_evts]].p4.z.flatten(),
            energy= jets[np.arange(len(wjb_inds))[valid_evts], wjb_inds[valid_evts]].p4.energy.flatten(),
            jetIdx=wjb_inds[valid_evts],
        )
    else:
        best_WJb = JaggedCandidateArray.candidatesfromcounts(
            counts=valid_evts.astype(int),
            px    = np.zeros(valid_evts.sum()),
            py    = np.zeros(valid_evts.sum()),
            pz    = np.zeros(valid_evts.sum()),
            energy= np.zeros(valid_evts.sum()),
            jetIdx= np.ones(valid_evts.sum())*(-10),
        )

        ## Nu
    best_Nu = JaggedCandidateArray.candidatesfromcounts(
        counts = valid_evts.astype(int),
        px = bp_nus[:, 0][valid_evts],
        py = bp_nus[:, 1][valid_evts],
        pz = bp_nus[:, 2][valid_evts],
        mass = np.zeros(valid_evts.sum()),
        chi2 = bp_nus[:, 3][valid_evts],
    )

        ## WHad
    whad_p4 = best_WJa.p4.flatten() + best_WJb.p4.flatten()
    best_WHad = JaggedCandidateArray.candidatesfromcounts(
        counts=valid_evts.astype(int),
        px    = whad_p4.x,
        py    = whad_p4.y,
        pz    = whad_p4.z,
        energy= whad_p4.energy,
    )

        ## THad
    thad_p4 = best_BHad.p4.flatten() + best_WHad.p4.flatten()
    best_THad = JaggedCandidateArray.candidatesfromcounts(
        counts=valid_evts.astype(int),
        px    = thad_p4.x,
        py    = thad_p4.y,
        pz    = thad_p4.z,
        energy= thad_p4.energy,
    )

        ## WLep
    wlep_p4 = leptons[valid_evts].p4.flatten() + best_Nu.p4.flatten()
    best_WLep = JaggedCandidateArray.candidatesfromcounts(
        counts=valid_evts.astype(int),
        px    = wlep_p4.x,
        py    = wlep_p4.y,
        pz    = wlep_p4.z,
        energy= wlep_p4.energy,
    )

        ## TLep
    tlep_p4 = best_BLep.p4.flatten() + best_WLep.p4.flatten()
    best_TLep = JaggedCandidateArray.candidatesfromcounts(
        counts=valid_evts.astype(int),
        px    = tlep_p4.x,
        py    = tlep_p4.y,
        pz    = tlep_p4.z,
        energy= tlep_p4.energy,
    )

        ## TTbar
    tt_p4 = best_THad.p4.flatten() + best_TLep.p4.flatten()
    best_TTbar = JaggedCandidateArray.candidatesfromcounts(
        counts=valid_evts.astype(int),
        px    = tt_p4.x,
        py    = tt_p4.y,
        pz    = tt_p4.z,
        energy= tt_p4.energy,
    )

    remove_fast = lambda x : x.split('fast_')[-1]
        ## Lepton
    dict_vars = {'Lepton' : {remove_fast(key) : leptons[valid_evts][key].flatten() for key in leptons.columns if key != 'p4'}}
    best_Lep = JaggedCandidateArray.candidatesfromcounts(
        counts=valid_evts.astype(int),
        **dict_vars['Lepton']
    )

        ## MET
    dict_vars.update({'MET' : {remove_fast(key) : MET[valid_evts][key].flatten() for key in MET.columns if key != 'p4'}})
    best_MET = JaggedCandidateArray.candidatesfromcounts(
        counts=valid_evts.astype(int),
        **dict_vars['MET']
    )

        ## Combine everything into a single table, all objects are JaggedArrays
    best_permutations = awkward.Table(
        BHad = best_BHad,
        BLep = best_BLep,
        WJa = best_WJa,
        WJb = best_WJb,
        Lepton = best_Lep,
        MET = best_MET,
        Nu = best_Nu,
        WLep = best_WLep,
        TLep = best_TLep,
        WHad = best_WHad,
        THad = best_THad,
        TTbar= best_TTbar,
        Prob = awkward.JaggedArray.fromcounts(valid_evts.astype(int), bp_probs[valid_evts][:,0]),
        MassDiscr = awkward.JaggedArray.fromcounts(valid_evts.astype(int), bp_probs[valid_evts][:,1]),
        NuDiscr = awkward.JaggedArray.fromcounts(valid_evts.astype(int), bp_probs[valid_evts][:,2]),
    )

    #set_trace()
    return best_permutations

