from numba import njit, objmode
import numpy as np
from pdb import set_trace
import compiled.pynusolver as pynusolver
import awkward as ak

@njit()
def get_test_permutations(njets_array, jets, leptons, met, btag_req=True):
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

    perms_ordering = []
    perms_nu = []

    #set_trace()
    for njets in njets_array:
        #print('evt idx: ', evt_idx, ', njets: ', njets)
        #print('evt idx: ', evt_idx)
        stop += njets

        ### initialize best_perm lists for event
        evt_ordering = []
        evt_nu = []

        for j0 in range(start, stop):
            ## require btagging
            if btag_req:
                if jets[j0, 4] < 0.5: continue

            ## run NS to get nschi2 to be used solver
            nu = np.zeros(4)
            with objmode():
                pynusolver.run_nu_solver(leptons[evt_idx], jets[j0], met[evt_idx], nu) # input order for nusolver is (lepton, jet, met, nu)
            #print("nu: ", nu)

            ## lepton and met info won't change inside njets for loop
            ## jet info will change for each value of j0
                ## get best perms for 3 jet events
            if njets == 3:
                for j1 in range(start, stop):
                    if j1 == j0: continue
                    if btag_req:
                        if jets[j1, 4] < 0.5: continue
                    for j2 in range(j1+1, stop):
                        if j2 == j0: continue
                        evt_ordering.append([j0-start, j1-start, j2-start, -999])
                        evt_nu.append(nu)

                ## get best perms for 4+ jet events
            else:
                for j1 in range(start, stop):
                    if j1 == j0: continue
                    if btag_req:
                        if jets[j1, 4] < 0.5: continue
                    for j2 in range(start, stop):
                        if j2 == j0 or j2 == j1: continue
                        for j3 in range(j2+1, stop):
                            if j3 == j0 or j3 == j1: continue
                            evt_ordering.append([j0-start, j1-start, j2-start, j3-start])
                            evt_nu.append(nu)

        perms_ordering.append(evt_ordering)
        perms_nu.append(evt_nu)

        start += njets
        evt_idx += 1
    return perms_ordering, perms_nu

#@njit()
def find_permutations(jets, leptons, MET, btagWP):
    '''
    Inputs:
        Jets, leptons, MET, and if jets pass btag WP
    Returns:
        List of (jet assignment ordering, associated neutrino solutions)
    '''

    jets_inputs = np.stack((ak.to_numpy(ak.ak.flatten(jets.px)), ak.to_numpy(ak.flatten(jets.py)), ak.to_numpy(ak.flatten(jets.pz)), ak.to_numpy(ak.flatten(jets.energy)), ak.to_numpy(ak.flatten(jets[btagWP]))), axis=1).astype('float64') # one row has (px, py, pyz, E)
    lepton_inputs = np.stack((ak.to_numpy(ak.flatten(leptons.px)), ak.to_numpy(ak.flatten(leptons.py)), ak.to_numpy(ak.flatten(leptons.pz)), ak.to_numpy(ak.flatten(leptons.energy))), axis=1).astype('float64') # one row has (px, py, pyz, E)
    met_inputs = np.stack((ak.to_numpy(MET.px), ak.to_numpy(MET.py)), axis=1).astype('float64') # one row has (px, py)
    p_ordering, p_nu = get_test_permutations(njets_array=ak.num(jets), jets=jets_inputs, leptons=lepton_inputs, met=met_inputs)

    #set_trace()
    test_perms = ak.Array({
        'blepIdx' : ak.from_iter(p_ordering)[:, :, 0],
        'bhadIdx' : ak.from_iter(p_ordering)[:, :, 1],
        'wjaIdx' : ak.from_iter(p_ordering)[:, :, 2],
        'wjbIdx' : ak.from_iter(p_ordering)[:, :, 3],
        'Nu' : ak.Array({
            'px' : ak.from_iter(p_nu)[:, :, 0],
            'py' : ak.from_iter(p_nu)[:, :, 1],
            'pz' : ak.from_iter(p_nu)[:, :, 2],
            'chi2' : ak.from_iter(p_nu)[:, :, 3],
        })
    })
    return test_perms
