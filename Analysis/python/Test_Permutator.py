from numba import njit, objmode
import numpy as np
from pdb import set_trace
import compiled.pynusolver as pynusolver

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

    perms_ordering_nu = []

    #set_trace()
    for njets in njets_array:
        #print('evt idx: ', evt_idx, ', njets: ', njets)
        #print('evt idx: ', evt_idx)
        stop += njets

        ### initialize best_perm lists for event
        evt_ordering_nu = []

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
                        evt_ordering_nu.append(([j0-start, j1-start, j2-start], nu))

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
                            evt_ordering_nu.append(([j0-start, j1-start, j2-start, j3-start], nu))

        perms_ordering_nu.append(evt_ordering_nu)

        start += njets
        evt_idx += 1
    return perms_ordering_nu

#@njit()
def find_permutations(jets, leptons, MET, btagWP):
    '''
    Inputs:
        Jets, leptons, MET, and event weights
    Returns:
        awkward Table containing
            1: Jet objects (BLeps, BHads, WJas, WJbs)
            2: Leptons, Neutrinos
            3: Calculated probabilities (Total -> Prob, mass discriminant -> MassDiscr, neutrino discrminant -> NuDiscr)
            4: Number of jets in events (njets)
            5: Event weights (evt_wts)
    '''

    jets_inputs = np.stack((jets.p4.x.flatten(), jets.p4.y.flatten(), jets.p4.z.flatten(), jets.p4.energy.flatten(), jets[btagWP].flatten()), axis=1).astype('float64') # one row has (px, py, pyz, E)
    lepton_inputs = np.stack((leptons.p4.x.flatten(), leptons.p4.y.flatten(), leptons.p4.z.flatten(), leptons.p4.energy.flatten()), axis=1).astype('float64') # one row has (px, py, pyz, E)
    met_inputs = np.stack((MET.p4.x.flatten(), MET.p4.y.flatten()), axis=1).astype('float64') # one row has (px, py)
    p_ordering_nu = get_test_permutations(njets_array=jets.counts, jets=jets_inputs, leptons=lepton_inputs, met=met_inputs)
    #set_trace()

    return np.array(p_ordering_nu).astype('O')
