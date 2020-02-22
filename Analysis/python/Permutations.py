from numba import njit
import numpy as np
from pdb import set_trace
import python.TTBarSolver as solver


@njit()
def get_permutations_4j(njets_array, jets, leptons, met):
    'Inputs:\n\t1D array with number of jets per event\n\t2D array of jets (px, py, pz, E, btag pass)\n\t2D array of leptons (px, py, pz, E)\n\t2D array of MET (px, py)'
    start = 0
    stop = 0
    evt_idx = 0

    for njets in njets_array:
        stop += njets
        for j0 in range(start, stop):
            ## require btagging
            if jets[j0, 4] < 0.5: continue

            ## run NS to get nschi2 to be used solver
            #run_nu_solver(leptons[evt_idx], jets[j0], met[evt_idx]) ## passes all info for only one lepton (px, py, pz, E), one jet (px, py, pz, E), and met (px, py) from the event at a time

            ## lepton and met info won't change inside njets for loop
            ## jet info will change for each value of j0 (4 separate for event with 4 jets)
            for j1 in range(start, stop):
                if j1 == j0: continue
                for j2 in range(start, stop):
                    if j2 == j0 or j2 == j1: continue
                    for j3 in range(j2+1, stop):
                        if j3 == j0 or j3 == j1: continue
                        #set_trace()
                            ## advanced indexing doesn't work with numba at this point, can only add px, py, pz separately
                        mthad = np.sqrt( np.sum(jets[:, 3].take([j1, j2, j3]))**2 - (np.sum(jets[:, 0].take([j1, j2, j3]))**2 + np.sum(jets[:, 1].take([j1, j2, j3]))**2 + np.sum(jets[:, 2].take([j1, j2, j3]))**2) ) ## sqrt(E2-p2) of combined j1+j2+j3 4-vector
                        mwhad = np.sqrt( np.sum(jets[:, 3].take([j2, j3]))**2 - (np.sum(jets[:, 0].take([j2, j3]))**2 + np.sum(jets[:, 1].take([j2, j3]))**2 + np.sum(jets[:, 2].take([j2, j3]))**2) ) ## sqrt(E2-p2) of combined j2+j3 4-vector
                        print('mthad = ', mthad, ', mwhad = ', mwhad)
                        #nschi = 400. # test value
                        #nudiscr, massdiscr, prob = solver.solve_4PJ(mthad, mwhad, nschi)
                        #print('mthad = ', mthad, ', mwhad = ', mwhad, ', nschi = ', nschi, ', Prob = ', prob, ', MassDiscr = ', massdiscr, ', NuDiscr = ', nudiscr)
                        print(j0, j1, j2, j3)

        start += njets
        evt_idx += 1
        print()

@njit()
def get_permutations_3j(njets_array, jets, leptons, met, use_merged=False):
    'Inputs:\n\t1D array with number of jets per event\n\t2D array of jets (px, py, pz, E, btag pass)\n\t2D array of leptons (px, py, pz, E)\n\t2D array of MET (px, py)'
    start = 0
    stop = 0
    evt_idx = 0

    #set_trace()
    for njets in njets_array:
        stop += njets
        for j0 in range(start, stop):
            ## require btagging
            if jets[j0, 4] < 0.5: continue

            ## run NS to get nschi2 to be used solver
            #run_nu_solver(leptons[evt_idx], jets[j0], met[evt_idx]) ## passes all info for only one lepton (px, py, pz, E), one jet (px, py, pz, E), and met (px, py) from the event at a time

            ## lepton and met info won't change inside njets for loop
            ## jet info will change for each value of j0 (3 separate for event with 3 jets)
            for j1 in range(start, stop):
                if j1 == j0: continue
                for j2 in range(j1+1, stop):
                    if j2 == j0: continue
                    #set_trace()
                        ## advanced indexing doesn't work with numba at this point, can only add px, py, pz separately
                    mbpjet = np.sqrt( np.sum(jets[:, 3].take([j2, j1]))**2 - (np.sum(jets[:, 0].take([j2, j1]))**2 + np.sum(jets[:, 1].take([j2, j1]))**2 + np.sum(jets[:, 2].take([j2, j1]))**2) ) ## sqrt(E2-p2) of combined j1+j2 4-vector
                    #nschi = 400. # test value
                    if use_merged:
                        #set_trace()
                        maxmjet = max( np.sqrt(jets[j1, 3]**2 - np.sum(jets[j1, 0:3]**2)), np.sqrt(jets[j2, 3]**2 - np.sum(jets[j2, 0:3]**2)) ) ## get max mass between j1 and j2
                        #nudiscr, massdiscr, prob = solver.solve_3J_merged(maxmjet, mbpjet, nschi)
                        #print('max(mjet) = ', maxmjet, ', m(b+j) = ', mbpjet, ', nschi = ', nschi, ', Prob = ', prob, ', MassDiscr = ', massdiscr, ', NuDiscr = ', nudiscr)
                        print('max(mjet) = ', maxmjet, ', m(b+j) = ', mbpjet)
                    else:
                        #nudiscr, massdiscr, prob = solver.solve_3J_lost(mbpjet, nschi)
                        #print('m(b+j) = ', mbpjet, ', nschi = ', nschi, ', Prob = ', prob, ', MassDiscr = ', massdiscr, ', NuDiscr = ', nudiscr)
                        print('m(b+j) = ', mbpjet)
                    print(j0, j1, j2)

        start += njets
        evt_idx += 1
        print()
        #set_trace()

#@njit()
def make_permutations(jets, leptons, MET):
    nj = jets.counts[0:5]
        ## make jets matrix
    jpx = jets.p4.x.flatten()[0:nj[0]*5]
    jpy = jets.p4.y.flatten()[0:nj[0]*5]
    jpz = jets.p4.z.flatten()[0:nj[0]*5]
    jE = jets.p4.energy.flatten()[0:nj[0]*5]
    jBtag = jets.BTAG_DEEPCSVMEDIUM.flatten()[0:nj[0]*5]
    #set_trace()
    jets_inputs = np.stack((jpx, jpy, jpz, jE, jBtag), axis=1) # one row has (px, py, pyz, E, btag Pass)
    #jets_inputs = np.stack((jets.p4.x.flatten(), jets.p4.y.flatten(), jets.p4.z.flatten(), jets.p4.energy.flatten()), axis=1) # one row has (px, py, pyz, E)

        ## make leptons matrix
    leppx = leptons.p4.x.flatten()[0:5]
    leppy = leptons.p4.y.flatten()[0:5]
    leppz = leptons.p4.z.flatten()[0:5]
    lepE = leptons.p4.energy.flatten()[0:5]
    #lepm = leptons.p4.mass.flatten()[0:5]
    leptons_inputs = np.stack((leppx, leppy, leppz, lepE), axis=1) # one row has (px, py, pyz, E)
    #leptons_inputs = np.stack((leptons.p4.x.flatten(), leptons.p4.y.flatten(), leptons.p4.z.flatten(), leptons.p4.energy.flatten()), axis=1) # one row has (px, py, pyz, E)

        ## make MET matrix
    met_px = MET.x[0:5]
    met_py = MET.y[0:5]
    met_inputs = np.stack((met_px, met_py), axis=1) # one row has (px, py)
    #met_inputs = np.stack((MET.x, MET.y), axis=1) # one row has (px, py)

    set_trace()
    get_permutations_4j(njets_array=nj, jets=jets_inputs, leptons=leptons_inputs, met=met_inputs)
    get_permutations_3j(njets_array=nj, jets=jets_inputs, leptons=leptons_inputs, met=met_inputs)

