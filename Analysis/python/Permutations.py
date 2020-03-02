from numba import njit, objmode
from numba.typed import List
import numpy as np
from pdb import set_trace
import python.TTBarSolver as solver
import compiled.pynusolver as pynusolver
from coffea.arrays import Initialize
import awkward

@njit()
def get_permutations(njets_array, jets, leptons, met, use_merged=False):
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
    print('Finding best perms')
    start = 0
    stop = 0
    evt_idx = 0

    #best_perms = List()
    best_perms_ordering = List()
    best_perms_nus = List()
    best_perms_probs = List()

    for njets in njets_array:
        #print('evt idx: ', evt_idx, ', njets: ', njets)
        stop += njets

        ## initialize best_perm lists for event
        best_perm_ordering = np.array([-10, -10, -10, -10])
        best_perm_nu = np.zeros(4)
        best_perm_probs = np.array([np.inf, np.inf, np.inf])

        for j0 in range(start, stop):
            ## require btagging
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
                    for j2 in range(j1+1, stop):
                        if j2 == j0: continue
                            ## advanced indexing doesn't work with numba at this point, can only add px, py, pz separately
                        mbpjet = np.sqrt( np.sum(jets[:, 3].take([j2, j1]))**2 - (np.sum(jets[:, 0].take([j2, j1]))**2 + np.sum(jets[:, 1].take([j2, j1]))**2 + np.sum(jets[:, 2].take([j2, j1]))**2) ) ## sqrt(E2-p2) of combined j1+j2 4-vector
                        if use_merged:
                            maxmjet = max( np.sqrt(jets[j1, 3]**2 - np.sum(jets[j1, 0:3]**2)), np.sqrt(jets[j2, 3]**2 - np.sum(jets[j2, 0:3]**2)) ) ## get max mass between j1 and j2
                            nudiscr, massdiscr, prob = solver.solve_3J_merged(maxmjet, mbpjet, nu[3]) ## input order is (maxmjet, mbpjet, nschi2 value)
                            #print('max(mjet) = ', maxmjet, ', m(b+j) = ', mbpjet, ', nschi = ', nu[3], ', Prob = ', prob, ', MassDiscr = ', massdiscr, ', NuDiscr = ', nudiscr)
                        else:
                            nudiscr, massdiscr, prob = solver.solve_3J_lost(mbpjet, nu[3]) ## input order is (mbpjet, nschi2 value)
                            #print('m(b+j) = ', mbpjet, ', nschi = ', nu[3], ', Prob = ', prob, ', MassDiscr = ', massdiscr, ', NuDiscr = ', nudiscr)
                        if prob < best_perm_probs[0]: ## have to be careful with indexing!!
                            best_perm_ordering = np.array([j0-start, j1-start, j2-start, -10]) # set last element to -10 in order to make array have the same shape as 4+ jets
                            best_perm_nu = nu
                            best_perm_probs = np.array([prob, massdiscr, nudiscr])


                ## get best perms for 4+ jet events
            else:
                #print('4+ jets')
                for j1 in range(start, stop):
                    if j1 == j0: continue
                    for j2 in range(start, stop):
                        if j2 == j0 or j2 == j1: continue
                        for j3 in range(j2+1, stop):
                            if j3 == j0 or j3 == j1: continue
                                ## advanced indexing doesn't work with numba at this point, can only add px, py, pz separately
                            mthad = np.sqrt( np.sum(jets[:, 3].take([j1, j2, j3]))**2 - (np.sum(jets[:, 0].take([j1, j2, j3]))**2 + np.sum(jets[:, 1].take([j1, j2, j3]))**2 + np.sum(jets[:, 2].take([j1, j2, j3]))**2) ) ## sqrt(E2-p2) of combined j1+j2+j3 4-vector
                            mwhad = np.sqrt( np.sum(jets[:, 3].take([j2, j3]))**2 - (np.sum(jets[:, 0].take([j2, j3]))**2 + np.sum(jets[:, 1].take([j2, j3]))**2 + np.sum(jets[:, 2].take([j2, j3]))**2) ) ## sqrt(E2-p2) of combined j2+j3 4-vector
                            nudiscr, massdiscr, prob = solver.solve_4PJ(mthad, mwhad, nu[3]) ## input order is (mthad, mwhad, nschi2 value)
                            #print('mthad = ', mthad, ', mwhad = ', mwhad, ', nschi = ', nu[3], ', Prob = ', prob, ', MassDiscr = ', massdiscr, ', NuDiscr = ', nudiscr)
                            if prob < best_perm_probs[0]: ## have to be careful with indexing!!
                                best_perm_ordering = np.array([j0-start, j1-start, j2-start, j3-start])
                                best_perm_nu = nu
                                best_perm_probs = np.array([prob, massdiscr, nudiscr])

        best_perms_ordering.append(best_perm_ordering)
        best_perms_nus.append(best_perm_nu)
        best_perms_probs.append(best_perm_probs)

        start += njets
        if np.mod(evt_idx, 10) == 0:
            print('\t', evt_idx+1,'/',len(njets_array), 'events processed')
        evt_idx += 1
    print('Finished finding best perms')
    return best_perms_ordering, best_perms_nus, best_perms_probs


def make_perm_table(jets, leptons, nus, jets_orderings, probs, evt_weights):
    '''
    Inputs:
        Jets, leptons, neutrinos, jet assignment object ordering, array of disrminiant probabilities (Total, mass discriminant, neutrino discrminant), and associated event weights
    Returns:
        awkward Table containing
            1: Jet objects (BLeps, BHads, WJas, WJbs)
            2: Leptons, Neutrinos
            3: Calculated probabilities (Total -> Prob, mass discriminant -> MassDiscr, neutrino discrminant -> NuDiscr)
            4: Number of jets in events (njets)
            5: Event weights (evt_wts)
    '''
    neutrinos = Initialize({
        'px' : nus[:, 0:1].flatten(),
        'py' : nus[:, 1:2].flatten(),
        'pz' : nus[:, 2:3].flatten(),
        'mass' : 0,
        'chi2' : nus[:, 3:4].flatten()
    })

        ## columns of jet_orderings correspond to blep, bhad, wja wjb
    ## blep objects jets[np.arange(len(jets_orderings)), jets_orderings[:,0]]
    BLeps = Initialize({
        'px' : jets[np.arange(len(jets_orderings)), jets_orderings[:,0]].p4.x.flatten(),
        'py' : jets[np.arange(len(jets_orderings)), jets_orderings[:,0]].p4.y.flatten(),
        'pz' : jets[np.arange(len(jets_orderings)), jets_orderings[:,0]].p4.z.flatten(),
        'mass' : jets[np.arange(len(jets_orderings)), jets_orderings[:,0]].p4.mass.flatten(),
    })
    BHads = Initialize({
        'px' : jets[np.arange(len(jets_orderings)), jets_orderings[:,1]].p4.x.flatten(),
        'py' : jets[np.arange(len(jets_orderings)), jets_orderings[:,1]].p4.y.flatten(),
        'pz' : jets[np.arange(len(jets_orderings)), jets_orderings[:,1]].p4.z.flatten(),
        'mass' : jets[np.arange(len(jets_orderings)), jets_orderings[:,1]].p4.mass.flatten(),
    })
    WJas = Initialize({
        'px' : jets[np.arange(len(jets_orderings)), jets_orderings[:,2]].p4.x.flatten(),
        'py' : jets[np.arange(len(jets_orderings)), jets_orderings[:,2]].p4.y.flatten(),
        'pz' : jets[np.arange(len(jets_orderings)), jets_orderings[:,2]].p4.z.flatten(),
        'mass' : jets[np.arange(len(jets_orderings)), jets_orderings[:,2]].p4.mass.flatten(),
    })
        ## initialize WJb to zeros
    WJbs = Initialize({
        'px' : np.zeros(len(jets_orderings)),
        'py' : np.zeros(len(jets_orderings)),
        'pz' : np.zeros(len(jets_orderings)),
        'mass' : np.zeros(len(jets_orderings)),
    })
        ## fill values only for events with 4+ jets
    WJbs.x[jets_orderings[:,3] >= 0] = jets[np.arange(len(jets_orderings))[jets_orderings[:,3] >= 0], jets_orderings[:,3][jets_orderings[:,3] >= 0]].p4.x.flatten()
    WJbs.y[jets_orderings[:,3] >= 0] = jets[np.arange(len(jets_orderings))[jets_orderings[:,3] >= 0], jets_orderings[:,3][jets_orderings[:,3] >= 0]].p4.y.flatten()
    WJbs.z[jets_orderings[:,3] >= 0] = jets[np.arange(len(jets_orderings))[jets_orderings[:,3] >= 0], jets_orderings[:,3][jets_orderings[:,3] >= 0]].p4.z.flatten()
    WJbs.mass[jets_orderings[:,3] >= 0] = jets[np.arange(len(jets_orderings))[jets_orderings[:,3] >= 0], jets_orderings[:,3][jets_orderings[:,3] >= 0]].p4.mass.flatten()
    WJbs.E[jets_orderings[:,3] >= 0] = jets[np.arange(len(jets_orderings))[jets_orderings[:,3] >= 0], jets_orderings[:,3][jets_orderings[:,3] >= 0]].p4.E.flatten()

        ## Combine everything into a single table
    perm_table = awkward.Table(
        BLeps=BLeps,
        BHads=BHads,
        WJas=WJas,
        WJbs=WJbs,
        Leptons=leptons,
        Nus=neutrinos,
        Prob=probs[:,0],
        MassDiscr=probs[:,1],
        NuDiscr=probs[:,2],
        njets=jets.counts,
        evt_wts=evt_weights,
    )
    print('Table of permutation objects created')
    return perm_table


#@njit()
def find_best_permutations(jets, leptons, MET, evt_weights):
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
    print('Making permutations')
        ## for testing
    ##set_trace()
    #nj = jets.counts[0:25]
    #    ## make jets matrix
    #jpx = jets.p4.x.flatten()[0:sum(nj)]
    #jpy = jets.p4.y.flatten()[0:sum(nj)]
    #jpz = jets.p4.z.flatten()[0:sum(nj)]
    #jE = jets.p4.energy.flatten()[0:sum(nj)]
    #jBtag = jets.BTAG_DEEPCSVMEDIUM.flatten()[0:sum(nj)]
    ##set_trace()
    #jets_inputs = np.stack((jpx, jpy, jpz, jE, jBtag), axis=1) # one row has (px, py, pyz, E, btag Pass)

    #    ## make leptons matrix
    #leppx = leptons.p4.x.flatten()[0:len(nj)]
    #leppy = leptons.p4.y.flatten()[0:len(nj)]
    #leppz = leptons.p4.z.flatten()[0:len(nj)]
    #lepE = leptons.p4.energy.flatten()[0:len(nj)]
    ##lepm = leptons.p4.mass.flatten()[0:len(nj)]
    #leptons_inputs = np.stack((leppx, leppy, leppz, lepE), axis=1) # one row has (px, py, pyz, E)

    #    ## make MET matrix
    #met_px = MET.x[0:len(nj)]
    #met_py = MET.y[0:len(nj)]
    #met_inputs = np.stack((met_px, met_py), axis=1) # one row has (px, py)

    ##bp = get_permutations(njets_array=nj, jets=jets_inputs, leptons=leptons_inputs, met=met_inputs)
    #bp_ordering, bp_nus, bp_probs = get_permutations(njets_array=nj, jets=jets_inputs, leptons=leptons_inputs, met=met_inputs)
    ##set_trace()
        ##

    jets_inputs = np.stack((jets.p4.x.flatten(), jets.p4.y.flatten(), jets.p4.z.flatten(), jets.p4.energy.flatten()), axis=1) # one row has (px, py, pyz, E)
    leptons_inputs = np.stack((leptons.p4.x.flatten(), leptons.p4.y.flatten(), leptons.p4.z.flatten(), leptons.p4.energy.flatten()), axis=1) # one row has (px, py, pyz, E)
    met_inputs = np.stack((MET.x, MET.y), axis=1) # one row has (px, py)
    bp_ordering, bp_nus, bp_probs = get_permutations(njets_array=jets.counts, jets=jets_inputs, leptons=leptons_inputs, met=met_inputs)

    ##set_trace()
    bp_ordering = np.asarray(bp_ordering)
    bp_nus = np.asarray(bp_nus)
    bp_probs = np.asarray(bp_probs)

        ## only keep permutations with some sort of solution (prob != infinity)
    #set_trace()
    mask = (bp_probs[:,0] != np.inf)
    valid_jets = jets[np.arange(len(bp_nus))][(mask)] 
    valid_leptons = leptons[np.arange(len(bp_nus))][(mask)]
    valid_nus = bp_nus[(mask)] 
    valid_ordering = bp_ordering[(mask)] 
    valid_probs = bp_probs[(mask)] 
    bp_table = make_perm_table(valid_jets, valid_leptons, valid_nus, valid_ordering, valid_probs, evt_weights[(mask)])
    #set_trace()

    return bp_table

