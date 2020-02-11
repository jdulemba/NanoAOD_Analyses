from numba import njit
import numpy as np
from pdb import set_trace
import python.TTBarSolver as ttsolver

solver = ttsolver.TTBarSolver()
'''
The following numpy arrays are taken directly from the first five events after selecting tight leptons and 4 jets.
Only for debugging purposes and testing at this point.
'''

test_4jets_array = np.array([4, 4, 4, 4, 4])
    ## (px, py, pz, E)
test_4jets = np.array([
    [ 1.46588218e+00, -9.68639069e+01, -4.20600510e+01, 1.05743782e+02], [ 8.11790390e+01,  3.44528503e+01,  4.08130493e+01, 9.74636612e+01],
    [-3.15172863e+01,  4.83540573e+01, -6.89234390e+01, 8.99276505e+01], [-7.96223402e+00, -2.11776638e+01, -5.77623100e+01, 6.22156715e+01],

    [ 3.90159683e+01,  1.20661369e+02, -3.57599678e+01, 1.32530609e+02], [-1.08127327e+01, -8.90333252e+01,  2.68187771e+01, 9.39489670e+01],
    [-3.50565262e+01, -6.98180161e+01,  1.05574127e+02, 1.31688095e+02], [-3.01058578e+01,  5.52911520e+00, -6.20904198e+01, 6.95760498e+01],

    [-8.41904984e+01,  4.62597466e+01, -3.28814163e+01, 1.02658005e+02], [ 3.70532799e+01, -5.96098938e+01, -8.36687946e+00, 7.08345718e+01],
    [ 4.14425392e+01,  3.66331635e+01, -1.06966200e+01, 5.65677986e+01], [ 1.84616814e+01,  1.45654039e+01,  3.60267448e+01, 4.33081932e+01],

    [ 5.94906197e+01,  2.19912910e+00, -4.02021790e+01, 7.27764816e+01], [-1.38686590e+01, -3.26109886e+01, -1.00762281e+01, 3.69552574e+01],
    [ 1.11924324e+01, -3.21049767e+01, -3.61392059e+01, 4.96597328e+01], [-1.13249636e+01,  2.74256325e+01,  3.71226239e+00, 3.06056423e+01],

    [-4.54414253e+01, -1.78687485e+02,  4.08158073e+01, 1.89382706e+02], [ 3.12490836e-02,  1.30375000e+02,  6.01012695e+02, 6.15328918e+02],
    [ 1.57233524e+01,  1.07481003e+02,  2.08606018e+02, 2.35406708e+02], [ 3.91699829e+01, -6.18972969e+01,  7.21763992e+01, 1.03180359e+02]
])
    ## (px, py, pz, E)
test_4jets_leptons = np.array([
    [-28.915369 ,  44.079334 , -62.7505   ,  81.95561  ],
    [ -9.069429 , -82.03243  ,  24.337626 ,  86.04596  ],
    [ 34.898754 , -55.507015 ,  -7.9238267,  66.04355  ],
    [ 10.259361 , -29.684729 , -33.34455  ,  45.80729  ],
    [ 12.323418 ,  88.61602  , 172.76591  , 194.5578   ]
])
    ## (px, py)
test_4jets_met = np.array([
    [-54.914692 ,  49.528954 ],
    [ 39.98818  ,   6.676497 ],
    [ -6.4915404, -43.05015  ],
    [-77.96277  ,  33.763004 ],
    [-10.805938 ,  22.055433 ]
])

'''
The following numpy arrays are taken directly from the first five events after selecting tight leptons and 3 jets.
Only for debugging purposes and testing at this point.
'''

test_3jets_array = np.array([3, 3, 3, 3, 3])
    ## (px, py, pz, E)
test_3jets = np.array([
    [  10.838933 ,   70.48397  , -233.93109  ,  244.79097  ], [ -50.59421  ,   42.6733   ,   56.524982 ,   87.87881  ], [  -8.887823 ,   27.07871  ,  -10.066564 ,   30.332691 ],
    [ 134.49626  ,    1.002671 , -134.5415   ,  190.8868   ], [ -75.891464 ,  -25.504702 ,  139.24947  ,  161.12848  ], [  41.0424   ,   18.980433 ,   -1.345666 ,   45.53924  ],
    [ -22.50881  ,  -69.37746  ,  -20.23488  ,   76.86495  ], [  54.648956 ,   19.614931 ,   27.129189 ,   64.31939  ], [   2.4663734,  -20.477003 ,  -11.23544  ,   23.759773 ],
    [  23.022575 ,  -68.61561  ,  124.99501  ,  144.5074   ], [  21.343796 ,   39.04731  ,  187.74185  ,  193.03026  ], [ -33.189785 ,   25.327862 ,   72.32822  ,   83.95407  ],
    [  60.489555 ,   18.26674  ,  195.12999  ,  205.40303  ], [  -4.60639  ,   30.497839 ,    5.30353  ,   32.046825 ], [  10.017437 ,  -28.443811 ,   -9.153948 ,   31.59562  ]
])
    ## (px, py, pz, E)
test_3jets_leptons = np.array([
    [ -8.560334 ,  25.910042 ,  -9.73551  ,  28.972418 ],
    [ 34.05833  ,  14.846232 ,  -0.8565457,  37.163494 ],
    [ 47.136646 ,  17.292995 ,  23.182486 ,  55.302353 ],
    [ 20.004333 , -59.668453 , 108.167015 , 125.14237  ],
    [  9.489876 , -27.157267 ,  -8.444733 ,  29.981653 ]
])
    ## (px, py)
test_3jets_met = np.array([
    [  10.396741 , -126.84293  ],
    [ -92.50585  ,   -3.2074502],
    [   3.3195245,   44.79466  ],
    [  -8.249034 ,  -70.06654  ],
    [ -67.42234  ,  -28.88607  ]
])


#@njit()
def get_permutations_4j(njets_array, jets, leptons, met):
    start = 0
    stop = 0
    evt_idx = 0

    #set_trace()
    for njets in njets_array:
        stop += njets
        for j0 in range(start, stop):
            ## run NS to get nschi2 to be used solver
            #run_nu_solver(leptons[evt_idx], jets[j0], met[evt_idx]) ## passes all info for only one lepton (px, py, pz, E), one jet (px, py, pz, E), and met (px, py) from the event at a time

            ## lepton and met info won't change inside njets for loop
            ## jet info will change for each value of j0 (4 separate for event with 4 jets)
            for j1 in range(start, stop):
                if j1 == j0: continue
                for j2 in range(start, stop):
                    if j2 == j0 or j2 == j1: continue
                    for j3 in range(start, stop):
                        if j3 == j0 or j3 == j1 or j3 == j2: continue
                        #set_trace()
                        mthad = np.sqrt( (jets[j1] + jets[j2] + jets[j3])[3]**2 -( (jets[j1]+jets[j2]+jets[j3])[0]**2 + (jets[j1]+jets[j2]+jets[j3])[1]**2 + (jets[j1]+jets[j2]+jets[j3])[2]**2 ) ) ## sqrt(E2-p2) of combined j1+j2+j3 4-vector
                        mwhad = np.sqrt( (jets[j2] + jets[j3])[3]**2 -( (jets[j2]+jets[j3])[0]**2 + (jets[j2]+jets[j3])[1]**2 + (jets[j2]+jets[j3])[2]**2 ) ) ## sqrt(E2-p2) of combined j2+j3 4-vector
                        nschi = 400. # test value
                        solver.solve_4PJ_(mthad, mwhad, nschi)
                        print('mthad = ', mthad, ', mwhad = ', mwhad, ', Prob = ', solver.Prob_, ', MassDiscr = ', solver.MassDiscr_, ', NuDiscr = ', solver.NuDiscr_)

        start += njets
        evt_idx += 1
        print()
        #set_trace()

#@njit()
def get_permutations_3j(njets_array, jets, leptons, met, use_merged=False): # jets(px, py, pz, E), leptons(px, py, pz, E), met(px, py)
    start = 0
    stop = 0
    evt_idx = 0

    #set_trace()
    for njets in njets_array:
        stop += njets
        for j0 in range(start, stop):
            ## run NS to get nschi2 to be used solver
            #run_nu_solver(leptons[evt_idx], jets[j0], met[evt_idx]) ## passes all info for only one lepton (px, py, pz, E), one jet (px, py, pz, E), and met (px, py) from the event at a time

            ## lepton and met info won't change inside njets for loop
            ## jet info will change for each value of j0 (3 separate for event with 3 jets)
            for j1 in range(start, stop):
                if j1 == j0: continue
                for j2 in range(start, stop):
                    if j2 == j0 or j2 == j1: continue
                    #set_trace()
                    mbpjet = np.sqrt((jets[j1]+jets[j2])[3]**2-((jets[j1]+jets[j2])[0]**2+(jets[j1]+jets[j2])[1]**2+(jets[j1]+jets[j2])[2]**2)) ## sqrt(E2-p2) of combined j1+j2 4-vector
                    nschi = 400. # test value
                    if use_merged:
                        maxmjet = max(np.sqrt(jets[j1][3]**2-(jets[j1][0]**2+jets[j1][1]**2+jets[j1][2]**2)), np.sqrt(jets[j2][3]**2-(jets[j2][0]**2+jets[j2][1]**2+jets[j2][2]**2)))
                        solver.solve_3J_merged_(maxmjet, mbpjet, nschi)
                        print('max(mjet) = ', maxmjet, ', m(b+j) = ', mbpjet, ', Prob = ', solver.Prob_, ', MassDiscr = ', solver.MassDiscr_, ', NuDiscr = ', solver.NuDiscr_)
                    else:
                        solver.solve_3J_lost_(mbpjet, nschi)
                        print('m(b+j) = ', mbpjet, ', Prob = ', solver.Prob_, ', MassDiscr = ', solver.MassDiscr_, ', NuDiscr = ', solver.NuDiscr_)

        start += njets
        evt_idx += 1
        print()
        #set_trace()

get_permutations_3j(test_3jets_array, test_3jets, test_3jets_leptons, test_3jets_met)
#get_permutations_3j(test_3jets_array, test_3jets, test_3jets_leptons, test_3jets_met, use_merged=True)
#get_permutations_4j(test_4jets_array, test_4jets, test_4jets_leptons, test_4jets_met)
#set_trace()

#@njit()
def make_permutations(jets, leptons, MET):
    nj = jets.counts[0:5]
        ## make jets matrix
    jpx = jets.p4.x.flatten()[0:nj[0]*5]
    jpy = jets.p4.y.flatten()[0:nj[0]*5]
    jpz = jets.p4.z.flatten()[0:nj[0]*5]
    jE = jets.p4.energy.flatten()[0:nj[0]*5]
    #set_trace()
    jets_inputs = np.stack((jpx, jpy, jpz, jE), axis=1) # one row has (px, py, pyz, E)
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

