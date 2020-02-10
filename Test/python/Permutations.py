from numba import njit
import numpy as np
from pdb import set_trace

'''
The following numpy arrays are taken directly from the first five events after selecting tight leptons and 4 jets.
Only for debugging purposes and testing at this point.
'''

test_njets_array = np.array([4, 4, 4, 4, 4])
    ## (px, py, pz, E)
test_jets = np.array([
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
test_leptons = np.array([
    [-28.915369 ,  44.079334 , -62.7505   ,  81.95561  ],
    [ -9.069429 , -82.03243  ,  24.337626 ,  86.04596  ],
    [ 34.898754 , -55.507015 ,  -7.9238267,  66.04355  ],
    [ 10.259361 , -29.684729 , -33.34455  ,  45.80729  ],
    [ 12.323418 ,  88.61602  , 172.76591  , 194.5578   ]
])
    ## (px, py)
test_met = np.array([
    [-54.914692 ,  49.528954 ],
    [ 39.98818  ,   6.676497 ],
    [ -6.4915404, -43.05015  ],
    [-77.96277  ,  33.763004 ],
    [-10.805938 ,  22.055433 ]
])

@njit()
def get_permutations_4j(njets_array, jets, leptons, met):
    start = 0
    stop = 0
    evt_idx = 0

    #set_trace()
    for njets in njets_array:
        stop += njets
        print('\tEvt index: ', evt_idx)
        for j0 in range(start, stop):
            print('j0: ', jets[j0], ', lep: ', leptons[evt_idx], ', met: ', met[evt_idx])
            ## run NS
            #run_nu_solver(leptons[evt_idx], jets[j0], met[evt_idx]) ## passes all info for only one lepton (px, py, pz, E), one jet (px, py, pz, E), and met (px, py) from the event at a time
            ## lepton and met info won't change inside njets for loop
            ## jet info will change for each value of j0 (4 separate for event with 4 jets)
            for j1 in range(start, stop):
                if j1 == j0: continue
                for j2 in range(start, stop):
                    if j2 == j0 or j2 == j1: continue
                    for j3 in range(start, stop):
                        if j3 == j0 or j3 == j1 or j3 == j2: continue
                        #print(j0, j1, j2, j3)
                        #print('jets px: ', jets_px[j0], jets_px[j1], jets_px[j2], jets_px[j3])
        start += njets
        evt_idx += 1
        #set_trace()

get_permutations_4j(test_njets_array, test_jets, test_leptons, test_met)
#set_trace()

##@njit()
def make_permutations(jets, leptons, MET):
    nj = jets.counts[0:5]
        ## make jets matrix
    jpx = jets.p4.x.flatten()[0:20]
    jpy = jets.p4.y.flatten()[0:20]
    jpz = jets.p4.z.flatten()[0:20]
    jE = jets.p4.energy.flatten()[0:20]
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


#    for i in range(jets.size):
#        evt_perms = list(itertools.permutations(jets[i]))
#        for idx, perm in enumerate(evt_perms):
#            set_trace()
#            if len(perm) == 3:
#                max_mjet = max(perm[0].__fast_mass, perm[1].__fast_mass)
#                mass_dijet = (perm[0].p4 + perm[1].p4).mass
#
#                print('Perm index ', idx, 'max(mjet)=', max_mjet, ', m(jj)=', mass_dijet)
#            else:
#                mThad = (perm[0].p4 + perm[1].p4 + perm[2].p4).mass
#                mWhad = (perm[0].p4 + perm[1].p4).mass
#
#                ## check matching for event categorization
#
#                print('Perm index ', idx, 'm(thad)=', mThad, ', m(Whad)=', mWhad)
#        #if len(evt_perms[0]) == 3:
#        #    perm_3jets(evt_perms)
#        #else:
#        #    perm_4pjets(evt_perms)
#
#        ## create function to evaluate best perm
#        #best_evt_perm = best_perm(evt_perms)
#
#        set_trace()
