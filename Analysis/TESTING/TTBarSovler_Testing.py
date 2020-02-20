from pdb import set_trace
import numpy as np
import python.TTBarSolver as solver
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('solver', choices=['4+', 'merged', 'lost'], help='Choose which solver to test (merged and lost are for 3 jets)')
args = parser.parse_args()

## test inputs (mthad, mwhad, nschi) for same outputs (prob, nsdisc, massdisc) from c++ FW (converted numeric_limits == 1.79769e+308 to np.inf)
## [mthad, mwhad, nschi, MassDiscr, NuDiscr, Prob]

    ## taken from rake 'test[bin/ttbar_post_alpha_reco.cc, ttJets$, cfg_files/htt_baseline_j20_l50_MT40.cfg, l 10]' in 9410 from Solve_4PJ_Analytical
test_4jets = np.array([
    [127.718, 88.807, np.inf, 10.4354, np.inf, np.inf],
    [127.718, 38.4145, np.inf, 11.1134, np.inf, np.inf],
    [127.718, 84.0735, np.inf, 10.5469, np.inf, np.inf],
    [182.667, 88.807, 201400, 7.43935, 7.31189, 14.7512],
    [182.667, 139.808, 201400, 12.5266, 7.31189, 19.8385],
    [182.667, 79.3406, 201400, 7.78639, 7.31189, 15.0983],
    [176.222, 38.4145, np.inf, 12.3892, np.inf, np.inf],
    [176.222, 139.808, np.inf, np.inf, np.inf, np.inf],
    [176.222, 101.846, np.inf, 8.88217, np.inf, np.inf],
    [153.013, 84.0735, 19.4144, 7.91432, 4.64635, 12.5607],
    [153.013, 79.3406, 19.4144, 7.52266, 4.64635, 12.169],
    [153.013, 101.846, 19.4144, 10.4274, 4.64635, 15.0738],
    
    [130.208, 90.4691, np.inf, 10.8719, np.inf, np.inf],
    [130.208, 39.3227, np.inf, 11.1361, np.inf, np.inf],
    [130.208, 85.7135, np.inf, 10.0688, np.inf, np.inf],
    [185.642, 90.4691, 185216, 7.70249, 7.3307, 15.0332],
    [185.642, 142.115, 185216, 13.5392, 7.3307, 20.8699],
    [185.642, 80.3226, 185216, 7.98642, 7.3307, 15.3171],
    [179.238, 39.3227, np.inf, 12.7611, np.inf, np.inf],
    [179.238, 142.115, np.inf, np.inf, np.inf, np.inf],
    [179.238, 103.599, np.inf, 9.16131, np.inf, np.inf],
    [155.566, 85.7135, 40.1674, 7.91698, 4.85834, 12.7753],
    [155.566, 80.3226, 40.1674, 7.47648, 4.85834, 12.3348],
    [155.566, 103.599, 40.1674, 10.8116, 4.85834, 15.6699],
    
    [125.226, 87.1434, np.inf, 10.691, np.inf, np.inf],
    [125.226, 37.5063, np.inf, 11.4569, np.inf, np.inf],
    [125.226, 82.4315, np.inf, 10.1536, np.inf, np.inf],
    [179.687, 87.1434, 219189, 7.35923, 7.36811, 14.7273],
    [179.687, 137.496, 219189, 13.5635, 7.36811, 20.9316],
    [179.687, 78.3583, 219189, 7.70542, 7.36811, 15.0735],
    [173.2, 37.5063, np.inf, 12.7873, np.inf, np.inf],
    [173.2, 137.496, np.inf, 12.4446, np.inf, np.inf],
    [173.2, 100.088, np.inf, 8.97838, np.inf, np.inf],
    [150.457, 82.4315, 5.46926, 7.99238, 4.3101, 12.3025],
    [150.457, 78.3583, 5.46926, 7.758, 4.3101, 12.0681],
    [150.457, 100.088, 5.46926, 10.4241, 4.3101, 14.7342]
])

    ## taken from 'test[bin/ttbar_post_alpha_reco.cc, ttJets$, cfg_files/htt_baseline_j20_l50_MT40.cfg, l 10]' in 9410 from Solve_3J_Lost_Analytical
# [Mbpjet, nschi, MassDiscr, NuDiscr, Prob]
test_3j_lost = np.array([
    [150.079, 422.248, 3.12498, 5.42828, 8.55326],
    [150.079, 422.248, 3.12498, 5.42828, 8.55326],
    [65.2468, 81.9903, 3.3043, 4.9988, 8.3031],
    [65.2468, 81.9903, 3.3043, 4.9988, 8.3031],
    [151.925, 424.878, 3.12498, 5.42828, 8.55326],
    [151.925, 424.878, 3.12498, 5.42828, 8.55326],
    
    [151.997, 479.689, 3.12498, 5.55867, 8.68365],
    [151.997, 479.689, 3.12498, 5.55867, 8.68365],
    [66.1589, 36.974, 3.3043, 4.77979, 8.08409],
    [66.1589, 36.974, 3.3043, 4.77979, 8.08409],
    [153.786, 479.092, 3.12498, 5.55867, 8.68365],
    [153.786, 479.092, 3.12498, 5.55867, 8.68365],
    
    [148.161, 367.954, 2.54095, 5.34663, 7.88758],
    [148.161, 367.954, 2.54095, 5.34663, 7.88758],
    [64.3347, 146.42, 3.3043, 5.21927, 8.52357],
    [64.3347, 146.42, 3.3043, 5.21927, 8.52357],
    [150.064, 373.29, 3.12498, 5.34663, 8.47161],
    [150.064, 373.29, 3.12498, 5.34663, 8.47161]
])

    ## taken from 'test[bin/ttbar_post_alpha_reco.cc, ttJets$, cfg_files/htt_baseline_j20_l50_MT40.cfg, l 10]' in 9410 from Solve_3J_Merged_Analytical
# [Mbpjet, MaxMjet, nschi, MassDiscr, NuDiscr, Prob]
test_3j_merged = np.array([
    [150.079, 14.307, 422.248, 7.25512, 5.60423, 12.8594],
    [150.079, 14.307, 422.248, 7.25512, 5.60423, 12.8594],
    [65.2468, 13.2822, 81.9903, 12.725, 4.79218, 17.5172],
    [65.2468, 13.2822, 81.9903, 12.725, 4.79218, 17.5172],
    [151.925, 14.307, 424.878, 7.25512, 5.60423, 12.8594],
    [151.925, 14.307, 424.878, 7.25512, 5.60423, 12.8594],
    
    [151.997, 14.4651, 479.689, 7.25512, 5.68832, 12.9434],
    [151.997, 14.4651, 479.689, 7.25512, 5.68832, 12.9434],
    [66.1589, 13.4608, 36.974, 12.725, 4.96922, 17.6943],
    [66.1589, 13.4608, 36.974, 12.725, 4.96922, 17.6943],
    [153.786, 14.4651, 479.092, 7.25512, 5.68832, 12.9434],
    [153.786, 14.4651, 479.092, 7.25512, 5.68832, 12.9434],
    
    [148.161, 14.1489, 367.954, 7.30527, 5.56269, 12.868],
    [148.161, 14.1489, 367.954, 7.30527, 5.56269, 12.868],
    [64.3347, 13.1035, 146.42, 12.725, 5.17114, 17.8962],
    [64.3347, 13.1035, 146.42, 12.725, 5.17114, 17.8962],
    [150.064, 14.1489, 373.29, 7.25512, 5.56269, 12.8178],
    [150.064, 14.1489, 373.29, 7.25512, 5.56269, 12.8178]
])

def same_output(computed=[], expected=[]):
    mass_discrs = [computed[0], expected[0]]
    nu_discrs = [computed[1], expected[1]]
    probs = [computed[2], expected[2]]

    #set_trace()
    if not (len(np.unique(np.around(mass_discrs, decimals=4))) == 1):
        print('MassDiscr values not the same: ', mass_discrs)
    if not (len(np.unique(np.around(nu_discrs, decimals=4))) == 1):
        print('NSDiscr values not the same: ', nu_discrs)
    if not (len(np.unique(np.around(probs, decimals=4))) == 1):
        print('Prob values not the same: ', probs)


if args.solver == '4+':
    for perm_idx in range(len(test_4jets[:])):
        mthad = test_4jets[perm_idx][0]
        mwhad = test_4jets[perm_idx][1]
        nschi = test_4jets[perm_idx][2]
        massD = test_4jets[perm_idx][3]
        nuD = test_4jets[perm_idx][4]
        prob = test_4jets[perm_idx][5]
    
        nudiscr, massdiscr, prob = solver.solve_4PJ(mthad, mwhad, nschi)
        #set_trace()
        #print('NSchi: ', nschi, ', NuDiscr: ' , nudiscr)
        #print('mthad: ', mthad, ', mwhad: ' , mwhad, ', MassDiscr: ', massdiscr)
        print('Computed MassDiscr = ', massdiscr, ', NuDiscr = ', nudiscr, ', Prob = ', prob)
        #same_output(computed=[massdiscr, nudiscr, prob], expected=[massD, nuD, prob])

if args.solver == 'lost':
    for perm_idx in range(len(test_3j_lost[:])):
        mbpjet = test_3j_lost[perm_idx][0]
        nschi = test_3j_lost[perm_idx][1]
        massD = test_3j_lost[perm_idx][2]
        nuD = test_3j_lost[perm_idx][3]
        prob = test_3j_lost[perm_idx][4]
    
        nudiscr, massdiscr, prob = solver.solve_3J_lost(mbpjet, nschi)
        #print('NSchi: ', nschi, ', NuDiscr: ' , nudiscr)
        print('m(b+jet): ', mbpjet, ', MassDiscr: ', massdiscr)
        #print('Computed MassDiscr = ', massdiscr, ', NuDiscr = ', nudiscr, ', Prob = ', prob)
    
        #same_output(computed=[massdiscr, nudiscr, prob], expected=[massD, nuD, prob])

if args.solver == 'merged':
    for perm_idx in range(len(test_3j_merged[:])):
        mbpjet = test_3j_merged[perm_idx][0]
        maxmjet = test_3j_merged[perm_idx][1]
        nschi = test_3j_merged[perm_idx][2]
        massD = test_3j_merged[perm_idx][3]
        nuD = test_3j_merged[perm_idx][4]
        prob = test_3j_merged[perm_idx][5]
    
        nudiscr, massdiscr, prob = solver.solve_3J_merged(maxmjet, mbpjet, nschi)
        #print('NSchi: ', nschi, ', NuDiscr: ' , nudiscr)
        print('max m(jet): ', maxmjet, ', m(b+jet): ', mbpjet, ', MassDiscr: ', massdiscr)
        #print('Computed MassDiscr = ', massdiscr, ', NuDiscr = ', nudiscr, ', Prob = ', prob)
    
        #same_output(computed=[massdiscr, nudiscr, prob], expected=[massD, nuD, prob])
    
