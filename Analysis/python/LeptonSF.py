from coffea.lookup_tools.root_converters import convert_histo_root_file
from coffea.lookup_tools.dense_lookup import dense_lookup
from pdb import set_trace
import os
import numpy as np

proj_dir = os.environ['PROJECT_DIR']

    ## hardcoded right now for testing with 2016
year = "2016"
lep_info = {}
lep_info["2016"] = {
    'Muons' : {
        'filename': 'muon_2016Legacy_ID_ISO_oldTrg_SFs.root',
        'SFs' : {
            'trg' : {
                'available' : 1,
                'abs_eta' : 1
            },
            'trk' : {
                'available' : 0,
                #'abs_eta' : 
            },
            'iso' : {
                'available' : 1,
                'abs_eta' : 0
            },
            'id'  : {
                'available' : 1,
                'abs_eta' : 0
            }
        },
        'pt_as_x' : 0,
    },
    'Electrons' : {
        'filename': 'electron_sf_2016Legacy_TightCutID_Preliminary.root',
        'SFs' : {
            'trg' : {
                'available' : 1,
                'abs_eta' : 0,
            },
            'trk' : {
                'available' : 0
            },
            'iso' : {
                'available' : 0,
                'abs_eta' : 0
            },
            'id' : {
                'available' : 1,
                'abs_eta' : 0
            },
        },
        'pt_as_x' : 0,
    }
}
    ##

class LeptonSF(object):
    def __init__(self):

        self.lepSFs_ = {}
        for lepton in lep_info[year].keys():
            SFfile = convert_histo_root_file('%s/inputs/data/%s' % (proj_dir, lep_info[year][lepton]['filename']))
            for sf_type in lep_info[year][lepton]['SFs'].keys():
                if lep_info[year][lepton]['SFs'][sf_type]['available']:
                    self.lepSFs_['%s_%s' % (lepton, sf_type)] = dense_lookup(*SFfile[(sf_type, 'dense_lookup')])

    def get_2d_weights_(self, lepton, sf_type, pt_array, eta_array, shift=None): ## pt and eta are numpy arrays
        eta_array = np.abs(eta_array) if lep_info[year][lepton]['SFs'][sf_type]['abs_eta'] else eta_array
        if lep_info[year][lepton]['pt_as_x']:
            xvar = pt_array
            yvar = eta_array
        else:
            xvar = eta_array
            yvar = pt_array
        #weights = self.lepSFs_['%s_%s' % (lepton, sf_type)](xvar, yvar)
        #print('%s %s SFs: ' % (lepton, sf_type), weights)
        #return weights
        return self.lepSFs_['%s_%s' % (lepton, sf_type)](xvar, yvar)

    def get_sf_(self, lepton=None, pt_array=None, eta_array=None): ## pt and eta are separate 1D numpy arrays atm
        tot_sf = np.ones(pt_array.size)
        for sf in self.lepSFs_.keys():
            if lepton not in sf: continue
            tot_sf = np.multiply(tot_sf, self.get_2d_weights_(lepton=sf.split('_')[0], sf_type=sf.split('_')[1], pt_array=pt_array, eta_array=eta_array))

        return tot_sf


'''
Test inputs taken from rake 'test[bin/ttbar_post_alpha_reco.cc, ttJets$, cfg_files/htt_baseline_j20_l50_MT40.cfg, l 50]' in 9410

Mu [pt, eta, expected SF]  -> output is the same as expected 

rake 'test[bin/ttbar_post_alpha_reco.cc, ttJets$, cfg_files/htt_baseline_j20_l50_MT40.cfg, l 100]' in 9410
Electron [pt, etaSC, expected SF] -> output is the same as expected
'''

#lepSF = LeptonSF()
#
#test_muons = np.array([
#    [32.477, -0.382367, 0.969724],
#    [102.845, -0.509618, 0.979812],
#    [44.8985, 0.466773, 0.972154],
#    [122.024, 0.499411, 0.977203],
#    [154.409, -0.167465, 0.968942],
#])
#
#test_electrons = np.array([
#    [44.1984, 1.41288, 0.957728],
#    [41.2064, -1.24926, 0.960154],
#    [69.1812, -0.531815, 0.94979],
#    [38.232, -0.15318, 0.931888],
#    [43.6675, 0.231901, 0.975533],
#])
#
#mu_pt = test_muons[:, 0]
#mu_eta = test_muons[:, 1]
#mu_expected_sfs = test_muons[:, 2]
#
#el_pt = test_electrons[:, 0]
#el_eta = test_electrons[:, 1]
#el_expected_sfs = test_electrons[:, 2]
#
##print()
#mu_sfs = lepSF.get_sf_(lepton='Muons', pt_array=mu_pt, eta_array=mu_eta)
#print('Muons sf: ', mu_sfs)
#print('Expected: ', mu_expected_sfs)
#el_sfs = lepSF.get_sf_(lepton='Electrons', pt_array=el_pt, eta_array=el_eta)
#print('Computed: ', el_sfs)
#print('Expected: ', el_expected_sfs)
##set_trace()
