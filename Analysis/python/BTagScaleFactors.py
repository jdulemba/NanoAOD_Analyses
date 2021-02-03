from coffea.lookup_tools.csv_converters import convert_btag_csv_file
import numpy as np
from coffea.lookup_tools.dense_evaluated_lookup import dense_evaluated_lookup
from collections import defaultdict
nested_dict = lambda: defaultdict(nested_dict)
from pdb import set_trace
import Utilities.prettyjson as prettyjson
import os
from coffea.util import load

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
base_jobid = os.environ['base_jobid']

if base_jobid == 'NanoAODv6':
    btag_csvFiles = {
        '2016' : {
            'DeepJet' : 'DeepJet_2016LegacySF_V1_used.csv',
            'DeepCSV' : 'DeepCSV_2016LegacySF_V1_used.csv',
        },
        '2017' : {
            'DeepJet' : 'DeepJet_2017SF_V4_B_F_used.csv',
            'DeepCSV' : 'DeepCSV_2017SF_V5_B_F_used.csv',
        },
        '2018' : {
            'DeepJet' : 'DeepJet_2018SF_V1_used.csv',
            'DeepCSV' : 'DeepCSV_2018SF_V1_used.csv',
        },
    }
elif base_jobid == 'ULnanoAOD':
    btag_csvFiles = {
        #'2016' : {
        #    'DeepJet' : 'DeepJet_2016LegacySF_V1_used.csv',
        #    'DeepCSV' : 'DeepCSV_2016LegacySF_V1_used.csv',
        #},
        '2017' : {
            'DeepJet' : 'DeepJet_106XUL17SF_used.csv',
            'DeepCSV' : 'DeepCSV_106XUL17SF_used.csv',
        },
        '2018' : {
            'DeepJet' : 'DeepJet_106XUL18SF_WPonly_used.csv',
            'DeepCSV' : 'DeepCSV_106XUL18SF_WPonly_used.csv',
        },
    }
else:
    raise ValueError("%s not currently supported" % base_jobid)

wp_lookup = [
    'Loose',
    'Medium',
    'Tight',
    'Reshape'
]

flav_2_name = ['B', 'C', 'UDSG']

import awkward
def reshuffle_sf_dict(sf_dict, label_modifier = lambda x: x):
    retval = nested_dict()
    for key, val in sf_dict.items():
        label, check = key
        if check != 'dense_evaluated_lookup':
            raise ValueError(f'Value for label {label} is not a dense_evaluated_lookup')
        
        label = label_modifier(label)
        split = label.split('_')
        if len(split) == 5:
            algo, wp, source, sys, flav = tuple(split)
        elif len(split) == 6:
            algo, wp, source, sys, up_down, flav = tuple(split)
            sys = '_'.join([sys, up_down])
        else:
            raise RuntimeError(f'This should not happen {label}')
            
        retval[algo][source][wp_lookup[int(wp)]]['_'.join([flav_2_name[int(flav)], sys])] = val
    return retval

def recursive_compile(sf_dict):
    '''compiles the values of the dict into coffea functions'''
    retval = {}
    for key, val in sf_dict.items():
        if isinstance(val, dict):
            retval[key] = recursive_compile(val)
        else:
            retval[key] = dense_evaluated_lookup(*val)
    return retval

from copy import deepcopy
class BTagSF(object):
    def __init__(self, csv = None, wp_key = None, effs = None):
        '''SF computation according to method 1a of 
        https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods
        Inputs: csv, wp_key, eff_file, pattern
        csv: path to a b-tagging CSV file
        wp_key: a tuple of three elements containing (Algo name, SF method, WP name) 
        effs: dictionary containing the efficiencies for each flavour as dense_lookups'''
        parsed_csv = reshuffle_sf_dict(
            convert_btag_csv_file(csv)
            )
        self.sf_ = recursive_compile(parsed_csv[wp_key[0]][wp_key[1]][wp_key[2]])
        # FIXME: move to correlated/uncorrelated
        # Define, by hand, the proper correlation among taggers, 
        # somewhere unfortunately needs to be hardcoded by hand
        # tuple of names for UDSG, B, C
        self.schema_ = { 
            'central' : ('UDSG_central', 'C_central', 'B_central'),
            'bc_up' : ('UDSG_central', 'C_up', 'B_up'),
            'bc_down' : ('UDSG_central', 'C_down', 'B_down'),
            'udsg_up' : ('UDSG_up', 'C_central', 'B_central'),
            'udsg_down' : ('UDSG_down', 'C_central', 'B_central'),
            }

        self.eff_ = {
            'B'    : effs['bottom'],
            'C'    : effs['charm' ],
            'UDSG' : effs['light' ],
        }

    def match_flav_(self, light, charm, bottom, flav):
        '''returns a np.array with the correct output matched
        according to the flavour'''
        ret = deepcopy(light)
        is_c = (flav == 4)
        is_b = (flav == 5)
        ret[is_c] = charm[is_c]
        ret[is_b] = bottom[is_b]
        return ret

    def efficiency_(self, pt, eta, flav):
        ''''computes the efficiency under each 
        flavour assumption and then matches it'''
        eff_l = self.eff_['UDSG'](pt, eta)
        eff_c = self.eff_['C'](pt, eta)
        eff_b = self.eff_['B'](pt, eta)
        return self.match_flav_(eff_l, eff_c, eff_b, flav)

    def get_scale_factor(self, jets, passing_cut):
        '''Starting from a jet collection and a string pointing to 
        the flag defining if the jet is b-tagged or not computes the 
        per-jet weight to be used. Supports only a single WP for the 
        moment'''
        # First of all flatten everything to make it easier to handle
        pt = jets.pt.flatten()
        eta = jets.eta.flatten()
        flav = jets.hadronFlav.flatten()
        pass_wp = jets[passing_cut].flatten()

        # Get the MC efficiency
        eff = self.efficiency_(pt, eta, flav)
        # for each systematic/central value compute the proper SF
        # cache the SF values as sometimes they are repeated, there 
        # might also be systematic combinations that are never accessed
        # but pruning them at the beginning can be hard
        # use schema to define combinations, lcb is a tuple with the sf keys
        # for light, charm, bottom for each systematic
        flavour_sf_cache = {}
        scale_factors = {} # our final product
        for key, lcb in self.schema_.items(): 
            # populate cache if needed
            for i in range(3):
                flavour_sf_cache[lcb[i]] = flavour_sf_cache.get(
                    # for some reason there is an additional dimension, pass_wp has no effect
                    lcb[i], self.sf_[lcb[i]](eta, pt, pass_wp) 
                )
            scale_factors[key] = eff * self.match_flav_(
                flavour_sf_cache[lcb[0]],
                flavour_sf_cache[lcb[1]],
                flavour_sf_cache[lcb[2]],
                flav
            )
        
        # use SF and eff to compute p(data) and p(MC)
        p_data = {key : np.where(pass_wp, val, 1 - val) 
                  for key, val in scale_factors.items()}
        p_mc = np.where(pass_wp, eff, 1 - eff)

        # return the jagged version of the ratio
        return {key : awkward.JaggedArray.fromcounts(jets.pt.counts, i/p_mc)
                for key, i in p_data.items()}


## evts = NanoEvents.from_file('/afs/cern.ch/work/j/jdulemba/public/ttJets2016Nano_0.root')
## evts['Jet']['DeepCSVMedium'] = evts['Jet']['btagDeepB'] > 0.8
## 
## sf_computer = BTagSF(
##     csv = 'DeepJet_2016LegacySF_V1.csv', 
##     wp_key = ('DeepJet', 'used', 'Medium'), 
##     eff_file = 'htt_DeepJet_2016Legacy_j20l50MT40_1lep_DEEPJETMEDIUM_DEEPJETMEDIUM_3PJets_efficiencies.root', 
##     pattern = '{0}/DEEPJETMEDIUM_eff_3Jets'
## )
## sf_weight = sf_computer.get_scale_factor(evts['Jet'], 'DeepCSVMedium')
