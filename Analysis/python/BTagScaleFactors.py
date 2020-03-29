from coffea.lookup_tools.csv_converters import convert_btag_csv_file
import numpy as np
from coffea.lookup_tools.dense_evaluated_lookup import dense_evaluated_lookup
from collections import defaultdict
nested_dict = lambda: defaultdict(nested_dict)
from coffea.lookup_tools.root_converters import convert_histo_root_file
from coffea.lookup_tools.dense_lookup import dense_lookup
from pdb import set_trace
import Utilities.prettyjson as prettyjson
import os

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']

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
    def __init__(self, csv = None, wp_key = None, eff_file = None, pattern = None):
        '''SF computation according to method 1a of 
        https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods
        Inputs: csv, wp_key, eff_file, pattern
        csv: path to a b-tagging CSV file
        wp_key: a tuple of three elements containing (Algo name, SF method, WP name) 
        eff_file: root file containing the efficiency histograms
        pattern: formattable string that accepts one parameter for flavour definition'''
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
            'udgs_up' : ('UDSG_up', 'C_central', 'B_central'),
            'udsg_down' : ('UDSG_down', 'C_central', 'B_central'),
            }

        effs = convert_histo_root_file(eff_file)
        self.eff_ = {
            'B'    : dense_lookup(*effs[(pattern.format('bottom'), 'dense_lookup')]),
            'C'    : dense_lookup(*effs[(pattern.format('charm' ), 'dense_lookup')]),
            'UDSG' : dense_lookup(*effs[(pattern.format('light' ), 'dense_lookup')]),
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
                    # for some reason there is an additional dimension
                    lcb[i], self.sf_[lcb[i]](pt, eta, pass_wp) 
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


def create_btag_sf_computer(year, njets):
    '''
    Method to create object that computes btag scale factors for a single jet multiplicity.

    Inputs:
        number of jets category, must be string (3 or 4+)
    Returns:
        Object that computes btag scale factors
    '''
    cfg_file = prettyjson.loads(open('%s/cfg_files/cfg_pars.json' % proj_dir).read())
    
    btag_wps = cfg_file['Jets']
    btagger = 'DeepJet' if btag_wps['btagger'] == 'DEEPJET' else 'DeepCSV' ## name in csv file
    wps = list(set([btag_wps['permutations']['tightb'], btag_wps['permutations']['looseb']]))
    if len(wps) > 1:
        raise IOError("Only single working point supported right now.")
    wp = wps[0].lower().capitalize()
    pattern_wp = (btagger+wp).upper()
    
    btag_files = cfg_file['BTagging'][year][btagger]
    csv_path = '/'.join([proj_dir, 'inputs', 'data', btag_files['csv_file']])
    eff_path = '/'.join([proj_dir, 'inputs', '%s_%s' % (year, jobid), 'INPUTS', btag_files['eff_file']])
    if not os.path.isfile(csv_path):
        raise IOError('BTagging csv file %s not found.' % csv_path)
    if not os.path.isfile(eff_path):
        raise IOError('BTagging efficiencies file %s not found.' % eff_path)

    if not (njets == '3' or njets == '4+'):
        raise IOError("Number of jets can only be 3 or 4+")

    njets_cat = '3Jets' if njets == '3' else '4PJets'
    sf_computer = BTagSF(
        csv = csv_path,
        wp_key = (btagger, 'used', wp),
        eff_file = eff_path,
        pattern = '{0}/%s_eff_%s' %  (pattern_wp, njets_cat)
    )

    print('BTag SF constructed for %s jets' % njets)
    return sf_computer

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
