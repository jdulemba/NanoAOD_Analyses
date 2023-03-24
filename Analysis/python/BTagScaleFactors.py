from Utilities.convert_btag_csv_file import convert_btag_csv_file
import numpy as np
from coffea.lookup_tools.dense_evaluated_lookup import dense_evaluated_lookup
from collections import defaultdict
nested_dict = lambda: defaultdict(nested_dict)
from pdb import set_trace
import Utilities.prettyjson as prettyjson
import os
from coffea.util import load
import awkward as ak

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
base_jobid = os.environ["base_jobid"]

btag_csvFiles = {
    "2016APV" : {
        #"DeepJet" : "DeepJet_106XUL16SF_mujets_used.csv",
        #"DeepCSV" : "DeepCSV_106XUL16SF_mujets_used.csv",
        "DeepJet" : "hig22013_wp_deepJet_UL2016preVFP_used.csv",
        #"DeepJet" : "DeepJet_106XUL16preVFPSF_v1_used.csv",
        "DeepCSV" : "DeepCSV_106XUL16preVFPSF_v1_used.csv",
    },
    "2016" : {
        #"DeepJet" : "DeepJet_106XUL16SF_mujets_used.csv",
        #"DeepCSV" : "DeepCSV_106XUL16SF_mujets_used.csv",
        "DeepJet" : "hig22013_wp_deepJet_UL2016postVFP_used.csv",
        #"DeepJet" : "DeepJet_106XUL16postVFPSF_v2_used.csv",
        "DeepCSV" : "DeepCSV_106XUL16postVFPSF_v2_used.csv",
    },
    "2017" : {
        #"DeepJet" : "DeepJet_106XUL17SF_WPonly_V2p1_mujets_used.csv",
        #"DeepCSV" : "DeepCSV_106XUL17SF_WPonly_V2p1_mujets_used.csv",
        "DeepJet" : "hig22013_wp_deepJet_UL2017_used.csv",
        #"DeepJet" : "wp_deepJet_106XUL17_v3_used.csv",
        "DeepCSV" : "wp_deepCSV_106XUL17_v3_used.csv",
    },
    "2018" : {
        #"DeepJet" : "DeepJet_106XUL18SF_WPonly_mujets_used.csv",
        #"DeepCSV" : "DeepCSV_106XUL18SF_WPonly_mujets_used.csv",
        "DeepJet" : "hig22013_wp_deepJet_UL2018_used.csv",
        #"DeepJet" : "wp_deepJet_106XUL18_v2_used.csv",
        "DeepCSV" : "wp_deepCSV_106XUL18_v2_used.csv",
    },
}

wp_lookup_dict = {
    "L" : "Loose",
    "M" : "Medium",
    "T" : "Tight",
    "R" : "Reshape"
}
wp_lookup = [
    "Loose",
    "Medium",
    "Tight",
    "Reshape"
]

flav_2_name_dict = {"5" : "B", "4" : "C", "0" : "UDSG"}
flav_2_name = ["B", "C", "UDSG"]

import awkward
def reshuffle_sf_dict(sf_dict, label_modifier = lambda x: x, algorithm = None):
    #set_trace()
    retval = nested_dict()
    for key, val in sf_dict.items():
        label, check = key
        if check != "dense_evaluated_lookup":
            raise ValueError(f"Value for label {label} is not a dense_evaluated_lookup")
        
        label = label_modifier(label)
        split = label.split("_")
        if len(split) == 5:
            algo, wp, source, sys, flav = tuple(split)
        elif len(split) == 6:
            algo, wp, source, sys, up_down, flav = tuple(split)
            sys = "_".join([sys, up_down])
        else:
            raise RuntimeError(f"This should not happen {label}")

        if (algo == "btagsf") and (algorithm is not None):
            algo = algorithm

        #set_trace()            
        retval[algo][source][wp_lookup_dict[str(wp)]]["_".join([flav_2_name_dict[str(flav)], sys])] = val
        #retval[algo][source][wp_lookup[int(wp)]]["_".join([flav_2_name[int(flav)], sys])] = val
    #set_trace()
    return retval

def recursive_compile(sf_dict):
    """compiles the values of the dict into coffea functions"""
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
        """SF computation according to method 1a of 
        https://twiki.cern.ch/twiki/bin/view/CMS/BTagSFMethods
        Inputs: csv, wp_key, eff_file, pattern
        csv: path to a b-tagging CSV file
        wp_key: a tuple of three elements containing (Algo name, SF method, WP name) 
        effs: dictionary containing the efficiencies for each flavour as dense_lookups"""
        #set_trace()
        parsed_csv = reshuffle_sf_dict(
            convert_btag_csv_file(csv),
            algorithm = wp_key[0]
            )
        self.sf_ = recursive_compile(parsed_csv[wp_key[0]][wp_key[1]][wp_key[2]])
        #set_trace()
        # FIXME: move to correlated/uncorrelated
        # Define, by hand, the proper correlation among taggers, 
        # somewhere unfortunately needs to be hardcoded by hand
        # tuple of names for UDSG, B, C
        self.schema_ = { 
            "central" : ("UDSG_central", "C_central", "B_central"),
            "bc_up" : ("UDSG_central", "C_up", "B_up"),
            "bc_up_correlated" : ("UDSG_central", "C_up_correlated", "B_up_correlated"),
            "bc_up_uncorrelated" : ("UDSG_central", "C_up_uncorrelated", "B_up_uncorrelated"),
            "bc_down" : ("UDSG_central", "C_down", "B_down"),
            "bc_down_correlated" : ("UDSG_central", "C_down_correlated", "B_down_correlated"),
            "bc_down_uncorrelated" : ("UDSG_central", "C_down_uncorrelated", "B_down_uncorrelated"),
            "l_up" : ("UDSG_up", "C_central", "B_central"),
            "l_up_correlated" : ("UDSG_up_correlated", "C_central", "B_central"),
            "l_up_uncorrelated" : ("UDSG_up_uncorrelated", "C_central", "B_central"),
            "l_down" : ("UDSG_down", "C_central", "B_central"),
            "l_down_correlated" : ("UDSG_down_correlated", "C_central", "B_central"),
            "l_down_uncorrelated" : ("UDSG_down_uncorrelated", "C_central", "B_central"),
            #"udsg_up" : ("UDSG_up", "C_central", "B_central"),
            #"udsg_up_correlated" : ("UDSG_up_correlated", "C_central", "B_central"),
            #"udsg_up_uncorrelated" : ("UDSG_up_uncorrelated", "C_central", "B_central"),
            #"udsg_down" : ("UDSG_down", "C_central", "B_central"),
            #"udsg_down_correlated" : ("UDSG_down_correlated", "C_central", "B_central"),
            #"udsg_down_uncorrelated" : ("UDSG_down_uncorrelated", "C_central", "B_central"),
            "bc_jes_up" : ("UDSG_central", "C_up_jes", "B_up_jes"),
            "bc_pileup_up" : ("UDSG_central", "C_up_pileup", "B_up_pileup"),
            "bc_statistic_up" : ("UDSG_central", "C_up_statistic", "B_up_statistic"),
            "bc_bfragmentation_up" : ("UDSG_central", "C_up_bfragmentation", "B_up_bfragmentation"),
            "bc_btempcorr_up" : ("UDSG_central", "C_up_btempcorr", "B_up_btempcorr"),
            "bc_cb_up" : ("UDSG_central", "C_up_cb", "B_up_cb"),
            "bc_cfragmentation_up" : ("UDSG_central", "C_up_cfragmentation", "B_up_cfragmentation"),
            "bc_cjets_up" : ("UDSG_central", "C_up_cjets", "B_up_cjets"),
            "bc_dmux_up" : ("UDSG_central", "C_up_dmux", "B_up_dmux"),
            "bc_gluonsplitting_up" : ("UDSG_central", "C_up_gluonsplitting", "B_up_gluonsplitting"),
            "bc_jetaway_up" : ("UDSG_central", "C_up_jetaway", "B_up_jetaway"),
            "bc_ksl_up" : ("UDSG_central", "C_up_ksl", "B_up_ksl"),
            "bc_l2c_up" : ("UDSG_central", "C_up_l2c", "B_up_l2c"),
            "bc_ltothers_up" : ("UDSG_central", "C_up_ltothers", "B_up_ltothers"),
            "bc_mudr_up" : ("UDSG_central", "C_up_mudr", "B_up_mudr"),
            "bc_mupt_up" : ("UDSG_central", "C_up_mupt", "B_up_mupt"),
            "bc_ptrel_up" : ("UDSG_central", "C_up_ptrel", "B_up_ptrel"),
            #"bc_type3_up" : ("UDSG_central", "C_up_type3", "B_up_type3"),

            "bc_jes_down" : ("UDSG_central", "C_down_jes", "B_down_jes"),
            "bc_pileup_down" : ("UDSG_central", "C_down_pileup", "B_down_pileup"),
            "bc_statistic_down" : ("UDSG_central", "C_down_statistic", "B_down_statistic"),
            "bc_bfragmentation_down" : ("UDSG_central", "C_down_bfragmentation", "B_down_bfragmentation"),
            "bc_btempcorr_down" : ("UDSG_central", "C_down_btempcorr", "B_down_btempcorr"),
            "bc_cb_down" : ("UDSG_central", "C_down_cb", "B_down_cb"),
            "bc_cfragmentation_down" : ("UDSG_central", "C_down_cfragmentation", "B_down_cfragmentation"),
            "bc_cjets_down" : ("UDSG_central", "C_down_cjets", "B_down_cjets"),
            "bc_dmux_down" : ("UDSG_central", "C_down_dmux", "B_down_dmux"),
            "bc_gluonsplitting_down" : ("UDSG_central", "C_down_gluonsplitting", "B_down_gluonsplitting"),
            "bc_jetaway_down" : ("UDSG_central", "C_down_jetaway", "B_down_jetaway"),
            "bc_ksl_down" : ("UDSG_central", "C_down_ksl", "B_down_ksl"),
            "bc_l2c_down" : ("UDSG_central", "C_down_l2c", "B_down_l2c"),
            "bc_ltothers_down" : ("UDSG_central", "C_down_ltothers", "B_down_ltothers"),
            "bc_mudr_down" : ("UDSG_central", "C_down_mudr", "B_down_mudr"),
            "bc_mupt_down" : ("UDSG_central", "C_down_mupt", "B_down_mupt"),
            "bc_ptrel_down" : ("UDSG_central", "C_down_ptrel", "B_down_ptrel"),
            #"bc_type3_down" : ("UDSG_central", "C_down_type3", "B_down_type3"),
            }

        self.eff_ = {
            "B"    : effs["bottom"],
            "C"    : effs["charm" ],
            "UDSG" : effs["light" ],
        }

    def match_flav_(self, light, charm, bottom, flav):
        """returns a np.array with the correct output matched
        according to the flavour"""
        ret = deepcopy(light)
        is_c = (flav == 4)
        is_b = (flav == 5)
        ret[is_c] = charm[is_c]
        ret[is_b] = bottom[is_b]
        return ret

    def efficiency_(self, pt, eta, flav):
        """"computes the efficiency under each 
        flavour assumption and then matches it"""
        eff_l = self.eff_["UDSG"](pt, eta)
        eff_c = self.eff_["C"](pt, eta)
        eff_b = self.eff_["B"](pt, eta)
        return self.match_flav_(eff_l, eff_c, eff_b, flav)

    def get_scale_factor(self, jets, passing_cut):
        """Starting from a jet collection and a string pointing to 
        the flag defining if the jet is b-tagged or not computes the 
        per-jet weight to be used. Supports only a single WP for the 
        moment"""
        # First of all flatten everything to make it easier to handle
        pt = ak.to_numpy(ak.flatten(jets.pt))
        eta = ak.to_numpy(ak.flatten(jets.eta))
        flav = ak.to_numpy(ak.flatten(jets.hadronFlavour))
        pass_wp = ak.to_numpy(ak.flatten(jets[passing_cut]))

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
                # protect against using prelim csv files that don't have UDSG wps
                flavour_sf_cache[lcb[i]] = np.ones(eta.size) if lcb[i] not in self.sf_.keys() else flavour_sf_cache.get(lcb[i], self.sf_[lcb[i]](eta, pt, pass_wp))
                
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

        # return the unflattened version of the ratio
        return {key : ak.unflatten(i/p_mc, ak.num(jets.pt))
                for key, i in p_data.items()}


## evts = NanoEvents.from_file("/afs/cern.ch/work/j/jdulemba/public/ttJets2016Nano_0.root")
## evts["Jet"]["DeepCSVMedium"] = evts["Jet"]["btagDeepB"] > 0.8
## 
## sf_computer = BTagSF(
##     csv = "DeepJet_2016LegacySF_V1.csv", 
##     wp_key = ("DeepJet", "used", "Medium"), 
##     eff_file = "htt_DeepJet_2016Legacy_j20l50MT40_1lep_DEEPJETMEDIUM_DEEPJETMEDIUM_3PJets_efficiencies.root", 
##     pattern = "{0}/DEEPJETMEDIUM_eff_3Jets"
## )
## sf_weight = sf_computer.get_scale_factor(evts["Jet"], "DeepCSVMedium")
