import awkward as ak
from pdb import set_trace
import Utilities.prettyjson as prettyjson
import numpy as np
import os

btag_values = {}
if os.environ['base_jobid'] == 'NanoAODv6':
    btag_values["2016"] = {
        'btagDeepB' : {
            'DeepCSVLoose' : 0.2217,
            'DeepCSVMedium': 0.6321,
            'DeepCSVTight' : 0.8953,
        },
        'btagDeepFlavB' : {
            'DeepJetLoose' : 0.0614,
            'DeepJetMedium': 0.3093,
            'DeepJetTight' : 0.7221,
        }
    }
    btag_values["2017"] = {
        'btagDeepB' : {
            'DeepCSVLoose' : 0.1522,
            'DeepCSVMedium': 0.4941,
            'DeepCSVTight' : 0.8001,
        },
        'btagDeepFlavB' : {
            'DeepJetLoose' : 0.0521,
            'DeepJetMedium': 0.3033,
            'DeepJetTight' : 0.7489,
        }
    }
    btag_values["2018"] = {
        'btagDeepB' : {
            'DeepCSVLoose' : 0.1241,
            'DeepCSVMedium': 0.4184,
            'DeepCSVTight' : 0.7527,
        },
        'btagDeepFlavB' : {
            'DeepJetLoose' : 0.0494,
            'DeepJetMedium': 0.2770,
            'DeepJetTight' : 0.7264,
        }
    }
    
elif os.environ['base_jobid'] == 'ULnanoAOD':
        # 2016(APV) values are still old ones
    btag_values["2016APV"] = {
        'btagDeepB' : {
            'DeepCSVLoose' : 0.2217,
            'DeepCSVMedium': 0.6321,
            'DeepCSVTight' : 0.8953,
        },
        'btagDeepFlavB' : {
            'DeepJetLoose' : 0.0614,
            'DeepJetMedium': 0.3093,
            'DeepJetTight' : 0.7221
        }
    }
    btag_values["2016"] = {
        'btagDeepB' : {
            'DeepCSVLoose' : 0.2217,
            'DeepCSVMedium': 0.6321,
            'DeepCSVTight' : 0.8953,
        },
        'btagDeepFlavB' : {
            'DeepJetLoose' : 0.0614,
            'DeepJetMedium': 0.3093,
            'DeepJetTight' : 0.7221
        }
    }
    btag_values["2017"] = {
        'btagDeepB' : {
            'DeepCSVLoose' : 0.1355,
            'DeepCSVMedium': 0.4506,
            'DeepCSVTight' : 0.7738,
        },
        'btagDeepFlavB' : {
            'DeepJetLoose' : 0.0532,
            'DeepJetMedium': 0.3040,
            'DeepJetTight' : 0.7476,
        }
    }
    btag_values["2018"] = {
        'btagDeepB' : {
            'DeepCSVLoose' : 0.1208,
            'DeepCSVMedium': 0.4168,
            'DeepCSVTight' : 0.7665,
        },
        'btagDeepFlavB' : {
            'DeepJetLoose' : 0.0490,
            'DeepJetMedium': 0.2783,
            'DeepJetTight' : 0.7100,
        }
    }

else:
    raise ValueError("base_jobid not set")
    
jet_pars = prettyjson.loads(open(os.path.join(os.environ['PROJECT_DIR'], 'cfg_files', 'cfg_pars_%s.json' % os.environ['jobid'])).read())['Jets']

valid_taggers = ['DeepCSV', 'DeepJet']
valid_WPs = ['Loose', 'Medium', 'Tight']

if jet_pars['btagger'] not in valid_taggers:
    raise IOError("%s is not a supported b-tagger" % jet_pars['btagger'])
if jet_pars['permutations']['tightb'] not in valid_WPs:
    raise IOError("%s is not a valid working point" % jet_pars['permutations']['tightb'])
if jet_pars['permutations']['looseb'] not in valid_WPs:
    raise IOError("%s is not a valid working point" % jet_pars['permutations']['looseb'])


def make_pt_eta_cuts(jets):
    pt_cut = (jets['pt'] >= jet_pars['ptmin'])
    eta_cut = (np.abs(jets['eta']) <= jet_pars['etamax'])
    return (pt_cut & eta_cut)

def make_leadjet_pt_cut(jets):
    leadpt_cut = (ak.max(jets['pt'], axis=1) >= jet_pars['lead_ptmin'])
    return leadpt_cut

def process_jets(events, year, corrections=None):

    jets = events['Jet']
    jets['pt_raw'] = (1 - jets['rawFactor']) * jets['pt']
    jets['mass_raw'] = (1 - jets['rawFactor']) * jets['mass']
    jets['rho'] = ak.broadcast_arrays(events.fixedGridRhoFastjetAll, jets.pt)[0]
    if not events.metadata['dataset'].startswith('data_Single'): jets['pt_gen'] = ak.values_astype(ak.fill_none(jets.matched_gen.pt, 0), np.float32)

       ## add btag wps
    for bdiscr in btag_values[year].keys():
        for wp in btag_values[year][bdiscr].keys():
            jets[wp] = (jets[bdiscr] > btag_values[year][bdiscr][wp])

        ## apply jet corrections
    if (jet_pars['applyJER'] == 1) and corrections is not None:
        if events.metadata['dataset'].startswith('data_Single'):
            era = [key for key in corrections['DATA'].keys() if events.metadata['dataset'].split(year)[-1] in key]
            if ('2016' in year) and (('Bv2' in events.metadata['dataset']) or ('C' in events.metadata['dataset']) or ('D' in events.metadata['dataset'])): era = ['BCD']
            if ('2016' in year) and (('E' in events.metadata['dataset']) or ('F' in events.metadata['dataset'])): era = ['EF']
            if ('2016' in year) and (('G' in events.metadata['dataset']) or ('H' in events.metadata['dataset'])): era = ['GH']
            #if ('2016' in year) and ('Bv2' in events.metadata['dataset']): era = ['BCD']
            if len(era) != 1: raise ValueError("Only one era should be used for %s" % events.metadata['dataset'])
            jet_factory = corrections['DATA'][era[0]]['JetsFactory']
            met_factory = corrections['DATA'][era[0]]['METFactory']

        else:
            jet_factory = corrections['MC']['JetsFactory']
            met_factory = corrections['MC']['METFactory']

        events_cache = events.caches[0]
        corrected_jets = jet_factory.build(jets, lazy_cache=events_cache)
        corrected_met = met_factory.build(events['MET'], corrected_jets, lazy_cache=events_cache)

    return corrected_jets, corrected_met
