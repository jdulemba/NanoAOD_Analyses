import coffea.processor.dataframe
import awkward
from pdb import set_trace
import Utilities.prettyjson as prettyjson
import numpy as np
import os

year = "2016"
btag_values = {}
btag_values["2016"] = {
    'btagDeepB' : {
        'DEEPCSVLOOSE' : 0.2217,
        'DEEPCSVMEDIUM': 0.6321,
        'DEEPCSVTIGHT' : 0.8953,
    },
    'btagDeepFlavB' : {
        'DEEPJETLOOSE' : 0.0614,
        'DEEPJETMEDIUM': 0.3093,
        'DEEPJETTIGHT' : 0.7221
    }
}
btag_values["2017"] = {
    'btagDeepB' : {
        'DEEPCSVLOOSE' : 0.1522,
        'DEEPCSVMEDIUM': 0.4941,
        'DEEPCSVTIGHT' : 0.8001,
    },
    'btagDeepFlavB' : {
        'DEEPJETLOOSE' : 0.0521,
        'DEEPJETMEDIUM': 0.3033,
        'DEEPJETTIGHT' : 0.7489
    }
}
btag_values["2018"] = {
    'btagDeepB' : {
        'DEEPCSVLOOSE' : 0.1241,
        'DEEPCSVMEDIUM': 0.4184,
        'DEEPCSVTIGHT' : 0.7527,
    },
    'btagDeepFlavB' : {
        'DEEPJETLOOSE' : 0.0494,
        'DEEPJETMEDIUM': 0.2770,
        'DEEPJETTIGHT' : 0.7264
    }
}

jet_pars = prettyjson.loads(open('%s/cfg_files/cfg_pars.json' % os.environ['PROJECT_DIR']).read())['Jets']

valid_taggers = ['DEEPCSV', 'DEEPJET']
valid_WPs = ['LOOSE', 'MEDIUM', 'TIGHT']

if jet_pars['btagger'] not in valid_taggers:
    raise IOError("%s is not a supported b-tagger" % jet_pars['btagger'])
if jet_pars['permutations']['tightb'] not in valid_WPs:
    raise IOError("%s is not a valid working point" % jet_pars['permutations']['tightb'])
if jet_pars['permutations']['looseb'] not in valid_WPs:
    raise IOError("%s is not a valid working point" % jet_pars['permutations']['looseb'])

bdiscr = 'btagDeepB' if jet_pars['btagger'] == 'DEEPCSV' else 'btagDeepFlavB'
wps = list(set([jet_pars['btagger']+jet_pars['permutations']['tightb'], jet_pars['btagger']+jet_pars['permutations']['looseb']]))



def make_kin_cuts(jets):
    if isinstance(jets, awkward.array.base.AwkwardArray):
        pt_cut = (jets.pt >= jet_pars['ptmin'])
        leadpt_cut = (jets.pt.max() >= jet_pars['lead_ptmin'])
        eta_cut = (np.abs(jets.eta) <= jet_pars['etamax'])
        kin_cuts = (pt_cut) & (leadpt_cut) & (eta_cut)

    else:
        raise ValueError("Only AwkwardArrays are supported")

    return kin_cuts


def add_btag_wps(jets, btagger, wps=[]):
    if btagger not in valid_taggers:
        raise IOError("%s is not a supported b-tagger" % btagger)
    bdiscr = 'btagDeepB' if btagger == 'DEEPCSV' else 'btagDeepFlavB'

    if isinstance(jets, awkward.array.base.AwkwardArray):
        for wp in wps:
            if 'BTAG_%s' % wp in jets.columns:
                continue
            if wp not in valid_WPs:
                raise IOError("%s is not a valid working point" % wp)
            jets['BTAG_%s' % wp] = (jets[bdiscr] > btag_values[year][bdiscr][wp])

    else:
        raise ValueError("Only AwkwardArrays are supported")

    return jets


def process_jets(df):

    if not isinstance(df, coffea.processor.dataframe.LazyDataFrame):
        raise IOError("This function only works for LazyDataFrame objects")

    from coffea.analysis_objects import JaggedCandidateArray
    
    Jet = JaggedCandidateArray.candidatesfromcounts(
        df['nJet'],
        pt=df['Jet_pt'],
        eta=df['Jet_eta'],
        phi=df['Jet_phi'],
        mass=df['Jet_mass'],
        btagDeepB=df['Jet_btagDeepB'],
        btagDeepFlavB=df['Jet_btagDeepFlavB'],
        Id=df['Jet_jetId'],
        cleanmask=df['Jet_cleanmask'],
        hadronFlav=df['Jet_hadronFlavour'],
    )

    #set_trace()
        ## add kinematic cuts
    Jet['Kin_Cuts'] = make_kin_cuts(Jet)

        ## add btag wps
    if len(wps) > 1:
        raise IOError("Only one btag wp supported right now")
    for wp in wps:
        Jet['BTAG_%s' % wp] = (Jet[bdiscr] > btag_values[year][bdiscr][wp])
    
    return Jet

