import coffea.nanoaod.nanoevents
import coffea.processor.dataframe
import awkward
from pdb import set_trace
import Utilities.prettyjson as prettyjson
import numpy as np
import os

btag_values = {
    'btagDeepB' : {
        'DEEPCSVLOOSE' : 0.2217,
        'DEEPCSVMEDIUM' : 0.6321,
        'DEEPCSVTIGHT' : 0.8953,
    },
    'btagDeepFlavB' : {
        'DEEPJETLOOSE' : 0.0614,
        'DEEPJETMEDIUM' : 0.3093,
        'DEEPJETTIGHT' : 0.7221
    }
}

jet_pars = prettyjson.loads(open('%s/cfg_files/cfg_pars.json' % os.environ['PROJECT_DIR']).read())['Jets']

valid_taggers = ['DEEPCSV', 'DEEPJET']
valid_WPs = ['LOOSE', 'MEDIUM', 'TIGHT']

def make_kin_cuts(jets):
    if isinstance(jets, awkward.array.base.AwkwardArray):
        pt_cut = (jets.pt >= jet_pars['ptmin'])
        leadpt_cut = (jets.pt.max() >= jet_pars['lead_ptmin'])
        eta_cut = (np.abs(jets.eta) <= jet_pars['etamax'])
        jets['kin_cuts'] = (pt_cut) & (leadpt_cut) & (eta_cut)

    else:
        raise ValueError("Only AwkwardArrays from NanoEvents are supported right now")

    return jets

def add_btag_disc(jets):

    if jet_pars['btagger'] not in valid_taggers:
        raise IOError("%s is not a supported b-tagger" % jet_pars['btagger'])
    if jet_pars['permutations']['tightb'] not in valid_WPs:
        raise IOError("%s is not a valid working point" % jet_pars['permutations']['tightb'])
    if jet_pars['permutations']['looseb'] not in valid_WPs:
        raise IOError("%s is not a valid working point" % jet_pars['permutations']['looseb'])

    #set_trace()
    if isinstance(jets, awkward.array.base.AwkwardArray):
        bdiscr = 'btagDeepB' if jet_pars['btagger'] == 'DEEPCSV' else 'btagDeepFlavB'
        wps = list(set([jet_pars['btagger']+jet_pars['permutations']['tightb'], jet_pars['btagger']+jet_pars['permutations']['looseb']]))
        for wp in wps:
            jets['BTAG_%s' % wp] = (jets[bdiscr] > btag_values[bdiscr][wp])

    else:
        raise ValueError("Only AwkwardArrays from NanoEvents are supported right now")

    return jets


def build_jets(jets):
    jets = make_kin_cuts(jets)
    jets = add_btag_disc(jets)

    return jets


def process_jets(df):
    print('Inside process_jets')

    if isinstance(df, coffea.nanoaod.nanoevents.NanoEvents):
        bdiscr = 'btagDeepB' if btagger == 'DEEPCSV' else 'btagDeepFlavB'
        wps = list(set([btagger+tightb, btagger+looseb]))
        set_trace()
        for wp in wps:
            df['Jet'][wp] = (df['Jet'][bdiscr] > btag_values[bdiscr][wp])
        set_trace()

    elif isinstance(df, coffea.processor.dataframe.LazyDataFrame):
        from coffea.analysis_objects import JaggedCandidateArray
    
        jets = JaggedCandidateArray.candidatesfromcounts(
            df['nJet'],
            pt=df['Jet_pt'],
            eta=df['Jet_eta'],
            phi=df['Jet_phi'],
            mass=df['Jet_mass'],
            DeepCSVb=df['Jet_btagDeepB'],
            DeepJetb=df['Jet_btagDeepFlavB'],
            Id=df['Jet_jetId'],
            cleanmask=df['Jet_cleanmask'],
        )
    
        #tight_id = (jets.Id == 3) # 3 for passing loose+tight for 2016, other years are different https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD
    
            ## add btag wps
        deepcsvL = (jets.DeepCSVb > 0.2217)
        deepcsvM = (jets.DeepCSVb > 0.6321)
        deepcsvT = (jets.DeepCSVb > 0.8953)
    
        deepjetL = (jets.DeepJetb > 0.0614)
        deepjetM = (jets.DeepJetb > 0.3093)
        deepjetT = (jets.DeepJetb > 0.7221)
    
        jets.add_attributes(
            #tightID = tight_id,
            deepcsvL = deepcsvL,
            deepcsvM = deepcsvM,
            deepcsvT = deepcsvT,
            deepjetL = deepjetL,
            deepjetM = deepjetM,
            deepjetT = deepjetT,
        )
    
        return jets

