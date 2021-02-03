import coffea.processor.dataframe
import awkward
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
    #btag_values["2016"] = {
    #    'btagDeepB' : {
    #        'DeepCSVLoose' : 0.2217,
    #        'DeepCSVMedium': 0.6321,
    #        'DeepCSVTight' : 0.8953,
    #    },
    #    'btagDeepFlavB' : {
    #        'DeepJetLoose' : 0.0614,
    #        'DeepJetMedium': 0.3093,
    #        'DeepJetTight' : 0.7221
    #    }
    #}
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
    
jet_pars = prettyjson.loads(open('%s/cfg_files/cfg_pars_%s.json' % (os.environ['PROJECT_DIR'], os.environ['jobid'])).read())['Jets']

valid_taggers = ['DeepCSV', 'DeepJet']
valid_WPs = ['Loose', 'Medium', 'Tight']

if jet_pars['btagger'] not in valid_taggers:
    raise IOError("%s is not a supported b-tagger" % jet_pars['btagger'])
if jet_pars['permutations']['tightb'] not in valid_WPs:
    raise IOError("%s is not a valid working point" % jet_pars['permutations']['tightb'])
if jet_pars['permutations']['looseb'] not in valid_WPs:
    raise IOError("%s is not a valid working point" % jet_pars['permutations']['looseb'])


_signature_map = {
    'JetPt': 'pt',
    'JetEta': 'eta',
    'Rho': 'rho',
    'JetA': 'area'
}

#bdiscr = 'btagDeepB' if jet_pars['btagger'] == 'DeepCSV' else 'btagDeepFlavB'
#wps = list(set([jet_pars['btagger']+jet_pars['permutations']['tightb'], jet_pars['btagger']+jet_pars['permutations']['looseb']]))
#wps = [''.join(wp) for wp in itertools.product(valid_taggers, valid_WPs)]


def make_pt_eta_cuts(jets):
    if isinstance(jets, awkward.array.base.AwkwardArray):
        pt_cut = (jets.pt >= jet_pars['ptmin'])
        eta_cut = (np.abs(jets.eta) <= jet_pars['etamax'])
        kin_cuts = (pt_cut & eta_cut)

    else:
        raise ValueError("Only AwkwardArrays are supported")

    return kin_cuts

def make_leadjet_pt_cut(jets):
    if isinstance(jets, awkward.array.base.AwkwardArray):
        leadpt_cut = (jets.pt.max() >= jet_pars['lead_ptmin'])
    else:
        raise ValueError("Only AwkwardArrays are supported")

    return leadpt_cut

def HEM_15_16_issue(jets):
    if isinstance(jets, awkward.array.base.AwkwardArray):
        hem_region = ((jets.eta > -3.2) & (jets.eta < -1.3) & (jets.phi > -1.57) & (jets.phi < -0.87))
    else:
        raise ValueError("Only AwkwardArrays are supported")

    return hem_region


def get_ptGenJet(jets, genjets, dr_max, pt_max_factor):
    '''
    requirements defined here: https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution#Smearing_procedures
    Procedure based off this: https://github.com/cms-sw/cmssw/blob/CMSSW_8_0_25/PhysicsTools/PatUtils/interface/SmearedJetProducerT.h#L59-L87
    '''
        # At each jet 1 location is the index of the genJet that it matched best with
        # <<<<important>>>> selves without a match will get a -1 to preserve counts structure
    matched_genJet_inds = jets.argmatch(genjets, deltaRCut=dr_max/2, deltaPtCut=(pt_max_factor*jets.JER))
        # initialize ptGenJet to zeros, must use jets not genjets for shape!
    ptGenJet = jets.pt.zeros_like()
        # find negation of events in which njets > 0 but ngenjets == 0 (results in empty sections of jagged arrays)
    valid_inds = ~((jets.counts > 0)  & (genjets.counts == 0))
        # set ptGenJet for jets that have corresponding genjet matches
    ptGenJet[valid_inds][matched_genJet_inds[valid_inds] != -1] = genjets[valid_inds][matched_genJet_inds[valid_inds][matched_genJet_inds[valid_inds] != -1]].pt
    jets['ptGenJet'] =  ptGenJet


def process_jets(df, year, corrections=None, shift=None):

    if not isinstance(df, coffea.processor.dataframe.LazyDataFrame):
        raise IOError("This function only works for LazyDataFrame objects")

    from coffea.analysis_objects import JaggedCandidateArray
    
    Jet = JaggedCandidateArray.candidatesfromcounts(
        counts=df['nJet'],
        pt=df['Jet_pt'],
        eta=df['Jet_eta'],
        phi=df['Jet_phi'],
        mass=df['Jet_mass'],
        btagDeepB=df['Jet_btagDeepB'],
        btagDeepFlavB=df['Jet_btagDeepFlavB'],
        Id=df['Jet_jetId'],
        #cleanmask=df['Jet_cleanmask'],
        #genJetIdx=df['Jet_genJetIdx'],
        area=df['Jet_area'],
        rawFactor=df['Jet_rawFactor'],
    )

    Jet['rho'] = Jet.pt.ones_like()*df['fixedGridRhoFastjetAll']
    Jet['ptRaw'] = Jet.pt*(1.-Jet['rawFactor'])
    Jet['massRaw'] = Jet.mass*(1.-Jet['rawFactor'])

       ## add btag wps
    for bdiscr in btag_values[year].keys():
        for wp in btag_values[year][bdiscr].keys():
            Jet[wp] = (Jet[bdiscr] > btag_values[year][bdiscr][wp])

        ## apply jet corrections
    if (jet_pars['applyJER'] == 1) and corrections is not None:
        if df.dataset.startswith('data_Single'):
            era = [key for key in corrections['DATA'].keys() if df.dataset.split(year)[-1] in key]
            if (year == '2016') and ('Bv2' in df.dataset): era = ['BCD']
            if len(era) != 1: raise ValueError("Only one era should be used for %s" % df.dataset)
            JER = corrections['DATA'][era[0]]['JER']
            Jet['JER'] = JER.getResolution(JetEta=Jet.eta, JetPt=Jet.pt, Rho=Jet.rho)
            Jet_transformer = corrections['DATA'][era[0]]['JT']

        else:
            Jet['hadronFlav'] = awkward.JaggedArray.fromcounts(Jet.counts, df['Jet_hadronFlavour'])
                ## apply JER
                    ## create gen jets for matching
            import python.GenParticleSelector as genpsel
            df['genJets'] = genpsel.process_genJets(df)

                # match jets to genJets to get ptGenJet
            JER = corrections['MC']['JER']
            Jet['JER'] = JER.getResolution(JetEta=Jet.eta, JetPt=Jet.pt, Rho=Jet.rho)
            get_ptGenJet(Jet, df['genJets'], dr_max=0.4, pt_max_factor=3)
            Jet_transformer = corrections['MC']['JT']
            #set_trace()

        import python.JetMET_corrections as JETMET_corr
        JETMET_corr.transform(Jet_transformer, jet=Jet, met=df['MET'], shift=shift)
        #Jet_transformer.transform(Jet, met=df['MET'])

    return Jet

