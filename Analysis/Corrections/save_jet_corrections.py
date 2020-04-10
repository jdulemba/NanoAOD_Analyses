from coffea.lookup_tools import extractor
from coffea.jetmet_tools import FactorizedJetCorrector, JetCorrectionUncertainty, JetTransformer, JetResolution, JetResolutionScaleFactor
import os
from pdb import set_trace
from coffea.util import save

proj_dir = os.environ['PROJECT_DIR']

def tag_to_DataTag(tag, era):
    campaign, reco_ver, version = tag.split('_')
    reco_ver = reco_ver+era
    return '_'.join([campaign, reco_ver, version])

jet_type = 'AK4PFchs'
jec_levels_MC = ['L1FastJet', 'L2Relative', 'L3Absolute']
jec_levels_Data = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']

jecfiles = {
    'Unc' : 'Uncertainty', # or UncertaintySources for splitting
    '2016' : {
        'tag' : 'Summer16_07Aug2017_V11',
        'DATA' : ['BCD', 'EF', 'GH']
    },
    '2017' : {
        'tag' : 'Fall17_17Nov2017_V32',
        'DATA' : ['B', 'C', 'DE', 'F'],
    },
    '2018' : {
        'tag' : 'Autumn18_V19',
        'DATA' : ['A', 'B', 'C', 'D']
    },
}

jerfiles = {
    '2016' : {
        'tag' : 'Summer16_25nsV1',
        'JER' : 'PtResolution',
        'JERSF' : 'SF',
    },
    '2017' : {
        'tag' : 'Fall17_V3',
        'JER' : 'PtResolution',
        'JERSF' : 'SF',
    },
    '2018' : {
        'tag' : 'Autumn18_V7',
        'JER' : 'PtResolution',
        'JERSF' : 'SF',
    },
}

Jetext = extractor()
for dirid in ['jec', 'junc', 'jr', 'jersf']:
    for dtype in ['MC']:
    #for dtype in ['DATA', 'MC']:
        print(dirid, dtype)
        #if dirid == 'jersf' and dtype == 'DATA': set_trace()#continue
        directory = '%s/inputs/data/%s/%s' % (proj_dir, dirid, dtype)
        for filename in os.listdir(directory):
            if 'AK4PFchs' not in filename: continue
            filename = directory+'/'+filename
            Jetext.add_weight_sets(['* * '+filename])
            print('%s added to weights' % filename)
Jetext.finalize()
Jetevaluator = Jetext.make_evaluator()

jet_corrections = {
    '2016' : {
        'MC' : {},
#        'DATA': {},
    },
    '2017' : {
        'MC' : {},
#        'DATA': {},
    },
    '2018' : {
        'MC' : {},
#        'DATA': {},
    },
}

for year in ['2016', '2017', '2018']:
    jec_tag = jecfiles[year]['tag']
    jer_tag = jerfiles[year]['tag']

    #for era in jecfiles[year]['DATA']:
    #    data_tag = tag_to_DataTag(jec_tag, era)
    #    JECcorrector = FactorizedJetCorrector(**{name: Jetevaluator[name] for name in ['%s_DATA_%s_%s' % (data_tag, level, jet_type) for level in jec_levels_Data]})
    #    JECuncertainties = JetCorrectionUncertainty(**{name:Jetevaluator[name] for name in ['%s_DATA_%s_%s' % (data_tag, jecfiles['Unc'], jet_type)]})
    #    JER = JetResolution(**{name:Jetevaluator[name] for name in ['%s_DATA_%s_%s' % (jer_tag, jerfiles[year]['JER'], jet_type)]})
    #    JERsf = JetResolutionScaleFactor(**{name:Jetevaluator[name] for name in ['%s_DATA_%s_%s' % (jer_tag, jerfiles[year]['JERSF'], jet_type)]})
    #    Jet_transformer = JetTransformer(jec=JECcorrector, junc=JECuncertainties, jer=JER, jersf=JERsf)
    #    jec_corrections[year]['DATA'][era] = {
    #        'JEC' : JECcorrector,
    #        'JECUnc' : JECuncertainties,
    #        'JER' : JER,
    #        'JERsf' : JERsf,
    #        'JT' : Jet_transformer,
    #    }

    JECcorrector = FactorizedJetCorrector(**{name: Jetevaluator[name] for name in ['%s_MC_%s_%s' % (jec_tag, level, jet_type) for level in jec_levels_MC]})
    JECuncertainties = JetCorrectionUncertainty(**{name:Jetevaluator[name] for name in ['%s_MC_%s_%s' % (jec_tag, jecfiles['Unc'], jet_type)]})
    JER = JetResolution(**{name:Jetevaluator[name] for name in ['%s_MC_%s_%s' % (jer_tag, jerfiles[year]['JER'], jet_type)]})
    JERsf = JetResolutionScaleFactor(**{name:Jetevaluator[name] for name in ['%s_MC_%s_%s' % (jer_tag, jerfiles[year]['JERSF'], jet_type)]})
    Jet_transformer = JetTransformer(jec=JECcorrector, junc=JECuncertainties, jer=JER, jersf=JERsf)
    jet_corrections[year]['MC'] = {
        'JEC' : JECcorrector,
        'JECUnc' : JECuncertainties,
        'JER' : JER,
        'JERsf' : JERsf,
        'JT' : Jet_transformer,
    }

#set_trace()
save(jet_corrections, '%s/Corrections/JetCorrections.coffea' % proj_dir)
print('\n%s/Corrections/JetCorrections.coffea written' % proj_dir)
