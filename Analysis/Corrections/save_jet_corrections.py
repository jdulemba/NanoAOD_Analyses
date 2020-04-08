from coffea.lookup_tools import extractor
from coffea.jetmet_tools import FactorizedJetCorrector, JetCorrectionUncertainty, JetTransformer, JetResolution, JetResolutionScaleFactor
import os
from pdb import set_trace
from coffea.util import save

proj_dir = os.environ['PROJECT_DIR']

jerfiles = {
    '2016' : {
        'JER' : ['Summer16_25nsV1_MC_PtResolution_AK4PFchs'],
        'JERSF' : ['Summer16_25nsV1_MC_SF_AK4PFchs'],
    },
    '2017' : {
        'JER' : ['Fall17_V3_MC_PtResolution_AK4PFchs'],
        'JERSF' : ['Fall17_V3_MC_SF_AK4PFchs'],
    },
    '2018' : {
        'JER' : ['Autumn18_V7_MC_PtResolution_AK4PFchs'],
        'JERSF' : ['Autumn18_V7_MC_SF_AK4PFchs'],
    },
}

jecfiles = {
    '2016' : {
        'Corr' : [
            'Summer16_07Aug2017_V11_MC_L1FastJet_AK4PFchs',
            'Summer16_07Aug2017_V11_MC_L2L3Residual_AK4PFchs',
            'Summer16_07Aug2017_V11_MC_L2Relative_AK4PFchs',
            'Summer16_07Aug2017_V11_MC_L2Residual_AK4PFchs',
            'Summer16_07Aug2017_V11_MC_L3Absolute_AK4PFchs',
        ],
        'Unc' : ['Summer16_07Aug2017_V11_MC_Uncertainty_AK4PFchs'],
    },
    '2017' : {
        'Corr' : [
            'Fall17_17Nov2017_V32_MC_L1FastJet_AK4PFchs',
            'Fall17_17Nov2017_V32_MC_L2L3Residual_AK4PFchs',
            'Fall17_17Nov2017_V32_MC_L2Relative_AK4PFchs',
            'Fall17_17Nov2017_V32_MC_L2Residual_AK4PFchs',
            'Fall17_17Nov2017_V32_MC_L3Absolute_AK4PFchs',
        ],
        'Unc' : ['Fall17_17Nov2017_V32_MC_Uncertainty_AK4PFchs'],
    },
    '2018' : {
        'Corr' : [
            'Autumn18_V19_MC_L1FastJet_AK4PFchs', 
            'Autumn18_V19_MC_L2L3Residual_AK4PFchs',
            'Autumn18_V19_MC_L2Relative_AK4PFchs',
            'Autumn18_V19_MC_L2Residual_AK4PFchs',
            'Autumn18_V19_MC_L3Absolute_AK4PFchs',
       ],
        'Unc' : ['Autumn18_V19_MC_Uncertainty_AK4PFchs'],
    },
}


Jetext = extractor()
for dirid in ['jec', 'junc', 'jr', 'jersf']:
    print(dirid)
    directory = '%s/inputs/data/%s' % (proj_dir, dirid)
    for filename in os.listdir(directory):
        if 'AK4PFchs' not in filename: continue
        filename = directory+'/'+filename
        Jetext.add_weight_sets(['* * '+filename])
        print('%s added to weights' % filename)
Jetext.finalize()
Jetevaluator = Jetext.make_evaluator()

jet_corrections = {
    '2016' : {},
    '2017' : {},
    '2018' : {},
}

for year in ['2016', '2017', '2018']:
    JECcorrector = FactorizedJetCorrector(**{name: Jetevaluator[name] for name in jecfiles[year]['Corr']})
    JECuncertainties = JetCorrectionUncertainty(**{name:Jetevaluator[name] for name in jecfiles[year]['Unc']})
    JER = JetResolution(**{name:Jetevaluator[name] for name in jerfiles[year]['JER']})
    JERsf = JetResolutionScaleFactor(**{name:Jetevaluator[name] for name in jerfiles[year]['JERSF']})
    Jet_transformer = JetTransformer(jec=JECcorrector,junc=JECuncertainties, jer = JER, jersf = JERsf)
    jet_corrections[year]['JEC'] = JECcorrector
    jet_corrections[year]['JECUnc'] = JECuncertainties
    jet_corrections[year]['JER'] = JER
    jet_corrections[year]['JERsf'] = JERsf
    jet_corrections[year]['JT'] = Jet_transformer

save(jet_corrections, '%s/Corrections/JetCorrections.coffea' % proj_dir)
print('\n%s/Corrections/JetCorrections.coffea written' % proj_dir)
