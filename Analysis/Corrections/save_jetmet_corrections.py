from coffea.lookup_tools import extractor
from coffea.jetmet_tools import FactorizedJetCorrector, JetCorrectionUncertainty
from coffea.jetmet_tools import JECStack
from coffea.jetmet_tools import JetResolution, JetResolutionScaleFactor
from coffea.jetmet_tools import CorrectedJetsFactory, CorrectedMETFactory
import os
from pdb import set_trace
from coffea.util import save

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('--split_uncs', action='store_true', help='Use individual jec uncertainty sources file')
parser.add_argument('--test', action='store_true', help="Output txt file named 'test_jetmet.txt'")

args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
base_jobid = os.environ['base_jobid']

jet_type = 'AK4PFchs'
jec_levels_MC = ['L1FastJet', 'L2Relative', 'L3Absolute']
jec_levels_Data = ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']

split_by_era = {
    'jec' : True,
    'junc' : True,
    'jr' : False,
    'jersf' : False,
}

jecfiles = {
    'Unc' : 'UncertaintySources' if args.split_uncs else 'Uncertainty', # or UncertaintySources for splitting
    '2016APV' : {
        'tag' : 'Summer19UL16APV',
        'v' : 'V7',
        'DATA' : ['RunBCD', 'RunEF'],
        'eras' : ['BCD', 'EF'],
    },
    '2016' : {
        'tag' : 'Summer16_07Aug2017' if base_jobid == 'NanoAODv6' else 'Summer19UL16',
        'v' : 'V11' if base_jobid == 'NanoAODv6' else 'V7',
        'DATA' : ['BCD', 'EF', 'GH'] if base_jobid == 'NanoAODv6' else ['RunFGH'],
        'eras' : ['BCD', 'EF', 'GH'] if base_jobid == 'NanoAODv6' else ['FGH'],
    },
    '2017' : {
        'tag' : 'Fall17_17Nov2017' if base_jobid == 'NanoAODv6' else 'Summer19UL17',
        'v' : 'V32' if base_jobid == 'NanoAODv6' else 'V5',
        'DATA' : ['B', 'C', 'DE', 'F'] if base_jobid == 'NanoAODv6' else ['RunB', 'RunC', 'RunD', 'RunE', 'RunF'],
        'eras' : ['B', 'C', 'DE', 'F'] if base_jobid == 'NanoAODv6' else ['B', 'C', 'D', 'E', 'F'],
    },
    '2018' : {
        'tag' : 'Autumn18' if base_jobid == 'NanoAODv6' else 'Summer19UL18',
        'v' : 'V19' if base_jobid == 'NanoAODv6' else 'V5',
        'DATA' : ['RunA', 'RunB', 'RunC', 'RunD'],
        'eras' : ['A', 'B', 'C', 'D'],
    },
}

jerfiles = {
    '2016APV' : {
        'tag' : 'Summer20UL16APV_JRV3',
        'JER' : 'PtResolution',
        'JERSF' : 'SF',
    },
    '2016' : {
        'tag' : 'Summer16_25nsV1' if base_jobid == 'NanoAODv6' else 'Summer20UL16_JRV3',
        'JER' : 'PtResolution',
        'JERSF' : 'SF',
    },
    '2017' : {
        'tag' : 'Fall17_V3' if base_jobid == 'NanoAODv6' else 'Summer19UL17_JRV2',
        'JER' : 'PtResolution',
        'JERSF' : 'SF',
    },
    '2018' : {
        'tag' : 'Autumn18_V7' if base_jobid == 'NanoAODv6' else 'Summer19UL18_JRV2',
        'JER' : 'PtResolution',
        'JERSF' : 'SF',
    },
}


def make_name_map(jstack, isMC):
    # values are hardcoded based on names in IDJet
    name_map = jstack.blank_name_map
        # set name_map keys needed for CorrectedJetsFactory
    name_map['JetPt'] = 'pt'
    name_map['JetMass'] = 'mass'
    name_map['JetEta'] = 'eta'
    name_map['JetA'] = 'area'
    if isMC: name_map['ptGenJet'] = 'pt_gen'
    name_map['ptRaw'] = 'pt_raw'
    name_map['massRaw'] = 'mass_raw'
    name_map['Rho'] = 'rho'
        # set name_map keys needed for CorrectedMETFactory
    name_map['METpt'] = 'pt'
    name_map['METphi'] = 'phi'
    name_map['JetPhi'] = 'phi'
    name_map['UnClusteredEnergyDeltaX'] = 'MetUnclustEnUpDeltaX'
    name_map['UnClusteredEnergyDeltaY'] = 'MetUnclustEnUpDeltaY'

    return name_map


years_to_run = ['2016APV', '2016', '2017', '2018'] if base_jobid == 'ULnanoAOD' else ['2016', '2017', '2018']
if args.test: years_to_run = ['2018']

Jetext = extractor()
for dirid in ['jec', 'junc', 'jr', 'jersf']:
    for dtype in ['DATA', 'MC']:
        for year in years_to_run:
            runs = jecfiles[year][dtype] if ( (dtype == 'DATA') and (split_by_era[dirid]) ) else ['']
            for run in runs:
                print(dirid, dtype, year, run)
                directory = os.path.join(proj_dir, 'inputs', 'data', base_jobid, dirid, dtype, year, run)
                for filename in os.listdir(directory):
                    if jet_type not in filename: continue
                    fname = os.path.join(directory, filename)
                    Jetext.add_weight_sets(['* * '+fname])
                    print('%s added to weights' % fname)
Jetext.finalize()
Jetevaluator = Jetext.make_evaluator()


jet_corrections = {year:{'MC' : {}, 'DATA' : {}} for year in years_to_run}
for year in years_to_run:
    jec_mc_tag = '_'.join([jecfiles[year]['tag'], jecfiles[year]['v'], 'MC'])
    jer_tag = jerfiles[year]['tag']

        # DATA
    data_JER = JetResolution(**{name:Jetevaluator[name] for name in ['%s_DATA_%s_%s' % (jer_tag, jerfiles[year]['JER'], jet_type)]})
    data_JERsf = JetResolutionScaleFactor(**{name:Jetevaluator[name] for name in ['%s_DATA_%s_%s' % (jer_tag, jerfiles[year]['JERSF'], jet_type)]})
    for idx, era in enumerate(jecfiles[year]['eras']):
        jec_data_tag = '_'.join([jecfiles[year]['tag']+jecfiles[year]['DATA'][idx], jecfiles[year]['v'], 'DATA']) if '_' in jecfiles[year]['tag'] else '_'.join([jecfiles[year]['tag'], jecfiles[year]['DATA'][idx], jecfiles[year]['v'], 'DATA'])
            # get jec, junc, jr, jersf
        data_JECcorrector = FactorizedJetCorrector(**{name: Jetevaluator[name] for name in ['%s_%s_%s' % (jec_data_tag, level, jet_type) for level in jec_levels_Data]})
        data_JECuncertainties = JetCorrectionUncertainty(**{name:Jetevaluator[name] for name in Jetevaluator.keys() if name.startswith('%s_%s_%s' % (jec_data_tag, jecfiles['Unc'], jet_type))})
            # make JEC stack of all corrections
        data_JECStack = JECStack({}, jec=data_JECcorrector, junc=data_JECuncertainties)
            # make jet and met factory
        data_name_map = make_name_map(data_JECStack, isMC=False)
        data_jet_factory = CorrectedJetsFactory(data_name_map, data_JECStack)
        data_met_factory = CorrectedMETFactory(data_name_map)
        jet_corrections[year]['DATA'][era] = {
            'JetsFactory' : data_jet_factory,
            'METFactory' : data_met_factory,
        }

        # MC
            # get jec, junc, jr, jersf
    MC_JECcorrector = FactorizedJetCorrector(**{name: Jetevaluator[name] for name in ['%s_%s_%s' % (jec_mc_tag, level, jet_type) for level in jec_levels_MC]})
    MC_JECuncertainties = JetCorrectionUncertainty(**{name:Jetevaluator[name] for name in Jetevaluator.keys() if name.startswith('%s_%s_%s' % (jec_mc_tag, jecfiles['Unc'], jet_type))})
    MC_JER = JetResolution(**{name:Jetevaluator[name] for name in ['%s_MC_%s_%s' % (jer_tag, jerfiles[year]['JER'], jet_type)]})
    MC_JERsf = JetResolutionScaleFactor(**{name:Jetevaluator[name] for name in ['%s_MC_%s_%s' % (jer_tag, jerfiles[year]['JERSF'], jet_type)]})
        # make JEC stack of all corrections
    MC_JECStack = JECStack({}, jec=MC_JECcorrector, junc=MC_JECuncertainties, jer=MC_JER, jersf=MC_JERsf)
        # make jet and met factory
    MC_name_map = make_name_map(MC_JECStack, isMC=True)
    MC_jet_factory = CorrectedJetsFactory(MC_name_map, MC_JECStack)
    MC_met_factory = CorrectedMETFactory(MC_name_map)
    jet_corrections[year]['MC'] = {
        'JetsFactory' : MC_jet_factory,
        'METFactory' : MC_met_factory,
    }
    print('Jet corrections for %s saved' % year)


fname = os.path.join(proj_dir, 'Corrections', base_jobid, 'JetMETCorrections_UncSources.coffea') if args.split_uncs else os.path.join(proj_dir, 'Corrections', base_jobid, 'JetMETCorrections.coffea')
if args.test: fname = os.path.join(proj_dir, 'test_jetmet.coffea')

save(jet_corrections, fname)
print('\n%s written' % fname)
