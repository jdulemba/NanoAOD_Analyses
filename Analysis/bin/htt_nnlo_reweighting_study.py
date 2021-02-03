#!/usr/bin/env python

from coffea import hist
from coffea.util import save, load
import coffea.processor as processor
from pdb import set_trace
import os, sys
import python.ObjectSelection as objsel
import Utilities.plot_tools as plt_tools
import python.MCWeights as MCWeights
import numpy as np
import Utilities.prettyjson as prettyjson
import coffea.lumi_tools.lumi_tools as lumi_tools
import Utilities.make_variables as make_vars
from python.IDJet import btag_values as btag_values
import python.GenParticleSelector as genpsel
import python.TTGenMatcher as ttmatcher
import python.TTPermutator as ttpermutator
from python.Permutations import compare_matched_best_perms
import Utilities.systematics as systematics

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
base_jobid = os.environ['base_jobid']
analyzer = 'htt_nnlo_reweighting_study'

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('fset', type=str, help='Fileset dictionary (in string form) to be used for the processor')
parser.add_argument('year', choices=['2016', '2017', '2018'], help='Specify which year to run over')
parser.add_argument('outfname', type=str, help='Specify output filename, including directory and file extension')
parser.add_argument('--debug', action='store_true', help='Uses iterative_executor for debugging purposes, otherwise futures_excutor will be used (faster)')

args = parser.parse_args()

# convert input string of fileset dictionary to actual dictionary
fdict = (args.fset).replace("\'", "\"")
fileset = prettyjson.loads(fdict)

if len(fileset.keys()) > 1:
    raise ValueError("Only one topology run at a time in order to determine which corrections and systematics to run")
samplename = list(fileset.keys())[0]

# get dataset classification, used for corrections/systematics
    ## specify ttJets samples
if args.year == '2016':
    Nominal_ttJets = ['ttJets_PS']#, 'ttJets']
else:
    Nominal_ttJets = ['ttJetsSL', 'ttJetsHad', 'ttJetsDiLep']

isTTbar_ = samplename.startswith('ttJets')
isNominal_ttJets_ = samplename in Nominal_ttJets
if not isNominal_ttJets_:
    raise ValueError("Only nominal ttbar events should be run")

## init tt probs for likelihoods
ttpermutator.year_to_run(year=args.year)

## load lumimask for data and corrections for event weights
cfg_pars = prettyjson.loads(open(os.path.join(proj_dir, 'cfg_files', 'cfg_pars_%s.json' % jobid)).read())
pu_correction = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'MC_PU_Weights.coffea'))
lepSF_correction = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'leptonSFs.coffea'))
jet_corrections = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'JetCorrections.coffea'))[args.year]
alpha_fname = os.path.join(proj_dir, 'Corrections', base_jobid, 'alpha_correction_%s.coffea' % jobid) if os.path.isfile(os.path.join(proj_dir, 'Corrections', base_jobid, 'alpha_correction_%s.coffea' % jobid)) else os.path.join(proj_dir, 'Corrections', base_jobid, 'alpha_correction_%s.coffea' % base_jobid)
alpha_corrections = load(alpha_fname)[args.year]['E']['All_2D'] # E, All_2D determined by post application plots
nnlo_reweighting = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'NNLO_to_Tune_Orig_Interp_Ratios_%s.coffea' % base_jobid))[args.year]
corrections = {
    'Pileup' : pu_correction,
    'Prefire' : True,
    'LeptonSF' : lepSF_correction,
    'BTagSF' : True,
    'JetCor' : jet_corrections,
    'Alpha' : alpha_corrections,
    'Kin_Rewt' : nnlo_reweighting,
}
##set_trace()
#if cfg_pars['NNLO']['Var'] is not None:
#    corrections.update({'Kin_Rewt' : nnlo_reweighting})
#    #corrections.update({
#    #    'Kin_Rewt' : {'Var' : cfg_pars['NNLO']['Var'], 'Correction' : nnlo_reweighting[cfg_pars['NNLO']['Var']]},
#    #    #'Kin_Rewt' : {'Var' : 'thad_pt', 'Correction' : nnlo_reweighting},
#    #    #'Kin_Rewt' : {'Var' : 'mtt_vs_thad_ctstar', 'Correction' : nnlo_reweighting},
#    #    #'Kin_Rewt' : {'Var' : 'mtt_vs_thad_ctstar_Interp', 'Correction' : nnlo_reweighting},
#    #})

    ## parameters for b-tagging
jet_pars = cfg_pars['Jets']
btaggers = ['DeepCSV']

wps_to_use = list(set([jet_pars['permutations']['tightb'],jet_pars['permutations']['looseb']]))
if not( len(wps_to_use) == 1):
    raise IOError("Only 1 unique btag working point supported now")
btag_wps = ['DeepCSVMedium']

if corrections['BTagSF'] == True:
    sf_file = os.path.join(proj_dir, 'Corrections', base_jobid, jet_pars['btagging']['btagSF_file'])
    if not os.path.isfile(sf_file):
        raise IOError("BTag SF file %s doesn't exist" % sf_file)

    btag_sfs = load(sf_file)
    corrections.update({'BTag_Constructors' : {}})
    for btagger in btaggers:
        threeJets = btag_sfs[args.year][btagger]['3Jets'][wps_to_use[0]]
        fourPJets = btag_sfs[args.year][btagger]['4PJets'][wps_to_use[0]]
        corrections['BTag_Constructors'].update({btagger : {'3Jets' : threeJets, '4PJets' : fourPJets} })
#set_trace()

MTcut = jet_pars['MT']

# 0 == '' (no gen matching), 1 == 'right', 2 == 'matchable', 3 == 'unmatchable', 4 == 'sl_tau', 5 == 'other' (not semilep)
perm_cats = {
    0 : '',
    1 : 'right',
    2 : 'matchable',
    3 : 'unmatchable',
    4 : 'sl_tau',
    5 : 'other',
}

    

# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class htt_nnlo_reweighting_study(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        #self.sys_axis = hist.Cat("sys", "Systematic")
        self.evt_type_axis = hist.Cat("evt_type", "Event Type")
        self.nnlo_reweighting_axis = hist.Cat("rewt", "NNLO Reweighting Type")
        self.jetmult_axis = hist.Cat("jmult", "nJets")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")
        self.lepcat_axis = hist.Cat("lepcat", "Lepton Category")
        self.btag_axis = hist.Cat("btag", "BTag Region")
        self.pt_axis = hist.Bin("pt", "p_{T} [GeV]", 200, 0, 1000)
        self.eta_axis = hist.Bin("eta", r"$\eta$", 200, -5., 5.)
        self.mtop_axis = hist.Bin("mtop", "m(top) [GeV]", 300, 0, 300)
        self.mtt_axis = hist.Bin("mtt", "m($t\overline{t}$) [GeV]", 360, 200, 2000)
        self.ctstar_axis = hist.Bin("ctstar", "cos($\\theta^{*}$)", 200, -1., 1.)
        #self.ctstar_axis = hist.Bin("ctstar", "cos($\\theta^{*}$)", 40, -1., 1.)
        self.ctstar_abs_axis = hist.Bin("ctstar_abs", "|cos($\\theta^{*}$)|", 200, 0., 1.)
        #self.ctstar_abs_axis = hist.Bin("ctstar_abs", "|cos($\\theta^{*}$)|", 20, 0., 1.)
        self.reso_pt_axis         = hist.Bin("reso_pt", "p_{T} [GeV]", 200, -200, 200)
        self.reso_eta_axis        = hist.Bin("reso_eta", r"$\eta$", 200, -5., 5.)
        self.reso_mtop_axis       = hist.Bin("reso_mtop", "m(top) [GeV]", 400, -200, 200)
        self.reso_mtt_axis        = hist.Bin("reso_mtt", "m($t\overline{t}$) [GeV]", 400, -200, 200)
        self.reso_ctstar_axis     = hist.Bin("reso_ctstar", "cos($\\theta^{*}$)", 400, -2., 2.)
        self.reso_ctstar_abs_axis = hist.Bin("reso_ctstar_abs", "|cos($\\theta^{*}$)|", 200, -1., 1.)
        self.nnlo_wts_axis = hist.Bin("nnlo_wts", "NNLO Weight", 200, 0., 2.)

            ## make dictionary of hists
        histo_dict = {}
                ## make semilep plots
        semilep_hists = self.make_semilep_hists()
        histo_dict.update(semilep_hists)        
                ## make dilep/had plots
        dl_had_hists = self.make_dl_had_hists()
        histo_dict.update(dl_had_hists)        

        self.sample_name = ''
        self.corrections = corrections

            ## make dict of cutflow for each systematic variation
        histo_dict['cutflow_nosys'] = processor.defaultdict_accumulator(int)
    
        self._accumulator = processor.dict_accumulator(histo_dict)


    @property
    def accumulator(self):
        return self._accumulator

    def make_semilep_hists(self):
        histo_dict = {}
            # gen/reco hists
                ## plots for cross checking
        histo_dict['Nominal_mtt']        = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mtt_axis)
        histo_dict['mtt_vs_nnlo_weights']= hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mtt_axis, self.nnlo_wts_axis)

        ##histo_dict['mtt']      = hist.Hist("Events", self.dataset_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mtt_axis)
        histo_dict['mtt']      = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mtt_axis)
        histo_dict['mthad']    = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mtop_axis)
        histo_dict['pt_thad']  = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.pt_axis)
        histo_dict['pt_tlep']  = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.pt_axis)
        histo_dict['pt_tt']    = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.pt_axis)
        histo_dict['eta_thad'] = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.eta_axis)
        histo_dict['eta_tlep'] = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.eta_axis)
        histo_dict['eta_tt']   = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.eta_axis)

        histo_dict['tlep_ctstar']     = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.ctstar_axis)
        histo_dict['tlep_ctstar_abs'] = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.ctstar_abs_axis)

        histo_dict['mtt_vs_tlep_ctstar_abs'] = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mtt_axis, self.ctstar_abs_axis)

            # gen-reco resolution hists
        histo_dict['Reso_mtt']      = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_mtt_axis)
        histo_dict['Reso_mthad']    = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_mtop_axis)
        histo_dict['Reso_pt_thad']  = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_pt_axis)
        histo_dict['Reso_pt_tlep']  = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_pt_axis)
        histo_dict['Reso_pt_tt']    = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_pt_axis)
        histo_dict['Reso_eta_thad'] = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_eta_axis)
        histo_dict['Reso_eta_tlep'] = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_eta_axis)
        histo_dict['Reso_eta_tt']   = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_eta_axis)

        histo_dict['Reso_tlep_ctstar']     = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_ctstar_axis)
        histo_dict['Reso_tlep_ctstar_abs'] = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_ctstar_abs_axis)

        histo_dict['Reso_mtt_vs_tlep_ctstar_abs'] = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_mtt_axis, self.reso_ctstar_abs_axis)

        return histo_dict

    def make_dl_had_hists(self):
        histo_dict = {}
            # gen/reco hists
        histo_dict['mtop']    = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mtop_axis)
        histo_dict['mtbar']   = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mtop_axis)
        histo_dict['pt_top']  = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.pt_axis)
        histo_dict['pt_tbar'] = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.pt_axis)
        histo_dict['eta_top'] = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.eta_axis)
        histo_dict['eta_tbar']= hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.eta_axis)

        histo_dict['top_ctstar']     = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.ctstar_axis)
        histo_dict['tbar_ctstar']    = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.ctstar_axis)
        histo_dict['top_ctstar_abs'] = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.ctstar_abs_axis)
        histo_dict['tbar_ctstar_abs']= hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.ctstar_abs_axis)

        histo_dict['mtt_vs_top_ctstar_abs'] = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mtt_axis, self.ctstar_abs_axis)
        histo_dict['mtt_vs_tbar_ctstar_abs']= hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mtt_axis, self.ctstar_abs_axis)

            # gen-reco resolution hists
        histo_dict['Reso_mtop']    = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_mtop_axis)
        histo_dict['Reso_mtbar']   = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_mtop_axis)
        histo_dict['Reso_pt_top']  = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_pt_axis)
        histo_dict['Reso_pt_tbar'] = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_pt_axis)
        histo_dict['Reso_eta_top'] = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_eta_axis)
        histo_dict['Reso_eta_tbar']= hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_eta_axis)

        histo_dict['Reso_top_ctstar']     = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_ctstar_axis)
        histo_dict['Reso_tbar_ctstar']    = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_ctstar_axis)
        histo_dict['Reso_top_ctstar_abs'] = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_ctstar_abs_axis)
        histo_dict['Reso_tbar_ctstar_abs']= hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_ctstar_abs_axis)

        histo_dict['Reso_mtt_vs_top_ctstar_abs'] = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_mtt_axis, self.reso_ctstar_abs_axis)
        histo_dict['Reso_mtt_vs_tbar_ctstar_abs']= hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_mtt_axis, self.reso_ctstar_abs_axis)

        return histo_dict

    def process(self, df):
        np.random.seed(10) # sets seed so values from random distributions are reproducible (JER corrections)
        output = self.accumulator.identity()

        self.sample_name = df.dataset

                ## initialize regions
        regions = {
            'Muon' : {
                'Tight' : {
                    'btagPass' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'tight_MU', 'DeepCSV_pass'},
                        '4PJets' : {'objselection', 'jets_4p', 'tight_MU', 'DeepCSV_pass'},
                    },
                    'btagFail' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'tight_MU', 'DeepCSV_fail'},
                        '4PJets' : {'objselection', 'jets_4p', 'tight_MU', 'DeepCSV_fail'},
                    },
                },
                'Loose' : {
                    'btagPass' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'loose_MU', 'DeepCSV_pass'},
                        '4PJets' : {'objselection', 'jets_4p', 'loose_MU', 'DeepCSV_pass'},
                    },
                    'btagFail' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'loose_MU', 'DeepCSV_fail'},
                        '4PJets' : {'objselection', 'jets_4p', 'loose_MU', 'DeepCSV_fail'},
                    },
                },
            },
            'Electron' : {
                'Tight' : {
                    'btagPass' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'tight_EL', 'DeepCSV_pass'},
                        '4PJets' : {'objselection', 'jets_4p', 'tight_EL', 'DeepCSV_pass'},
                    },
                    'btagFail' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'tight_EL', 'DeepCSV_fail'},
                        '4PJets' : {'objselection', 'jets_4p', 'tight_EL', 'DeepCSV_fail'},
                    },
                },
                'Loose' : {
                    'btagPass' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'loose_EL', 'DeepCSV_pass'},
                        '4PJets' : {'objselection', 'jets_4p', 'loose_EL', 'DeepCSV_pass'},
                    },
                    'btagFail' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'loose_EL', 'DeepCSV_fail'},
                        '4PJets' : {'objselection', 'jets_4p', 'loose_EL', 'DeepCSV_fail'},
                    },
                },
            },
        }

            ## make event weights
                # data or MC distinction made internally
        mu_evt_weights = MCWeights.get_event_weights(df, year=args.year, corrections=self.corrections, BTagSFs=btaggers)#, isTTbar=isNominal_ttJets_)
        el_evt_weights = MCWeights.get_event_weights(df, year=args.year, corrections=self.corrections, BTagSFs=btaggers)#, isTTbar=isNominal_ttJets_)

            ## initialize selections
        selection = processor.PackedSelection()

            ## object selection
        objsel_evts = objsel.select(df, year=args.year, corrections=self.corrections, cutflow=output['cutflow_nosys'])
        output['cutflow_nosys']['nEvts passing jet and lepton obj selection'] += objsel_evts.sum()
        selection.add('jets_3', df['Jet'].counts == 3)
        selection.add('jets_4p', df['Jet'].counts > 3)
        selection.add('objselection', objsel_evts)
        #selection.add('DeepJet_pass', df['Jet']['DeepJet'+wps_to_use[0]].sum() >= 2)            
        selection.add('DeepCSV_pass', df['Jet']['DeepCSV'+wps_to_use[0]].sum() >= 2)
        #selection.add('DeepCSV_fail', (df['Jet']['btagDeepB'].max() < 0.2) & (df['Jet']['btagDeepB'].max() > 0.1) ) # max DeepCSV value is between 0.1 and 0.2

            # sort jets by btag value
        df['Jet'] = df['Jet'][df['Jet']['btagDeepB'].argsort(ascending=False)] if btaggers[0] == 'DeepCSV' else df['Jet'][df['Jet']['btagDeepFlavB'].argsort(ascending=False)]

            # btag fail sideband
        deepcsv_sorted = df['Jet'][df['Jet']['btagDeepB'].argsort(ascending=False)]['btagDeepB']
        valid_counts_inds = np.where(df['Jet'].counts > 1)[0]
        deepcsv_fail = np.zeros(df.size).astype(bool)
        deepcsv_fail[valid_counts_inds] = (deepcsv_sorted[valid_counts_inds][:, 0] < btag_values[args.year]['btagDeepB']['DeepCSV'+wps_to_use[0]]) & (deepcsv_sorted[valid_counts_inds][:, 1] < btag_values[args.year]['btagDeepB']['DeepCSV'+wps_to_use[0]])
        selection.add('DeepCSV_fail', deepcsv_fail ) # highest and second highest DeepCSV values don't pass tight and loose WPs

            ## add different selections
                ## muons
        selection.add('tight_MU', df['Muon']['TIGHTMU'].sum() == 1) # one muon passing TIGHT criteria
        selection.add('loose_MU', df['Muon']['LOOSEMU'].sum() == 1) # one muon passing LOOSE criteria
                ## electrons
        selection.add('tight_EL', df['Electron']['TIGHTEL'].sum() == 1) # one electron passing TIGHT criteria
        selection.add('loose_EL', df['Electron']['LOOSEEL'].sum() == 1) # one electron passing LOOSE criteria

        #set_trace()
        ### apply lepton SFs to MC (only applicable to tight leptons)
        if 'LeptonSF' in corrections.keys():
            tight_mu_cut = selection.require(objselection=True, tight_MU=True) # find events passing muon object selection with one tight muon
            tight_muons = df['Muon'][tight_mu_cut][(df['Muon'][tight_mu_cut]['TIGHTMU'] == True)]
            muSFs_dict =  MCWeights.get_lepton_sf(year=args.year, lepton='Muons', corrections=lepSF_correction,
                pt=tight_muons.pt.flatten(), eta=tight_muons.eta.flatten())
            mu_reco_cen = np.ones(df.size)
            mu_reco_err = np.zeros(df.size)
            mu_trig_cen = np.ones(df.size)
            mu_trig_err = np.zeros(df.size)
            mu_reco_cen[tight_mu_cut] = muSFs_dict['RECO_CEN']
            mu_reco_err[tight_mu_cut] = muSFs_dict['RECO_ERR']
            mu_trig_cen[tight_mu_cut] = muSFs_dict['TRIG_CEN']
            mu_trig_err[tight_mu_cut] = muSFs_dict['TRIG_ERR']
            mu_evt_weights.add('Lep_RECO', mu_reco_cen, mu_reco_err, mu_reco_err, shift=True)
            mu_evt_weights.add('Lep_TRIG', mu_trig_cen, mu_trig_err, mu_trig_err, shift=True)

            tight_el_cut = selection.require(objselection=True, tight_EL=True) # find events passing electron object selection with one tight electron
            tight_electrons = df['Electron'][tight_el_cut][(df['Electron'][tight_el_cut]['TIGHTEL'] == True)]
            elSFs_dict = MCWeights.get_lepton_sf(year=args.year, lepton='Electrons', corrections=lepSF_correction,
                pt=tight_electrons.pt.flatten(), eta=tight_electrons.etaSC.flatten())
            el_reco_cen = np.ones(df.size)
            el_reco_err = np.zeros(df.size)
            el_trig_cen = np.ones(df.size)
            el_trig_err = np.zeros(df.size)
            el_reco_cen[tight_el_cut] = elSFs_dict['RECO_CEN']
            el_reco_err[tight_el_cut] = elSFs_dict['RECO_ERR']
            el_trig_cen[tight_el_cut] = elSFs_dict['TRIG_CEN']
            el_trig_err[tight_el_cut] = elSFs_dict['TRIG_ERR']
            el_evt_weights.add('Lep_RECO', el_reco_cen, el_reco_err, el_reco_err, shift=True)
            el_evt_weights.add('Lep_TRIG', el_trig_cen, el_trig_err, el_trig_err, shift=True)

            ## apply btagging SFs to MC
        if corrections['BTagSF'] == True:
            #set_trace()
            deepcsv_cen = np.ones(df.size)
            deepcsv_bc_up = np.ones(df.size)
            deepcsv_bc_dw = np.ones(df.size)
            deepcsv_l_up = np.ones(df.size)
            deepcsv_l_dw = np.ones(df.size)
            threeJets_cut = selection.require(objselection=True, jets_3=True)
            deepcsv_3j_wts = self.corrections['BTag_Constructors']['DeepCSV']['3Jets'].get_scale_factor(jets=df['Jet'][threeJets_cut], passing_cut='DeepCSV'+wps_to_use[0])
            deepcsv_cen[threeJets_cut] = deepcsv_3j_wts['central'].prod()
            deepcsv_bc_up[threeJets_cut] = deepcsv_3j_wts['bc_up'].prod()
            deepcsv_bc_dw[threeJets_cut] = deepcsv_3j_wts['bc_down'].prod()
            deepcsv_l_up[threeJets_cut] = deepcsv_3j_wts['udsg_up'].prod()
            deepcsv_l_dw[threeJets_cut] = deepcsv_3j_wts['udsg_down'].prod()

            fourplusJets_cut = selection.require(objselection=True, jets_4p=True)
            deepcsv_4pj_wts = self.corrections['BTag_Constructors']['DeepCSV']['4PJets'].get_scale_factor(jets=df['Jet'][fourplusJets_cut], passing_cut='DeepCSV'+wps_to_use[0])
            deepcsv_cen[fourplusJets_cut] = deepcsv_4pj_wts['central'].prod()
            deepcsv_bc_up[fourplusJets_cut] = deepcsv_4pj_wts['bc_up'].prod()
            deepcsv_bc_dw[fourplusJets_cut] = deepcsv_4pj_wts['bc_down'].prod()
            deepcsv_l_up[fourplusJets_cut] = deepcsv_4pj_wts['udsg_up'].prod()
            deepcsv_l_dw[fourplusJets_cut] = deepcsv_4pj_wts['udsg_down'].prod()

                # make dict of btag weights
            btag_weights = {
                'DeepCSV_CEN' : deepcsv_cen,
                'DeepCSV_bc_UP' : deepcsv_bc_up,
                'DeepCSV_bc_DW' : deepcsv_bc_dw,
                'DeepCSV_l_UP' : deepcsv_l_up,
                'DeepCSV_l_DW' : deepcsv_l_dw,
            }

        GenTTbar = genpsel.select(df, mode='NORMAL')
        if 'Kin_Rewt' in self.corrections.keys():
            ##kin_wts = MCWeights.get_kin_weights(self.corrections['Kin_Rewt'], GenTTbar)
            #mu_evt_weights.add('%s_reweighting' % corrections['Kin_Rewt']['Var'], kin_wts)
            #el_evt_weights.add('%s_reweighting' % corrections['Kin_Rewt']['Var'], kin_wts)
            #set_trace()
            nnlo_weights = {var : MCWeights.get_kin_weights({'Var': var, 'Correction' : self.corrections['Kin_Rewt'][var]}, GenTTbar) for var in self.corrections['Kin_Rewt'].keys()}
            nnlo_weights.update({'Nominal' : np.ones(df.size)})


        #set_trace()
        ## fill hists for each region
        for lepton in regions.keys():
            evt_weights = mu_evt_weights if lepton == 'Muon' else el_evt_weights
            for leptype in regions[lepton].keys():
                for btagregion in regions[lepton][leptype].keys():
                    for jmult in regions[lepton][leptype][btagregion].keys():
                        cut = selection.all(*regions[lepton][leptype][btagregion][jmult])
                        #set_trace()

                        output['cutflow_nosys']['nEvts %s' % ', '.join([lepton, leptype, btagregion, jmult])] += cut.sum()

                        if args.debug: print(lepton, leptype, btagregion, jmult)
                        if cut.sum() > 0:
                            ltype = 'MU' if lepton == 'Muon' else 'EL'
                            if 'loose_or_tight_%s' % ltype in regions[lepton][leptype][btagregion][jmult]:
                                lep_mask = ((df[lepton][cut]['TIGHT%s' % ltype] == True) | (df[lepton][cut]['LOOSE%s' % ltype] == True))
                            elif 'tight_%s' % ltype in regions[lepton][leptype][btagregion][jmult]:
                                lep_mask = (df[lepton][cut]['TIGHT%s' % ltype] == True)
                            elif 'loose_%s' % ltype in regions[lepton][leptype][btagregion][jmult]:
                                lep_mask = (df[lepton][cut]['LOOSE%s' % ltype] == True)
                            else:
                                raise ValueError("Not sure what lepton type to choose for event")


                                # find best permutations
                            best_perms = ttpermutator.find_best_permutations(jets=df['Jet'][cut], leptons=df[lepton][cut][lep_mask], MET=df['MET'][cut], btagWP=btag_wps[0], btag_req=False if btagregion == 'btagFail' else True)
                            valid_perms = best_perms['TTbar'].counts > 0
                            output['cutflow_nosys']['nEvts %s: valid perms' % ', '.join([lepton, leptype, btagregion, jmult])] += valid_perms.sum()

                            bp_status = np.zeros(cut.size, dtype=int) # 0 == '' (no gen matching), 1 == 'right', 2 == 'matchable', 3 == 'unmatchable', 4 == 'noslep'
                                # get matched permutation (semilep ttbar only)
                            if isTTbar_:
                                semilep_evts = GenTTbar['SL']['TTbar'].counts > 0
                                bp_status[~semilep_evts] = 5
                                #if (cut & semilep_evts).sum() == 0:
                                #    continue
                                if 'loose_or_tight_%s' % ltype in regions[lepton][leptype][btagregion][jmult]:
                                    sl_lep_mask = ((df[lepton][(cut & semilep_evts)]['TIGHT%s' % ltype] == True) | (df[lepton][(cut & semilep_evts)]['LOOSE%s' % ltype] == True))
                                elif 'tight_%s' % ltype in regions[lepton][leptype][btagregion][jmult]:
                                    sl_lep_mask = (df[lepton][(cut & semilep_evts)]['TIGHT%s' % ltype] == True)
                                elif 'loose_%s' % ltype in regions[lepton][leptype][btagregion][jmult]:
                                    sl_lep_mask = (df[lepton][(cut & semilep_evts)]['LOOSE%s' % ltype] == True)
                                else:
                                    raise ValueError("Not sure what lepton type to choose for event")
                                mp = ttmatcher.best_match(gen_hyp=GenTTbar[(cut & semilep_evts)], jets=df['Jet'][(cut & semilep_evts)], leptons=df[lepton][(cut & semilep_evts)][sl_lep_mask], met=df['MET'][(cut & semilep_evts)])
                                perm_cat_array = compare_matched_best_perms(mp, best_perms, njets=jmult, bp_mask=semilep_evts[cut])
                                bp_status[(cut & semilep_evts)] = perm_cat_array
                                    # consider tau events as separate
                                sl_tau_evts = (np.abs(GenTTbar['SL']['Lepton'].pdgId) == 15).pad(1).fillna(False).flatten() # fillna to make sure array is same size as df.size
                                bp_status[sl_tau_evts] = 4

                                ## calculate MT
                            MT = make_vars.MT(df[lepton][cut][lep_mask], df['MET'][cut])
                            MTHigh = (MT[valid_perms] >= MTcut).flatten()
                            output['cutflow_nosys']['nEvts %s: pass MT cut' % ', '.join([lepton, leptype, btagregion, jmult])] += MTHigh.sum()

                                # make final objects/values after applying all cuts
                            evt_wts = (evt_weights.weight()*btag_weights['%s_CEN' % btaggers[0]])[cut][valid_perms][MTHigh]
                            #final_evt_wts = (evt_weights.weight()*btag_weights['%s_CEN' % btaggers[0]])[cut][valid_perms][MTHigh]
                            final_gentt = GenTTbar[cut][valid_perms][MTHigh]
                            final_bps = best_perms[valid_perms][MTHigh]
                            final_bp_status = bp_status[cut][valid_perms][MTHigh]
                            
                            sl_evts = final_gentt['SL']['TTbar'].counts > 0
                            dl_evts = final_gentt['DL']['TTbar'].counts > 0
                            had_evts = final_gentt['Had']['TTbar'].counts > 0

                            #set_trace()
                            for rewt_type in nnlo_weights.keys():
                                #set_trace()
                                final_nnlo_wts = nnlo_weights[rewt_type][cut][valid_perms][MTHigh]
                                final_evt_wts = evt_wts*final_nnlo_wts
                                    # fill hists
                                if sl_evts.sum() > 0:
                                    output = self.fill_semilep_hists(acc=output, rewt_type=rewt_type, jetmult=jmult, leptype=lepton, lepcat=leptype, btagregion=btagregion, permarray=final_bp_status[sl_evts], perm=final_bps[sl_evts], gentt=final_gentt[sl_evts], evt_weights=final_evt_wts[sl_evts], nnlo_weights=final_nnlo_wts[sl_evts])
                                    #output = self.fill_semilep_hists(acc=output, jetmult=jmult, leptype=lepton, lepcat=leptype, btagregion=btagregion, permarray=final_bp_status[sl_evts], perm=final_bps[sl_evts], gentt=final_gentt[sl_evts], evt_weights=final_evt_wts[sl_evts])
                                if dl_evts.sum() > 0:
                                    output = self.fill_dilep_had_hists(acc=output, rewt_type=rewt_type, jetmult=jmult, leptype=lepton, lepcat=leptype, btagregion=btagregion, ttdecay='DL', permarray=final_bp_status[dl_evts], perm=final_bps[dl_evts], gentt=final_gentt[dl_evts], evt_weights=final_evt_wts[dl_evts], nnlo_weights=final_nnlo_wts[dl_evts])
                                    #output = self.fill_dilep_had_hists(acc=output, jetmult=jmult, leptype=lepton, lepcat=leptype, btagregion=btagregion, ttdecay='DL', permarray=final_bp_status[dl_evts], perm=final_bps[dl_evts], gentt=final_gentt[dl_evts], evt_weights=final_evt_wts[dl_evts])
                                if had_evts.sum() > 0:
                                    output = self.fill_dilep_had_hists(acc=output, rewt_type=rewt_type, jetmult=jmult, leptype=lepton, lepcat=leptype, btagregion=btagregion, ttdecay='Had', permarray=final_bp_status[had_evts], perm=final_bps[had_evts], gentt=final_gentt[had_evts], evt_weights=final_evt_wts[had_evts], nnlo_weights=final_nnlo_wts[had_evts])
                                    #output = self.fill_dilep_had_hists(acc=output, jetmult=jmult, leptype=lepton, lepcat=leptype, btagregion=btagregion, ttdecay='Had', permarray=final_bp_status[had_evts], perm=final_bps[had_evts], gentt=final_gentt[had_evts], evt_weights=final_evt_wts[had_evts])


        return output



    def fill_semilep_hists(self, acc, rewt_type, jetmult, leptype, lepcat, btagregion, permarray, perm, gentt, evt_weights, nnlo_weights):
    #def fill_semilep_hists(self, acc, jetmult, leptype, lepcat, btagregion, permarray, perm, gentt, evt_weights):
            ## apply alpha correction for 3Jets
        if jetmult == '3Jets':
            orig_thad_p4, reco_tlep_p4 = perm['THad'].p4.flatten(), perm['TLep'].p4.flatten()
            alpha_corr = self.corrections['Alpha'](172.5/orig_thad_p4.mass)
            reco_thad_p4 = orig_thad_p4*alpha_corr

        else:
            reco_thad_p4, reco_tlep_p4 = perm['THad'].p4.flatten(), perm['TLep'].p4.flatten()
            
        reco_ttbar_p4 = (reco_thad_p4 + reco_tlep_p4)
        reco_thad_ctstar, reco_tlep_ctstar = make_vars.ctstar_flat(reco_thad_p4, reco_tlep_p4)

        gen_ttbar_p4, gen_thad_p4, gen_tlep_p4 = gentt['SL']['TTbar'].p4.flatten(), gentt['SL']['THad'].p4.flatten(), gentt['SL']['TLep'].p4.flatten()
        gen_thad_ctstar, gen_tlep_ctstar = make_vars.ctstar_flat(gen_thad_p4, gen_tlep_p4)

        #set_trace()
        for permval in np.unique(permarray).tolist():
            #set_trace()
            perm_inds = np.where(permarray == permval)
            dataset_name = '%s_%s' % (self.sample_name, perm_cats[permval]) if permval != 0 else self.sample_name

                ## fill reco-level quantities
            #acc['mtt'].fill(     dataset=dataset_name, evt_type='RECO', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtt=reco_ttbar_p4.mass[perm_inds], weight=evt_weights[perm_inds])
            acc['mtt'].fill(     dataset=dataset_name, rewt=rewt_type, evt_type='RECO', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtt=reco_ttbar_p4.mass[perm_inds], weight=evt_weights[perm_inds])
            acc['mthad'].fill(   dataset=dataset_name, rewt=rewt_type, evt_type='RECO', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtop=reco_thad_p4.mass[perm_inds], weight=evt_weights[perm_inds])
            acc['pt_thad'].fill( dataset=dataset_name, rewt=rewt_type, evt_type='RECO', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=reco_thad_p4.pt[perm_inds], weight=evt_weights[perm_inds])
            acc['pt_tlep'].fill( dataset=dataset_name, rewt=rewt_type, evt_type='RECO', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=reco_tlep_p4.pt[perm_inds], weight=evt_weights[perm_inds])
            acc['pt_tt'].fill(   dataset=dataset_name, rewt=rewt_type, evt_type='RECO', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=reco_ttbar_p4.pt[perm_inds], weight=evt_weights[perm_inds])
            acc['eta_thad'].fill(dataset=dataset_name, rewt=rewt_type, evt_type='RECO', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=reco_thad_p4.eta[perm_inds], weight=evt_weights[perm_inds])
            acc['eta_tlep'].fill(dataset=dataset_name, rewt=rewt_type, evt_type='RECO', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=reco_tlep_p4.eta[perm_inds], weight=evt_weights[perm_inds])
            acc['eta_tt'].fill(  dataset=dataset_name, rewt=rewt_type, evt_type='RECO', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=reco_ttbar_p4.eta[perm_inds], weight=evt_weights[perm_inds])

            acc['tlep_ctstar'].fill(    dataset=dataset_name, rewt=rewt_type, evt_type='RECO', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, ctstar=reco_tlep_ctstar[perm_inds], weight=evt_weights[perm_inds])
            acc['tlep_ctstar_abs'].fill(dataset=dataset_name, rewt=rewt_type, evt_type='RECO', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, ctstar_abs=np.abs(reco_tlep_ctstar[perm_inds]), weight=evt_weights[perm_inds])

            acc['mtt_vs_tlep_ctstar_abs'].fill(dataset=dataset_name, rewt=rewt_type, evt_type='RECO', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtt=reco_ttbar_p4.mass[perm_inds], ctstar_abs=np.abs(reco_tlep_ctstar[perm_inds]), weight=evt_weights[perm_inds])
            #set_trace()
            acc['mtt_vs_nnlo_weights'].fill(dataset=dataset_name, rewt=rewt_type, evt_type='RECO', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtt=reco_ttbar_p4.mass[perm_inds], nnlo_wts=nnlo_weights[perm_inds])
            acc['Nominal_mtt'].fill(        dataset=dataset_name, rewt=rewt_type, evt_type='RECO', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtt=reco_ttbar_p4.mass[perm_inds], weight=(evt_weights/nnlo_weights)[perm_inds])

                ## fill gen-level quantities
            #acc['mtt'].fill(     dataset=dataset_name, evt_type='GEN', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtt=gen_ttbar_p4.mass[perm_inds], weight=evt_weights[perm_inds])
            acc['mtt'].fill(     dataset=dataset_name, rewt=rewt_type, evt_type='GEN', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtt=gen_ttbar_p4.mass[perm_inds], weight=evt_weights[perm_inds])
            acc['mthad'].fill(   dataset=dataset_name, rewt=rewt_type, evt_type='GEN', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtop=gen_thad_p4.mass[perm_inds], weight=evt_weights[perm_inds])
            acc['pt_thad'].fill( dataset=dataset_name, rewt=rewt_type, evt_type='GEN', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=gen_thad_p4.pt[perm_inds], weight=evt_weights[perm_inds])
            acc['pt_tlep'].fill( dataset=dataset_name, rewt=rewt_type, evt_type='GEN', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=gen_tlep_p4.pt[perm_inds], weight=evt_weights[perm_inds])
            acc['pt_tt'].fill(   dataset=dataset_name, rewt=rewt_type, evt_type='GEN', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=gen_ttbar_p4.pt[perm_inds], weight=evt_weights[perm_inds])
            acc['eta_thad'].fill(dataset=dataset_name, rewt=rewt_type, evt_type='GEN', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=gen_thad_p4.eta[perm_inds], weight=evt_weights[perm_inds])
            acc['eta_tlep'].fill(dataset=dataset_name, rewt=rewt_type, evt_type='GEN', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=gen_tlep_p4.eta[perm_inds], weight=evt_weights[perm_inds])
            acc['eta_tt'].fill(  dataset=dataset_name, rewt=rewt_type, evt_type='GEN', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=gen_ttbar_p4.eta[perm_inds], weight=evt_weights[perm_inds])

            acc['tlep_ctstar'].fill(    dataset=dataset_name, rewt=rewt_type, evt_type='GEN', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, ctstar=gen_tlep_ctstar[perm_inds], weight=evt_weights[perm_inds])
            acc['tlep_ctstar_abs'].fill(dataset=dataset_name, rewt=rewt_type, evt_type='GEN', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, ctstar_abs=np.abs(gen_tlep_ctstar[perm_inds]), weight=evt_weights[perm_inds])

            acc['mtt_vs_tlep_ctstar_abs'].fill(dataset=dataset_name, rewt=rewt_type, evt_type='GEN', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtt=gen_ttbar_p4.mass[perm_inds], ctstar_abs=np.abs(gen_tlep_ctstar[perm_inds]), weight=evt_weights[perm_inds])

            acc['mtt_vs_nnlo_weights'].fill(dataset=dataset_name, rewt=rewt_type, evt_type='GEN', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtt=gen_ttbar_p4.mass[perm_inds], nnlo_wts=nnlo_weights[perm_inds])
            acc['Nominal_mtt'].fill(        dataset=dataset_name, rewt=rewt_type, evt_type='GEN', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtt=gen_ttbar_p4.mass[perm_inds], weight=(evt_weights/nnlo_weights)[perm_inds])


                ## fill gen-reco quantities
            acc['Reso_mtt'].fill(     dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, reso_mtt=gen_ttbar_p4.mass[perm_inds]-reco_ttbar_p4.mass[perm_inds], weight=evt_weights[perm_inds])
            acc['Reso_mthad'].fill(   dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, reso_mtop=gen_thad_p4.mass[perm_inds]-reco_thad_p4.mass[perm_inds], weight=evt_weights[perm_inds])
            acc['Reso_pt_thad'].fill( dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, reso_pt=gen_thad_p4.pt[perm_inds]-reco_thad_p4.pt[perm_inds], weight=evt_weights[perm_inds])
            acc['Reso_pt_tlep'].fill( dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, reso_pt=gen_tlep_p4.pt[perm_inds]-reco_tlep_p4.pt[perm_inds], weight=evt_weights[perm_inds])
            acc['Reso_pt_tt'].fill(   dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, reso_pt=gen_ttbar_p4.pt[perm_inds]-reco_ttbar_p4.pt[perm_inds], weight=evt_weights[perm_inds])
            acc['Reso_eta_thad'].fill(dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, reso_eta=gen_thad_p4.eta[perm_inds]-reco_thad_p4.eta[perm_inds], weight=evt_weights[perm_inds])
            acc['Reso_eta_tlep'].fill(dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, reso_eta=gen_tlep_p4.eta[perm_inds]-reco_tlep_p4.eta[perm_inds], weight=evt_weights[perm_inds])
            acc['Reso_eta_tt'].fill(  dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, reso_eta=gen_ttbar_p4.eta[perm_inds]-reco_ttbar_p4.eta[perm_inds], weight=evt_weights[perm_inds])

            acc['Reso_tlep_ctstar'].fill(    dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, reso_ctstar=gen_tlep_ctstar[perm_inds]-reco_tlep_ctstar[perm_inds], weight=evt_weights[perm_inds])
            acc['Reso_tlep_ctstar_abs'].fill(dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, reso_ctstar_abs=np.abs(gen_tlep_ctstar[perm_inds])-np.abs(reco_tlep_ctstar[perm_inds]), weight=evt_weights[perm_inds])

            acc['Reso_mtt_vs_tlep_ctstar_abs'].fill(dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, reso_mtt=gen_ttbar_p4.mass[perm_inds]-reco_ttbar_p4.mass[perm_inds], reso_ctstar_abs=np.abs(gen_tlep_ctstar[perm_inds])-np.abs(reco_tlep_ctstar[perm_inds]), weight=evt_weights[perm_inds])

        return acc        

    def fill_dilep_had_hists(self, acc, rewt_type, jetmult, leptype, lepcat, btagregion, ttdecay, permarray, perm, gentt, evt_weights, nnlo_weights):
    #def fill_dilep_had_hists(self, acc, jetmult, leptype, lepcat, btagregion, ttdecay, permarray, perm, gentt, evt_weights):
            ## apply alpha correction for 3Jets
        if jetmult == '3Jets':
            orig_thad_p4, reco_tlep_p4 = perm['THad'].p4.flatten(), perm['TLep'].p4.flatten()
            alpha_corr = self.corrections['Alpha'](172.5/orig_thad_p4.mass)
            reco_thad_p4 = orig_thad_p4*alpha_corr

        else:
            reco_thad_p4, reco_tlep_p4 = perm['THad'].p4.flatten(), perm['TLep'].p4.flatten()

        #set_trace()
        reco_ttbar_p4 = (reco_thad_p4 + reco_tlep_p4)
        reco_thad_ctstar, reco_tlep_ctstar = make_vars.ctstar_flat(reco_thad_p4, reco_tlep_p4)
            # convert thad/tlep into top/tbar for gen comparison
        reco_top_mass = np.where(perm['Lepton'].charge.flatten() == 1, reco_tlep_p4.mass, reco_thad_p4.mass)            
        reco_tbar_mass = np.where(perm['Lepton'].charge.flatten() == -1, reco_tlep_p4.mass, reco_thad_p4.mass)            
        reco_top_pt = np.where(perm['Lepton'].charge.flatten() == 1, reco_tlep_p4.pt, reco_thad_p4.pt)            
        reco_tbar_pt = np.where(perm['Lepton'].charge.flatten() == -1, reco_tlep_p4.pt, reco_thad_p4.pt)            
        reco_top_eta = np.where(perm['Lepton'].charge.flatten() == 1, reco_tlep_p4.eta, reco_thad_p4.eta)            
        reco_tbar_eta = np.where(perm['Lepton'].charge.flatten() == -1, reco_tlep_p4.eta, reco_thad_p4.eta)            
        reco_top_ctstar = np.where(perm['Lepton'].charge.flatten() == 1, reco_tlep_ctstar, reco_thad_ctstar)            
        reco_tbar_ctstar = np.where(perm['Lepton'].charge.flatten() == -1, reco_tlep_ctstar, reco_thad_ctstar)            

        #set_trace()
        gen_ttbar_p4, gen_top_p4, gen_tbar_p4 = gentt[ttdecay]['TTbar'].p4.flatten(), gentt[ttdecay]['Top'].p4.flatten(), gentt[ttdecay]['Tbar'].p4.flatten()
        gen_top_ctstar, gen_tbar_ctstar = make_vars.ctstar_flat(gen_top_p4, gen_tbar_p4)

        for permval in np.unique(permarray).tolist():
            #set_trace()
            perm_inds = np.where(permarray == permval)
            dataset_name = '%s_%s' % (self.sample_name, perm_cats[permval]) if permval != 0 else self.sample_name

                ## fill reco-level quantities
            #acc['mtt'].fill(     dataset=dataset_name, evt_type='RECO', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtt=reco_ttbar_p4.mass[perm_inds], weight=evt_weights[perm_inds])
            acc['mtt'].fill(     dataset=dataset_name, rewt=rewt_type, evt_type='RECO', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtt=reco_ttbar_p4.mass[perm_inds], weight=evt_weights[perm_inds])
            acc['mtop'].fill(    dataset=dataset_name, rewt=rewt_type, evt_type='RECO', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtop=reco_top_mass[perm_inds], weight=evt_weights[perm_inds])
            acc['mtbar'].fill(   dataset=dataset_name, rewt=rewt_type, evt_type='RECO', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtop=reco_tbar_mass[perm_inds], weight=evt_weights[perm_inds])
            acc['pt_top'].fill(  dataset=dataset_name, rewt=rewt_type, evt_type='RECO', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=reco_top_pt[perm_inds], weight=evt_weights[perm_inds])
            acc['pt_tbar'].fill( dataset=dataset_name, rewt=rewt_type, evt_type='RECO', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=reco_tbar_pt[perm_inds], weight=evt_weights[perm_inds])
            acc['pt_tt'].fill(   dataset=dataset_name, rewt=rewt_type, evt_type='RECO', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=reco_ttbar_p4.pt[perm_inds], weight=evt_weights[perm_inds])
            acc['eta_top'].fill( dataset=dataset_name, rewt=rewt_type, evt_type='RECO', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=reco_top_eta[perm_inds], weight=evt_weights[perm_inds])
            acc['eta_tbar'].fill(dataset=dataset_name, rewt=rewt_type, evt_type='RECO', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=reco_tbar_eta[perm_inds], weight=evt_weights[perm_inds])
            acc['eta_tt'].fill(  dataset=dataset_name, rewt=rewt_type, evt_type='RECO', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=reco_ttbar_p4.eta[perm_inds], weight=evt_weights[perm_inds])

            acc['top_ctstar'].fill(     dataset=dataset_name, rewt=rewt_type, evt_type='RECO',  jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, ctstar=reco_top_ctstar[perm_inds], weight=evt_weights[perm_inds])
            acc['tbar_ctstar'].fill(    dataset=dataset_name, rewt=rewt_type, evt_type='RECO',  jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, ctstar=reco_tbar_ctstar[perm_inds], weight=evt_weights[perm_inds])
            acc['top_ctstar_abs'].fill( dataset=dataset_name, rewt=rewt_type, evt_type='RECO',  jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, ctstar_abs=np.abs(reco_top_ctstar[perm_inds]), weight=evt_weights[perm_inds])
            acc['tbar_ctstar_abs'].fill(dataset=dataset_name, rewt=rewt_type, evt_type='RECO',  jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, ctstar_abs=np.abs(reco_tbar_ctstar[perm_inds]), weight=evt_weights[perm_inds])

            acc['mtt_vs_top_ctstar_abs'].fill( dataset=dataset_name, rewt=rewt_type, evt_type='RECO',  jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtt=reco_ttbar_p4.mass[perm_inds], ctstar_abs=np.abs(reco_top_ctstar[perm_inds]), weight=evt_weights[perm_inds])
            acc['mtt_vs_tbar_ctstar_abs'].fill(dataset=dataset_name, rewt=rewt_type, evt_type='RECO',  jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtt=reco_ttbar_p4.mass[perm_inds], ctstar_abs=np.abs(reco_tbar_ctstar[perm_inds]), weight=evt_weights[perm_inds])

            acc['mtt_vs_nnlo_weights'].fill(dataset=dataset_name, rewt=rewt_type, evt_type='RECO', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtt=reco_ttbar_p4.mass[perm_inds], nnlo_wts=nnlo_weights[perm_inds])
            acc['Nominal_mtt'].fill(        dataset=dataset_name, rewt=rewt_type, evt_type='RECO', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtt=reco_ttbar_p4.mass[perm_inds], weight=(evt_weights/nnlo_weights)[perm_inds])

                ## fill gen-level quantities
            acc['mtt'].fill(     dataset=dataset_name, rewt=rewt_type, evt_type='GEN', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtt=gen_ttbar_p4.mass[perm_inds], weight=evt_weights[perm_inds])
            acc['mtop'].fill(    dataset=dataset_name, rewt=rewt_type, evt_type='GEN', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtop=gen_top_p4.mass[perm_inds], weight=evt_weights[perm_inds])
            acc['mtbar'].fill(   dataset=dataset_name, rewt=rewt_type, evt_type='GEN', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtop=gen_tbar_p4.mass[perm_inds], weight=evt_weights[perm_inds])
            acc['pt_top'].fill(  dataset=dataset_name, rewt=rewt_type, evt_type='GEN', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=gen_top_p4.pt[perm_inds], weight=evt_weights[perm_inds])
            acc['pt_tbar'].fill( dataset=dataset_name, rewt=rewt_type, evt_type='GEN', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=gen_tbar_p4.pt[perm_inds], weight=evt_weights[perm_inds])
            acc['pt_tt'].fill(   dataset=dataset_name, rewt=rewt_type, evt_type='GEN', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=gen_ttbar_p4.pt[perm_inds], weight=evt_weights[perm_inds])
            acc['eta_top'].fill( dataset=dataset_name, rewt=rewt_type, evt_type='GEN', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=gen_top_p4.eta[perm_inds], weight=evt_weights[perm_inds])
            acc['eta_tbar'].fill(dataset=dataset_name, rewt=rewt_type, evt_type='GEN', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=gen_tbar_p4.eta[perm_inds], weight=evt_weights[perm_inds])
            acc['eta_tt'].fill(  dataset=dataset_name, rewt=rewt_type, evt_type='GEN', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=gen_ttbar_p4.eta[perm_inds], weight=evt_weights[perm_inds])

            acc['top_ctstar'].fill(     dataset=dataset_name, rewt=rewt_type, evt_type='GEN', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, ctstar=gen_top_ctstar[perm_inds], weight=evt_weights[perm_inds])
            acc['tbar_ctstar'].fill(    dataset=dataset_name, rewt=rewt_type, evt_type='GEN', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, ctstar=gen_tbar_ctstar[perm_inds], weight=evt_weights[perm_inds])
            acc['top_ctstar_abs'].fill( dataset=dataset_name, rewt=rewt_type, evt_type='GEN', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, ctstar_abs=np.abs(gen_top_ctstar[perm_inds]), weight=evt_weights[perm_inds])
            acc['tbar_ctstar_abs'].fill(dataset=dataset_name, rewt=rewt_type, evt_type='GEN', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, ctstar_abs=np.abs(gen_tbar_ctstar[perm_inds]), weight=evt_weights[perm_inds])

            acc['mtt_vs_top_ctstar_abs'].fill( dataset=dataset_name, rewt=rewt_type, evt_type='GEN', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtt=gen_ttbar_p4.mass[perm_inds], ctstar_abs=np.abs(gen_top_ctstar[perm_inds]), weight=evt_weights[perm_inds])
            acc['mtt_vs_tbar_ctstar_abs'].fill(dataset=dataset_name, rewt=rewt_type, evt_type='GEN', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtt=gen_ttbar_p4.mass[perm_inds], ctstar_abs=np.abs(gen_tbar_ctstar[perm_inds]), weight=evt_weights[perm_inds])

            acc['mtt_vs_nnlo_weights'].fill(dataset=dataset_name, rewt=rewt_type, evt_type='GEN', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtt=gen_ttbar_p4.mass[perm_inds], nnlo_wts=nnlo_weights[perm_inds])
            acc['Nominal_mtt'].fill(        dataset=dataset_name, rewt=rewt_type, evt_type='GEN', jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtt=gen_ttbar_p4.mass[perm_inds], weight=(evt_weights/nnlo_weights)[perm_inds])

                ## fill gen-reco quantities
            acc['Reso_mtt'].fill(     dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, reso_mtt=gen_ttbar_p4.mass[perm_inds]-reco_ttbar_p4.mass[perm_inds], weight=evt_weights[perm_inds])
            acc['Reso_mtop'].fill(    dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, reso_mtop=gen_top_p4.mass[perm_inds]-reco_top_mass[perm_inds], weight=evt_weights[perm_inds])
            acc['Reso_mtbar'].fill(   dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, reso_mtop=gen_tbar_p4.mass[perm_inds]-reco_tbar_mass[perm_inds], weight=evt_weights[perm_inds])
            acc['Reso_pt_top'].fill(  dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, reso_pt=gen_top_p4.pt[perm_inds]-reco_top_pt[perm_inds], weight=evt_weights[perm_inds])
            acc['Reso_pt_tbar'].fill( dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, reso_pt=gen_tbar_p4.pt[perm_inds]-reco_tbar_pt[perm_inds], weight=evt_weights[perm_inds])
            acc['Reso_pt_tt'].fill(   dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, reso_pt=gen_ttbar_p4.pt[perm_inds]-reco_ttbar_p4.pt[perm_inds], weight=evt_weights[perm_inds])
            acc['Reso_eta_top'].fill( dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, reso_eta=gen_top_p4.eta[perm_inds]-reco_top_eta[perm_inds], weight=evt_weights[perm_inds])
            acc['Reso_eta_tbar'].fill(dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, reso_eta=gen_tbar_p4.eta[perm_inds]-reco_tbar_eta[perm_inds], weight=evt_weights[perm_inds])
            acc['Reso_eta_tt'].fill(  dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, reso_eta=gen_ttbar_p4.eta[perm_inds]-reco_ttbar_p4.eta[perm_inds], weight=evt_weights[perm_inds])

            acc['Reso_top_ctstar'].fill(     dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, reso_ctstar=gen_top_ctstar[perm_inds]-reco_top_ctstar[perm_inds], weight=evt_weights[perm_inds])
            acc['Reso_tbar_ctstar'].fill(    dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, reso_ctstar=gen_tbar_ctstar[perm_inds]-reco_tbar_ctstar[perm_inds], weight=evt_weights[perm_inds])
            acc['Reso_top_ctstar_abs'].fill( dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, reso_ctstar_abs=np.abs(gen_top_ctstar[perm_inds])-np.abs(reco_top_ctstar[perm_inds]), weight=evt_weights[perm_inds])
            acc['Reso_tbar_ctstar_abs'].fill(dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, reso_ctstar_abs=np.abs(gen_tbar_ctstar[perm_inds])-np.abs(reco_tbar_ctstar[perm_inds]), weight=evt_weights[perm_inds])

            acc['Reso_mtt_vs_top_ctstar_abs'].fill( dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, reso_mtt=gen_ttbar_p4.mass[perm_inds]-reco_ttbar_p4.mass[perm_inds], reso_ctstar_abs=np.abs(gen_top_ctstar[perm_inds])-np.abs(reco_top_ctstar[perm_inds]), weight=evt_weights[perm_inds])
            acc['Reso_mtt_vs_tbar_ctstar_abs'].fill(dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, reso_mtt=gen_ttbar_p4.mass[perm_inds]-reco_ttbar_p4.mass[perm_inds], reso_ctstar_abs=np.abs(gen_tbar_ctstar[perm_inds])-np.abs(reco_tbar_ctstar[perm_inds]), weight=evt_weights[perm_inds])

        return acc        




    def postprocess(self, accumulator):
        return accumulator

proc_executor = processor.iterative_executor if args.debug else processor.futures_executor

output = processor.run_uproot_job(fileset,
    treename='Events',
    processor_instance=htt_nnlo_reweighting_study(),
    executor=proc_executor,
    executor_args={
        'workers': 8,
        'flatten' : True,
        'compression': 8,
        #'compression': 5,
    },
    chunksize=10000 if args.debug else 100000,
    #chunksize=10000 if args.debug else 50000,
    #chunksize=50000,
)


if args.debug:
    print(output)
#set_trace()
#print(output['cutflow'])

save(output, args.outfname)
print('%s has been written' % args.outfname)
