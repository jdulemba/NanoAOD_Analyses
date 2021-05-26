#!/usr/bin/env python

from coffea import hist, processor
from coffea.nanoevents import NanoAODSchema
from coffea.analysis_tools import PackedSelection
import awkward as ak
from coffea.nanoevents.methods import vector
ak.behavior.update(vector.behavior)

from coffea.util import save, load
from pdb import set_trace
import os#, sys
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
import python.IDJet as IDJet
from copy import deepcopy

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
base_jobid = os.environ['base_jobid']
analyzer = 'evtWeights'

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('fset', type=str, help='Fileset dictionary (in string form) to be used for the processor')
parser.add_argument('year', choices=['2017', '2018'] if base_jobid == 'ULnanoAOD' else ['2016', '2017', '2018'], help='Specify which year to run over')
parser.add_argument('outfname', type=str, help='Specify output filename, including directory and file extension')
parser.add_argument('--debug', action='store_true', help='Uses iterative_executor for debugging purposes, otherwise futures_excutor will be used (faster)')
parser.add_argument('--evt_sys', type=str, help='Specify event systematics to run. Default is (NONE,NONE) and all opts are capitalized through run_analyzer')
parser.add_argument('--rewt_sys', type=str, help='Specify reweighting systematics to run. Default is (NONE,NONE) and all opts are capitalized through run_analyzer')
parser.add_argument('--only_sys', type=int, help='Only run specified systematics and not nominal weights (nosys)')

args = parser.parse_args()

# convert input string of fileset dictionary to actual dictionary
fdict = (args.fset).replace("\'", "\"")
fileset = prettyjson.loads(fdict)

if len(fileset.keys()) > 1:
    raise ValueError("Only one topology run at a time in order to determine which corrections and systematics to run")
samplename = list(fileset.keys())[0]

# get dataset classification, used for corrections/systematics
isTTbar_ = samplename.startswith('ttJets')
isTTSL_ = samplename.startswith('ttJetsSL')
isSignal_ = (samplename.startswith('AtoTT') or samplename.startswith('HtoTT'))
isInt_ = isSignal_ and ('Int' in samplename)
isData_ = samplename.startswith('data')
if isData_:
    raise ValueError("This should only be run on MC")

Nominal_ttJets = ['ttJets_PS'] if ((args.year == '2016') and (base_jobid == 'NanoAODv6')) else ['ttJetsSL', 'ttJetsHad', 'ttJetsDiLep']
isNominal_ttJets_ = samplename in Nominal_ttJets
isTTShift_ = isTTbar_ and (samplename.endswith('UP') or samplename.endswith('DOWN') or 'mtop' in samplename)

## init tt probs for likelihoods
ttpermutator.year_to_run(year=args.year)

## load lumimask for data and corrections for event weights
pu_correction = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'MC_PU_Weights.coffea'))[args.year]
lepSF_correction = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'leptonSFs.coffea'))[args.year]
jet_corrections = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'JetMETCorrections.coffea'))[args.year]
#nnlo_reweighting = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'NNLO_to_Tune_Orig_Interp_Ratios_%s.coffea' % base_jobid))[args.year]
corrections = {
    'Pileup' : pu_correction,
    'Prefire' : True,
    'LeptonSF' : lepSF_correction,
    'BTagSF' : True,
    'JetCor' : jet_corrections,
    #'Kin_Rewt' : {'Var' : 'mtt_vs_thad_ctstar_Interp', 'Correction' : nnlo_reweighting},
}

cfg_pars = prettyjson.loads(open(os.path.join(proj_dir, 'cfg_files', 'cfg_pars_%s.json' % jobid)).read())
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

# get systematics to run
if isTTShift_: # data or separate systematics dataset
    event_systematics_to_run = ['nosys']
    reweight_systematics_to_run = ['nosys']
else: # MC samples
    import fnmatch
    if bool(args.only_sys): # don't run 'nosys'
        if (args.evt_sys == 'NONE') and (args.rewt_sys == 'NONE'):
            raise ValueError("At least one systematic must be specified in order to run only on systematics!")
    
        event_systematics_to_run = ['nosys'] if (args.rewt_sys != 'NONE') else []
        reweight_systematics_to_run = []
    
    else:
        event_systematics_to_run = ['nosys']
        reweight_systematics_to_run = ['nosys']

    event_systematics_to_run += [systematics.event_sys_opts[args.year][name] for name in systematics.event_sys_opts[args.year].keys() if fnmatch.fnmatch(name, args.evt_sys)]
    if isSignal_:
        reweight_systematics_to_run += [systematics.signal_reweight_opts[args.year][name] for name in systematics.signal_reweight_opts[args.year].keys() if fnmatch.fnmatch(name, args.rewt_sys)]
    else:
        reweight_systematics_to_run += [systematics.reweight_sys_opts[args.year][name] for name in systematics.reweight_sys_opts[args.year].keys() if fnmatch.fnmatch(name, args.rewt_sys)]

        ## check that systematics only related to ttbar events aren't used for non-ttbar events
    if not isTTbar_:
        reweight_systematics_to_run = [sys for sys in reweight_systematics_to_run if sys not in systematics.ttJets_sys.values()]
    
print('Running with event systematics:', *sorted(set(event_systematics_to_run).difference(set(['nosys']))), sep=', ') if 'nosys' in event_systematics_to_run else print('Running with event systematics:', *sorted(event_systematics_to_run), sep=', ')
print('  and reweight systematics:', *sorted(reweight_systematics_to_run), sep=', ')
#set_trace()


# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class evtWeights(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.sys_axis = hist.Cat("sys", "Systematic")
        self.jetmult_axis = hist.Cat("jmult", "nJets")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")
        self.lepcat_axis = hist.Cat("lepcat", "Lepton Category")
        self.btag_axis = hist.Cat("btag", "BTag Region")
        self.genwt_axis = hist.Bin("genwt", "Gen Weight", 1000, -500, 500)
        self.totwt_axis = hist.Bin("totwt", "Total Weight", 1000, -500, 500)

            ## make dictionary of hists
        histo_dict = {}
                ## make hists
        evtWeight_hists = self.make_hists()
        histo_dict.update(evtWeight_hists)

        self.sample_name = ''
        self.corrections = corrections
        self.reweight_systematics_to_run = reweight_systematics_to_run
        self.event_systematics_to_run = event_systematics_to_run

            ## make dict of cutflow for each systematic variation
        for sys in self.event_systematics_to_run:
            histo_dict['cutflow_%s' % sys] = processor.defaultdict_accumulator(int)
    
        self._accumulator = processor.dict_accumulator(histo_dict)

            ## define regions depending on data or MC, systematic or nosys
        base_regions = {
            'Muon' : {
                'Tight' : {
                    'btagPass' : {
                        '3Jets'  : {'lep_and_filter_pass', 'passing_jets', 'jets_3' , 'tight_MU', 'DeepCSV_pass'},
                        '4PJets' : {'lep_and_filter_pass', 'passing_jets', 'jets_4p', 'tight_MU', 'DeepCSV_pass'},
                    },
                    'btagFail' : {
                        '3Jets'  : {'lep_and_filter_pass', 'passing_jets', 'jets_3' , 'tight_MU', 'DeepCSV_fail'},
                        '4PJets' : {'lep_and_filter_pass', 'passing_jets', 'jets_4p', 'tight_MU', 'DeepCSV_fail'},
                    },
                },
                'Loose' : {
                    'btagPass' : {
                        '3Jets'  : {'lep_and_filter_pass', 'passing_jets', 'jets_3' , 'loose_MU', 'DeepCSV_pass'},
                        '4PJets' : {'lep_and_filter_pass', 'passing_jets', 'jets_4p', 'loose_MU', 'DeepCSV_pass'},
                    },
                    'btagFail' : {
                        '3Jets'  : {'lep_and_filter_pass', 'passing_jets', 'jets_3' , 'loose_MU', 'DeepCSV_fail'},
                        '4PJets' : {'lep_and_filter_pass', 'passing_jets', 'jets_4p', 'loose_MU', 'DeepCSV_fail'},
                    },
                },
            },
            'Electron' : {
                'Tight' : {
                    'btagPass' : {
                        '3Jets'  : {'lep_and_filter_pass', 'passing_jets', 'jets_3' , 'tight_EL', 'DeepCSV_pass'},
                        '4PJets' : {'lep_and_filter_pass', 'passing_jets', 'jets_4p', 'tight_EL', 'DeepCSV_pass'},
                    },
                    'btagFail' : {
                        '3Jets'  : {'lep_and_filter_pass', 'passing_jets', 'jets_3' , 'tight_EL', 'DeepCSV_fail'},
                        '4PJets' : {'lep_and_filter_pass', 'passing_jets', 'jets_4p', 'tight_EL', 'DeepCSV_fail'},
                    },
                },
                'Loose' : {
                    'btagPass' : {
                        '3Jets'  : {'lep_and_filter_pass', 'passing_jets', 'jets_3' , 'loose_EL', 'DeepCSV_pass'},
                        '4PJets' : {'lep_and_filter_pass', 'passing_jets', 'jets_4p', 'loose_EL', 'DeepCSV_pass'},
                    },
                    'btagFail' : {
                        '3Jets'  : {'lep_and_filter_pass', 'passing_jets', 'jets_3' , 'loose_EL', 'DeepCSV_fail'},
                        '4PJets' : {'lep_and_filter_pass', 'passing_jets', 'jets_4p', 'loose_EL', 'DeepCSV_fail'},
                    },
                },
            },
        }
        self.regions = {sys:deepcopy(base_regions) for sys in self.event_systematics_to_run}
            # remove sideband regions for systematics
        for sys in self.regions.keys():
            if (sys != 'nosys') or isTTShift_ or isSignal_:
                del self.regions[sys]['Muon']['Loose']
                del self.regions[sys]['Electron']['Loose']
                del self.regions[sys]['Muon']['Tight']['btagFail']
                del self.regions[sys]['Electron']['Tight']['btagFail']

        if isTTSL_:
            for sys in self.regions.keys():
                for lepton in self.regions[sys].keys():
                    for leptype in self.regions[sys][lepton].keys():
                        for btagregion in self.regions[sys][lepton][leptype].keys():
                            for jmult in self.regions[sys][lepton][leptype][btagregion].keys():
                                self.regions[sys][lepton][leptype][btagregion][jmult].update({'semilep'})


    @property
    def accumulator(self):
        return self._accumulator


    def make_hists(self):
        histo_dict = {}
        histo_dict['GenWeights'] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.genwt_axis)
        histo_dict['TotWeights'] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.totwt_axis)

        return histo_dict


    def process(self, events):
        np.random.seed(10) # sets seed so values from random distributions are reproducible (JER corrections)
        output = self.accumulator.identity()

        self.sample_name = events.metadata['dataset']

            ## make event weights
                # data or MC distinction made internally
        mu_evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections, isTTbar=isNominal_ttJets_)
        el_evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections, isTTbar=isNominal_ttJets_)

            ## initialize selections
        selection = {evt_sys: PackedSelection() for evt_sys in self.event_systematics_to_run}

            # get all passing leptons
        lep_and_filter_pass = objsel.select_leptons(events, year=args.year)#, cutflow=output['cutflow'])
        {selection[sys].add('lep_and_filter_pass', lep_and_filter_pass) for sys in selection.keys()} # add passing leptons requirement to all systematics

            ## build corrected jets and MET
        events['Jet'], events['MET'] = IDJet.process_jets(events, args.year, self.corrections['JetCor'])

        ## add different selections
                ## muons
        tight_mu_sel, loose_mu_sel = ak.sum(events['Muon']['TIGHTMU'], axis=1) == 1, ak.sum(events['Muon']['LOOSEMU'], axis=1) == 1
        {selection[sys].add('tight_MU', tight_mu_sel) for sys in selection.keys()} # one muon passing TIGHT criteria
        {selection[sys].add('loose_MU', loose_mu_sel) for sys in selection.keys()} # one muon passing LOOSE criteria
                ## electrons
        tight_el_sel, loose_el_sel = ak.sum(events['Electron']['TIGHTEL'], axis=1) == 1, ak.sum(events['Electron']['LOOSEEL'], axis=1) == 1
        {selection[sys].add('tight_EL', tight_el_sel) for sys in selection.keys()} # one electron passing TIGHT criteria
        {selection[sys].add('loose_EL', loose_el_sel) for sys in selection.keys()} # one electron passing LOOSE criteria

        ### apply lepton SFs to MC (only applicable to tight leptons)
        if 'LeptonSF' in corrections.keys():
            tight_muons = events['Muon'][tight_mu_sel][(events['Muon'][tight_mu_sel]['TIGHTMU'] == True)]
            muSFs_dict =  MCWeights.get_lepton_sf(year=args.year, lepton='Muons', corrections=self.corrections['LeptonSF'],
                pt=ak.flatten(tight_muons['pt']), eta=ak.flatten(tight_muons['eta']))
            mu_reco_cen = np.ones(len(events))
            mu_reco_err = np.zeros(len(events))
            mu_trig_cen = np.ones(len(events))
            mu_trig_err = np.zeros(len(events))
            mu_reco_cen[tight_mu_sel] = muSFs_dict['RECO_CEN']
            mu_reco_err[tight_mu_sel] = muSFs_dict['RECO_ERR']
            mu_trig_cen[tight_mu_sel] = muSFs_dict['TRIG_CEN']
            mu_trig_err[tight_mu_sel] = muSFs_dict['TRIG_ERR']
            mu_evt_weights.add('Lep_RECO', mu_reco_cen, mu_reco_err, mu_reco_err, shift=True)
            mu_evt_weights.add('Lep_TRIG', mu_trig_cen, mu_trig_err, mu_trig_err, shift=True)
    
            tight_electrons = events['Electron'][tight_el_sel][(events['Electron'][tight_el_sel]['TIGHTEL'] == True)]
            elSFs_dict = MCWeights.get_lepton_sf(year=args.year, lepton='Electrons', corrections=self.corrections['LeptonSF'],
                pt=ak.flatten(tight_electrons['pt']), eta=ak.flatten(tight_electrons['etaSC']))
            el_reco_cen = np.ones(len(events))
            el_reco_err = np.zeros(len(events))
            el_trig_cen = np.ones(len(events))
            el_trig_err = np.zeros(len(events))
            el_reco_cen[tight_el_sel] = elSFs_dict['RECO_CEN']
            el_reco_err[tight_el_sel] = elSFs_dict['RECO_ERR']
            el_trig_cen[tight_el_sel] = elSFs_dict['TRIG_CEN']
            el_trig_err[tight_el_sel] = elSFs_dict['TRIG_ERR']
            el_evt_weights.add('Lep_RECO', el_reco_cen, el_reco_err, el_reco_err, shift=True)
            el_evt_weights.add('Lep_TRIG', el_trig_cen, el_trig_err, el_trig_err, shift=True)

            # find gen level particles for ttbar system and other ttbar corrections
        if isTTbar_:
            genpsel.select(events, mode='NORMAL')
            {selection[sys].add('semilep', ak.num(events['SL']) > 0) for sys in selection.keys()}
            if 'Kin_Rewt' in self.corrections.keys():
                kin_wts = MCWeights.get_kin_weights(self.corrections['Kin_Rewt'], GenTTbar)
                mu_evt_weights.add('%s_reweighting' % corrections['Kin_Rewt']['Var'], kin_wts)
                el_evt_weights.add('%s_reweighting' % corrections['Kin_Rewt']['Var'], kin_wts)


        #if args.debug: set_trace()
            # run over systematics that require changes to event objects (jets+MET)
        for evt_sys in self.event_systematics_to_run:
            output['cutflow_%s' % evt_sys]['lep_and_filter_pass'] += ak.sum(lep_and_filter_pass)
                # jet selection
            passing_jets = objsel.jets_selection(events, year=args.year, cutflow=output['cutflow_%s' % evt_sys], shift=evt_sys)
            output['cutflow_%s' % evt_sys]['nEvts passing jet and lepton obj selection'] += ak.sum(passing_jets & lep_and_filter_pass)
            selection[evt_sys].add('passing_jets', passing_jets)
            selection[evt_sys].add('jets_3',  ak.num(events['SelectedJets']) == 3)
            selection[evt_sys].add('jets_4p',  ak.num(events['SelectedJets']) > 3) # only for getting btag weights
            selection[evt_sys].add('DeepCSV_pass', ak.sum(events['SelectedJets'][btag_wps[0]], axis=1) >= 2)

                # sort jets by btag value
            events['SelectedJets'] = events['SelectedJets'][ak.argsort(events['SelectedJets']['btagDeepB'], ascending=False)] if btaggers[0] == 'DeepCSV' else events['SelectedJets'][ak.argsort(events['SelectedJets']['btagDeepFlavB'], ascending=False)]

                # btag fail sideband
            deepcsv_sorted = events['SelectedJets'][ak.argsort(events['SelectedJets']['btagDeepB'], ascending=False)]['btagDeepB']
            valid_counts_inds = ak.where(ak.num(events['SelectedJets']) > 1)[0]
            deepcsv_fail = np.zeros(len(events)).astype(bool)
            deepcsv_fail[valid_counts_inds] = (deepcsv_sorted[valid_counts_inds][:, 0] < btag_values[args.year]['btagDeepB']['DeepCSV'+wps_to_use[0]]) & (deepcsv_sorted[valid_counts_inds][:, 1] < btag_values[args.year]['btagDeepB']['DeepCSV'+wps_to_use[0]])
            selection[evt_sys].add('DeepCSV_fail', deepcsv_fail ) # highest and second highest DeepCSV values don't pass tight and loose WPs

                ## apply btagging SFs to MC
            if (not isData_) and (corrections['BTagSF'] == True):
                deepcsv_cen   = np.ones(len(events))
                deepcsv_bc_up = np.ones(len(events))
                deepcsv_bc_dw = np.ones(len(events))
                deepcsv_l_up = np.ones(len(events))
                deepcsv_l_dw = np.ones(len(events))

                threeJets_cut = selection[evt_sys].require(lep_and_filter_pass=True, passing_jets=True, jets_3=True)
                deepcsv_3j_wts = self.corrections['BTag_Constructors']['DeepCSV']['3Jets'].get_scale_factor(jets=events['SelectedJets'][threeJets_cut], passing_cut='DeepCSV'+wps_to_use[0])
                deepcsv_cen[threeJets_cut] = ak.prod(deepcsv_3j_wts['central'], axis=1)
                deepcsv_bc_up[threeJets_cut] = ak.prod(deepcsv_3j_wts['bc_up'], axis=1)
                deepcsv_bc_dw[threeJets_cut] = ak.prod(deepcsv_3j_wts['bc_down'], axis=1)
                deepcsv_l_up[threeJets_cut] = ak.prod(deepcsv_3j_wts['udsg_up'], axis=1)
                deepcsv_l_dw[threeJets_cut] = ak.prod(deepcsv_3j_wts['udsg_down'], axis=1)
    
                fourplusJets_cut = selection[evt_sys].require(lep_and_filter_pass=True, passing_jets=True, jets_4p=True)
                deepcsv_4pj_wts = self.corrections['BTag_Constructors']['DeepCSV']['4PJets'].get_scale_factor(jets=events['SelectedJets'][fourplusJets_cut], passing_cut='DeepCSV'+wps_to_use[0])
                deepcsv_cen[fourplusJets_cut] = ak.prod(deepcsv_4pj_wts['central'], axis=1)
                deepcsv_bc_up[fourplusJets_cut] = ak.prod(deepcsv_4pj_wts['bc_up'], axis=1)
                deepcsv_bc_dw[fourplusJets_cut] = ak.prod(deepcsv_4pj_wts['bc_down'], axis=1)
                deepcsv_l_up[fourplusJets_cut] = ak.prod(deepcsv_4pj_wts['udsg_up'], axis=1)
                deepcsv_l_dw[fourplusJets_cut] = ak.prod(deepcsv_4pj_wts['udsg_down'], axis=1)
                    # make dict of btag weights
                btag_weights = {
                    'DeepCSV_CEN' : deepcsv_cen,
                    'DeepCSV_bc_UP' : deepcsv_bc_up,
                    'DeepCSV_bc_DW' : deepcsv_bc_dw,
                    'DeepCSV_l_UP' : deepcsv_l_up,
                    'DeepCSV_l_DW' : deepcsv_l_dw,
                }

            #set_trace()
            ## fill hists for each region
            for lepton in self.regions[evt_sys].keys():
                evt_weights = mu_evt_weights if lepton == 'Muon' else el_evt_weights
                for leptype in self.regions[evt_sys][lepton].keys():
                    for btagregion in self.regions[evt_sys][lepton][leptype].keys():
                        for jmult in self.regions[evt_sys][lepton][leptype][btagregion].keys():
                            cut = selection[evt_sys].all(*self.regions[evt_sys][lepton][leptype][btagregion][jmult])
                            #set_trace()

                            output['cutflow_%s' % evt_sys]['nEvts %s' % ', '.join([lepton, leptype, btagregion, jmult])] += cut.sum()

                            if args.debug: print(lepton, leptype, btagregion, jmult)
                            if cut.sum() > 0:
                                ltype = 'MU' if lepton == 'Muon' else 'EL'
                                if 'loose_or_tight_%s' % ltype in self.regions[evt_sys][lepton][leptype][btagregion][jmult]:
                                    leptons = events[lepton][cut][((events[lepton][cut]['TIGHT%s' % ltype] == True) | (events[lepton][cut]['LOOSE%s' % ltype] == True))]
                                elif 'tight_%s' % ltype in self.regions[evt_sys][lepton][leptype][btagregion][jmult]:
                                    leptons = events[lepton][cut][(events[lepton][cut]['TIGHT%s' % ltype] == True)]
                                elif 'loose_%s' % ltype in self.regions[evt_sys][lepton][leptype][btagregion][jmult]:
                                    leptons = events[lepton][cut][(events[lepton][cut]['LOOSE%s' % ltype] == True)]
                                else:
                                    raise ValueError("Not sure what lepton type to choose for event")

                                    # get jets and MET
                                jets, met = events['SelectedJets'][cut], events['SelectedMET'][cut]

                                    # find best permutations
                                best_perms = ttpermutator.find_best_permutations(jets=jets, leptons=leptons, MET=met, btagWP=btag_wps[0], btag_req=False if btagregion == 'btagFail' else True)
                                valid_perms = np.zeros(cut.size, dtype=bool)
                                valid_perms[cut] = ak.num(best_perms['TTbar'].pt) > 0
                                output['cutflow_%s' % evt_sys]['nEvts %s: valid perms' % ', '.join([lepton, leptype, btagregion, jmult])] += ak.sum(valid_perms)

                                bp_status = np.zeros(cut.size, dtype=int) # 0 == '' (no gen matching), 1 == 'right', 2 == 'matchable', 3 == 'unmatchable', 4 == 'sl_tau', 5 == 'noslep'
                                    # get matched permutation (semilep ttbar only)
                                if isTTbar_:
                                    semilep_evts = selection[evt_sys].require(semilep=True)
                                    bp_status[~semilep_evts] = 5
                                    if semilep_evts.sum() > 0:
                                            # find matched permutations
                                        mp = ttmatcher.best_match(gen_hyp=events['SL'][cut], jets=jets, leptons=leptons, met=met)
                                        perm_cat_array = compare_matched_best_perms(mp, best_perms, njets=jmult)
                                        bp_status[cut] = perm_cat_array
                                        sl_tau_evts = ak.where(ak.fill_none(ak.pad_none(np.abs(events['SL']['Lepton'].pdgId) == 15, 1), False) == True)[0]
                                        bp_status[sl_tau_evts] = 4

                                    ## create MT regions
                                MTHigh = np.zeros(cut.size, dtype=bool)
                                MT = make_vars.MT(leptons, met)
                                MTHigh[valid_perms] = ak.flatten(MT[valid_perms[cut]] >= MTcut) 
                                output['cutflow_%s' % evt_sys]['nEvts %s: pass MT cut' % ', '.join([lepton, leptype, btagregion, jmult])] += ak.sum(MTHigh)

                                    # fill hists for each systematic
                                if args.debug: print('  evt sys:', evt_sys)
                                #set_trace()
                                if evt_sys == 'nosys':
                                    for rewt_sys in self.reweight_systematics_to_run:
                                        #if args.debug: set_trace()
                                            ## only fill plots in signal region if systematic variation being used
                                        if (rewt_sys != 'nosys') and ('%s_%s' % (leptype, btagregion) != 'Tight_btagPass'): continue

                                        sysname = events.metadata['dataset'].split('_')[-1] if isTTShift_ else rewt_sys
                                        if args.debug: print('    sysname:', sysname)

                                        #set_trace()
                                        if rewt_sys == 'nosys':
                                            wts = (evt_weights.weight()*btag_weights['%s_CEN' % btaggers[0]])[MTHigh]
                                        elif rewt_sys.startswith('btag'):
                                            wts = (evt_weights.weight()*btag_weights[rewt_sys.replace('btag', btaggers[0])])[MTHigh]
                                        else:
                                            wts = (evt_weights.weight(rewt_sys)*btag_weights['%s_CEN' % btaggers[0]])[MTHigh]

                                            # fill hists for interference samples
                                        if isInt_:
                                            #set_trace()
                                                # fill hists for positive weights
                                            pos_evts = np.where(wts > 0)
                                            self.sample_name = '%s_pos' % events.metadata['dataset']
                                            output = self.fill_hists(acc=output, sys=sysname, jetmult=jmult, leptype=lepton, lepcat=leptype, btagregion=btagregion, permarray=bp_status[MTHigh][pos_evts],
                                                gen_wts=events['genWeight'][MTHigh][pos_evts], tot_wts=wts[pos_evts])
                                                # fill hists for negative weights
                                            neg_evts = np.where(wts < 0)
                                            self.sample_name = '%s_neg' % events.metadata['dataset']
                                            output = self.fill_hists(acc=output, sys=sysname, jetmult=jmult, leptype=lepton, lepcat=leptype, btagregion=btagregion, permarray=bp_status[MTHigh][neg_evts],
                                                gen_wts=events['genWeight'][MTHigh][neg_evts], tot_wts=wts[neg_evts])
                                        else:
                                            output = self.fill_hists(acc=output, sys=sysname, jetmult=jmult, leptype=lepton, lepcat=leptype, btagregion=btagregion, permarray=bp_status[MTHigh],
                                                gen_wts=events['genWeight'][MTHigh], tot_wts=wts)

                                else:
                                    if args.debug: print('    sysname:', evt_sys)
                                    wts = (evt_weights.weight()*btag_weights['%s_CEN' % btaggers[0]])[MTHigh]
                                            # fill hists for interference samples
                                    if isInt_:
                                        #set_trace()
                                            # fill hists for positive weights
                                        pos_evts = np.where(wts > 0)
                                        self.sample_name = '%s_pos' % events.metadata['dataset']
                                        output = self.fill_hists(acc=output, sys=evt_sys, jetmult=jmult, leptype=lepton, lepcat=leptype, btagregion=btagregion, permarray=bp_status[MTHigh][pos_evts],
                                            gen_wts=events['genWeight'][MTHigh][pos_evts], tot_wts=wts[pos_evts])
                                            # fill hists for negative weights
                                        neg_evts = np.where(wts < 0)
                                        self.sample_name = '%s_neg' % events.metadata['dataset']
                                        output = self.fill_hists(acc=output, sys=evt_sys, jetmult=jmult, leptype=lepton, lepcat=leptype, btagregion=btagregion, permarray=bp_status[MTHigh][neg_evts],
                                            gen_wts=events['genWeight'][MTHigh][neg_evts], tot_wts=wts[neg_evts])
                                    else:
                                        output = self.fill_hists(acc=output, sys=evt_sys, jetmult=jmult, leptype=lepton, lepcat=leptype, btagregion=btagregion, permarray=bp_status[MTHigh],
                                            gen_wts=events['genWeight'][MTHigh], tot_wts=wts)

        return output

    def fill_hists(self, acc, sys, jetmult, leptype, lepcat, btagregion, permarray, gen_wts, tot_wts):
        #set_trace()
        for permval in np.unique(permarray).tolist():
            perm_inds = np.where(permarray == permval)
            dataset_name = '%s_%s' % (self.sample_name, perm_cats[permval]) if permval != 0 else self.sample_name
            acc['GenWeights'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, genwt=gen_wts[perm_inds])
            acc['TotWeights'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, totwt=tot_wts[perm_inds])

        return acc        

    def postprocess(self, accumulator):
        return accumulator


proc_executor = processor.iterative_executor if args.debug else processor.futures_executor
proc_exec_args = {"schema": NanoAODSchema} if args.debug else {"schema": NanoAODSchema, "workers": 8}
output = processor.run_uproot_job(
    fileset,
    treename="Events",
    processor_instance=evtWeights(),
    executor=proc_executor,
    executor_args=proc_exec_args,
    #chunksize=10000 if args.debug else 100000,
    chunksize=100000,
)

save(output, args.outfname)
print('%s has been written' % args.outfname)
