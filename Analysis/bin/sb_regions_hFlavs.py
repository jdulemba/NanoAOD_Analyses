#!/usr/bin/env python

import time
tic = time.time()

from coffea import hist, processor
from coffea.nanoevents import NanoAODSchema
from coffea.analysis_tools import PackedSelection
import awkward as ak
from coffea.nanoevents.methods import vector
ak.behavior.update(vector.behavior)

from pdb import set_trace
from coffea.util import save, load
import os
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
analyzer = 'sb_regions_hFlavs'

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('fset', type=str, help='Fileset dictionary (in string form) to be used for the processor')
parser.add_argument('year', choices=['2016APV', '2016', '2017', '2018'] if base_jobid == 'ULnanoAOD' else ['2016', '2017', '2018'], help='Specify which year to run over')
parser.add_argument('outfname', type=str, help='Specify output filename, including directory and file extension')
parser.add_argument('opts', type=str, help='Fileset dictionary (in string form) to be used for the processor')
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
isSE_Data_ = samplename.startswith('data_SingleElectron')
isSM_Data_ = samplename.startswith('data_SingleMuon')
if isData_:
    lumiMask_path = os.path.join(proj_dir, 'inputs', 'data', base_jobid, 'LumiMasks', '%s_GoldenJson_%s.txt' % (args.year, base_jobid))

Nominal_ttJets = ['ttJets_PS'] if ((args.year == '2016') and (base_jobid == 'NanoAODv6')) else ['ttJetsSL', 'ttJetsHad', 'ttJetsDiLep']
isNominal_ttJets_ = samplename in Nominal_ttJets
isTTShift_ = isTTbar_ and (samplename.endswith('UP') or samplename.endswith('DOWN') or 'mtop' in samplename)

# convert input string of options dictionary to actual dictionary
odict = (args.opts).replace("\'", "\"")
opts_dict = prettyjson.loads(odict)

    ## set config options passed through argparse
import ast
to_debug = ast.literal_eval(opts_dict.get('debug', 'False'))
apply_hem = ast.literal_eval(opts_dict.get('apply_hem', 'True'))

#set_trace()
## init tt probs for likelihoods
ttpermutator.year_to_run(year=args.year)

## load lumimask for data and corrections for event weights
pu_correction = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'MC_PU_Weights.coffea'))[args.year]
lepSF_correction = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'leptonSFs.coffea'))[args.year]
jet_corrections = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'JetMETCorrections.coffea'))[args.year]
alpha_corrections = load(os.path.join(proj_dir, 'Corrections', base_jobid,'alpha_correction_%s.coffea' % jobid))[args.year]['E']['All_2D'] # E, All_2D determined by post application plots
nnlo_var = 'mtt_vs_thad_ctstar_Interp'
nnlo_reweighting = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'NNLO_to_Tune_Orig_Interp_Ratios_%s.coffea' % base_jobid))[args.year]
corrections = {
    'Pileup' : pu_correction,
    'Prefire' : True,
    'LeptonSF' : lepSF_correction,
    'BTagSF' : not isData_,
    'JetCor' : jet_corrections,
    'Alpha' : alpha_corrections,
    'NNLO_Rewt' : {'Var' : nnlo_var, 'Correction' : nnlo_reweighting[nnlo_var]},
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
evt_sys_to_run = opts_dict.get('evt_sys', 'NONE').upper()
rewt_sys_to_run = opts_dict.get('rewt_sys', 'NONE').upper()
only_sys = ast.literal_eval(opts_dict.get('only_sys', 'False'))
if isData_ or isTTShift_: # data or separate systematics dataset
    event_systematics_to_run = ['nosys']
    reweight_systematics_to_run = ['nosys']
else: # MC samples
    import fnmatch
    if only_sys: # don't run 'nosys'
        if (evt_sys_to_run == 'NONE') and (rewt_sys_to_run == 'NONE'):
            raise ValueError("At least one systematic must be specified in order to run only on systematics!")
    
        event_systematics_to_run = ['nosys'] if (rewt_sys_to_run != 'NONE') else []
        reweight_systematics_to_run = []
    
    else:
        event_systematics_to_run = ['nosys']
        reweight_systematics_to_run = ['nosys']

    event_systematics_to_run += [systematics.event_sys_opts[args.year][name] for name in systematics.event_sys_opts[args.year].keys() if fnmatch.fnmatch(name, evt_sys_to_run)]
    if isSignal_:
        reweight_systematics_to_run += [systematics.signal_reweight_opts[args.year][name] for name in systematics.signal_reweight_opts[args.year].keys() if fnmatch.fnmatch(name, rewt_sys_to_run)]
    else:
        reweight_systematics_to_run += [systematics.reweight_sys_opts[args.year][name] for name in systematics.reweight_sys_opts[args.year].keys() if fnmatch.fnmatch(name, rewt_sys_to_run)]

        ## check that systematics only related to ttbar events aren't used for non-ttbar events
    if not isTTbar_:
        reweight_systematics_to_run = [sys for sys in reweight_systematics_to_run if sys not in systematics.ttJets_sys.values()]
    
print('Running with event systematics:', *sorted(set(event_systematics_to_run).difference(set(['nosys']))), sep=', ') if 'nosys' in event_systematics_to_run else print('Running with event systematics:', *sorted(event_systematics_to_run), sep=', ')
print('  and reweight systematics:', *sorted(reweight_systematics_to_run), sep=', ')
#set_trace()

#btag_regions = {
#    '0p05' : (0.0, 0.05),
#    '0p10' : (0.05, 0.10),
#    '0p15' : (0.10, 0.15),
#    '0p20' : (0.15, 0.20),
#    '0p25' : (0.20, 0.25),
#    '0p30' : (0.25, 0.30),
#    '0p35' : (0.30, 0.35),
#    '0p40' : (0.35, 0.40),
#    '0p45' : (0.40, 0.45),
#    '0p50' : (0.45, 0.50),
#}
btag_regions = {
    'p00p15' : (0.00, 0.15),
    'p15p30' : (0.15, 0.30),
    'p30p45' : (0.30, 0.45),
    #'p00p05' : (0.00, 0.05),
    #'p05p10' : (0.05, 0.10),
    #'p10p15' : (0.10, 0.15),
    #'p15p20' : (0.15, 0.20),
    #'p20p25' : (0.20, 0.25),
    #'p25p30' : (0.25, 0.30),
    #'p30p35' : (0.30, 0.35),
    #'p35p40' : (0.35, 0.40),
    #'p40p45' : (0.40, 0.45),
}

# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class sb_regions_hFlavs(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.sys_axis = hist.Cat("sys", "Systematic")
        self.jetmult_axis = hist.Cat("jmult", "nJets")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")
        self.btag_axis = hist.Cat("btag", "BTag Region")
        self.hflav_axis = hist.Cat("hFlav", "Hadron Flavour")
        self.deepcsv_discr_axis = hist.Bin("deepcsv_bdisc", "DeepCSV bDiscr", 100, 0., 1.)
        #self.deepjet_discr_axis = hist.Bin("deepjet_bdisc", "DeepJet bDiscr", 100, 0., 1.)
        self.pt_axis = hist.Bin("pt", "p_{T} [GeV]", 100, 0, 1000)
        self.eta_axis = hist.Bin("eta", r"$\eta$", 100, -5., 5.)

            ## make dictionary of hists
        histo_dict = {}
                ## make jet hists
        jet_hists = self.make_jet_hists()
        histo_dict.update(jet_hists)

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
                'btagPass' : {
                    '3Jets'  : {'lep_and_filter_pass', 'passing_jets', 'jets_3' , 'tight_MU', 'DeepCSV_pass'},
                    '4PJets' : {'lep_and_filter_pass', 'passing_jets', 'jets_4p', 'tight_MU', 'DeepCSV_pass'},
                },
            },
            'Electron' : {
                'btagPass' : {
                    '3Jets'  : {'lep_and_filter_pass', 'passing_jets', 'jets_3' , 'tight_EL', 'DeepCSV_pass'},
                    '4PJets' : {'lep_and_filter_pass', 'passing_jets', 'jets_4p', 'tight_EL', 'DeepCSV_pass'},
                },
            },
        }
        self.regions = {sys:deepcopy(base_regions) for sys in self.event_systematics_to_run}
            # add sideband regions for nosys
        self.regions['nosys']['Muon'].update({
            key : {
                '3Jets'  : {'lep_and_filter_pass', 'passing_jets', 'jets_3' , 'tight_MU', 'DeepCSV_%s' % key},
                '4PJets' : {'lep_and_filter_pass', 'passing_jets', 'jets_4p', 'tight_MU', 'DeepCSV_%s' % key},
            } for key in btag_regions.keys()
        })
        self.regions['nosys']['Electron'].update({
            key : {
                '3Jets'  : {'lep_and_filter_pass', 'passing_jets', 'jets_3' , 'tight_EL', 'DeepCSV_%s' % key},
                '4PJets' : {'lep_and_filter_pass', 'passing_jets', 'jets_4p', 'tight_EL', 'DeepCSV_%s' % key},
            } for key in btag_regions.keys()
        })

        if isSM_Data_:
            if 'Electron' in self.regions['nosys'].keys(): del self.regions['nosys']['Electron']
        if isSE_Data_:
            if 'Muon' in self.regions['nosys'].keys(): del self.regions['nosys']['Muon']
        if isData_:
            for lepton in self.regions['nosys'].keys():
                for btagregion in self.regions['nosys'][lepton].keys():
                    for jmult in self.regions['nosys'][lepton][btagregion].keys():
                        self.regions['nosys'][lepton][btagregion][jmult].update({'lumimask'})

        if isTTSL_:
            for sys in self.regions.keys():
                for lepton in self.regions[sys].keys():
                    for btagregion in self.regions[sys][lepton].keys():
                        for jmult in self.regions[sys][lepton][btagregion].keys():
                            self.regions[sys][lepton][btagregion][jmult].update({'semilep'})


    @property
    def accumulator(self):
        return self._accumulator


    def make_jet_hists(self):
        histo_dict = {}
                # disc vs pt (eta) for all jets
        histo_dict['DeepCSV_bDisc_vs_Jet_pt']  = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.hflav_axis, self.pt_axis, self.deepcsv_discr_axis)
        histo_dict['DeepCSV_bDisc_vs_Jet_eta'] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.hflav_axis, self.eta_axis, self.deepcsv_discr_axis)
                # disc vs pt (eta) for jet with 1st highest disc value
        histo_dict['DeepCSV_bDisc_vs_Jet1_pt']  = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.hflav_axis, self.pt_axis, self.deepcsv_discr_axis)
        histo_dict['DeepCSV_bDisc_vs_Jet1_eta'] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.hflav_axis, self.eta_axis, self.deepcsv_discr_axis)
                # disc vs pt (eta) for jet with 2nd highest disc value
        histo_dict['DeepCSV_bDisc_vs_Jet2_pt']  = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.hflav_axis, self.pt_axis, self.deepcsv_discr_axis)
        histo_dict['DeepCSV_bDisc_vs_Jet2_eta'] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.hflav_axis, self.eta_axis, self.deepcsv_discr_axis)
                # disc vs pt (eta) for jet with 3rd highest disc value
        histo_dict['DeepCSV_bDisc_vs_Jet3_pt']  = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.hflav_axis, self.pt_axis, self.deepcsv_discr_axis)
        histo_dict['DeepCSV_bDisc_vs_Jet3_eta'] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.hflav_axis, self.eta_axis, self.deepcsv_discr_axis)

            # plots involving DeepCSV and DeepJet b-tagging discriminant values
        histo_dict['DeepCSV_bDisc'] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.hflav_axis, self.deepcsv_discr_axis)

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
        lep_and_filter_pass = objsel.select_leptons(events, year=args.year, hem_15_16=apply_hem)
        {selection[sys].add('lep_and_filter_pass', lep_and_filter_pass) for sys in selection.keys()} # add passing leptons requirement to all systematics

            ## build corrected jets and MET
        events['Jet'], events['MET'] = IDJet.process_jets(events, args.year, self.corrections['JetCor'])

        if isData_:
            runs = events.run
            lumis = events.luminosityBlock
            Golden_Json_LumiMask = lumi_tools.LumiMask(lumiMask_path)
            LumiMask = Golden_Json_LumiMask.__call__(runs, lumis) ## returns array of valid events
            {selection[sys].add('lumimask', LumiMask) for sys in selection.keys()}
                  ## object selection and add different selections
            if isSM_Data_:
                        ## muons
                {selection[sys].add('tight_MU', ak.sum(events['Muon']['TIGHTMU'], axis=1) == 1) for sys in selection.keys()} # one muon passing TIGHT criteria
                #{selection[sys].add('loose_MU', ak.sum(events['Muon']['LOOSEMU'], axis=1) == 1) for sys in selection.keys()} # one muon passing LOOSE criteria
            if isSE_Data_:
                        ## electrons
                {selection[sys].add('tight_EL', ak.sum(events['Electron']['TIGHTEL'], axis=1) == 1) for sys in selection.keys()} # one electron passing TIGHT criteria
                #{selection[sys].add('loose_EL', ak.sum(events['Electron']['LOOSEEL'], axis=1) == 1) for sys in selection.keys()} # one electron passing LOOSE criteria

        else:
            ## add different selections
                    ## muons
            tight_mu_sel, loose_mu_sel = ak.sum(events['Muon']['TIGHTMU'], axis=1) == 1, ak.sum(events['Muon']['LOOSEMU'], axis=1) == 1
            {selection[sys].add('tight_MU', tight_mu_sel) for sys in selection.keys()} # one muon passing TIGHT criteria
            #{selection[sys].add('loose_MU', loose_mu_sel) for sys in selection.keys()} # one muon passing LOOSE criteria
                    ## electrons
            tight_el_sel, loose_el_sel = ak.sum(events['Electron']['TIGHTEL'], axis=1) == 1, ak.sum(events['Electron']['LOOSEEL'], axis=1) == 1
            {selection[sys].add('tight_EL', tight_el_sel) for sys in selection.keys()} # one electron passing TIGHT criteria
            #{selection[sys].add('loose_EL', loose_el_sel) for sys in selection.keys()} # one electron passing LOOSE criteria

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
            if isTTSL_:
                genpsel.select(events, mode='NORMAL')
                {selection[sys].add('semilep', ak.num(events['SL']) > 0) for sys in selection.keys()}
            else:
                {selection[sys].add('semilep', np.zeros(len(events), dtype=bool)) for sys in selection.keys()}
            if 'NNLO_Rewt' in self.corrections.keys():
                if not isTTSL_:
                    genpsel.select(events, mode='NORMAL')
                    {selection[sys].add('semilep', ak.num(events['SL']) > 0) for sys in selection.keys()}
                nnlo_wts = MCWeights.get_nnlo_weights(self.corrections['NNLO_Rewt'], events)
                mu_evt_weights.add('%s_reweighting' % corrections['NNLO_Rewt']['Var'], nnlo_wts)
                el_evt_weights.add('%s_reweighting' % corrections['NNLO_Rewt']['Var'], nnlo_wts)


        #if to_debug: set_trace()
            # run over systematics that require changes to event objects (jets+MET)
        for evt_sys in self.event_systematics_to_run:
            output['cutflow_%s' % evt_sys]['lep_and_filter_pass'] += ak.sum(lep_and_filter_pass)
                # jet selection
            passing_jets = objsel.jets_selection(events, year=args.year, cutflow=output['cutflow_%s' % evt_sys], shift=evt_sys, hem_15_16=apply_hem)
            output['cutflow_%s' % evt_sys]['nEvts passing jet and lepton obj selection'] += ak.sum(passing_jets & lep_and_filter_pass)
            selection[evt_sys].add('passing_jets', passing_jets)
            selection[evt_sys].add('jets_3',  ak.num(events['SelectedJets']) == 3)
            selection[evt_sys].add('jets_4p',  ak.num(events['SelectedJets']) > 3) # only for getting btag weights
            selection[evt_sys].add('DeepCSV_pass', ak.sum(events['SelectedJets'][btag_wps[0]], axis=1) >= 2)

                # sort jets by btag value
            events['SelectedJets'] = events['SelectedJets'][ak.argsort(events['SelectedJets']['btagDeepB'], ascending=False)] if btaggers[0] == 'DeepCSV' else events['SelectedJets'][ak.argsort(events['SelectedJets']['btagDeepFlavB'], ascending=False)]

                # btag sidebands
            deepcsv_sorted = events['SelectedJets'][ak.argsort(events['SelectedJets']['btagDeepB'], ascending=False)]['btagDeepB']
            for btag_reg, (low_bound, up_bound) in btag_regions.items():
                selection[evt_sys].add('DeepCSV_%s' % btag_reg, (ak.max(deepcsv_sorted, axis=1) > low_bound) & (ak.max(deepcsv_sorted, axis=1) <= up_bound))

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
                for btagregion in self.regions[evt_sys][lepton].keys():
                    for jmult in self.regions[evt_sys][lepton][btagregion].keys():
                        #set_trace()
                        cut = selection[evt_sys].all(*self.regions[evt_sys][lepton][btagregion][jmult])

                        output['cutflow_%s' % evt_sys]['nEvts %s' % ', '.join([lepton, btagregion, jmult])] += cut.sum()

                        if to_debug: print(lepton, btagregion, jmult)
                        if cut.sum() > 0:
                            ltype = 'MU' if lepton == 'Muon' else 'EL'
                            if 'loose_or_tight_%s' % ltype in self.regions[evt_sys][lepton][btagregion][jmult]:
                                leptons = events[lepton][cut][((events[lepton][cut]['TIGHT%s' % ltype] == True) | (events[lepton][cut]['LOOSE%s' % ltype] == True))]
                            elif 'tight_%s' % ltype in self.regions[evt_sys][lepton][btagregion][jmult]:
                                leptons = events[lepton][cut][(events[lepton][cut]['TIGHT%s' % ltype] == True)]
                            elif 'loose_%s' % ltype in self.regions[evt_sys][lepton][btagregion][jmult]:
                                leptons = events[lepton][cut][(events[lepton][cut]['LOOSE%s' % ltype] == True)]
                            else:
                                raise ValueError("Not sure what lepton type to choose for event")

                                # get jets and MET
                            jets, met = events['SelectedJets'][cut], events['SelectedMET'][cut]

                                # find best permutations
                            best_perms = ttpermutator.find_best_permutations(jets=jets, leptons=leptons, MET=met, btagWP=btag_wps[0], btag_req=True if btagregion == 'btagPass' else False)
                            valid_perms = ak.num(best_perms['TTbar'].pt) > 0
                            output['cutflow_%s' % evt_sys]['nEvts %s: valid perms' % ', '.join([lepton, btagregion, jmult])] += ak.sum(valid_perms)

                            bp_status = np.zeros(cut.size, dtype=int) # 0 == '' (no gen matching), 1 == 'right', 2 == 'matchable', 3 == 'unmatchable', 4 == 'sl_tau', 5 == 'noslep'
                                # get matched permutation (semilep ttbar only)
                            if isTTbar_:
                                ## find gen level particles for ttbar system and other ttbar corrections
                                #set_trace()
                                #if isTTSL_:
                                #    genpsel.select(events, mode='NORMAL')
                                #    {selection[sys].add('semilep', ak.num(events['SL']) > 0) for sys in selection.keys()}
                                #else:
                                #    {selection[sys].add('semilep', np.zeros(len(events), dtype=bool)) for sys in selection.keys()}
                                #if 'NNLO_Rewt' in self.corrections.keys():
                                #    if not isTTSL_:
                                #        genpsel.select(events, mode='NORMAL')
                                #        {selection[sys].add('semilep', ak.num(events['SL']) > 0) for sys in selection.keys()}
                                #    nnlo_wts = MCWeights.get_nnlo_weights(self.corrections['NNLO_Rewt'], events)
                                #    mu_evt_weights.add('%s_reweighting' % corrections['NNLO_Rewt']['Var'], nnlo_wts)
                                #    el_evt_weights.add('%s_reweighting' % corrections['NNLO_Rewt']['Var'], nnlo_wts)
                                ###

                                semilep_evts = selection[evt_sys].require(semilep=True)
                                bp_status[~semilep_evts] = 5
                                if semilep_evts.sum() > 0:
                                        # find matched permutations
                                    mp = ttmatcher.best_match(gen_hyp=events['SL'][cut], jets=jets, leptons=leptons, met=met)
                                    #set_trace()
                                    perm_cat_array = compare_matched_best_perms(mp, best_perms, njets=jmult)
                                    bp_status[cut] = perm_cat_array
                                    if ak.any(ak.num(events['SL']['Lepton'].pdgId) != 1): raise ValueError("Number of leptons is incorrect for classifying tau+jets events")
                                    sl_tau_evts = ak.where(np.abs(events['SL']['Lepton'].pdgId) == 15)[0]
                                    bp_status[sl_tau_evts] = 4

                                ## create MT regions
                            MT = make_vars.MT(leptons, met)
                            MTHigh = ak.flatten(MT[valid_perms] >= MTcut)
                            output['cutflow_%s' % evt_sys]['nEvts %s: pass MT cut' % ', '.join([lepton, btagregion, jmult])] += ak.sum(MTHigh)

                                # fill hists for each systematic
                            if to_debug: print('  evt sys:', evt_sys)
                            #set_trace()
                            if evt_sys == 'nosys':
                                for rewt_sys in self.reweight_systematics_to_run:
                                    #if to_debug: set_trace()
                                        ## only fill plots in signal region if systematic variation being used
                                    if (rewt_sys != 'nosys') and (btagregion != 'btagPass'): continue

                                    if to_debug: print('    sysname:', rewt_sys)

                                    if rewt_sys == 'nosys':
                                        wts = evt_weights.weight()[cut][valid_perms][MTHigh] if isData_ else (evt_weights.weight()*btag_weights['%s_CEN' % btaggers[0]])[cut][valid_perms][MTHigh]
                                    elif rewt_sys.startswith('btag'):
                                        wts = (evt_weights.weight()*btag_weights[rewt_sys.replace('btag', btaggers[0])])[cut][valid_perms][MTHigh]
                                    else:
                                        wts = (evt_weights.weight(rewt_sys)*btag_weights['%s_CEN' % btaggers[0]])[cut][valid_perms][MTHigh]

                                    sysname = events.metadata['dataset'].split('_')[-1] if isTTShift_ else rewt_sys

                                        # fill hists for interference samples
                                    if isInt_:
                                        #set_trace()
                                            # fill hists for positive weights
                                        pos_evts = np.where(wts > 0)
                                        self.sample_name = '%s_pos' % events.metadata['dataset']
                                        if not isData_:
                                            output = self.fill_mc_hists(acc=output, sys=sysname, jetmult=jmult, leptype=lepton, btagregion=btagregion, permarray=bp_status[cut][valid_perms][MTHigh][pos_evts],
                                                perm=best_perms[valid_perms][MTHigh][pos_evts], jets=jets[valid_perms][MTHigh][pos_evts], leptons=leptons[valid_perms][MTHigh][pos_evts], MTvals=MT[valid_perms][MTHigh][pos_evts], evt_wts=wts[pos_evts])
                                        else:
                                            output = self.fill_data_hists(acc=output, sys=sysname, jetmult=jmult, leptype=lepton, btagregion=btagregion, permarray=bp_status[cut][valid_perms][MTHigh][pos_evts],
                                                perm=best_perms[valid_perms][MTHigh][pos_evts], jets=jets[valid_perms][MTHigh][pos_evts], leptons=leptons[valid_perms][MTHigh][pos_evts], MTvals=MT[valid_perms][MTHigh][pos_evts], evt_wts=wts[pos_evts])
                                            # fill hists for negative weights
                                        neg_evts = np.where(wts < 0)
                                        self.sample_name = '%s_neg' % events.metadata['dataset']
                                        if not isData_:
                                            output = self.fill_mc_hists(acc=output, sys=sysname, jetmult=jmult, leptype=lepton, btagregion=btagregion, permarray=bp_status[cut][valid_perms][MTHigh][neg_evts],
                                                perm=best_perms[valid_perms][MTHigh][neg_evts], jets=jets[valid_perms][MTHigh][neg_evts], leptons=leptons[valid_perms][MTHigh][neg_evts], MTvals=MT[valid_perms][MTHigh][neg_evts], evt_wts=wts[neg_evts])
                                        else:
                                            output = self.fill_data_hists(acc=output, sys=sysname, jetmult=jmult, leptype=lepton, btagregion=btagregion, permarray=bp_status[cut][valid_perms][MTHigh][neg_evts],
                                                perm=best_perms[valid_perms][MTHigh][neg_evts], jets=jets[valid_perms][MTHigh][neg_evts], leptons=leptons[valid_perms][MTHigh][neg_evts], MTvals=MT[valid_perms][MTHigh][neg_evts], evt_wts=wts[neg_evts])
                                    else:
                                        if not isData_:
                                            output = self.fill_mc_hists(acc=output, sys=sysname, jetmult=jmult, leptype=lepton, btagregion=btagregion, permarray=bp_status[cut][valid_perms][MTHigh],
                                                perm=best_perms[valid_perms][MTHigh], jets=jets[valid_perms][MTHigh], leptons=leptons[valid_perms][MTHigh], MTvals=MT[valid_perms][MTHigh], evt_wts=wts)
                                        else:
                                            output = self.fill_data_hists(acc=output, sys=sysname, jetmult=jmult, leptype=lepton, btagregion=btagregion, permarray=bp_status[cut][valid_perms][MTHigh],
                                                perm=best_perms[valid_perms][MTHigh], jets=jets[valid_perms][MTHigh], leptons=leptons[valid_perms][MTHigh], MTvals=MT[valid_perms][MTHigh], evt_wts=wts)

                            else:
                                if to_debug: print('    sysname:', evt_sys)
                                #if to_debug: set_trace()
                                wts = (evt_weights.weight()*btag_weights['%s_CEN' % btaggers[0]])[cut][valid_perms][MTHigh]
                                        # fill hists for interference samples
                                if isInt_:
                                    #set_trace()
                                        # fill hists for positive weights
                                    pos_evts = np.where(wts > 0)
                                    self.sample_name = '%s_pos' % events.metadata['dataset']
                                    if not isData_:
                                        output = self.fill_mc_hists(acc=output, sys=evt_sys, jetmult=jmult, leptype=lepton, btagregion=btagregion, permarray=bp_status[cut][valid_perms][MTHigh][pos_evts],
                                            perm=best_perms[valid_perms][MTHigh][pos_evts], jets=jets[valid_perms][MTHigh][pos_evts], leptons=leptons[valid_perms][MTHigh][pos_evts], MTvals=MT[valid_perms][MTHigh][pos_evts], evt_wts=wts[pos_evts])
                                    else:
                                        output = self.fill_data_hists(acc=output, sys=evt_sys, jetmult=jmult, leptype=lepton, btagregion=btagregion, permarray=bp_status[cut][valid_perms][MTHigh][pos_evts],
                                            perm=best_perms[valid_perms][MTHigh][pos_evts], jets=jets[valid_perms][MTHigh][pos_evts], leptons=leptons[valid_perms][MTHigh][pos_evts], MTvals=MT[valid_perms][MTHigh][pos_evts], evt_wts=wts[pos_evts])
                                        # fill hists for negative weights
                                    neg_evts = np.where(wts < 0)
                                    self.sample_name = '%s_neg' % events.metadata['dataset']
                                    if not isData_:
                                        output = self.fill_mc_hists(acc=output, sys=evt_sys, jetmult=jmult, leptype=lepton, btagregion=btagregion, permarray=bp_status[cut][valid_perms][MTHigh][neg_evts],
                                            perm=best_perms[valid_perms][MTHigh][neg_evts], jets=jets[valid_perms][MTHigh][neg_evts], leptons=leptons[valid_perms][MTHigh][neg_evts], MTvals=MT[valid_perms][MTHigh][neg_evts], evt_wts=wts[neg_evts])
                                    else:
                                        output = self.fill_data_hists(acc=output, sys=evt_sys, jetmult=jmult, leptype=lepton, btagregion=btagregion, permarray=bp_status[cut][valid_perms][MTHigh][neg_evts],
                                            perm=best_perms[valid_perms][MTHigh][neg_evts], jets=jets[valid_perms][MTHigh][neg_evts], leptons=leptons[valid_perms][MTHigh][neg_evts], MTvals=MT[valid_perms][MTHigh][neg_evts], evt_wts=wts[neg_evts])
                                else:
                                    if not isData_:
                                        output = self.fill_mc_hists(acc=output, sys=evt_sys, jetmult=jmult, leptype=lepton, btagregion=btagregion, permarray=bp_status[cut][valid_perms][MTHigh],
                                            perm=best_perms[valid_perms][MTHigh], jets=jets[valid_perms][MTHigh], leptons=leptons[valid_perms][MTHigh], MTvals=MT[valid_perms][MTHigh], evt_wts=wts)
                                    else:
                                        output = self.fill_data_hists(acc=output, sys=evt_sys, jetmult=jmult, leptype=lepton, btagregion=btagregion, permarray=bp_status[cut][valid_perms][MTHigh],
                                            perm=best_perms[valid_perms][MTHigh], jets=jets[valid_perms][MTHigh], leptons=leptons[valid_perms][MTHigh], MTvals=MT[valid_perms][MTHigh], evt_wts=wts)

        return output


    def fill_data_hists(self, acc, sys, jetmult, leptype, btagregion, permarray, perm, jets, leptons, MTvals, evt_wts):
            ## apply alpha correction for 3Jets
        if (jetmult == '3Jets') and ('Alpha' in self.corrections):
            alpha_corr = self.corrections['Alpha'](172.5/perm['THad'].mass)
            perm['THad'] = perm['THad'].multiply(alpha_corr) # correct thad
            perm['TTbar'] = ak.flatten(perm['THad']+perm['TLep']) # correct ttbar


        pt_sorted_jets = jets[ak.argsort(jets.pt, ascending=False)]
        deepcsv_sorted_jets = jets[ak.argsort(jets.btagDeepB, ascending=False)]
        deepjet_sorted_jets = jets[ak.argsort(jets.btagDeepFlavB, ascending=False)]

        for permval in np.unique(permarray).tolist():
            perm_inds = np.where(permarray == permval)
            dataset_name = '%s_%s' % (self.sample_name, perm_cats[permval]) if permval != 0 else self.sample_name

            #set_trace()
                # jets are b jet
            acc['DeepCSV_bDisc_vs_Jet_pt'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='bjet',
                pt=ak.flatten(jets[perm_inds].pt), deepcsv_bdisc=ak.flatten(jets[perm_inds].btagDeepB), weight=ak.flatten((ak.ones_like(jets.pt)*evt_wts)[perm_inds]))
            acc['DeepCSV_bDisc_vs_Jet_eta'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='bjet',
                eta=ak.flatten(jets[perm_inds].eta), deepcsv_bdisc=ak.flatten(jets[perm_inds].btagDeepB), weight=ak.flatten((ak.ones_like(jets.pt)*evt_wts)[perm_inds]))

                # jet with highest DeepCSV value (jet 1) is b jet
            acc['DeepCSV_bDisc_vs_Jet1_pt'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='bjet',
                pt=deepcsv_sorted_jets[perm_inds].pt[:, 0], deepcsv_bdisc=deepcsv_sorted_jets[perm_inds].btagDeepB[:, 0], weight=evt_wts[perm_inds])
            acc['DeepCSV_bDisc_vs_Jet1_eta'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='bjet',
                eta=deepcsv_sorted_jets[perm_inds].eta[:, 0], deepcsv_bdisc=deepcsv_sorted_jets[perm_inds].btagDeepB[:, 0], weight=evt_wts[perm_inds])

                # jet with second highest DeepCSV value (jet 2) is b jet
            acc['DeepCSV_bDisc_vs_Jet2_pt'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='bjet',
                pt=deepcsv_sorted_jets[perm_inds].pt[:, 1], deepcsv_bdisc=deepcsv_sorted_jets[perm_inds].btagDeepB[:, 1], weight=evt_wts[perm_inds])
            acc['DeepCSV_bDisc_vs_Jet2_eta'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='bjet',
                eta=deepcsv_sorted_jets[perm_inds].eta[:, 1], deepcsv_bdisc=deepcsv_sorted_jets[perm_inds].btagDeepB[:, 1], weight=evt_wts[perm_inds])

                # jet with third highest DeepCSV value (jet 3) is b jet
            acc['DeepCSV_bDisc_vs_Jet3_pt'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='bjet',
                pt=deepcsv_sorted_jets[perm_inds].pt[:, 2], deepcsv_bdisc=deepcsv_sorted_jets[perm_inds].btagDeepB[:, 2], weight=evt_wts[perm_inds])
            acc['DeepCSV_bDisc_vs_Jet3_eta'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='bjet',
                eta=deepcsv_sorted_jets[perm_inds].eta[:, 2], deepcsv_bdisc=deepcsv_sorted_jets[perm_inds].btagDeepB[:, 2], weight=evt_wts[perm_inds])

                # split DeepCSV value based on b, c, and light jets
            acc['DeepCSV_bDisc'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='bjet',
                deepcsv_bdisc=ak.flatten(deepcsv_sorted_jets[perm_inds].btagDeepB), weight=ak.flatten((ak.ones_like(jets.pt)*evt_wts)[perm_inds]))

        return acc        

    def fill_mc_hists(self, acc, sys, jetmult, leptype, btagregion, permarray, perm, jets, leptons, MTvals, evt_wts):
            ## apply alpha correction for 3Jets
        if (jetmult == '3Jets') and ('Alpha' in self.corrections):
            alpha_corr = self.corrections['Alpha'](172.5/perm['THad'].mass)
            perm['THad'] = perm['THad'].multiply(alpha_corr) # correct thad
            perm['TTbar'] = ak.flatten(perm['THad']+perm['TLep']) # correct ttbar


        pt_sorted_jets = jets[ak.argsort(jets.pt, ascending=False)]
        deepcsv_sorted_jets = jets[ak.argsort(jets.btagDeepB, ascending=False)]
        deepjet_sorted_jets = jets[ak.argsort(jets.btagDeepFlavB, ascending=False)]

            # split deepcsv sorted jets into b, c, and light jets
        deepcsv_sorted_bjets = deepcsv_sorted_jets[deepcsv_sorted_jets['hadronFlavour'] == 5]
        deepcsv_sorted_cjets = deepcsv_sorted_jets[deepcsv_sorted_jets['hadronFlavour'] == 4]
        deepcsv_sorted_ljets = deepcsv_sorted_jets[deepcsv_sorted_jets['hadronFlavour'] == 0]


        for permval in np.unique(permarray).tolist():
            perm_inds = np.where(permarray == permval)
            dataset_name = '%s_%s' % (self.sample_name, perm_cats[permval]) if permval != 0 else self.sample_name

            #set_trace()
                # jets are b jet
            acc['DeepCSV_bDisc_vs_Jet_pt'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='bjet',
                pt=ak.flatten(jets[perm_inds][jets[perm_inds]['hadronFlavour'] == 5].pt),
                deepcsv_bdisc=ak.flatten(jets[perm_inds][jets[perm_inds]['hadronFlavour'] == 5].btagDeepB),
                weight=ak.flatten((ak.ones_like(jets[jets['hadronFlavour'] == 5].pt)*evt_wts)[perm_inds]))
            acc['DeepCSV_bDisc_vs_Jet_eta'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='bjet',
                eta=ak.flatten(jets[perm_inds][jets[perm_inds]['hadronFlavour'] == 5].eta),
                deepcsv_bdisc=ak.flatten(jets[perm_inds][jets[perm_inds]['hadronFlavour'] == 5].btagDeepB),
                weight=ak.flatten((ak.ones_like(jets[jets['hadronFlavour'] == 5].eta)*evt_wts)[perm_inds]))
                # jets are c jet
            acc['DeepCSV_bDisc_vs_Jet_pt'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='cjet',
                pt=ak.flatten(jets[perm_inds][jets[perm_inds]['hadronFlavour'] == 4].pt),
                deepcsv_bdisc=ak.flatten(jets[perm_inds][jets[perm_inds]['hadronFlavour'] == 4].btagDeepB),
                weight=ak.flatten((ak.ones_like(jets[jets['hadronFlavour'] == 4].pt)*evt_wts)[perm_inds]))
            acc['DeepCSV_bDisc_vs_Jet_eta'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='cjet',
                eta=ak.flatten(jets[perm_inds][jets[perm_inds]['hadronFlavour'] == 4].eta),
                deepcsv_bdisc=ak.flatten(jets[perm_inds][jets[perm_inds]['hadronFlavour'] == 4].btagDeepB),
                weight=ak.flatten((ak.ones_like(jets[jets['hadronFlavour'] == 4].eta)*evt_wts)[perm_inds]))
                # jets are l jet
            acc['DeepCSV_bDisc_vs_Jet_pt'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='ljet',
                pt=ak.flatten(jets[perm_inds][jets[perm_inds]['hadronFlavour'] == 0].pt),
                deepcsv_bdisc=ak.flatten(jets[perm_inds][jets[perm_inds]['hadronFlavour'] == 0].btagDeepB),
                weight=ak.flatten((ak.ones_like(jets[jets['hadronFlavour'] == 0].pt)*evt_wts)[perm_inds]))
            acc['DeepCSV_bDisc_vs_Jet_eta'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='ljet',
                eta=ak.flatten(jets[perm_inds][jets[perm_inds]['hadronFlavour'] == 0].eta),
                deepcsv_bdisc=ak.flatten(jets[perm_inds][jets[perm_inds]['hadronFlavour'] == 0].btagDeepB),
                weight=ak.flatten((ak.ones_like(jets[jets['hadronFlavour'] == 0].eta)*evt_wts)[perm_inds]))

                # jet with highest DeepCSV value (jet 1) is b jet
            acc['DeepCSV_bDisc_vs_Jet1_pt'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='bjet',
                pt=deepcsv_sorted_jets[:, 0][perm_inds][deepcsv_sorted_jets[:, 0][perm_inds]['hadronFlavour'] == 5].pt,
                deepcsv_bdisc=deepcsv_sorted_jets[:, 0][perm_inds][deepcsv_sorted_jets[:, 0][perm_inds]['hadronFlavour'] == 5].btagDeepB,
                weight=evt_wts[perm_inds][deepcsv_sorted_jets[:, 0][perm_inds]['hadronFlavour'] == 5])
            acc['DeepCSV_bDisc_vs_Jet1_eta'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='bjet',
                eta=deepcsv_sorted_jets[:, 0][perm_inds][deepcsv_sorted_jets[:, 0][perm_inds]['hadronFlavour'] == 5].eta,
                deepcsv_bdisc=deepcsv_sorted_jets[:, 0][perm_inds][deepcsv_sorted_jets[:, 0][perm_inds]['hadronFlavour'] == 5].btagDeepB,
                weight=evt_wts[perm_inds][deepcsv_sorted_jets[:, 0][perm_inds]['hadronFlavour'] == 5])
                # jet with highest DeepCSV value (jet 1) is c jet
            acc['DeepCSV_bDisc_vs_Jet1_pt'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='cjet',
                pt=deepcsv_sorted_jets[:, 0][perm_inds][deepcsv_sorted_jets[:, 0][perm_inds]['hadronFlavour'] == 4].pt,
                deepcsv_bdisc=deepcsv_sorted_jets[:, 0][perm_inds][deepcsv_sorted_jets[:, 0][perm_inds]['hadronFlavour'] == 4].btagDeepB,
                weight=evt_wts[perm_inds][deepcsv_sorted_jets[:, 0][perm_inds]['hadronFlavour'] == 4])
            acc['DeepCSV_bDisc_vs_Jet1_eta'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='cjet',
                eta=deepcsv_sorted_jets[:, 0][perm_inds][deepcsv_sorted_jets[:, 0][perm_inds]['hadronFlavour'] == 4].eta,
                deepcsv_bdisc=deepcsv_sorted_jets[:, 0][perm_inds][deepcsv_sorted_jets[:, 0][perm_inds]['hadronFlavour'] == 4].btagDeepB,
                weight=evt_wts[perm_inds][deepcsv_sorted_jets[:, 0][perm_inds]['hadronFlavour'] == 4])
                # jet with highest DeepCSV value (jet 1) is l jet
            acc['DeepCSV_bDisc_vs_Jet1_pt'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='ljet',
                pt=deepcsv_sorted_jets[:, 0][perm_inds][deepcsv_sorted_jets[:, 0][perm_inds]['hadronFlavour'] == 0].pt,
                deepcsv_bdisc=deepcsv_sorted_jets[:, 0][perm_inds][deepcsv_sorted_jets[:, 0][perm_inds]['hadronFlavour'] == 0].btagDeepB,
                weight=evt_wts[perm_inds][deepcsv_sorted_jets[:, 0][perm_inds]['hadronFlavour'] == 0])
            acc['DeepCSV_bDisc_vs_Jet1_eta'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='ljet',
                eta=deepcsv_sorted_jets[:, 0][perm_inds][deepcsv_sorted_jets[:, 0][perm_inds]['hadronFlavour'] == 0].eta,
                deepcsv_bdisc=deepcsv_sorted_jets[:, 0][perm_inds][deepcsv_sorted_jets[:, 0][perm_inds]['hadronFlavour'] == 0].btagDeepB,
                weight=evt_wts[perm_inds][deepcsv_sorted_jets[:, 0][perm_inds]['hadronFlavour'] == 0])

                # jet with second highest DeepCSV value (jet 2) is b jet
            acc['DeepCSV_bDisc_vs_Jet2_pt'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='bjet',
                pt=deepcsv_sorted_jets[:, 1][perm_inds][deepcsv_sorted_jets[:, 1][perm_inds]['hadronFlavour'] == 5].pt,
                deepcsv_bdisc=deepcsv_sorted_jets[:, 1][perm_inds][deepcsv_sorted_jets[:, 1][perm_inds]['hadronFlavour'] == 5].btagDeepB,
                weight=evt_wts[perm_inds][deepcsv_sorted_jets[:, 1][perm_inds]['hadronFlavour'] == 5])
            acc['DeepCSV_bDisc_vs_Jet2_eta'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='bjet',
                eta=deepcsv_sorted_jets[:, 1][perm_inds][deepcsv_sorted_jets[:, 1][perm_inds]['hadronFlavour'] == 5].eta,
                deepcsv_bdisc=deepcsv_sorted_jets[:, 1][perm_inds][deepcsv_sorted_jets[:, 1][perm_inds]['hadronFlavour'] == 5].btagDeepB,
                weight=evt_wts[perm_inds][deepcsv_sorted_jets[:, 1][perm_inds]['hadronFlavour'] == 5])
                # jet with second highest DeepCSV value (jet 2) is c jet
            acc['DeepCSV_bDisc_vs_Jet2_pt'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='cjet',
                pt=deepcsv_sorted_jets[:, 1][perm_inds][deepcsv_sorted_jets[:, 1][perm_inds]['hadronFlavour'] == 4].pt,
                deepcsv_bdisc=deepcsv_sorted_jets[:, 1][perm_inds][deepcsv_sorted_jets[:, 1][perm_inds]['hadronFlavour'] == 4].btagDeepB,
                weight=evt_wts[perm_inds][deepcsv_sorted_jets[:, 1][perm_inds]['hadronFlavour'] == 4])
            acc['DeepCSV_bDisc_vs_Jet2_eta'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='cjet',
                eta=deepcsv_sorted_jets[:, 1][perm_inds][deepcsv_sorted_jets[:, 1][perm_inds]['hadronFlavour'] == 4].eta,
                deepcsv_bdisc=deepcsv_sorted_jets[:, 1][perm_inds][deepcsv_sorted_jets[:, 1][perm_inds]['hadronFlavour'] == 4].btagDeepB,
                weight=evt_wts[perm_inds][deepcsv_sorted_jets[:, 1][perm_inds]['hadronFlavour'] == 4])
                # jet with second highest DeepCSV value (jet 2) is l jet
            acc['DeepCSV_bDisc_vs_Jet2_pt'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='ljet',
                pt=deepcsv_sorted_jets[:, 1][perm_inds][deepcsv_sorted_jets[:, 1][perm_inds]['hadronFlavour'] == 0].pt,
                deepcsv_bdisc=deepcsv_sorted_jets[:, 1][perm_inds][deepcsv_sorted_jets[:, 1][perm_inds]['hadronFlavour'] == 0].btagDeepB,
                weight=evt_wts[perm_inds][deepcsv_sorted_jets[:, 1][perm_inds]['hadronFlavour'] == 0])
            acc['DeepCSV_bDisc_vs_Jet2_eta'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='ljet',
                eta=deepcsv_sorted_jets[:, 1][perm_inds][deepcsv_sorted_jets[:, 1][perm_inds]['hadronFlavour'] == 0].eta,
                deepcsv_bdisc=deepcsv_sorted_jets[:, 1][perm_inds][deepcsv_sorted_jets[:, 1][perm_inds]['hadronFlavour'] == 0].btagDeepB,
                weight=evt_wts[perm_inds][deepcsv_sorted_jets[:, 1][perm_inds]['hadronFlavour'] == 0])

                # jet with third highest DeepCSV value (jet 3) is b jet
            acc['DeepCSV_bDisc_vs_Jet3_pt'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='bjet',
                pt=deepcsv_sorted_jets[:, 2][perm_inds][deepcsv_sorted_jets[:, 2][perm_inds]['hadronFlavour'] == 5].pt,
                deepcsv_bdisc=deepcsv_sorted_jets[:, 2][perm_inds][deepcsv_sorted_jets[:, 2][perm_inds]['hadronFlavour'] == 5].btagDeepB,
                weight=evt_wts[perm_inds][deepcsv_sorted_jets[:, 2][perm_inds]['hadronFlavour'] == 5])
            acc['DeepCSV_bDisc_vs_Jet3_eta'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='bjet',
                eta=deepcsv_sorted_jets[:, 2][perm_inds][deepcsv_sorted_jets[:, 2][perm_inds]['hadronFlavour'] == 5].eta,
                deepcsv_bdisc=deepcsv_sorted_jets[:, 2][perm_inds][deepcsv_sorted_jets[:, 2][perm_inds]['hadronFlavour'] == 5].btagDeepB,
                weight=evt_wts[perm_inds][deepcsv_sorted_jets[:, 2][perm_inds]['hadronFlavour'] == 5])
                # jet with third highest DeepCSV value (jet 3) is c jet
            acc['DeepCSV_bDisc_vs_Jet3_pt'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='cjet',
                pt=deepcsv_sorted_jets[:, 2][perm_inds][deepcsv_sorted_jets[:, 2][perm_inds]['hadronFlavour'] == 4].pt,
                deepcsv_bdisc=deepcsv_sorted_jets[:, 2][perm_inds][deepcsv_sorted_jets[:, 2][perm_inds]['hadronFlavour'] == 4].btagDeepB,
                weight=evt_wts[perm_inds][deepcsv_sorted_jets[:, 2][perm_inds]['hadronFlavour'] == 4])
            acc['DeepCSV_bDisc_vs_Jet3_eta'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='cjet',
                eta=deepcsv_sorted_jets[:, 2][perm_inds][deepcsv_sorted_jets[:, 2][perm_inds]['hadronFlavour'] == 4].eta,
                deepcsv_bdisc=deepcsv_sorted_jets[:, 2][perm_inds][deepcsv_sorted_jets[:, 2][perm_inds]['hadronFlavour'] == 4].btagDeepB,
                weight=evt_wts[perm_inds][deepcsv_sorted_jets[:, 2][perm_inds]['hadronFlavour'] == 4])
                # jet with third highest DeepCSV value (jet 3) is l jet
            acc['DeepCSV_bDisc_vs_Jet3_pt'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='ljet',
                pt=deepcsv_sorted_jets[:, 2][perm_inds][deepcsv_sorted_jets[:, 2][perm_inds]['hadronFlavour'] == 0].pt,
                deepcsv_bdisc=deepcsv_sorted_jets[:, 2][perm_inds][deepcsv_sorted_jets[:, 2][perm_inds]['hadronFlavour'] == 0].btagDeepB,
                weight=evt_wts[perm_inds][deepcsv_sorted_jets[:, 2][perm_inds]['hadronFlavour'] == 0])
            acc['DeepCSV_bDisc_vs_Jet3_eta'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='ljet',
                eta=deepcsv_sorted_jets[:, 2][perm_inds][deepcsv_sorted_jets[:, 2][perm_inds]['hadronFlavour'] == 0].eta,
                deepcsv_bdisc=deepcsv_sorted_jets[:, 2][perm_inds][deepcsv_sorted_jets[:, 2][perm_inds]['hadronFlavour'] == 0].btagDeepB,
                weight=evt_wts[perm_inds][deepcsv_sorted_jets[:, 2][perm_inds]['hadronFlavour'] == 0])


                # split DeepCSV value based on b, c, and light jets
            acc['DeepCSV_bDisc'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='bjet', 
                deepcsv_bdisc=ak.flatten(deepcsv_sorted_bjets[perm_inds].btagDeepB), weight=ak.flatten((ak.ones_like(deepcsv_sorted_bjets.pt)*evt_wts)[perm_inds]))
            acc['DeepCSV_bDisc'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='cjet', 
                deepcsv_bdisc=ak.flatten(deepcsv_sorted_cjets[perm_inds].btagDeepB), weight=ak.flatten((ak.ones_like(deepcsv_sorted_cjets.pt)*evt_wts)[perm_inds]))
            acc['DeepCSV_bDisc'].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, hFlav='ljet', 
                deepcsv_bdisc=ak.flatten(deepcsv_sorted_ljets[perm_inds].btagDeepB), weight=ak.flatten((ak.ones_like(deepcsv_sorted_ljets.pt)*evt_wts)[perm_inds]))

        return acc        


    def postprocess(self, accumulator):
        return accumulator


proc_executor = processor.iterative_executor if to_debug else processor.futures_executor
proc_exec_args = {"schema": NanoAODSchema} if to_debug else {"schema": NanoAODSchema, "workers": 8}
output = processor.run_uproot_job(
    fileset,
    treename="Events",
    processor_instance=sb_regions_hFlavs(),
    executor=proc_executor,
    executor_args=proc_exec_args,
    #chunksize=5000,
    chunksize=100000,
)

save(output, args.outfname)
print('%s has been written' % args.outfname)

toc = time.time()
print("Total time: %.1f" % (toc - tic))
