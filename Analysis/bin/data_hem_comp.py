#!/usr/bin/env python

from coffea import hist, processor
from coffea.nanoevents import NanoAODSchema
from coffea.analysis_tools import PackedSelection
import awkward as ak
from coffea.nanoevents.methods import vector
ak.behavior.update(vector.behavior)

from coffea.util import save, load
from pdb import set_trace
import os
import python.ObjectSelection as objsel
import python.MCWeights as MCWeights
import numpy as np
import Utilities.prettyjson as prettyjson
import coffea.lumi_tools.lumi_tools as lumi_tools
import Utilities.make_variables as make_vars
from python.IDJet import btag_values as btag_values
import python.TTPermutator as ttpermutator
import python.IDJet as IDJet

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
base_jobid = os.environ['base_jobid']
analyzer = 'data_hem_comp'

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('fset', type=str, help='Fileset dictionary (in string form) to be used for the processor')
parser.add_argument('year', choices=['2018'], help='Specify which year to run over')
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
isData_ = samplename.startswith('data')
isSE_Data_ = samplename.startswith('data_SingleElectron')
isSM_Data_ = samplename.startswith('data_SingleMuon')
if isData_:
    lumiMask_path = os.path.join(proj_dir, 'inputs', 'data', base_jobid, 'LumiMasks', '%s_GoldenJson_%s.txt' % (args.year, base_jobid))
else:
    raise ValueError("Only data should be run with this analyzer.")

# convert input string of options dictionary to actual dictionary
odict = (args.opts).replace("\'", "\"")
opts_dict = prettyjson.loads(odict)

    ## set config options passed through argparse
import ast
to_debug = ast.literal_eval(opts_dict.get('debug', 'False'))
#apply_hem = ast.literal_eval(opts_dict.get('apply_hem', 'False')) and (not isData_)

## init tt probs for likelihoods
ttpermutator.year_to_run(year=args.year)

## load lumimask for data and corrections for event weights
jet_corrections = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'JetMETCorrections.coffea'))[args.year]
alpha_corrections = load(os.path.join(proj_dir, 'Corrections', base_jobid,'alpha_correction_%s.coffea' % jobid))[args.year]['E']['All_2D'] # E, All_2D determined by post application plots
corrections = {
    'JetCor' : jet_corrections,
    'Alpha' : alpha_corrections,
}

cfg_pars = prettyjson.loads(open(os.path.join(proj_dir, 'cfg_files', 'cfg_pars_%s.json' % jobid)).read())
    ## parameters for b-tagging
jet_pars = cfg_pars['Jets']
btaggers = ['DeepCSV']

wps_to_use = list(set([jet_pars['permutations']['tightb'],jet_pars['permutations']['looseb']]))
if not( len(wps_to_use) == 1):
    raise IOError("Only 1 unique btag working point supported now")
btag_wps = ['DeepCSVMedium']

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

## define run number that corresponds to before/after HEM issue https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/2000.html
hem_run_after = 319077 # runs >= 319077


# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class data_hem_comp(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.hem_axis = hist.Cat("hem", "Before/After HEM")
        self.jetmult_axis = hist.Cat("jmult", "nJets")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")
        self.lepcat_axis = hist.Cat("lepcat", "Lepton Category")
        self.btag_axis = hist.Cat("btag", "BTag Region")
        self.pt_axis = hist.Bin("pt", "p_{T} [GeV]", 200, 0, 1000)
        self.eta_axis = hist.Bin("eta", r"$\eta$", 200, -5., 5.)
        self.phi_axis = hist.Bin("phi", r"$\phi$", 160, -4, 4)
        #self.energy_axis = hist.Bin("energy", "E [GeV]", 200, 0, 1000)
        self.njets_axis = hist.Bin("njets", "n_{jets}", 20, 0, 20)
        self.lepIso_axis = hist.Bin("iso", "pfRelIso", 100, 0., 1.)
        self.mt_axis = hist.Bin("mt", "M_{T}", 200, 0., 1000.)
        self.mtop_axis = hist.Bin("mtop", "m(top) [GeV]", 300, 0, 300)
        self.mtt_axis = hist.Bin("mtt", "m($t\overline{t}$) [GeV]", 360, 200, 2000)
        self.probDisc_axis = hist.Bin("prob", "$\lambda_{C}$", 300, 0, 30)
        self.massDisc_axis = hist.Bin("massdisc", "$\lambda_{M}$", 300, 0, 30)
        self.nsDisc_axis = hist.Bin("nsdisc", "$\lambda_{NS}$", 300, 0, 30)
        self.nu_dist_axis = hist.Bin("nu_dist", r"$D_{\nu, min}$", 150, 0., 150.)
        self.ctstar_axis = hist.Bin("ctstar", "cos($\\theta^{*}$)", 200, -1., 1.)
        self.ctstar_abs_axis = hist.Bin("ctstar_abs", "|cos($\\theta^{*}$)|", 200, 0., 1.)

            ## make dictionary of hists
        histo_dict = {}
                ## make jet hists
        jet_hists = self.make_jet_hists()
        histo_dict.update(jet_hists)
                ## make lepton hists
        lep_hists = self.make_lep_hists()
        histo_dict.update(lep_hists)        
                ## make selection plots
        selection_hists = self.make_selection_hists()
        histo_dict.update(selection_hists)        

        self.sample_name = ''
        self.corrections = corrections

            ## make dict of cutflow for each systematic variation
        histo_dict['cutflow'] = processor.defaultdict_accumulator(int)
    
        self._accumulator = processor.dict_accumulator(histo_dict)

            ## define regions depending on data or MC, systematic or nosys
        self.regions = {
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

        if isSM_Data_:
            del self.regions['Electron']
        if isSE_Data_:
            del self.regions['Muon']
        for lepton in self.regions.keys():
            for leptype in self.regions[lepton].keys():
                for btagregion in self.regions[lepton][leptype].keys():
                    for jmult in self.regions[lepton][leptype][btagregion].keys():
                        self.regions[lepton][leptype][btagregion][jmult].update({'lumimask'})


    @property
    def accumulator(self):
        return self._accumulator


    def make_jet_hists(self):
        histo_dict = {}
        histo_dict['Jets_pt']    = hist.Hist("Events", self.dataset_axis, self.hem_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.pt_axis)
        histo_dict['Jets_eta']   = hist.Hist("Events", self.dataset_axis, self.hem_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.eta_axis)
        histo_dict['Jets_phi']   = hist.Hist("Events", self.dataset_axis, self.hem_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.phi_axis)
        histo_dict['Jets_njets'] = hist.Hist("Events", self.dataset_axis, self.hem_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.njets_axis)
        histo_dict['Jets_LeadJet_pt']    = hist.Hist("Events", self.dataset_axis, self.hem_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.pt_axis)
        histo_dict['Jets_LeadJet_eta']   = hist.Hist("Events", self.dataset_axis, self.hem_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.eta_axis)
        histo_dict['Jets_phi_vs_eta'] = hist.Hist("Events", self.dataset_axis, self.hem_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.phi_axis, self.eta_axis)

        return histo_dict

    def make_lep_hists(self):
        histo_dict = {}
        histo_dict['Lep_pt']    = hist.Hist("Events", self.dataset_axis, self.hem_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.pt_axis)
        histo_dict['Lep_eta']   = hist.Hist("Events", self.dataset_axis, self.hem_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.eta_axis)
        histo_dict['Lep_iso']   = hist.Hist("Events", self.dataset_axis, self.hem_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.lepIso_axis)
        histo_dict['Lep_phi']   = hist.Hist("Events", self.dataset_axis, self.hem_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.phi_axis)
        histo_dict['Lep_phi_vs_eta'] = hist.Hist("Events", self.dataset_axis, self.hem_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.phi_axis, self.eta_axis)
    #    histo_dict['Lep_etaSC'] = hist.Hist("Events", self.dataset_axis, self.hem_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.eta_axis)
    #    histo_dict['Lep_energy']= hist.Hist("Events", self.dataset_axis, self.hem_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.energy_axis)

        return histo_dict

    def make_selection_hists(self):
        histo_dict = {}
        histo_dict['mtt']      = hist.Hist("Events", self.dataset_axis, self.hem_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mtt_axis)
        histo_dict['mthad']    = hist.Hist("Events", self.dataset_axis, self.hem_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mtop_axis)
        histo_dict['pt_thad']  = hist.Hist("Events", self.dataset_axis, self.hem_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.pt_axis)
        histo_dict['pt_tlep']  = hist.Hist("Events", self.dataset_axis, self.hem_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.pt_axis)
        histo_dict['pt_tt']    = hist.Hist("Events", self.dataset_axis, self.hem_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.pt_axis)
        histo_dict['eta_thad'] = hist.Hist("Events", self.dataset_axis, self.hem_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.eta_axis)
        histo_dict['eta_tlep'] = hist.Hist("Events", self.dataset_axis, self.hem_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.eta_axis)
        histo_dict['eta_tt']   = hist.Hist("Events", self.dataset_axis, self.hem_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.eta_axis)

        histo_dict['tlep_ctstar']     = hist.Hist("Events", self.dataset_axis, self.hem_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.ctstar_axis)
        histo_dict['tlep_ctstar_abs'] = hist.Hist("Events", self.dataset_axis, self.hem_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.ctstar_abs_axis)

        histo_dict['full_disc'] = hist.Hist("Events", self.dataset_axis, self.hem_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.probDisc_axis)
        histo_dict['mass_disc'] = hist.Hist("Events", self.dataset_axis, self.hem_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.massDisc_axis)
        histo_dict['ns_disc']   = hist.Hist("Events", self.dataset_axis, self.hem_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.nsDisc_axis)
        histo_dict['ns_dist']   = hist.Hist("Events", self.dataset_axis, self.hem_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.nu_dist_axis)

        histo_dict['MT'] = hist.Hist("Events", self.dataset_axis, self.hem_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mt_axis)

        histo_dict['MET_pt'] = hist.Hist("Events", self.dataset_axis, self.hem_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.pt_axis)
        histo_dict['MET_phi']= hist.Hist("Events", self.dataset_axis, self.hem_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.phi_axis)

        histo_dict['mtt_vs_tlep_ctstar_abs'] = hist.Hist("Events", self.dataset_axis, self.hem_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mtt_axis, self.ctstar_abs_axis)

        return histo_dict


    def process(self, events):
        np.random.seed(10) # sets seed so values from random distributions are reproducible (JER corrections)
        output = self.accumulator.identity()

        self.sample_name = events.metadata['dataset']

            ## make event weights
                # data or MC distinction made internally
        mu_evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections)
        el_evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections)

            ## initialize selections
        selection = PackedSelection()

            # get all passing leptons
        lep_and_filter_pass = objsel.select_leptons(events, year=args.year)#, cutflow=output['cutflow'])
        selection.add('lep_and_filter_pass', lep_and_filter_pass) # add passing leptons requirement to all systematics

            ## build corrected jets and MET
        events['Jet'], events['MET'] = IDJet.process_jets(events, args.year, self.corrections['JetCor'])

        runs = events.run
        lumis = events.luminosityBlock
        Golden_Json_LumiMask = lumi_tools.LumiMask(lumiMask_path)
        LumiMask = Golden_Json_LumiMask.__call__(runs, lumis) ## returns array of valid events
        selection.add('lumimask', LumiMask)
              ## object selection and add different selections
        if isSM_Data_:
                    ## muons
            selection.add('tight_MU', ak.sum(events['Muon']['TIGHTMU'], axis=1) == 1) # one muon passing TIGHT criteria
            selection.add('loose_MU', ak.sum(events['Muon']['LOOSEMU'], axis=1) == 1) # one muon passing LOOSE criteria
        if isSE_Data_:
                    ## electrons
            selection.add('tight_EL', ak.sum(events['Electron']['TIGHTEL'], axis=1) == 1) # one electron passing TIGHT criteria
            selection.add('loose_EL', ak.sum(events['Electron']['LOOSEEL'], axis=1) == 1) # one electron passing LOOSE criteria

        output['cutflow']['lep_and_filter_pass'] += ak.sum(lep_and_filter_pass)
            # jet selection
        passing_jets = objsel.jets_selection(events, year=args.year, cutflow=output['cutflow'])
        output['cutflow']['nEvts passing jet and lepton obj selection'] += ak.sum(passing_jets & lep_and_filter_pass)
        selection.add('passing_jets', passing_jets)
        selection.add('jets_3',  ak.num(events['SelectedJets']) == 3)
        selection.add('jets_4p',  ak.num(events['SelectedJets']) > 3) # only for getting btag weights
        selection.add('DeepCSV_pass', ak.sum(events['SelectedJets'][btag_wps[0]], axis=1) >= 2)

            # sort jets by btag value
        events['SelectedJets'] = events['SelectedJets'][ak.argsort(events['SelectedJets']['btagDeepB'], ascending=False)] if btaggers[0] == 'DeepCSV' else events['SelectedJets'][ak.argsort(events['SelectedJets']['btagDeepFlavB'], ascending=False)]

            # btag fail sideband
        deepcsv_sorted = events['SelectedJets'][ak.argsort(events['SelectedJets']['btagDeepB'], ascending=False)]['btagDeepB']
        valid_counts_inds = ak.where(ak.num(events['SelectedJets']) > 1)[0]
        deepcsv_fail = np.zeros(len(events)).astype(bool)
        deepcsv_fail[valid_counts_inds] = (deepcsv_sorted[valid_counts_inds][:, 0] < btag_values[args.year]['btagDeepB']['DeepCSV'+wps_to_use[0]]) & (deepcsv_sorted[valid_counts_inds][:, 1] < btag_values[args.year]['btagDeepB']['DeepCSV'+wps_to_use[0]])
        selection.add('DeepCSV_fail', deepcsv_fail) # highest and second highest DeepCSV values don't pass tight and loose WPs

        ## fill hists for each region
        for lepton in self.regions.keys():
            evt_weights = mu_evt_weights if lepton == 'Muon' else el_evt_weights
            for leptype in self.regions[lepton].keys():
                for btagregion in self.regions[lepton][leptype].keys():
                    for jmult in self.regions[lepton][leptype][btagregion].keys():
                        cut = selection.all(*self.regions[lepton][leptype][btagregion][jmult])
                        #set_trace()

                        output['cutflow']['nEvts %s' % ', '.join([lepton, leptype, btagregion, jmult])] += cut.sum()

                        if to_debug: print(lepton, leptype, btagregion, jmult)
                        if cut.sum() > 0:
                            ltype = 'MU' if lepton == 'Muon' else 'EL'
                            if 'loose_or_tight_%s' % ltype in self.regions[lepton][leptype][btagregion][jmult]:
                                leptons = events[lepton][cut][((events[lepton][cut]['TIGHT%s' % ltype] == True) | (events[lepton][cut]['LOOSE%s' % ltype] == True))]
                            elif 'tight_%s' % ltype in self.regions[lepton][leptype][btagregion][jmult]:
                                leptons = events[lepton][cut][(events[lepton][cut]['TIGHT%s' % ltype] == True)]
                            elif 'loose_%s' % ltype in self.regions[lepton][leptype][btagregion][jmult]:
                                leptons = events[lepton][cut][(events[lepton][cut]['LOOSE%s' % ltype] == True)]
                            else:
                                raise ValueError("Not sure what lepton type to choose for event")

                                # get jets and MET
                            jets, met = events['SelectedJets'][cut], events['SelectedMET'][cut]

                                # find best permutations
                            best_perms = ttpermutator.find_best_permutations(jets=jets, leptons=leptons, MET=met, btagWP=btag_wps[0], btag_req=False if btagregion == 'btagFail' else True)
                            valid_perms = ak.num(best_perms['TTbar'].pt) > 0
                            output['cutflow']['nEvts %s: valid perms' % ', '.join([lepton, leptype, btagregion, jmult])] += ak.sum(valid_perms)

                            bp_status = np.zeros(cut.size, dtype=int) # 0 == '' (no gen matching), 1 == 'right', 2 == 'matchable', 3 == 'unmatchable', 4 == 'sl_tau', 5 == 'noslep'

                                ## create MT regions
                            MT = make_vars.MT(leptons, met)
                            MTHigh = ak.flatten(MT[valid_perms] >= MTcut)
                            output['cutflow']['nEvts %s: pass MT cut' % ', '.join([lepton, leptype, btagregion, jmult])] += ak.sum(MTHigh)

                                # fill hists for before/after HEM issue
                            wts = evt_weights.weight()[cut][valid_perms][MTHigh]
                            valid_runs = runs[cut][valid_perms][MTHigh]
                            runs_before = ak.where(valid_runs < hem_run_after)
                            runs_after = ak.where(valid_runs >= hem_run_after)
                            if len(runs_before[0]) > 0:
                                output = self.fill_hists(acc=output, hem_reg='Before', jetmult=jmult, leptype=lepton, lepcat=leptype, btagregion=btagregion, permarray=bp_status[cut][valid_perms][MTHigh][runs_before],
                                  perm=best_perms[valid_perms][MTHigh][runs_before], jets=jets[valid_perms][MTHigh][runs_before], leptons=leptons[valid_perms][MTHigh][runs_before], MTvals=MT[valid_perms][MTHigh][runs_before], evt_wts=wts[runs_before])
                            if len(runs_after[0]) > 0:
                                output = self.fill_hists(acc=output, hem_reg='After', jetmult=jmult, leptype=lepton, lepcat=leptype, btagregion=btagregion, permarray=bp_status[cut][valid_perms][MTHigh][runs_after],
                                  perm=best_perms[valid_perms][MTHigh][runs_after], jets=jets[valid_perms][MTHigh][runs_after], leptons=leptons[valid_perms][MTHigh][runs_after], MTvals=MT[valid_perms][MTHigh][runs_after], evt_wts=wts[runs_after])

        return output

    def fill_hists(self, acc, hem_reg, jetmult, leptype, lepcat, btagregion, permarray, perm, jets, leptons, MTvals, evt_wts):
            ## apply alpha correction for 3Jets
        if (jetmult == '3Jets') and ('Alpha' in self.corrections):
            alpha_corr = self.corrections['Alpha'](172.5/perm['THad'].mass)
            perm['THad'] = perm['THad'].multiply(alpha_corr) # correct thad
            perm['TTbar'] = ak.flatten(perm['THad']+perm['TLep']) # correct ttbar

        thad_ctstar, tlep_ctstar = make_vars.ctstar(perm['THad'], perm['TLep'])
        thad_ctstar, tlep_ctstar = ak.flatten(thad_ctstar, axis=None), ak.flatten(tlep_ctstar, axis=None)

        pt_sorted_jets = jets[ak.argsort(jets.pt, ascending=False)]

        for permval in np.unique(permarray).tolist():
            perm_inds = np.where(permarray == permval)
            dataset_name = '%s_%s' % (self.sample_name, perm_cats[permval]) if permval != 0 else self.sample_name
            acc['mtt'].fill(     dataset=dataset_name, hem=hem_reg, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtt=ak.flatten(perm['TTbar'].mass)[perm_inds], weight=evt_wts[perm_inds])
            acc['mthad'].fill(   dataset=dataset_name, hem=hem_reg, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtop=ak.flatten(perm['THad'].mass)[perm_inds], weight=evt_wts[perm_inds])
            acc['pt_thad'].fill( dataset=dataset_name, hem=hem_reg, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=ak.flatten(perm['THad'].pt)[perm_inds], weight=evt_wts[perm_inds])
            acc['pt_tlep'].fill( dataset=dataset_name, hem=hem_reg, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=ak.flatten(perm['TLep'].pt)[perm_inds], weight=evt_wts[perm_inds])
            acc['pt_tt'].fill(   dataset=dataset_name, hem=hem_reg, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=ak.flatten(perm['TTbar'].pt)[perm_inds], weight=evt_wts[perm_inds])
            acc['eta_thad'].fill(dataset=dataset_name, hem=hem_reg, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=ak.flatten(perm['THad'].eta)[perm_inds], weight=evt_wts[perm_inds])
            acc['eta_tlep'].fill(dataset=dataset_name, hem=hem_reg, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=ak.flatten(perm['TLep'].eta)[perm_inds], weight=evt_wts[perm_inds])
            acc['eta_tt'].fill(  dataset=dataset_name, hem=hem_reg, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=ak.flatten(perm['TTbar'].eta)[perm_inds], weight=evt_wts[perm_inds])

            acc['tlep_ctstar'].fill(    dataset=dataset_name, hem=hem_reg, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, ctstar=tlep_ctstar[perm_inds], weight=evt_wts[perm_inds])
            acc['tlep_ctstar_abs'].fill(dataset=dataset_name, hem=hem_reg, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, ctstar_abs=np.abs(tlep_ctstar[perm_inds]), weight=evt_wts[perm_inds])

            acc['full_disc'].fill(dataset=dataset_name, hem=hem_reg, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, prob=ak.flatten(perm['Prob'])[perm_inds], weight=evt_wts[perm_inds])
            acc['mass_disc'].fill(dataset=dataset_name, hem=hem_reg, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, massdisc=ak.flatten(perm['MassDiscr'])[perm_inds], weight=evt_wts[perm_inds])
            acc['ns_disc'].fill(  dataset=dataset_name, hem=hem_reg, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, nsdisc=ak.flatten(perm['NuDiscr'])[perm_inds], weight=evt_wts[perm_inds])
            acc['ns_dist'].fill(  dataset=dataset_name, hem=hem_reg, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, nu_dist=ak.flatten(np.sqrt(perm['Nu'].chi2))[perm_inds], weight=evt_wts[perm_inds])

            acc['MET_pt'].fill( dataset=dataset_name, hem=hem_reg, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=ak.flatten(perm['MET'].pt)[perm_inds], weight=evt_wts[perm_inds])
            acc['MET_phi'].fill(dataset=dataset_name, hem=hem_reg, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, phi=ak.flatten(perm['MET'].phi)[perm_inds], weight=evt_wts[perm_inds])

            acc['mtt_vs_tlep_ctstar_abs'].fill(dataset=dataset_name, hem=hem_reg,  jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion,
                mtt=ak.flatten(perm['TTbar'].mass)[perm_inds], ctstar_abs=np.abs(tlep_ctstar[perm_inds]), weight=evt_wts[perm_inds])

            acc['Jets_pt'].fill(   dataset=dataset_name, hem=hem_reg, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=ak.flatten(jets.pt[perm_inds]), weight=ak.flatten((ak.ones_like(jets.pt)*evt_wts)[perm_inds]))
            acc['Jets_eta'].fill(  dataset=dataset_name, hem=hem_reg, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=ak.flatten(jets.eta[perm_inds]), weight=ak.flatten((ak.ones_like(jets.eta)*evt_wts)[perm_inds]))
            acc['Jets_phi'].fill(  dataset=dataset_name, hem=hem_reg, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, phi=ak.flatten(jets.phi[perm_inds]), weight=ak.flatten((ak.ones_like(jets.phi)*evt_wts)[perm_inds]))
            acc['Jets_njets'].fill(dataset=dataset_name, hem=hem_reg, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, njets=ak.num(jets)[perm_inds], weight=evt_wts[perm_inds])
            acc['Jets_phi_vs_eta'].fill(dataset=dataset_name, hem=hem_reg, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion,
                phi=ak.flatten(jets.phi[perm_inds]), eta=ak.flatten(jets.eta[perm_inds]), weight=ak.flatten((ak.ones_like(jets.phi)*evt_wts)[perm_inds]))

            acc['Jets_LeadJet_pt'].fill( dataset=dataset_name, hem=hem_reg, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=pt_sorted_jets[perm_inds].pt[:, 0], weight=evt_wts[perm_inds])
            acc['Jets_LeadJet_eta'].fill(dataset=dataset_name, hem=hem_reg, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=pt_sorted_jets[perm_inds].eta[:, 0], weight=evt_wts[perm_inds])

            acc['Lep_pt'].fill(    dataset=dataset_name, hem=hem_reg, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=ak.flatten(leptons.pt[perm_inds]), weight=evt_wts[perm_inds])
            acc['Lep_eta'].fill(   dataset=dataset_name, hem=hem_reg, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=ak.flatten(leptons.eta[perm_inds]), weight=evt_wts[perm_inds])
            acc['Lep_phi'].fill(   dataset=dataset_name, hem=hem_reg, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, phi=ak.flatten(leptons.phi[perm_inds]), weight=evt_wts[perm_inds])
            acc['Lep_iso'].fill(   dataset=dataset_name, hem=hem_reg, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion,
                iso=ak.flatten(leptons['pfRelIso04_all'][perm_inds]) if leptype == 'Muon' else ak.flatten(leptons['pfRelIso03_all'][perm_inds]), weight=evt_wts[perm_inds])
            acc['Lep_phi_vs_eta'].fill(dataset=dataset_name, hem=hem_reg, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion,
                phi=ak.flatten(leptons.phi[perm_inds]), eta=ak.flatten(leptons.eta[perm_inds]), weight=evt_wts[perm_inds])

            acc['MT'].fill(dataset=dataset_name, hem=hem_reg, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mt=ak.flatten(MTvals)[perm_inds], weight=evt_wts[perm_inds])

        return acc        

    def postprocess(self, accumulator):
        return accumulator


proc_executor = processor.iterative_executor if to_debug else processor.futures_executor
proc_exec_args = {"schema": NanoAODSchema} if to_debug else {"schema": NanoAODSchema, "workers": 8}
output = processor.run_uproot_job(
    fileset,
    treename="Events",
    processor_instance=data_hem_comp(),
    executor=proc_executor,
    executor_args=proc_exec_args,
    #chunksize=10000 if to_debug else 100000,
    chunksize=100000,
)

save(output, args.outfname)
print('%s has been written' % args.outfname)
