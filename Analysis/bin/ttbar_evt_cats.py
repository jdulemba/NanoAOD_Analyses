#!/usr/bin/env python

import time
tic = time.time()

from coffea import hist, processor
from coffea.nanoevents import NanoAODSchema
from coffea.analysis_tools import PackedSelection
import awkward as ak
from coffea.nanoevents.methods import vector
ak.behavior.update(vector.behavior)

from coffea.util import save, load
from pdb import set_trace
import os, sys
import python.ObjectSelection as objsel
import Utilities.plot_tools as plt_tools
import python.MCWeights as MCWeights
import Utilities.make_variables as make_vars
import numpy as np
import Utilities.prettyjson as prettyjson
import coffea.lumi_tools.lumi_tools as lumi_tools
import python.GenParticleSelector as genpsel
import python.TTGenMatcher as ttmatcher
import python.IDJet as IDJet

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
base_jobid = os.environ['base_jobid']
analyzer = 'ttbar_evt_cats'

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('fset', type=str, help='Fileset dictionary (in string form) to be used for the processor')
parser.add_argument('year', choices=['2016', '2017', '2018'] if base_jobid == 'NanoAODv6' else ['2016APV', '2016', '2017', '2018'], help='Specify which year to run over')
parser.add_argument('outfname', type=str, help='Specify output filename, including directory and file extension')
parser.add_argument('opts', type=str, help='Fileset dictionary (in string form) to be used for the processor')
args = parser.parse_args()

# convert input string of fileset dictionary to actual dictionary
fdict = (args.fset).replace("\'", "\"")
fileset = prettyjson.loads(fdict)
if len(fileset.keys()) > 1:
    raise ValueError("Only one topology run at a time in order to determine which corrections and systematics to run")
samplename = list(fileset.keys())[0]

    ## specify ttJets samples
allowed_samples = ['ttJets_PS'] if ((args.year == '2016') and (base_jobid == 'NanoAODv6')) else ['ttJetsSL']
isAllowed = np.array([(key in allowed_samples) for key in fileset.keys()]).all()
if not isAllowed:
    raise ValueError("Not a valid dataset to run on")

# convert input string of options dictionary to actual dictionary
odict = (args.opts).replace("\'", "\"")
opts_dict = prettyjson.loads(odict)

    ## set config options passed through argparse
import ast
to_debug = ast.literal_eval(opts_dict.get('debug', 'False'))
apply_hem = ast.literal_eval(opts_dict.get('apply_hem', 'True'))


## load corrections for event weights
pu_correction = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'MC_PU_Weights.coffea'))[args.year]
lepSF_correction = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'leptonSFs.coffea'))[args.year]
jet_corrections = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'JetMETCorrections.coffea'))[args.year]
nnlo_var = 'mtt_vs_thad_ctstar_Interp'
nnlo_reweighting = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'NNLO_to_Tune_Orig_Interp_Ratios_%s.coffea' % base_jobid))[args.year]
corrections = {
    'Pileup' : pu_correction,
    'Prefire' : True,
    'LeptonSF' : lepSF_correction,
    'BTagSF' : True,
    #'BTagSF' : False,
    'JetCor' : jet_corrections,
    'NNLO_Rewt' : {'Var' : nnlo_var, 'Correction' : nnlo_reweighting[nnlo_var]},
}

cfg_pars = prettyjson.loads(open(os.path.join(proj_dir, 'cfg_files', 'cfg_pars_%s.json' % jobid)).read())
    ## parameters for b-tagging
jet_pars = cfg_pars['Jets']
btaggers = [jet_pars['btagger']]

wps_to_use = list(set([jet_pars['permutations']['tightb'],jet_pars['permutations']['looseb']]))
if not( len(wps_to_use) == 1):
    raise IOError("Only 1 unique btag working point supported now")
btag_wps = ['DeepCSVMedium']

if corrections['BTagSF'] == True:
    sf_file = os.path.join(proj_dir, 'Corrections', jobid, jet_pars['btagging']['btagSF_file'])
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


# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class ttbar_evt_cats(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")
        self.jetmult_axis = hist.Cat("jmult", "nJets")
        self.object_axis = hist.Cat("objtype", "Gen Object")
        self.cat_axis = hist.Cat("cat", "Evt Cat")
        self.pt_axis = hist.Bin("pt", "p_{T} [GeV]", 200, 0, 1000)
        self.eta_axis = hist.Bin("eta", r"$\eta$", 100, -5, 5)
        self.mtt_axis = hist.Bin("mtt", "Gen m($t\overline{t}$) [GeV]", 360, 200, 2000)

            ## make dictionary of hists
        histo_dict = {}
                ## make jet hists
        gen_hists = self.make_hists()
        histo_dict.update(gen_hists)

        histo_dict['cutflow'] = processor.defaultdict_accumulator(int)

        self._accumulator = processor.dict_accumulator(histo_dict)
        self.sample_name = ''
        self.corrections = corrections

        self.regions = {
            'Muon' : {
                'Tight' : {
                    '3Jets'  : {'lep_and_filter_pass', 'passing_jets', 'jets_3' , 'tight_MU', 'DeepCSV_pass', 'semilep'},
                    #'4PJets'  : {'lep_and_filter_pass', 'passing_jets', 'jets_4p' , 'tight_MU', 'DeepCSV_pass', 'semilep'},
                },
            },
            'Electron' : {
                'Tight' : {
                    '3Jets'  : {'lep_and_filter_pass', 'passing_jets', 'jets_3' , 'tight_EL', 'DeepCSV_pass', 'semilep'},
                    #'4PJets'  : {'lep_and_filter_pass', 'passing_jets', 'jets_4p' , 'tight_EL', 'DeepCSV_pass', 'semilep'},
                },
            },
        }
    
    @property
    def accumulator(self):
        return self._accumulator


    def make_hists(self):
        histo_dict = {}
        histo_dict['mtt']      = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.cat_axis, self.mtt_axis)
        histo_dict['pass_mtt'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.cat_axis, self.object_axis, self.mtt_axis)
        histo_dict['fail_mtt'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.cat_axis, self.object_axis, self.mtt_axis)

        return histo_dict



    def process(self, events):
        np.random.seed(10) # sets seed so values from random distributions are reproducible (JER corrections)
        output = self.accumulator.identity()

        self.sample_name = events.metadata['dataset']

            ## make event weights
                # data or MC distinction made internally
        mu_evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections, isTTbar=True)
        el_evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections, isTTbar=True)

            ## initialize selections and regions
        selection = PackedSelection()

            # get all passing leptons
        lep_and_filter_pass = objsel.select_leptons(events, year=args.year, hem_15_16=apply_hem)
        output['cutflow']['lep_and_filter_pass'] += ak.sum(lep_and_filter_pass)
        selection.add('lep_and_filter_pass', lep_and_filter_pass)

            ## build corrected jets and MET
        events['Jet'], events['MET'] = IDJet.process_jets(events, args.year, self.corrections['JetCor'])

        ## add different selections
                ## muons
        tight_mu_sel = ak.sum(events['Muon']['TIGHTMU'], axis=1) == 1
        selection.add('tight_MU', tight_mu_sel) # one muon passing TIGHT criteria
                ## electrons
        tight_el_sel = ak.sum(events['Electron']['TIGHTEL'], axis=1) == 1
        selection.add('tight_EL', tight_el_sel) # one electron passing TIGHT criteria

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


            # jet selection
        passing_jets = objsel.jets_selection(events, year=args.year, cutflow=output['cutflow'], hem_15_16=apply_hem)
        output['cutflow']['nEvts passing jet and lepton obj selection'] += ak.sum(passing_jets & lep_and_filter_pass)
        selection.add('passing_jets', passing_jets)
        selection.add('jets_3',  ak.num(events['SelectedJets']) == 3)
        selection.add('jets_4p',  ak.num(events['SelectedJets']) > 3) # only for getting btag weights
        selection.add('DeepCSV_pass', ak.sum(events['SelectedJets'][btag_wps[0]], axis=1) >= 2)

            # sort jets by btag value
        events['SelectedJets'] = events['SelectedJets'][ak.argsort(events['SelectedJets']['btagDeepB'], ascending=False)] if btaggers[0] == 'DeepCSV' else events['SelectedJets'][ak.argsort(events['SelectedJets']['btagDeepFlavB'], ascending=False)]

            ## apply btagging SFs to MC
        if (corrections['BTagSF'] == True):
            deepcsv_cen   = np.ones(len(events))
            #deepcsv_bc_up = np.ones(len(events))
            #deepcsv_bc_dw = np.ones(len(events))
            #deepcsv_l_up = np.ones(len(events))
            #deepcsv_l_dw = np.ones(len(events))

            threeJets_cut = selection.require(lep_and_filter_pass=True, passing_jets=True, jets_3=True)
            deepcsv_3j_wts = self.corrections['BTag_Constructors']['DeepCSV']['3Jets'].get_scale_factor(jets=events['SelectedJets'][threeJets_cut], passing_cut='DeepCSV'+wps_to_use[0])
            deepcsv_cen[threeJets_cut] = ak.prod(deepcsv_3j_wts['central'], axis=1)
            #deepcsv_bc_up[threeJets_cut] = ak.prod(deepcsv_3j_wts['bc_up'], axis=1)
            #deepcsv_bc_dw[threeJets_cut] = ak.prod(deepcsv_3j_wts['bc_down'], axis=1)
            #deepcsv_l_up[threeJets_cut] = ak.prod(deepcsv_3j_wts['udsg_up'], axis=1)
            #deepcsv_l_dw[threeJets_cut] = ak.prod(deepcsv_3j_wts['udsg_down'], axis=1)

            fourplusJets_cut = selection.require(lep_and_filter_pass=True, passing_jets=True, jets_4p=True)
            deepcsv_4pj_wts = self.corrections['BTag_Constructors']['DeepCSV']['4PJets'].get_scale_factor(jets=events['SelectedJets'][fourplusJets_cut], passing_cut='DeepCSV'+wps_to_use[0])
            deepcsv_cen[fourplusJets_cut] = ak.prod(deepcsv_4pj_wts['central'], axis=1)
            #deepcsv_bc_up[fourplusJets_cut] = ak.prod(deepcsv_4pj_wts['bc_up'], axis=1)
            #deepcsv_bc_dw[fourplusJets_cut] = ak.prod(deepcsv_4pj_wts['bc_down'], axis=1)
            #deepcsv_l_up[fourplusJets_cut] = ak.prod(deepcsv_4pj_wts['udsg_up'], axis=1)
            #deepcsv_l_dw[fourplusJets_cut] = ak.prod(deepcsv_4pj_wts['udsg_down'], axis=1)
                # make dict of btag weights
            btag_weights = {
                'DeepCSV_CEN' : deepcsv_cen,
                #'DeepCSV_bc_UP' : deepcsv_bc_up,
                #'DeepCSV_bc_DW' : deepcsv_bc_dw,
                #'DeepCSV_l_UP' : deepcsv_l_up,
                #'DeepCSV_l_DW' : deepcsv_l_dw,
            }
        else:
            btag_weights = {
                'DeepCSV_CEN'   : np.ones(len(events)),
                #'DeepCSV_bc_UP' : np.ones(len(events)),
                #'DeepCSV_bc_DW' : np.ones(len(events)),
                #'DeepCSV_l_UP'  : np.ones(len(events)),
                #'DeepCSV_l_DW'  : np.ones(len(events)),
            }

        # find gen level particles for ttbar system
        genpsel.select(events, mode='NORMAL')
        selection.add('semilep', ak.num(events['SL']) > 0)
        if 'NNLO_Rewt' in self.corrections.keys():
            nnlo_wts = MCWeights.get_nnlo_weights(self.corrections['NNLO_Rewt'], events)
            mu_evt_weights.add('%s_reweighting' % corrections['NNLO_Rewt']['Var'], nnlo_wts)
            el_evt_weights.add('%s_reweighting' % corrections['NNLO_Rewt']['Var'], nnlo_wts)


        ## fill hists for each region
        for lepton in self.regions.keys():
            evt_weights = mu_evt_weights if lepton == 'Muon' else el_evt_weights
            for leptype in self.regions[lepton].keys():
                for jmult in self.regions[lepton][leptype].keys():

                    if to_debug: print(lepton, leptype, 'btagPass', jmult)
                    cut = selection.all(*self.regions[lepton][leptype][jmult])
                    output['cutflow']['nEvts %s' % ', '.join([lepton, leptype, 'btagPass', jmult])] += cut.sum()
                    if cut.sum() == 0: continue

                    ltype = 'MU' if lepton == 'Muon' else 'EL'
                    if 'loose_or_tight_%s' % ltype in self.regions[lepton][leptype][jmult]:
                        leptons = events[lepton][cut][((events[lepton][cut]['TIGHT%s' % ltype] == True) | (events[lepton][cut]['LOOSE%s' % ltype] == True))]
                    elif 'tight_%s' % ltype in self.regions[lepton][leptype][jmult]:
                        leptons = events[lepton][cut][(events[lepton][cut]['TIGHT%s' % ltype] == True)]
                    elif 'loose_%s' % ltype in self.regions[lepton][leptype][jmult]:
                        leptons = events[lepton][cut][(events[lepton][cut]['LOOSE%s' % ltype] == True)]
                    else:
                        raise ValueError("Not sure what lepton type to choose for event")

                        # get jets and MET
                    jets, met = events['SelectedJets'][cut], events['SelectedMET'][cut]

                        ## create MT regions
                    MT = make_vars.MT(leptons, met)
                    MTHigh = ak.flatten(MT >= MTcut)
                    output['cutflow']['nEvts %s: pass MT cut' % ', '.join([lepton, leptype, 'btagPass', jmult])] += ak.sum(MTHigh)

                        # find matched permutations
                    mp = ttmatcher.best_match(gen_hyp=events['SL'][cut], jets=jets, leptons=leptons, met=met)

                    wts = (evt_weights.weight()*btag_weights['%s_CEN' % btaggers[0]])[cut][MTHigh]
                    valid_gen_objs = events['SL'][cut][MTHigh]
                    valid_mp = mp[MTHigh]

                    ## fill hists of gen mttbar for all events
                    output['mtt'].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, cat='All', mtt=ak.flatten(valid_gen_objs['TTbar'].mass, axis=None), weight=wts)
                        # for each gen parton find ones that pass and fail jet pt and eta cuts
                    for gen_obj in ['BHad', 'BLep', 'WJa', 'WJb']:
                        pt_eta_mask = ak.flatten((valid_gen_objs[gen_obj].pt >= jet_pars['ptmin']) & (np.abs(valid_gen_objs[gen_obj].eta) <= jet_pars['etamax']))
                        output['pass_mtt'].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, cat='All', objtype=gen_obj, mtt=ak.flatten(valid_gen_objs[pt_eta_mask]['TTbar'].mass, axis=None), weight=wts[pt_eta_mask])
                        output['fail_mtt'].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, cat='All', objtype=gen_obj, mtt=ak.flatten(valid_gen_objs[~pt_eta_mask]['TTbar'].mass, axis=None), weight=wts[~pt_eta_mask])

                    ## fill hists of gen mttbar for each lost jet/partially merged event category
                    evt_cats = ['Merged_BHadBLep', 'Merged_BHadWJa', 'Merged_BHadWJb', 'Merged_BLepWJa', 'Merged_BLepWJb', 'Merged_WJets', 'Lost_BHad', 'Lost_BLep', 'Lost_WJa', 'Lost_WJb']
                    for evt_cat in evt_cats:
                        evt_cat_mask = ak.flatten(valid_mp[evt_cat])
                        output['mtt'].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, cat=evt_cat, mtt=ak.flatten(valid_gen_objs[evt_cat_mask]['TTbar'].mass, axis=None), weight=wts[evt_cat_mask])

                            # for each gen parton find ones that pass and fail jet pt and eta cuts
                        for gen_obj in ['BHad', 'BLep', 'WJa', 'WJb']:
                            pt_eta_mask = ak.flatten((valid_gen_objs[evt_cat_mask][gen_obj].pt >= jet_pars['ptmin']) & (np.abs(valid_gen_objs[evt_cat_mask][gen_obj].eta) <= jet_pars['etamax']))
                            output['pass_mtt'].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, cat=evt_cat, objtype=gen_obj, mtt=ak.flatten(valid_gen_objs[evt_cat_mask][pt_eta_mask]['TTbar'].mass, axis=None), weight=wts[evt_cat_mask][pt_eta_mask])
                            output['fail_mtt'].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, cat=evt_cat, objtype=gen_obj, mtt=ak.flatten(valid_gen_objs[evt_cat_mask][~pt_eta_mask]['TTbar'].mass, axis=None), weight=wts[evt_cat_mask][~pt_eta_mask])
                        
                    ## fill hists of gen mttbar for 'other' (not merged/lost) events
                    other_cat_mask = ak.flatten(~(valid_mp['Merged_Event'] | valid_mp['Lost_Event']))
                    output['mtt'].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, cat='Other', mtt=ak.flatten(valid_gen_objs[other_cat_mask]['TTbar'].mass, axis=None), weight=wts[other_cat_mask])
                        # for each gen parton find ones that pass and fail jet pt and eta cuts
                    for gen_obj in ['BHad', 'BLep', 'WJa', 'WJb']:
                        pt_eta_mask = ak.flatten((valid_gen_objs[other_cat_mask][gen_obj].pt >= jet_pars['ptmin']) & (np.abs(valid_gen_objs[other_cat_mask][gen_obj].eta) <= jet_pars['etamax']))
                        output['pass_mtt'].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, cat='Other', objtype=gen_obj, mtt=ak.flatten(valid_gen_objs[other_cat_mask][pt_eta_mask]['TTbar'].mass, axis=None), weight=wts[other_cat_mask][pt_eta_mask])
                        output['fail_mtt'].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, cat='Other', objtype=gen_obj, mtt=ak.flatten(valid_gen_objs[other_cat_mask][~pt_eta_mask]['TTbar'].mass, axis=None), weight=wts[other_cat_mask][~pt_eta_mask])



        return output



    def postprocess(self, accumulator):
        return accumulator

proc_executor = processor.iterative_executor if to_debug else processor.futures_executor
proc_exec_args = {"schema": NanoAODSchema} if to_debug else {"schema": NanoAODSchema, "workers": 8}
output = processor.run_uproot_job(
    fileset,
    treename='Events',
    processor_instance=ttbar_evt_cats(),
    executor=proc_executor,
    executor_args=proc_exec_args,
    chunksize=100000,
)

save(output, args.outfname)
print('%s has been written' % args.outfname)

toc = time.time()
print("Total time: %.1f" % (toc - tic))
