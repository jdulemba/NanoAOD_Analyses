#!/usr/bin/env python

import time
tic = time.time()

from coffea.util import save, load
from coffea import hist, processor
from coffea.nanoevents import NanoAODSchema
from coffea.analysis_tools import PackedSelection
import awkward as ak
from coffea.nanoevents.methods import vector
ak.behavior.update(vector.behavior)
from pdb import set_trace
import os
import python.ObjectSelection as objsel
import Utilities.plot_tools as plt_tools
import python.MCWeights as MCWeights
import numpy as np
import Utilities.prettyjson as prettyjson
import python.GenParticleSelector as genpsel
import python.TTGenMatcher as ttmatcher
import Utilities.make_variables as make_vars
import python.Test_Permutator as tperm
import python.IDJet as IDJet

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
base_jobid = os.environ['base_jobid']
analyzer = 'permProbComputer'

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

# convert input string of options dictionary to actual dictionary
odict = (args.opts).replace("\'", "\"")
opts_dict = prettyjson.loads(odict)

    ## set config options passed through argparse
import ast
to_debug = ast.literal_eval(opts_dict.get('debug', 'False'))

    ## run on correct samples
Nominal_ttJets = ['ttJetsSL'] if base_jobid == 'ULnanoAOD' else ['ttJetsSL', 'ttJets_PS']
isTTbar = np.array([(key in Nominal_ttJets) for key in fileset.keys()]).all()
if not isTTbar:
    raise ValueError("This should only be run on nominal ttbar events!")

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
    'BTagSF' : False,
    'JetCor' : jet_corrections,
    'NNLO_Rewt' : {'Var' : nnlo_var, 'Correction' : nnlo_reweighting[nnlo_var]},
}

cfg_pars = prettyjson.loads(open(os.path.join(proj_dir, 'cfg_files', 'cfg_pars_%s.json' % jobid)).read())
    ## parameters for b-tagging
jet_pars = cfg_pars['Jets']
btagger = jet_pars['btagger']
wps_to_use = list(set([jet_pars['permutations']['tightb'],jet_pars['permutations']['looseb']]))
if not( len(wps_to_use) == 1):
    raise ValueError("Only 1 unique btag working point supported now")
btag_wp = btagger+wps_to_use[0]

MTcut = jet_pars['MT']

# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class permProbComputer(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.jetmult_axis = hist.Cat("jmult", "nJets")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")
        self.lepcat_axis = hist.Cat("lepcat", "Lepton Category")
        self.mtregion_axis = hist.Cat("mtregion", "MT Category")
        self.permcat_axis = hist.Cat("permcat", "permutation Category")
        self.tMass_axis = hist.Bin("topmass", "m(t_{had}) [GeV]", 500, 0., 500.)
        self.wMass_axis = hist.Bin("wmass", "m(W_{had}) [GeV]", 500, 0., 500.)
        self.nu_chi2_axis = hist.Bin("nu_chi2", r"$\chi_{\nu}^{2}$", 1000, 0., 1000.)
        self.nu_dist_axis = hist.Bin("nu_dist", r"$D_{\nu, min}$", 150, 0., 150.)
        self.maxmjet_axis = hist.Bin("maxmjet", "max m(jet) [GeV]", 500, 0., 500.)
        self.mThadProxy_axis = hist.Bin("mThadProxy", "m(b+j) [GeV]", 1000, 0., 2000.)

            ## make dictionary of hists
        histo_dict = {}
            ## make jet hists
                # 3j merged
        merged_3j_hists = self.make_3j_merged_hists()
        histo_dict.update(merged_3j_hists)
                # 3j lost 
        lost_3j_hists = self.make_3j_lost_hists()
        histo_dict.update(lost_3j_hists)
                # 4+ jets
        perm_4pj_hists = self.make_4pj_hists()
        histo_dict.update(perm_4pj_hists)

        histo_dict['cutflow'] = processor.defaultdict_accumulator(int)

        self._accumulator = processor.dict_accumulator(histo_dict)
        self.sample_name = ''
        self.corrections = corrections

    
    @property
    def accumulator(self):
        return self._accumulator

    def make_3j_merged_hists(self):
        histo_dict = {}
        histo_dict['Merged_mTHadProxy_vs_maxmjet'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepcat_axis, self.mtregion_axis, self.permcat_axis, self.maxmjet_axis, self.mThadProxy_axis)
        histo_dict['Merged_nusolver_chi2']         = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepcat_axis, self.mtregion_axis, self.permcat_axis, self.nu_chi2_axis)
        histo_dict['Merged_nusolver_dist']         = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepcat_axis, self.mtregion_axis, self.permcat_axis, self.nu_dist_axis)

        return histo_dict

    def make_3j_lost_hists(self):
        histo_dict = {}
        histo_dict['Lost_mTHadProxy']    = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepcat_axis, self.mtregion_axis, self.permcat_axis, self.mThadProxy_axis)
        histo_dict['Lost_nusolver_chi2'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepcat_axis, self.mtregion_axis, self.permcat_axis, self.nu_chi2_axis)
        histo_dict['Lost_nusolver_dist'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepcat_axis, self.mtregion_axis, self.permcat_axis, self.nu_dist_axis)

        return histo_dict


    def make_4pj_hists(self):
        histo_dict = {}
        histo_dict['mWHad_vs_mTHad'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepcat_axis, self.mtregion_axis, self.permcat_axis, self.tMass_axis, self.wMass_axis)
        histo_dict['nusolver_chi2']  = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepcat_axis, self.mtregion_axis, self.permcat_axis, self.nu_chi2_axis)
        histo_dict['nusolver_dist']  = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepcat_axis, self.mtregion_axis, self.permcat_axis, self.nu_dist_axis)

        return histo_dict

    
    def process(self, events):
        np.random.seed(10) # sets seed so values from random distributions are reproducible (JER corrections)
        output = self.accumulator.identity()

        self.sample_name = events.metadata['dataset']

            ## make event weights
                # data or MC distinction made internally
        mu_evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections)
        el_evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections)

            ## initialize selections and regions
        selection = PackedSelection()
        regions = {
            'Muon' : {
                #'LoT' : {
                #    '3Jets'  : {'lep_and_filter_pass', 'passing_jets', 'jets_3' , 'loose_or_tight_MU', 'btag_pass', 'semilep'},
                #    '4PJets' : {'lep_and_filter_pass', 'passing_jets', 'jets_4p', 'loose_or_tight_MU', 'btag_pass', 'semilep'},
                #},
                'Tight' : {
                    '3Jets'  : {'lep_and_filter_pass', 'passing_jets', 'jets_3' , 'tight_MU', 'btag_pass', 'semilep'},
                    '4PJets' : {'lep_and_filter_pass', 'passing_jets', 'jets_4p', 'tight_MU', 'btag_pass', 'semilep'},
                },
            },
            'Electron' : {
                #'LoT' : {
                #    '3Jets'  : {'lep_and_filter_pass', 'passing_jets', 'jets_3' , 'loose_or_tight_EL', 'btag_pass', 'semilep'},
                #    '4PJets' : {'lep_and_filter_pass', 'passing_jets', 'jets_4p', 'loose_or_tight_EL', 'btag_pass', 'semilep'},
                #},
                'Tight' : {
                    '3Jets'  : {'lep_and_filter_pass', 'passing_jets', 'jets_3' , 'tight_EL', 'btag_pass', 'semilep'},
                    '4PJets' : {'lep_and_filter_pass', 'passing_jets', 'jets_4p', 'tight_EL', 'btag_pass', 'semilep'},
                },
            },
        }

            # get all passing leptons
        lep_and_filter_pass = objsel.select_leptons(events, year=args.year, cutflow=output['cutflow'])
        selection.add('lep_and_filter_pass', lep_and_filter_pass)
                ## muons
        selection.add('tight_MU', ak.sum(events['Muon']['TIGHTMU'], axis=1) == 1) # one muon passing TIGHT criteria
        selection.add('loose_or_tight_MU', ak.sum(events['Muon']['LOOSEMU'] | events['Muon']['TIGHTMU'], axis=1) == 1) # one muon passing LOOSE or TIGHT criteria
                ## electrons
        selection.add('tight_EL', ak.sum(events['Electron']['TIGHTEL'], axis=1) == 1) # one muon passing TIGHT criteria
        selection.add('loose_or_tight_EL', ak.sum(events['Electron']['LOOSEEL'] | events['Electron']['TIGHTEL'], axis=1) == 1) # one muon passing LOOSE or TIGHT criteria

            ## build corrected jets and MET
        events['Jet'], events['MET'] = IDJet.process_jets(events, args.year, self.corrections['JetCor'])

            # jet selection
        passing_jets = objsel.jets_selection(events, year=args.year, cutflow=output['cutflow'])
        selection.add('passing_jets', passing_jets)
        output['cutflow']['nEvts passing jet and lepton obj selection'] += ak.sum(passing_jets & lep_and_filter_pass)
        selection.add('jets_3', ak.num(events['SelectedJets']) == 3)
        selection.add('jets_4p', ak.num(events['SelectedJets']) > 3)
        selection.add('btag_pass', ak.sum(events['SelectedJets'][btag_wp], axis=1) >= 2)

            # sort jets by btag value
        events['SelectedJets'] = events['SelectedJets'][ak.argsort(events['SelectedJets']['btagDeepB'], ascending=False)] if btagger == 'DeepCSV' else events['SelectedJets'][ak.argsort(events['SelectedJets']['btagDeepFlavB'], ascending=False)]

        ## apply lepton SFs to MC (only applicable to tight leptons)
        if 'LeptonSF' in corrections.keys():
            tight_mu_cut = selection.require(tight_MU=True) # find events passing muon object selection with one tight muon
            tight_muons = events['Muon'][tight_mu_cut][(events['Muon'][tight_mu_cut]['TIGHTMU'] == True)]
            muSFs_dict =  MCWeights.get_lepton_sf(year=args.year, lepton='Muons', corrections=self.corrections['LeptonSF'],
                pt=ak.flatten(tight_muons['pt']), eta=ak.flatten(tight_muons['eta']))
            mu_reco_cen = np.ones(len(events))
            mu_reco_err = np.zeros(len(events))
            mu_trig_cen = np.ones(len(events))
            mu_trig_err = np.zeros(len(events))
            mu_reco_cen[tight_mu_cut] = muSFs_dict['RECO_CEN']
            mu_reco_err[tight_mu_cut] = muSFs_dict['RECO_ERR']
            mu_trig_cen[tight_mu_cut] = muSFs_dict['TRIG_CEN']
            mu_trig_err[tight_mu_cut] = muSFs_dict['TRIG_ERR']
            mu_evt_weights.add('RECO', mu_reco_cen, mu_reco_err, mu_reco_err, shift=True)
            mu_evt_weights.add('TRIG', mu_trig_cen, mu_trig_err, mu_trig_err, shift=True)

            tight_el_cut = selection.require(tight_EL=True) # find events passing electron object selection with one tight electron
            tight_electrons = events['Electron'][tight_el_cut][(events['Electron'][tight_el_cut]['TIGHTEL'] == True)]
            elSFs_dict = MCWeights.get_lepton_sf(year=args.year, lepton='Electrons', corrections=self.corrections['LeptonSF'],
                pt=ak.flatten(tight_electrons['pt']), eta=ak.flatten(tight_electrons['etaSC']))
            el_reco_cen = np.ones(len(events))
            el_reco_err = np.zeros(len(events))
            el_trig_cen = np.ones(len(events))
            el_trig_err = np.zeros(len(events))
            el_reco_cen[tight_el_cut] = elSFs_dict['RECO_CEN']
            el_reco_err[tight_el_cut] = elSFs_dict['RECO_ERR']
            el_trig_cen[tight_el_cut] = elSFs_dict['TRIG_CEN']
            el_trig_err[tight_el_cut] = elSFs_dict['TRIG_ERR']
            el_evt_weights.add('RECO', el_reco_cen, el_reco_err, el_reco_err, shift=True)
            el_evt_weights.add('TRIG', el_trig_cen, el_trig_err, el_trig_err, shift=True)

            # find gen level particles for ttbar system
        genpsel.select(events, mode='NORMAL')
        selection.add('semilep', ak.num(events['SL']) > 0)
        if 'NNLO_Rewt' in self.corrections.keys():
            nnlo_wts = MCWeights.get_nnlo_weights(self.corrections['NNLO_Rewt'], events)
            mu_evt_weights.add('%s_reweighting' % corrections['NNLO_Rewt']['Var'], nnlo_wts)
            el_evt_weights.add('%s_reweighting' % corrections['NNLO_Rewt']['Var'], nnlo_wts)

        ## fill hists for each region
        for lepton in regions.keys():
            evt_weights = mu_evt_weights if lepton == 'Muon' else el_evt_weights
            for lepcat in regions[lepton].keys():
                for jmult in regions[lepton][lepcat].keys():
                    cut = selection.all(*regions[lepton][lepcat][jmult])
                    wts = evt_weights.weight()[cut] # get event weights

                    if to_debug: print(lepton, lepcat, jmult)
                    #print(lepton, lepcat, jmult)
                    #set_trace()
                    if ak.sum(cut) > 0:
                        leptype = 'MU' if lepton == 'Muon' else 'EL'
                        if 'loose_or_tight_%s' % leptype in regions[lepton][lepcat][jmult]:
                            leptons = events[lepton][cut][((events[lepton][cut]['TIGHT%s' % leptype] == True) | (events[lepton][cut]['LOOSE%s' % leptype] == True))]
                        elif 'tight_%s' % leptype in regions[lepton][lepcat][jmult]:
                            leptons = events[lepton][cut][(events[lepton][cut]['TIGHT%s' % leptype] == True)]
                        elif 'loose_%s' % leptype in regions[lepton][lepcat][jmult]:
                            leptons = events[lepton][cut][(events[lepton][cut]['LOOSE%s' % leptype] == True)]
                        else:
                            raise ValueError("Not sure what lepton type to choose for event")

                            # get jets and MET
                        jets, met = events['SelectedJets'][cut], events['SelectedMET'][cut]

                            ## create MT regions
                        MT = make_vars.MT(leptons, events['MET'][cut])
                        MTHigh = ak.flatten(MT >= MTcut)


                        test_perms = tperm.find_permutations(jets=jets, leptons=leptons, MET=met, btagWP=btag_wp)
                        matched_perm = ttmatcher.best_match(gen_hyp=events['SL'][cut], jets=jets, leptons=leptons, met=met)

                        #set_trace()
                        if jmult == '3Jets':
                            output = self.make_3j_categories(accumulator=output, jetmult=jmult, leptype=lepton, lepcat=lepcat, mtregion='MTHigh', jets=jets[MTHigh],  tps=test_perms[MTHigh],  mp=matched_perm[MTHigh],  evt_weights=wts[MTHigh])
                            #output = self.make_3j_categories(accumulator=output, jetmult=jmult, leptype=lepton, lepcat=lepcat, mtregion='MTLow',  jets=jets[~MTHigh], tps=test_perms[~MTHigh], mp=matched_perm[~MTHigh], evt_weights=wts[~MTHigh])
                        else:
                            output = self.fill_4pj_hists(accumulator=output, jetmult=jmult, leptype=lepton, lepcat=lepcat, mtregion='MTHigh', jets=jets[MTHigh],  tps=test_perms[MTHigh],  mp=matched_perm[MTHigh],  evt_weights=wts[MTHigh])
                            #output = self.fill_4pj_hists(accumulator=output, jetmult=jmult, leptype=lepton, lepcat=lepcat, mtregion='MTLow',  jets=jets[~MTHigh], tps=test_perms[~MTHigh], mp=matched_perm[~MTHigh], evt_weights=wts[~MTHigh])

        return output

    def fill_4pj_hists(self, accumulator, jetmult, leptype, lepcat, mtregion, jets, tps, mp, evt_weights):
        #set_trace()
        valid_evts = (ak.num(mp['TTbar'].pt) > 0) & (ak.flatten(mp['unique_matches'] == 4))
        if ak.sum(valid_evts) > 0:
            evt_wts = ak.ones_like(tps[valid_evts].blepIdx)*evt_weights[valid_evts] # make event weights same shape as test perms

            tp_blepInds, mp_blepInds = ak.unzip(ak.cartesian( [tps[valid_evts].blepIdx, mp[valid_evts]['BLep'].jetIdx] )) # make mp inds same size as test perms, then split
            tp_bhadInds, mp_bhadInds = ak.unzip(ak.cartesian( [tps[valid_evts].bhadIdx, mp[valid_evts]['BHad'].jetIdx] )) # make mp inds same size as test perms, then split
            tp_wjaInds, mp_wjaInds   = ak.unzip(ak.cartesian( [tps[valid_evts].wjaIdx, mp[valid_evts]['WJa'].jetIdx] )) # make mp inds same size as test perms, then split
            tp_wjbInds, mp_wjbInds   = ak.unzip(ak.cartesian( [tps[valid_evts].wjbIdx, mp[valid_evts]['WJb'].jetIdx] )) # make mp inds same size as test perms, then split

            mW = (jets[valid_evts][tp_wjaInds] + jets[valid_evts][tp_wjbInds]).mass
            mT = (jets[valid_evts][tp_wjaInds] + jets[valid_evts][tp_wjbInds] + jets[valid_evts][tp_bhadInds]).mass

                # define when tp objects are correctly matched
            correct_blep_matched, correct_bhad_matched = (tp_blepInds == mp_blepInds), (tp_bhadInds == mp_bhadInds)
            correct_whad_matched = ( ((tp_wjaInds == mp_wjaInds) & (tp_wjbInds == mp_wjbInds)) | ((tp_wjaInds == mp_wjbInds) & (tp_wjbInds == mp_wjaInds)) )

                # all tp objects correctly assigned 
            all_correct_tpInds = correct_blep_matched & correct_bhad_matched & correct_whad_matched
                # only thad objects correctly assigned for tp
            only_right_thad_tpInds = (~correct_blep_matched) & correct_bhad_matched & correct_whad_matched
                # only blep (tlep) objects correctly assigned for tp
            only_right_tlep_tpInds = correct_blep_matched & (~correct_bhad_matched) & (~correct_whad_matched)
                # all other assignment categories
            all_wrong_tpInds =  ~(all_correct_tpInds | only_right_thad_tpInds | only_right_tlep_tpInds)

                # find Nu chi2 solutions for all tp categories
            all_correct_tp_NuChi2 = ak.flatten(tps[valid_evts]['Nu']['chi2'][all_correct_tpInds], axis=None)
            only_right_thad_tp_NuChi2 = ak.flatten(tps[valid_evts]['Nu']['chi2'][only_right_thad_tpInds], axis=None)
            only_right_tlep_tp_NuChi2 = ak.flatten(tps[valid_evts]['Nu']['chi2'][only_right_tlep_tpInds], axis=None)
            all_wrong_tp_NuChi2 = ak.flatten(tps[valid_evts]['Nu']['chi2'][all_wrong_tpInds], axis=None)

            # fill hists
                # Correct
            accumulator['nusolver_chi2'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, mtregion=mtregion, permcat='Correct', nu_chi2=all_correct_tp_NuChi2, weight=ak.flatten(evt_wts[all_correct_tpInds], axis=None))
            accumulator['nusolver_dist'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, mtregion=mtregion, permcat='Correct', nu_dist=np.sqrt(all_correct_tp_NuChi2), weight=ak.flatten(evt_wts[all_correct_tpInds], axis=None))
            accumulator['mWHad_vs_mTHad'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, mtregion=mtregion, permcat='Correct', topmass=ak.flatten(mT[all_correct_tpInds], axis=None), wmass=ak.flatten(mW[all_correct_tpInds], axis=None), weight=ak.flatten(evt_wts[all_correct_tpInds], axis=None))
                # Right THad
            accumulator['nusolver_chi2'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, mtregion=mtregion, permcat='Right_THad', nu_chi2=only_right_thad_tp_NuChi2, weight=ak.flatten(evt_wts[only_right_thad_tpInds], axis=None))
            accumulator['nusolver_dist'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, mtregion=mtregion, permcat='Right_THad', nu_dist=np.sqrt(only_right_thad_tp_NuChi2), weight=ak.flatten(evt_wts[only_right_thad_tpInds], axis=None))
            accumulator['mWHad_vs_mTHad'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, mtregion=mtregion, permcat='Right_THad', topmass=ak.flatten(mT[only_right_thad_tpInds], axis=None), wmass=ak.flatten(mW[only_right_thad_tpInds], axis=None), weight=ak.flatten(evt_wts[only_right_thad_tpInds], axis=None))
                # Right TLep
            accumulator['nusolver_chi2'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, mtregion=mtregion, permcat='Right_TLep', nu_chi2=only_right_tlep_tp_NuChi2, weight=ak.flatten(evt_wts[only_right_tlep_tpInds], axis=None))
            accumulator['nusolver_dist'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, mtregion=mtregion, permcat='Right_TLep', nu_dist=np.sqrt(only_right_tlep_tp_NuChi2), weight=ak.flatten(evt_wts[only_right_tlep_tpInds], axis=None))
            accumulator['mWHad_vs_mTHad'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, mtregion=mtregion, permcat='Right_TLep', topmass=ak.flatten(mT[only_right_tlep_tpInds], axis=None), wmass=ak.flatten(mW[only_right_tlep_tpInds], axis=None), weight=ak.flatten(evt_wts[only_right_tlep_tpInds], axis=None))
                # Wrong
            accumulator['nusolver_chi2'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, mtregion=mtregion, permcat='Wrong', nu_chi2=all_wrong_tp_NuChi2, weight=ak.flatten(evt_wts[all_wrong_tpInds], axis=None))
            accumulator['nusolver_dist'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, mtregion=mtregion, permcat='Wrong', nu_dist=np.sqrt(all_wrong_tp_NuChi2), weight=ak.flatten(evt_wts[all_wrong_tpInds], axis=None))
            accumulator['mWHad_vs_mTHad'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, mtregion=mtregion, permcat='Wrong', topmass=ak.flatten(mT[all_wrong_tpInds], axis=None), wmass=ak.flatten(mW[all_wrong_tpInds], axis=None), weight=ak.flatten(evt_wts[all_wrong_tpInds], axis=None))

        return accumulator        


    def make_3j_categories(self, accumulator, jetmult, leptype, lepcat, mtregion, jets, tps, mp, evt_weights):
        #set_trace()
        valid_evts = (ak.num(mp['TTbar'].pt) > 0) & (ak.flatten(mp['unique_matches'] >= 3))

        correct_merged = valid_evts & ak.flatten(mp['Merged_Event'] & (mp['Merged_BHadWJa'] | mp['Merged_BHadWJb'] | mp['Merged_WJets']))
        correct_lost   = valid_evts & ak.flatten(mp['Lost_Event'] & (mp['Lost_WJa'] | mp['Lost_WJb']))

            ## Lost Events
        if ak.sum(correct_lost) > 0:
            if to_debug: print("\tAnalyzing Lost")
            lost_evt_wts = ak.ones_like(tps[correct_lost].blepIdx)*evt_weights[correct_lost] # make event weights same shape as test perms

            lost_tp_blepInds, lost_mp_blepInds = ak.unzip(ak.cartesian( [tps[correct_lost].blepIdx, mp[correct_lost]['BLep'].jetIdx] )) # make mp inds same size as test perms, then split
            lost_tp_bhadInds, lost_mp_bhadInds = ak.unzip(ak.cartesian( [tps[correct_lost].bhadIdx, mp[correct_lost]['BHad'].jetIdx] )) # make mp inds same size as test perms, then split
                        # find inds where wja or wjb exist for mp
            lost_mp_wjaInds, lost_mp_wjbInds = ak.flatten(ak.fill_none(ak.pad_none(mp[correct_lost]['WJa'].jetIdx, 1), -999)), ak.flatten(ak.fill_none(ak.pad_none(mp[correct_lost]['WJb'].jetIdx, 1), -999))
            lost_tp_wjaInds,  lost_mp_wjetInds = ak.unzip(ak.cartesian( [tps[correct_lost].wjaIdx, ak.where(lost_mp_wjaInds > 0, lost_mp_wjaInds, lost_mp_wjbInds)] )) # make mp inds same size as test perms, then split
                    # find correctly assigned test perm jet assignments
            correct_lost_tpInds = (lost_tp_blepInds == lost_mp_blepInds) & (lost_tp_bhadInds == lost_mp_bhadInds) & (lost_tp_wjaInds == lost_mp_wjetInds)
            correct_lost_tp_NuChi2 = ak.flatten(tps[correct_lost]['Nu']['chi2'][correct_lost_tpInds], axis=None)
            wrong_lost_tpInds = ~correct_lost_tpInds
            wrong_lost_tp_NuChi2 = ak.flatten(tps[correct_lost]['Nu']['chi2'][wrong_lost_tpInds], axis=None)
            lost_tp_thad_proxy_mass = (jets[correct_lost][lost_tp_bhadInds]+jets[correct_lost][lost_tp_wjaInds]).mass

                    # fill hists
                        # Correct
            accumulator['Lost_nusolver_chi2'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, mtregion=mtregion, permcat='Correct', nu_chi2=correct_lost_tp_NuChi2, weight=ak.flatten(lost_evt_wts[correct_lost_tpInds], axis=None))
            accumulator['Lost_nusolver_dist'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, mtregion=mtregion, permcat='Correct', nu_dist=np.sqrt(correct_lost_tp_NuChi2), weight=ak.flatten(lost_evt_wts[correct_lost_tpInds], axis=None))
            accumulator['Lost_mTHadProxy'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, mtregion=mtregion, permcat='Correct', mThadProxy=ak.flatten(lost_tp_thad_proxy_mass[correct_lost_tpInds], axis=None), weight=ak.flatten(lost_evt_wts[correct_lost_tpInds], axis=None))
                        # Wrong
            accumulator['Lost_nusolver_chi2'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, mtregion=mtregion, permcat='Wrong', nu_chi2=wrong_lost_tp_NuChi2, weight=ak.flatten(lost_evt_wts[wrong_lost_tpInds], axis=None))
            accumulator['Lost_nusolver_dist'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, mtregion=mtregion, permcat='Wrong', nu_dist=np.sqrt(wrong_lost_tp_NuChi2), weight=ak.flatten(lost_evt_wts[wrong_lost_tpInds], axis=None))
            accumulator['Lost_mTHadProxy'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, mtregion=mtregion, permcat='Wrong', mThadProxy=ak.flatten(lost_tp_thad_proxy_mass[wrong_lost_tpInds], axis=None), weight=ak.flatten(lost_evt_wts[wrong_lost_tpInds], axis=None))

            ## Merged Events
        if ak.sum(correct_merged) > 0:
            if to_debug: print("\tAnalyzing Merged")
            merged_evt_wts = ak.ones_like(tps[correct_merged].blepIdx)*evt_weights[correct_merged] # make event weights same shape as test perms

            merged_tp_blepInds, merged_mp_blepInds = ak.unzip(ak.cartesian( [tps[correct_merged].blepIdx, mp[correct_merged]['BLep'].jetIdx] )) # make mp inds same size as test perms, then split
            merged_tp_bhadInds, merged_mp_bhadInds = ak.unzip(ak.cartesian( [tps[correct_merged].bhadIdx, mp[correct_merged]['BHad'].jetIdx] )) # make mp inds same size as test perms, then split
            merged_tp_wjaInds, merged_mp_wjaInds   = ak.unzip(ak.cartesian( [tps[correct_merged].wjaIdx, mp[correct_merged]['WJa'].jetIdx] )) # make mp inds same size as test perms, then split
            merged_tp_wjbInds, merged_mp_wjbInds   = ak.unzip(ak.cartesian( [tps[correct_merged].wjbIdx, mp[correct_merged]['WJb'].jetIdx] )) # make mp inds same size as test perms, then split
                    # find correctly assigned test perm jet assignments
            correct_merged_tpInds = (merged_tp_blepInds == merged_mp_blepInds) & (merged_tp_bhadInds == merged_mp_bhadInds) & ( (merged_tp_wjaInds == merged_mp_wjaInds) | (merged_tp_wjaInds == merged_mp_wjbInds) )
            correct_merged_tp_NuChi2 = ak.flatten(tps[correct_merged]['Nu']['chi2'][correct_merged_tpInds], axis=None)
            wrong_merged_tpInds = ~correct_merged_tpInds
            wrong_merged_tp_NuChi2 = ak.flatten(tps[correct_merged]['Nu']['chi2'][wrong_merged_tpInds], axis=None)
            merged_tp_thad_proxy_mass = (jets[correct_merged][merged_tp_bhadInds]+jets[correct_merged][merged_tp_wjaInds]).mass
            merged_tp_maxmjet = ak.where(jets[correct_merged][merged_tp_bhadInds].mass > jets[correct_merged][merged_tp_wjaInds].mass, jets[correct_merged][merged_tp_bhadInds].mass, jets[correct_merged][merged_tp_wjaInds].mass)

                    # fill hists
                        # Correct
            accumulator['Merged_nusolver_chi2'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, mtregion=mtregion, permcat='Correct', nu_chi2=correct_merged_tp_NuChi2, weight=ak.flatten(merged_evt_wts[correct_merged_tpInds], axis=None))
            accumulator['Merged_nusolver_dist'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, mtregion=mtregion, permcat='Correct', nu_dist=np.sqrt(correct_merged_tp_NuChi2), weight=ak.flatten(merged_evt_wts[correct_merged_tpInds], axis=None))
            accumulator['Merged_mTHadProxy_vs_maxmjet'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, mtregion=mtregion, permcat='Correct', maxmjet=ak.flatten(merged_tp_maxmjet[correct_merged_tpInds], axis=None), mThadProxy=ak.flatten(merged_tp_thad_proxy_mass[correct_merged_tpInds], axis=None), weight=ak.flatten(merged_evt_wts[correct_merged_tpInds], axis=None))
                        # Wrong
            accumulator['Merged_nusolver_chi2'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, mtregion=mtregion, permcat='Wrong', nu_chi2=wrong_merged_tp_NuChi2, weight=ak.flatten(merged_evt_wts[wrong_merged_tpInds], axis=None))
            accumulator['Merged_nusolver_dist'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, mtregion=mtregion, permcat='Wrong', nu_dist=np.sqrt(wrong_merged_tp_NuChi2), weight=ak.flatten(merged_evt_wts[wrong_merged_tpInds], axis=None))
            accumulator['Merged_mTHadProxy_vs_maxmjet'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, mtregion=mtregion, permcat='Wrong', maxmjet=ak.flatten(merged_tp_maxmjet[wrong_merged_tpInds], axis=None), mThadProxy=ak.flatten(merged_tp_thad_proxy_mass[wrong_merged_tpInds], axis=None), weight=ak.flatten(merged_evt_wts[wrong_merged_tpInds], axis=None))

        return accumulator


    def postprocess(self, accumulator):
        return accumulator


proc_executor = processor.iterative_executor if to_debug else processor.futures_executor
proc_exec_args = {"schema": NanoAODSchema} if to_debug else {"schema": NanoAODSchema, "workers": 8}
output = processor.run_uproot_job(
    fileset,
    treename="Events",
    processor_instance=permProbComputer(),
    executor=proc_executor,
    executor_args=proc_exec_args,
    #chunksize=10000 if to_debug else 100000,
    chunksize=100000,
)

save(output, args.outfname)
print('%s has been written' % args.outfname)

toc = time.time()
print("Total time: %.1f" % (toc - tic))
