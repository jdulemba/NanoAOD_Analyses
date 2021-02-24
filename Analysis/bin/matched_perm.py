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
import numpy as np
import Utilities.prettyjson as prettyjson
import python.GenParticleSelector as genpsel
import python.TTGenMatcher as ttmatcher
import python.TTPermutator as ttpermutator
import Utilities.make_variables as make_vars
import python.ObjectSelection as objsel
import python.MCWeights as MCWeights
from python.IDJet import btag_values as btag_values

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
base_jobid = os.environ['base_jobid']
analyzer = 'matched_perm'

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('fset', type=str, help='Fileset dictionary (in string form) to be used for the processor')
parser.add_argument('year', choices=['2016APV', '2016', '2017', '2018'] if base_jobid == 'ULnanoAOD' else ['2016', '2017', '2018'], help='Specify which year to run over')
parser.add_argument('outfname', type=str, help='Specify output filename, including directory and file extension')
parser.add_argument('--debug', action='store_true', help='Uses iterative_executor for debugging purposes, otherwise futures_excutor will be used (faster)')

args = parser.parse_args()

# convert input string of fileset dictionary to actual dictionary
fdict = (args.fset).replace("\'", "\"")
fileset = prettyjson.loads(fdict)

## init tt probs for likelihoods
ttpermutator.year_to_run(year=args.year)

    ## run on correct samples
Nominal_ttJets = ['ttJets_PS', 'ttJets'] if ((args.year == '2016') and (base_jobid == 'NanoAODv6')) else ['ttJetsSL']
isTTbar = np.array([(key in Nominal_ttJets) for key in fileset.keys()]).all()
if not isTTbar:
    raise ValueError("This should only be run on nominal ttbar events!")

init_btag = ~(np.array([key.startswith('data') for key in fileset.keys()]).all())

## load corrections for event weights
pu_correction = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'MC_PU_Weights.coffea'))[args.year]
lepSF_correction = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'leptonSFs.coffea'))[args.year]
jet_corrections = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'JetMETCorrections.coffea'))[args.year]
corrections = {
    'Pileup' : pu_correction,
    'Prefire' : True,
    'LeptonSF' : lepSF_correction,
    'BTagSF' : init_btag,
    'JetCor' : jet_corrections,
}


cfg_pars = prettyjson.loads(open(os.path.join(proj_dir, 'cfg_files', 'cfg_pars_%s.json' % jobid)).read())
    ## parameters for b-tagging
jet_pars = cfg_pars['Jets']
#btagger = jet_pars['btagger']
#btaggers = ['DeepJet', 'DeepCSV']
btaggers = ['DeepCSV']

wps_to_use = list(set([jet_pars['permutations']['tightb'],jet_pars['permutations']['looseb']]))
if not( len(wps_to_use) == 1):
    raise IOError("Only 1 unique btag working point supported now")
#btag_wp = btagger+wps_to_use[0]
#btag_wps = ['DeepJetMedium', 'DeepCSVMedium']
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
class matched_perm(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.objtype_axis = hist.Cat("objtype", "Object Type")
        self.jetmult_axis = hist.Cat("jmult", "nJets")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")
        self.lepcat_axis = hist.Cat("lepcat", "Lepton Category")
        self.btag_axis = hist.Cat("btag", "BTag Region")
        self.pt_axis     = hist.Bin("pt", "p_{T} [GeV]", 200, 0, 1000)
        self.eta_axis    = hist.Bin("eta", r"$\eta$", 60, -3, 3)
        self.phi_axis    = hist.Bin("phi", r"$\phi$", 160, -4, 4)
        self.mass_axis   = hist.Bin("mass", "m [GeV]", 1000, 0., 2000.)
        self.energy_axis = hist.Bin("energy", "E [GeV]", 200, 0, 1000)
        self.reso_pt_axis     = hist.Bin("reso_pt", "p_{T} [GeV]", 200, -200, 200)
        self.reso_eta_axis    = hist.Bin("reso_eta", r"$\eta$", 100, -5, 5)
        self.reso_phi_axis    = hist.Bin("reso_phi", r"$\phi$", 100, -5, 5)
        self.reso_mass_axis   = hist.Bin("reso_mass", "m [GeV]", 1000, -200., 200.)
        self.reso_energy_axis = hist.Bin("reso_energy", "E [GeV]", 1000, -200, 200)

            ## make dictionary of hists
        histo_dict = {}
            ## make jet hists
                # Gen
        gen_hists = self.make_gen_hists()
        histo_dict.update(gen_hists)
                # Reco
        reco_hists = self.make_reco_hists()
        histo_dict.update(reco_hists)
                # Reso
        reso_hists = self.make_reso_hists()
        histo_dict.update(reso_hists)

        histo_dict['cutflow'] = processor.defaultdict_accumulator(int)

        self._accumulator = processor.dict_accumulator(histo_dict)
        self.sample_name = ''
        self.corrections = corrections

    
    @property
    def accumulator(self):
        return self._accumulator

    def make_gen_hists(self):
        histo_dict = {}
        histo_dict['Gen_pt']     = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.objtype_axis, self.pt_axis)
        histo_dict['Gen_eta']    = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.objtype_axis, self.eta_axis)
        histo_dict['Gen_phi']    = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.objtype_axis, self.phi_axis)
        histo_dict['Gen_mass']   = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.objtype_axis, self.mass_axis)
        histo_dict['Gen_energy'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.objtype_axis, self.energy_axis)

        return histo_dict

    def make_reco_hists(self):
        histo_dict = {}
        histo_dict['Reco_pt']     = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.objtype_axis, self.pt_axis)
        histo_dict['Reco_eta']    = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.objtype_axis, self.eta_axis)
        histo_dict['Reco_phi']    = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.objtype_axis, self.phi_axis)
        histo_dict['Reco_mass']   = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.objtype_axis, self.mass_axis)
        histo_dict['Reco_energy'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.objtype_axis, self.energy_axis)

        return histo_dict

    def make_reso_hists(self):
        histo_dict = {}
        histo_dict['Reso_pt']     = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.objtype_axis, self.reso_pt_axis)
        histo_dict['Reso_eta']    = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.objtype_axis, self.reso_eta_axis)
        histo_dict['Reso_phi']    = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.objtype_axis, self.reso_phi_axis)
        histo_dict['Reso_mass']   = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.objtype_axis, self.reso_mass_axis)
        histo_dict['Reso_energy'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.objtype_axis, self.reso_energy_axis)

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
                'Tight' : {
                    'btagPass' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'tight_MU', 'DeepCSV_pass', 'semilep'},
                        '4PJets' : {'objselection', 'jets_4p', 'tight_MU', 'DeepCSV_pass', 'semilep'},
                    },
                    'btagFail' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'tight_MU', 'DeepCSV_fail', 'semilep'},
                        '4PJets' : {'objselection', 'jets_4p', 'tight_MU', 'DeepCSV_fail', 'semilep'},
                    },
                },
                'Loose' : {
                    'btagPass' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'loose_MU', 'DeepCSV_pass', 'semilep'},
                        '4PJets' : {'objselection', 'jets_4p', 'loose_MU', 'DeepCSV_pass', 'semilep'},
                    },
                    'btagFail' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'loose_MU', 'DeepCSV_fail', 'semilep'},
                        '4PJets' : {'objselection', 'jets_4p', 'loose_MU', 'DeepCSV_fail', 'semilep'},
                    },
                },
            },
            'Electron' : {
                'Tight' : {
                    'btagPass' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'tight_EL', 'DeepCSV_pass', 'semilep'},
                        '4PJets' : {'objselection', 'jets_4p', 'tight_EL', 'DeepCSV_pass', 'semilep'},
                    },
                    'btagFail' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'tight_EL', 'DeepCSV_fail', 'semilep'},
                        '4PJets' : {'objselection', 'jets_4p', 'tight_EL', 'DeepCSV_fail', 'semilep'},
                    },
                },
                'Loose' : {
                    'btagPass' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'loose_EL', 'DeepCSV_pass', 'semilep'},
                        '4PJets' : {'objselection', 'jets_4p', 'loose_EL', 'DeepCSV_pass', 'semilep'},
                    },
                    'btagFail' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'loose_EL', 'DeepCSV_fail', 'semilep'},
                        '4PJets' : {'objselection', 'jets_4p', 'loose_EL', 'DeepCSV_fail', 'semilep'},
                    },
                },
            },
        }

            ## object selection
        objsel_evts = objsel.select(events, year=args.year, corrections=self.corrections, cutflow=output['cutflow'])
        output['cutflow']['nEvts passing jet and lepton obj selection'] += ak.sum(objsel_evts)
        selection.add('objselection', objsel_evts)
        selection.add('jets_3',  ak.num(events['Jet']) == 3)
        selection.add('jets_4p', ak.num(events['Jet']) > 3)
        selection.add('DeepCSV_pass', ak.sum(events['Jet']['DeepCSV'+wps_to_use[0]], axis=1) >= 2)

            # sort jets by btag value
        events['Jet'] = events['Jet'][ak.argsort(events['Jet']['btagDeepB'], ascending=False)] if btagger == 'DeepCSV' else events['Jet'][ak.argsort(events['Jet']['btagDeepFlavB'], ascending=False)]

            # btag fail sideband
        deepcsv_sorted = events['Jet'][ak.argsort(events['Jet']['btagDeepB'], ascending=False)]['btagDeepB']
        valid_counts_inds = ak.where(ak.num(events['Jet']) > 1)[0]
        deepcsv_fail = np.zeros(len(events)).astype(bool)
        deepcsv_fail[valid_counts_inds] = (deepcsv_sorted[valid_counts_inds][:, 0] < btag_values[args.year]['btagDeepB']['DeepCSV'+wps_to_use[0]]) & (deepcsv_sorted[valid_counts_inds][:, 1] < btag_values[args.year]['btagDeepB']['DeepCSV'+wps_to_use[0]])
        selection.add('DeepCSV_fail', deepcsv_fail ) # highest and second highest DeepCSV values don't pass tight and loose WPs

            ## add different selections
                ## muons
        selection.add('tight_MU', ak.sum(events['Muon']['TIGHTMU'], axis=1) == 1) # one muon passing TIGHT criteria
        selection.add('loose_MU', ak.sum(events['Muon']['LOOSEMU'], axis=1) == 1) # one muon passing LOOSE criteria

                ## electrons
        selection.add('tight_EL', ak.sum(events['Electron']['TIGHTEL'], axis=1) == 1) # one electron passing TIGHT criteria
        selection.add('loose_EL', ak.sum(events['Electron']['LOOSEEL'], axis=1) == 1) # one electron passing LOOSE criteria

        ### apply lepton SFs to MC (only applicable to tight leptons)
        if 'LeptonSF' in corrections.keys():
            tight_mu_cut = selection.require(objselection=True, tight_MU=True) # find events passing muon object selection with one tight muon
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

            tight_el_cut = selection.require(objselection=True, tight_EL=True) # find events passing electron object selection with one tight electron
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

        if corrections['BTagSF'] == True:
            deepcsv_cen   = np.ones(len(events))
            threeJets_cut = selection.require(objselection=True, jets_3=True)
            deepcsv_3j_wts = self.corrections['BTag_Constructors']['DeepCSV']['3Jets'].get_scale_factor(jets=events['Jet'][threeJets_cut], passing_cut='DeepCSV'+wps_to_use[0])
            deepcsv_cen[threeJets_cut] = ak.prod(deepcsv_3j_wts['central'], axis=1)

            fourplusJets_cut = selection.require(objselection=True, jets_4p=True)
            deepcsv_4pj_wts = self.corrections['BTag_Constructors']['DeepCSV']['4PJets'].get_scale_factor(jets=events['Jet'][fourplusJets_cut], passing_cut='DeepCSV'+wps_to_use[0])
            deepcsv_cen[fourplusJets_cut] = ak.prod(deepcsv_4pj_wts['central'], axis=1)

                # make dict of btag weights
            btag_weights = {
                'DeepCSV_CEN' : deepcsv_cen,
            }


            # find gen level particles for ttbar system
        genpsel.select(events, mode='NORMAL')
        selection.add('semilep', ak.num(events['SL']) > 0)

        ## fill hists for each region
        for lepton in regions.keys():
            evt_weights = mu_evt_weights if lepton == 'Muon' else el_evt_weights
            for leptype in regions[lepton].keys():
                for btagregion in regions[lepton][leptype].keys():
                    for jmult in regions[lepton][leptype][btagregion].keys():
                        cut = selection.all(*regions[lepton][leptype][btagregion][jmult])

                        if args.debug: print(lepton, leptype, btagregion, jmult)
                        if cut.sum() > 0:
                            ltype = 'MU' if lepton == 'Muon' else 'EL'
                            if 'loose_or_tight_%s' % ltype in regions[lepton][leptype][btagregion][jmult]:
                                lep_mask = ((events[lepton][cut]['TIGHT%s' % ltype] == True) | (events[lepton][cut]['LOOSE%s' % ltype] == True))
                            elif 'tight_%s' % ltype in regions[lepton][leptype][btagregion][jmult]:
                                lep_mask = (events[lepton][cut]['TIGHT%s' % ltype] == True)
                            elif 'loose_%s' % ltype in regions[lepton][leptype][btagregion][jmult]:
                                lep_mask = (events[lepton][cut]['LOOSE%s' % ltype] == True)
                            else:
                                raise ValueError("Not sure what lepton type to choose for event")

                                # find matched perm
                            mp = ttmatcher.best_match(gen_hyp=events['SL'][cut], jets=events['Jet'][cut], leptons=events[lepton][cut][lep_mask], met=events['MET'][cut])

                                ## create MT regions
                            MT = make_vars.MT(events[lepton][cut][lep_mask], events['MET'][cut])
                            MTHigh = ak.flatten(MT >= MTcut)


                            wts = (evt_weights.weight()*btag_weights['%s_CEN' % btaggers[0]])[cut][MTHigh]

                            mp_status = np.zeros(cut.sum(), dtype=int) # 0 == '' (no gen matching), 1 == 'right', 2 == 'matchable', 3 == 'unmatchable', 4 == 'sl_tau', 5 == 'noslep'
                            if jmult == '3Jets':
                                valid_evts = (ak.num(mp['TTbar'].pt) > 0) & (ak.flatten(mp['unique_matches'] >= 3))
                        
                                    # merged events
                                merged_evts = valid_evts & ak.flatten(mp['Merged_Event'])
                                correct_merged = merged_evts & ak.flatten(mp['Merged_BHadWJa'] | mp['Merged_BHadWJb'] | mp['Merged_WJets'])
                                wrong_merged = merged_evts & ~(ak.flatten(mp['Merged_BHadWJa'] | mp['Merged_BHadWJb'] | mp['Merged_WJets']))
                        
                                    # lost events
                                lost_evts = valid_evts & ak.flatten(mp['Lost_Event'])
                                correct_lost = lost_evts & ak.flatten((mp['Lost_WJa'] | mp['Lost_WJb']))
                                wrong_lost = lost_evts & ~(ak.flatten(mp['Lost_WJa'] | mp['Lost_WJb']))
                        
                                # event categorization
                                    # unmatchable events
                                unmatchable_evts = (~valid_evts | wrong_merged | wrong_lost)
                                    # matchable events
                                matchable_evts = (correct_lost | correct_merged) # matched perm is correct event type irregardless of object matching

                            else:
                                valid_evts = (ak.num(mp['TTbar'].pt) > 0) & (ak.flatten(mp['unique_matches'] == 4))
                                    # unmatchable events
                                unmatchable_evts = ~valid_evts
                                    # matchable events
                                matchable_evts = valid_evts

                            mp_status[matchable_evts] = 2
                            mp_status[unmatchable_evts] = 3

                            output = self.make_categories(acc=output, jetmult=jmult, leptype=lepton, lepcat=leptype, btagregion=btagregion, permarray=mp_status[MTHigh], genttbar=events['SL'][cut][MTHigh], matched_perm=mp[MTHigh], evt_weights=wts)

        return output


    def make_categories(self, acc, jetmult, leptype, lepcat, btagregion, permarray, genttbar, matched_perm, evt_weights):
        for permval in np.unique(permarray).tolist():
            perm_inds = np.where(permarray == permval)
            dataset_name = '%s_%s' % (self.sample_name, perm_cats[permval]) if permval != 0 else self.sample_name

            for obj in ['TTbar', 'THad', 'TLep', 'BHad', 'BLep', 'WHad', 'WLep', 'Lepton']:
                mp_exists = ak.where(ak.num(matched_perm[obj][perm_inds].pt) > 0)[0]
                gen_tt = genttbar[obj][perm_inds][mp_exists]
                mp = matched_perm[obj][perm_inds][mp_exists]
                e_weights = evt_weights[perm_inds][mp_exists]

                acc['Gen_pt'    ].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, pt=ak.flatten(gen_tt.pt, axis=None), weight=e_weights)
                acc['Gen_eta'   ].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, eta=ak.flatten(gen_tt.eta, axis=None), weight=e_weights)
                acc['Gen_phi'   ].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, phi=ak.flatten(gen_tt.phi, axis=None), weight=e_weights)
                acc['Gen_mass'  ].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, mass=ak.flatten(gen_tt.mass, axis=None), weight=e_weights)
                acc['Gen_energy'].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, energy=ak.flatten(gen_tt.energy, axis=None), weight=e_weights)

                acc['Reco_pt'    ].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, pt=ak.flatten(mp.pt, axis=None), weight=e_weights)
                acc['Reco_eta'   ].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, eta=ak.flatten(mp.eta, axis=None), weight=e_weights)
                acc['Reco_phi'   ].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, phi=ak.flatten(mp.phi, axis=None), weight=e_weights)
                acc['Reco_mass'  ].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, mass=ak.flatten(mp.mass, axis=None), weight=e_weights)
                acc['Reco_energy'].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, energy=ak.flatten(mp.energy, axis=None), weight=e_weights)

                acc['Reso_pt'    ].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, reso_pt=ak.flatten(gen_tt.pt, axis=None)-ak.flatten(mp.pt, axis=None), weight=e_weights)
                acc['Reso_eta'   ].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, reso_eta=ak.flatten(gen_tt.eta, axis=None)-ak.flatten(mp.eta, axis=None), weight=e_weights)
                acc['Reso_phi'   ].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, reso_phi=ak.flatten(gen_tt.phi, axis=None)-ak.flatten(mp.phi, axis=None), weight=e_weights)
                acc['Reso_mass'  ].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, reso_mass=ak.flatten(gen_tt.mass, axis=None)-ak.flatten(mp.mass, axis=None), weight=e_weights)
                acc['Reso_energy'].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, reso_energy=ak.flatten(gen_tt.energy, axis=None)-ak.flatten(mp.energy, axis=None), weight=e_weights)

        return acc

    def postprocess(self, accumulator):
        return accumulator


proc_executor = processor.iterative_executor if args.debug else processor.futures_executor
proc_exec_args = {"schema": NanoAODSchema} if args.debug else {"schema": NanoAODSchema, "workers": 8}
output = processor.run_uproot_job(
    fileset,
    treename="Events",
    processor_instance=matched_perm(),
    executor=proc_executor,
    executor_args=proc_exec_args,
    #chunksize=10000 if args.debug else 100000,
    chunksize=100000,
    #chunksize=10000,
)

save(output, args.outfname)
print('%s has been written' % args.outfname)
