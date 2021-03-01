#!/usr/bin/env python

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
import numpy as np
import Utilities.prettyjson as prettyjson
import coffea.lumi_tools.lumi_tools as lumi_tools
import Utilities.make_variables as make_vars
from python.IDJet import btag_values as btag_values
import python.GenParticleSelector as genpsel
import python.TTGenMatcher as ttmatcher
import python.TTPermutator as ttpermutator
from python.Permutations import compare_matched_best_perms

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
base_jobid = os.environ['base_jobid']
analyzer = 'best_perms'

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

isTTbar_ = samplename.startswith('ttJets')
isTTSL_ = samplename.startswith('ttJetsSL')
isData_ = samplename.startswith('data')
if isData_:
    lumiMask_path = os.path.join(proj_dir, 'inputs', 'data', base_jobid, 'LumiMasks', '%s_GoldenJson_%s.txt' % (args.year, base_jobid))

## init tt probs for likelihoods
ttpermutator.year_to_run(year=args.year)

## load lumimask for data and corrections for event weights
pu_correction = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'MC_PU_Weights.coffea'))[args.year]
lepSF_correction = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'leptonSFs.coffea'))[args.year]
jet_corrections = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'JetMETCorrections.coffea'))[args.year]
alpha_corrections = load(os.path.join(proj_dir, 'Corrections', jobid,'alpha_correction_%s.coffea' % jobid))[args.year]['E']['All_2D'] # E, All_2D determined by post application plots
corrections = {
    'Pileup' : pu_correction,
    'Prefire' : True,
    'LeptonSF' : lepSF_correction,
    'BTagSF' : not isData_,
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

if corrections['BTagSF'] == True:
    sf_file = '%s/Corrections/%s/%s' % (proj_dir, jobid, jet_pars['btagging']['btagSF_file'])
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
class best_perms(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.jetmult_axis = hist.Cat("jmult", "nJets")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")
        self.lepcat_axis = hist.Cat("lepcat", "Lepton Category")
        self.btag_axis = hist.Cat("btag", "BTag Region")
        self.objtype_axis = hist.Cat("objtype", "Object Type")
        self.pt_axis = hist.Bin("pt", "p_{T} [GeV]", 200, 0, 1000)
        self.eta_axis = hist.Bin("eta", r"$\eta$", 200, -5., 5.)
        self.phi_axis = hist.Bin("phi", r"$\phi$", 160, -4, 4)
        self.energy_axis = hist.Bin("energy", "E [GeV]", 200, 0, 1000)
        self.mass_axis = hist.Bin("mass", "m [GeV]", 1000, 0., 2000.)
        self.njets_axis = hist.Bin("njets", "n_{jets}", 20, 0, 20)
        self.lepIso_axis = hist.Bin("iso", "pfRelIso", 100, 0., 1.)
        self.mt_axis = hist.Bin("mt", "M_{T}", 200, 0., 1000.)
        self.mtop_axis = hist.Bin("mtop", "m(top) [GeV]", 300, 0, 300)
        self.mtt_axis = hist.Bin("mtt", "m($t\overline{t}$) [GeV]", 180, 200, 2000)
        self.probDisc_axis = hist.Bin("prob", "$\lambda_{C}$", 300, 0, 30)
        self.massDisc_axis = hist.Bin("massdisc", "$\lambda_{M}$", 300, 0, 30)
        self.nsDisc_axis = hist.Bin("nsdisc", "$\lambda_{NS}$", 300, 0, 30)
        self.ctstar_axis = hist.Bin("ctstar", "cos($\\theta^{*}$)", 200, -1., 1.)
        self.ctstar_abs_axis = hist.Bin("ctstar_abs", "|cos($\\theta^{*}$)|", 200, 0., 1.)

            ## make dictionary of hists
        histo_dict = {}
                ## make jet hists
        jet_hists = self.make_jet_hists()
        histo_dict.update(jet_hists)
                ## make selection plots
        selection_hists = self.make_selection_hists()
        histo_dict.update(selection_hists)        
                ## make bp plots
        bp_hists = self.make_bp_hists()
        histo_dict.update(bp_hists)        

        histo_dict['cutflow'] = processor.defaultdict_accumulator(int)

        self._accumulator = processor.dict_accumulator(histo_dict)
        self.sample_name = ''
        self.corrections = corrections
        self.isData = True

    
    @property
    def accumulator(self):
        return self._accumulator


    def make_bp_hists(self):
        histo_dict = {}
        histo_dict['BP_pt']     = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.objtype_axis, self.pt_axis)
        histo_dict['BP_eta']    = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.objtype_axis, self.eta_axis)
        histo_dict['BP_phi']    = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.objtype_axis, self.phi_axis)
        histo_dict['BP_mass']   = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.objtype_axis, self.mass_axis)
        histo_dict['BP_energy'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.objtype_axis, self.energy_axis)

        return histo_dict


    def make_jet_hists(self):
        histo_dict = {}
        histo_dict['Jets_pt']    = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.pt_axis)
        histo_dict['Jets_eta']   = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.eta_axis)
    #    histo_dict['Jets_phi']   = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.phi_axis)
    #    histo_dict['Jets_energy']= hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.energy_axis)
        histo_dict['Jets_njets'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.njets_axis)
        histo_dict['Jets_LeadJet_pt']    = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.pt_axis)
        histo_dict['Jets_LeadJet_eta']   = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.eta_axis)
    #    histo_dict['Jets_LeadJet_phi']   = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.phi_axis)
    #    histo_dict['Jets_LeadJet_energy']= hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.energy_axis)

        return histo_dict

    def make_selection_hists(self):
        histo_dict = {}
        histo_dict['mtt']      = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mtt_axis)
        histo_dict['mthad']    = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mtop_axis)
        histo_dict['pt_thad']  = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.pt_axis)
        histo_dict['pt_tt']    = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.pt_axis)
        histo_dict['eta_thad'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.eta_axis)
        histo_dict['eta_tt']   = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.eta_axis)

        histo_dict['full_disc'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.probDisc_axis)
        histo_dict['mass_disc'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.massDisc_axis)
        histo_dict['ns_disc']   = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.nsDisc_axis)

        histo_dict['tlep_ctstar']     = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.ctstar_axis)
        histo_dict['tlep_ctstar_abs'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.ctstar_abs_axis)

        histo_dict['MT'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mt_axis)

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
                        '3Jets'  : {'objselection', 'jets_3' , 'tight_MU', 'DeepCSV_pass'},
                        '4PJets' : {'objselection', 'jets_4p', 'tight_MU', 'DeepCSV_pass'},
                    },
                    #'btagFail' : {
                    #    '3Jets'  : {'objselection', 'jets_3' , 'tight_MU', 'DeepCSV_fail'},
                    #    '4PJets' : {'objselection', 'jets_4p', 'tight_MU', 'DeepCSV_fail'},
                    #},
                },
                #'Loose' : {
                #    'btagPass' : {
                #        '3Jets'  : {'objselection', 'jets_3' , 'loose_MU', 'DeepCSV_pass'},
                #        '4PJets' : {'objselection', 'jets_4p', 'loose_MU', 'DeepCSV_pass'},
                #    },
                #    'btagFail' : {
                #        '3Jets'  : {'objselection', 'jets_3' , 'loose_MU', 'DeepCSV_fail'},
                #        '4PJets' : {'objselection', 'jets_4p', 'loose_MU', 'DeepCSV_fail'},
                #    },
                #},
            },
            'Electron' : {
                'Tight' : {
                    'btagPass' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'tight_EL', 'DeepCSV_pass'},
                        '4PJets' : {'objselection', 'jets_4p', 'tight_EL', 'DeepCSV_pass'},
                    },
                    #'btagFail' : {
                    #    '3Jets'  : {'objselection', 'jets_3' , 'tight_EL', 'DeepCSV_fail'},
                    #    '4PJets' : {'objselection', 'jets_4p', 'tight_EL', 'DeepCSV_fail'},
                    #},
                },
                #'Loose' : {
                #    'btagPass' : {
                #        '3Jets'  : {'objselection', 'jets_3' , 'loose_EL', 'DeepCSV_pass'},
                #        '4PJets' : {'objselection', 'jets_4p', 'loose_EL', 'DeepCSV_pass'},
                #    },
                #    'btagFail' : {
                #        '3Jets'  : {'objselection', 'jets_3' , 'loose_EL', 'DeepCSV_fail'},
                #        '4PJets' : {'objselection', 'jets_4p', 'loose_EL', 'DeepCSV_fail'},
                #    },
                #},
            },
        }

            ## object selection
        objsel_evts = objsel.select(events, year=args.year, corrections=self.corrections, cutflow=output['cutflow'])
        #set_trace()
        output['cutflow']['nEvts passing jet and lepton obj selection'] += ak.sum(objsel_evts)
        selection.add('objselection', objsel_evts)
        selection.add('jets_3',  ak.num(events['Jet']) == 3)
        selection.add('jets_4p', ak.num(events['Jet']) > 3)
        selection.add('DeepCSV_pass', ak.sum(events['Jet'][btag_wps[0]], axis=1) >= 2)

            # sort jets by btag value
        events['Jet'] = events['Jet'][ak.argsort(events['Jet']['btagDeepB'], ascending=False)] if btaggers[0] == 'DeepCSV' else events['Jet'][ak.argsort(events['Jet']['btagDeepFlavB'], ascending=False)]

            # btag fail sideband
        deepcsv_sorted = events['Jet'][ak.argsort(events['Jet']['btagDeepB'], ascending=False)]['btagDeepB']
        valid_counts_inds = ak.where(ak.num(events['Jet']) > 1)[0]
        deepcsv_fail = np.zeros(len(events)).astype(bool)
        deepcsv_fail[valid_counts_inds] = (deepcsv_sorted[valid_counts_inds][:, 0] < btag_values[args.year]['btagDeepB']['DeepCSV'+wps_to_use[0]]) & (deepcsv_sorted[valid_counts_inds][:, 1] < btag_values[args.year]['btagDeepB']['DeepCSV'+wps_to_use[0]])
        selection.add('DeepCSV_fail', deepcsv_fail ) # highest and second highest DeepCSV values don't pass tight and loose WPs

        #set_trace()
        if isData_:
            isSE_Data = self.sample_name.startswith('data_SingleElectron')
            isSM_Data = self.sample_name.startswith('data_SingleMuon')
            runs = events.run
            lumis = events.luminosityBlock
            Golden_Json_LumiMask = lumi_tools.LumiMask(lumiMask_path)
            LumiMask = Golden_Json_LumiMask.__call__(runs, lumis) ## returns array of valid events
            selection.add('lumimask', LumiMask)
   
                ## object selection and add different selections
            if isSM_Data:
                del regions['Electron']
                        ## muons
                selection.add('tight_MU', ak.sum(events['Muon']['TIGHTMU'], axis=1) == 1) # one muon passing TIGHT criteria
                selection.add('loose_MU', ak.sum(events['Muon']['LOOSEMU'], axis=1) == 1) # one muon passing LOOSE criteria
            if isSE_Data:
                del regions['Muon']
                        ## electrons
                selection.add('tight_EL', ak.sum(events['Electron']['TIGHTEL'], axis=1) == 1) # one electron passing TIGHT criteria
                selection.add('loose_EL', ak.sum(events['Electron']['LOOSEEL'], axis=1) == 1) # one electron passing LOOSE criteria

            for lepton in regions.keys():
                for leptype in regions[lepton].keys():
                    for btagregion in regions[lepton][leptype].keys():
                        for jmult in regions[lepton][leptype][btagregion].keys():
                            regions[lepton][leptype][btagregion][jmult].update({'lumimask'})

        else:
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

                ## apply btagging SFs to MC
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


        if isTTbar_:
                # find gen level particles for ttbar system
            genpsel.select(events, mode='NORMAL')
            selection.add('semilep', ak.num(events['SL']) > 0)
            if isTTSL_:
                for lepton in regions.keys():
                    for leptype in regions[lepton].keys():
                        for btagregion in regions[lepton][leptype].keys():
                            for jmult in regions[lepton][leptype][btagregion].keys():
                                regions[lepton][leptype][btagregion][jmult].update({'semilep'})


        ## fill hists for each region
        for lepton in regions.keys():
            evt_weights = mu_evt_weights if lepton == 'Muon' else el_evt_weights
            for leptype in regions[lepton].keys():
                for btagregion in regions[lepton][leptype].keys():
                    for jmult in regions[lepton][leptype][btagregion].keys():
                        cut = selection.all(*regions[lepton][leptype][btagregion][jmult])
                        #set_trace()

                        #if args.debug: print(lepton, leptype, btagregion, jmult)
                        if cut.sum() > 0:
                                # get leptons
                            ltype = 'MU' if lepton == 'Muon' else 'EL'
                            if 'loose_or_tight_%s' % ltype in regions[lepton][leptype][btagregion][jmult]:
                                leptons = events[lepton][cut][((events[lepton][cut]['TIGHT%s' % ltype] == True) | (events[lepton][cut]['LOOSE%s' % ltype] == True))]
                            elif 'tight_%s' % ltype in regions[lepton][leptype][btagregion][jmult]:
                                leptons = events[lepton][cut][(events[lepton][cut]['TIGHT%s' % ltype] == True)]
                            elif 'loose_%s' % ltype in regions[lepton][leptype][btagregion][jmult]:
                                leptons = events[lepton][cut][(events[lepton][cut]['LOOSE%s' % ltype] == True)]
                            else:
                                raise ValueError("Not sure what lepton type to choose for event")

                                # get jets and MET
                            jets, met = events['Jet'][cut], events['MET'][cut]

                                # find best permutations
                            best_perms = ttpermutator.find_best_permutations(jets=jets, leptons=leptons, MET=met, btagWP=btag_wps[0], btag_req=False if btagregion == 'btagFail' else True)
                            valid_perms = ak.num(best_perms['TTbar'].pt) > 0

                            bp_status = np.zeros(cut.size, dtype=int) # 0 == '' (no gen matching), 1 == 'right', 2 == 'matchable', 3 == 'unmatchable', 4 == 'sl_tau', 5 == 'noslep'
                                # get matched permutation (semilep ttbar only)
                            if isTTbar_:
                                semilep_evts = selection.require(semilep=True)
                                bp_status[~semilep_evts] = 5
                                if semilep_evts.sum() > 0:
                                        # find matched permutations
                                    mp = ttmatcher.best_match(gen_hyp=events['SL'][cut], jets=jets, leptons=leptons, met=met)
                                    perm_cat_array = compare_matched_best_perms(mp, best_perms, njets=jmult)
                                    bp_status[cut] = perm_cat_array
                                    sl_tau_evts = ak.where(ak.fill_none(ak.pad_none(np.abs(events['SL']['Lepton'].pdgId) == 15, 1), False) == True)[0]
                                    bp_status[sl_tau_evts] = 4


                                ## create MT regions
                            MT = make_vars.MT(leptons, met)
                            MTHigh = ak.flatten(MT[valid_perms] >= MTcut)

                            wts = evt_weights.weight()[cut][valid_perms][MTHigh] if self.isData else (evt_weights.weight()*btag_weights['%s_CEN' % btaggers[0]])[cut][valid_perms][MTHigh]
                                # fill hists
                            output = self.fill_selection_hists(acc=output, jetmult=jmult, leptype=lepton, lepcat=leptype, btagregion=btagregion, permarray=bp_status[cut][valid_perms][MTHigh],
                                perm=best_perms[valid_perms][MTHigh], MTvals=MT[valid_perms][MTHigh], evt_wts=wts)
                            #set_trace()
                            output = self.fill_jet_hists(acc=output, jetmult=jmult, leptype=lepton, lepcat=leptype, btagregion=btagregion,
                                permarray=bp_status[cut][valid_perms][MTHigh], obj=jets[valid_perms][MTHigh], evt_wts=wts)


        return output



    def fill_selection_hists(self, acc, jetmult, leptype, lepcat, btagregion, permarray, perm, MTvals, evt_wts):
        #set_trace()
            ## apply alpha correction for 3Jets
        if jetmult == '3Jets':
            alpha_corr = self.corrections['Alpha'](172.5/perm['THad'].mass)
            perm['THad'] = perm['THad'].multiply(alpha_corr) # correct thad
            perm['TTbar'] = ak.flatten(perm['THad']+perm['TLep']) # correct ttbar

        thad_ctstar, tlep_ctstar = make_vars.ctstar(perm['THad'], perm['TLep'])
        thad_ctstar, tlep_ctstar = ak.flatten(thad_ctstar, axis=None), ak.flatten(tlep_ctstar, axis=None)

        for permval in np.unique(permarray).tolist():
            perm_inds = np.where(permarray == permval)
            dataset_name = '%s_%s' % (self.sample_name, perm_cats[permval]) if permval != 0 else self.sample_name

            for obj in ['TTbar', 'THad', 'TLep', 'BHad', 'BLep', 'WHad', 'WLep', 'Lepton']:
                perm_obj = perm[obj]
                acc['BP_pt'    ].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, pt=ak.flatten(perm_obj.pt[perm_inds]), weight=evt_wts[perm_inds])
                acc['BP_eta'   ].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, eta=ak.flatten(perm_obj.eta[perm_inds]), weight=evt_wts[perm_inds])
                acc['BP_phi'   ].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, phi=ak.flatten(perm_obj.phi[perm_inds]), weight=evt_wts[perm_inds])
                acc['BP_mass'  ].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, mass=ak.flatten(perm_obj.mass[perm_inds]), weight=evt_wts[perm_inds])
                acc['BP_energy'].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, energy=ak.flatten(perm_obj.energy[perm_inds]), weight=evt_wts[perm_inds])

            acc['mtt'].fill(     dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtt=ak.flatten(perm['TTbar'].mass)[perm_inds], weight=evt_wts[perm_inds])
            acc['mthad'].fill(   dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtop=ak.flatten(perm['THad'].mass)[perm_inds], weight=evt_wts[perm_inds])
            acc['pt_thad'].fill( dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=ak.flatten(perm['THad'].pt)[perm_inds], weight=evt_wts[perm_inds])
            acc['pt_tt'].fill(   dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=ak.flatten(perm['TTbar'].pt)[perm_inds], weight=evt_wts[perm_inds])
            acc['eta_thad'].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=ak.flatten(perm['THad'].eta)[perm_inds], weight=evt_wts[perm_inds])
            acc['eta_tt'].fill(  dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=ak.flatten(perm['TTbar'].eta)[perm_inds], weight=evt_wts[perm_inds])

            acc['tlep_ctstar'].fill(    dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, ctstar=tlep_ctstar[perm_inds], weight=evt_wts[perm_inds])
            acc['tlep_ctstar_abs'].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, ctstar_abs=np.abs(tlep_ctstar[perm_inds]), weight=evt_wts[perm_inds])

            acc['full_disc'].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, prob=ak.flatten(perm['Prob'])[perm_inds], weight=evt_wts[perm_inds])
            acc['mass_disc'].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, massdisc=ak.flatten(perm['MassDiscr'])[perm_inds], weight=evt_wts[perm_inds])
            acc['ns_disc'].fill(  dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, nsdisc=ak.flatten(perm['NuDiscr'])[perm_inds], weight=evt_wts[perm_inds])
            acc['MT'].fill(       dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mt=ak.flatten(MTvals)[perm_inds], weight=evt_wts[perm_inds])

        return acc        

    def fill_jet_hists(self, acc, jetmult, leptype, lepcat, btagregion, permarray, obj, evt_wts):
        #set_trace()
        pt_sorted_jets = obj[ak.argsort(obj.pt, ascending=False)]
        for permval in np.unique(permarray).tolist():
            perm_inds = np.where(permarray == permval)
            dataset_name = '%s_%s' % (self.sample_name, perm_cats[permval]) if permval != 0 else self.sample_name

            acc['Jets_pt'].fill(    dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=ak.flatten(obj.pt[perm_inds]), weight=ak.flatten((ak.ones_like(obj.pt)*evt_wts)[perm_inds]))
            acc['Jets_eta'].fill(   dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=ak.flatten(obj.eta[perm_inds]), weight=ak.flatten((ak.ones_like(obj.eta)*evt_wts)[perm_inds]))
            acc['Jets_njets'].fill( dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, njets=ak.num(obj)[perm_inds], weight=evt_wts[perm_inds])
            #acc['Jets_phi'].fill(   dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, phi=obj.phi.flatten(), weight=(obj.pt.ones_like()*evt_wts).flatten())
            #acc['Jets_energy'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, energy=obj.p4.E.flatten(), weight=(obj.pt.ones_like()*evt_wts).flatten())

            acc['Jets_LeadJet_pt'].fill(    dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=pt_sorted_jets[perm_inds].pt[:, 0], weight=evt_wts[perm_inds])
            acc['Jets_LeadJet_eta'].fill(   dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=pt_sorted_jets[perm_inds].eta[:, 0], weight=evt_wts[perm_inds])
            #acc['Jets_LeadJet_phi'].fill(   dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, phi=pt_sorted_jets.phi[:, 0], weight=evt_wts)
            #acc['Jets_LeadJet_energy'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, energy=pt_sorted_jets.p4.E[:, 0], weight=evt_wts)

        return acc        



    def postprocess(self, accumulator):
        return accumulator


proc_executor = processor.iterative_executor if args.debug else processor.futures_executor
proc_exec_args = {"schema": NanoAODSchema} if args.debug else {"schema": NanoAODSchema, "workers": 8}
output = processor.run_uproot_job(
    fileset,
    treename="Events",
    processor_instance=best_perms(),
    executor=proc_executor,
    executor_args=proc_exec_args,
    chunksize=10000 if args.debug else 100000,
    #chunksize=100000,
)

save(output, args.outfname)
print('%s has been written' % args.outfname)
