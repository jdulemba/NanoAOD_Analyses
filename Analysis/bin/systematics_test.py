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
analyzer = 'systematics_test'

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

init_btag = ~(np.array([key.startswith('data') for key in fileset.keys()]).all())

## init tt probs for likelihoods
ttpermutator.year_to_run(year=args.year)

## load lumimask for data and corrections for event weights
pu_correction = load('%s/Corrections/%s/MC_PU_Weights.coffea' % (proj_dir, jobid))
lepSF_correction = load('%s/Corrections/leptonSFs.coffea' % proj_dir)
jet_corrections = load('%s/Corrections/JetCorrections.coffea' % proj_dir)[args.year]
alpha_corrections = load('%s/Corrections/%s/alpha_correction_%s.coffea' % (proj_dir, jobid, jobid))[args.year]['E']['All'] # E, All determined by post application plots
corrections = {
    'Pileup' : pu_correction,
    'Prefire' : True,
    'LeptonSF' : lepSF_correction,
    'BTagSF' : init_btag,
    'JetCor' : jet_corrections,
    'Alpha' : alpha_corrections,
}

cfg_pars = prettyjson.loads(open('%s/cfg_files/cfg_pars_%s.json' % (proj_dir, jobid)).read())
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

    ## specify ttJets samples
if args.year == '2016':
    Nominal_ttJets = ['ttJets_PS', 'ttJets']
else:
    Nominal_ttJets = ['ttJetsSL', 'ttJetsHad', 'ttJetsDiLep']


# 0 == '' (no gen matching), 1 == 'right', 2 == 'matchable', 3 == 'unmatchable', 4 == 'other' (not semilep)
perm_cats = {
    0 : '',
    1 : 'right',
    2 : 'matchable',
    3 : 'unmatchable',
    4 : 'other'
}

allowed_systematics = ['nosys']
if init_btag: # MC
    allowed_systematics += systematics.lep_sys + systematics.misc_sys + systematics.btag_sys



# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class systematics_test(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.sys_axis = hist.Cat("sys", "Systematic")
        self.jetmult_axis = hist.Cat("jmult", "nJets")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")
        self.lepcat_axis = hist.Cat("lepcat", "Lepton Category")
        self.btag_axis = hist.Cat("btag", "BTag Region")
        #self.mtregion_axis = hist.Cat("mtregion", "MT Region")
        self.pt_axis = hist.Bin("pt", "p_{T} [GeV]", 200, 0, 1000)
        self.eta_axis = hist.Bin("eta", r"$\eta$", 200, -5., 5.)
        #self.phi_axis = hist.Bin("phi", r"$\phi$", 160, -4, 4)
        #self.energy_axis = hist.Bin("energy", "E [GeV]", 200, 0, 1000)
        self.njets_axis = hist.Bin("njets", "n_{jets}", 20, 0, 20)
        self.lepIso_axis = hist.Bin("iso", "pfRelIso", 100, 0., 1.)
        self.mt_axis = hist.Bin("mt", "M_{T}", 200, 0., 1000.)
        #self.btagSF_axis = hist.Bin("btagSF", "SF_{btag}", 100, 0., 5.)
        #self.lepSF_axis = hist.Bin("lepSF", "SF_{lep}", 100, 0., 2.)
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
                ## make lepton hists
        lep_hists = self.make_lep_hists()
        histo_dict.update(lep_hists)        
        #        ## make other hists
        #other_hists = self.make_other_hists()
        #histo_dict.update(other_hists)        
        #        ## make selection plots
        #selection_hists = self.make_selection_hists()
        #histo_dict.update(selection_hists)        

        histo_dict['cutflow'] = processor.defaultdict_accumulator(int)

        self._accumulator = processor.dict_accumulator(histo_dict)
        self.sample_name = ''
        self.corrections = corrections
        self.isData = True

    
    @property
    def accumulator(self):
        return self._accumulator


    def make_jet_hists(self):
        histo_dict = {}
        histo_dict['Jets_pt']        = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.pt_axis)
        histo_dict['Jets_njets']     = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.njets_axis)
        histo_dict['Jets_LeadJet_pt']= hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.pt_axis)

        return histo_dict

    
    def make_lep_hists(self):
        histo_dict = {}
        histo_dict['Lep_pt']    = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.pt_axis)

        return histo_dict

    #def make_other_hists(self):
    #    histo_dict = {}
    #    for btagger in btaggers:
    #        histo_dict['%s_SF' % btagger] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mtregion_axis, self.btagSF_axis)
    #    histo_dict['Lep_SF']  = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mtregion_axis, self.lepSF_axis)
    #    if corrections['Prefire'] == True:
    #        histo_dict['prefire_weight']  = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mtregion_axis, self.lepSF_axis)
    #    if corrections['Pileup'] is not None:
    #        histo_dict['pileup_weight']  = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mtregion_axis, self.lepSF_axis)

    #    return histo_dict

    #def make_selection_hists(self):
    #    histo_dict = {}
    #    histo_dict['mtt']      = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mtt_axis)
    #    histo_dict['mthad']    = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mtop_axis)
    #    histo_dict['pt_thad']  = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.pt_axis)
    #    histo_dict['pt_tlep']  = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.pt_axis)
    #    histo_dict['pt_tt']    = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.pt_axis)
    #    histo_dict['eta_thad'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.eta_axis)
    #    histo_dict['eta_tlep'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.eta_axis)
    #    histo_dict['eta_tt']   = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.eta_axis)

    #    histo_dict['full_disc'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.probDisc_axis)
    #    histo_dict['mass_disc'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.massDisc_axis)
    #    histo_dict['ns_disc']   = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.nsDisc_axis)

    #    histo_dict['tlep_ctstar']     = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.ctstar_axis)
    #    histo_dict['tlep_ctstar_abs'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.ctstar_abs_axis)

    #    histo_dict['MT'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mt_axis)

    #    return histo_dict

    def process(self, df):
        np.random.seed(10) # sets seed so values from random distributions are reproducible (JER corrections)
        output = self.accumulator.identity()

        self.sample_name = df.dataset
        isTTbar = self.sample_name.startswith('ttJets')

            ## make event weights
                # data or MC distinction made internally
        mu_evt_weights = MCWeights.get_event_weights(df, year=args.year, corrections=self.corrections, BTagSFs=btaggers)
        el_evt_weights = MCWeights.get_event_weights(df, year=args.year, corrections=self.corrections, BTagSFs=btaggers)

            ## initialize selections and regions
        selection = processor.PackedSelection()
        regions = {
            'Muon' : {
                'Tight' : {
                    'btagPass' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'tight_MU', 'DeepCSV_pass'},
                        '4PJets' : {'objselection', 'jets_4p', 'tight_MU', 'DeepCSV_pass'},
                    },
                },
            },
            'Electron' : {
                'Tight' : {
                    'btagPass' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'tight_EL', 'DeepCSV_pass'},
                        '4PJets' : {'objselection', 'jets_4p', 'tight_EL', 'DeepCSV_pass'},
                    },
                },
            },
        }

            ## object selection
        objsel_evts = objsel.select(df, year=args.year, corrections=self.corrections, accumulator=output)
        output['cutflow']['nEvts passing jet and lepton obj selection'] += objsel_evts.sum()
        selection.add('jets_3', df['Jet'].counts == 3)
        selection.add('jets_4p', df['Jet'].counts > 3)
        selection.add('objselection', objsel_evts)
        selection.add('DeepCSV_pass', df['Jet']['DeepCSV'+wps_to_use[0]].sum() >= 2)

            # sort jets by btag value
        df['Jet'] = df['Jet'][df['Jet']['btagDeepB'].argsort(ascending=False)] if btaggers[0] == 'DeepCSV' else df['Jet'][df['Jet']['btagDeepFlavB'].argsort(ascending=False)]

            # btag fail sideband
        deepcsv_sorted = df['Jet'][df['Jet']['btagDeepB'].argsort(ascending=False)]['btagDeepB']
        valid_counts_inds = np.where(df['Jet'].counts > 1)[0]
        deepcsv_fail = np.zeros(df.size).astype(bool)
        deepcsv_fail[valid_counts_inds] = (deepcsv_sorted[valid_counts_inds][:, 0] < btag_values[args.year]['btagDeepB']['DeepCSV'+wps_to_use[0]]) & (deepcsv_sorted[valid_counts_inds][:, 1] < btag_values[args.year]['btagDeepB']['DeepCSV'+wps_to_use[0]])
        selection.add('DeepCSV_fail', deepcsv_fail ) # highest and second highest DeepCSV values don't pass tight and loose WPs

        #set_trace()
        self.isData = self.sample_name.startswith('data_Single')
        if self.isData:
            isSE_Data = self.sample_name.startswith('data_SingleElectron')
            isSM_Data = self.sample_name.startswith('data_SingleMuon')
            runs = df.run
            lumis = df.luminosityBlock
            Golden_Json_LumiMask = lumi_tools.LumiMask('%s/inputs/data/LumiMasks/%s_GoldenJson.txt' % (proj_dir, args.year))
            LumiMask = Golden_Json_LumiMask.__call__(runs, lumis) ## returns array of valid events
            selection.add('lumimask', LumiMask)
   
                ## object selection and add different selections
            if isSM_Data:
                del regions['Electron']
                        ## muons
                selection.add('tight_MU', df['Muon']['TIGHTMU'].sum() == 1) # one muon passing TIGHT criteria
                selection.add('loose_MU', df['Muon']['LOOSEMU'].sum() == 1) # one muon passing LOOSE criteria
            if isSE_Data:
                del regions['Muon']
                        ## electrons
                selection.add('tight_EL', df['Electron']['TIGHTEL'].sum() == 1) # one electron passing TIGHT criteria
                selection.add('loose_EL', df['Electron']['LOOSEEL'].sum() == 1) # one electron passing LOOSE criteria

            for lepton in regions.keys():
                #for btagger in regions[lepton].keys():
                for leptype in regions[lepton].keys():
                    for btagregion in regions[lepton][leptype].keys():
                        for jmult in regions[lepton][leptype][btagregion].keys():
                            regions[lepton][leptype][btagregion][jmult].update({'lumimask'})

        if not self.isData:
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
                #set_trace()

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

            # don't use ttbar events with indices % 10 == 0, 1, 2
            if self.sample_name in Nominal_ttJets:
                events = df.event
                selection.add('keep_ttbar', ~np.stack([((events % 10) == idx) for idx in [0, 1, 2]], axis=1).any(axis=1))
                for lepton in regions.keys():
                    for leptype in regions[lepton].keys():
                        for btagregion in regions[lepton][leptype].keys():
                            for jmult in regions[lepton][leptype][btagregion].keys():
                                sel = regions[lepton][leptype][btagregion][jmult]
                                sel.update({'keep_ttbar'})

            if isTTbar:
                genp_mode = 'NORMAL'
                GenTTbar = genpsel.select(df, systype='FINAL', mode=genp_mode)


        #set_trace()
        ## fill hists for each region
        for lepton in regions.keys():
            evt_weights = mu_evt_weights if lepton == 'Muon' else el_evt_weights
            for leptype in regions[lepton].keys():
                for btagregion in regions[lepton][leptype].keys():
                    for jmult in regions[lepton][leptype][btagregion].keys():
                        cut = selection.all(*regions[lepton][leptype][btagregion][jmult])
                        #set_trace()

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

                            bp_status = np.zeros(cut.size, dtype=int) # 0 == '' (no gen matching), 1 == 'right', 2 == 'matchable', 3 == 'unmatchable', 4 == 'noslep'
                                # get matched permutation (semilep ttbar only)
                            if isTTbar:
                                semilep_evts = GenTTbar['SL']['TTbar'].counts > 0
                                bp_status[~semilep_evts] = 4
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


                                ## calculate MT
                            MT = make_vars.MT(df[lepton][cut][lep_mask], df['MET'][cut])
                            MTHigh = (MT[valid_perms] >= MTcut).flatten()

                                # fill hists for each systematic
                            #set_trace()
                            for sys in allowed_systematics:
                                if args.debug: print('    ', sys)
                            #for sys in ['nosys']+sorted(evt_weights.variations):
                            #    if sys not in allowed_systematics:
                            #        raise ValueError("Systematic %s not in list of allowed systematics to run" % sys)

                                if sys == 'nosys':
                                    evt_weights_to_use = evt_weights.weight() if self.isData else evt_weights.weight()*btag_weights['%s_CEN' % btaggers[0]]
                                elif sys.startswith('btag'):
                                    evt_weights_to_use = evt_weights.weight()*btag_weights[sys.replace('btag', btaggers[0])]
                                else:
                                    evt_weights_to_use = evt_weights.weight(sys)*btag_weights['%s_CEN' % btaggers[0]]

                                #output = self.fill_selection_hists(acc=output, jetmult=jmult, leptype=lepton, lepcat=leptype, btagregion=btagregion, permarray=bp_status[cut][valid_perms][MTHigh], MTcut=MTHigh, perm=best_perms, evt_weights=evt_weights_to_use[cut][valid_perms])
                                output = self.fill_jet_hists(acc=output, sys=sys, jetmult=jmult, leptype=lepton, lepcat=leptype, btagregion=btagregion, permarray=bp_status[cut][valid_perms][MTHigh], obj=df['Jet'][cut][valid_perms][MTHigh], evt_weights=evt_weights_to_use[cut][valid_perms][MTHigh])
                                output = self.fill_lep_hists(acc=output, sys=sys, jetmult=jmult, leptype=lepton, lepcat=leptype, btagregion=btagregion, permarray=bp_status[cut][valid_perms][MTHigh], obj=df[lepton][cut][lep_mask][valid_perms][MTHigh], evt_weights=evt_weights_to_use[cut][valid_perms][MTHigh])
                                ##output = self.fill_sf_hists(acc=output, jetmult=jmult, leptype=lepton, lepcat=leptype, btagregion=btagregion, evt_weights=evt_weights._weights if not self.isData else evt_weights_to_use, mask=cut , mtcut=MTHigh)

                                #set_trace()
                                #for permval in np.unique(bp_status[cut][valid_perms][MTHigh]).tolist():
                                #    #set_trace()
                                #    perm_inds = np.where(bp_status[cut][valid_perms][MTHigh] == permval)
                                #    dataset_name = '%s_%s' % (self.sample_name, perm_cats[permval]) if permval != 0 else self.sample_name
                                #    output['MT'].fill(dataset=dataset_name, jmult=jmult, leptype=lepton, lepcat=leptype, btag=btagregion, mt=MT[valid_perms][MTHigh].flatten()[perm_inds], weight=evt_weights_to_use[cut][valid_perms][MTHigh][perm_inds])


        return output



    #def fill_selection_hists(self, acc, jetmult, leptype, lepcat, btagregion, permarray, MTcut, perm, evt_weights):
    #        ## apply alpha correction for 3Jets
    #    if jetmult == '3Jets':
    #        orig_thad_p4, tlep_p4 = perm['THad'].p4.flatten()[MTcut], perm['TLep'].p4.flatten()[MTcut]
    #        alpha_corr = self.corrections['Alpha'](172.5/orig_thad_p4.mass)
    #        thad_p4 = orig_thad_p4*alpha_corr

    #    else:
    #        thad_p4, tlep_p4 = perm['THad'].p4.flatten()[MTcut], perm['TLep'].p4.flatten()[MTcut]
    #        
    #    ttbar_p4 = (thad_p4 + tlep_p4)
    #    thad_ctstar, tlep_ctstar = make_vars.ctstar_flat(thad_p4, tlep_p4)

    #    for permval in np.unique(permarray).tolist():
    #        #set_trace()
    #        perm_inds = np.where(permarray == permval)
    #        dataset_name = '%s_%s' % (self.sample_name, perm_cats[permval]) if permval != 0 else self.sample_name

    #        acc['mtt'].fill(     dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtt=ttbar_p4.mass[perm_inds], weight=evt_weights[MTcut][perm_inds])
    #        acc['mthad'].fill(   dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtop=thad_p4.mass[perm_inds], weight=evt_weights[MTcut][perm_inds])
    #        acc['pt_thad'].fill( dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=thad_p4.pt[perm_inds], weight=evt_weights[MTcut][perm_inds])
    #        acc['pt_tlep'].fill( dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=tlep_p4.pt[perm_inds], weight=evt_weights[MTcut][perm_inds])
    #        acc['pt_tt'].fill(   dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=ttbar_p4.pt[perm_inds], weight=evt_weights[MTcut][perm_inds])

    #        acc['tlep_ctstar'].fill(    dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, ctstar=tlep_ctstar[perm_inds], weight=evt_weights[MTcut][perm_inds])
    #        acc['tlep_ctstar_abs'].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, ctstar_abs=np.abs(tlep_ctstar[perm_inds]), weight=evt_weights[MTcut][perm_inds])

    #        acc['full_disc'].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, prob=perm['Prob'].flatten()[MTcut][perm_inds], weight=evt_weights[MTcut][perm_inds])
    #        acc['mass_disc'].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, massdisc=perm['MassDiscr'].flatten()[MTcut][perm_inds], weight=evt_weights[MTcut][perm_inds])
    #        acc['ns_disc'].fill(  dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, nsdisc=perm['NuDiscr'].flatten()[MTcut][perm_inds], weight=evt_weights[MTcut][perm_inds])

    #    return acc        

    def fill_jet_hists(self, acc, sys, jetmult, leptype, lepcat, btagregion, permarray, obj, evt_weights):
        #set_trace()
        for permval in np.unique(permarray).tolist():
            perm_inds = np.where(permarray == permval)
            dataset_name = '%s_%s' % (self.sample_name, perm_cats[permval]) if permval != 0 else self.sample_name
            acc['Jets_pt'].fill(    dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=obj.pt.flatten()[perm_inds], weight=(obj.pt.ones_like()*evt_weights).flatten()[perm_inds])
            acc['Jets_njets'].fill( dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, njets=obj.counts[perm_inds], weight=evt_weights[perm_inds])

            pt_sorted_jets = obj[obj.pt.argsort(ascending=False)]
            acc['Jets_LeadJet_pt'].fill(    dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=pt_sorted_jets.pt[:, 0][perm_inds], weight=evt_weights[perm_inds])

        return acc        

    def fill_lep_hists(self, acc, sys, jetmult, leptype, lepcat, btagregion, permarray, obj, evt_weights):
        #set_trace()
        for permval in np.unique(permarray).tolist():
            perm_inds = np.where(permarray == permval)
            dataset_name = '%s_%s' % (self.sample_name, perm_cats[permval]) if permval != 0 else self.sample_name
            acc['Lep_pt'].fill(    dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=obj.pt.flatten()[perm_inds], weight=evt_weights[perm_inds])

        return acc        

    #def fill_sf_hists(self, acc, jetmult, leptype, lepcat, btagregion, mtregion, evt_weights, mask, mtcut):
    #    #set_trace()
    #    if not self.isData:
    #        for btagger in btaggers:
    #            acc['%s_SF' % btagger].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtregion=mtregion, btagSF=evt_weights[btagger][mask][mtcut])
    #        acc['Lep_SF'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtregion=mtregion, lepSF=evt_weights['%s_SF' % leptype][mask][mtcut])
    #        if 'pileup_weight' in acc.keys(): 
    #            acc['pileup_weight'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtregion=mtregion, lepSF=evt_weights['pileup_weight'][mask][mtcut])
    #        if 'prefire_weight' in acc.keys(): 
    #            acc['prefire_weight'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtregion=mtregion, lepSF=evt_weights['prefire_weight'][mask][mtcut])
    #    else:
    #        for btagger in btaggers:
    #            acc['%s_SF' % btagger].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtregion=mtregion, btagSF=evt_weights[mask][mtcut])
    #        acc['Lep_SF' ].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtregion=mtregion, lepSF=evt_weights[mask][mtcut])
    #        if 'pileup_weight' in acc.keys(): 
    #            acc['pileup_weight'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtregion=mtregion, lepSF=evt_weights[mask][mtcut])
    #        if 'prefire_weight' in acc.keys(): 
    #            acc['prefire_weight'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtregion=mtregion, lepSF=evt_weights[mask][mtcut])

    #    return acc


    def postprocess(self, accumulator):
        return accumulator

proc_executor = processor.iterative_executor if args.debug else processor.futures_executor

output = processor.run_uproot_job(fileset,
    treename='Events',
    processor_instance=systematics_test(),
    executor=proc_executor,
    executor_args={
        'workers': 8,
        'flatten' : True,
        'compression': 5,
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
