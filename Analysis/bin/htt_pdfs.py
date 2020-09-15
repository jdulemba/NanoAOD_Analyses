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
analyzer = 'htt_pdfs'

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('fset', type=str, help='Fileset dictionary (in string form) to be used for the processor')
parser.add_argument('year', choices=['2016', '2017', '2018'], help='Specify which year to run over')
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
    ## specify ttJets samples
if args.year == '2016':
    Nominal_ttJets = ['ttJets_PS']#, 'ttJets']
else:
    Nominal_ttJets = ['ttJetsSL', 'ttJetsHad', 'ttJetsDiLep']

isTTbar_ = samplename.startswith('ttJets')
isNominal_ttJets_ = samplename in Nominal_ttJets
if not isNominal_ttJets_:
    raise ValueError("Only nominal ttJets samples are to be run with this analyzer.")
#isTTShift_ = isTTbar_ and (samplename.endswith('UP') or samplename.endswith('DOWN') or 'mtop' in samplename)
#isData_ = samplename.startswith('data')

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
    'BTagSF' : True,
    'JetCor' : jet_corrections,
    'Alpha' : alpha_corrections,
}

cfg_pars = prettyjson.loads(open('%s/cfg_files/cfg_pars_%s.json' % (proj_dir, jobid)).read())
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

# 0 == '' (no gen matching), 1 == 'right', 2 == 'matchable', 3 == 'unmatchable', 4 == 'other' (not semilep)
perm_cats = {
    0 : '',
    1 : 'right',
    2 : 'matchable',
    3 : 'unmatchable',
    4 : 'other'
}

# get systematics to run
event_systematics_to_run = ['nosys']
reweight_systematics_to_run = ['nosys']
print('Running with event systematics:', *sorted(set(event_systematics_to_run).difference(set(['nosys']))), sep=', ') if 'nosys' in event_systematics_to_run else print('Running with event systematics:', *sorted(event_systematics_to_run), sep=', ')
print('  and reweight systematics:', *sorted(reweight_systematics_to_run), sep=', ')
#set_trace()


# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class htt_pdfs(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.sys_axis = hist.Cat("sys", "Systematic")
        self.jetmult_axis = hist.Cat("jmult", "nJets")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")
        self.lepcat_axis = hist.Cat("lepcat", "Lepton Category")
        self.btag_axis = hist.Cat("btag", "BTag Region")
        self.pdf_axis = hist.Cat("pdfset", "PDF Set")
        #self.mtregion_axis = hist.Cat("mtregion", "MT Region")
        #self.pt_axis = hist.Bin("pt", "p_{T} [GeV]", 200, 0, 1000)
        #self.eta_axis = hist.Bin("eta", r"$\eta$", 200, -5., 5.)
        #self.phi_axis = hist.Bin("phi", r"$\phi$", 160, -4, 4)
        #self.energy_axis = hist.Bin("energy", "E [GeV]", 200, 0, 1000)
        #self.njets_axis = hist.Bin("njets", "n_{jets}", 20, 0, 20)
        #self.lepIso_axis = hist.Bin("iso", "pfRelIso", 100, 0., 1.)
        #self.mt_axis = hist.Bin("mt", "M_{T}", 200, 0., 1000.)
        #self.btagSF_axis = hist.Bin("btagSF", "SF_{btag}", 100, 0., 5.)
        #self.lepSF_axis = hist.Bin("lepSF", "SF_{lep}", 100, 0., 2.)
        #self.mtop_axis = hist.Bin("mtop", "m(top) [GeV]", 300, 0, 300)
        self.mtt_axis = hist.Bin("mtt", "m($t\overline{t}$) [GeV]", 360, 200, 2000)
        #self.probDisc_axis = hist.Bin("prob", "$\lambda_{C}$", 300, 0, 30)
        #self.massDisc_axis = hist.Bin("massdisc", "$\lambda_{M}$", 300, 0, 30)
        #self.nsDisc_axis = hist.Bin("nsdisc", "$\lambda_{NS}$", 300, 0, 30)
        #self.ctstar_axis = hist.Bin("ctstar", "cos($\\theta^{*}$)", 40, -1., 1.)
        self.ctstar_abs_axis = hist.Bin("ctstar_abs", "|cos($\\theta^{*}$)|", 20, 0., 1.)

            ## make dictionary of hists
        histo_dict = {}
                ## make selection plots
        selection_hists = self.make_selection_hists()
        histo_dict.update(selection_hists)        

        self.sample_name = ''
        self.corrections = corrections
        self.reweight_systematics_to_run = reweight_systematics_to_run
        self.event_systematics_to_run = event_systematics_to_run

            ## make dict of cutflow for each systematic variation
        for sys in self.event_systematics_to_run:
            histo_dict['cutflow_%s' % sys] = processor.defaultdict_accumulator(int)
    
        self._accumulator = processor.dict_accumulator(histo_dict)


    @property
    def accumulator(self):
        return self._accumulator


    def make_selection_hists(self):
        histo_dict = {}
        #histo_dict['mtt'] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mtt_axis)

        #histo_dict['tlep_ctstar']     = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.ctstar_axis)
        #histo_dict['tlep_ctstar_abs'] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.ctstar_abs_axis)

        histo_dict['mtt_vs_tlep_ctstar_abs'] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.pdf_axis, self.mtt_axis, self.ctstar_abs_axis)

        return histo_dict

    def process(self, df):
        np.random.seed(10) # sets seed so values from random distributions are reproducible (JER corrections)
        output = self.accumulator.identity()

        self.sample_name = df.dataset

            # run over systematics that require changes to event objects (jets+MET)
        for evt_sys in self.event_systematics_to_run:

                ## initialize regions
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

            if (evt_sys != 'nosys'):
                ## only keep signal region when event systematic isn't 'nosys'
                del regions['Muon']['Loose']
                del regions['Electron']['Loose']
                del regions['Muon']['Tight']['btagFail']
                del regions['Electron']['Tight']['btagFail']

                ## make event weights
                    # data or MC distinction made internally
            mu_evt_weights = MCWeights.get_event_weights(df, year=args.year, corrections=self.corrections, BTagSFs=btaggers, isTTbar=isNominal_ttJets_)
            el_evt_weights = MCWeights.get_event_weights(df, year=args.year, corrections=self.corrections, BTagSFs=btaggers, isTTbar=isNominal_ttJets_)

                ## initialize selections
            selection = processor.PackedSelection()

                ## object selection
            objsel_evts = objsel.select(df, year=args.year, corrections=self.corrections, cutflow=output['cutflow_%s' % evt_sys]) if evt_sys == 'nosys' else objsel.select(df, year=args.year, corrections=self.corrections, cutflow=output['cutflow_%s' % evt_sys], shift=evt_sys)
            output['cutflow_%s' % evt_sys]['nEvts passing jet and lepton obj selection'] += objsel_evts.sum()
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

            #set_trace()
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

            # don't use ttbar events with indices % 10 == 0, 1, 2
            if isNominal_ttJets_:
                events = df.event
                selection.add('keep_ttbar', ~np.stack([((events % 10) == idx) for idx in [0, 1, 2]], axis=1).any(axis=1))
                for lepton in regions.keys():
                    for leptype in regions[lepton].keys():
                        for btagregion in regions[lepton][leptype].keys():
                            for jmult in regions[lepton][leptype][btagregion].keys():
                                sel = regions[lepton][leptype][btagregion][jmult]
                                sel.update({'keep_ttbar'})

            if isTTbar_:
                genp_mode = 'NORMAL'
                GenTTbar = genpsel.select(df, systype='FINAL', mode=genp_mode)
                MCWeights.get_pdf_weights(df)
                    

            #set_trace()
            ## fill hists for each region
            for lepton in regions.keys():
                evt_weights = mu_evt_weights if lepton == 'Muon' else el_evt_weights
                for leptype in regions[lepton].keys():
                    for btagregion in regions[lepton][leptype].keys():
                        for jmult in regions[lepton][leptype][btagregion].keys():
                            cut = selection.all(*regions[lepton][leptype][btagregion][jmult])
                            #set_trace()

                            output['cutflow_%s' % evt_sys]['nEvts %s' % ', '.join([lepton, leptype, btagregion, jmult])] += cut.sum()

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
                                output['cutflow_%s' % evt_sys]['nEvts %s: valid perms' % ', '.join([lepton, leptype, btagregion, jmult])] += valid_perms.sum()

                                    ## calculate MT
                                MT = make_vars.MT(df[lepton][cut][lep_mask], df['MET'][cut])
                                MTHigh = (MT[valid_perms] >= MTcut).flatten()
                                output['cutflow_%s' % evt_sys]['nEvts %s: pass MT cut' % ', '.join([lepton, leptype, btagregion, jmult])] += MTHigh.sum()
                                if MTHigh.sum() == 0: continue

                                    # fill hists for each systematic
                                if args.debug: print('  evt sys:', evt_sys)
                                #set_trace()
                                if evt_sys == 'nosys':
                                    for rewt_sys in self.reweight_systematics_to_run:
                                        if args.debug: print('    sysname:', rewt_sys)

                                        #set_trace()
                                            ## only fill plots in signal region if systematic variation being used
                                        if (rewt_sys != 'nosys') and ('%s_%s' % (leptype, btagregion) != 'Tight_btagPass'): continue

                                        if rewt_sys == 'nosys':
                                            evt_weights_to_use = evt_weights.weight()*btag_weights['%s_CEN' % btaggers[0]]
                                        elif rewt_sys.startswith('btag'):
                                            evt_weights_to_use = evt_weights.weight()*btag_weights[rewt_sys.replace('btag', btaggers[0])]
                                        else:
                                            evt_weights_to_use = evt_weights.weight(rewt_sys)*btag_weights['%s_CEN' % btaggers[0]]

                                        sysname = rewt_sys
                                        output = self.fill_selection_hists(acc=output, sys=sysname, jetmult=jmult, leptype=lepton, lepcat=leptype, btagregion=btagregion, MTcut=MTHigh, perm=best_perms, pdf_wts=df['PDFWeights'][cut][valid_perms][MTHigh], evt_weights=evt_weights_to_use[cut][valid_perms])


                                else:
                                    if args.debug: print('    sysname:', evt_sys)
                                    evt_weights_to_use = evt_weights.weight()*btag_weights['%s_CEN' % btaggers[0]]

                                    output = self.fill_selection_hists(acc=output, sys=evt_sys, jetmult=jmult, leptype=lepton, lepcat=leptype, btagregion=btagregion, MTcut=MTHigh, perm=best_perms, pdf_wts=df['PDFWeights'][cut][valid_perms][MTHigh], evt_weights=evt_weights_to_use[cut][valid_perms])


        return output



    def fill_selection_hists(self, acc, sys, jetmult, leptype, lepcat, btagregion, MTcut, perm, pdf_wts, evt_weights):
            ## apply alpha correction for 3Jets
        if jetmult == '3Jets':
            orig_thad_p4, tlep_p4 = perm['THad'].p4.flatten()[MTcut], perm['TLep'].p4.flatten()[MTcut]
            alpha_corr = self.corrections['Alpha'](172.5/orig_thad_p4.mass)
            thad_p4 = orig_thad_p4*alpha_corr

        else:
            thad_p4, tlep_p4 = perm['THad'].p4.flatten()[MTcut], perm['TLep'].p4.flatten()[MTcut]
            
        ttbar_p4 = (thad_p4 + tlep_p4)
        thad_ctstar, tlep_ctstar = make_vars.ctstar_flat(thad_p4, tlep_p4)

        for idx in range(sorted(set(pdf_wts.counts))[0]):
            acc['mtt_vs_tlep_ctstar_abs'].fill(dataset=self.sample_name, sys=sys,  jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pdfset=str(idx), mtt=ttbar_p4.mass, ctstar_abs=np.abs(tlep_ctstar), weight=evt_weights[MTcut]*pdf_wts[:,idx])

        return acc        



    def postprocess(self, accumulator):
        return accumulator

proc_executor = processor.iterative_executor if args.debug else processor.futures_executor

output = processor.run_uproot_job(fileset,
    treename='Events',
    processor_instance=htt_pdfs(),
    executor=proc_executor,
    executor_args={
        'workers': 8,
        'flatten' : True,
        #'compression': 8,
        'compression': 5,
    },
    chunksize=10000 if args.debug else 100000,
    #chunksize=10000 if args.debug else 500000,
)


if args.debug:
    print(output)

save(output, args.outfname)
print('%s has been written' % args.outfname)
