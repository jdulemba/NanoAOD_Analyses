#!/usr/bin/env python

from coffea import hist
from coffea.util import save, load
import coffea.processor as processor
from pdb import set_trace
import os, sys
import python.ObjectSelection as objsel
import python.MCWeights as MCWeights
import numpy as np
import Utilities.prettyjson as prettyjson
import coffea.lumi_tools.lumi_tools as lumi_tools
import Utilities.make_variables as make_vars
from python.IDJet import btag_values as btag_values
#import python.GenParticleSelector as genpsel
#import python.TTGenMatcher as ttmatcher
import python.TTPermutator as ttpermutator
import Utilities.make_variables as make_vars

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'htt_test_of_homogeneity'

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
corrections = {
    'Pileup' : pu_correction,
    'Prefire' : True,
    'LeptonSF' : lepSF_correction,
    'BTagSF' : init_btag,
    'JetCor' : jet_corrections,
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


# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class htt_test_of_homogeneity(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.jetmult_axis = hist.Cat("jmult", "nJets")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")
        self.lepcat_axis = hist.Cat("lepcat", "Lepton Category")
        self.btag_axis = hist.Cat("btag", "btagging Category")
        self.pt_axis = hist.Bin("pt", "p_{T} [GeV]", 200, 0, 1000)
        self.njets_axis = hist.Bin("njets", "n_{jets}", 20, 0, 20)
        self.lepIso_axis = hist.Bin("iso", "pfRelIso", 2000, 0., 20.)
        self.mtt_axis = hist.Bin("mtt", "m($t\overline{t}$) [GeV]", 180, 200, 2000)
        self.ctstar_axis = hist.Bin("ctstar", "cos($\\theta^{*}$)", 200, -1., 1.)

            ## make dictionary of hists
        histo_dict = {}
                ## make jet hists
        hists = self.make_hists()
        histo_dict.update(hists)

        histo_dict['cutflow'] = processor.defaultdict_accumulator(int)

        self._accumulator = processor.dict_accumulator(histo_dict)
        self.sample_name = ''
        self.corrections = corrections
        self.isData = True

    
    @property
    def accumulator(self):
        return self._accumulator


    def make_hists(self):
        histo_dict = {}
        histo_dict['mtt_vs_iso']= hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepcat_axis, self.btag_axis, self.mtt_axis, self.lepIso_axis)
        histo_dict['tlep_ctstar_vs_iso'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepcat_axis, self.btag_axis, self.ctstar_axis, self.lepIso_axis)

        return histo_dict


    def process(self, df):
        np.random.seed(10) # sets seed so values from random distributions are reproducible (JER corrections)
        output = self.accumulator.identity()

        self.sample_name = df.dataset

            ## make event weights
                # data or MC distinction made internally
        evt_weights = MCWeights.get_event_weights(df, year=args.year, corrections=self.corrections, BTagSFs=btaggers)

            ## initialize selections and regions
        selection = processor.PackedSelection()
        regions = {
            'Muon' : {
                'Loose' : {
                    'btagPass' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'loose_MU', 'DeepCSV_pass'},
                        '4PJets' : {'objselection', 'jets_4p', 'loose_MU', 'DeepCSV_pass'},
                    },
                    #'btagFail' : {
                    #    '3Jets'  : {'objselection', 'jets_3' , 'loose_MU', 'DeepCSV_fail'},
                    #    '4PJets' : {'objselection', 'jets_4p', 'loose_MU', 'DeepCSV_fail'},
                    #},
                },
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
            },
            'Electron' : {
                'Loose' : {
                    'btagPass' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'loose_EL', 'DeepCSV_pass'},
                        '4PJets' : {'objselection', 'jets_4p', 'loose_EL', 'DeepCSV_pass'},
                    },
                    #'btagFail' : {
                    #    '3Jets'  : {'objselection', 'jets_3' , 'loose_EL', 'DeepCSV_fail'},
                    #    '4PJets' : {'objselection', 'jets_4p', 'loose_EL', 'DeepCSV_fail'},
                    #},
                },
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
            },
        }

            ## object selection
        objsel_evts = objsel.select(df, year=args.year, corrections=self.corrections, accumulator=output)
        output['cutflow']['nEvts passing jet and lepton obj selection'] += objsel_evts.sum()
        selection.add('jets_3', df['Jet'].counts == 3)
        selection.add('jets_4p', df['Jet'].counts > 3)
        selection.add('objselection', objsel_evts)
        #selection.add('DeepJet_pass', df['Jet']['DeepJet'+wps_to_use[0]].sum() >= 2)            
        selection.add('DeepCSV_pass', df['Jet']['DeepCSV'+wps_to_use[0]].sum() >= 2)            

        #set_trace()
            # sort jets by btag value, needed for whenever permutations are made
        df['Jet'] = df['Jet'][df['Jet']['btagDeepB'].argsort(ascending=False)] if btaggers[0] == 'DeepCSV' else df['Jet'][df['Jet']['btagDeepFlavB'].argsort(ascending=False)]

        #    # btag fail sideband
        #deepcsv_sorted = df['Jet'][df['Jet']['btagDeepB'].argsort(ascending=False)]['btagDeepB']
        #valid_counts_inds = np.where(df['Jet'].counts > 1)[0]
        #deepcsv_fail = np.zeros(df.size).astype(bool)
        #deepcsv_fail[valid_counts_inds] = (deepcsv_sorted[valid_counts_inds][:, 0] < btag_values[args.year]['btagDeepB']['DeepCSV'+wps_to_use[0]]) & (deepcsv_sorted[valid_counts_inds][:, 1] < btag_values[args.year]['btagDeepB']['DeepCSV'+wps_to_use[0]])
        #selection.add('DeepCSV_fail', deepcsv_fail ) # highest and second highest DeepCSV values don't pass tight and loose WPs


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
                for lepcat in regions[lepton].keys():
                    for btagregion in regions[lepton][lepcat].keys():
                        for jmult in regions[lepton][lepcat][btagregion].keys():
                            regions[lepton][lepcat][btagregion][jmult].update({'lumimask'})

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
                evt_weights._weights['Muon_SF'][tight_mu_cut] = MCWeights.get_lepton_sf(year=args.year, lepton='Muons', corrections=lepSF_correction,
                    pt=tight_muons.pt.flatten(), eta=tight_muons.eta.flatten())
                tight_el_cut = selection.require(objselection=True, tight_EL=True) # find events passing electron object selection with one tight electron
                tight_electrons = df['Electron'][tight_el_cut][(df['Electron'][tight_el_cut]['TIGHTEL'] == True)]
                evt_weights._weights['Electron_SF'][tight_el_cut] = MCWeights.get_lepton_sf(year=args.year, lepton='Electrons', corrections=lepSF_correction,
                    pt=tight_electrons.pt.flatten(), eta=tight_electrons.etaSC.flatten())

                ## apply btagging SFs to MC
            if corrections['BTagSF'] == True:
                #set_trace()
                threeJets_cut = selection.require(objselection=True, jets_3=True)
                #deepjet_3j_wts = self.corrections['BTag_Constructors']['DeepJet']['3Jets'].get_scale_factor(jets=df['Jet'][threeJets_cut], passing_cut='DeepJet'+wps_to_use[0])
                #evt_weights._weights['DeepJet'][threeJets_cut] = deepjet_3j_wts['central'].prod()
                deepcsv_3j_wts = self.corrections['BTag_Constructors']['DeepCSV']['3Jets'].get_scale_factor(jets=df['Jet'][threeJets_cut], passing_cut='DeepCSV'+wps_to_use[0])
                evt_weights._weights['DeepCSV'][threeJets_cut] = deepcsv_3j_wts['central'].prod()

                fourplusJets_cut = selection.require(objselection=True, jets_4p=True)
                #deepjet_4pj_wts = self.corrections['BTag_Constructors']['DeepJet']['4PJets'].get_scale_factor(jets=df['Jet'][fourplusJets_cut], passing_cut='DeepJet'+wps_to_use[0])
                #evt_weights._weights['DeepJet'][fourplusJets_cut] = deepjet_4pj_wts['central'].prod()
                deepcsv_4pj_wts = self.corrections['BTag_Constructors']['DeepCSV']['4PJets'].get_scale_factor(jets=df['Jet'][fourplusJets_cut], passing_cut='DeepCSV'+wps_to_use[0])
                evt_weights._weights['DeepCSV'][fourplusJets_cut] = deepcsv_4pj_wts['central'].prod()

            # don't use ttbar events with indices % 10 == 0, 1, 2
            if self.sample_name in Nominal_ttJets:
                events = df.event
                selection.add('keep_ttbar', ~np.stack([((events % 10) == idx) for idx in [0, 1, 2]], axis=1).any(axis=1))
                for lepton in regions.keys():
                    for lepcat in regions[lepton].keys():
                        for btagregion in regions[lepton][lepcat].keys():
                            for jmult in regions[lepton][lepcat][btagregion].keys():
                                sel = regions[lepton][lepcat][btagregion][jmult]
                                sel.update({'keep_ttbar'})


        #set_trace()
        ## fill hists for each region
        for lepton in regions.keys():
            lepSF_to_exclude = 'Electron_SF' if lepton == 'Muon' else 'Muon_SF'
            btagSF_to_exclude = 'DeepJet'
            for lepcat in regions[lepton].keys():
                for btagregion in regions[lepton][lepcat].keys():
                    for jmult in regions[lepton][lepcat][btagregion].keys():
                        cut = selection.all(*regions[lepton][lepcat][btagregion][jmult])
                        #set_trace()

                        if cut.sum() > 0:
                            ltype = 'MU' if lepton == 'Muon' else 'EL'
                            if 'loose_or_tight_%s' % ltype in regions[lepton][lepcat][btagregion][jmult]:
                                lep_mask = ((df[lepton][cut]['TIGHT%s' % ltype] == True) | (df[lepton][cut]['LOOSE%s' % ltype] == True))
                            elif 'tight_%s' % ltype in regions[lepton][lepcat][btagregion][jmult]:
                                lep_mask = (df[lepton][cut]['TIGHT%s' % ltype] == True)
                            elif 'loose_%s' % ltype in regions[lepton][lepcat][btagregion][jmult]:
                                lep_mask = (df[lepton][cut]['LOOSE%s' % ltype] == True)
                            else:
                                raise ValueError("Not sure what lepton type to choose for event")

                                # find best permutations
                            best_perms = ttpermutator.find_best_permutations(jets=df['Jet'][cut], leptons=df[lepton][cut][lep_mask], MET=df['MET'][cut], btagWP=btag_wps[0])
                            valid_perms = best_perms['TTbar'].counts > 0

                                ## calculate MT
                            MT = make_vars.MT(df[lepton][cut][lep_mask], df['MET'][cut])
                            MTHigh = (MT[valid_perms] >= MTcut).flatten()

                            evt_weights_to_use = evt_weights.weight()
                            if not self.isData:
                                evt_weights_to_use = evt_weights.partial_weight(exclude=[lepSF_to_exclude, btagSF_to_exclude])

                            #set_trace()
                            lepIso = df[lepton][cut][lep_mask][valid_perms][MTHigh].pfRelIso.flatten()
                            mtt = best_perms['TTbar'].mass.flatten()
                            thad_ctstar, tlep_ctstar = make_vars.ctstar(best_perms['THad'].p4, best_perms['TLep'].p4)

                            e_weights = evt_weights_to_use[cut][valid_perms][MTHigh]
                            output['mtt_vs_iso'].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, lepcat=lepcat, btag=btagregion, mtt=mtt[MTHigh], iso=lepIso, weight=e_weights)
                            output['tlep_ctstar_vs_iso'].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, lepcat=lepcat, btag=btagregion, ctstar=tlep_ctstar[MTHigh], iso=lepIso, weight=e_weights)
                            #set_trace()

        return output



    def postprocess(self, accumulator):
        return accumulator

proc_executor = processor.iterative_executor if args.debug else processor.futures_executor

output = processor.run_uproot_job(fileset,
    treename='Events',
    processor_instance=htt_test_of_homogeneity(),
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
