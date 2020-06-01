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

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'btagging_invest'

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
class btagging_invest(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.jetmult_axis = hist.Cat("jmult", "nJets")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")
        self.lepcat_axis = hist.Cat("lepcat", "Lepton Category")
        self.btag_axis = hist.Cat("btag", "btagging Category")
        self.pt_axis = hist.Bin("pt", "p_{T} [GeV]", 200, 0, 1000)
        self.njets_axis = hist.Bin("njets", "n_{jets}", 20, 0, 20)
        self.mt_axis = hist.Bin("mt", "M_{T}", 200, 0., 1000.)
        self.SF_axis = hist.Bin("sf", "SF", 500, 0., 5.)
        self.deepcsv_axis = hist.Bin("deepcsv", "SF", 100, 0., 1.)
        self.pu_axis = hist.Bin("pu", "nTrueInt", 100, 0., 100.)

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
        histo_dict['Jets_LeadJet_pt']= hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepcat_axis, self.btag_axis, self.pt_axis)
        histo_dict['Jets_pt']        = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepcat_axis, self.btag_axis, self.pt_axis)
        histo_dict['Jets_njets']     = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepcat_axis, self.btag_axis, self.njets_axis)
        histo_dict['Jets_DeepCSV']   = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepcat_axis, self.btag_axis, self.deepcsv_axis)
        histo_dict['Lep_pt']         = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepcat_axis, self.btag_axis, self.pt_axis)
        histo_dict['MT']             = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepcat_axis, self.btag_axis, self.mt_axis)

        histo_dict['nTrueInt_puweight'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepcat_axis, self.btag_axis, self.pu_axis)
        histo_dict['nTrueInt_noweight'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepcat_axis, self.btag_axis, self.pu_axis)
        histo_dict['BTagSF'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepcat_axis, self.btag_axis, self.SF_axis)
        histo_dict['LepSF'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepcat_axis, self.btag_axis, self.SF_axis)
        histo_dict['PileupWeight'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepcat_axis, self.btag_axis, self.SF_axis)
        histo_dict['EvtWeight'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepcat_axis, self.btag_axis, self.SF_axis)

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
                    'zero_btags' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'loose_MU', 'DeepCSV_0'},
                        '4PJets' : {'objselection', 'jets_4p', 'loose_MU', 'DeepCSV_0'},
                    },
                    'one_btag' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'loose_MU', 'DeepCSV_1'},
                        '4PJets' : {'objselection', 'jets_4p', 'loose_MU', 'DeepCSV_1'},
                    },
                    'two_btags' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'loose_MU', 'DeepCSV_2'},
                        '4PJets' : {'objselection', 'jets_4p', 'loose_MU', 'DeepCSV_2'},
                    },
                    'threePlus_btags' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'loose_MU', 'DeepCSV_3p'},
                        '4PJets' : {'objselection', 'jets_4p', 'loose_MU', 'DeepCSV_3p'},
                    },
                },
                'Tight' : {
                    'zero_btags' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'tight_MU', 'DeepCSV_0'},
                        '4PJets' : {'objselection', 'jets_4p', 'tight_MU', 'DeepCSV_0'},
                    },
                    'one_btag' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'tight_MU', 'DeepCSV_1'},
                        '4PJets' : {'objselection', 'jets_4p', 'tight_MU', 'DeepCSV_1'},
                    },
                    'two_btags' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'tight_MU', 'DeepCSV_2'},
                        '4PJets' : {'objselection', 'jets_4p', 'tight_MU', 'DeepCSV_2'},
                    },
                    'threePlus_btags' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'tight_MU', 'DeepCSV_3p'},
                        '4PJets' : {'objselection', 'jets_4p', 'tight_MU', 'DeepCSV_3p'},
                    },
                },
            },
            'Electron' : {
                'Loose' : {
                    'zero_btags' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'loose_EL', 'DeepCSV_0'},
                        '4PJets' : {'objselection', 'jets_4p', 'loose_EL', 'DeepCSV_0'},
                    },
                    'one_btag' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'loose_EL', 'DeepCSV_1'},
                        '4PJets' : {'objselection', 'jets_4p', 'loose_EL', 'DeepCSV_1'},
                    },
                    'two_btags' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'loose_EL', 'DeepCSV_2'},
                        '4PJets' : {'objselection', 'jets_4p', 'loose_EL', 'DeepCSV_2'},
                    },
                    'threePlus_btags' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'loose_EL', 'DeepCSV_3p'},
                        '4PJets' : {'objselection', 'jets_4p', 'loose_EL', 'DeepCSV_3p'},
                    },
                },
                'Tight' : {
                    'zero_btags' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'tight_EL', 'DeepCSV_0'},
                        '4PJets' : {'objselection', 'jets_4p', 'tight_EL', 'DeepCSV_0'},
                    },
                    'one_btag' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'tight_EL', 'DeepCSV_1'},
                        '4PJets' : {'objselection', 'jets_4p', 'tight_EL', 'DeepCSV_1'},
                    },
                    'two_btags' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'tight_EL', 'DeepCSV_2'},
                        '4PJets' : {'objselection', 'jets_4p', 'tight_EL', 'DeepCSV_2'},
                    },
                    'threePlus_btags' : {
                        '3Jets'  : {'objselection', 'jets_3' , 'tight_EL', 'DeepCSV_3p'},
                        '4PJets' : {'objselection', 'jets_4p', 'tight_EL', 'DeepCSV_3p'},
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
        #selection.add('DeepJet_pass', df['Jet']['DeepJet'+wps_to_use[0]].sum() >= 2)            
        #selection.add('DeepCSV_pass', df['Jet']['DeepCSV'+wps_to_use[0]].sum() >= 2)            
        selection.add('DeepCSV_0', df['Jet']['DeepCSV'+wps_to_use[0]].sum() == 0)            
        selection.add('DeepCSV_1', df['Jet']['DeepCSV'+wps_to_use[0]].sum() == 1)            
        selection.add('DeepCSV_2', df['Jet']['DeepCSV'+wps_to_use[0]].sum() == 2)            
        selection.add('DeepCSV_3p', df['Jet']['DeepCSV'+wps_to_use[0]].sum() >= 3)            

        #set_trace()
            # sort jets by btag value, needed when making permutations
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
                #selection.add('loose_or_tight_MU', (df['Muon']['LOOSEMU'] | df['Muon']['TIGHTMU']).sum() == 1) # one muon passing LOOSE or TIGHT criteria
            if isSE_Data:
                del regions['Muon']
                        ## electrons
                selection.add('tight_EL', df['Electron']['TIGHTEL'].sum() == 1) # one electron passing TIGHT criteria
                selection.add('loose_EL', df['Electron']['LOOSEEL'].sum() == 1) # one electron passing LOOSE criteria
                #selection.add('loose_or_tight_EL', (df['Electron']['LOOSEEL'] | df['Electron']['TIGHTEL']).sum() == 1) # one electron passing LOOSE or TIGHT criteria

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
            #selection.add('loose_or_tight_MU', (df['Muon']['LOOSEMU'] | df['Muon']['TIGHTMU']).sum() == 1) # one muon passing LOOSE or TIGHT criteria
                    ## electrons
            selection.add('tight_EL', df['Electron']['TIGHTEL'].sum() == 1) # one electron passing TIGHT criteria
            selection.add('loose_EL', df['Electron']['LOOSEEL'].sum() == 1) # one electron passing LOOSE criteria
            #selection.add('loose_or_tight_EL', (df['Electron']['LOOSEEL'] | df['Electron']['TIGHTEL']).sum() == 1) # one electron passing LOOSE or TIGHT criteria

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

                                ## calculate MT
                            MT = make_vars.MT(df[lepton][cut][lep_mask], df['MET'][cut])
                            MTHigh = (MT >= MTcut).flatten()

                            evt_weights_to_use = evt_weights.weight()
                            if not self.isData:
                                evt_weights_to_use = evt_weights.partial_weight(exclude=[lepSF_to_exclude, btagSF_to_exclude])

                            jets = df['Jet'][cut][MTHigh]
                            leptons = df[lepton][cut][lep_mask][MTHigh]

                            btagSF = np.ones(MTHigh.size) if self.isData else evt_weights._weights[btaggers[0]][cut][MTHigh].flatten()
                            lepSF = np.ones(MTHigh.size) if self.isData else evt_weights._weights['%s_SF' % lepton][cut][MTHigh].flatten()
                            pu_weight = np.ones(MTHigh.size) if self.isData else evt_weights._weights['pileup_weight'][cut][MTHigh].flatten()
                            tot_weight = evt_weights_to_use[cut][MTHigh].flatten()
                            #set_trace()
                            if not self.isData:
                                output['nTrueInt_puweight'].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, lepcat=lepcat, btag=btagregion, pu=df['Pileup_nTrueInt'][cut][MTHigh], weight=pu_weight)
                                output['nTrueInt_noweight'].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, lepcat=lepcat, btag=btagregion, pu=df['Pileup_nTrueInt'][cut][MTHigh])
                            output['BTagSF'].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, lepcat=lepcat, btag=btagregion, sf=btagSF)
                            output['LepSF'].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, lepcat=lepcat, btag=btagregion, sf=lepSF)
                            output['PileupWeight'].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, lepcat=lepcat, btag=btagregion, sf=pu_weight)
                            output['EvtWeight'].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, lepcat=lepcat, btag=btagregion, sf=tot_weight)

                            output = self.fill_hists(acc=output, jetmult=jmult, leptype=lepton, lepcat=lepcat, btag=btagregion, jets=jets, leptons=leptons, MT=MT[MTHigh].flatten(), evt_weights=tot_weight)


        return output

    def fill_hists(self, acc, jetmult, leptype, lepcat, btag, jets, leptons, MT, evt_weights):
        #set_trace()
        acc['Jets_njets'].fill( dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btag, njets=jets.counts, weight=evt_weights)
        acc['Jets_LeadJet_pt'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btag, pt=jets.pt.max(), weight=evt_weights)
        acc['Jets_pt'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btag, pt=jets.pt.flatten(), weight=(jets.pt.ones_like()*evt_weights).flatten())
        acc['Jets_DeepCSV'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btag, deepcsv=jets.btagDeepB.flatten(), weight=(jets.pt.ones_like()*evt_weights).flatten())

        acc['Lep_pt'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btag, pt=leptons.pt.flatten(), weight=evt_weights)

        acc['MT'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btag, mt=MT, weight=evt_weights)

        return acc        


    def postprocess(self, accumulator):
        return accumulator

proc_executor = processor.iterative_executor if args.debug else processor.futures_executor

output = processor.run_uproot_job(fileset,
    treename='Events',
    processor_instance=btagging_invest(),
    executor=proc_executor,
    executor_args={
        'workers': 8,
        'flatten' : True,
        'compression': 5,
    },
    chunksize=10000 if args.debug else 100000,
    #chunksize=10000 if args.debug else 50000,
)


if args.debug:
    print(output)
#set_trace()
#print(output['cutflow'])

save(output, args.outfname)
print('%s has been written' % args.outfname)
