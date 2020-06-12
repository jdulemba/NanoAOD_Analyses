#!/usr/bin/env python

from coffea import hist
from coffea.util import save, load
import coffea.processor as processor
from pdb import set_trace
import os, sys
import python.ObjectSelection as objsel
#import coffea.processor.dataframe
import Utilities.plot_tools as plt_tools
import python.MCWeights as MCWeights
import numpy as np
import Utilities.prettyjson as prettyjson
import python.GenParticleSelector as genpsel
import python.TTGenMatcher as ttmatcher
import python.TTPermutator as ttpermutator
import Utilities.make_variables as make_vars
from python.IDJet import btag_values as btag_values

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'matched_perm'

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

## init tt probs for likelihoods
ttpermutator.year_to_run(year=args.year)

    ## run on correct samples
if args.year == '2016':
    Nominal_ttJets = ['ttJets_PS', 'ttJets']
else:
    Nominal_ttJets = ['ttJetsSL', 'ttJetsHad', 'ttJetsDiLep']
isTTbar = np.array([(key in Nominal_ttJets) for key in fileset.keys()]).all()
if not isTTbar:
    raise ValueError("This should only be run on nominal ttbar events!")

init_btag = ~(np.array([key.startswith('data') for key in fileset.keys()]).all())

## load corrections for event weights
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

# 0 == '' (no gen matching), 1 == 'right', 2 == 'matchable', 3 == 'unmatchable', 4 == 'other' (not semilep)
perm_cats = {
    0 : '',
    1 : 'right',
    2 : 'matchable',
    3 : 'unmatchable',
    4 : 'other'
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

        #set_trace()
            ## object selection
        objsel_evts = objsel.select(df, year=args.year, corrections=self.corrections, accumulator=output)
        output['cutflow']['nEvts passing jet and lepton obj selection'] += objsel_evts.sum()
        selection.add('objselection', objsel_evts)
        selection.add('jets_3', df['Jet'].counts == 3)
        selection.add('jets_4p', df['Jet'].counts > 3)
        selection.add('DeepCSV_pass', df['Jet']['DeepCSV'+wps_to_use[0]].sum() >= 2)

            # sort jets by btag value
        df['Jet'] = df['Jet'][df['Jet']['btagDeepB'].argsort(ascending=False)] if btaggers[0] == 'DeepCSV' else df['Jet'][df['Jet']['btagDeepFlavB'].argsort(ascending=False)]

            # btag fail sideband
        deepcsv_sorted = df['Jet'][df['Jet']['btagDeepB'].argsort(ascending=False)]['btagDeepB']
        valid_counts_inds = np.where(df['Jet'].counts > 1)[0]
        deepcsv_fail = np.zeros(df.size).astype(bool)
        deepcsv_fail[valid_counts_inds] = (deepcsv_sorted[valid_counts_inds][:, 0] < btag_values[args.year]['btagDeepB']['DeepCSV'+wps_to_use[0]]) & (deepcsv_sorted[valid_counts_inds][:, 1] < btag_values[args.year]['btagDeepB']['DeepCSV'+wps_to_use[0]])
        selection.add('DeepCSV_fail', deepcsv_fail ) # highest and second highest DeepCSV values don't pass tight and loose WPs

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

        if corrections['BTagSF'] == True:
            #set_trace()
            threeJets_cut = selection.require(objselection=True, jets_3=True)
            deepcsv_3j_wts = self.corrections['BTag_Constructors']['DeepCSV']['3Jets'].get_scale_factor(jets=df['Jet'][threeJets_cut], passing_cut='DeepCSV'+wps_to_use[0])
            evt_weights._weights['DeepCSV'][threeJets_cut] = deepcsv_3j_wts['central'].prod()

            fourplusJets_cut = selection.require(objselection=True, jets_4p=True)
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


            # find gen level particles for ttbar system
        #set_trace()
        genp_mode = 'NORMAL'
        GenTTbar = genpsel.select(df, systype='FINAL', mode=genp_mode)
        #genpsel.select(df, mode='LHE')
        #set_trace()
        selection.add('semilep', GenTTbar['SL']['TTbar'].counts > 0)

        ## fill hists for each region
        for lepton in regions.keys():
            lepSF_to_exclude = 'Electron_SF' if lepton == 'Muon' else 'Muon_SF'
            for leptype in regions[lepton].keys():
            #for btagger in regions[lepton].keys():
                #btagSF_to_exclude = 'DeepJet' if btagger == 'DeepCSV' else 'DeepCSV'
                btagSF_to_exclude = 'DeepJet'
                for btagregion in regions[lepton][leptype].keys():
                    for jmult in regions[lepton][leptype][btagregion].keys():
                        cut = selection.all(*regions[lepton][leptype][btagregion][jmult])

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

                                # find matched perm
                            mp = ttmatcher.best_match(gen_hyp=GenTTbar[cut], jets=df['Jet'][cut], leptons=df[lepton][cut][lep_mask], met=df['MET'][cut])

                                ## create MT regions
                            MT = make_vars.MT(df[lepton][cut][lep_mask], df['MET'][cut])
                            MTHigh = (MT >= MTcut).flatten()
                            #MTHigh = (MT[valid_bps] >= MTcut).flatten()

                            evt_weights_to_use = evt_weights.weight()
                            ## apply lepton SFs to MC (only applicable to tight leptons)
                            if 'LeptonSF' in corrections.keys():
                                evt_weights_to_use = evt_weights.partial_weight(exclude=['Electron_SF']) if lepton == 'Muon' else evt_weights.partial_weight(exclude=['Muon_SF']) # exclude SF from other lepton

                            mp_status = np.zeros(cut.sum(), dtype=int) # 0 == '' (no gen matching), 1 == 'right', 2 == 'matchable', 3 == 'unmatchable', 4 == 'noslep'
                            if jmult == '3Jets':
                                valid_evts = ((mp['TTbar'].counts > 0) & ((mp['unique_matches'] >= 3).flatten()))
                        
                                    # merged events
                                merged_evts = valid_evts & mp['Merged_Event'].flatten()
                                correct_merged = merged_evts & (mp['Merged_BHadWJa'] | mp['Merged_BHadWJb'] | mp['Merged_WJets']).flatten()
                                wrong_merged = merged_evts & ~(mp['Merged_BHadWJa'] | mp['Merged_BHadWJb'] | mp['Merged_WJets']).flatten()
                        
                                    # lost events
                                lost_evts = valid_evts & mp['Lost_Event'].flatten()
                                correct_lost = lost_evts & ((mp['Lost_WJa'] | mp['Lost_WJb'])).flatten()
                                wrong_lost = lost_evts & ~((mp['Lost_WJa'] | mp['Lost_WJb'])).flatten()
                        
                                # event categorization
                                    # unmatchable events
                                unmatchable_evts = (~valid_evts | wrong_merged | wrong_lost)
                                    # matchable events
                                matchable_evts = (correct_lost | correct_merged) # matched perm is correct event type irregardless of object matching

                            else:
                                valid_evts = ((mp['TTbar'].counts > 0) & ((mp['unique_matches'] == 4).flatten()))
                                    # unmatchable events
                                unmatchable_evts = ~valid_evts
                                    # matchable events
                                matchable_evts = valid_evts

                            #set_trace()
                            mp_status[matchable_evts] = 2
                            mp_status[unmatchable_evts] = 3

                            output = self.make_categories(acc=output, jetmult=jmult, leptype=lepton, lepcat=leptype, btagregion=btagregion, permarray=mp_status[MTHigh], genttbar=GenTTbar[cut][MTHigh], matched_perm=mp[MTHigh], evt_weights=evt_weights_to_use[cut][MTHigh])

        return output


    def make_categories(self, acc, jetmult, leptype, lepcat, btagregion, permarray, genttbar, matched_perm, evt_weights):
        for permval in np.unique(permarray).tolist():
            perm_inds = np.where(permarray == permval)
            dataset_name = '%s_%s' % (self.sample_name, perm_cats[permval]) if permval != 0 else self.sample_name
            #set_trace()

            for obj in ['TTbar', 'THad', 'TLep', 'BHad', 'BLep', 'WHad', 'WLep', 'Lepton']:
                mp_exists = np.where(matched_perm[obj][perm_inds].counts > 0)
                gen_tt = genttbar['SL'][obj][perm_inds][mp_exists]
                mp = matched_perm[obj][perm_inds][mp_exists]
                e_weights = evt_weights[perm_inds][mp_exists]
                acc['Gen_pt'    ].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, pt=gen_tt.p4.pt.flatten(), weight=e_weights)
                acc['Gen_eta'   ].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, eta=gen_tt.p4.eta.flatten(), weight=e_weights)
                acc['Gen_phi'   ].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, phi=gen_tt.p4.phi.flatten(), weight=e_weights)
                acc['Gen_mass'  ].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, mass=gen_tt.p4.mass.flatten(), weight=e_weights)
                acc['Gen_energy'].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, energy=gen_tt.p4.energy.flatten(), weight=e_weights)

                acc['Reco_pt'    ].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, pt=mp.p4.pt.flatten(), weight=e_weights)
                acc['Reco_eta'   ].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, eta=mp.p4.eta.flatten(), weight=e_weights)
                acc['Reco_phi'   ].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, phi=mp.p4.phi.flatten(), weight=e_weights)
                acc['Reco_mass'  ].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, mass=mp.p4.mass.flatten(), weight=e_weights)
                acc['Reco_energy'].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, energy=mp.p4.energy.flatten(), weight=e_weights)

                acc['Reso_pt'    ].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, reso_pt=gen_tt.p4.pt.flatten()-mp.p4.pt.flatten(), weight=e_weights)
                acc['Reso_eta'   ].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, reso_eta=gen_tt.p4.eta.flatten()-mp.p4.eta.flatten(), weight=e_weights)
                acc['Reso_phi'   ].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, reso_phi=gen_tt.p4.phi.flatten()-mp.p4.phi.flatten(), weight=e_weights)
                acc['Reso_mass'  ].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, reso_mass=gen_tt.p4.mass.flatten()-mp.p4.mass.flatten(), weight=e_weights)
                acc['Reso_energy'].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, objtype=obj, reso_energy=gen_tt.p4.energy.flatten()-mp.p4.energy.flatten(), weight=e_weights)

        return acc

    def postprocess(self, accumulator):
        return accumulator

proc_executor = processor.iterative_executor if args.debug else processor.futures_executor

output = processor.run_uproot_job(fileset,
    treename='Events',
    processor_instance=matched_perm(),
    executor=proc_executor,
    executor_args={
        'workers': 8,
        'flatten' : True,
        'compression': 5,
    },
    chunksize=10000 if args.debug else 100000,
    #chunksize=10000 if args.debug else 50000,
    #chunksize=10000,
)


if args.debug:
    print(output)
#set_trace()
#print(output['cutflow'])

save(output, args.outfname)
print('%s has been written' % args.outfname)
