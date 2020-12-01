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
import python.GenParticleSelector as genpsel
import python.TTGenMatcher as ttmatcher
import Utilities.make_variables as make_vars
import python.Test_Permutator as tperm

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'permProbComputer'

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('fset', type=str, help='Fileset dictionary (in string form) to be used for the processor')
parser.add_argument('year', choices=['2016APV', '2016', '2017', '2018'], help='Specify which year to run over')
parser.add_argument('outfname', type=str, help='Specify output filename, including directory and file extension')
parser.add_argument('--debug', action='store_true', help='Uses iterative_executor for debugging purposes, otherwise futures_excutor will be used (faster)')

args = parser.parse_args()

# convert input string of fileset dictionary to actual dictionary
fdict = (args.fset).replace("\'", "\"")
fileset = prettyjson.loads(fdict)

    ## run on correct samples
Nominal_ttJets = ['ttJetsSL']
isTTbar = np.array([(key in Nominal_ttJets) for key in fileset.keys()]).all()
if not isTTbar:
    raise ValueError("This should only be run on nominal ttbar events!")

## load corrections for event weights
pu_correction = load(os.path.join(proj_dir, 'Corrections', jobid, 'MC_PU_Weights.coffea'))
lepSF_correction = load(os.path.join(proj_dir, 'Corrections', jobid, 'leptonSFs.coffea'))
jet_corrections = load(os.path.join(proj_dir, 'Corrections', jobid, 'JetCorrections.coffea'))[args.year]
corrections = {
    'Pileup' : pu_correction,
    'Prefire' : True,
    'LeptonSF' : lepSF_correction,
    'BTagSF' : False,
    'JetCor' : jet_corrections,
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
        self.mbpjet_axis = hist.Bin("mbpjet", "m(b+j) [GeV]", 1000, 0., 2000.)

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
        histo_dict['Merged_mbpjet_vs_maxmjet'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepcat_axis, self.mtregion_axis, self.permcat_axis, self.maxmjet_axis, self.mbpjet_axis)
        histo_dict['Merged_nusolver_chi2']     = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepcat_axis, self.mtregion_axis, self.permcat_axis, self.nu_chi2_axis)
        histo_dict['Merged_nusolver_dist']     = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepcat_axis, self.mtregion_axis, self.permcat_axis, self.nu_dist_axis)

        return histo_dict

    def make_3j_lost_hists(self):
        histo_dict = {}
        histo_dict['Lost_mbpjet']        = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepcat_axis, self.mtregion_axis, self.permcat_axis, self.mbpjet_axis)
        histo_dict['Lost_nusolver_chi2'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepcat_axis, self.mtregion_axis, self.permcat_axis, self.nu_chi2_axis)
        histo_dict['Lost_nusolver_dist'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepcat_axis, self.mtregion_axis, self.permcat_axis, self.nu_dist_axis)

        return histo_dict


    def make_4pj_hists(self):
        histo_dict = {}
        histo_dict['mWHad_vs_mTHad'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepcat_axis, self.mtregion_axis, self.permcat_axis, self.tMass_axis, self.wMass_axis)
        histo_dict['nusolver_chi2']  = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepcat_axis, self.mtregion_axis, self.permcat_axis, self.nu_chi2_axis)
        histo_dict['nusolver_dist']  = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepcat_axis, self.mtregion_axis, self.permcat_axis, self.nu_dist_axis)

        return histo_dict

    
    def process(self, df):
        np.random.seed(10) # sets seed so values from random distributions are reproducible (JER corrections)
        output = self.accumulator.identity()

        self.sample_name = df.dataset

            ## make event weights
                # data or MC distinction made internally
        evt_weights = MCWeights.get_event_weights(df, year=args.year, corrections=self.corrections)

            ## initialize selections and regions
        selection = processor.PackedSelection()
        regions = {
            'Muon' : {
                'LoT' : {
                    '3Jets'  : {'objselection', 'jets_3' , 'loose_or_tight_MU', 'btag_pass', 'semilep'},
                    '4PJets' : {'objselection', 'jets_4p', 'loose_or_tight_MU', 'btag_pass', 'semilep'},
                    #'4Jets' :  {'objselection', 'jets_4',  'loose_or_tight_MU', 'btag_pass', 'semilep'},
                    #'5Jets' :  {'objselection', 'jets_5',  'loose_or_tight_MU', 'btag_pass', 'semilep'},
                    #'6PJets' : {'objselection', 'jets_6p', 'loose_or_tight_MU', 'btag_pass', 'semilep'},
                },
                'Tight' : {
                    '3Jets'  : {'objselection', 'jets_3' , 'tight_MU', 'btag_pass', 'semilep'},
                    '4PJets' : {'objselection', 'jets_4p', 'tight_MU', 'btag_pass', 'semilep'},
                #    '4Jets' :  {'objselection', 'jets_4',  'tight_MU', 'btag_pass', 'semilep'},
                #    '5Jets' :  {'objselection', 'jets_5',  'tight_MU', 'btag_pass', 'semilep'},
                #    '6PJets' : {'objselection', 'jets_6p', 'tight_MU', 'btag_pass', 'semilep'},
                },
            },
            'Electron' : {
                'LoT' : {
                    '3Jets'  : {'objselection', 'jets_3' , 'loose_or_tight_EL', 'btag_pass', 'semilep'},
                    '4PJets' : {'objselection', 'jets_4p', 'loose_or_tight_EL', 'btag_pass', 'semilep'},
                    #'4Jets' :  {'objselection', 'jets_4',  'loose_or_tight_EL', 'btag_pass', 'semilep'},
                    #'5Jets' :  {'objselection', 'jets_5',  'loose_or_tight_EL', 'btag_pass', 'semilep'},
                    #'6PJets' : {'objselection', 'jets_6p', 'loose_or_tight_EL', 'btag_pass', 'semilep'},
                },
                'Tight' : {
                    '3Jets'  : {'objselection', 'jets_3' , 'tight_EL', 'btag_pass', 'semilep'},
                    '4PJets' : {'objselection', 'jets_4p', 'tight_EL', 'btag_pass', 'semilep'},
                #    '4Jets' :  {'objselection', 'jets_4',  'tight_EL', 'btag_pass', 'semilep'},
                #    '5Jets' :  {'objselection', 'jets_5',  'tight_EL', 'btag_pass', 'semilep'},
                #    '6PJets' : {'objselection', 'jets_6p', 'tight_EL', 'btag_pass', 'semilep'},
                },
            },
        }

            ## object selection
        objsel_evts = objsel.select(df, year=args.year, corrections=self.corrections, cutflow=output['cutflow'])
        output['cutflow']['nEvts passing jet and lepton obj selection'] += objsel_evts.sum()
        selection.add('jets_3', df['Jet'].counts == 3)
        selection.add('jets_4p', df['Jet'].counts > 3)
        selection.add('jets_4', df['Jet'].counts == 4)
        selection.add('jets_5', df['Jet'].counts == 5)
        selection.add('jets_6p', df['Jet'].counts >= 6)
        selection.add('objselection', objsel_evts)
        selection.add('btag_pass', df['Jet'][btag_wp].sum() >= 2)

            # sort jets by btag value
        df['Jet'] = df['Jet'][df['Jet']['btagDeepB'].argsort(ascending=False)] if btagger == 'DeepCSV' else df['Jet'][df['Jet']['btagDeepFlavB'].argsort(ascending=False)]

            ## add different selections
                ## muons
        selection.add('tight_MU', df['Muon']['TIGHTMU'].sum() == 1) # one muon passing TIGHT criteria
        #selection.add('loose_MU', df['Muon']['LOOSEMU'].sum() == 1) # one muon passing LOOSE criteria
        selection.add('loose_or_tight_MU', (df['Muon']['LOOSEMU'] | df['Muon']['TIGHTMU']).sum() == 1) # one muon passing LOOSE or TIGHT criteria
                ## electrons
        selection.add('tight_EL', df['Electron']['TIGHTEL'].sum() == 1) # one electron passing TIGHT criteria
        #selection.add('loose_EL', df['Electron']['LOOSEEL'].sum() == 1) # one electron passing LOOSE criteria
        selection.add('loose_or_tight_EL', (df['Electron']['LOOSEEL'] | df['Electron']['TIGHTEL']).sum() == 1) # one electron passing LOOSE or TIGHT criteria

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
            evt_weights.add('Muon_RECO', mu_reco_cen, mu_reco_err, mu_reco_err, shift=True)
            evt_weights.add('Muon_TRIG', mu_trig_cen, mu_trig_err, mu_trig_err, shift=True)

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
            evt_weights.add('Electron_RECO', el_reco_cen, el_reco_err, el_reco_err, shift=True)
            evt_weights.add('Electron_TRIG', el_trig_cen, el_trig_err, el_trig_err, shift=True)

            # find gen level particles for ttbar system
        #set_trace()
        genp_mode = 'NORMAL'
        GenTTbar = genpsel.select(df, mode=genp_mode)
        #set_trace()
        selection.add('semilep', GenTTbar['SL']['TTbar'].counts > 0)

        ## fill hists for each region
        for lepton in regions.keys():
            for lepcat in regions[lepton].keys():
                for jmult in regions[lepton][lepcat].keys():
                    cut = selection.all(*regions[lepton][lepcat][jmult])

                    #set_trace()
                    if cut.sum() > 0:
                        leptype = 'MU' if lepton == 'Muon' else 'EL'
                        if 'loose_or_tight_%s' % leptype in regions[lepton][lepcat][jmult]:
                            lep_mask = ((df[lepton][cut]['TIGHT%s' % leptype]) | (df[lepton][cut]['LOOSE%s' % leptype]))
                        elif 'tight_%s' % leptype in regions[lepton][lepcat][jmult]:
                            lep_mask = (df[lepton][cut]['TIGHT%s' % leptype])
                        elif 'loose_%s' % leptype in regions[lepton][lepcat][jmult]:
                            lep_mask = (df[lepton][cut]['LOOSE%s' % leptype])
                        else:
                            raise ValueError("Not sure what lepton type to choose for event")

                            ## create MT regions
                        MT = make_vars.MT(df[lepton][cut][lep_mask], df['MET'][cut])
                        MTHigh = (MT >= MTcut).flatten()

                        #print(jmult)

                        test_perms = tperm.find_permutations(jets=df['Jet'][cut], leptons=df[lepton][cut][lep_mask], MET=df['MET'][cut], btagWP=btag_wp)
                        matched_perm = ttmatcher.best_match(gen_hyp=GenTTbar[cut], jets=df['Jet'][cut], leptons=df[lepton][cut][lep_mask], met=df['MET'][cut])

                        ## apply lepton SFs to MC (only applicable to tight leptons)
                        if 'LeptonSF' in corrections.keys():
                            wts = evt_weights.partial_weight(exclude=['Electron_RECO', 'Electron_TRIG']) if lepton == 'Muon' else evt_weights.partial_weight(exclude=['Muon_RECO', 'Muon_TRIG']) # exclude SF from other lepton

                        #set_trace()
                        if jmult == '3Jets':
                            output = self.make_3j_categories(accumulator=output, jetmult=jmult, leptype=lepton, lepcat=lepcat, mtregion='MTHigh', jets=df['Jet'][cut][MTHigh], tps=test_perms[MTHigh], mp=matched_perm[MTHigh], evt_weights=wts[cut][MTHigh])
                            output = self.make_3j_categories(accumulator=output, jetmult=jmult, leptype=lepton, lepcat=lepcat, mtregion='MTLow',  jets=df['Jet'][cut][~MTHigh], tps=test_perms[~MTHigh], mp=matched_perm[~MTHigh], evt_weights=wts[cut][~MTHigh])
                        else:
                            output = self.fill_4pj_hists(accumulator=output, jetmult=jmult, leptype=lepton, lepcat=lepcat, mtregion='MTHigh', jets=df['Jet'][cut][MTHigh], tps=test_perms[MTHigh], mp=matched_perm[MTHigh], evt_weights=wts[cut][MTHigh])
                            output = self.fill_4pj_hists(accumulator=output, jetmult=jmult, leptype=lepton, lepcat=lepcat, mtregion='MTLow', jets=df['Jet'][cut][~MTHigh], tps=test_perms[~MTHigh], mp=matched_perm[~MTHigh], evt_weights=wts[cut][~MTHigh])

        return output

    def fill_4pj_hists(self, accumulator, jetmult, leptype, lepcat, mtregion, jets, tps, mp, evt_weights):
        #set_trace()

        valid_evts = (mp['TTbar'].counts > 0) & ((mp['unique_matches'] == 4).flatten())
        for evt in np.where(valid_evts)[0]:
            #print('mp order:', mp[evt]['BLep'].jetIdx[0], mp[evt]['BHad'].jetIdx[0], mp[evt]['WJa'].jetIdx[0], mp[evt]['WJb'].jetIdx[0])
            for tp in tps[evt]:
                tp_blep_ind, tp_bhad_ind, tp_wja_ind, tp_wjb_ind = tp[0].astype(int) if not isinstance(tp[0], list) else tp[0]
                tp_nu = tp[1]
                if (tp_blep_ind == mp[evt]['BLep'].jetIdx[0]) and (tp_bhad_ind == mp[evt]['BHad'].jetIdx[0]) and ( ((tp_wja_ind == mp[evt]['WJa'].jetIdx[0]) and (tp_wjb_ind == mp[evt]['WJb'].jetIdx[0])) or ((tp_wja_ind == mp[evt]['WJb'].jetIdx[0]) and (tp_wjb_ind == mp[evt]['WJa'].jetIdx[0])) ):
                    tp_cat = 'Correct'
                elif (tp_bhad_ind == mp[evt]['BHad'].jetIdx[0]) and ( ((tp_wja_ind == mp[evt]['WJa'].jetIdx[0]) and (tp_wjb_ind == mp[evt]['WJb'].jetIdx[0])) or ((tp_wja_ind == mp[evt]['WJb'].jetIdx[0]) and (tp_wjb_ind == mp[evt]['WJa'].jetIdx[0])) ):
                    tp_cat = 'Right_THad'
                elif (tp_blep_ind == mp[evt]['BLep'].jetIdx[0]):
                    tp_cat = 'Right_TLep'
                else:
                    tp_cat = 'Wrong'

                #print('mp nu:', np.array(mp[evt]['Nu'].chi2).flatten(), ', tp nu:', tp_nu[3], ', tp order:', tp[0], ', cat:', tp_cat)
                #set_trace()
                mW = (jets[evt][tp_wja_ind].p4 + jets[evt][tp_wjb_ind].p4).mass
                mT = (jets[evt][tp_bhad_ind].p4 + jets[evt][tp_wja_ind].p4 + jets[evt][tp_wjb_ind].p4).mass

                accumulator['nusolver_chi2'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, mtregion=mtregion, permcat=tp_cat, nu_chi2=np.array([tp_nu[3]]), weight=np.array([evt_weights[evt]]))
                accumulator['nusolver_dist'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, mtregion=mtregion, permcat=tp_cat, nu_dist=np.array([np.sqrt(tp_nu[3])]), weight=np.array([evt_weights[evt]]))
                accumulator['mWHad_vs_mTHad'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, mtregion=mtregion, permcat=tp_cat, topmass=np.array([mT]), wmass=np.array([mW]), weight=np.array([evt_weights[evt]]))
            #set_trace()

        return accumulator        


    def make_3j_categories(self, accumulator, jetmult, leptype, lepcat, mtregion, jets, tps, mp, evt_weights):
        #set_trace()
        valid_evts = (mp['TTbar'].counts > 0) & ((mp['unique_matches'] >= 3).flatten())

        correct_merged = valid_evts & (mp['Merged_Event'] & (mp['Merged_BHadWJa'] | mp['Merged_BHadWJb'] | mp['Merged_WJets'])).flatten()
        correct_lost   = valid_evts & (mp['Lost_Event'] & (mp['Lost_WJa'] | mp['Lost_WJb'])).flatten()

            ## Lost Events
        for lost_evt in np.where(correct_lost)[0]:
            for tp in tps[lost_evt]:
                tp_blep_ind, tp_bhad_ind, tp_wja_ind = tp[0].astype(int) if not isinstance(tp[0], list) else tp[0]
                tp_nu = tp[1]
                mp_wjet = 'WJa' if mp[lost_evt]['WJb'].size == 0 else 'WJb' # find which wjet exists in matched perm for comparison
                tp_cat = 'Correct' if (tp_blep_ind == mp[lost_evt]['BLep'].jetIdx[0]) and (tp_bhad_ind == mp[lost_evt]['BHad'].jetIdx[0]) and (tp_wja_ind == mp[lost_evt][mp_wjet].jetIdx[0]) else 'Wrong'

                accumulator['Lost_nusolver_chi2'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, mtregion=mtregion, permcat=tp_cat, nu_chi2=np.array([tp_nu[3]]), weight=np.array([evt_weights[lost_evt]]))
                accumulator['Lost_nusolver_dist'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, mtregion=mtregion, permcat=tp_cat, nu_dist=np.array([np.sqrt(tp_nu[3])]), weight=np.array([evt_weights[lost_evt]]))
                mbpjet = (jets[lost_evt][tp_bhad_ind].p4 + jets[lost_evt][tp_wja_ind].p4).mass
                accumulator['Lost_mbpjet'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, mtregion=mtregion, permcat=tp_cat, mbpjet=np.array([mbpjet]), weight=np.array([evt_weights[lost_evt]]))

            ## Merged Events
        for merged_evt in np.where(correct_merged)[0]:
            for tp in tps[merged_evt]:
                tp_blep_ind, tp_bhad_ind, tp_wja_ind = tp[0].astype(int) if not isinstance(tp[0], list) else tp[0]
                tp_nu = tp[1]
                tp_cat = 'Correct' if (tp_blep_ind == mp[merged_evt]['BLep'].jetIdx[0]) and (tp_bhad_ind == mp[merged_evt]['BHad'].jetIdx[0]) and ((tp_wja_ind == mp[merged_evt]['WJa'].jetIdx[0]) or (tp_wja_ind == mp[merged_evt]['WJb'].jetIdx[0])) else 'Wrong'

                accumulator['Merged_nusolver_chi2'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, mtregion=mtregion, permcat=tp_cat, nu_chi2=np.array([tp_nu[3]]), weight=np.array([evt_weights[merged_evt]]))
                accumulator['Merged_nusolver_dist'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, mtregion=mtregion, permcat=tp_cat, nu_dist=np.array([np.sqrt(tp_nu[3])]), weight=np.array([evt_weights[merged_evt]]))
                maxmjet = np.maximum(jets[merged_evt][tp_bhad_ind].p4.mass, jets[merged_evt][tp_wja_ind].p4.mass)
                mbpjet = (jets[merged_evt][tp_bhad_ind].p4 + jets[merged_evt][tp_wja_ind].p4).mass
                accumulator['Merged_mbpjet_vs_maxmjet'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, lepcat=lepcat, mtregion=mtregion, permcat=tp_cat, maxmjet=np.array([maxmjet]), mbpjet=np.array([mbpjet]), weight=np.array([evt_weights[merged_evt]]))

        return accumulator

    def postprocess(self, accumulator):
        return accumulator

proc_executor = processor.iterative_executor if args.debug else processor.futures_executor
output = processor.run_uproot_job(fileset,
    treename='Events',
    processor_instance=permProbComputer(),
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

save(output, args.outfname)
print('%s has been written' % args.outfname)
