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

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'ttbar_alpha_reco'

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

    ## parameters for b-tagging
jet_pars = prettyjson.loads(open('%s/cfg_files/cfg_pars_%s.json' % (proj_dir, jobid)).read())['Jets']
btagger = jet_pars['btagger']
wps_to_use = list(set([jet_pars['permutations']['tightb'],jet_pars['permutations']['looseb']]))
if not( len(wps_to_use) == 1):
    raise ValueError("Only 1 unique btag working point supported now")
btag_wp = btagger+wps_to_use[0]

if corrections['BTagSF'] == True:
    sf_file = '%s/Corrections/%s/%s' % (proj_dir, jobid, jet_pars['btagging']['btagSF_file'])
    if not os.path.isfile(sf_file):
        raise ValueError("BTag SF file %s doesn't exist" % sf_file)

    btag_sfs = load(sf_file)
    threeJets = btag_sfs[args.year][btagger]['3Jets'][wps_to_use[0]]
    fourPJets = btag_sfs[args.year][btagger]['4PJets'][wps_to_use[0]]
    corrections.update({'BTag_Constructors' : {'3Jets' : threeJets, '4PJets' : fourPJets} })


# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class ttbar_alpha_reco(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")
        self.bp_cat_axis  = hist.Cat("bp_cat", "Best Perm Categoryization")
        self.bp_mtt_axis  = hist.Bin("bp_mtt", "m($t\overline{t}$)", 36, 200., 2000.)
        self.alpha_P_axis   = hist.Bin("alpha_p", "Gen/Reco P($t_{h}$)", 500, 0., 10.)
        self.alpha_E_axis   = hist.Bin("alpha_e", "Gen/Reco E($t_{h}$)", 500, 0., 10.)
        self.norm_mthad_axis= hist.Bin("norm_mthad", "172.5/m($t_{had}$)", 500, 0., 5.)

            ## make dictionary of hists
        histo_dict = {}
            ## make jet hists
                # Gen
        alpha_hists = self.make_alpha_hists()
        histo_dict.update(alpha_hists)

        #histo_dict['cutflow'] = processor.defaultdict_accumulator(int)

        self._accumulator = processor.dict_accumulator(histo_dict)
        self.sample_name = ''
        self.corrections = corrections

    
    @property
    def accumulator(self):
        return self._accumulator

    def make_alpha_hists(self):
        histo_dict = {}
        histo_dict['Alpha_THad_P'] = hist.Hist("Events", self.dataset_axis, self.leptype_axis, self.bp_cat_axis, self.bp_mtt_axis, self.norm_mthad_axis, self.alpha_P_axis)
        histo_dict['Alpha_THad_E'] = hist.Hist("Events", self.dataset_axis, self.leptype_axis, self.bp_cat_axis, self.bp_mtt_axis, self.norm_mthad_axis, self.alpha_E_axis)

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
                'objselection', 'jets_3' , 'tight_MU', 'btag_pass', 'semilep'
            },
            'Electron' : {
                'objselection', 'jets_3' , 'tight_EL', 'btag_pass', 'semilep'
            },
        }

            ## object selection
        objsel_evts = objsel.select(df, year=args.year, corrections=self.corrections)#, accumulator=output)
        #output['cutflow']['nEvts passing jet and lepton obj selection'] += objsel_evts.sum()
        selection.add('jets_3', df['Jet'].counts == 3)
        selection.add('objselection', objsel_evts)
        selection.add('btag_pass', df['Jet'][btag_wp].sum() >= 2)

            ## add different selections
                ## muons
        selection.add('tight_MU', df['Muon']['TIGHTMU'].sum() == 1) # one muon passing TIGHT criteria
        #selection.add('loose_or_tight_MU', (df['Muon']['LOOSEMU'] | df['Muon']['TIGHTMU']).sum() == 1) # one muon passing LOOSE or TIGHT criteria

                ## electrons
        selection.add('tight_EL', df['Electron']['TIGHTEL'].sum() == 1) # one electron passing TIGHT criteria
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

        if corrections['BTagSF'] == True:
            #set_trace()
            threeJets_cut = selection.require(objselection=True, jets_3=True)
            threeJets_btagwts = self.corrections['BTag_Constructors']['3Jets'].get_scale_factor(jets=df['Jet'][threeJets_cut], passing_cut=btag_wp)
            evt_weights._weights['Btag_SF'][threeJets_cut] = threeJets_btagwts['central'].prod()


            # find gen level particles for ttbar system
        #set_trace()
        genp_mode = 'NORMAL'
        GenTTbar = genpsel.select(df, systype='FINAL', mode=genp_mode)
        #genpsel.select(df, mode='LHE')
        #set_trace()
        selection.add('semilep', GenTTbar['SL']['TTbar'].counts > 0)

        ## fill hists for each region
        for lepton in regions.keys():
            cut = selection.all(*regions[lepton])

            #set_trace()
            if cut.sum() > 0:
                leptype = 'MU' if lepton == 'Muon' else 'EL'
                if 'loose_or_tight_%s' % leptype in regions[lepton]:
                    lep_mask = ((df[lepton][cut]['TIGHT%s' % leptype]) | (df[lepton][cut]['LOOSE%s' % leptype]))
                elif 'tight_%s' % leptype in regions[lepton]:
                    lep_mask = (df[lepton][cut]['TIGHT%s' % leptype])
                elif 'loose_%s' % leptype in regions[lepton]:
                    lep_mask = (df[lepton][cut]['LOOSE%s' % leptype])
                else:
                    raise ValueError("Not sure what lepton type to choose for event")

                matched_perm = ttmatcher.best_match(gen_hyp=GenTTbar[cut], jets=df['Jet'][cut], leptons=df[lepton][cut][lep_mask], met=df['MET'][cut])

                    # find best permutations
                best_perms = ttpermutator.find_best_permutations(jets=df['Jet'][cut], leptons=df[lepton][cut][lep_mask], MET=df['MET'][cut], btagWP=btag_wp)
                valid_perms = best_perms['TTbar'].counts > 0

                    ## create MT regions
                MT = make_vars.MT(df[lepton][cut][lep_mask], df['MET'][cut])
                MTHigh = (MT[valid_perms] >= 40.).flatten()

                evt_weights_to_use = evt_weights.weight()
                ## apply lepton SFs to MC (only applicable to tight leptons)
                if 'LeptonSF' in corrections.keys():
                    evt_weights_to_use = evt_weights.partial_weight(exclude=['Electron_SF']) if lepton == 'Muon' else evt_weights.partial_weight(exclude=['Muon_SF']) # exclude SF from other lepton

                #set_trace()
                output = self.make_3j_categories(acc=output, leptype=lepton, genttbar=GenTTbar[cut][valid_perms][MTHigh], mp=matched_perm[valid_perms][MTHigh], bp=best_perms, mtcut=MTHigh, evt_weights=evt_weights_to_use[cut][valid_perms][MTHigh])

        return output


    def make_3j_categories(self, acc, leptype, genttbar, mp, bp, mtcut, evt_weights):
        #set_trace()

        valid_bp = bp['TTbar'].counts > 0
        valid_evts = (mp['TTbar'].counts > 0) & ((mp['unique_matches'] >= 3).flatten())

        correct_mp_merged = valid_evts & (mp['Merged_Event'] & (mp['Merged_BHadWJa'] | mp['Merged_BHadWJb'] | mp['Merged_WJets'])).flatten()
        correct_mp_lost   = valid_evts & (mp['Lost_Event'] & (mp['Lost_WJa'] | mp['Lost_WJb'])).flatten()

            ## initialize arrays for best perm categories
        bp_right = np.zeros(valid_evts.size).astype(bool)
        bp_matchable = np.zeros(valid_evts.size).astype(bool)
        bp_unmatchable = ~(correct_mp_merged | correct_mp_lost)

        ## classify best perm based on matched perm objects
            # merged events
        bp_merged_right_matchable = np.array(((mp[correct_mp_merged]['BLep'].jetIdx == bp['BLep'][valid_bp][mtcut][correct_mp_merged].jetIdx) & (mp[correct_mp_merged]['BHad'].jetIdx == bp['BHad'][valid_bp][mtcut][correct_mp_merged].jetIdx)\
            & ( (mp[correct_mp_merged]['WJa'].jetIdx == bp['WJa'][valid_bp][mtcut][correct_mp_merged].jetIdx) | (mp[correct_mp_merged]['WJb'].jetIdx == bp['WJa'][valid_bp][mtcut][correct_mp_merged].jetIdx) )).flatten()) # right are True, matchable False

        bp_right[correct_mp_merged] = bp_merged_right_matchable
        bp_matchable[correct_mp_merged] = ~bp_merged_right_matchable

            # lost events
        lost_same_wjet = np.zeros(bp['WJa'][valid_bp][mtcut][correct_mp_lost].size).astype(bool)
        lost_wja = np.where(mp[correct_mp_lost]['WJa'].counts == 0)[0]
        lost_wjb = np.where(mp[correct_mp_lost]['WJb'].counts == 0)[0]
        lost_same_wjet[lost_wjb] = np.array((mp[correct_mp_lost]['WJa'][lost_wjb].jetIdx == bp['WJa'][valid_bp][mtcut][correct_mp_lost][lost_wjb].jetIdx).flatten())
        lost_same_wjet[lost_wja] = np.array((mp[correct_mp_lost]['WJb'][lost_wja].jetIdx == bp['WJa'][valid_bp][mtcut][correct_mp_lost][lost_wja].jetIdx).flatten())
        bp_lost_right_matchable = np.array(((mp[correct_mp_lost]['BLep'].jetIdx == bp['BLep'][valid_bp][mtcut][correct_mp_lost].jetIdx) & (mp[correct_mp_lost]['BHad'].jetIdx == bp['BHad'][valid_bp][mtcut][correct_mp_lost].jetIdx)).flatten()) & lost_same_wjet

        bp_right[correct_mp_lost] = bp_lost_right_matchable
        bp_matchable[correct_mp_lost] = ~bp_lost_right_matchable

            ## get variables
        bp_thad_E = bp['THad'][valid_bp][mtcut].p4.energy.flatten()
        bp_thad_P = np.sqrt((bp['THad'][valid_bp][mtcut].p4.x**2) + (bp['THad'][valid_bp][mtcut].p4.y**2) + (bp['THad'][valid_bp][mtcut].p4.z**2)).flatten()
        bp_thad_M = bp['THad'][valid_bp][mtcut].p4.mass.flatten()
        bp_mtt = bp['TTbar'][valid_bp][mtcut].p4.mass.flatten()

        gen_thad_E = genttbar['SL']['THad'].p4.energy.flatten()
        gen_thad_P = np.sqrt((genttbar['SL']['THad'].p4.x**2) + (genttbar['SL']['THad'].p4.y**2) + (genttbar['SL']['THad'].p4.z**2)).flatten()

            # right matching
        acc['Alpha_THad_P'].fill(dataset=self.sample_name, leptype=leptype, bp_cat='Right', bp_mtt=bp_mtt[bp_right], norm_mthad=172.5/bp_thad_M[bp_right], alpha_p=gen_thad_P[bp_right]/bp_thad_P[bp_right], weight=evt_weights[bp_right])
        acc['Alpha_THad_E'].fill(dataset=self.sample_name, leptype=leptype, bp_cat='Right', bp_mtt=bp_mtt[bp_right], norm_mthad=172.5/bp_thad_M[bp_right], alpha_e=gen_thad_E[bp_right]/bp_thad_E[bp_right], weight=evt_weights[bp_right])

            # matchable matching
        acc['Alpha_THad_P'].fill(dataset=self.sample_name, leptype=leptype, bp_cat='Matchable', bp_mtt=bp_mtt[bp_matchable], norm_mthad=172.5/bp_thad_M[bp_matchable], alpha_p=gen_thad_P[bp_matchable]/bp_thad_P[bp_matchable], weight=evt_weights[bp_matchable])
        acc['Alpha_THad_E'].fill(dataset=self.sample_name, leptype=leptype, bp_cat='Matchable', bp_mtt=bp_mtt[bp_matchable], norm_mthad=172.5/bp_thad_M[bp_matchable], alpha_e=gen_thad_E[bp_matchable]/bp_thad_E[bp_matchable], weight=evt_weights[bp_matchable])

            # unmatchable matching
        acc['Alpha_THad_P'].fill(dataset=self.sample_name, leptype=leptype, bp_cat='Unmatchable', bp_mtt=bp_mtt[bp_unmatchable], norm_mthad=172.5/bp_thad_M[bp_unmatchable], alpha_p=gen_thad_P[bp_unmatchable]/bp_thad_P[bp_unmatchable], weight=evt_weights[bp_unmatchable])
        acc['Alpha_THad_E'].fill(dataset=self.sample_name, leptype=leptype, bp_cat='Unmatchable', bp_mtt=bp_mtt[bp_unmatchable], norm_mthad=172.5/bp_thad_M[bp_unmatchable], alpha_e=gen_thad_E[bp_unmatchable]/bp_thad_E[bp_unmatchable], weight=evt_weights[bp_unmatchable])

        return acc

    def postprocess(self, accumulator):
        return accumulator

proc_executor = processor.iterative_executor if args.debug else processor.futures_executor

output = processor.run_uproot_job(fileset,
    treename='Events',
    processor_instance=ttbar_alpha_reco(),
    executor=proc_executor,
    executor_args={
        'workers': 8,
        'flatten' : True,
        'compression': 5,
    },
    chunksize=10000 if args.debug else 50000,
    #chunksize=50000,
)


if args.debug:
    print(output)
#set_trace()
#print(output['cutflow'])

save(output, args.outfname)
print('%s has been written' % args.outfname)
