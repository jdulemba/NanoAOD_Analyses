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
import python.TTPermutator as ttpermutator
import Utilities.make_variables as make_vars
from python.Permutations import compare_matched_best_perms

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
base_jobid = os.environ['base_jobid']
analyzer = 'ttbar_alpha_reco'

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

    ## specify ttJets samples
if (args.year == '2016') and (base_jobid == 'NanoAODv6'):
    Nominal_ttJets = ['ttJets_PS']
#if args.year == '2016':
#    Nominal_ttJets = ['ttJets_PS', 'ttJets']
else:
    Nominal_ttJets = ['ttJetsSL']#, 'ttJetsHad', 'ttJetsDiLep']
if not Nominal_ttJets:
    raise ValueError("This should only be run on nominal ttbar events!")


## init tt probs for likelihoods
ttpermutator.year_to_run(year=args.year)

## load corrections for event weights
pu_correction = load(os.path.join(proj_dir, 'Corrections', jobid, 'MC_PU_Weights.coffea'))
lepSF_correction = load(os.path.join(proj_dir, 'Corrections', jobid, 'leptonSFs.coffea'))
jet_corrections = load(os.path.join(proj_dir, 'Corrections', jobid, 'JetCorrections.coffea'))[args.year]
corrections = {
    'Pileup' : pu_correction,
    'Prefire' : True,
    'LeptonSF' : lepSF_correction,
    'BTagSF' : True,
    'JetCor' : jet_corrections,
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
class ttbar_alpha_reco(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")
        self.bp_mtt_axis  = hist.Bin("bp_mtt", "m($t\overline{t}$)", 36, 200., 2000.)
        #self.alpha_P_axis   = hist.Bin("alpha_p", "Gen/Reco P($t_{h}$)", 500, 0., 10.)
        #self.alpha_E_axis   = hist.Bin("alpha_e", "Gen/Reco E($t_{h}$)", 500, 0., 10.)
        self.alpha_P_axis   = hist.Bin("alpha_p", "Gen/Reco P($t_{h}$)", 1000, 0., 5.)
        self.alpha_E_axis   = hist.Bin("alpha_e", "Gen/Reco E($t_{h}$)", 1000, 0., 5.)
        self.norm_mthad_axis= hist.Bin("norm_mthad", "172.5/m($t_{had}$)", 500, 0., 5.)

            ## make dictionary of hists
        histo_dict = {}
            ## make jet hists
                # Gen
        alpha_hists = self.make_alpha_hists()
        histo_dict.update(alpha_hists)

        self._accumulator = processor.dict_accumulator(histo_dict)
        self.sample_name = ''
        self.corrections = corrections

    
    @property
    def accumulator(self):
        return self._accumulator

    def make_alpha_hists(self):
        histo_dict = {}
        histo_dict['Alpha_THad_P'] = hist.Hist("Events", self.dataset_axis, self.leptype_axis, self.bp_mtt_axis, self.norm_mthad_axis, self.alpha_P_axis)
        histo_dict['Alpha_THad_E'] = hist.Hist("Events", self.dataset_axis, self.leptype_axis, self.bp_mtt_axis, self.norm_mthad_axis, self.alpha_E_axis)
        #histo_dict['Alpha_THad_E'] = hist.Hist("Events", self.dataset_axis, self.leptype_axis, self.bp_mtt_axis, self.norm_mthad_axis)

        return histo_dict

    
    def process(self, df):
        np.random.seed(10) # sets seed so values from random distributions are reproducible (JER corrections)
        output = self.accumulator.identity()

        self.sample_name = df.dataset

            ## make event weights
                # data or MC distinction made internally
        mu_evt_weights = MCWeights.get_event_weights(df, year=args.year, corrections=self.corrections, BTagSFs=btaggers)
        el_evt_weights = MCWeights.get_event_weights(df, year=args.year, corrections=self.corrections, BTagSFs=btaggers)

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
        selection.add('jets_3', df['Jet'].counts == 3)
        selection.add('objselection', objsel_evts)
        selection.add('btag_pass', df['Jet'][btag_wps[0]].sum() >= 2)

            # sort jets by btag value
        df['Jet'] = df['Jet'][df['Jet']['btagDeepB'].argsort(ascending=False)] if btagger == 'DeepCSV' else df['Jet'][df['Jet']['btagDeepFlavB'].argsort(ascending=False)]

            ## add different selections
                ## muons
        selection.add('tight_MU', df['Muon']['TIGHTMU'].sum() == 1) # one muon passing TIGHT criteria

                ## electrons
        selection.add('tight_EL', df['Electron']['TIGHTEL'].sum() == 1) # one electron passing TIGHT criteria

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

        if corrections['BTagSF'] == True:
            #set_trace()
            deepcsv_cen = np.ones(df.size)
            threeJets_cut = selection.require(objselection=True, jets_3=True)
            deepcsv_3j_wts = self.corrections['BTag_Constructors']['DeepCSV']['3Jets'].get_scale_factor(jets=df['Jet'][threeJets_cut], passing_cut='DeepCSV'+wps_to_use[0])
            deepcsv_cen[threeJets_cut] = deepcsv_3j_wts['central'].prod()
        
                # make dict of btag weights
            btag_weights = {
                'DeepCSV_CEN' : deepcsv_cen,
            }


            # find gen level particles for ttbar system
        GenTTbar = genpsel.select(df, mode='NORMAL')
        selection.add('semilep', GenTTbar['SL']['TTbar'].counts > 0)

        ## fill hists for each region
        for lepton in regions.keys():
            evt_weights = mu_evt_weights if lepton == 'Muon' else el_evt_weights
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

                #set_trace()
                    # find matched permutations
                mp = ttmatcher.best_match(gen_hyp=GenTTbar[cut], jets=df['Jet'][cut], leptons=df[lepton][cut][lep_mask], met=df['MET'][cut])

                    # find best permutations
                best_perms = ttpermutator.find_best_permutations(jets=df['Jet'][cut], leptons=df[lepton][cut][lep_mask], MET=df['MET'][cut], btagWP=btag_wps[0])
                valid_perms = best_perms['TTbar'].counts > 0

                    # compare matched per to best perm
                bp_status = np.zeros(cut.size, dtype=int) # 0 == '' (no gen matching), 1 == 'right', 2 == 'matchable', 3 == 'unmatchable', 4 == 'sl_tau', 5 == 'noslep'
                perm_cat_array = compare_matched_best_perms(mp, best_perms, njets='3Jets')
                bp_status[(cut)] = perm_cat_array
                sl_tau_evts = (np.abs(GenTTbar['SL']['Lepton'].pdgId) == 15).pad(1).fillna(False).flatten() # fillna to make sure array is same size as df.size
                bp_status[sl_tau_evts] = 4

                    ## create MT regions
                MT = make_vars.MT(df[lepton][cut][lep_mask], df['MET'][cut])
                MTHigh = (MT[valid_perms] >= MTcut).flatten()

                evt_weights_to_use = (evt_weights.weight()*btag_weights['%s_CEN' % btaggers[0]])[cut][valid_perms][MTHigh]

                #set_trace()
                output = self.make_3j_categories(acc=output, leptype=lepton, permarray=bp_status[cut][valid_perms][MTHigh], genttbar=GenTTbar[cut][valid_perms][MTHigh], bp=best_perms, mtcut=MTHigh, evt_weights=evt_weights_to_use)

        return output


    def make_3j_categories(self, acc, leptype, permarray, genttbar, bp, mtcut, evt_weights):
        #set_trace()
            ## get variables
        bp_thad_E = bp['THad'].p4.energy.flatten()[mtcut]
        bp_thad_P = np.sqrt((bp['THad'].p4.x**2) + (bp['THad'].p4.y**2) + (bp['THad'].p4.z**2)).flatten()[mtcut]
        bp_thad_M = bp['THad'].p4.mass.flatten()[mtcut]
        bp_mtt = bp['TTbar'].p4.mass.flatten()[mtcut]

        gen_thad_E = genttbar['SL']['THad'].p4.energy.flatten()
        gen_thad_P = np.sqrt((genttbar['SL']['THad'].p4.x**2) + (genttbar['SL']['THad'].p4.y**2) + (genttbar['SL']['THad'].p4.z**2)).flatten()

        for permval in np.unique(permarray).tolist():
            #set_trace()
            if permval != 1: continue # only fill hists for 'right' events
            perm_inds = np.where(permarray == permval)
            dataset_name = '%s_%s' % (self.sample_name, perm_cats[permval]) if permval != 0 else self.sample_name

            acc['Alpha_THad_P'].fill(dataset=dataset_name, leptype=leptype, bp_mtt=bp_mtt[perm_inds], norm_mthad=172.5/bp_thad_M[perm_inds], alpha_p=gen_thad_P[perm_inds]/bp_thad_P[perm_inds], weight=evt_weights[perm_inds])
            acc['Alpha_THad_E'].fill(dataset=dataset_name, leptype=leptype, bp_mtt=bp_mtt[perm_inds], norm_mthad=172.5/bp_thad_M[perm_inds], alpha_e=gen_thad_E[perm_inds]/bp_thad_E[perm_inds], weight=evt_weights[perm_inds])


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
