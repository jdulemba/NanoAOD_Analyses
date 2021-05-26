#!/usr/bin/env python

from coffea import hist
from coffea.util import save, load
import coffea.processor as processor
from pdb import set_trace
import os, sys
import python.ObjectSelection as objsel
import coffea.processor.dataframe
import Utilities.plot_tools as plt_tools
import python.MCWeights as MCWeights
import Utilities.make_variables as make_vars
import numpy as np
import Utilities.prettyjson as prettyjson
import coffea.lumi_tools.lumi_tools as lumi_tools
import python.GenParticleSelector as genpsel
import python.TTGenMatcher as ttmatcher

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'ttbar_evt_cats'

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

    ## specify ttJets samples
allowed_samples = ['ttJets_PS'] if args.year == '2016' else ['ttJetsSL']
isAllowed = np.array([(key in allowed_samples) for key in fileset.keys()]).all()
if not isAllowed:
    raise ValueError("Not a valid dataset to run on")

## load corrections for event weights
pu_correction = load(os.path.join(proj_dir, 'Corrections', jobid, 'MC_PU_Weights.coffea'))
lepSF_correction = load(os.path.join(proj_dir, 'Corrections', jobid, 'leptonSFs.coffea'))
jet_corrections = load(os.path.join(proj_dir, 'Corrections', jobid, 'JetCorrections.coffea'))[args.year]
corrections = {
    'Pileup' : pu_correction,
    'Prefire' : True,
    'LeptonSF' : lepSF_correction,
    'BTagSF' : True,
    #'BTagSF' : False,
    'JetCor' : jet_corrections,
}

cfg_pars = prettyjson.loads(open('%s/cfg_files/cfg_pars_%s.json' % (proj_dir, jobid)).read())
    ## parameters for b-tagging
jet_pars = cfg_pars['Jets']
btagger = jet_pars['btagger']

wps_to_use = list(set([jet_pars['permutations']['tightb'],jet_pars['permutations']['looseb']]))
if not( len(wps_to_use) == 1):
    raise IOError("Only 1 unique btag working point supported now")
btag_wp = btagger+wps_to_use[0]

if corrections['BTagSF'] == True:
    sf_file = os.path.join(proj_dir, 'Corrections', jobid, jet_pars['btagging']['btagSF_file'])
    if not os.path.isfile(sf_file):
        raise IOError("BTag SF file %s doesn't exist" % sf_file)

    btag_sfs = load(sf_file)
    corrections.update({'BTag_Constructors' : {}})
    threeJets = btag_sfs[args.year][btagger]['3Jets'][wps_to_use[0]]
    fourPJets = btag_sfs[args.year][btagger]['4PJets'][wps_to_use[0]]
    corrections['BTag_Constructors'].update({btagger : {'3Jets' : threeJets, '4PJets' : fourPJets} })
#set_trace()

MTcut = jet_pars['MT']


# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class ttbar_evt_cats(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")
        self.jetmult_axis = hist.Cat("jmult", "nJets")
        self.mtreg_axis = hist.Cat("mtreg", "MT region")
        self.object_axis = hist.Cat("objtype", "Gen Object")
        self.cat_axis = hist.Cat("cat", "Evt Cat")
        self.pt_axis = hist.Bin("pt", "p_{T} [GeV]", 200, 0, 1000)
        self.eta_axis = hist.Bin("eta", r"$\eta$", 100, -5, 5)
        self.mtt_axis = hist.Bin("mtt", "Gen m($t\overline{t}$) [GeV]", 360, 200, 2000)

            ## make dictionary of hists
        histo_dict = {}
                ## make jet hists
        gen_hists = self.make_hists()
        histo_dict.update(gen_hists)

        histo_dict['cutflow'] = processor.defaultdict_accumulator(int)

        self._accumulator = processor.dict_accumulator(histo_dict)
        self.sample_name = ''
        self.corrections = corrections

    
    @property
    def accumulator(self):
        return self._accumulator


    def make_hists(self):
        histo_dict = {}
        histo_dict['mtt']      = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.mtreg_axis, self.cat_axis, self.mtt_axis)
        histo_dict['pass_mtt'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.mtreg_axis, self.cat_axis, self.object_axis, self.mtt_axis)
        histo_dict['fail_mtt'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.mtreg_axis, self.cat_axis, self.object_axis, self.mtt_axis)

        return histo_dict



    def process(self, df):
        np.random.seed(10) # sets seed so values from random distributions are reproducible (JER corrections)
        output = self.accumulator.identity()

        if not isinstance(df, coffea.processor.dataframe.LazyDataFrame):
            raise IOError("This function only works for LazyDataFrame objects")

        self.sample_name = df.dataset

            ## make event weights
                # data or MC distinction made internally
        mu_evt_weights = MCWeights.get_event_weights(df, year=args.year, corrections=self.corrections, BTagSFs=[btagger], isTTbar=True)
        el_evt_weights = MCWeights.get_event_weights(df, year=args.year, corrections=self.corrections, BTagSFs=[btagger], isTTbar=True)

            ## initialize selections and regions
        selection = processor.PackedSelection()
        regions = {
            'Muon' : {
                'Tight' : {
                    '3Jets'  : {'objselection', 'jets_3' , 'tight_MU', 'btag_pass', 'semilep'},
                    #'4PJets' : {'objselection', 'jets_4p', 'tight_MU', 'btag_pass', 'semilep'},
                },
            },
            'Electron' : {
                'Tight' : {
                    '3Jets'  : {'objselection', 'jets_3' , 'tight_EL', 'btag_pass', 'semilep'},
                    #'4PJets' : {'objselection', 'jets_4p', 'tight_EL', 'btag_pass', 'semilep'},
                },
            },
        }

            ## object selection
        objsel_evts = objsel.select(df, year=args.year, corrections=self.corrections, cutflow=output['cutflow'])
        output['cutflow']['nEvts passing jet and lepton obj selection'] += objsel_evts.sum()
        selection.add('jets_3', df['Jet'].counts == 3)
        selection.add('jets_4p', df['Jet'].counts > 3)
        selection.add('objselection', objsel_evts)
        selection.add('btag_pass', df['Jet'][btag_wp].sum() >= 2)

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
            mu_evt_weights.add('Muon_RECO', mu_reco_cen, mu_reco_err, mu_reco_err, shift=True)
            mu_evt_weights.add('Muon_TRIG', mu_trig_cen, mu_trig_err, mu_trig_err, shift=True)

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
            el_evt_weights.add('Electron_RECO', el_reco_cen, el_reco_err, el_reco_err, shift=True)
            el_evt_weights.add('Electron_TRIG', el_trig_cen, el_trig_err, el_trig_err, shift=True)

                 ## apply btagging SFs to MC
            if corrections['BTagSF'] == True:
                #set_trace()
                deepcsv_cen = np.ones(df.size)
                deepcsv_bc_up = np.ones(df.size)
                deepcsv_bc_dw = np.ones(df.size)
                deepcsv_l_up = np.ones(df.size)
                deepcsv_l_dw = np.ones(df.size)
                threeJets_cut = selection.require(objselection=True, jets_3=True)
                deepcsv_3j_wts = self.corrections['BTag_Constructors'][btagger]['3Jets'].get_scale_factor(jets=df['Jet'][threeJets_cut], passing_cut=btag_wp)
                deepcsv_cen[threeJets_cut] = deepcsv_3j_wts['central'].prod()
                deepcsv_bc_up[threeJets_cut] = deepcsv_3j_wts['bc_up'].prod()
                deepcsv_bc_dw[threeJets_cut] = deepcsv_3j_wts['bc_down'].prod()
                deepcsv_l_up[threeJets_cut] = deepcsv_3j_wts['udsg_up'].prod()
                deepcsv_l_dw[threeJets_cut] = deepcsv_3j_wts['udsg_down'].prod()

                fourplusJets_cut = selection.require(objselection=True, jets_4p=True)
                deepcsv_4pj_wts = self.corrections['BTag_Constructors'][btagger]['4PJets'].get_scale_factor(jets=df['Jet'][fourplusJets_cut], passing_cut=btag_wp)
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

            # get gen ttbar objects
        GenTTbar = genpsel.select(df, mode='NORMAL')
        selection.add('semilep', GenTTbar['SL']['TTbar'].counts > 0)


        ## fill hists for each region
        for lepton in regions.keys():
            evt_weights = mu_evt_weights if lepton == 'Muon' else el_evt_weights
            for leptype in regions[lepton].keys():
                for jmult in regions[lepton][leptype].keys():
                    if args.debug: print(lepton, leptype, 'btagPass', jmult)
                    cut = selection.all(*regions[lepton][leptype][jmult])
                    output['cutflow']['nEvts %s' % ', '.join([lepton, leptype, 'btagPass', jmult])] += cut.sum()
                    if cut.sum() == 0: continue

                    ltype = 'MU' if lepton == 'Muon' else 'EL'
                    if 'loose_or_tight_%s' % ltype in regions[lepton][leptype][jmult]:
                        lep_mask = ((df[lepton][cut]['TIGHT%s' % ltype] == True) | (df[lepton][cut]['LOOSE%s' % ltype] == True))
                    elif 'tight_%s' % ltype in regions[lepton][leptype][jmult]:
                        lep_mask = (df[lepton][cut]['TIGHT%s' % ltype] == True)
                    elif 'loose_%s' % ltype in regions[lepton][leptype][jmult]:
                        lep_mask = (df[lepton][cut]['LOOSE%s' % ltype] == True)
                    else:
                        raise ValueError("Not sure what lepton type to choose for event")

                        ## calculate MT
                    MT = make_vars.MT(df[lepton][cut][lep_mask], df['MET'][cut])
                    MTHigh = (MT >= MTcut).flatten()
                    output['cutflow']['nEvts %s: pass MT cut' % ', '.join([lepton, leptype, 'btagPass', jmult])] += MTHigh.sum()
                    MTLow = ~MTHigh
                    output['cutflow']['nEvts %s: fail MT cut' % ', '.join([lepton, leptype, 'btagPass', jmult])] += MTLow.sum()

                        # get matched perm
                    gen_objs = GenTTbar[cut]
                    mp = ttmatcher.best_match(gen_hyp=gen_objs, jets=df['Jet'][cut], leptons=df[lepton][cut][lep_mask], met=df['MET'][cut])

                    #set_trace()
                    for mt_region in ['MTHigh', 'MTLow']:
                        mt_cut = MTHigh if mt_region == 'MTHigh' else MTLow
                        evt_weights_to_use = (evt_weights.weight()*btag_weights['%s_CEN' % btagger])[cut][mt_cut]

                        valid_gen_objs = gen_objs[mt_cut]
                        valid_mp = mp[mt_cut]

                        ## fill hists of gen mttbar for all events
                        output['mtt'].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, mtreg=mt_region, cat='All', mtt=valid_gen_objs['SL']['TTbar'].p4.mass.flatten(), weight=evt_weights_to_use)
                            # for each gen parton find ones that pass and fail jet pt and eta cuts
                        for gen_obj in ['BHad', 'BLep', 'WJa', 'WJb']:
                            pt_eta_mask = ((valid_gen_objs['SL'][gen_obj].p4.pt >= jet_pars['ptmin']) & (np.abs(valid_gen_objs['SL'][gen_obj].p4.eta) <= jet_pars['etamax'])).flatten()
                            output['pass_mtt'].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, mtreg=mt_region, cat='All', objtype=gen_obj, mtt=valid_gen_objs['SL'][pt_eta_mask]['TTbar'].p4.mass.flatten(), weight=evt_weights_to_use[pt_eta_mask])
                            output['fail_mtt'].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, mtreg=mt_region, cat='All', objtype=gen_obj, mtt=valid_gen_objs['SL'][~pt_eta_mask]['TTbar'].p4.mass.flatten(), weight=evt_weights_to_use[~pt_eta_mask])


                        ## fill hists of gen mttbar for each lost jet/partially merged event category
                        evt_cats = ['Merged_BHadBLep', 'Merged_BHadWJa', 'Merged_BHadWJb', 'Merged_BLepWJa', 'Merged_BLepWJb', 'Merged_WJets', 'Lost_BHad', 'Lost_BLep', 'Lost_WJa', 'Lost_WJb']
                        for evt_cat in evt_cats:
                            evt_cat_mask = valid_mp[evt_cat].flatten()
                            output['mtt'].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, mtreg=mt_region, cat=evt_cat, mtt=valid_gen_objs['SL'][evt_cat_mask]['TTbar'].p4.mass.flatten(), weight=evt_weights_to_use[evt_cat_mask])

                                # for each gen parton find ones that pass and fail jet pt and eta cuts
                            for gen_obj in ['BHad', 'BLep', 'WJa', 'WJb']:
                                pt_eta_mask = ((valid_gen_objs[evt_cat_mask]['SL'][gen_obj].p4.pt >= jet_pars['ptmin']) & (np.abs(valid_gen_objs[evt_cat_mask]['SL'][gen_obj].p4.eta) <= jet_pars['etamax'])).flatten()
                                output['pass_mtt'].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, mtreg=mt_region, cat=evt_cat, objtype=gen_obj, mtt=valid_gen_objs[evt_cat_mask]['SL'][pt_eta_mask]['TTbar'].p4.mass.flatten(), weight=evt_weights_to_use[evt_cat_mask][pt_eta_mask])
                                output['fail_mtt'].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, mtreg=mt_region, cat=evt_cat, objtype=gen_obj, mtt=valid_gen_objs[evt_cat_mask]['SL'][~pt_eta_mask]['TTbar'].p4.mass.flatten(), weight=evt_weights_to_use[evt_cat_mask][~pt_eta_mask])
                            
                        ## fill hists of gen mttbar for 'other' (not merged/lost) events
                        other_cat_mask = ~((valid_mp['Merged_Event'] | valid_mp['Lost_Event']).flatten())
                        output['mtt'].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, mtreg=mt_region, cat='Other', mtt=valid_gen_objs['SL'][other_cat_mask]['TTbar'].p4.mass.flatten(), weight=evt_weights_to_use[other_cat_mask])
                            # for each gen parton find ones that pass and fail jet pt and eta cuts
                        for gen_obj in ['BHad', 'BLep', 'WJa', 'WJb']:
                            pt_eta_mask = ((valid_gen_objs[other_cat_mask]['SL'][gen_obj].p4.pt >= jet_pars['ptmin']) & (np.abs(valid_gen_objs[other_cat_mask]['SL'][gen_obj].p4.eta) <= jet_pars['etamax'])).flatten()
                            output['pass_mtt'].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, mtreg=mt_region, cat='Other', objtype=gen_obj, mtt=valid_gen_objs[other_cat_mask]['SL'][pt_eta_mask]['TTbar'].p4.mass.flatten(), weight=evt_weights_to_use[other_cat_mask][pt_eta_mask])
                            output['fail_mtt'].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, mtreg=mt_region, cat='Other', objtype=gen_obj, mtt=valid_gen_objs[other_cat_mask]['SL'][~pt_eta_mask]['TTbar'].p4.mass.flatten(), weight=evt_weights_to_use[other_cat_mask][~pt_eta_mask])

                

        return output



    def postprocess(self, accumulator):
        return accumulator

proc_executor = processor.iterative_executor if args.debug else processor.futures_executor

output = processor.run_uproot_job(fileset,
    treename='Events',
    processor_instance=ttbar_evt_cats(),
    executor=proc_executor,
    executor_args={
        'workers': 8,
        'flatten' : True,
        'compression': 5,
    },
    chunksize=10000 if args.debug else 100000,
)


if args.debug:
    print(output)

save(output, args.outfname)
print('%s has been written' % args.outfname)
