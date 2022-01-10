#!/usr/bin/env python

import time
tic = time.time()

from coffea import hist, processor
from coffea.analysis_tools import PackedSelection
import awkward as ak
from coffea.nanoevents.methods import vector
ak.behavior.update(vector.behavior)

from pdb import set_trace
from coffea.util import save, load
import os
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
import python.IDJet as IDJet
from copy import deepcopy

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
base_jobid = os.environ["base_jobid"]
analyzer = "htt_btag_sb_regions"

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("fset", type=str, help="Fileset dictionary (in string form) to be used for the processor")
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="Specify which year to run over")
parser.add_argument("outfname", type=str, help="Specify output filename, including directory and file extension")
parser.add_argument("opts", type=str, help="Fileset dictionary (in string form) to be used for the processor")
args = parser.parse_args()

# convert input string of fileset dictionary to actual dictionary
fdict = (args.fset).replace("\'", "\"")
fileset = prettyjson.loads(fdict)

if len(fileset.keys()) > 1:
    raise ValueError("Only one topology run at a time in order to determine which corrections and systematics to run")
samplename = list(fileset.keys())[0]

# get dataset classification, used for corrections/systematics
isTTbar_ = samplename.startswith("ttJets")
isTTSL_ = samplename.startswith("ttJetsSL")
Nominal_ttJets = ["ttJetsSL", "ttJetsHad", "ttJetsDiLep"]
isNominal_ttJets_ = samplename in Nominal_ttJets
if not isNominal_ttJets_:
    raise ValueError("This should only be run with nominal ttJets samples")

# convert input string of options dictionary to actual dictionary
odict = (args.opts).replace("\'", "\"")
opts_dict = prettyjson.loads(odict)

    ## set config options passed through argparse
import ast
to_debug = ast.literal_eval(opts_dict.get("debug", "False"))
apply_hem = ast.literal_eval(opts_dict.get("apply_hem", "True"))

## init tt probs for likelihoods
ttpermutator.year_to_run(year=args.year)

#set_trace()
cfg_pars = prettyjson.loads(open(os.path.join(proj_dir, "cfg_files", f"cfg_pars_{jobid}.json")).read())
## load corrections for event weights
pu_correction = load(os.path.join(proj_dir, "Corrections", base_jobid, cfg_pars["corrections"]["pu"]))[args.year]
lepSF_correction = load(os.path.join(proj_dir, "Corrections", base_jobid, cfg_pars["corrections"]["lepton"]))[args.year]
jet_corrections = load(os.path.join(proj_dir, "Corrections", base_jobid, cfg_pars["corrections"]["jetmet"]["sources"]))[args.year]
alpha_corrections = load(os.path.join(proj_dir, "Corrections", base_jobid, cfg_pars["corrections"]["alpha"]))[args.year]["E"]["All_2D"]
nnlo_reweighting = load(os.path.join(proj_dir, "Corrections", base_jobid, cfg_pars["corrections"]["nnlo"]["filename"]))
ewk_reweighting = load(os.path.join(proj_dir, "Corrections", base_jobid, cfg_pars["corrections"]["ewk"]["file"]))
corrections = {
    "Pileup" : pu_correction,
    "Prefire" : True,
    "LeptonSF" : lepSF_correction,
    "BTagSF" : True,
    "JetCor" : jet_corrections,
    "Alpha" : alpha_corrections,
    "NNLO_Rewt" : {"Var" : cfg_pars["corrections"]["nnlo"]["var"], "Correction" : nnlo_reweighting[cfg_pars["corrections"]["nnlo"]["var"]]},
    "EWK_Rewt" : {"Correction" : ewk_reweighting, "wt" : cfg_pars["corrections"]["ewk"]["wt"]},
}

    ## parameters for b-tagging
jet_pars = cfg_pars["Jets"]
btaggers = ["DeepCSV"]

wps_to_use = list(set([jet_pars["permutations"]["tightb"],jet_pars["permutations"]["looseb"]]))
if not( len(wps_to_use) == 1):
    raise IOError("Only 1 unique btag working point supported now")
btag_wps = ["DeepCSVMedium"]

if corrections["BTagSF"] == True:
    sf_file = os.path.join(proj_dir, "Corrections", base_jobid, jet_pars["btagging"]["btagSF_file"])
    if not os.path.isfile(sf_file):
        raise IOError(f"BTag SF file {sf_file} doesn't exist")

    btag_sfs = load(sf_file)
    corrections.update({"BTag_Constructors" : {}})
    for btagger in btaggers:
        threeJets = btag_sfs[args.year][btagger]["3Jets"][wps_to_use[0]]
        fourPJets = btag_sfs[args.year][btagger]["4PJets"][wps_to_use[0]]
        corrections["BTag_Constructors"].update({btagger : {"3Jets" : threeJets, "4PJets" : fourPJets} })
#set_trace()

MTcut = jet_pars["MT"]

# 0 == "" (no gen matching), 1 == "right", 2 == "matchable", 3 == "unmatchable", 4 == "sl_tau", 5 == "other" (not semilep)
perm_cats = {
    0 : "",
    1 : "right",
    2 : "matchable",
    3 : "unmatchable",
    4 : "sl_tau",
    5 : "other",
}

# get systematics to run
event_systematics_to_run = ["nosys"]

# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class htt_btag_sb_regions(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.sys_axis = hist.Cat("sys", "Systematic")
        self.jetmult_axis = hist.Cat("jmult", "nJets")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")
        self.btag_axis = hist.Cat("btag", "BTag Region")

        self.mtt_axis = hist.Bin("mtt", "m($t\overline{t}$) [GeV]", np.around(np.linspace(200., 2000., 361), decimals=0))
        self.ctstar_abs_axis = hist.Bin("ctstar_abs", "|cos($\\theta^{*}$)|", np.around(np.linspace(0., 1., 21), decimals=2))

            ## make dictionary of hists
        histo_dict = {}
            ## make selection plots
        selection_hists = self.make_selection_hists()
        histo_dict.update(selection_hists)        

        self.sample_name = ""
        self.corrections = corrections
        self.event_systematics_to_run = event_systematics_to_run

            ## make dict of cutflow for each systematic variation
        for sys in self.event_systematics_to_run:
            histo_dict["cutflow_%s" % sys] = processor.defaultdict_accumulator(int)
    
        self._accumulator = processor.dict_accumulator(histo_dict)

            ## define regions depending on data or MC, systematic or nosys
        base_regions = {
            "Muon" : {
                    "3Jets"  : {"lep_and_filter_pass", "passing_jets", "jets_3" , "tight_MU", "DeepCSV_pass"},
                    "4PJets" : {"lep_and_filter_pass", "passing_jets", "jets_4p", "tight_MU", "DeepCSV_pass"},
                #},
            },
            "Electron" : {
                    "3Jets"  : {"lep_and_filter_pass", "passing_jets", "jets_3" , "tight_EL", "DeepCSV_pass"},
                    "4PJets" : {"lep_and_filter_pass", "passing_jets", "jets_4p", "tight_EL", "DeepCSV_pass"},
                #},
            },
        }
        self.regions = {sys:deepcopy(base_regions) for sys in self.event_systematics_to_run}
        if isTTSL_:
            for sys in self.regions.keys():
                for lepton in self.regions[sys].keys():
                    for jmult in self.regions[sys][lepton].keys():
                        self.regions[sys][lepton][jmult].update({"semilep"})


    @property
    def accumulator(self):
        return self._accumulator

    def make_selection_hists(self):
        histo_dict = {}

        histo_dict["mtt_vs_tlep_ctstar_abs"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.mtt_axis, self.ctstar_abs_axis)

        return histo_dict


    def process(self, events):
        output = self.accumulator.identity()
        self.sample_name = events.metadata["dataset"]

            ## make event weights
                # data or MC distinction made internally
        #set_trace()
        mu_evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections, isTTbar=isNominal_ttJets_)
        el_evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections, isTTbar=isNominal_ttJets_)
        pdfWeights = events["LHEPdfWeight"]
        #set_trace()

            ## initialize selections
        selection = {evt_sys: PackedSelection() for evt_sys in self.event_systematics_to_run}

            # get all passing leptons
        lep_and_filter_pass = objsel.select_leptons(events, year=args.year, hem_15_16=apply_hem)
        {selection[sys].add("lep_and_filter_pass", lep_and_filter_pass) for sys in selection.keys()} # add passing leptons requirement to all systematics

            ## build corrected jets and MET
        events["Jet"], events["MET"] = IDJet.process_jets(events, args.year, self.corrections["JetCor"])

        ## add different selections
                ## muons
        tight_mu_sel, loose_mu_sel = ak.sum(events["Muon"]["TIGHTMU"], axis=1) == 1, ak.sum(events["Muon"]["LOOSEMU"], axis=1) == 1
        {selection[sys].add("tight_MU", tight_mu_sel) for sys in selection.keys()} # one muon passing TIGHT criteria
                ## electrons
        tight_el_sel, loose_el_sel = ak.sum(events["Electron"]["TIGHTEL"], axis=1) == 1, ak.sum(events["Electron"]["LOOSEEL"], axis=1) == 1
        {selection[sys].add("tight_EL", tight_el_sel) for sys in selection.keys()} # one electron passing TIGHT criteria

        ### apply lepton SFs to MC (only applicable to tight leptons)
        if "LeptonSF" in corrections.keys():
            #set_trace()
            tight_muons = events["Muon"][tight_mu_sel][(events["Muon"][tight_mu_sel]["TIGHTMU"] == True)]
            muSFs_dict =  MCWeights.get_lepton_sf(year=args.year, lepton="Muons", corrections=self.corrections["LeptonSF"],
                pt=ak.flatten(tight_muons["pt"]), eta=ak.flatten(tight_muons["eta"]))
            mu_reco_cen = np.ones(len(events))
            mu_reco_err = np.zeros(len(events))
            mu_trig_cen = np.ones(len(events))
            mu_trig_err = np.zeros(len(events))
            mu_reco_cen[tight_mu_sel] = muSFs_dict["RECO_CEN"]
            mu_reco_err[tight_mu_sel] = muSFs_dict["RECO_ERR"]
            mu_trig_cen[tight_mu_sel] = muSFs_dict["TRIG_CEN"]
            mu_trig_err[tight_mu_sel] = muSFs_dict["TRIG_ERR"]
            mu_evt_weights.add("Lep_RECO", np.copy(mu_reco_cen), np.copy(mu_reco_cen+mu_reco_err), np.copy(mu_reco_cen-mu_reco_err))
            mu_evt_weights.add("Lep_TRIG", np.copy(mu_trig_cen), np.copy(mu_trig_cen+mu_trig_err), np.copy(mu_trig_cen-mu_trig_err))
    
            tight_electrons = events["Electron"][tight_el_sel][(events["Electron"][tight_el_sel]["TIGHTEL"] == True)]
            elSFs_dict = MCWeights.get_lepton_sf(year=args.year, lepton="Electrons", corrections=self.corrections["LeptonSF"],
                pt=ak.flatten(tight_electrons["pt"]), eta=ak.flatten(tight_electrons["etaSC"]))
            el_reco_cen = np.ones(len(events))
            el_reco_err = np.zeros(len(events))
            el_trig_cen = np.ones(len(events))
            el_trig_err = np.zeros(len(events))
            el_reco_cen[tight_el_sel] = elSFs_dict["RECO_CEN"]
            el_reco_err[tight_el_sel] = elSFs_dict["RECO_ERR"]
            el_trig_cen[tight_el_sel] = elSFs_dict["TRIG_CEN"]
            el_trig_err[tight_el_sel] = elSFs_dict["TRIG_ERR"]
            el_evt_weights.add("Lep_RECO", np.copy(el_reco_cen), np.copy(el_reco_cen+el_reco_err), np.copy(el_reco_cen-el_reco_err))
            el_evt_weights.add("Lep_TRIG", np.copy(el_trig_cen), np.copy(el_trig_cen+el_trig_err), np.copy(el_trig_cen-el_trig_err))
            #set_trace()

            # find gen level particles for ttbar system and other ttbar corrections
        if isTTbar_:
            if isTTSL_:
                genpsel.select(events, mode="NORMAL")
                {selection[sys].add("semilep", ak.num(events["SL"]) > 0) for sys in selection.keys()}
            else:
                {selection[sys].add("semilep", np.zeros(len(events), dtype=bool)) for sys in selection.keys()}

            EWcorr_cen = np.ones(len(events))
            EWcorr_unc = np.ones(len(events))
            if "NNLO_Rewt" in self.corrections.keys():
                    # find gen level particles for ttbar system
                nnlo_wts = MCWeights.get_nnlo_weights(self.corrections["NNLO_Rewt"], events)
                EWcorr_cen *= np.copy(ak.to_numpy(nnlo_wts))

            if "EWK_Rewt" in self.corrections.keys():
                    ## NLO EW weights
                if self.corrections["EWK_Rewt"]["wt"] == "Otto":
                    ewk_wts_dict = MCWeights.get_Otto_ewk_weights(self.corrections["EWK_Rewt"]["Correction"], events)
                    EWcorr_cen *= np.copy(ak.to_numpy(ewk_wts_dict["Rebinned_KFactor_1.0"]))
                    EWcorr_unc = ewk_wts_dict["DeltaQCD"]*ewk_wts_dict["Rebinned_DeltaEW_1.0"]

                    mu_evt_weights.add("EWcorr",
                        np.copy(EWcorr_cen),
                        np.copy(EWcorr_cen*EWcorr_unc+1)
                    )
                    el_evt_weights.add("EWcorr",
                        np.copy(EWcorr_cen),
                        np.copy(EWcorr_cen*EWcorr_unc+1)
                    )

                if self.corrections["EWK_Rewt"]["wt"] == "Afiq":
                    ewk_wts_dict = MCWeights.get_ewk_weights(self.corrections["EWK_Rewt"]["Correction"], events)
                    EWcorr_cen *= np.copy(ewk_wts_dict["yt_1"])

                    mu_evt_weights.add("Yukawa",  # really just varying value of Yt
                        np.copy(EWcorr_cen),
                        np.copy(nnlo_wts*np.copy(ewk_wts_dict["yt_1.5"])),
                        np.copy(nnlo_wts*np.copy(ewk_wts_dict["yt_0.5"])),
                    )
                    el_evt_weights.add("Yukawa",
                        np.copy(EWcorr_cen),
                        np.copy(nnlo_wts*np.copy(ewk_wts_dict["yt_1.5"])),
                        np.copy(nnlo_wts*np.copy(ewk_wts_dict["yt_0.5"])),
                    )


            # run over systematics that require changes to event objects (jets+MET)
        for evt_sys in self.event_systematics_to_run:
            output["cutflow_%s" % evt_sys]["lep_and_filter_pass"] += ak.sum(lep_and_filter_pass)
                # jet selection
            passing_jets = objsel.jets_selection(events, year=args.year, cutflow=output["cutflow_%s" % evt_sys], shift=evt_sys, hem_15_16=apply_hem)
            output["cutflow_%s" % evt_sys]["nEvts passing jet and lepton obj selection"] += ak.sum(passing_jets & lep_and_filter_pass)
            selection[evt_sys].add("passing_jets", passing_jets)
            selection[evt_sys].add("jets_3",  ak.num(events["SelectedJets"]) == 3)
            selection[evt_sys].add("jets_4p",  ak.num(events["SelectedJets"]) > 3) # only for getting btag weights
            selection[evt_sys].add("DeepCSV_pass", ak.sum(events["SelectedJets"][btag_wps[0]], axis=1) >= 2)

                # sort jets by btag value
            events["SelectedJets"] = events["SelectedJets"][ak.argsort(events["SelectedJets"]["btagDeepB"], ascending=False)] if btaggers[0] == "DeepCSV" else events["SelectedJets"][ak.argsort(events["SelectedJets"]["btagDeepFlavB"], ascending=False)]

                # btag sidebands
            deepcsv_sorted = events["SelectedJets"][ak.argsort(events["SelectedJets"]["btagDeepB"], ascending=False)]["btagDeepB"]

                ## apply btagging SFs to MC
            if corrections["BTagSF"] == True:
                deepcsv_cen   = np.ones(len(events))
                deepcsv_bc_up = np.ones(len(events))
                deepcsv_bc_dw = np.ones(len(events))
                deepcsv_l_up = np.ones(len(events))
                deepcsv_l_dw = np.ones(len(events))

                threeJets_cut = selection[evt_sys].require(lep_and_filter_pass=True, passing_jets=True, jets_3=True)
                deepcsv_3j_wts = self.corrections["BTag_Constructors"]["DeepCSV"]["3Jets"].get_scale_factor(jets=events["SelectedJets"][threeJets_cut], passing_cut="DeepCSV"+wps_to_use[0])
                deepcsv_cen[threeJets_cut] = ak.prod(deepcsv_3j_wts["central"], axis=1)
                deepcsv_bc_up[threeJets_cut] = ak.prod(deepcsv_3j_wts["bc_up"], axis=1)
                deepcsv_bc_dw[threeJets_cut] = ak.prod(deepcsv_3j_wts["bc_down"], axis=1)
                deepcsv_l_up[threeJets_cut] = ak.prod(deepcsv_3j_wts["udsg_up"], axis=1)
                deepcsv_l_dw[threeJets_cut] = ak.prod(deepcsv_3j_wts["udsg_down"], axis=1)
    
                fourplusJets_cut = selection[evt_sys].require(lep_and_filter_pass=True, passing_jets=True, jets_4p=True)
                deepcsv_4pj_wts = self.corrections["BTag_Constructors"]["DeepCSV"]["4PJets"].get_scale_factor(jets=events["SelectedJets"][fourplusJets_cut], passing_cut="DeepCSV"+wps_to_use[0])
                deepcsv_cen[fourplusJets_cut] = ak.prod(deepcsv_4pj_wts["central"], axis=1)
                deepcsv_bc_up[fourplusJets_cut] = ak.prod(deepcsv_4pj_wts["bc_up"], axis=1)
                deepcsv_bc_dw[fourplusJets_cut] = ak.prod(deepcsv_4pj_wts["bc_down"], axis=1)
                deepcsv_l_up[fourplusJets_cut] = ak.prod(deepcsv_4pj_wts["udsg_up"], axis=1)
                deepcsv_l_dw[fourplusJets_cut] = ak.prod(deepcsv_4pj_wts["udsg_down"], axis=1)
                    # make dict of btag weights
                btag_weights = {
                    "DeepCSV_CEN" : np.copy(deepcsv_cen),
                    "DeepCSV_bc_UP" : np.copy(deepcsv_bc_up),
                    "DeepCSV_bc_DW" : np.copy(deepcsv_bc_dw),
                    "DeepCSV_l_UP" : np.copy(deepcsv_l_up),
                    "DeepCSV_l_DW" : np.copy(deepcsv_l_dw),
                }

            #set_trace()
            ## fill hists for each region
            for lepton in self.regions[evt_sys].keys():
                evt_weights = mu_evt_weights if lepton == "Muon" else el_evt_weights
                for jmult in self.regions[evt_sys][lepton].keys():
                    cut = selection[evt_sys].all(*self.regions[evt_sys][lepton][jmult])

                    output["cutflow_%s" % evt_sys]["nEvts %s" % ", ".join([lepton, jmult])] += cut.sum()

                    if to_debug: print(lepton, jmult)
                    if cut.sum() > 0:
                        ltype = "MU" if lepton == "Muon" else "EL"
                        if "loose_or_tight_%s" % ltype in self.regions[evt_sys][lepton][jmult]:
                            leptons = events[lepton][cut][((events[lepton][cut]["TIGHT%s" % ltype] == True) | (events[lepton][cut]["LOOSE%s" % ltype] == True))]
                        elif "tight_%s" % ltype in self.regions[evt_sys][lepton][jmult]:
                            leptons = events[lepton][cut][(events[lepton][cut]["TIGHT%s" % ltype] == True)]
                        elif "loose_%s" % ltype in self.regions[evt_sys][lepton][jmult]:
                            leptons = events[lepton][cut][(events[lepton][cut]["LOOSE%s" % ltype] == True)]
                        else:
                            raise ValueError("Not sure what lepton type to choose for event")

                            # get jets and MET
                        jets, met = events["SelectedJets"][cut], events["SelectedMET"][cut]

                            # find best permutations
                        best_perms = ttpermutator.find_best_permutations(jets=jets, leptons=leptons, MET=met, btagWP=btag_wps[0], btag_req=True)
                        valid_perms = ak.num(best_perms["TTbar"].pt) > 0
                        output["cutflow_%s" % evt_sys]["nEvts %s: valid perms" % ", ".join([lepton, jmult])] += ak.sum(valid_perms)

                        bp_status = np.zeros(cut.size, dtype=int) # 0 == "" (no gen matching), 1 == "right", 2 == "matchable", 3 == "unmatchable", 4 == "sl_tau", 5 == "noslep"
                            # get matched permutation (semilep ttbar only)
                        if isTTbar_:
                            semilep_evts = selection[evt_sys].require(semilep=True)
                            bp_status[~semilep_evts] = 5
                            if semilep_evts.sum() > 0:
                                    # find matched permutations
                                mp = ttmatcher.best_match(gen_hyp=events["SL"][cut], jets=jets, leptons=leptons, met=met)
                                #set_trace()
                                perm_cat_array = compare_matched_best_perms(mp, best_perms, njets=jmult)
                                bp_status[cut] = perm_cat_array
                                if ak.any(ak.num(events["SL"]["Lepton"].pdgId) != 1): raise ValueError("Number of leptons is incorrect for classifying tau+jets events")
                                sl_tau_evts = ak.where(np.abs(events["SL"]["Lepton"].pdgId) == 15)[0]
                                bp_status[sl_tau_evts] = 4

                            ## create MT regions
                        MT = make_vars.MT(leptons, met)
                        MTHigh = ak.flatten(MT[valid_perms] >= MTcut)
                        output["cutflow_%s" % evt_sys]["nEvts %s: pass MT cut" % ", ".join([lepton, jmult])] += ak.sum(MTHigh)

                            # fill hists for each systematic
                        #set_trace()
                        if to_debug: print("  evt sys:", evt_sys)
                        if evt_sys == "nosys":
                            nom_wts = (evt_weights.weight()*btag_weights["%s_CEN" % btaggers[0]])[cut][valid_perms][MTHigh]
                                # fill hists for nominal
                            output = self.fill_hists(acc=output, sys="nosys", jetmult=jmult, leptype=lepton, permarray=bp_status[cut][valid_perms][MTHigh],
                                perm=best_perms[valid_perms][MTHigh], evt_wts=nom_wts)

                                # fill hists for each pdf variation
                            for pdf_idx in range(ak.size(pdfWeights[0])):
                                if pdf_idx == 101: pdf_name = "alphaSDown"
                                elif pdf_idx == 102: pdf_name = "alphaSUp"
                                else: pdf_name = f"pdf_{pdf_idx}"

                                if to_debug: print(f"\tsysname: {pdf_name}")
                                pdf_wts = pdfWeights[:, pdf_idx][cut][valid_perms][MTHigh]
                                wts = nom_wts*pdf_wts
                                output = self.fill_hists(acc=output, sys=pdf_name, jetmult=jmult, leptype=lepton, permarray=bp_status[cut][valid_perms][MTHigh],
                                    perm=best_perms[valid_perms][MTHigh], evt_wts=wts)


        return output

    def fill_hists(self, acc, sys, jetmult, leptype, permarray, perm, evt_wts):
            ## apply alpha correction for 3Jets
        if (jetmult == "3Jets") and ("Alpha" in self.corrections):
            alpha_corr = self.corrections["Alpha"](172.5/perm["THad"].mass)
            perm["THad"] = perm["THad"].multiply(alpha_corr) # correct thad
            perm["TTbar"] = ak.flatten(perm["THad"]+perm["TLep"]) # correct ttbar

        thad_ctstar, tlep_ctstar = make_vars.ctstar(perm["THad"], perm["TLep"])
        thad_ctstar, tlep_ctstar = ak.flatten(thad_ctstar, axis=None), ak.flatten(tlep_ctstar, axis=None)

        #set_trace()
        for permval in np.unique(permarray).tolist():
            perm_inds = np.where(permarray == permval)
            dataset_name = "%s_%s" % (self.sample_name, perm_cats[permval]) if permval != 0 else self.sample_name

            acc["mtt_vs_tlep_ctstar_abs"].fill(dataset=dataset_name, sys=sys,  jmult=jetmult, leptype=leptype,
                mtt=ak.flatten(perm["TTbar"].mass)[perm_inds], ctstar_abs=np.abs(tlep_ctstar[perm_inds]), weight=evt_wts[perm_inds])

        return acc        

    def postprocess(self, accumulator):
        return accumulator


proc_executor = processor.iterative_executor if to_debug else processor.futures_executor
proc_exec_args = {"schema": processor.NanoAODSchema} if to_debug else {"schema": processor.NanoAODSchema, "workers": 8}
output = processor.run_uproot_job(
    fileset,
    treename="Events",
    processor_instance=htt_btag_sb_regions(),
    executor=proc_executor,
    executor_args=proc_exec_args,
    chunksize=100000,
)

save(output, args.outfname)
print(f"{args.outfname} has been written")

toc = time.time()
print("Total time: %.1f" % (toc - tic))
