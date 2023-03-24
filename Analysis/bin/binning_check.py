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
analyzer = "binning_check"

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("fset", type=str, help="Fileset dictionary (in string form) to be used for the processor")
parser.add_argument("year", choices=["2016", "2017", "2018"] if base_jobid == "NanoAODv6" else ["2016APV", "2016", "2017", "2018"], help="Specify which year to run over")
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
isSignal_ = (samplename.startswith("AtoTT") or samplename.startswith("HtoTT"))
isInt_ = isSignal_ and ("Int" in samplename)

Nominal_ttJets = ["ttJetsSL", "ttJetsHad", "ttJetsDiLep"]
isNominal_ttJets_ = samplename in Nominal_ttJets
if (not isTTSL_) and (not isSignal_):
    raise ValueError(f"{samplename} not allowed for this analyzer. Only nominal semileptonic ttJets samples or signal samples are allowed.")

# convert input string of options dictionary to actual dictionary
odict = (args.opts).replace("\'", "\"")
opts_dict = prettyjson.loads(odict)

    ## set config options passed through argparse
import ast
to_debug = ast.literal_eval(opts_dict.get("debug", "False"))
apply_hem = ast.literal_eval(opts_dict.get("apply_hem", "True"))

## init tt probs for likelihoods
ttpermutator.year_to_run(year=args.year)

        # copy fileset root files to local condor node if running on condor
if "isCondor" in opts_dict.keys():
    if ast.literal_eval(opts_dict["isCondor"]):
        from subprocess import check_output, STDOUT
        sites_to_try = ["root://xrootd-cms.infn.it/", "root://cmsxrootd.fnal.gov/"]
        n_retries = len(sites_to_try) + 1
        for idx, rfile in enumerate(fileset[samplename]):
            cp_success = False
            for cp_attempt in range(n_retries):
                if cp_success: continue
                cp_rfile = rfile if cp_attempt == 0 else "/".join([sites_to_try[cp_attempt-1], rfile.split("//")[-1]]) # replace whatever redirector is used to regional Bari one
                print(f"Attempt {cp_attempt+1} to copy {cp_rfile} to /tmp")
                try:
                    output = check_output(["xrdcp", "-f", f"{cp_rfile}", "/tmp"], timeout=1200, stderr=STDOUT)
                    #output = check_output(["xrdcp", "-f", f"{cp_rfile}", "/tmp"], timeout=600, stderr=STDOUT)
                    #output = check_output(["xrdcp", "-f", f"{cp_rfile}", "/tmp"], timeout=300, stderr=STDOUT)
                    cp_success = True
                except:
                    cp_success = False
                    continue
            if not cp_success:
                raise ValueError(f"{cp_rfile} not found")
            fileset[samplename][idx] = f"/tmp/{rfile.split('/')[-1]}"

#set_trace()
cfg_pars = prettyjson.loads(open(os.path.join(proj_dir, "cfg_files", f"cfg_pars_{jobid}.json")).read())
## load corrections for event weights
pu_correction = load(os.path.join(proj_dir, "Corrections", base_jobid, cfg_pars["corrections"]["pu"]))[args.year]
lepSF_correction = load(os.path.join(proj_dir, "Corrections", base_jobid, cfg_pars["corrections"]["lepton"]))[args.year]
jet_corrections = load(os.path.join(proj_dir, "Corrections", base_jobid, cfg_pars["corrections"]["jetmet"][cfg_pars["corrections"]["jetmet"]["to_use"]]))[args.year]
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

jet_pars = cfg_pars["Jets"]
btaggers = [jet_pars["btagger"]]
#btaggers = ["DeepCSV"]

wps_to_use = list(set([jet_pars["permutations"]["tightb"],jet_pars["permutations"]["looseb"]]))
if not( len(wps_to_use) == 1):
    raise IOError("Only 1 unique btag working point supported now")
btag_wp = btaggers[0]+wps_to_use[0]
#btag_wps = ["DeepCSVMedium"]

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
evt_sys_to_run = opts_dict.get("evt_sys", "NONE").upper()
rewt_sys_to_run = opts_dict.get("rewt_sys", "NONE").upper()
only_sys = ast.literal_eval(opts_dict.get("only_sys", "False"))
#set_trace()
event_systematics_to_run = ["nosys"]
reweight_systematics_to_run = ["nosys"]
#set_trace()


# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class binning_check(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.jetmult_axis = hist.Cat("jmult", "nJets")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")

        self.mtt_axis        = hist.Bin("mtt", "m($t\overline{t}$) [GeV]", np.around(np.linspace(200., 2000., 361), decimals=0))
        self.ctstar_axis     = hist.Bin("ctstar", "cos($\\theta^{*}$)", np.around(np.linspace(-1., 1., 41), decimals=2))
        self.ctstar_abs_axis = hist.Bin("ctstar_abs", "|cos($\\theta^{*}$)|", np.around(np.linspace(0., 1., 21), decimals=2))
        self.top_rap_axis    = hist.Bin("top_rap", "yt", np.around(np.linspace(-5., 5., 201), decimals=2))
        self.deltaYtt_axis   = hist.Bin("deltaYtt", "$\\delta Y_{t \overline{t}}$", np.around(np.linspace(-10., 10., 401), decimals=2))
        self.deltaYtt_abs_axis    = hist.Bin("deltaYtt_abs", "|$\\delta Y_{t \overline{t}}$|", np.around(np.linspace(0., 5., 101), decimals=2))
        self.reso_mtt_axis        = hist.Bin("reso_mtt", "m($t\overline{t}$) [GeV]", np.around(np.linspace(-500., 500., 1001), decimals=0))
        self.reso_ctstar_axis     = hist.Bin("reso_ctstar", "cos($\\theta^{*}$)", np.around(np.linspace(-2., 2., 401), decimals=2))
        self.reso_ctstar_abs_axis = hist.Bin("reso_ctstar_abs", "|cos($\\theta^{*}$)|", np.around(np.linspace(-1., 1., 201), decimals=2))
        self.rel_reso_mtt_axis    = hist.Bin("rel_reso_mtt", "m($t\overline{t}$) [GeV]", np.around(np.linspace(-5., 5., 201), decimals=2))

            ## make dictionary of hists
        histo_dict = {}
                ## make reco plots
        reco_hists = self.reco_hists()
        histo_dict.update(reco_hists)
        if isTTSL_:
                ## make gen plots
            gen_hists = self.gen_hists()
            histo_dict.update(gen_hists)
                ## make reso plots
            reso_hists = self.reso_hists()
            histo_dict.update(reso_hists)

        self.sample_name = ""
        self.corrections = corrections
        self.reweight_systematics_to_run = reweight_systematics_to_run
        self.event_systematics_to_run = event_systematics_to_run

            ## make dict of cutflow for each systematic variation
        for sys in self.event_systematics_to_run:
            histo_dict[f"cutflow_{sys}"] = processor.defaultdict_accumulator(int)
    
        self._accumulator = processor.dict_accumulator(histo_dict)

            ## define regions depending on data or MC, systematic or nosys
        base_regions = {
            "Muon" : {
                "btagPass" : {
                    "3Jets"  : {"lep_and_filter_pass", "passing_jets", "jets_3" , "tight_MU", f"{btaggers[0]}_pass"},
                    "4PJets" : {"lep_and_filter_pass", "passing_jets", "jets_4p", "tight_MU", f"{btaggers[0]}_pass"},
                    #"3Jets"  : {"lep_and_filter_pass", "passing_jets", "jets_3" , "tight_MU", "DeepCSV_pass"},
                    #"4PJets" : {"lep_and_filter_pass", "passing_jets", "jets_4p", "tight_MU", "DeepCSV_pass"},
                },
            },
            "Electron" : {
                "btagPass" : {
                    "3Jets"  : {"lep_and_filter_pass", "passing_jets", "jets_3" , "tight_EL", f"{btaggers[0]}_pass"},
                    "4PJets" : {"lep_and_filter_pass", "passing_jets", "jets_4p", "tight_EL", f"{btaggers[0]}_pass"},
                    #"3Jets"  : {"lep_and_filter_pass", "passing_jets", "jets_3" , "tight_EL", "DeepCSV_pass"},
                    #"4PJets" : {"lep_and_filter_pass", "passing_jets", "jets_4p", "tight_EL", "DeepCSV_pass"},
                },
            },
        }
        self.regions = {sys:deepcopy(base_regions) for sys in self.event_systematics_to_run}

        if isTTSL_:
            for sys in self.regions.keys():
                for lepton in self.regions[sys].keys():
                    for btagregion in self.regions[sys][lepton].keys():
                        for jmult in self.regions[sys][lepton][btagregion].keys():
                            self.regions[sys][lepton][btagregion][jmult].update({"semilep"})


    @property
    def accumulator(self):
        return self._accumulator


    def reco_hists(self):
        histo_dict = {}
        histo_dict["RECO_mtt_vs_tlep_ctstar"]     = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.mtt_axis, self.ctstar_axis)
        histo_dict["RECO_mtt_vs_tlep_ctstar_abs"] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.mtt_axis, self.ctstar_abs_axis)
        histo_dict["RECO_mtt_vs_topRapidity"]     = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.mtt_axis, self.top_rap_axis)
        histo_dict["RECO_mtt_vs_deltaYtt"]        = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.mtt_axis, self.deltaYtt_axis)
        histo_dict["RECO_mtt_vs_deltaYtt_abs"]    = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.mtt_axis, self.deltaYtt_abs_axis)

        return histo_dict

    def gen_hists(self):
        histo_dict = {}
        histo_dict["GEN_mtt_vs_tlep_ctstar"]     = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.mtt_axis, self.ctstar_axis)
        histo_dict["GEN_mtt_vs_tlep_ctstar_abs"] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.mtt_axis, self.ctstar_abs_axis)
        histo_dict["GEN_mtt_vs_topRapidity"]     = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.mtt_axis, self.top_rap_axis)
        histo_dict["GEN_mtt_vs_deltaYtt"]        = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.mtt_axis, self.deltaYtt_axis)
        histo_dict["GEN_mtt_vs_deltaYtt_abs"]    = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.mtt_axis, self.deltaYtt_abs_axis)

        return histo_dict

    def reso_hists(self):
        histo_dict = {}
        histo_dict["RESO_mtt_vs_tlep_ctstar"]     = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.reso_mtt_axis, self.reso_ctstar_axis)
        histo_dict["RESO_mtt_vs_tlep_ctstar_abs"] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.reso_mtt_axis, self.reso_ctstar_abs_axis)
        histo_dict["RESO_mtt_vs_mtt"]       = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.mtt_axis, self.reso_mtt_axis)
        histo_dict["RESO_ctstar_vs_mtt"]    = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.mtt_axis, self.reso_ctstar_axis)
        histo_dict["RESO_ctstar_abs_vs_mtt"]= hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.mtt_axis, self.reso_ctstar_abs_axis)
        histo_dict["REL_RESO_mtt_vs_mtt"]   = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.mtt_axis, self.rel_reso_mtt_axis)

        return histo_dict


    def process(self, events):
        output = self.accumulator.identity()

        self.sample_name = events.metadata["dataset"]

            ## make event weights
                # data or MC distinction made internally
        #set_trace()
        mu_evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections, isTTbar=isNominal_ttJets_, isSignal=isSignal_)
        el_evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections, isTTbar=isNominal_ttJets_, isSignal=isSignal_)
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
            tight_muons = events["Muon"][tight_mu_sel][(events["Muon"][tight_mu_sel]["TIGHTMU"] == True)]
            muSFs_dict =  MCWeights.get_lepton_sf(sf_dict=self.corrections["LeptonSF"]["Muons"],
                pt=ak.flatten(tight_muons["pt"]), eta=ak.flatten(tight_muons["eta"]), tight_lep_mask=tight_mu_sel, leptype="Muons")

            tight_electrons = events["Electron"][tight_el_sel][(events["Electron"][tight_el_sel]["TIGHTEL"] == True)]
            elSFs_dict = MCWeights.get_lepton_sf(sf_dict=self.corrections["LeptonSF"]["Electrons"],
                pt=ak.flatten(tight_electrons["pt"]), eta=ak.flatten(tight_electrons["etaSC"]), tight_lep_mask=tight_el_sel, leptype="Electrons")

            # find gen level particles for ttbar system and other ttbar corrections
        if isTTbar_:
            if "NNLO_Rewt" in self.corrections.keys():
                    # find gen level particles for ttbar system
                nnlo_wts = MCWeights.get_nnlo_weights(self.corrections["NNLO_Rewt"], events)
                mu_evt_weights.add("NNLOqcd",
                    np.copy(nnlo_wts),
                )
                el_evt_weights.add("NNLOqcd",
                    np.copy(nnlo_wts),
                )

            if "EWK_Rewt" in self.corrections.keys():
                    ## NLO EW weights
                if self.corrections["EWK_Rewt"]["wt"] == "Otto":
                    ewk_wts_dict = MCWeights.get_Otto_ewk_weights(self.corrections["EWK_Rewt"]["Correction"], events)
                        # add Yukawa coupling variation
                    mu_evt_weights.add("Yukawa",  # really just varying value of Yt
                        np.copy(ewk_wts_dict["Rebinned_KFactor_1.0"]),
                        #np.copy(ewk_wts_dict["Rebinned_KFactor_1.1"]),
                        #np.copy(ewk_wts_dict["Rebinned_KFactor_0.9"]),
                    )
                    el_evt_weights.add("Yukawa",  # really just varying value of Yt
                        np.copy(ewk_wts_dict["Rebinned_KFactor_1.0"]),
                        #np.copy(ewk_wts_dict["Rebinned_KFactor_1.1"]),
                        #np.copy(ewk_wts_dict["Rebinned_KFactor_0.9"]),
                    )

            if isTTSL_:
                genpsel.select(events, mode="NORMAL")
                {selection[sys].add("semilep", ak.num(events["SL"]) > 0) for sys in selection.keys()}
            else:
                {selection[sys].add("semilep", np.zeros(len(events), dtype=bool)) for sys in selection.keys()}


        #if to_debug: set_trace()
            # run over systematics that require changes to event objects (jets+MET)
        for evt_sys in self.event_systematics_to_run:
            output[f"cutflow_{evt_sys}"]["lep_and_filter_pass"] += ak.sum(lep_and_filter_pass)
                # jet selection
            passing_jets = objsel.jets_selection(events, year=args.year, cutflow=output[f"cutflow_{evt_sys}"], shift=evt_sys, hem_15_16=apply_hem)
            output[f"cutflow_{evt_sys}"]["nEvts passing jet and lepton obj selection"] += ak.sum(passing_jets & lep_and_filter_pass)
            selection[evt_sys].add("passing_jets", passing_jets)
            selection[evt_sys].add("jets_3",  ak.num(events["SelectedJets"]) == 3)
            selection[evt_sys].add("jets_4p",  ak.num(events["SelectedJets"]) > 3) # only for getting btag weights
            selection[evt_sys].add(f"{btaggers[0]}_pass", ak.sum(events["SelectedJets"][btag_wp], axis=1) >= 2)
            #selection[evt_sys].add("DeepCSV_pass", ak.sum(events["SelectedJets"][btag_wps[0]], axis=1) >= 2)

                # sort jets by btag value
            events["SelectedJets"] = events["SelectedJets"][ak.argsort(events["SelectedJets"][IDJet.btag_tagger_to_disc_name[btaggers[0]]], ascending=False)]
            #events["SelectedJets"] = events["SelectedJets"][ak.argsort(events["SelectedJets"]["btagDeepB"], ascending=False)] if btaggers[0] == "DeepCSV" else events["SelectedJets"][ak.argsort(events["SelectedJets"]["btagDeepFlavB"], ascending=False)]

                # btag sidebands
            #btagger_sorted = events["SelectedJets"][ak.argsort(events["SelectedJets"][IDJet.btag_tagger_to_disc_name[btaggers[0]]], ascending=False)][IDJet.btag_tagger_to_disc_name[btaggers[0]]]
            #deepcsv_sorted = events["SelectedJets"][ak.argsort(events["SelectedJets"]["btagDeepB"], ascending=False)]["btagDeepB"]

                ## apply btagging SFs to MC
            if (corrections["BTagSF"] == True):
                btag_weights = {key : np.ones(len(events)) for key in self.corrections["BTag_Constructors"][btaggers[0]]["3Jets"].schema_.keys()}

                threeJets_cut = selection[evt_sys].require(lep_and_filter_pass=True, passing_jets=True, jets_3=True)
                btagger_3j_wts = self.corrections["BTag_Constructors"][btaggers[0]]["3Jets"].get_scale_factor(jets=events["SelectedJets"][threeJets_cut], passing_cut=btag_wp)
                fourplusJets_cut = selection[evt_sys].require(lep_and_filter_pass=True, passing_jets=True, jets_4p=True)
                btagger_4pj_wts = self.corrections["BTag_Constructors"][btaggers[0]]["4PJets"].get_scale_factor(jets=events["SelectedJets"][fourplusJets_cut], passing_cut=btag_wp)

                    # fll dict of btag weights
                for wt_name in btagger_3j_wts.keys():
                    btag_weights[wt_name][threeJets_cut] = ak.prod(btagger_3j_wts[wt_name], axis=1)
                    btag_weights[wt_name][fourplusJets_cut] = ak.prod(btagger_4pj_wts[wt_name], axis=1)

                #btag_weights = {key : np.ones(len(events)) for key in self.corrections["BTag_Constructors"]["DeepCSV"]["3Jets"].schema_.keys()}

                #threeJets_cut = selection[evt_sys].require(lep_and_filter_pass=True, passing_jets=True, jets_3=True)
                #deepcsv_3j_wts = self.corrections["BTag_Constructors"]["DeepCSV"]["3Jets"].get_scale_factor(jets=events["SelectedJets"][threeJets_cut], passing_cut="DeepCSV"+wps_to_use[0])
                #fourplusJets_cut = selection[evt_sys].require(lep_and_filter_pass=True, passing_jets=True, jets_4p=True)
                #deepcsv_4pj_wts = self.corrections["BTag_Constructors"]["DeepCSV"]["4PJets"].get_scale_factor(jets=events["SelectedJets"][fourplusJets_cut], passing_cut="DeepCSV"+wps_to_use[0])

                #    # fll dict of btag weights
                #for wt_name in deepcsv_3j_wts.keys():
                #    btag_weights[wt_name][threeJets_cut] = ak.prod(deepcsv_3j_wts[wt_name], axis=1)
                #    btag_weights[wt_name][fourplusJets_cut] = ak.prod(deepcsv_4pj_wts[wt_name], axis=1)

            elif (corrections["BTagSF"] == False):
                raise ValueError("BTag SFs not applied to MC")

            #set_trace()
            ## fill hists for each region
            for lepton in self.regions[evt_sys].keys():
                evt_weights = mu_evt_weights if lepton == "Muon" else el_evt_weights
                lep_SFs = muSFs_dict if lepton == "Muon" else elSFs_dict
                for btagregion in self.regions[evt_sys][lepton].keys():
                    for jmult in self.regions[evt_sys][lepton][btagregion].keys():
                        #set_trace()
                        cut = selection[evt_sys].all(*self.regions[evt_sys][lepton][btagregion][jmult])

                        output[f"cutflow_{evt_sys}"]["nEvts %s" % ", ".join([lepton, btagregion, jmult])] += cut.sum()

                        if to_debug: print(lepton, btagregion, jmult)
                        if cut.sum() > 0:
                            ltype = "MU" if lepton == "Muon" else "EL"
                            if f"loose_or_tight_{ltype}" in self.regions[evt_sys][lepton][btagregion][jmult]:
                                leptons = events[lepton][cut][((events[lepton][cut][f"TIGHT{ltype}"] == True) | (events[lepton][cut][f"LOOSE{ltype}"] == True))]
                            elif f"tight_{ltype}" in self.regions[evt_sys][lepton][btagregion][jmult]:
                                leptons = events[lepton][cut][(events[lepton][cut][f"TIGHT{ltype}"] == True)]
                            elif f"loose_{ltype}" in self.regions[evt_sys][lepton][btagregion][jmult]:
                                leptons = events[lepton][cut][(events[lepton][cut][f"LOOSE{ltype}"] == True)]
                            else:
                                raise ValueError("Not sure what lepton type to choose for event")

                                # get jets and MET
                            jets, met = events["SelectedJets"][cut], events["SelectedMET"][cut]

                                # find best permutations
                            best_perms = ttpermutator.find_best_permutations(jets=jets, leptons=leptons, MET=met, btagWP=btag_wp, btag_req=True if btagregion == "btagPass" else False)
                            #best_perms = ttpermutator.find_best_permutations(jets=jets, leptons=leptons, MET=met, btagWP=btag_wps[0], btag_req=True if btagregion == "btagPass" else False)
                            valid_perms = ak.num(best_perms["TTbar"].pt) > 0
                            output[f"cutflow_{evt_sys}"]["nEvts %s: valid perms" % ", ".join([lepton, btagregion, jmult])] += ak.sum(valid_perms)

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
                            output[f"cutflow_{evt_sys}"]["nEvts %s: pass MT cut" % ", ".join([lepton, btagregion, jmult])] += ak.sum(MTHigh)

                                # fill hists for each systematic
                            if to_debug: print("  evt sys:", evt_sys)
                            if evt_sys == "nosys":
                                for rewt_sys in self.reweight_systematics_to_run:
                                    #if to_debug: set_trace()
                                    if to_debug: print("    sysname:", rewt_sys)

                                    wts = (evt_weights.weight() * btag_weights["central"] * lep_SFs["central"])[cut][valid_perms][MTHigh]
                                    #wts = (evt_weights.weight()*btag_weights["%s_CEN" % btaggers[0]])[cut][valid_perms][MTHigh]

                                        # fill hists for interference samples
                                    if isInt_:
                                        #set_trace()
                                            # fill hists for positive weights
                                        pos_evts = np.where(wts > 0)
                                        self.sample_name = "%s_pos" % events.metadata["dataset"]
                                        output = self.fill_signal_hists(acc=output, jetmult=jmult, leptype=lepton, permarray=bp_status[cut][valid_perms][MTHigh][pos_evts],
                                            perm=best_perms[valid_perms][MTHigh][pos_evts], evt_wts=wts[pos_evts])
                                            # fill hists for negative weights
                                        neg_evts = np.where(wts < 0)
                                        self.sample_name = "%s_neg" % events.metadata["dataset"]
                                        output = self.fill_signal_hists(acc=output, jetmult=jmult, leptype=lepton, permarray=bp_status[cut][valid_perms][MTHigh][neg_evts],
                                            perm=best_perms[valid_perms][MTHigh][neg_evts], evt_wts=wts[neg_evts])
                                    else:
                                        if isTTSL_:
                                            output = self.fill_SL_hists(acc=output, jetmult=jmult, leptype=lepton, permarray=bp_status[cut][valid_perms][MTHigh],
                                                perm=best_perms[valid_perms][MTHigh], gentt=events["SL"][cut][valid_perms][MTHigh], evt_wts=wts)
                                        else:
                                            output = self.fill_signal_hists(acc=output, jetmult=jmult, leptype=lepton, permarray=bp_status[cut][valid_perms][MTHigh],
                                                perm=best_perms[valid_perms][MTHigh], evt_wts=wts)

        return output

    def fill_SL_hists(self, acc, jetmult, leptype, permarray, perm, gentt, evt_wts):
            ## apply alpha correction for 3Jets
        if (jetmult == "3Jets") and ("Alpha" in self.corrections):
            alpha_corr = self.corrections["Alpha"](172.5/perm["THad"].mass)
            perm["THad"] = perm["THad"].multiply(alpha_corr) # correct thad
            perm["TTbar"] = ak.flatten(perm["THad"]+perm["TLep"]) # correct ttbar

        reco_thad_ctstar, reco_tlep_ctstar = make_vars.ctstar(perm["THad"], perm["TLep"], flatten=True)
        gen_thad_ctstar, gen_tlep_ctstar = make_vars.ctstar(gentt["THad"], gentt["TLep"], flatten=True)

        reco_top_rap, reco_tbar_rap, reco_deltaYtt = make_vars.deltaYtt(perm["Top"], perm["Tbar"], return_indiv=True, flatten=True)
        gen_top_rap, gen_tbar_rap, gen_deltaYtt = make_vars.deltaYtt(gentt["Top"], gentt["Tbar"], return_indiv=True, flatten=True)

        #set_trace()
        for permval in np.unique(permarray).tolist():
            perm_inds = np.where(permarray == permval)
            dataset_name = "%s_%s" % (self.sample_name, perm_cats[permval]) if permval != 0 else self.sample_name

            #set_trace()
                # reco plots
            acc["RECO_mtt_vs_tlep_ctstar"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype,
                mtt=ak.flatten(perm["TTbar"].mass, axis=None)[perm_inds], ctstar=reco_tlep_ctstar[perm_inds], weight=evt_wts[perm_inds])
            acc["RECO_mtt_vs_tlep_ctstar_abs"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype,
                mtt=ak.flatten(perm["TTbar"].mass, axis=None)[perm_inds], ctstar_abs=np.abs(reco_tlep_ctstar[perm_inds]), weight=evt_wts[perm_inds])
            acc["RECO_mtt_vs_topRapidity"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype,
                mtt=ak.flatten(perm["TTbar"].mass, axis=None)[perm_inds], top_rap=reco_top_rap[perm_inds], weight=evt_wts[perm_inds])
            acc["RECO_mtt_vs_deltaYtt"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype,
                mtt=ak.flatten(perm["TTbar"].mass, axis=None)[perm_inds], deltaYtt=reco_deltaYtt[perm_inds], weight=evt_wts[perm_inds])
            acc["RECO_mtt_vs_deltaYtt_abs"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype,
                mtt=ak.flatten(perm["TTbar"].mass, axis=None)[perm_inds], deltaYtt_abs=np.abs(reco_deltaYtt[perm_inds]), weight=evt_wts[perm_inds])

                # gen plots
            acc["GEN_mtt_vs_tlep_ctstar"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype,
                mtt=ak.flatten(gentt["TTbar"].mass, axis=None)[perm_inds], ctstar=gen_tlep_ctstar[perm_inds], weight=evt_wts[perm_inds])
            acc["GEN_mtt_vs_tlep_ctstar_abs"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype,
                mtt=ak.flatten(gentt["TTbar"].mass, axis=None)[perm_inds], ctstar_abs=np.abs(gen_tlep_ctstar[perm_inds]), weight=evt_wts[perm_inds])
            acc["GEN_mtt_vs_topRapidity"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype,
                mtt=ak.flatten(gentt["TTbar"].mass, axis=None)[perm_inds], top_rap=gen_top_rap[perm_inds], weight=evt_wts[perm_inds])
            acc["GEN_mtt_vs_deltaYtt"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype,
                mtt=ak.flatten(gentt["TTbar"].mass, axis=None)[perm_inds], deltaYtt=gen_deltaYtt[perm_inds], weight=evt_wts[perm_inds])
            acc["GEN_mtt_vs_deltaYtt_abs"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype,
                mtt=ak.flatten(gentt["TTbar"].mass, axis=None)[perm_inds], deltaYtt_abs=np.abs(gen_deltaYtt[perm_inds]), weight=evt_wts[perm_inds])

                # gen-reco reso plots as a function of reco mtt
            acc["RESO_mtt_vs_mtt"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype,
                mtt=ak.flatten(perm["TTbar"].mass, axis=None)[perm_inds],
                reso_mtt=ak.flatten(gentt["TTbar"].mass, axis=None)[perm_inds] - ak.flatten(perm["TTbar"].mass, axis=None)[perm_inds],
                weight=evt_wts[perm_inds])
            acc["RESO_ctstar_vs_mtt"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype,
                mtt=ak.flatten(perm["TTbar"].mass, axis=None)[perm_inds],
                reso_ctstar=gen_tlep_ctstar[perm_inds] - reco_tlep_ctstar[perm_inds], weight=evt_wts[perm_inds])
            acc["RESO_ctstar_abs_vs_mtt"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype,
                mtt=ak.flatten(perm["TTbar"].mass, axis=None)[perm_inds],
                reso_ctstar_abs=np.abs(gen_tlep_ctstar[perm_inds]) - np.abs(reco_tlep_ctstar[perm_inds]), weight=evt_wts[perm_inds])
            acc["RESO_mtt_vs_tlep_ctstar"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype,
                reso_mtt=ak.flatten(gentt["TTbar"].mass, axis=None)[perm_inds] - ak.flatten(perm["TTbar"].mass)[perm_inds],
                reso_ctstar=gen_tlep_ctstar[perm_inds] - reco_tlep_ctstar[perm_inds], weight=evt_wts[perm_inds])
            acc["RESO_mtt_vs_tlep_ctstar_abs"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype,
                reso_mtt=ak.flatten(gentt["TTbar"].mass, axis=None)[perm_inds] - ak.flatten(perm["TTbar"].mass)[perm_inds],
                reso_ctstar_abs=gen_tlep_ctstar[perm_inds] - reco_tlep_ctstar[perm_inds], weight=evt_wts[perm_inds])

                # (reco-gen)/gen reso plots
            acc["REL_RESO_mtt_vs_mtt"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype,
                mtt=ak.flatten(perm["TTbar"].mass, axis=None)[perm_inds],
                rel_reso_mtt=(ak.flatten(perm["TTbar"].mass, axis=None)[perm_inds] - ak.flatten(gentt["TTbar"].mass, axis=None)[perm_inds])/ak.flatten(gentt["TTbar"].mass, axis=None)[perm_inds],
                weight=evt_wts[perm_inds])

        return acc        


    def fill_signal_hists(self, acc, jetmult, leptype, permarray, perm, evt_wts):
            ## apply alpha correction for 3Jets
        if (jetmult == "3Jets") and ("Alpha" in self.corrections):
            alpha_corr = self.corrections["Alpha"](172.5/perm["THad"].mass)
            perm["THad"] = perm["THad"].multiply(alpha_corr) # correct thad
            perm["TTbar"] = ak.flatten(perm["THad"]+perm["TLep"]) # correct ttbar

        thad_ctstar, tlep_ctstar = make_vars.ctstar(perm["THad"], perm["TLep"], flatten=True)

        reco_top_rap, reco_tbar_rap, reco_deltaYtt = make_vars.deltaYtt(perm["Top"], perm["Tbar"], return_indiv=True, flatten=True)

        #set_trace()
        for permval in np.unique(permarray).tolist():
            perm_inds = np.where(permarray == permval)
            dataset_name = "%s_%s" % (self.sample_name, perm_cats[permval]) if permval != 0 else self.sample_name

            acc["RECO_mtt_vs_tlep_ctstar"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype,
                mtt=ak.flatten(perm["TTbar"].mass, axis=None)[perm_inds], ctstar=tlep_ctstar[perm_inds], weight=evt_wts[perm_inds])
            acc["RECO_mtt_vs_tlep_ctstar_abs"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype,
                mtt=ak.flatten(perm["TTbar"].mass, axis=None)[perm_inds], ctstar_abs=np.abs(tlep_ctstar[perm_inds]), weight=evt_wts[perm_inds])
            acc["RECO_mtt_vs_topRapidity"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype,
                mtt=ak.flatten(perm["TTbar"].mass, axis=None)[perm_inds], top_rap=reco_top_rap[perm_inds], weight=evt_wts[perm_inds])
            acc["RECO_mtt_vs_deltaYtt"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype,
                mtt=ak.flatten(perm["TTbar"].mass, axis=None)[perm_inds], deltaYtt=reco_deltaYtt[perm_inds], weight=evt_wts[perm_inds])
            acc["RECO_mtt_vs_deltaYtt_abs"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype,
                mtt=ak.flatten(perm["TTbar"].mass, axis=None)[perm_inds], deltaYtt_abs=np.abs(reco_deltaYtt[perm_inds]), weight=evt_wts[perm_inds])

        return acc        

    def postprocess(self, accumulator):
        return accumulator


proc_executor = processor.iterative_executor if to_debug else processor.futures_executor
proc_exec_args = {"schema": processor.NanoAODSchema} if to_debug else {"schema": processor.NanoAODSchema, "workers": 8}
output = processor.run_uproot_job(
    fileset,
    treename="Events",
    processor_instance=binning_check(),
    executor=proc_executor,
    executor_args=proc_exec_args,
    chunksize=100000,
)

save(output, args.outfname)
print(f"{args.outfname} has been written")

if "isCondor" in opts_dict.keys():
    if ast.literal_eval(opts_dict["isCondor"]):
        print(f"Deleting files from /tmp")
        os.system(f"rm {' '.join(fileset[samplename])}")

toc = time.time()
print("Total time: %.1f" % (toc - tic))
