#!/usr/bin/env python

"""
This analyzer compares the effect of applying the pileup corrections to simulation
"""

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
import python.GenParticleSelector as genpsel
import python.TTPermutator as ttpermutator
import python.IDJet as IDJet
from copy import deepcopy

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
base_jobid = os.environ["base_jobid"]

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
isSingleTop_ = samplename.startswith("singlet")
isTopSample_ = isTTbar_ or isSingleTop_ # ttbar or single top samples
isSignal_ = (samplename.startswith("AtoTT") or samplename.startswith("HtoTT"))
if isSignal_: raise ValueError("This should only be run with data and SM background")
isInt_ = isSignal_ and ("Int" in samplename)
isData_ = samplename.startswith("data")
isSE_Data_ = samplename.startswith("data_SingleElectron")
isSM_Data_ = samplename.startswith("data_SingleMuon")
if isData_:
    lumiMask_path = os.path.join(proj_dir, "inputs", "data", base_jobid, "LumiMasks", f"{args.year}_GoldenJson_{base_jobid}.txt")

Nominal_ttJets = ["ttJetsSL", "ttJetsHad", "ttJetsDiLep"]
isNominal_ttJets_ = samplename in Nominal_ttJets
isTTShift_ = isTTbar_ and (samplename.endswith("UP") or samplename.endswith("DOWN") or "mtop" in samplename)

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
        import subprocess
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
                    output = check_output(["xrdcp", "-f", f"{cp_rfile}", "/tmp"], timeout=600, stderr=STDOUT)
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
    "BTagSF" : not isData_,
    "JetCor" : jet_corrections,
    "Alpha" : alpha_corrections,
    "NNLO_Rewt" : {"Var" : cfg_pars["corrections"]["nnlo"]["var"], "Correction" : nnlo_reweighting[cfg_pars["corrections"]["nnlo"]["var"]]},
    "EWK_Rewt" : {"Correction" : ewk_reweighting, "wt" : cfg_pars["corrections"]["ewk"]["wt"]},
}

    ## parameters for b-tagging
jet_pars = cfg_pars["Jets"]
btaggers = [jet_pars["btagger"]]

wps_to_use = list(set([jet_pars["permutations"]["tightb"],jet_pars["permutations"]["looseb"]]))
if not( len(wps_to_use) == 1):
    raise IOError("Only 1 unique btag working point supported now")
btag_wp = btaggers[0]+wps_to_use[0]

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
reweight_systematics_to_run = ["nosys"] if isData_ else ["nosys", "PileupDown", "PileupUp", "No_Pileup"]
    
print("Running with event systematics:", *sorted(set(event_systematics_to_run).difference(set(["nosys"]))), sep=", ") if "nosys" in event_systematics_to_run else print("Running with event systematics:", *sorted(event_systematics_to_run), sep=", ")
print("\t\tand reweight systematics:", *sorted(reweight_systematics_to_run), sep=", ")
#set_trace()

    # sideband regions are determined by dividing deepcsv medium wp values by 3 for each year
btag_regions = {}

# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class MyAnalyzer(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.sys_axis = hist.Cat("sys", "Systematic")
        self.jetmult_axis = hist.Cat("jmult", "nJets")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")

        self.nTI_axis = hist.Bin("nTI", "nTI", np.around(np.linspace(0., 100., 101), decimals=1))
        self.nPV_axis = hist.Bin("nPV", "nPV", np.around(np.linspace(0., 100., 101), decimals=1))
        self.rho_axis = hist.Bin("rho", "rho", np.around(np.linspace(0., 100., 1001), decimals=1))
        self.wt_axis = hist.Bin("wt", "wt", np.around(np.linspace(0., 5., 101), decimals=2))

            ## make dictionary of hists
        histo_dict = {}
                ## make selection plots
        selection_hists = self.make_selection_hists()
        histo_dict.update(selection_hists)        

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
                },
            },
            "Electron" : {
                "btagPass" : {
                    "3Jets"  : {"lep_and_filter_pass", "passing_jets", "jets_3" , "tight_EL", f"{btaggers[0]}_pass"},
                    "4PJets" : {"lep_and_filter_pass", "passing_jets", "jets_4p", "tight_EL", f"{btaggers[0]}_pass"},
                },
            },
        }
        self.regions = {sys:deepcopy(base_regions) for sys in self.event_systematics_to_run}
            # add sideband regions for nosys
        if "nosys" in self.event_systematics_to_run:
            self.regions["nosys"]["Muon"].update({
                key : {
                    "3Jets"  : {"lep_and_filter_pass", "passing_jets", "jets_3" , "tight_MU", f"{btaggers[0]}_{key}"},
                    "4PJets" : {"lep_and_filter_pass", "passing_jets", "jets_4p", "tight_MU", f"{btaggers[0]}_{key}"},
                } for key in btag_regions.keys()
            })
            self.regions["nosys"]["Electron"].update({
                key : {
                    "3Jets"  : {"lep_and_filter_pass", "passing_jets", "jets_3" , "tight_EL", f"{btaggers[0]}_{key}"},
                    "4PJets" : {"lep_and_filter_pass", "passing_jets", "jets_4p", "tight_EL", f"{btaggers[0]}_{key}"},
                } for key in btag_regions.keys()
            })

        if isSM_Data_:
            if "Electron" in self.regions["nosys"].keys(): del self.regions["nosys"]["Electron"]
        if isSE_Data_:
            if "Muon" in self.regions["nosys"].keys(): del self.regions["nosys"]["Muon"]
        if isData_:
            for lepton in self.regions["nosys"].keys():
                for btagregion in self.regions["nosys"][lepton].keys():
                    for jmult in self.regions["nosys"][lepton][btagregion].keys():
                        self.regions["nosys"][lepton][btagregion][jmult].update({"lumimask"})

        if isTTSL_:
            for sys in self.regions.keys():
                for lepton in self.regions[sys].keys():
                    for btagregion in self.regions[sys][lepton].keys():
                        for jmult in self.regions[sys][lepton][btagregion].keys():
                            self.regions[sys][lepton][btagregion][jmult].update({"semilep"})


    @property
    def accumulator(self):
        return self._accumulator

    def make_selection_hists(self):
        histo_dict = {}
        #histo_dict["nTrueInt"]  = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.nTI_axis)
        histo_dict["nPV"]      = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.nPV_axis)
        histo_dict["rho"]      = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.rho_axis)
        #histo_dict["pu_wt"]      = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.wt_axis)
        #histo_dict["TOT_nTrueInt"]  = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.nTI_axis)
        histo_dict["TOT_nPV"]  = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.nPV_axis)
        histo_dict["TOT_rho"]  = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.rho_axis)
        #histo_dict["TOT_pu_wt"]  = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.wt_axis)

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

        for rewt in self.reweight_systematics_to_run:
            if rewt == "nosys":
                pu_wt = ak.ones_like(events["PV"]["npvs"]) if isData_ else ak.copy(self.corrections["Pileup"][self.sample_name]["central"](events["Pileup"]["nTrueInt"]))
            elif rewt == "PileupUp":
                pu_wt = ak.copy(self.corrections["Pileup"][self.sample_name]["up"](events["Pileup"]["nTrueInt"]))
            elif rewt == "PileupDown":
                pu_wt = ak.copy(self.corrections["Pileup"][self.sample_name]["down"](events["Pileup"]["nTrueInt"]))
            elif rewt == "No_Pileup":
                pu_wt = ak.ones_like(events["PV"]["npvs"])

            evt_wt = ak.ones_like(events["PV"]["npvs"]) if isData_ else ak.copy(events["genWeight"]) * pu_wt

            #output["TOT_pu_wt"].fill(dataset=self.sample_name, sys=rewt, wt=ak.flatten(pu_wt, axis=None))
            #output["TOT_nTrueInt"].fill(dataset=self.sample_name, sys=rewt, nTI=ak.flatten(events["Pileup"]["nTrueInt"], axis=None), weight=pu_wt)
            #output["TOT_nPV"].fill(dataset=self.sample_name, sys=rewt, nPV=ak.flatten(events["PV"]["npvs"], axis=None), weight=pu_wt)
            #output["TOT_rho"].fill(dataset=self.sample_name, sys=rewt, rho=ak.flatten(events["fixedGridRhoFastjetAll"], axis=None), weight=pu_wt)
            #output["TOT_nTrueInt"].fill(dataset=self.sample_name, sys=rewt, nTI=ak.flatten(events["Pileup"]["nTrueInt"], axis=None), weight=evt_wt)
            output["TOT_nPV"].fill(dataset=self.sample_name, sys=rewt, nPV=ak.flatten(events["PV"]["npvs"], axis=None), weight=evt_wt)
            output["TOT_rho"].fill(dataset=self.sample_name, sys=rewt, rho=ak.flatten(events["fixedGridRhoFastjetAll"], axis=None), weight=evt_wt)


            ## initialize selections
        selection = {evt_sys: PackedSelection() for evt_sys in self.event_systematics_to_run}

            # get all passing leptons
        lep_and_filter_pass = objsel.select_leptons(events, year=args.year, hem_15_16=apply_hem)
        {selection[sys].add("lep_and_filter_pass", lep_and_filter_pass) for sys in selection.keys()} # add passing leptons requirement to all systematics

            ## build corrected jets and MET
        events["Jet"], events["MET"] = IDJet.process_jets(events, args.year, self.corrections["JetCor"])

        if isData_:
            runs = events.run
            lumis = events.luminosityBlock
            Golden_Json_LumiMask = lumi_tools.LumiMask(lumiMask_path)
            LumiMask = Golden_Json_LumiMask.__call__(runs, lumis) ## returns array of valid events
            {selection[sys].add("lumimask", LumiMask) for sys in selection.keys()}
                  ## object selection and add different selections
            if isSM_Data_:
                        ## muons
                {selection[sys].add("tight_MU", ak.sum(events["Muon"]["TIGHTMU"], axis=1) == 1) for sys in selection.keys()} # one muon passing TIGHT criteria
            if isSE_Data_:
                        ## electrons
                {selection[sys].add("tight_EL", ak.sum(events["Electron"]["TIGHTEL"], axis=1) == 1) for sys in selection.keys()} # one electron passing TIGHT criteria

        else:
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
                #set_trace()
                    ## NLO EW weights
                if self.corrections["EWK_Rewt"]["wt"] == "Otto":
                    ewk_wts_dict = MCWeights.get_Otto_ewk_weights(self.corrections["EWK_Rewt"]["Correction"], events)
                    mu_evt_weights.add("EWunc",
                        np.ones(len(events)),
                        np.copy(ewk_wts_dict["DeltaQCD"]*ewk_wts_dict["Rebinned_DeltaEW_1.0"]),
                        shift=True
                    )
                    el_evt_weights.add("EWunc",
                        np.ones(len(events)),
                        np.copy(ewk_wts_dict["DeltaQCD"]*ewk_wts_dict["Rebinned_DeltaEW_1.0"]),
                        shift=True
                    )

                        # add Yukawa coupling variation
                    mu_evt_weights.add("Yukawa",  # really just varying value of Yt
                        np.copy(ewk_wts_dict["Rebinned_KFactor_1.0"]),
                        np.copy(ewk_wts_dict["Rebinned_KFactor_1.11"]),
                        np.copy(ewk_wts_dict["Rebinned_KFactor_0.88"]),
                    )
                    el_evt_weights.add("Yukawa",  # really just varying value of Yt
                        np.copy(ewk_wts_dict["Rebinned_KFactor_1.0"]),
                        np.copy(ewk_wts_dict["Rebinned_KFactor_1.11"]),
                        np.copy(ewk_wts_dict["Rebinned_KFactor_0.88"]),
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

                # sort jets by btag value
            events["SelectedJets"] = events["SelectedJets"][ak.argsort(events["SelectedJets"][IDJet.btag_tagger_to_disc_name[btaggers[0]]], ascending=False)]
            btagger_sorted = events["SelectedJets"][ak.argsort(events["SelectedJets"][IDJet.btag_tagger_to_disc_name[btaggers[0]]], ascending=False)][IDJet.btag_tagger_to_disc_name[btaggers[0]]]
            for btag_reg, (low_bound, up_bound) in btag_regions.items():
                selection[evt_sys].add(f"{btaggers[0]}_{btag_reg}", (ak.max(btagger_sorted, axis=1) > low_bound) & (ak.max(btagger_sorted, axis=1) <= up_bound))

                ## apply btagging SFs to MC
            if (not isData_) and (corrections["BTagSF"] == True):
                btag_weights = {key : np.ones(len(events)) for key in self.corrections["BTag_Constructors"][btaggers[0]]["3Jets"].schema_.keys()}

                threeJets_cut = selection[evt_sys].require(lep_and_filter_pass=True, passing_jets=True, jets_3=True)
                btagger_3j_wts = self.corrections["BTag_Constructors"][btaggers[0]]["3Jets"].get_scale_factor(jets=events["SelectedJets"][threeJets_cut], passing_cut=btag_wp)
                fourplusJets_cut = selection[evt_sys].require(lep_and_filter_pass=True, passing_jets=True, jets_4p=True)
                btagger_4pj_wts = self.corrections["BTag_Constructors"][btaggers[0]]["4PJets"].get_scale_factor(jets=events["SelectedJets"][fourplusJets_cut], passing_cut=btag_wp)

                    # fll dict of btag weights
                for wt_name in btagger_3j_wts.keys():
                    btag_weights[wt_name][threeJets_cut] = ak.prod(btagger_3j_wts[wt_name], axis=1)
                    btag_weights[wt_name][fourplusJets_cut] = ak.prod(btagger_4pj_wts[wt_name], axis=1)

            elif (not isData_) and (corrections["BTagSF"] == False):
                raise ValueError("BTag SFs not applied to MC")

            #set_trace()
            ## fill hists for each region
            for lepton in self.regions[evt_sys].keys():
                evt_weights = mu_evt_weights if lepton == "Muon" else el_evt_weights
                if not isData_:
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
                            valid_perms = ak.num(best_perms["TTbar"].pt) > 0
                            output[f"cutflow_{evt_sys}"]["nEvts %s: valid perms" % ", ".join([lepton, btagregion, jmult])] += ak.sum(valid_perms)

                            bp_status = np.zeros(cut.size, dtype=int) # 0 == "" (no gen matching), 1 == "right", 2 == "matchable", 3 == "unmatchable", 4 == "sl_tau", 5 == "noslep"

                                ## create MT regions
                            MT = make_vars.MT(leptons, met)
                            MTHigh = ak.flatten(MT[valid_perms] >= MTcut)
                            output[f"cutflow_{evt_sys}"]["nEvts %s: pass MT cut" % ", ".join([lepton, btagregion, jmult])] += ak.sum(MTHigh)

                                # fill hists for each systematic
                            if to_debug: print(f"\t\tevt sys: {evt_sys}")
                            if evt_sys == "nosys":
                                #set_trace()
                                for rewt_sys in self.reweight_systematics_to_run:
                                    #if to_debug: set_trace()
                                        ## only fill plots in signal region if systematic variation being used
                                    if (rewt_sys != "nosys") and (btagregion != "btagPass"): continue

                                    if to_debug: print(f"\t\tsysname: {rewt_sys}")

                                    #set_trace()
                                    sysname = events.metadata["dataset"].split("_")[-1] if isTTShift_ else rewt_sys

                                    #if "Yukawa" in rewt_sys: set_trace()
                                    if rewt_sys == "nosys":
                                        pu_wt = ak.ones_like(events["PV"]["npvs"])[cut][valid_perms][MTHigh] if isData_ else evt_weights._weights["Pileup"][cut][valid_perms][MTHigh]
                                        wts = evt_weights.weight()[cut][valid_perms][MTHigh] if isData_ else (evt_weights.weight() * btag_weights["central"] * lep_SFs["central"])[cut][valid_perms][MTHigh]
                                    elif rewt_sys == "No_Pileup":
                                        pu_wt = ak.ones_like(events["PV"]["npvs"])[cut][valid_perms][MTHigh]
                                        wts = (evt_weights.partial_weight(exclude=["Pileup"]) * btag_weights["central"] * lep_SFs["central"])[cut][valid_perms][MTHigh]
                                    elif rewt_sys.startswith("btag"):
                                        wts = (evt_weights.weight() * btag_weights[rewt_sys.split("btag_")[-1]] * lep_SFs["central"])[cut][valid_perms][MTHigh]
                                    elif rewt_sys.startswith("Lep"):
                                        #set_trace()
                                        if rewt_sys.split("_")[-1] in lep_SFs.keys():
                                            wts = (evt_weights.weight() * btag_weights["central"] * lep_SFs[rewt_sys.split("_")[-1]])[cut][valid_perms][MTHigh]
                                        else:
                                            print(f"{rewt_sys.split('_')[-1]} not found in {lepton} SF dict. Skipping")
                                            continue
                                    else:
                                        if rewt_sys not in evt_weights.variations:
                                            print(f"{rewt_sys} not option in event weights. Skipping")
                                            continue
                                        pu_wt = (evt_weights.partial_weight(include=["Pileup"])*evt_weights._modifiers[rewt_sys])[cut][valid_perms][MTHigh]
                                        wts = (evt_weights.weight(rewt_sys) * btag_weights["central"] * lep_SFs["central"])[cut][valid_perms][MTHigh]

                                    #output["pu_wt"].fill(dataset=self.sample_name, sys=sysname, jmult=jmult, leptype=lepton, wt=ak.flatten(pu_wt, axis=None))
                                    output["nPV"].fill(dataset=self.sample_name, sys=sysname, jmult=jmult, leptype=lepton, nPV=ak.flatten(events["PV"]["npvs"], axis=None)[cut][valid_perms][MTHigh], weight=wts)
                                    output["rho"].fill(dataset=self.sample_name, sys=sysname, jmult=jmult, leptype=lepton, rho=ak.flatten(events["fixedGridRhoFastjetAll"], axis=None)[cut][valid_perms][MTHigh], weight=wts)

                                    

        return output


    def postprocess(self, accumulator):
        return accumulator

if to_debug:
    proc_executor = processor.iterative_executor
    proc_exec_args = {"schema": processor.NanoAODSchema}
else:
    proc_executor = processor.futures_executor
    proc_exec_args = {
        "schema": processor.NanoAODSchema,
        "workers": 8,
        "merging": True,
    }
output = processor.run_uproot_job(
    fileset,
    treename="Events",
    processor_instance=MyAnalyzer(),
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
