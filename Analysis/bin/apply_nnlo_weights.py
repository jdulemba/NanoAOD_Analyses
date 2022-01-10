#!/usr/bin/env python

import time
tic = time.time()

from coffea import hist, processor
#from coffea.nanoevents import NanoAODSchema
from coffea.analysis_tools import PackedSelection
import awkward as ak
from coffea.nanoevents.methods import vector
ak.behavior.update(vector.behavior)

from coffea.util import save, load
from pdb import set_trace
import os
import python.ObjectSelection as objsel
import python.MCWeights as MCWeights
import numpy as np
import Utilities.prettyjson as prettyjson
import Utilities.make_variables as make_vars
from python.IDJet import btag_values as btag_values
import python.GenParticleSelector as genpsel
import python.TTGenMatcher as ttmatcher
import python.TTPermutator as ttpermutator
from python.Permutations import compare_matched_best_perms
import python.IDJet as IDJet
from copy import deepcopy

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
base_jobid = os.environ["base_jobid"]
analyzer = "apply_nnlo_weights"

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("fset", type=str, help="Fileset dictionary (in string form) to be used for the processor")
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"] if base_jobid == "ULnanoAOD" else ["2016", "2017", "2018"], help="Specify which year to run over")
parser.add_argument("outfname", type=str, help="Specify output filename, including directory and file extension")
parser.add_argument("opts", type=str, help="Fileset dictionary (in string form) to be used for the processor")
args = parser.parse_args()

# convert input string of fileset dictionary to actual dictionary
fdict = (args.fset).replace("\'", "\"")
fileset = prettyjson.loads(fdict)

# convert input string of options dictionary to actual dictionary
odict = (args.opts).replace("\'", "\"")
opts_dict = prettyjson.loads(odict)

    ## set config options passed through argparse
import ast
to_debug = ast.literal_eval(opts_dict.get("debug", "False"))


if len(fileset.keys()) > 1:
    raise ValueError("Only one topology run at a time in order to determine which corrections and systematics to run")
samplename = list(fileset.keys())[0]

# get dataset classification, used for corrections/systematics
    ## specify ttJets samples
Nominal_ttJets = ["ttJets_PS"] if ((args.year == "2016") and (base_jobid == "NanoAODv6")) else ["ttJetsSL", "ttJetsHad", "ttJetsDiLep"]
isTTbar_ = samplename.startswith("ttJets")
isTTSL_ = samplename.startswith("ttJetsSL")
isNominal_ttJets_ = samplename in Nominal_ttJets
if not isNominal_ttJets_:
    raise ValueError("Only nominal ttbar events should be run")

## init tt probs for likelihoods
ttpermutator.year_to_run(year=args.year)

## load lumimask for data and corrections for event weights
pu_correction = load(os.path.join(proj_dir, "Corrections", base_jobid, "MC_PU_Weights.coffea"))[args.year]
lepSF_correction = load(os.path.join(proj_dir, "Corrections", base_jobid, "leptonSFs.coffea"))[args.year]
jet_corrections = load(os.path.join(proj_dir, "Corrections", base_jobid, "JetMETCorrections.coffea"))[args.year]
alpha_corrections = load(os.path.join(proj_dir, "Corrections", base_jobid,"alpha_correction_%s.coffea" % jobid))[args.year]["E"]["All_2D"] # E, All_2D determined by post application plots
nnlo_reweighting = load(os.path.join(proj_dir, "Corrections", base_jobid, f"NNLO_to_Tune_Ratios_{base_jobid}_Test.coffea"))[args.year]
#nnlo_reweighting = load(os.path.join(proj_dir, "Corrections", base_jobid, "NNLO_to_Tune_Orig_Interp_Ratios_%s.coffea" % base_jobid))[args.year]
corrections = {
    "Pileup" : pu_correction,
    "Prefire" : True,
    "LeptonSF" : lepSF_correction,
    "BTagSF" : True,
    "JetCor" : jet_corrections,
    "Alpha" : alpha_corrections,
    "NNLO_Rewt" : nnlo_reweighting,
}

cfg_pars = prettyjson.loads(open(os.path.join(proj_dir, "cfg_files", "cfg_pars_%s.json" % jobid)).read())
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
        raise IOError("BTag SF file %s doesn't exist" % sf_file)

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

    

# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class apply_nnlo_weights(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        #self.sys_axis = hist.Cat("sys", "Systematic")
        self.evt_type_axis = hist.Cat("evt_type", "Event Type")
        self.nnlo_reweighting_axis = hist.Cat("rewt", "NNLO Reweighting Type")
        self.jetmult_axis = hist.Cat("jmult", "nJets")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")
        self.lepcat_axis = hist.Cat("lepcat", "Lepton Category")
        self.btag_axis = hist.Cat("btag", "BTag Region")
        self.pt_axis = hist.Bin("pt", "p_{T} [GeV]", 200, 0, 1000)
        self.eta_axis = hist.Bin("eta", r"$\eta$", 200, -5., 5.)
        self.mtop_axis = hist.Bin("mtop", "m(top) [GeV]", 300, 0, 300)
        self.mtt_axis = hist.Bin("mtt", "m($t\overline{t}$) [GeV]", 360, 200, 2000)
        self.ctstar_axis = hist.Bin("ctstar", "cos($\\theta^{*}$)", 200, -1., 1.)
        #self.ctstar_axis = hist.Bin("ctstar", "cos($\\theta^{*}$)", 40, -1., 1.)
        self.ctstar_abs_axis = hist.Bin("ctstar_abs", "|cos($\\theta^{*}$)|", 200, 0., 1.)
        #self.ctstar_abs_axis = hist.Bin("ctstar_abs", "|cos($\\theta^{*}$)|", 20, 0., 1.)
        self.reso_pt_axis         = hist.Bin("reso_pt", "p_{T} [GeV]", 200, -200, 200)
        self.reso_eta_axis        = hist.Bin("reso_eta", r"$\eta$", 200, -5., 5.)
        self.reso_mtop_axis       = hist.Bin("reso_mtop", "m(top) [GeV]", 400, -200, 200)
        self.reso_mtt_axis        = hist.Bin("reso_mtt", "m($t\overline{t}$) [GeV]", 400, -200, 200)
        self.reso_ctstar_axis     = hist.Bin("reso_ctstar", "cos($\\theta^{*}$)", 400, -2., 2.)
        self.reso_ctstar_abs_axis = hist.Bin("reso_ctstar_abs", "|cos($\\theta^{*}$)|", 200, -1., 1.)
        self.nnlo_wts_axis = hist.Bin("nnlo_wts", "NNLO Weight", 200, 0., 2.)

            ## make dictionary of hists
        histo_dict = {}
                ## make semilep plots
        semilep_hists = self.make_semilep_hists()
        histo_dict.update(semilep_hists)        
                ## make dilep/had plots
        dl_had_hists = self.make_dl_had_hists()
        histo_dict.update(dl_had_hists)        

        self.sample_name = ""
        self.corrections = corrections

            ## make dict of cutflow for each systematic variation
        histo_dict["cutflow_nosys"] = processor.defaultdict_accumulator(int)
    
        self._accumulator = processor.dict_accumulator(histo_dict)

            ## define regions for MC
        base_regions = {
            "Muon" : {
                "Tight" : {
                    "btagPass" : {
                        "3Jets"  : {"lep_and_filter_pass", "passing_jets", "jets_3" , "tight_MU", "DeepCSV_pass"},
                        "4PJets" : {"lep_and_filter_pass", "passing_jets", "jets_4p", "tight_MU", "DeepCSV_pass"},
                    },
                },
            },
            "Electron" : {
                "Tight" : {
                    "btagPass" : {
                        "3Jets"  : {"lep_and_filter_pass", "passing_jets", "jets_3" , "tight_EL", "DeepCSV_pass"},
                        "4PJets" : {"lep_and_filter_pass", "passing_jets", "jets_4p", "tight_EL", "DeepCSV_pass"},
                    },
                },
            },
        }
        self.regions = deepcopy(base_regions)
        if isTTSL_:
            for lepton in self.regions.keys():
                for leptype in self.regions[lepton].keys():
                    for btagregion in self.regions[lepton][leptype].keys():
                        for jmult in self.regions[lepton][leptype][btagregion].keys():
                            self.regions[lepton][leptype][btagregion][jmult].update({"semilep"})

    @property
    def accumulator(self):
        return self._accumulator

    def make_semilep_hists(self):
        histo_dict = {}
            # gen/reco hists
                ## plots for cross checking
        histo_dict["mtt_vs_nnlo_weights"]= hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mtt_axis, self.nnlo_wts_axis)

        histo_dict["mtt"]      = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mtt_axis)
        histo_dict["mthad"]    = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mtop_axis)
        histo_dict["pt_thad"]  = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.pt_axis)
        histo_dict["pt_tlep"]  = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.pt_axis)
        histo_dict["pt_tt"]    = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.pt_axis)
        histo_dict["eta_thad"] = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.eta_axis)
        histo_dict["eta_tlep"] = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.eta_axis)
        histo_dict["eta_tt"]   = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.eta_axis)

        histo_dict["tlep_ctstar"]     = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.ctstar_axis)
        histo_dict["tlep_ctstar_abs"] = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.ctstar_abs_axis)

        histo_dict["mtt_vs_tlep_ctstar_abs"] = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mtt_axis, self.ctstar_abs_axis)

            # gen-reco resolution hists
        histo_dict["Reso_mtt"]      = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_mtt_axis)
        histo_dict["Reso_mthad"]    = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_mtop_axis)
        histo_dict["Reso_pt_thad"]  = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_pt_axis)
        histo_dict["Reso_pt_tlep"]  = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_pt_axis)
        histo_dict["Reso_pt_tt"]    = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_pt_axis)
        histo_dict["Reso_eta_thad"] = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_eta_axis)
        histo_dict["Reso_eta_tlep"] = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_eta_axis)
        histo_dict["Reso_eta_tt"]   = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_eta_axis)

        histo_dict["Reso_tlep_ctstar"]     = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_ctstar_axis)
        histo_dict["Reso_tlep_ctstar_abs"] = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_ctstar_abs_axis)

        histo_dict["Reso_mtt_vs_tlep_ctstar_abs"] = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_mtt_axis, self.reso_ctstar_abs_axis)

        return histo_dict

    def make_dl_had_hists(self):
        histo_dict = {}
            # gen/reco hists
        histo_dict["mtop"]    = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mtop_axis)
        histo_dict["mtbar"]   = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mtop_axis)
        histo_dict["pt_top"]  = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.pt_axis)
        histo_dict["pt_tbar"] = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.pt_axis)
        histo_dict["eta_top"] = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.eta_axis)
        histo_dict["eta_tbar"]= hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.eta_axis)

        histo_dict["top_ctstar"]     = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.ctstar_axis)
        histo_dict["tbar_ctstar"]    = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.ctstar_axis)
        histo_dict["top_ctstar_abs"] = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.ctstar_abs_axis)
        histo_dict["tbar_ctstar_abs"]= hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.ctstar_abs_axis)

        histo_dict["mtt_vs_top_ctstar_abs"] = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mtt_axis, self.ctstar_abs_axis)
        histo_dict["mtt_vs_tbar_ctstar_abs"]= hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.evt_type_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mtt_axis, self.ctstar_abs_axis)

            # gen-reco resolution hists
        histo_dict["Reso_mtop"]    = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_mtop_axis)
        histo_dict["Reso_mtbar"]   = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_mtop_axis)
        histo_dict["Reso_pt_top"]  = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_pt_axis)
        histo_dict["Reso_pt_tbar"] = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_pt_axis)
        histo_dict["Reso_eta_top"] = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_eta_axis)
        histo_dict["Reso_eta_tbar"]= hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_eta_axis)

        histo_dict["Reso_top_ctstar"]     = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_ctstar_axis)
        histo_dict["Reso_tbar_ctstar"]    = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_ctstar_axis)
        histo_dict["Reso_top_ctstar_abs"] = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_ctstar_abs_axis)
        histo_dict["Reso_tbar_ctstar_abs"]= hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_ctstar_abs_axis)

        histo_dict["Reso_mtt_vs_top_ctstar_abs"] = hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_mtt_axis, self.reso_ctstar_abs_axis)
        histo_dict["Reso_mtt_vs_tbar_ctstar_abs"]= hist.Hist("Events", self.dataset_axis, self.nnlo_reweighting_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.reso_mtt_axis, self.reso_ctstar_abs_axis)

        return histo_dict

    def process(self, events):
        np.random.seed(10) # sets seed so values from random distributions are reproducible (JER corrections)
        output = self.accumulator.identity()

        self.sample_name = events.metadata["dataset"]

            ## make event weights
                # data or MC distinction made internally
        mu_evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections)
        el_evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections)

            ## initialize selections
        selection = PackedSelection()

            # find gen level particles for ttbar system and other ttbar corrections
        genpsel.select(events, mode="NORMAL")
        selection.add("semilep", ak.num(events["SL"]) > 0)
        if "NNLO_Rewt" in self.corrections.keys():
            nnlo_weights = {var : MCWeights.get_nnlo_weights({"Var": var, "Correction" : self.corrections["NNLO_Rewt"][var]}, events) for var in self.corrections["NNLO_Rewt"].keys()}
            nnlo_weights.update({"Nominal" : np.ones(len(events))})

            # get all passing leptons
        lep_and_filter_pass = objsel.select_leptons(events, year=args.year)
        output["cutflow_nosys"]["lep_and_filter_pass"] += ak.sum(lep_and_filter_pass)
        selection.add("lep_and_filter_pass", lep_and_filter_pass)

        ## add different selections
                ## muons
        tight_mu_sel, loose_mu_sel = ak.sum(events["Muon"]["TIGHTMU"], axis=1) == 1, ak.sum(events["Muon"]["LOOSEMU"], axis=1) == 1
        selection.add("tight_MU", tight_mu_sel) # one muon passing TIGHT criteria
        selection.add("loose_MU", loose_mu_sel) # one muon passing LOOSE criteria
                ## electrons
        tight_el_sel, loose_el_sel = ak.sum(events["Electron"]["TIGHTEL"], axis=1) == 1, ak.sum(events["Electron"]["LOOSEEL"], axis=1) == 1
        selection.add("tight_EL", tight_el_sel) # one electron passing TIGHT criteria
        selection.add("loose_EL", loose_el_sel) # one electron passing LOOSE criteria

        ### apply lepton SFs to MC (only applicable to tight leptons)
        if "LeptonSF" in corrections.keys():
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

            ## build corrected jets and MET
        events["Jet"], events["MET"] = IDJet.process_jets(events, args.year, self.corrections["JetCor"])

            # jet selection
        passing_jets = objsel.jets_selection(events, year=args.year, cutflow=output["cutflow_nosys"])
        output["cutflow_nosys"]["nEvts passing jet and lepton obj selection"] += ak.sum(passing_jets & lep_and_filter_pass)
        selection.add("passing_jets", passing_jets)
        selection.add("jets_3",  ak.num(events["SelectedJets"]) == 3)
        selection.add("jets_4p",  ak.num(events["SelectedJets"]) > 3) # only for getting btag weights
        selection.add("DeepCSV_pass", ak.sum(events["SelectedJets"][btag_wps[0]], axis=1) >= 2)

            # sort jets by btag value
        events["SelectedJets"] = events["SelectedJets"][ak.argsort(events["SelectedJets"]["btagDeepB"], ascending=False)] if btaggers[0] == "DeepCSV" else events["SelectedJets"][ak.argsort(events["SelectedJets"]["btagDeepFlavB"], ascending=False)]

            # btag fail sideband
        deepcsv_sorted = events["SelectedJets"][ak.argsort(events["SelectedJets"]["btagDeepB"], ascending=False)]["btagDeepB"]
        valid_counts_inds = ak.where(ak.num(events["SelectedJets"]) > 1)[0]
        deepcsv_fail = np.zeros(len(events)).astype(bool)
        deepcsv_fail[valid_counts_inds] = (deepcsv_sorted[valid_counts_inds][:, 0] < btag_values[args.year]["btagDeepB"]["DeepCSV"+wps_to_use[0]]) & (deepcsv_sorted[valid_counts_inds][:, 1] < btag_values[args.year]["btagDeepB"]["DeepCSV"+wps_to_use[0]])
        selection.add("DeepCSV_fail", deepcsv_fail ) # highest and second highest DeepCSV values don"t pass tight and loose WPs

            ## apply btagging SFs to MC
        if corrections["BTagSF"] == True:
            deepcsv_cen   = np.ones(len(events))

            threeJets_cut = selection.require(lep_and_filter_pass=True, passing_jets=True, jets_3=True)
            deepcsv_3j_wts = self.corrections["BTag_Constructors"]["DeepCSV"]["3Jets"].get_scale_factor(jets=events["SelectedJets"][threeJets_cut], passing_cut="DeepCSV"+wps_to_use[0])
            deepcsv_cen[threeJets_cut] = ak.prod(deepcsv_3j_wts["central"], axis=1)

            fourplusJets_cut = selection.require(lep_and_filter_pass=True, passing_jets=True, jets_4p=True)
            deepcsv_4pj_wts = self.corrections["BTag_Constructors"]["DeepCSV"]["4PJets"].get_scale_factor(jets=events["SelectedJets"][fourplusJets_cut], passing_cut="DeepCSV"+wps_to_use[0])
            deepcsv_cen[fourplusJets_cut] = ak.prod(deepcsv_4pj_wts["central"], axis=1)


        #set_trace()
        ## fill hists for each region
        for lepton in self.regions.keys():
            evt_weights = mu_evt_weights if lepton == "Muon" else el_evt_weights
            for leptype in self.regions[lepton].keys():
                for btagregion in self.regions[lepton][leptype].keys():
                    for jmult in self.regions[lepton][leptype][btagregion].keys():
                        cut = selection.all(*self.regions[lepton][leptype][btagregion][jmult])
                        #set_trace()

                        output["cutflow_nosys"]["nEvts %s" % ", ".join([lepton, leptype, btagregion, jmult])] += cut.sum()

                        if to_debug: print(lepton, leptype, btagregion, jmult)
                        if cut.sum() > 0:
                            ltype = "MU" if lepton == "Muon" else "EL"
                            if "loose_or_tight_%s" % ltype in self.regions[lepton][leptype][btagregion][jmult]:
                                leptons = events[lepton][cut][((events[lepton][cut]["TIGHT%s" % ltype] == True) | (events[lepton][cut]["LOOSE%s" % ltype] == True))]
                            elif "tight_%s" % ltype in self.regions[lepton][leptype][btagregion][jmult]:
                                leptons = events[lepton][cut][(events[lepton][cut]["TIGHT%s" % ltype] == True)]
                            elif "loose_%s" % ltype in self.regions[lepton][leptype][btagregion][jmult]:
                                leptons = events[lepton][cut][(events[lepton][cut]["LOOSE%s" % ltype] == True)]
                            else:
                                raise ValueError("Not sure what lepton type to choose for event")

                                # get jets and MET
                            jets, met = events["SelectedJets"][cut], events["SelectedMET"][cut]

                                # find best permutations
                            best_perms = ttpermutator.find_best_permutations(jets=jets, leptons=leptons, MET=met, btagWP=btag_wps[0], btag_req=False if btagregion == "btagFail" else True)
                            valid_perms = ak.num(best_perms["TTbar"].pt) > 0
                            output["cutflow_nosys"]["nEvts %s: valid perms" % ", ".join([lepton, leptype, btagregion, jmult])] += ak.sum(valid_perms)

                            bp_status = np.zeros(cut.size, dtype=int) # 0 == "" (no gen matching), 1 == "right", 2 == "matchable", 3 == "unmatchable", 4 == "sl_tau", 5 == "noslep"
                                # get matched permutation (semilep ttbar only)
                            if isTTbar_:
                                semilep_evts = selection.require(semilep=True)
                                bp_status[~semilep_evts] = 5
                                if semilep_evts.sum() > 0:
                                        # find matched permutations
                                    mp = ttmatcher.best_match(gen_hyp=events["SL"][cut], jets=jets, leptons=leptons, met=met)
                                    perm_cat_array = compare_matched_best_perms(mp, best_perms, njets=jmult)
                                    bp_status[cut] = perm_cat_array
                                        # consider tau events as separate
                                    if ak.any(ak.num(events["SL"]["Lepton"].pdgId) != 1): raise ValueError("Number of leptons is incorrect for classifying tau+jets events")
                                    sl_tau_evts = ak.where(np.abs(events["SL"]["Lepton"].pdgId) == 15)[0]
                                    bp_status[sl_tau_evts] = 4

                                ## create MT regions
                            MT = make_vars.MT(leptons, met)
                            MTHigh = ak.flatten(MT[valid_perms] >= MTcut)
                            output["cutflow_nosys"]["nEvts %s: pass MT cut" % ", ".join([lepton, leptype, btagregion, jmult])] += ak.sum(MTHigh)

                                # make final objects/values after applying all cuts
                            evt_wts = (evt_weights.weight()*deepcsv_cen)[cut][valid_perms][MTHigh]
                            final_bps = best_perms[valid_perms][MTHigh]
                            final_bp_status = bp_status[cut][valid_perms][MTHigh]

                            sl_evts = ak.num(events["SL"][cut][valid_perms][MTHigh]["TTbar"]) > 0
                            dl_evts = ak.num(events["DL"][cut][valid_perms][MTHigh]["TTbar"]) > 0
                            had_evts = ak.num(events["Had"][cut][valid_perms][MTHigh]["TTbar"]) > 0

                            for rewt_type in nnlo_weights.keys():
                                final_nnlo_wts = nnlo_weights[rewt_type][cut][valid_perms][MTHigh]
                                final_evt_wts = evt_wts*final_nnlo_wts
                                    # fill hists
                                if ak.sum(sl_evts) > 0:
                                    output = self.fill_semilep_hists(acc=output, rewt_type=rewt_type, jetmult=jmult, leptype=lepton, lepcat=leptype, btagregion=btagregion, permarray=final_bp_status[sl_evts],
                                        perm=final_bps[sl_evts], gentt=events["SL"][cut][valid_perms][MTHigh][sl_evts], evt_weights=final_evt_wts[sl_evts], nnlo_weights=final_nnlo_wts[sl_evts])
                                if ak.sum(dl_evts) > 0:
                                    output = self.fill_dilep_had_hists(acc=output, rewt_type=rewt_type, jetmult=jmult, leptype=lepton, lepcat=leptype, btagregion=btagregion, permarray=final_bp_status[dl_evts],
                                        perm=final_bps[dl_evts], gentt=events["DL"][cut][valid_perms][MTHigh][dl_evts], evt_weights=final_evt_wts[dl_evts], nnlo_weights=final_nnlo_wts[dl_evts])
                                if ak.sum(had_evts) > 0:
                                    output = self.fill_dilep_had_hists(acc=output, rewt_type=rewt_type, jetmult=jmult, leptype=lepton, lepcat=leptype, btagregion=btagregion, permarray=final_bp_status[had_evts],
                                        perm=final_bps[had_evts], gentt=events["Had"][cut][valid_perms][MTHigh][had_evts], evt_weights=final_evt_wts[had_evts], nnlo_weights=final_nnlo_wts[had_evts])

        return output


    def fill_semilep_hists(self, acc, rewt_type, jetmult, leptype, lepcat, btagregion, permarray, perm, gentt, evt_weights, nnlo_weights):
            ## apply alpha correction for 3Jets
        #set_trace()
        if (jetmult == "3Jets") and ("Alpha" in self.corrections):
            alpha_corr = self.corrections["Alpha"](172.5/perm["THad"].mass)
            perm["THad"] = perm["THad"].multiply(alpha_corr) # correct thad
            perm["TTbar"] = ak.flatten(perm["THad"]+perm["TLep"]) # correct ttbar

        reco_thad_ctstar, reco_tlep_ctstar = make_vars.ctstar(perm["THad"], perm["TLep"])
        reco_thad_ctstar, reco_tlep_ctstar = ak.flatten(reco_thad_ctstar, axis=None), ak.flatten(reco_tlep_ctstar, axis=None)

        gen_thad_ctstar, gen_tlep_ctstar = make_vars.ctstar(gentt["THad"], gentt["TLep"])
        gen_thad_ctstar, gen_tlep_ctstar = ak.flatten(gen_thad_ctstar, axis=None), ak.flatten(gen_tlep_ctstar, axis=None)

        for permval in np.unique(permarray).tolist():
            perm_inds = np.where(permarray == permval)
            dataset_name = "%s_%s" % (self.sample_name, perm_cats[permval]) if permval != 0 else self.sample_name

                ## fill reco-level quantities
            acc["mtt"].fill(     dataset=dataset_name, rewt=rewt_type, evt_type="RECO", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtt=ak.flatten(perm["TTbar"].mass)[perm_inds], weight=evt_weights[perm_inds])
            acc["mthad"].fill(   dataset=dataset_name, rewt=rewt_type, evt_type="RECO", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtop=ak.flatten(perm["THad"].mass)[perm_inds], weight=evt_weights[perm_inds])
            acc["pt_thad"].fill( dataset=dataset_name, rewt=rewt_type, evt_type="RECO", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=ak.flatten(perm["THad"].pt)[perm_inds], weight=evt_weights[perm_inds])
            acc["pt_tlep"].fill( dataset=dataset_name, rewt=rewt_type, evt_type="RECO", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=ak.flatten(perm["TLep"].pt)[perm_inds], weight=evt_weights[perm_inds])
            acc["pt_tt"].fill(   dataset=dataset_name, rewt=rewt_type, evt_type="RECO", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=ak.flatten(perm["TTbar"].pt)[perm_inds], weight=evt_weights[perm_inds])
            acc["eta_thad"].fill(dataset=dataset_name, rewt=rewt_type, evt_type="RECO", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=ak.flatten(perm["THad"].eta)[perm_inds], weight=evt_weights[perm_inds])
            acc["eta_tlep"].fill(dataset=dataset_name, rewt=rewt_type, evt_type="RECO", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=ak.flatten(perm["TLep"].eta)[perm_inds], weight=evt_weights[perm_inds])
            acc["eta_tt"].fill(  dataset=dataset_name, rewt=rewt_type, evt_type="RECO", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=ak.flatten(perm["TTbar"].eta)[perm_inds], weight=evt_weights[perm_inds])

            acc["tlep_ctstar"].fill(    dataset=dataset_name, rewt=rewt_type, evt_type="RECO", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, ctstar=reco_tlep_ctstar[perm_inds], weight=evt_weights[perm_inds])
            acc["tlep_ctstar_abs"].fill(dataset=dataset_name, rewt=rewt_type, evt_type="RECO", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, ctstar_abs=np.abs(reco_tlep_ctstar[perm_inds]), weight=evt_weights[perm_inds])

            acc["mtt_vs_tlep_ctstar_abs"].fill(dataset=dataset_name, rewt=rewt_type, evt_type="RECO", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion,
                mtt=ak.flatten(perm["TTbar"].mass)[perm_inds], ctstar_abs=np.abs(reco_tlep_ctstar[perm_inds]), weight=evt_weights[perm_inds])

            acc["mtt_vs_nnlo_weights"].fill(dataset=dataset_name, rewt=rewt_type, evt_type="RECO", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtt=ak.flatten(perm["TTbar"].mass)[perm_inds], nnlo_wts=nnlo_weights[perm_inds])

                ## fill gen-level quantities
            acc["mtt"].fill(     dataset=dataset_name, rewt=rewt_type, evt_type="GEN", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtt=ak.flatten(gentt["TTbar"].mass, axis=None)[perm_inds], weight=evt_weights[perm_inds])
            acc["mthad"].fill(   dataset=dataset_name, rewt=rewt_type, evt_type="GEN", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtop=ak.flatten(gentt["THad"].mass, axis=None)[perm_inds], weight=evt_weights[perm_inds])
            acc["pt_thad"].fill( dataset=dataset_name, rewt=rewt_type, evt_type="GEN", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=ak.flatten(gentt["THad"].pt, axis=None)[perm_inds], weight=evt_weights[perm_inds])
            acc["pt_tlep"].fill( dataset=dataset_name, rewt=rewt_type, evt_type="GEN", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=ak.flatten(gentt["TLep"].pt, axis=None)[perm_inds], weight=evt_weights[perm_inds])
            acc["pt_tt"].fill(   dataset=dataset_name, rewt=rewt_type, evt_type="GEN", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=ak.flatten(gentt["TTbar"].pt, axis=None)[perm_inds], weight=evt_weights[perm_inds])
            acc["eta_thad"].fill(dataset=dataset_name, rewt=rewt_type, evt_type="GEN", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=ak.flatten(gentt["THad"].eta, axis=None)[perm_inds], weight=evt_weights[perm_inds])
            acc["eta_tlep"].fill(dataset=dataset_name, rewt=rewt_type, evt_type="GEN", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=ak.flatten(gentt["TLep"].eta, axis=None)[perm_inds], weight=evt_weights[perm_inds])
            acc["eta_tt"].fill(  dataset=dataset_name, rewt=rewt_type, evt_type="GEN", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=ak.flatten(gentt["TTbar"].eta, axis=None)[perm_inds], weight=evt_weights[perm_inds])

            acc["tlep_ctstar"].fill(    dataset=dataset_name, rewt=rewt_type, evt_type="GEN", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, ctstar=gen_tlep_ctstar[perm_inds], weight=evt_weights[perm_inds])
            acc["tlep_ctstar_abs"].fill(dataset=dataset_name, rewt=rewt_type, evt_type="GEN", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, ctstar_abs=np.abs(gen_tlep_ctstar[perm_inds]), weight=evt_weights[perm_inds])

            acc["mtt_vs_tlep_ctstar_abs"].fill(dataset=dataset_name, rewt=rewt_type, evt_type="GEN", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion,
                mtt=ak.flatten(gentt["TTbar"].mass, axis=None)[perm_inds], ctstar_abs=np.abs(gen_tlep_ctstar[perm_inds]), weight=evt_weights[perm_inds])

            acc["mtt_vs_nnlo_weights"].fill(dataset=dataset_name, rewt=rewt_type, evt_type="GEN", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtt=ak.flatten(gentt["TTbar"].mass, axis=None)[perm_inds], nnlo_wts=nnlo_weights[perm_inds])

                ## fill gen-reco quantities
            acc["Reso_mtt"].fill(dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion,
                reso_mtt=ak.flatten(gentt["TTbar"].mass, axis=None)[perm_inds]-ak.flatten(perm["TTbar"].mass)[perm_inds], weight=evt_weights[perm_inds])
            acc["Reso_mthad"].fill(dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion,
                reso_mtop=ak.flatten(gentt["THad"].mass, axis=None)[perm_inds]-ak.flatten(perm["THad"].mass)[perm_inds], weight=evt_weights[perm_inds])
            acc["Reso_pt_thad"].fill(dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion,
                reso_pt=ak.flatten(gentt["THad"].pt, axis=None)[perm_inds]-ak.flatten(perm["THad"].pt)[perm_inds], weight=evt_weights[perm_inds])
            acc["Reso_pt_tlep"].fill(dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion,
                reso_pt=ak.flatten(gentt["TLep"].pt, axis=None)[perm_inds]-ak.flatten(perm["TLep"].pt)[perm_inds], weight=evt_weights[perm_inds])
            acc["Reso_pt_tt"].fill(dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion,
                reso_pt=ak.flatten(gentt["TTbar"].pt, axis=None)[perm_inds]-ak.flatten(perm["TTbar"].pt)[perm_inds], weight=evt_weights[perm_inds])
            acc["Reso_eta_thad"].fill(dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion,
                reso_eta=ak.flatten(gentt["THad"].eta, axis=None)[perm_inds]-ak.flatten(perm["THad"].eta)[perm_inds], weight=evt_weights[perm_inds])
            acc["Reso_eta_tlep"].fill(dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion,
                reso_eta=ak.flatten(gentt["TLep"].eta, axis=None)[perm_inds]-ak.flatten(perm["TLep"].eta)[perm_inds], weight=evt_weights[perm_inds])
            acc["Reso_eta_tt"].fill(dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion,
                reso_eta=ak.flatten(gentt["TTbar"].eta, axis=None)[perm_inds]-ak.flatten(perm["TTbar"].eta)[perm_inds], weight=evt_weights[perm_inds])

            acc["Reso_tlep_ctstar"].fill(dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion,
                reso_ctstar=gen_tlep_ctstar[perm_inds]-reco_tlep_ctstar[perm_inds], weight=evt_weights[perm_inds])
            acc["Reso_tlep_ctstar_abs"].fill(dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion,
                reso_ctstar_abs=np.abs(gen_tlep_ctstar[perm_inds])-np.abs(reco_tlep_ctstar[perm_inds]), weight=evt_weights[perm_inds])

            acc["Reso_mtt_vs_tlep_ctstar_abs"].fill(dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion,
                reso_mtt=ak.flatten(gentt["TTbar"].mass, axis=None)[perm_inds]-ak.flatten(perm["TTbar"].mass)[perm_inds], reso_ctstar_abs=np.abs(gen_tlep_ctstar[perm_inds])-np.abs(reco_tlep_ctstar[perm_inds]), weight=evt_weights[perm_inds])

        return acc        

    def fill_dilep_had_hists(self, acc, rewt_type, jetmult, leptype, lepcat, btagregion, permarray, perm, gentt, evt_weights, nnlo_weights):
            ## apply alpha correction for 3Jets
        if (jetmult == "3Jets") and ("Alpha" in self.corrections):
            alpha_corr = self.corrections["Alpha"](172.5/perm["THad"].mass)
            perm["THad"] = perm["THad"].multiply(alpha_corr) # correct thad
            perm["TTbar"] = ak.flatten(perm["THad"]+perm["TLep"]) # correct ttbar

        reco_thad_ctstar, reco_tlep_ctstar = make_vars.ctstar(perm["THad"], perm["TLep"])
        reco_thad_ctstar, reco_tlep_ctstar = ak.flatten(reco_thad_ctstar, axis=None), ak.flatten(reco_tlep_ctstar, axis=None)
            # convert thad/tlep into top/tbar for gen comparison
        reco_top_mass = ak.flatten(ak.where(perm["Lepton"].charge == 1, perm["TLep"].mass, perm["THad"].mass), axis=None)
        reco_tbar_mass = ak.flatten(ak.where(perm["Lepton"].charge == -1, perm["TLep"].mass, perm["THad"].mass), axis=None)
        reco_top_pt = ak.flatten(ak.where(perm["Lepton"].charge == 1, perm["TLep"].pt, perm["THad"].pt), axis=None)
        reco_tbar_pt = ak.flatten(ak.where(perm["Lepton"].charge == -1, perm["TLep"].pt, perm["THad"].pt), axis=None)
        reco_top_eta = ak.flatten(ak.where(perm["Lepton"].charge == 1, perm["TLep"].eta, perm["THad"].eta), axis=None)
        reco_tbar_eta = ak.flatten(ak.where(perm["Lepton"].charge == -1, perm["TLep"].eta, perm["THad"].eta), axis=None)
        reco_top_ctstar = ak.flatten(ak.where(perm["Lepton"].charge == 1, reco_tlep_ctstar, reco_thad_ctstar))
        reco_tbar_ctstar = ak.flatten(ak.where(perm["Lepton"].charge == -1, reco_tlep_ctstar, reco_thad_ctstar))

        gen_top_ctstar, gen_tbar_ctstar = make_vars.ctstar(gentt["Top"], gentt["Tbar"])
        gen_top_ctstar, gen_tbar_ctstar = ak.flatten(gen_top_ctstar, axis=None), ak.flatten(gen_tbar_ctstar, axis=None)

        for permval in np.unique(permarray).tolist():
            perm_inds = np.where(permarray == permval)
            dataset_name = "%s_%s" % (self.sample_name, perm_cats[permval]) if permval != 0 else self.sample_name

                ## fill reco-level quantities
            acc["mtt"].fill(     dataset=dataset_name, rewt=rewt_type, evt_type="RECO", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtt=ak.flatten(perm["TTbar"].mass)[perm_inds], weight=evt_weights[perm_inds])
            acc["mtop"].fill(    dataset=dataset_name, rewt=rewt_type, evt_type="RECO", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtop=reco_top_mass[perm_inds], weight=evt_weights[perm_inds])
            acc["mtbar"].fill(   dataset=dataset_name, rewt=rewt_type, evt_type="RECO", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtop=reco_tbar_mass[perm_inds], weight=evt_weights[perm_inds])
            acc["pt_top"].fill(  dataset=dataset_name, rewt=rewt_type, evt_type="RECO", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=reco_top_pt[perm_inds], weight=evt_weights[perm_inds])
            acc["pt_tbar"].fill( dataset=dataset_name, rewt=rewt_type, evt_type="RECO", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=reco_tbar_pt[perm_inds], weight=evt_weights[perm_inds])
            acc["pt_tt"].fill(   dataset=dataset_name, rewt=rewt_type, evt_type="RECO", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=ak.flatten(perm["TTbar"].pt)[perm_inds], weight=evt_weights[perm_inds])
            acc["eta_top"].fill( dataset=dataset_name, rewt=rewt_type, evt_type="RECO", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=reco_top_eta[perm_inds], weight=evt_weights[perm_inds])
            acc["eta_tbar"].fill(dataset=dataset_name, rewt=rewt_type, evt_type="RECO", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=reco_tbar_eta[perm_inds], weight=evt_weights[perm_inds])
            acc["eta_tt"].fill(  dataset=dataset_name, rewt=rewt_type, evt_type="RECO", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=ak.flatten(perm["TTbar"].eta)[perm_inds], weight=evt_weights[perm_inds])

            acc["top_ctstar"].fill(     dataset=dataset_name, rewt=rewt_type, evt_type="RECO",  jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, ctstar=reco_top_ctstar[perm_inds], weight=evt_weights[perm_inds])
            acc["tbar_ctstar"].fill(    dataset=dataset_name, rewt=rewt_type, evt_type="RECO",  jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, ctstar=reco_tbar_ctstar[perm_inds], weight=evt_weights[perm_inds])
            acc["top_ctstar_abs"].fill( dataset=dataset_name, rewt=rewt_type, evt_type="RECO",  jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, ctstar_abs=np.abs(reco_top_ctstar[perm_inds]), weight=evt_weights[perm_inds])
            acc["tbar_ctstar_abs"].fill(dataset=dataset_name, rewt=rewt_type, evt_type="RECO",  jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, ctstar_abs=np.abs(reco_tbar_ctstar[perm_inds]), weight=evt_weights[perm_inds])

            acc["mtt_vs_top_ctstar_abs"].fill( dataset=dataset_name, rewt=rewt_type, evt_type="RECO",  jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion,
                mtt=ak.flatten(perm["TTbar"].mass)[perm_inds], ctstar_abs=np.abs(reco_top_ctstar[perm_inds]), weight=evt_weights[perm_inds])
            acc["mtt_vs_tbar_ctstar_abs"].fill(dataset=dataset_name, rewt=rewt_type, evt_type="RECO",  jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion,
                mtt=ak.flatten(perm["TTbar"].mass)[perm_inds], ctstar_abs=np.abs(reco_tbar_ctstar[perm_inds]), weight=evt_weights[perm_inds])

            acc["mtt_vs_nnlo_weights"].fill(dataset=dataset_name, rewt=rewt_type, evt_type="RECO", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtt=ak.flatten(perm["TTbar"].mass)[perm_inds], nnlo_wts=nnlo_weights[perm_inds])

                ## fill gen-level quantities
            acc["mtt"].fill(     dataset=dataset_name, rewt=rewt_type, evt_type="GEN", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtt=ak.flatten(gentt["TTbar"].mass, axis=None)[perm_inds], weight=evt_weights[perm_inds])
            acc["mtop"].fill(    dataset=dataset_name, rewt=rewt_type, evt_type="GEN", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtop=ak.flatten(gentt["Top"].mass, axis=None)[perm_inds], weight=evt_weights[perm_inds])
            acc["mtbar"].fill(   dataset=dataset_name, rewt=rewt_type, evt_type="GEN", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtop=ak.flatten(gentt["Tbar"].mass, axis=None)[perm_inds], weight=evt_weights[perm_inds])
            acc["pt_top"].fill(  dataset=dataset_name, rewt=rewt_type, evt_type="GEN", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=ak.flatten(gentt["Top"].pt, axis=None)[perm_inds], weight=evt_weights[perm_inds])
            acc["pt_tbar"].fill( dataset=dataset_name, rewt=rewt_type, evt_type="GEN", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=ak.flatten(gentt["Tbar"].pt, axis=None)[perm_inds], weight=evt_weights[perm_inds])
            acc["pt_tt"].fill(   dataset=dataset_name, rewt=rewt_type, evt_type="GEN", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, pt=ak.flatten(gentt["TTbar"].pt, axis=None)[perm_inds], weight=evt_weights[perm_inds])
            acc["eta_top"].fill( dataset=dataset_name, rewt=rewt_type, evt_type="GEN", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=ak.flatten(gentt["Top"].eta, axis=None)[perm_inds], weight=evt_weights[perm_inds])
            acc["eta_tbar"].fill(dataset=dataset_name, rewt=rewt_type, evt_type="GEN", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=ak.flatten(gentt["Tbar"].eta, axis=None)[perm_inds], weight=evt_weights[perm_inds])
            acc["eta_tt"].fill(  dataset=dataset_name, rewt=rewt_type, evt_type="GEN", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, eta=ak.flatten(gentt["TTbar"].eta, axis=None)[perm_inds], weight=evt_weights[perm_inds])

            acc["top_ctstar"].fill(     dataset=dataset_name, rewt=rewt_type, evt_type="GEN", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, ctstar=gen_top_ctstar[perm_inds], weight=evt_weights[perm_inds])
            acc["tbar_ctstar"].fill(    dataset=dataset_name, rewt=rewt_type, evt_type="GEN", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, ctstar=gen_tbar_ctstar[perm_inds], weight=evt_weights[perm_inds])
            acc["top_ctstar_abs"].fill( dataset=dataset_name, rewt=rewt_type, evt_type="GEN", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, ctstar_abs=np.abs(gen_top_ctstar[perm_inds]), weight=evt_weights[perm_inds])
            acc["tbar_ctstar_abs"].fill(dataset=dataset_name, rewt=rewt_type, evt_type="GEN", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, ctstar_abs=np.abs(gen_tbar_ctstar[perm_inds]), weight=evt_weights[perm_inds])

            acc["mtt_vs_top_ctstar_abs"].fill(dataset=dataset_name, rewt=rewt_type, evt_type="GEN", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion,
                mtt=ak.flatten(gentt["TTbar"].mass, axis=None)[perm_inds], ctstar_abs=np.abs(gen_top_ctstar[perm_inds]), weight=evt_weights[perm_inds])
            acc["mtt_vs_tbar_ctstar_abs"].fill(dataset=dataset_name, rewt=rewt_type, evt_type="GEN", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion,
                mtt=ak.flatten(gentt["TTbar"].mass, axis=None)[perm_inds], ctstar_abs=np.abs(gen_tbar_ctstar[perm_inds]), weight=evt_weights[perm_inds])

            acc["mtt_vs_nnlo_weights"].fill(dataset=dataset_name, rewt=rewt_type, evt_type="GEN", jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, mtt=ak.flatten(gentt["TTbar"].mass, axis=None)[perm_inds], nnlo_wts=nnlo_weights[perm_inds])

                ## fill gen-reco quantities
            acc["Reso_mtt"].fill(dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion,
                reso_mtt=ak.flatten(gentt["TTbar"].mass, axis=None)[perm_inds]-ak.flatten(perm["TTbar"].mass)[perm_inds], weight=evt_weights[perm_inds])
            acc["Reso_mtop"].fill(dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion,
                reso_mtop=ak.flatten(gentt["Top"].mass, axis=None)[perm_inds]-reco_top_mass[perm_inds], weight=evt_weights[perm_inds])
            acc["Reso_mtbar"].fill(dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion,
                reso_mtop=ak.flatten(gentt["Tbar"].mass, axis=None)[perm_inds]-reco_tbar_mass[perm_inds], weight=evt_weights[perm_inds])
            acc["Reso_pt_top"].fill(dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion,
                reso_pt=ak.flatten(gentt["Top"].pt, axis=None)[perm_inds]-reco_top_pt[perm_inds], weight=evt_weights[perm_inds])
            acc["Reso_pt_tbar"].fill(dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion,
                reso_pt=ak.flatten(gentt["Tbar"].pt, axis=None)[perm_inds]-reco_tbar_pt[perm_inds], weight=evt_weights[perm_inds])
            acc["Reso_pt_tt"].fill(dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion,
                reso_pt=ak.flatten(gentt["TTbar"].pt, axis=None)[perm_inds]-ak.flatten(perm["TTbar"].pt)[perm_inds], weight=evt_weights[perm_inds])
            acc["Reso_eta_top"].fill(dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion,
                reso_eta=ak.flatten(gentt["Top"].eta, axis=None)[perm_inds]-reco_top_eta[perm_inds], weight=evt_weights[perm_inds])
            acc["Reso_eta_tbar"].fill(dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion,
                reso_eta=ak.flatten(gentt["Tbar"].eta, axis=None)[perm_inds]-reco_tbar_eta[perm_inds], weight=evt_weights[perm_inds])
            acc["Reso_eta_tt"].fill(dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion,
                reso_eta=ak.flatten(gentt["TTbar"].eta, axis=None)[perm_inds]-ak.flatten(perm["TTbar"].eta)[perm_inds], weight=evt_weights[perm_inds])

            acc["Reso_top_ctstar"].fill(     dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, reso_ctstar=gen_top_ctstar[perm_inds]-reco_top_ctstar[perm_inds], weight=evt_weights[perm_inds])
            acc["Reso_tbar_ctstar"].fill(    dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, reso_ctstar=gen_tbar_ctstar[perm_inds]-reco_tbar_ctstar[perm_inds], weight=evt_weights[perm_inds])
            acc["Reso_top_ctstar_abs"].fill( dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, reso_ctstar_abs=np.abs(gen_top_ctstar[perm_inds])-np.abs(reco_top_ctstar[perm_inds]), weight=evt_weights[perm_inds])
            acc["Reso_tbar_ctstar_abs"].fill(dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion, reso_ctstar_abs=np.abs(gen_tbar_ctstar[perm_inds])-np.abs(reco_tbar_ctstar[perm_inds]), weight=evt_weights[perm_inds])

                ## fill gen-reco quantities

            acc["Reso_mtt_vs_top_ctstar_abs"].fill( dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion,
                reso_mtt=ak.flatten(gentt["TTbar"].mass, axis=None)[perm_inds]-ak.flatten(perm["TTbar"].mass)[perm_inds], reso_ctstar_abs=np.abs(gen_top_ctstar[perm_inds])-np.abs(reco_top_ctstar[perm_inds]), weight=evt_weights[perm_inds])
            acc["Reso_mtt_vs_tbar_ctstar_abs"].fill(dataset=dataset_name, rewt=rewt_type, jmult=jetmult, leptype=leptype, lepcat=lepcat, btag=btagregion,
                reso_mtt=ak.flatten(gentt["TTbar"].mass, axis=None)[perm_inds]-ak.flatten(perm["TTbar"].mass)[perm_inds], reso_ctstar_abs=np.abs(gen_tbar_ctstar[perm_inds])-np.abs(reco_tbar_ctstar[perm_inds]), weight=evt_weights[perm_inds])

        return acc        


    def postprocess(self, accumulator):
        return accumulator


proc_executor = processor.iterative_executor if to_debug else processor.futures_executor
#proc_exec_args = {"schema": NanoAODSchema} if to_debug else {"schema": NanoAODSchema, "workers": 8}
proc_exec_args = {"schema": processor.NanoAODSchema} if to_debug else {"schema": processor.NanoAODSchema, "workers": 8}
output = processor.run_uproot_job(
    fileset,
    treename="Events",
    processor_instance=apply_nnlo_weights(),
    executor=proc_executor,
    executor_args=proc_exec_args,
    chunksize=100000,
)

save(output, args.outfname)
print("%s has been written" % args.outfname)

toc = time.time()
print("Total time: %.1f" % (toc - tic))
