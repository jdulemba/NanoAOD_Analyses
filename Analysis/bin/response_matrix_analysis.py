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
import python.MCWeights as MCWeights
import numpy as np
import Utilities.prettyjson as prettyjson
import Utilities.make_variables as make_vars
import python.GenParticleSelector as genpsel
import python.TTPermutator as ttpermutator
import Utilities.systematics as systematics
import python.IDJet as IDJet
from copy import deepcopy
import Utilities.final_analysis_binning as final_binning

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

isTTSL_ = samplename == "ttJetsSL"
isSignal_ = (samplename.startswith("AtoTT") or samplename.startswith("HtoTT"))
isInt_ = isSignal_ and ("Int" in samplename)

if not (isTTSL_ or isSignal_):
    raise ValueError("This analyzer should only be run with tt l+jets or signal events")

# get dataset classification, used for corrections/systematics
isTTbar_ = samplename.startswith("ttJets")

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
#mtopcuts = {
#    "3Jets" : (jet_pars["mtop"]["3Jetsmin"], jet_pars["mtop"]["3Jetsmax"]),
#    "4PJets" : (jet_pars["mtop"]["4PJetsmin"], jet_pars["mtop"]["4PJetsmax"]),
#}

# get systematics to run
evt_sys_to_run = opts_dict.get("evt_sys", "NONE").upper()
rewt_sys_to_run = opts_dict.get("rewt_sys", "NONE").upper()
only_sys = ast.literal_eval(opts_dict.get("only_sys", "False"))

import fnmatch
if only_sys: # don't run 'nosys'
    if (evt_sys_to_run == "NONE") and (rewt_sys_to_run == "NONE"):
        raise ValueError("At least one systematic must be specified in order to run only on systematics!")

    event_systematics_to_run = ["nosys"] if (rewt_sys_to_run != "NONE") else []
    reweight_systematics_to_run = []

else:
    event_systematics_to_run = ["nosys"]
    reweight_systematics_to_run = ["nosys"]

event_systematics_to_run += [systematics.event_sys_opts[args.year][name] for name in systematics.event_sys_opts[args.year].keys() if fnmatch.fnmatch(name, evt_sys_to_run)]
if isSignal_:
    reweight_systematics_to_run += [systematics.signal_reweight_opts[args.year][name] for name in systematics.signal_reweight_opts[args.year].keys() if fnmatch.fnmatch(name, rewt_sys_to_run)]
else:
    reweight_systematics_to_run += [systematics.reweight_sys_opts[args.year][name] for name in systematics.reweight_sys_opts[args.year].keys() if fnmatch.fnmatch(name, rewt_sys_to_run)]

    ## check that systematics only related to ttbar events aren't used for non-ttbar events
if not isTTbar_:
    reweight_systematics_to_run = [sys for sys in reweight_systematics_to_run if sys not in systematics.ttJets_sys.values()]

print("\n\nRunning with event systematics:", *sorted(set(event_systematics_to_run).difference(set(["nosys"]))), sep=", ") if "nosys" in event_systematics_to_run else print("Running with event systematics:", *sorted(event_systematics_to_run), sep=", ")
print("\tand reweight systematics:", *sorted(reweight_systematics_to_run))
#set_trace()


# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class htt_btag_sb_regions(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.sys_axis = hist.Cat("sys", "Systematic")
        self.jetmult_axis = hist.Cat("jmult", "nJets")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")

        self.reco_st_mtt_ctstar_axis = hist.Bin("reco_st_mtt_ctstar_idx", "Reco st mtt ctstar",
            np.around(np.arange((final_binning.st_binning.size-1) * (final_binning.ctstar_abs_binning.size-1) * (final_binning.mtt_binning.size-1) + 1), decimals=0))
        self.gen_st_mtt_ctstar_axis = hist.Bin("gen_st_mtt_ctstar_idx", "Gen st mtt ctstar",
            np.around(np.arange((final_binning.st_binning.size-1) * (final_binning.ctstar_abs_binning.size-1) * (final_binning.mtt_binning.size-1) + 1), decimals=0))

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
            histo_dict["cutflow_%s" % sys] = processor.defaultdict_accumulator(int)
    
        self._accumulator = processor.dict_accumulator(histo_dict)

            ## define regions depending on data or MC, systematic or nosys
        base_regions = {
            "Muon" : {
                "3Jets"  : {"lep_and_filter_pass", "passing_jets", "jets_3" , "tight_MU", "DeepCSV_pass"},
                "4PJets" : {"lep_and_filter_pass", "passing_jets", "jets_4p", "tight_MU", "DeepCSV_pass"},
            },
            "Electron" : {
                "3Jets"  : {"lep_and_filter_pass", "passing_jets", "jets_3" , "tight_EL", "DeepCSV_pass"},
                "4PJets" : {"lep_and_filter_pass", "passing_jets", "jets_4p", "tight_EL", "DeepCSV_pass"},
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
        histo_dict["Reco_mtt_x_tlep_ctstar_abs_x_st_inds"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.reco_st_mtt_ctstar_axis)
        if isTTSL_:
            histo_dict["Reco_VS_Gen_mtt_x_tlep_ctstar_abs_x_st_inds"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis,
                self.gen_st_mtt_ctstar_axis, self.reco_st_mtt_ctstar_axis)

        return histo_dict


    def process(self, events):
        output = self.accumulator.identity()

        self.sample_name = events.metadata["dataset"]

            ## make event weights
                # data or MC distinction made internally
        mu_evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections, isTTbar=isTTSL_, isSignal=isSignal_)
        el_evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections, isTTbar=isTTSL_, isSignal=isSignal_)

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

            # find gen level particles for ttbar system and other ttbar corrections
        if isTTSL_:
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
                        np.copy(ewk_wts_dict["Rebinned_KFactor_1.1"]),
                        np.copy(ewk_wts_dict["Rebinned_KFactor_0.9"]),
                    )
                    el_evt_weights.add("Yukawa",  # really just varying value of Yt
                        np.copy(ewk_wts_dict["Rebinned_KFactor_1.0"]),
                        np.copy(ewk_wts_dict["Rebinned_KFactor_1.1"]),
                        np.copy(ewk_wts_dict["Rebinned_KFactor_0.9"]),
                    )

            genpsel.select(events, mode="NORMAL")
            {selection[sys].add("semilep", ak.num(events["SL"]) > 0) for sys in selection.keys()}


            # run over systematics that require changes to event objects (jets+MET)
        for evt_sys in self.event_systematics_to_run:
            output[f"cutflow_{evt_sys}"]["lep_and_filter_pass"] += ak.sum(lep_and_filter_pass)
                # jet selection
            passing_jets = objsel.jets_selection(events, year=args.year, cutflow=output[f"cutflow_{evt_sys}"], shift=evt_sys, hem_15_16=apply_hem)
            output[f"cutflow_{evt_sys}"]["nEvts passing jet and lepton obj selection"] += ak.sum(passing_jets & lep_and_filter_pass)
            selection[evt_sys].add("passing_jets", passing_jets)
            selection[evt_sys].add("jets_3",  ak.num(events["SelectedJets"]) == 3)
            selection[evt_sys].add("jets_4p",  ak.num(events["SelectedJets"]) > 3) # only for getting btag weights
            selection[evt_sys].add("DeepCSV_pass", ak.sum(events["SelectedJets"][btag_wps[0]], axis=1) >= 2)

                # sort jets by btag value
            events["SelectedJets"] = events["SelectedJets"][ak.argsort(events["SelectedJets"]["btagDeepB"], ascending=False)] if btaggers[0] == "DeepCSV" else events["SelectedJets"][ak.argsort(events["SelectedJets"]["btagDeepFlavB"], ascending=False)]

                ## apply btagging SFs to MC
            if corrections["BTagSF"] == True:
                btag_weights = {key : np.ones(len(events)) for key in self.corrections["BTag_Constructors"]["DeepCSV"]["3Jets"].schema_.keys()}

                threeJets_cut = selection[evt_sys].require(lep_and_filter_pass=True, passing_jets=True, jets_3=True)
                deepcsv_3j_wts = self.corrections["BTag_Constructors"]["DeepCSV"]["3Jets"].get_scale_factor(jets=events["SelectedJets"][threeJets_cut], passing_cut="DeepCSV"+wps_to_use[0])
                fourplusJets_cut = selection[evt_sys].require(lep_and_filter_pass=True, passing_jets=True, jets_4p=True)
                deepcsv_4pj_wts = self.corrections["BTag_Constructors"]["DeepCSV"]["4PJets"].get_scale_factor(jets=events["SelectedJets"][fourplusJets_cut], passing_cut="DeepCSV"+wps_to_use[0])

                    # fll dict of btag weights
                for wt_name in deepcsv_3j_wts.keys():
                    btag_weights[wt_name][threeJets_cut] = ak.prod(deepcsv_3j_wts[wt_name], axis=1)
                    btag_weights[wt_name][fourplusJets_cut] = ak.prod(deepcsv_4pj_wts[wt_name], axis=1)
            else:
                raise ValueError("BTag SFs not applied to MC")


            ## fill hists for each region
            for lepton in self.regions[evt_sys].keys():
                evt_weights = mu_evt_weights if lepton == "Muon" else el_evt_weights
                for jmult in self.regions[evt_sys][lepton].keys():
                    #set_trace()
                    cut = selection[evt_sys].all(*self.regions[evt_sys][lepton][jmult])

                    output[f"cutflow_{evt_sys}"]["nEvts %s" % ", ".join([lepton, jmult])] += cut.sum()

                    if to_debug: print(lepton, jmult)
                    if cut.sum() > 0:
                        ltype = "MU" if lepton == "Muon" else "EL"
                        if f"loose_or_tight_{ltype}" in self.regions[evt_sys][lepton][jmult]:
                            leptons = events[lepton][cut][((events[lepton][cut][f"TIGHT{ltype}"] == True) | (events[lepton][cut][f"LOOSE{ltype}"] == True))]
                        elif f"tight_{ltype}" in self.regions[evt_sys][lepton][jmult]:
                            leptons = events[lepton][cut][(events[lepton][cut][f"TIGHT{ltype}"] == True)]
                        elif f"loose_{ltype}" in self.regions[evt_sys][lepton][jmult]:
                            leptons = events[lepton][cut][(events[lepton][cut][f"LOOSE{ltype}"] == True)]
                        else:
                            raise ValueError("Not sure what lepton type to choose for event")

                            # get jets and MET
                        jets, met = events["SelectedJets"][cut], events["SelectedMET"][cut]

                            # find best permutations
                        best_perms = ttpermutator.find_best_permutations(jets=jets, leptons=leptons, MET=met, btagWP=btag_wps[0], btag_req=True)
                        valid_perms = ak.num(best_perms["TTbar"].pt) > 0
                        output[f"cutflow_{evt_sys}"]["nEvts %s: valid perms" % ", ".join([lepton, jmult])] += ak.sum(valid_perms)

                            ## create MT regions
                        MT = make_vars.MT(leptons, met)
                        MTHigh = ak.flatten(MT[valid_perms] >= MTcut)
                        output[f"cutflow_{evt_sys}"]["nEvts %s: pass MT cut" % ", ".join([lepton, jmult])] += ak.sum(MTHigh)

                            # fill hists for each systematic
                        if to_debug: print("  evt sys:", evt_sys)
                        if evt_sys == "nosys":
                            for rewt_sys in self.reweight_systematics_to_run:
                                #if to_debug: set_trace()
                                if to_debug: print("\tsysname:", rewt_sys)

                                if rewt_sys == "nosys":
                                    wts = (evt_weights.weight()*btag_weights["central"])[cut][valid_perms][MTHigh]
                                elif rewt_sys.startswith("btag"):
                                    wts = (evt_weights.weight()*btag_weights[rewt_sys.split("btag_")[-1]])[cut][valid_perms][MTHigh]
                                else:
                                    if rewt_sys not in evt_weights.variations:
                                    #if rewt_sys not in evt_weights._modifiers.keys():
                                        print(f"{rewt_sys} not option in event weights. Skipping")
                                        continue
                                    wts = (evt_weights.weight(rewt_sys)*btag_weights["central"])[cut][valid_perms][MTHigh]

                                #set_trace()
                                if isInt_:
                                        # fill hists for positive weights
                                    pos_evts = wts > 0
                                    self.sample_name = "%s_pos" % events.metadata["dataset"]
                                    output = self.fill_hists(acc=output, sys=rewt_sys, jetmult=jmult, leptype=lepton,
                                        reco_perm=best_perms[valid_perms][MTHigh][pos_evts], gentt=None, evt_wts=wts[pos_evts])
                                        # fill hists for negative weights
                                    neg_evts = wts < 0
                                    self.sample_name = "%s_neg" % events.metadata["dataset"]
                                    output = self.fill_hists(acc=output, sys=rewt_sys, jetmult=jmult, leptype=lepton,
                                        reco_perm=best_perms[valid_perms][MTHigh][neg_evts], gentt=None, evt_wts=wts[neg_evts])
                                else:
                                    output = self.fill_hists(acc=output, sys=rewt_sys, jetmult=jmult, leptype=lepton,
                                        reco_perm=best_perms[valid_perms][MTHigh], gentt=events["SL"][cut][valid_perms][MTHigh] if isTTSL_ else None, evt_wts=wts)

                        else:
                            if to_debug: print("\tsysname:", evt_sys)
                            wts = (evt_weights.weight()*btag_weights["central"])[cut][valid_perms][MTHigh]
                            if isInt_:
                                    # fill hists for positive weights
                                pos_evts = wts > 0
                                self.sample_name = "%s_pos" % events.metadata["dataset"]
                                output = self.fill_hists(acc=output, sys=evt_sys, jetmult=jmult, leptype=lepton,
                                    reco_perm=best_perms[valid_perms][MTHigh][pos_evts], gentt=None, evt_wts=wts[pos_evts])
                                    # fill hists for negative weights
                                neg_evts = wts < 0
                                self.sample_name = "%s_neg" % events.metadata["dataset"]
                                output = self.fill_hists(acc=output, sys=evt_sys, jetmult=jmult, leptype=lepton,
                                    reco_perm=best_perms[valid_perms][MTHigh][neg_evts], gentt=None, evt_wts=wts[neg_evts])
                            else:
                                output = self.fill_hists(acc=output, sys=evt_sys, jetmult=jmult, leptype=lepton,
                                    reco_perm=best_perms[valid_perms][MTHigh], gentt=events["SL"][cut][valid_perms][MTHigh] if isTTSL_ else None, evt_wts=wts)

        return output

    def fill_hists(self, acc, sys, jetmult, leptype, reco_perm, gentt, evt_wts):
            ## apply alpha correction for 3Jets
        if (jetmult == "3Jets") and ("Alpha" in self.corrections):
            alpha_corr = self.corrections["Alpha"](172.5/reco_perm["THad"].mass)
            reco_perm["THad"] = reco_perm["THad"].multiply(alpha_corr) # correct thad
            reco_perm["TTbar"] = ak.flatten(reco_perm["THad"]+reco_perm["TLep"]) # correct ttbar

        #    # keep events within top mass window
        #mtop_mask = (mtopcuts[jetmult][0] < reco_perm["THad"].mass) & (mtopcuts[jetmult][1] > reco_perm["THad"].mass)

            # find reco values
        reco_thad_ctstar, reco_tlep_ctstar = make_vars.ctstar(reco_perm["THad"], reco_perm["TLep"])
        reco_thad_ctstar, reco_tlep_ctstar = ak.flatten(reco_thad_ctstar, axis=None), ak.flatten(reco_tlep_ctstar, axis=None)
        reco_tlep_ctstar_abs = np.abs(reco_tlep_ctstar)
        reco_ST = ak.flatten(reco_perm["THad"].pt + reco_perm["TLep"].pt, axis=None)
        reco_mtt = ak.flatten(reco_perm["TTbar"].mass, axis=None)

        # Reco inds and masks
            # mtt
        reco_mtt_ind_vals = np.array([np.argmax(reco_mtt[idx] < final_binning.mtt_binning) for idx in range(len(reco_mtt))])
        reco_mtt_inds_mask = (reco_mtt > final_binning.mtt_binning[0]) & (reco_mtt <= final_binning.mtt_binning[-1])
            # ctstar
        reco_ctstar_ind_vals = np.array([np.argmax(reco_tlep_ctstar_abs[idx] < final_binning.ctstar_abs_binning) for idx in range(len(reco_tlep_ctstar_abs))])
        reco_ctstar_inds_mask = (reco_tlep_ctstar_abs > final_binning.ctstar_abs_binning[0]) & (reco_tlep_ctstar_abs <= final_binning.ctstar_abs_binning[-1])
            # ST
        reco_st_ind_vals = np.array([np.argmax(reco_ST[idx] < final_binning.st_binning) for idx in range(len(reco_ST))])
        reco_st_inds_mask = (reco_ST > final_binning.st_binning[0]) & (reco_ST <= final_binning.st_binning[-1])
            # ST x mtt x ctstar
        reco_st_mtt_ctstar_ind_vals = (reco_st_ind_vals - 1) * (final_binning.mtt_binning.size - 1) * (final_binning.ctstar_abs_binning.size - 1) + (reco_ctstar_ind_vals - 1) * (final_binning.mtt_binning.size - 1) + reco_mtt_ind_vals
        reco_st_mtt_ctstar_inds_mask = reco_mtt_inds_mask & reco_ctstar_inds_mask & reco_st_inds_mask #& ak.flatten(mtop_mask, axis=None)

            # fill 1D reco inds hist 
        acc["Reco_mtt_x_tlep_ctstar_abs_x_st_inds"].fill(dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype,
            reco_st_mtt_ctstar_idx=reco_st_mtt_ctstar_ind_vals[reco_st_mtt_ctstar_inds_mask]-0.5, weight=evt_wts[reco_st_mtt_ctstar_inds_mask])

        if isTTSL_:
            #set_trace()
            #    # get kinematic cuts on charged lepton and both b quarks
            #gen_kin_cuts =  ak.flatten((gentt["Lepton"]["pt"] >= 25.) & (np.abs(gentt["Lepton"]["eta"]) <= 2.5) & \
            #            (gentt["BHad"]["pt"] >= 25.) & (np.abs(gentt["BHad"]["eta"]) <= 2.5) & \
            #            (gentt["BLep"]["pt"] >= 25.) & (np.abs(gentt["BLep"]["eta"]) <= 2.5), axis=1)

                # find gen values
            gen_thad_ctstar, gen_tlep_ctstar = make_vars.ctstar(gentt["THad"], gentt["TLep"])
            gen_thad_ctstar, gen_tlep_ctstar = ak.flatten(gen_thad_ctstar, axis=None), ak.flatten(gen_tlep_ctstar, axis=None)
            gen_tlep_ctstar_abs = np.abs(gen_tlep_ctstar)
            gen_ST = ak.flatten(gentt["THad"].pt + gentt["TLep"].pt, axis=None)
            gen_mtt = ak.flatten(gentt["TTbar"].mass, axis=None)

            #set_trace()
            # Gen inds and masks
                # mtt
            gen_mtt_ind_vals = np.array([np.argmax(gen_mtt[idx] < final_binning.mtt_binning) for idx in range(len(gen_mtt))])
            gen_mtt_inds_mask = (gen_mtt > final_binning.mtt_binning[0]) & (gen_mtt <= final_binning.mtt_binning[-1])
                # ctstar
            gen_ctstar_ind_vals = np.array([np.argmax(gen_tlep_ctstar_abs[idx] < final_binning.ctstar_abs_binning) for idx in range(len(gen_tlep_ctstar_abs))])
            gen_ctstar_inds_mask = (gen_tlep_ctstar_abs > final_binning.ctstar_abs_binning[0]) & (gen_tlep_ctstar_abs <= final_binning.ctstar_abs_binning[-1])
                # ST
            gen_st_ind_vals = np.array([np.argmax(gen_ST[idx] < final_binning.st_binning) for idx in range(len(gen_ST))])
            gen_st_inds_mask = (gen_ST > final_binning.st_binning[0]) & (gen_ST <= final_binning.st_binning[-1])
                # ST x mtt x ctstar
            gen_st_mtt_ctstar_ind_vals = (gen_st_ind_vals - 1) * (final_binning.mtt_binning.size - 1) * (final_binning.ctstar_abs_binning.size - 1) + (gen_ctstar_ind_vals - 1) * (final_binning.mtt_binning.size - 1) + gen_mtt_ind_vals
            #gen_st_mtt_ctstar_inds_mask = gen_mtt_inds_mask & gen_ctstar_inds_mask & gen_st_inds_mask

                    # with kinematic cut
            gen_st_mtt_ctstar_inds_mask = gen_mtt_inds_mask & gen_ctstar_inds_mask & gen_st_inds_mask #& gen_kin_cuts

                # fill 2D gen (x axis) vs reco (y axis) inds hist
            gen_reco_inds_mask = reco_st_mtt_ctstar_inds_mask & gen_st_mtt_ctstar_inds_mask
            acc["Reco_VS_Gen_mtt_x_tlep_ctstar_abs_x_st_inds"].fill(dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype,
                gen_st_mtt_ctstar_idx=gen_st_mtt_ctstar_ind_vals[gen_reco_inds_mask]-0.5,
                reco_st_mtt_ctstar_idx=reco_st_mtt_ctstar_ind_vals[gen_reco_inds_mask]-0.5, weight=evt_wts[gen_reco_inds_mask])

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
