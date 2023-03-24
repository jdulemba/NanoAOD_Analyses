#!/usr/bin/env python

import time
tic = time.time()

from coffea import hist, processor
from coffea.nanoevents import NanoAODSchema
from coffea.analysis_tools import PackedSelection
import awkward as ak
from coffea.nanoevents.methods import vector
ak.behavior.update(vector.behavior)

from coffea.util import save, load
from pdb import set_trace
import os, sys
import python.ObjectSelection as objsel
import Utilities.plot_tools as plt_tools
import python.MCWeights as MCWeights
import Utilities.make_variables as make_vars
import numpy as np
import Utilities.prettyjson as prettyjson
import coffea.lumi_tools.lumi_tools as lumi_tools
import python.GenParticleSelector as genpsel
import python.TTGenMatcher as ttmatcher
import python.IDJet as IDJet

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
base_jobid = os.environ["base_jobid"]
analyzer = "ttbar_evt_cats"

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

    ## specify ttJets samples
allowed_samples = ["ttJets_PS"] if ((args.year == "2016") and (base_jobid == "NanoAODv6")) else ["ttJetsSL"]
isAllowed = np.array([(key in allowed_samples) for key in fileset.keys()]).all()
if not isAllowed:
    raise ValueError("Not a valid dataset to run on")

# convert input string of options dictionary to actual dictionary
odict = (args.opts).replace("\'", "\"")
opts_dict = prettyjson.loads(odict)

    ## set config options passed through argparse
import ast
to_debug = ast.literal_eval(opts_dict.get("debug", "False"))
apply_hem = ast.literal_eval(opts_dict.get("apply_hem", "True"))


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
btaggers = [jet_pars["btagger"]]

wps_to_use = list(set([jet_pars["permutations"]["tightb"],jet_pars["permutations"]["looseb"]]))
if not( len(wps_to_use) == 1):
    raise IOError("Only 1 unique btag working point supported now")
btag_wp = btaggers[0]+wps_to_use[0]
#btag_wps = ["DeepCSVMedium"]

if corrections["BTagSF"] == True:
    sf_file = os.path.join(proj_dir, "Corrections", jobid, jet_pars["btagging"]["btagSF_file"])
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


# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class ttbar_evt_cats(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")
        self.jetmult_axis = hist.Cat("jmult", "nJets")
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

        histo_dict["cutflow"] = processor.defaultdict_accumulator(int)

        self._accumulator = processor.dict_accumulator(histo_dict)
        self.sample_name = ""
        self.corrections = corrections

        self.regions = {
            "Muon" : {
                "Tight" : {
                    "3Jets"  : {"lep_and_filter_pass", "passing_jets", "jets_3" , "tight_MU", f"{btaggers[0]}_pass", "semilep"},
                    "4PJets"  : {"lep_and_filter_pass", "passing_jets", "jets_4p" , "tight_MU", f"{btaggers[0]}_pass", "semilep"},
                    #"3Jets"  : {"lep_and_filter_pass", "passing_jets", "jets_3" , "tight_MU", "DeepCSV_pass", "semilep"},
                    #"4PJets"  : {"lep_and_filter_pass", "passing_jets", "jets_4p" , "tight_MU", "DeepCSV_pass", "semilep"},
                },
            },
            "Electron" : {
                "Tight" : {
                    "3Jets"  : {"lep_and_filter_pass", "passing_jets", "jets_3" , "tight_EL", f"{btaggers[0]}_pass", "semilep"},
                    "4PJets"  : {"lep_and_filter_pass", "passing_jets", "jets_4p" , "tight_EL", f"{btaggers[0]}_pass", "semilep"},
                    #"3Jets"  : {"lep_and_filter_pass", "passing_jets", "jets_3" , "tight_EL", "DeepCSV_pass", "semilep"},
                    #"4PJets"  : {"lep_and_filter_pass", "passing_jets", "jets_4p" , "tight_EL", "DeepCSV_pass", "semilep"},
                },
            },
        }
    
    @property
    def accumulator(self):
        return self._accumulator


    def make_hists(self):
        histo_dict = {}
        histo_dict["mtt"]      = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.cat_axis, self.mtt_axis)
        histo_dict["pass_mtt"] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.cat_axis, self.object_axis, self.mtt_axis)
        histo_dict["fail_mtt"] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.cat_axis, self.object_axis, self.mtt_axis)

        return histo_dict



    def process(self, events):
        np.random.seed(10) # sets seed so values from random distributions are reproducible (JER corrections)
        output = self.accumulator.identity()

        self.sample_name = events.metadata["dataset"]

            ## make event weights
                # data or MC distinction made internally
        mu_evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections, isTTbar=True)
        el_evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections, isTTbar=True)

            ## initialize selections and regions
        selection = PackedSelection()

            # get all passing leptons
        lep_and_filter_pass = objsel.select_leptons(events, year=args.year, hem_15_16=apply_hem)
        output["cutflow"]["lep_and_filter_pass"] += ak.sum(lep_and_filter_pass)
        selection.add("lep_and_filter_pass", lep_and_filter_pass)

            ## build corrected jets and MET
        events["Jet"], events["MET"] = IDJet.process_jets(events, args.year, self.corrections["JetCor"])

        ## add different selections
                ## muons
        tight_mu_sel = ak.sum(events["Muon"]["TIGHTMU"], axis=1) == 1
        selection.add("tight_MU", tight_mu_sel) # one muon passing TIGHT criteria
                ## electrons
        tight_el_sel = ak.sum(events["Electron"]["TIGHTEL"], axis=1) == 1
        selection.add("tight_EL", tight_el_sel) # one electron passing TIGHT criteria

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

        ## find gen level particles for ttbar system
        genpsel.select(events, mode="NORMAL")
        selection.add("semilep", ak.num(events["SL"]) > 0)

            # jet selection
        passing_jets = objsel.jets_selection(events, year=args.year, cutflow=output["cutflow"], hem_15_16=apply_hem)
        output["cutflow"]["nEvts passing jet and lepton obj selection"] += ak.sum(passing_jets & lep_and_filter_pass)
        selection.add("passing_jets", passing_jets)
        selection.add("jets_3",  ak.num(events["SelectedJets"]) == 3)
        selection.add("jets_4p",  ak.num(events["SelectedJets"]) > 3) # only for getting btag weights
        selection.add(f"{btaggers[0]}_pass", ak.sum(events["SelectedJets"][btag_wp], axis=1) >= 2)

            # sort jets by btag value
        events["SelectedJets"] = events["SelectedJets"][ak.argsort(events["SelectedJets"][IDJet.btag_tagger_to_disc_name[btaggers[0]]], ascending=False)]

            ## apply btagging SFs to MC
        if (corrections["BTagSF"] == True):
            btag_weights = {key : np.ones(len(events)) for key in self.corrections["BTag_Constructors"][btaggers[0]]["3Jets"].schema_.keys()}

            threeJets_cut = selection.require(lep_and_filter_pass=True, passing_jets=True, jets_3=True)
            btagger_3j_wts = self.corrections["BTag_Constructors"][btaggers[0]]["3Jets"].get_scale_factor(jets=events["SelectedJets"][threeJets_cut], passing_cut=btag_wp)
            fourplusJets_cut = selection.require(lep_and_filter_pass=True, passing_jets=True, jets_4p=True)
            btagger_4pj_wts = self.corrections["BTag_Constructors"][btaggers[0]]["4PJets"].get_scale_factor(jets=events["SelectedJets"][fourplusJets_cut], passing_cut=btag_wp)

                # fll dict of btag weights
            for wt_name in btagger_3j_wts.keys():
                btag_weights[wt_name][threeJets_cut] = ak.prod(btagger_3j_wts[wt_name], axis=1)
                btag_weights[wt_name][fourplusJets_cut] = ak.prod(btagger_4pj_wts[wt_name], axis=1)
        elif (corrections["BTagSF"] == False):
            btag_weights = {"central" : np.ones(len(events))}
            print("BTag SFs not applied to MC")



        ## fill hists for each region
        for lepton in self.regions.keys():
            evt_weights = mu_evt_weights if lepton == "Muon" else el_evt_weights
            for leptype in self.regions[lepton].keys():
                for jmult in self.regions[lepton][leptype].keys():

                    if to_debug: print(lepton, leptype, "btagPass", jmult)
                    cut = selection.all(*self.regions[lepton][leptype][jmult])
                    output["cutflow"]["nEvts %s" % ", ".join([lepton, leptype, "btagPass", jmult])] += cut.sum()
                    if cut.sum() == 0: continue

                    ltype = "MU" if lepton == "Muon" else "EL"
                    if "loose_or_tight_%s" % ltype in self.regions[lepton][leptype][jmult]:
                        leptons = events[lepton][cut][((events[lepton][cut]["TIGHT%s" % ltype] == True) | (events[lepton][cut]["LOOSE%s" % ltype] == True))]
                    elif "tight_%s" % ltype in self.regions[lepton][leptype][jmult]:
                        leptons = events[lepton][cut][(events[lepton][cut]["TIGHT%s" % ltype] == True)]
                    elif "loose_%s" % ltype in self.regions[lepton][leptype][jmult]:
                        leptons = events[lepton][cut][(events[lepton][cut]["LOOSE%s" % ltype] == True)]
                    else:
                        raise ValueError("Not sure what lepton type to choose for event")

                        # get jets and MET
                    jets, met = events["SelectedJets"][cut], events["SelectedMET"][cut]

                        ## create MT regions
                    MT = make_vars.MT(leptons, met)
                    MTHigh = ak.flatten(MT >= MTcut)
                    output["cutflow"]["nEvts %s: pass MT cut" % ", ".join([lepton, leptype, "btagPass", jmult])] += ak.sum(MTHigh)

                        # find matched permutations
                    mp = ttmatcher.best_match(gen_hyp=events["SL"][cut], jets=jets, leptons=leptons, met=met)

                    wts = (evt_weights.weight()*btag_weights["central"])[cut][MTHigh]
                    #wts = (evt_weights.weight()*btag_weights["%s_CEN" % btaggers[0]])[cut][MTHigh]
                    valid_gen_objs = events["SL"][cut][MTHigh]
                    valid_mp = mp[MTHigh]

                    ## fill hists of gen mttbar for all events
                    output["mtt"].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, cat="All", mtt=ak.flatten(valid_gen_objs["TTbar"].mass, axis=None), weight=wts)
                        # for each gen parton find ones that pass and fail jet pt and eta cuts
                    for gen_obj in ["BHad", "BLep", "WJa", "WJb"]:
                        pt_eta_mask = ak.flatten((valid_gen_objs[gen_obj].pt >= jet_pars["ptmin"]) & (np.abs(valid_gen_objs[gen_obj].eta) <= jet_pars["etamax"]))
                        output["pass_mtt"].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, cat="All", objtype=gen_obj, mtt=ak.flatten(valid_gen_objs[pt_eta_mask]["TTbar"].mass, axis=None), weight=wts[pt_eta_mask])
                        output["fail_mtt"].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, cat="All", objtype=gen_obj, mtt=ak.flatten(valid_gen_objs[~pt_eta_mask]["TTbar"].mass, axis=None), weight=wts[~pt_eta_mask])

                    ## fill hists of gen mttbar for each lost jet/partially merged event category
                    evt_cats = ["Merged_BHadBLep", "Merged_BHadWJa", "Merged_BHadWJb", "Merged_BLepWJa", "Merged_BLepWJb", "Merged_WJets", "Lost_BHad", "Lost_BLep", "Lost_WJa", "Lost_WJb"]
                    for evt_cat in evt_cats:
                        evt_cat_mask = ak.flatten(valid_mp[evt_cat])
                        output["mtt"].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, cat=evt_cat, mtt=ak.flatten(valid_gen_objs[evt_cat_mask]["TTbar"].mass, axis=None), weight=wts[evt_cat_mask])

                            # for each gen parton find ones that pass and fail jet pt and eta cuts
                        for gen_obj in ["BHad", "BLep", "WJa", "WJb"]:
                            pt_eta_mask = ak.flatten((valid_gen_objs[evt_cat_mask][gen_obj].pt >= jet_pars["ptmin"]) & (np.abs(valid_gen_objs[evt_cat_mask][gen_obj].eta) <= jet_pars["etamax"]))
                            output["pass_mtt"].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, cat=evt_cat, objtype=gen_obj, mtt=ak.flatten(valid_gen_objs[evt_cat_mask][pt_eta_mask]["TTbar"].mass, axis=None), weight=wts[evt_cat_mask][pt_eta_mask])
                            output["fail_mtt"].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, cat=evt_cat, objtype=gen_obj, mtt=ak.flatten(valid_gen_objs[evt_cat_mask][~pt_eta_mask]["TTbar"].mass, axis=None), weight=wts[evt_cat_mask][~pt_eta_mask])
                        
                    ## fill hists of gen mttbar for "other" (not merged/lost) events
                    other_cat_mask = ak.flatten(~(valid_mp["Merged_Event"] | valid_mp["Lost_Event"]))
                    output["mtt"].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, cat="Other", mtt=ak.flatten(valid_gen_objs[other_cat_mask]["TTbar"].mass, axis=None), weight=wts[other_cat_mask])
                        # for each gen parton find ones that pass and fail jet pt and eta cuts
                    for gen_obj in ["BHad", "BLep", "WJa", "WJb"]:
                        pt_eta_mask = ak.flatten((valid_gen_objs[other_cat_mask][gen_obj].pt >= jet_pars["ptmin"]) & (np.abs(valid_gen_objs[other_cat_mask][gen_obj].eta) <= jet_pars["etamax"]))
                        output["pass_mtt"].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, cat="Other", objtype=gen_obj, mtt=ak.flatten(valid_gen_objs[other_cat_mask][pt_eta_mask]["TTbar"].mass, axis=None), weight=wts[other_cat_mask][pt_eta_mask])
                        output["fail_mtt"].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, cat="Other", objtype=gen_obj, mtt=ak.flatten(valid_gen_objs[other_cat_mask][~pt_eta_mask]["TTbar"].mass, axis=None), weight=wts[other_cat_mask][~pt_eta_mask])



        return output



    def postprocess(self, accumulator):
        return accumulator

proc_executor = processor.iterative_executor if to_debug else processor.futures_executor
proc_exec_args = {"schema": NanoAODSchema} if to_debug else {"schema": NanoAODSchema, "workers": 8}
output = processor.run_uproot_job(
    fileset,
    treename="Events",
    processor_instance=ttbar_evt_cats(),
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
