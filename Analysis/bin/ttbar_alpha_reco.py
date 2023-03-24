#!/usr/bin/env python

import time
tic = time.time()

from coffea import hist, processor
from coffea.analysis_tools import PackedSelection
import awkward as ak
from coffea.nanoevents.methods import vector
ak.behavior.update(vector.behavior)
from coffea.util import save, load
from pdb import set_trace
import os
import numpy as np
import Utilities.prettyjson as prettyjson
import python.GenParticleSelector as genpsel
import python.TTGenMatcher as ttmatcher
import python.TTPermutator as ttpermutator
import Utilities.make_variables as make_vars
import python.ObjectSelection as objsel
import python.MCWeights as MCWeights
from python.Permutations import compare_matched_best_perms
import python.IDJet as IDJet

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
base_jobid = os.environ["base_jobid"]
analyzer = "ttbar_alpha_reco"

print(f"Jobid: {jobid}")
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
samplename = list(fileset.keys())[0]

# convert input string of options dictionary to actual dictionary
odict = (args.opts).replace("\'", "\"")
opts_dict = prettyjson.loads(odict)

    ## set config options passed through argparse
import ast
to_debug = ast.literal_eval(opts_dict.get("debug", "False"))

    ## run on correct samples
Nominal_ttJets = ["ttJetsSL"]
isTTbar = np.array([(key in Nominal_ttJets) for key in fileset.keys()]).all()
if not isTTbar:
    raise ValueError("This should only be run on nominal ttbar events!")

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
                    #output = check_output(["xrdcp", "-f", f"{cp_rfile}", "/tmp"], timeout=600, stderr=STDOUT)
                    output = check_output(["xrdcp", "-f", f"{cp_rfile}", "/tmp"], timeout=300, stderr=STDOUT)
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
nnlo_reweighting = load(os.path.join(proj_dir, "Corrections", base_jobid, cfg_pars["corrections"]["nnlo"]["filename"]))
ewk_reweighting = load(os.path.join(proj_dir, "Corrections", base_jobid, cfg_pars["corrections"]["ewk"]["file"]))
corrections = {
    "Pileup" : pu_correction,
    "Prefire" : True,
    "LeptonSF" : lepSF_correction,
    "BTagSF" : True,
    "JetCor" : jet_corrections,
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
class ttbar_alpha_reco(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")
        self.bp_mtt_axis  = hist.Bin("bp_mtt", "m($t\overline{t}$)", np.around(np.array([200., 350., 400., 500., 700., 1000., 2000.]), decimals=0))
        #self.bp_mtt_axis  = hist.Bin("bp_mtt", "m($t\overline{t}$)", 36, 200., 2000.)
        self.alpha_P_axis   = hist.Bin("alpha_p", "Gen/Reco P($t_{h}$)", np.around(np.linspace(0., 5., 1001), decimals=3))
        self.alpha_E_axis   = hist.Bin("alpha_e", "Gen/Reco E($t_{h}$)", np.around(np.linspace(0., 5., 1001), decimals=3))
        self.norm_mthad_axis= hist.Bin("norm_mthad", "172.5/m($t_{had}$)", np.around(np.linspace(0., 5., 51), decimals=1))

            ## make dictionary of hists
        histo_dict = {}
            ## make jet hists
                # Gen
        alpha_hists = self.make_alpha_hists()
        histo_dict.update(alpha_hists)

        self._accumulator = processor.dict_accumulator(histo_dict)
        self.sample_name = ""
        self.corrections = corrections
        self.regions = {
            "Muon" : {
                "lep_and_filter_pass", "passing_jets", "jets_3" , "tight_MU", "btag_pass", "semilep"
            },
            "Electron" : {
                "lep_and_filter_pass", "passing_jets", "jets_3" , "tight_EL", "btag_pass", "semilep"
            },
        }

    
    @property
    def accumulator(self):
        return self._accumulator

    def make_alpha_hists(self):
        histo_dict = {}
        histo_dict["Alpha_THad_P"] = hist.Hist("Events", self.dataset_axis, self.leptype_axis, self.bp_mtt_axis, self.norm_mthad_axis, self.alpha_P_axis)
        histo_dict["Alpha_THad_E"] = hist.Hist("Events", self.dataset_axis, self.leptype_axis, self.bp_mtt_axis, self.norm_mthad_axis, self.alpha_E_axis)

        return histo_dict

    
    def process(self, events):
        output = self.accumulator.identity()

        self.sample_name = events.metadata["dataset"]

            ## make event weights
                # data or MC distinction made internally
        mu_evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections)
        el_evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections)

            ## initialize selections and regions
        selection = PackedSelection()

            # get all passing leptons
        lep_and_filter_pass = objsel.select_leptons(events, year=args.year)
        selection.add("lep_and_filter_pass", lep_and_filter_pass) # add passing leptons requirement to all systematics

        ## add different selections
                ## muons
        tight_mu_sel = ak.sum(events["Muon"]["TIGHTMU"], axis=1) == 1
        selection.add("tight_MU", tight_mu_sel) # one muon passing TIGHT criteria
                ## electrons
        tight_el_sel = ak.sum(events["Electron"]["TIGHTEL"], axis=1) == 1
        selection.add("tight_EL", tight_el_sel) # one electron passing TIGHT criteria

        ### apply lepton SFs to MC (only applicable to tight leptons)
        if "LeptonSF" in self.corrections.keys():
            tight_muons = events["Muon"][tight_mu_sel][(events["Muon"][tight_mu_sel]["TIGHTMU"] == True)]
            muSFs_dict =  MCWeights.get_lepton_sf(sf_dict=self.corrections["LeptonSF"]["Muons"],
                pt=ak.flatten(tight_muons["pt"]), eta=ak.flatten(tight_muons["eta"]), tight_lep_mask=tight_mu_sel, leptype="Muons")

            tight_electrons = events["Electron"][tight_el_sel][(events["Electron"][tight_el_sel]["TIGHTEL"] == True)]
            elSFs_dict = MCWeights.get_lepton_sf(sf_dict=self.corrections["LeptonSF"]["Electrons"],
                pt=ak.flatten(tight_electrons["pt"]), eta=ak.flatten(tight_electrons["etaSC"]), tight_lep_mask=tight_el_sel, leptype="Electrons")

            ## build corrected jets and MET
        events["Jet"], events["MET"] = IDJet.process_jets(events, args.year, self.corrections["JetCor"])

            # jet selection
        passing_jets = objsel.jets_selection(events, year=args.year)
        selection.add("passing_jets", passing_jets)
        selection.add("jets_3",  ak.num(events["SelectedJets"]) == 3)
        selection.add("btag_pass", ak.sum(events["SelectedJets"][btag_wp], axis=1) >= 2)

        events["SelectedJets"] = events["SelectedJets"][ak.argsort(events["SelectedJets"][IDJet.btag_tagger_to_disc_name[btaggers[0]]], ascending=False)]

        if self.corrections["BTagSF"] == True:
            btag_weights = {key : np.ones(len(events)) for key in self.corrections["BTag_Constructors"][btaggers[0]]["3Jets"].schema_.keys()}

            threeJets_cut = selection.require(lep_and_filter_pass=True, passing_jets=True, jets_3=True)
            btagger_3j_wts = self.corrections["BTag_Constructors"][btaggers[0]]["3Jets"].get_scale_factor(jets=events["SelectedJets"][threeJets_cut], passing_cut=btag_wp)

                # fll dict of btag weights
            for wt_name in btagger_3j_wts.keys():
                btag_weights[wt_name][threeJets_cut] = np.copy(ak.prod(btagger_3j_wts[wt_name], axis=1))


            # find gen level particles for ttbar system
        genpsel.select(events, mode="NORMAL")
        selection.add("semilep", ak.num(events["SL"]) > 0)

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
                    # add Yukawa coupling variation
                mu_evt_weights.add("Yukawa",  # really just varying value of Yt
                    np.copy(ewk_wts_dict["Rebinned_KFactor_1.0"]),
                    np.copy(ewk_wts_dict["Rebinned_KFactor_1.11"]),
                    np.copy(ewk_wts_dict["Rebinned_KFactor_0.88"]),
                    #np.copy(ewk_wts_dict["Rebinned_KFactor_1.1"]),
                    #np.copy(ewk_wts_dict["Rebinned_KFactor_0.9"]),
                )
                el_evt_weights.add("Yukawa",  # really just varying value of Yt
                    np.copy(ewk_wts_dict["Rebinned_KFactor_1.0"]),
                    np.copy(ewk_wts_dict["Rebinned_KFactor_1.11"]),
                    np.copy(ewk_wts_dict["Rebinned_KFactor_0.88"]),
                    #np.copy(ewk_wts_dict["Rebinned_KFactor_1.1"]),
                    #np.copy(ewk_wts_dict["Rebinned_KFactor_0.9"]),
                )

        ## fill hists for each region
        for lepton in self.regions.keys():
            evt_weights = mu_evt_weights if lepton == "Muon" else el_evt_weights
            lep_SFs = muSFs_dict if lepton == "Muon" else elSFs_dict
            cut = selection.all(*self.regions[lepton])

            #set_trace()
            if cut.sum() > 0:
                leptype = "MU" if lepton == "Muon" else "EL"
                if "loose_or_tight_%s" % leptype in self.regions[lepton]:
                    leptons = events[lepton][cut][((events[lepton][cut]["TIGHT%s" % leptype] == True) | (events[lepton][cut]["LOOSE%s" % leptype] == True))]
                elif "tight_%s" % leptype in self.regions[lepton]:
                    leptons = events[lepton][cut][(events[lepton][cut]["TIGHT%s" % leptype] == True)]
                elif "loose_%s" % leptype in self.regions[lepton]:
                    leptons = events[lepton][cut][(events[lepton][cut]["LOOSE%s" % leptype] == True)]
                else:
                    raise ValueError("Not sure what lepton type to choose for event")

                    # get jets and MET
                jets, met = events["SelectedJets"][cut], events["SelectedMET"][cut]

                    # find matched permutations
                mp = ttmatcher.best_match(gen_hyp=events["SL"][cut], jets=jets, leptons=leptons, met=met)

                    # find best permutations
                best_perms = ttpermutator.find_best_permutations(jets=jets, leptons=leptons, MET=met, btagWP=btag_wp)
                valid_perms = ak.num(best_perms["TTbar"].pt) > 0

                    # compare matched per to best perm
                bp_status = np.zeros(cut.size, dtype=int) # 0 == "" (no gen matching), 1 == "right", 2 == "matchable", 3 == "unmatchable", 4 == "sl_tau", 5 == "noslep"
                perm_cat_array = compare_matched_best_perms(mp, best_perms, njets="3Jets")
                bp_status[cut] = perm_cat_array
                if ak.any(ak.num(events["SL"]["Lepton"].pdgId) != 1): raise ValueError("Number of leptons is incorrect for classifying tau+jets events")
                sl_tau_evts = ak.where(np.abs(events["SL"]["Lepton"].pdgId) == 15)[0]
                bp_status[sl_tau_evts] = 4

                    ## create MT regions
                MT = make_vars.MT(leptons, met)
                MTHigh = ak.flatten(MT[valid_perms] >= MTcut)

                wts = (evt_weights.weight() * btag_weights["central"] * lep_SFs["central"])[cut][valid_perms][MTHigh]
                output = self.make_3j_categories(acc=output, leptype=lepton, permarray=bp_status[cut][valid_perms][MTHigh], genttbar=events["SL"][cut][valid_perms][MTHigh], bp=best_perms[valid_perms][MTHigh], evt_wts=wts)

        return output


    def make_3j_categories(self, acc, leptype, permarray, genttbar, bp, evt_wts):
            ## get variables
        bp_thad_E = ak.flatten(bp["THad"].energy)
        bp_thad_P = ak.flatten(bp["THad"].p)
        bp_thad_M = ak.flatten(bp["THad"].mass)
        bp_mtt    = ak.flatten(bp["TTbar"].mass)

        gen_thad_E = ak.flatten(genttbar["THad"].energy, axis=None)
        gen_thad_P = ak.flatten(genttbar["THad"].p, axis=None)

        for permval in np.unique(permarray).tolist():
            #set_trace()
            if permval != 1: continue # only fill hists for "right" events
            perm_inds = np.where(permarray == permval)
            dataset_name = "%s_%s" % (self.sample_name, perm_cats[permval]) if permval != 0 else self.sample_name

            acc["Alpha_THad_P"].fill(dataset=dataset_name, leptype=leptype, bp_mtt=bp_mtt[perm_inds], norm_mthad=172.5/bp_thad_M[perm_inds], alpha_p=gen_thad_P[perm_inds]/bp_thad_P[perm_inds], weight=evt_wts[perm_inds])
            acc["Alpha_THad_E"].fill(dataset=dataset_name, leptype=leptype, bp_mtt=bp_mtt[perm_inds], norm_mthad=172.5/bp_thad_M[perm_inds], alpha_e=gen_thad_E[perm_inds]/bp_thad_E[perm_inds], weight=evt_wts[perm_inds])


        return acc

    def postprocess(self, accumulator):
        return accumulator


proc_executor = processor.iterative_executor if to_debug else processor.futures_executor
proc_exec_args = {"schema": processor.NanoAODSchema} if to_debug else {"schema": processor.NanoAODSchema, "workers": 8}
output = processor.run_uproot_job(
    fileset,
    treename="Events",
    processor_instance=ttbar_alpha_reco(),
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
