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
import os#, sys
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
import fnmatch
import python.IDJet as IDJet

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
base_jobid = os.environ["base_jobid"]
analyzer = "ttbar_post_alpha_reco"

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

    ## specify ttJets samples
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
                    #output = check_output(["xrdcp", "-f", f"{cp_rfile}", "/tmp"], timeout=1200, stderr=STDOUT)
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
alpha_corrections = load(os.path.join(proj_dir, "Corrections", base_jobid, cfg_pars["corrections"]["alpha"]))[args.year]
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
if corrections["BTagSF"] == True:
    sf_file = os.path.join(proj_dir, "Corrections", jobid, jet_pars["btagging"]["btagSF_file"])
    if not os.path.isfile(sf_file):
        raise IOError(f"BTag SF file {sf_file} does not exist")

    btag_sfs = load(sf_file)
    corrections.update({"BTag_Constructors" : {}})
    for btagger in btaggers:
        threeJets = btag_sfs[args.year][btagger]["3Jets"][wps_to_use[0]]
        fourPJets = btag_sfs[args.year][btagger]["4PJets"][wps_to_use[0]]
        corrections["BTag_Constructors"].update({btagger : {"3Jets" : threeJets, "4PJets" : fourPJets} })

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

    ## determine which systematics to run
sys_opts = ["nosys", "JER_UP", "JER_DW", "JES_Total_UP", "JES_Total_DW"]
evt_sys_to_run = opts_dict.get("evt_sys", "NONE")
event_systematics_to_run = ["nosys"] if (evt_sys_to_run == "NONE") else [name for name in sys_opts if fnmatch.fnmatch(name, evt_sys_to_run)]
print("Running with event systematics:", *sorted(event_systematics_to_run), sep=", ")

# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class ttbar_post_alpha_reco(processor.ProcessorABC):
    def __init__(self):

        self.sample_name = ""
        self.corrections = corrections
        self.event_systematics_to_run = event_systematics_to_run

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.sys_axis = hist.Cat("sys", "Systematic")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")
        self.jetmult_axis = hist.Cat("jmult", "nJets")
        self.corrtype_axis = hist.Cat("corrtype", "Correction Type")
        self.mtt_axis  = hist.Bin("mtt", "m($t\overline{t}$) [GeV]", np.around(np.linspace(200., 2000., 181), decimals=0))
        self.mthad_axis= hist.Bin("mthad", "m($t_{had}$) [GeV]", np.around(np.linspace(0., 500., 501), decimals=0))
        self.ctstar_axis = hist.Bin("ctstar", "cos($\\theta^{*}$)", np.around(np.linspace(-1., 1., 201), decimals=2))
        self.ctstar_abs_axis = hist.Bin("ctstar_abs", "|cos($\\theta^{*}$)|", np.around(np.linspace(0., 1., 201), decimals=3))
        self.reso_mtt_axis  = hist.Bin("reso_mtt", "", np.around(np.linspace(-500., 500., 1001), decimals=0))
        self.reso_mthad_axis= hist.Bin("reso_mthad", "", np.around(np.linspace(-500., 500., 1001), decimals=0))
        self.reso_ctstar_axis = hist.Bin("reso_ctstar", "", np.around(np.linspace(-2., 2., 801), decimals=3))
        self.reso_ctstar_abs_axis = hist.Bin("reso_ctstar_abs", "", np.around(np.linspace(-1., 1., 801), decimals=4))

            ## make dictionary of hists
        histo_dict = {}
        hists = self.make_hists()
        histo_dict.update(hists)

            ## make dict of cutflow for each systematic variation
        for sys in self.event_systematics_to_run:
            histo_dict[f"cutflow_{sys}"] = processor.defaultdict_accumulator(int)

            # alpha correction choices
        self.alpha_corr_choices = [
            ("E", "All_1D"), ("E", "All_2D"), ("E", "Mtt"),
            ("P", "All_1D"), ("P", "All_2D"), ("P", "Mtt")
        ]
        self._accumulator = processor.dict_accumulator(histo_dict)
        self.regions = {
            "Muon" : {
                "3Jets"  : {"lep_and_filter_pass", "passing_jets", "jets_3"  , "tight_MU", "btag_pass", "semilep"},
                "4Jets"  : {"lep_and_filter_pass", "passing_jets", "jets_4"  , "tight_MU", "btag_pass", "semilep"},
                "5PJets" : {"lep_and_filter_pass", "passing_jets", "jets_5p" , "tight_MU", "btag_pass", "semilep"},
            },
            "Electron" : {
                "3Jets"  : {"lep_and_filter_pass", "passing_jets", "jets_3"  , "tight_EL", "btag_pass", "semilep"},
                "4Jets"  : {"lep_and_filter_pass", "passing_jets", "jets_4"  , "tight_EL", "btag_pass", "semilep"},
                "5PJets" : {"lep_and_filter_pass", "passing_jets", "jets_5p" , "tight_EL", "btag_pass", "semilep"},
            },
        }

    
    @property
    def accumulator(self):
        return self._accumulator

    def make_hists(self):
        histo_dict = {}
        histo_dict["Reco_mtt"]             = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.corrtype_axis, self.mtt_axis)
        histo_dict["Reco_mthad"]           = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.corrtype_axis, self.mthad_axis)
        histo_dict["Reco_thad_ctstar"]     = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.corrtype_axis, self.ctstar_axis)
        histo_dict["Reco_thad_ctstar_abs"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.corrtype_axis, self.ctstar_abs_axis)
        histo_dict["Reso_mtt"]             = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.corrtype_axis, self.reso_mtt_axis)
        histo_dict["Reso_mthad"]           = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.corrtype_axis, self.reso_mthad_axis)
        histo_dict["Reso_thad_ctstar"]     = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.corrtype_axis, self.reso_ctstar_axis)
        histo_dict["Reso_thad_ctstar_abs"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.corrtype_axis, self.reso_ctstar_abs_axis)

        return histo_dict

    
    def process(self, events):
        output = self.accumulator.identity()

        self.sample_name = events.metadata["dataset"]

            ## make event weights
                # data or MC distinction made internally
        mu_evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections)
        el_evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections)

                ## initialize selections and regions
        selection = {sys_name: PackedSelection() for sys_name in self.event_systematics_to_run}

            ## build corrected jets and MET
        events["Jet"], events["MET"] = IDJet.process_jets(events, args.year, self.corrections["JetCor"])

            # find gen level particles for ttbar system
        genpsel.select(events, mode="NORMAL")
        {selection[sys].add("semilep", ak.num(events["SL"]) > 0) for sys in selection.keys()} # add semilep requirement to all systematics

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
                    #np.copy(ewk_wts_dict["Rebinned_KFactor_1.1"]),
                    #np.copy(ewk_wts_dict["Rebinned_KFactor_0.9"]),
                )
                el_evt_weights.add("Yukawa",  # really just varying value of Yt
                    np.copy(ewk_wts_dict["Rebinned_KFactor_1.0"]),
                    #np.copy(ewk_wts_dict["Rebinned_KFactor_1.1"]),
                    #np.copy(ewk_wts_dict["Rebinned_KFactor_0.9"]),
                )

            # get all passing leptons
        lep_and_filter_pass = objsel.select_leptons(events, year=args.year)
        {selection[sys].add("lep_and_filter_pass", lep_and_filter_pass) for sys in selection.keys()} # add passing leptons requirement to all systematics

            # add tight e/mu selections
        tight_mu_sel = ak.sum(events["Muon"]["TIGHTMU"], axis=1) == 1
        tight_el_sel = ak.sum(events["Electron"]["TIGHTEL"], axis=1) == 1
        {selection[sys].add("tight_MU", tight_mu_sel) for sys in selection.keys()} # one muon passing TIGHT criteria
        {selection[sys].add("tight_EL", tight_el_sel) for sys in selection.keys()} # one electron passing TIGHT criteria

        ### apply lepton SFs to MC (only applicable to tight leptons)
        if "LeptonSF" in self.corrections.keys():
            tight_muons = events["Muon"][tight_mu_sel][(events["Muon"][tight_mu_sel]["TIGHTMU"] == True)]
            muSFs_dict =  MCWeights.get_lepton_sf(sf_dict=self.corrections["LeptonSF"]["Muons"],
                pt=ak.flatten(tight_muons["pt"]), eta=ak.flatten(tight_muons["eta"]), tight_lep_mask=tight_mu_sel, leptype="Muons")

            tight_electrons = events["Electron"][tight_el_sel][(events["Electron"][tight_el_sel]["TIGHTEL"] == True)]
            elSFs_dict = MCWeights.get_lepton_sf(sf_dict=self.corrections["LeptonSF"]["Electrons"],
                pt=ak.flatten(tight_electrons["pt"]), eta=ak.flatten(tight_electrons["etaSC"]), tight_lep_mask=tight_el_sel, leptype="Electrons")

            # loop over systematics involving jets/MET
        for evt_sys in self.event_systematics_to_run:
            output["cutflow_%s" % evt_sys]["lep_and_filter_pass"] += ak.sum(lep_and_filter_pass)
                # jet selection
            passing_jets = objsel.jets_selection(events, year=args.year, cutflow=output[f"cutflow_{evt_sys}"], shift=evt_sys)
            #set_trace()
            output["cutflow_%s" % evt_sys]["nEvts passing jet and lepton obj selection"] += ak.sum(passing_jets & lep_and_filter_pass)
            selection[evt_sys].add("passing_jets", passing_jets)
            selection[evt_sys].add("jets_3",  ak.num(events["SelectedJets"]) == 3)
            selection[evt_sys].add("jets_4",  ak.num(events["SelectedJets"]) == 4)
            selection[evt_sys].add("jets_5p",  ak.num(events["SelectedJets"]) > 4)
            selection[evt_sys].add("jets_4p",  ak.num(events["SelectedJets"]) > 3) # only for getting btag weights
            selection[evt_sys].add("btag_pass", ak.sum(events["SelectedJets"][btag_wp], axis=1) >= 2)

                # sort jets by btag value
            events["SelectedJets"] = events["SelectedJets"][ak.argsort(events["SelectedJets"][IDJet.btag_tagger_to_disc_name[btaggers[0]]], ascending=False)]

            if self.corrections["BTagSF"] == True:
                btag_weights = {key : np.ones(len(events)) for key in self.corrections["BTag_Constructors"][btaggers[0]]["3Jets"].schema_.keys()}

                threeJets_cut = selection[evt_sys].require(lep_and_filter_pass=True, passing_jets=True, jets_3=True)
                btagger_3j_wts = self.corrections["BTag_Constructors"][btaggers[0]]["3Jets"].get_scale_factor(jets=events["SelectedJets"][threeJets_cut], passing_cut=btag_wp)

                fourplusJets_cut = selection[evt_sys].require(lep_and_filter_pass=True, passing_jets=True, jets_4p=True)
                btagger_4pj_wts = self.corrections["BTag_Constructors"][btaggers[0]]["4PJets"].get_scale_factor(jets=events["SelectedJets"][fourplusJets_cut], passing_cut=btag_wp)

                    # fll dict of btag weights
                for wt_name in btagger_3j_wts.keys():
                    btag_weights[wt_name][threeJets_cut] = np.copy(ak.prod(btagger_3j_wts[wt_name], axis=1))
                    btag_weights[wt_name][fourplusJets_cut] = np.copy(ak.prod(btagger_4pj_wts[wt_name], axis=1))


            ## fill hists for each region
            for lepton in self.regions.keys():
                evt_weights = mu_evt_weights if lepton == "Muon" else el_evt_weights
                lep_SFs = muSFs_dict if lepton == "Muon" else elSFs_dict
                for jmult in self.regions[lepton].keys():
                    cut = selection[evt_sys].all(*self.regions[lepton][jmult])

                    #set_trace()
                    output[f"cutflow_{evt_sys}"]["nEvts %s" % ", ".join([lepton, jmult])] += cut.sum()        
                    if cut.sum() > 0:
                        ltype = "MU" if lepton == "Muon" else "EL"
                        if f"loose_or_tight_{ltype}" in self.regions[lepton][jmult]:
                            leptons = events[lepton][cut][((events[lepton][cut][f"TIGHT{ltype}"] == True) | (events[lepton][cut][f"LOOSE{ltype}"] == True))]
                        elif f"tight_{ltype}" in self.regions[lepton][jmult]:
                            leptons = events[lepton][cut][(events[lepton][cut][f"TIGHT{ltype}"] == True)]
                        elif f"loose_{ltype}" in self.regions[lepton][jmult]:
                            leptons = events[lepton][cut][(events[lepton][cut][f"LOOSE{ltype}"] == True)]
                        else:
                            raise ValueError("Not sure what lepton type to choose for event")

                        #set_trace()
                            # get jets and met
                        jets, met = events["SelectedJets"][cut], events["SelectedMET"][cut]

                            # find matched permutations
                        mp = ttmatcher.best_match(gen_hyp=events["SL"][cut], jets=jets, leptons=leptons, met=met)
        
                            # find best permutations
                        best_perms = ttpermutator.find_best_permutations(jets=jets, leptons=leptons, MET=met, btagWP=btag_wp)
                        valid_perms = ak.num(best_perms["TTbar"].pt) > 0
                        output[f"cutflow_{evt_sys}"]["nEvts %s: valid perms" % ", ".join([lepton, jmult])] += ak.sum(valid_perms)

                            # compare matched per to best perm
                        bp_status = np.zeros(cut.size, dtype=int) # 0 == "" (no gen matching), 1 == "right", 2 == "matchable", 3 == "unmatchable", 4 == "sl_tau", 5 == "noslep"
                        perm_cat_array = compare_matched_best_perms(mp, best_perms, njets="3Jets" if jmult == "3Jets" else "4PJets")
                        bp_status[cut] = perm_cat_array
                        sl_tau_evts = ak.where(np.abs(events["SL"]["Lepton"].pdgId) == 15)[0]
                        bp_status[sl_tau_evts] = 4
        
                            ## create MT regions
                        MT = make_vars.MT(leptons, met)
                        MTHigh = ak.flatten(MT[valid_perms] >= MTcut)
                        output[f"cutflow_{evt_sys}"]["nEvts %s: pass MT cut" % ", ".join([lepton, jmult])] += ak.sum(MTHigh)

                        wts = (evt_weights.weight() * btag_weights["central"] * lep_SFs["central"])[cut][valid_perms][MTHigh]   

                        if to_debug: print(*[lepton, jmult, evt_sys], sep=", ")
                        output = self.fill_hists(acc=output, sys=evt_sys, jetmult=jmult, leptype=lepton, permarray=bp_status[cut][valid_perms][MTHigh], genttbar=events["SL"][cut][valid_perms][MTHigh], bp=best_perms[valid_perms][MTHigh], evt_wts=wts)

        return output


    def fill_hists(self, acc, sys, jetmult, leptype, permarray, genttbar, bp, evt_wts):
        #set_trace()
            ## get gen variables that don"t change
        gen_mtt = ak.flatten(genttbar["TTbar"].mass, axis=None)
        gen_mthad = ak.flatten(genttbar["THad"].mass, axis=None)
        gen_thad_ctstar, gen_tlep_ctstar = make_vars.ctstar(genttbar["THad"], genttbar["TLep"])
        gen_thad_ctstar = ak.flatten(gen_thad_ctstar, axis=None)

        orig_bp_thad_ctstar, orig_bp_tlep_ctstar = make_vars.ctstar(bp["THad"], bp["TLep"])
        orig_bp_thad_ctstar = ak.flatten(orig_bp_thad_ctstar)
        orig_bp_mthad = ak.flatten(bp["THad"].mass)
        orig_bp_mtt = ak.flatten(bp["TTbar"].mass)

            # construct dict of alpha corrections
        alpha_corr_vals_dict = {"_".join([kvar, mtt_bin]) : self.corrections["Alpha"][kvar][mtt_bin](172.5/orig_bp_mthad, orig_bp_mtt) if mtt_bin == "Mtt" else self.corrections["Alpha"][kvar][mtt_bin](172.5/orig_bp_mthad) for kvar, mtt_bin in self.alpha_corr_choices}
            # construct dict of corrected thad/ttbar values
        corr_bp_vals = {}
        for key in alpha_corr_vals_dict.keys():
            corr_bp_thad_p4 = bp["THad"].multiply(alpha_corr_vals_dict[key])
            corr_bp_thad_ctstar, corr_bp_tlep_ctstar = make_vars.ctstar(corr_bp_thad_p4, bp["TLep"])
            corr_bp_vals[key] = {
                "mthad" : ak.flatten(corr_bp_thad_p4.mass),
                "thad_ctstar" : ak.flatten(corr_bp_thad_ctstar, axis=None),
                "mtt" : ak.flatten((corr_bp_thad_p4 + bp["TLep"]).mass, axis=None)
            }


        for permval in np.unique(permarray).tolist():
            #set_trace()
            perm_inds = np.where(permarray == permval)
            dataset_name = "%s_%s" % (self.sample_name, perm_cats[permval]) if permval != 0 else self.sample_name

                # apply no correction
                # Fill hists with no correction
            acc["Reco_mtt"].fill(            dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, corrtype="Uncorrected", mtt=orig_bp_mtt[perm_inds], weight=evt_wts[perm_inds])
            acc["Reco_mthad"].fill(          dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, corrtype="Uncorrected", mthad=orig_bp_mthad[perm_inds], weight=evt_wts[perm_inds])
            acc["Reco_thad_ctstar"].fill(    dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, corrtype="Uncorrected", ctstar=orig_bp_thad_ctstar[perm_inds], weight=evt_wts[perm_inds])
            acc["Reco_thad_ctstar_abs"].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, corrtype="Uncorrected", ctstar_abs=np.abs(orig_bp_thad_ctstar[perm_inds]), weight=evt_wts[perm_inds])
            acc["Reso_mtt"].fill(            dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, corrtype="Uncorrected", reso_mtt=gen_mtt[perm_inds] - orig_bp_mtt[perm_inds], weight=evt_wts[perm_inds])
            acc["Reso_mthad"].fill(          dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, corrtype="Uncorrected", reso_mthad=gen_mthad[perm_inds] - orig_bp_mthad[perm_inds], weight=evt_wts[perm_inds])
            acc["Reso_thad_ctstar"].fill(    dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, corrtype="Uncorrected", reso_ctstar=gen_thad_ctstar[perm_inds] - orig_bp_thad_ctstar[perm_inds], weight=evt_wts[perm_inds])
            acc["Reso_thad_ctstar_abs"].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, corrtype="Uncorrected", reso_ctstar_abs=np.abs(gen_thad_ctstar[perm_inds]) - np.abs(orig_bp_thad_ctstar[perm_inds]), weight=evt_wts[perm_inds])

            if jetmult != "3Jets": continue

            #set_trace()
                ## apply alpha corrections only for 3Jets
            for corrtype in corr_bp_vals.keys():
                corr_mthad = corr_bp_vals[corrtype]["mthad"][perm_inds]
                corr_ctstar = corr_bp_vals[corrtype]["thad_ctstar"][perm_inds]
                corr_mtt = corr_bp_vals[corrtype]["mtt"][perm_inds]
                acc["Reco_mtt"].fill(            dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, corrtype=corrtype, mtt=corr_mtt, weight=evt_wts[perm_inds])
                acc["Reco_mthad"].fill(          dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, corrtype=corrtype, mthad=corr_mthad, weight=evt_wts[perm_inds])
                acc["Reco_thad_ctstar"].fill(    dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, corrtype=corrtype, ctstar=corr_ctstar, weight=evt_wts[perm_inds])
                acc["Reco_thad_ctstar_abs"].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, corrtype=corrtype, ctstar_abs=np.abs(corr_ctstar), weight=evt_wts[perm_inds])
                acc["Reso_mtt"].fill(            dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, corrtype=corrtype, reso_mtt=gen_mtt[perm_inds] - corr_mtt, weight=evt_wts[perm_inds])
                acc["Reso_mthad"].fill(          dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, corrtype=corrtype, reso_mthad=gen_mthad[perm_inds] - corr_mthad, weight=evt_wts[perm_inds])
                acc["Reso_thad_ctstar"].fill(    dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, corrtype=corrtype, reso_ctstar=gen_thad_ctstar[perm_inds] - corr_ctstar, weight=evt_wts[perm_inds])
                acc["Reso_thad_ctstar_abs"].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, corrtype=corrtype, reso_ctstar_abs=np.abs(gen_thad_ctstar[perm_inds]) - np.abs(corr_ctstar), weight=evt_wts[perm_inds])

        return acc


    def postprocess(self, accumulator):
        return accumulator


proc_executor = processor.iterative_executor if to_debug else processor.futures_executor
proc_exec_args = {"schema": processor.NanoAODSchema} if to_debug else {"schema": processor.NanoAODSchema, "workers": 8}
output = processor.run_uproot_job(
    fileset,
    treename="Events",
    processor_instance=ttbar_post_alpha_reco(),
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
