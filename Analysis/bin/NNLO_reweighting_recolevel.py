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
import python.IDJet as IDJet

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

nominal_ttJets_ = ["ttJetsSL", "ttJetsHad", "ttJetsDiLep"]
isNominal_ttJets_ = samplename in nominal_ttJets_
if not isNominal_ttJets_:
    raise ValueError("Should only be run on ttbar datasets")

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
base_jobid = os.environ["base_jobid"]

## get dataset classification, used for corrections/systematics
isTTbar_ = samplename.startswith("ttJets")
isTTSL_ = samplename.startswith("ttJetsSL")

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

reweight_systematics_to_run = ["nosys", "TopPtUp", "TopPtDown", "Mtt_vs_top_ctstarUp", "Mtt_vs_top_ctstarDown"]
print("\tRunning reweight systematics: " + ", ".join(reweight_systematics_to_run))
#set_trace()


# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class htt_btag_sb_regions(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.sys_axis = hist.Cat("sys", "Systematic")
        self.jetmult_axis = hist.Cat("jmult", "nJets")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")

        self.btag_discr_axis = hist.Bin("bdisc", f"{btaggers[0]} bDiscr", np.around(np.linspace(0., 1., 101), decimals=2))
        self.pt_axis = hist.Bin("pt", "p_{T} [GeV]", np.around(np.linspace(0., 1000., 101), decimals=0))
        self.eta_axis = hist.Bin("eta", r"$\eta$", np.around(np.linspace(-5., 5., 101), decimals=1))
        self.eta_2d_axis = hist.Bin("eta_2d", r"$\eta$", np.around(np.array([-3.0, -2.5, -1.3, -0.7, 0., 0.7, 1.3, 2.5, 3.0]), decimals=1))
        self.phi_axis = hist.Bin("phi", r"$\phi$", np.around(np.linspace(-4., 4., 81), decimals=1))
        self.phi_2d_axis = hist.Bin("phi_2d", r"$\phi$", np.around(np.array([-3.2, -2.4, -1.57, -0.87, 0., 0.87, 1.57, 2.4, 3.2]), decimals=2))
        #self.energy_axis = hist.Bin("energy", "E [GeV]", 200, 0, 1000)
        self.njets_axis = hist.Bin("njets", "n_{jets}", np.around(np.linspace(0., 20., 21), decimals=0))
        self.lepIso_axis = hist.Bin("iso", "pfRelIso", np.around(np.linspace(0., 1., 101), decimals=2))
        self.mt_axis = hist.Bin("mt", "M_{T}", np.around(np.linspace(0., 1000., 201), decimals=0))
        self.mtop_axis = hist.Bin("mtop", "m(top) [GeV]", np.around(np.linspace(0., 300., 301), decimals=0))
        self.wmass_axis = hist.Bin("wmass", "m(W_{had}) [GeV]", np.around(np.linspace(0., 300., 301), decimals=0))
        self.mtt_axis = hist.Bin("mtt", "m($t\overline{t}$) [GeV]", np.around(np.linspace(200., 2000., 361), decimals=0))
        self.probDisc_axis = hist.Bin("prob", "$\lambda_{C}$", np.around(np.linspace(0., 30., 301), decimals=1))
        self.massDisc_axis = hist.Bin("massdisc", "$\lambda_{M}$", np.around(np.linspace(0., 30., 301), decimals=1))
        self.nsDisc_axis = hist.Bin("nsdisc", "$\lambda_{NS}$", np.around(np.linspace(0., 30., 301), decimals=1))
        self.nu_dist_axis = hist.Bin("nu_dist", r"$D_{\nu, min}$", np.around(np.linspace(0., 150., 151), decimals=0))
        self.ctstar_axis = hist.Bin("ctstar", "cos($\\theta^{*}$)", np.around(np.linspace(-1., 1., 41), decimals=2))
        self.ctstar_abs_axis = hist.Bin("ctstar_abs", "|cos($\\theta^{*}$)|", np.around(np.linspace(0., 1., 21), decimals=2))

            ## make dictionary of hists
        histo_dict = {}
                ## make jet hists
        jet_hists = self.make_jet_hists()
        histo_dict.update(jet_hists)
                ## make lepton hists
        lep_hists = self.make_lep_hists()
        histo_dict.update(lep_hists)        
                ## make selection plots
        selection_hists = self.make_selection_hists()
        histo_dict.update(selection_hists)        

        self.sample_name = ""
        self.corrections = corrections
        self.reweight_systematics_to_run = reweight_systematics_to_run

            ## make dict of cutflow for each systematic variation
        histo_dict["cutflow"] = processor.defaultdict_accumulator(int)
    
        self._accumulator = processor.dict_accumulator(histo_dict)

            ## define regions depending on data or MC, systematic or nosys
        self.regions = {
            "Muon" : {
                "3Jets"  : {"lep_and_filter_pass", "passing_jets", "jets_3" , "tight_MU", f"{btaggers[0]}_pass"},
                "4PJets" : {"lep_and_filter_pass", "passing_jets", "jets_4p", "tight_MU", f"{btaggers[0]}_pass"},
            },
            "Electron" : {
                "3Jets"  : {"lep_and_filter_pass", "passing_jets", "jets_3" , "tight_EL", f"{btaggers[0]}_pass"},
                "4PJets" : {"lep_and_filter_pass", "passing_jets", "jets_4p", "tight_EL", f"{btaggers[0]}_pass"},
            },
        }

        if isTTSL_:
            for lepton in self.regions.keys():
                for jmult in self.regions[lepton].keys():
                    self.regions[lepton][jmult].update({"semilep"})


    @property
    def accumulator(self):
        return self._accumulator


    def make_jet_hists(self):
        histo_dict = {}
        histo_dict["Jets_pt"]    = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.pt_axis)
        histo_dict["Jets_eta"]   = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.eta_axis)
        histo_dict["Jets_phi"]   = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.phi_axis)
        histo_dict["Jets_njets"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.njets_axis)
        histo_dict["Jets_LeadJet_pt"]    = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.pt_axis)
        histo_dict["Jets_LeadJet_eta"]   = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.eta_axis)
        histo_dict["Jets_phi_vs_eta"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.phi_2d_axis, self.eta_2d_axis)
        histo_dict[f"Jets_{btaggers[0]}_bDisc"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_discr_axis)

        return histo_dict

    def make_lep_hists(self):
        histo_dict = {}
        histo_dict["Lep_pt"]    = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.pt_axis)
        histo_dict["Lep_eta"]   = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.eta_axis)
        histo_dict["Lep_iso"]   = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.lepIso_axis)
        histo_dict["Lep_phi"]   = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.phi_axis)
        histo_dict["Lep_phi_vs_eta"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.phi_2d_axis, self.eta_2d_axis)

        return histo_dict

    def make_selection_hists(self):
        histo_dict = {}
        histo_dict["mtt"]      = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.mtt_axis)
        histo_dict["mthad"]    = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.mtop_axis)
        histo_dict["mWHad"]    = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.wmass_axis)
        histo_dict["mWLep"]    = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.wmass_axis)
        histo_dict["pt_thad"]  = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.pt_axis)
        histo_dict["pt_tlep"]  = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.pt_axis)
        histo_dict["pt_tt"]    = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.pt_axis)
        histo_dict["eta_thad"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.eta_axis)
        histo_dict["eta_tlep"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.eta_axis)
        histo_dict["eta_tt"]   = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.eta_axis)

        histo_dict["tlep_ctstar"]     = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.ctstar_axis)
        histo_dict["tlep_ctstar_abs"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.ctstar_abs_axis)

        histo_dict["full_disc"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.probDisc_axis)
        histo_dict["mass_disc"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.massDisc_axis)
        histo_dict["ns_disc"]   = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.nsDisc_axis)
        histo_dict["ns_dist"]   = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.nu_dist_axis)

        histo_dict["MT"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.mt_axis)

        histo_dict["MET_pt"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.pt_axis)
        histo_dict["MET_phi"]= hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.phi_axis)

        histo_dict["mtt_vs_tlep_ctstar_abs"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.mtt_axis, self.ctstar_abs_axis)

        return histo_dict


    def process(self, events):
        output = self.accumulator.identity()

        self.sample_name = events.metadata["dataset"]

            ## make event weights
                # data or MC distinction made internally
        mu_evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections, isTTbar=isNominal_ttJets_)
        el_evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections, isTTbar=isNominal_ttJets_)

            ## initialize selections
        selection = PackedSelection()

            # get all passing leptons
        lep_and_filter_pass = objsel.select_leptons(events, year=args.year, hem_15_16=apply_hem)
        selection.add("lep_and_filter_pass", lep_and_filter_pass) # add passing leptons requirement to all systematics

            ## build corrected jets and MET
        events["CorrectedJets"], events["CorrectedMET"] = IDJet.process_jets(events, args.year, self.corrections["JetCor"])

        ## add different selections
                ## muons
        tight_mu_sel, loose_mu_sel = ak.sum(events["Muon"]["TIGHTMU"], axis=1) == 1, ak.sum(events["Muon"]["LOOSEMU"], axis=1) == 1
        selection.add("tight_MU", tight_mu_sel) # one muon passing TIGHT criteria
                ## electrons
        tight_el_sel, loose_el_sel = ak.sum(events["Electron"]["TIGHTEL"], axis=1) == 1, ak.sum(events["Electron"]["LOOSEEL"], axis=1) == 1
        selection.add("tight_EL", tight_el_sel) # one electron passing TIGHT criteria

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
                # add mtt x ctstar reweighting
                nnlo_2d_wts = MCWeights.get_nnlo_weights(self.corrections["NNLO_Rewt"], events)
                mu_evt_weights.add("Mtt_vs_top_ctstar",
                    np.ones(len(events)),
                    np.copy(nnlo_2d_wts),
                )
                el_evt_weights.add("Mtt_vs_top_ctstar",
                    np.ones(len(events)),
                    np.copy(nnlo_2d_wts),
                )

            if "EWK_Rewt" in self.corrections.keys():
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

                # add top pt reweighting
            toppt_weights = ak.flatten(MCWeights.get_TopPt_weights(events), axis=None)
            mu_evt_weights.add("TopPt",
                np.ones(len(events)),
                np.copy(toppt_weights),
            )
            el_evt_weights.add("TopPt",
                np.ones(len(events)),
                np.copy(toppt_weights),
            )

            if isTTSL_:
                genpsel.select(events, mode="NORMAL")
                selection.add("semilep", ak.num(events["SL"]) > 0)
            else:
                selection.add("semilep", np.zeros(len(events), dtype=bool))

            # run over systematics that require changes to event objects (jets+MET)
        output["cutflow"]["lep_and_filter_pass"] += ak.sum(lep_and_filter_pass)
            # jet selection
        passing_jets = objsel.jets_selection(events, year=args.year, cutflow=output["cutflow"], shift="nosys", hem_15_16=apply_hem, met_factory=self.corrections["JetCor"]["MC"]["METFactory"])
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
            #raise ValueError("BTag SFs not applied to MC")

        ## fill hists for each region
        for lepton in self.regions.keys():
            evt_weights = mu_evt_weights if lepton == "Muon" else el_evt_weights
            lep_SFs = muSFs_dict if lepton == "Muon" else elSFs_dict
            for jmult in self.regions[lepton].keys():
                cut = selection.all(*self.regions[lepton][jmult])

                output["cutflow"]["nEvts %s" % ", ".join([lepton, jmult])] += cut.sum()

                if to_debug: print(lepton, jmult)
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

                        # get jets and MET
                    jets, met = events["SelectedJets"][cut], events["SelectedMET"][cut]

                        # find best permutations
                    best_perms = ttpermutator.find_best_permutations(jets=jets, leptons=leptons, MET=met, btagWP=btag_wp, btag_req=True)
                    valid_perms = ak.num(best_perms["TTbar"].pt) > 0
                    output["cutflow"]["nEvts %s: valid perms" % ", ".join([lepton, jmult])] += ak.sum(valid_perms)

                        ## create MT regions
                    MT = make_vars.MT(leptons, met)
                    MTHigh = ak.flatten(MT[valid_perms] >= MTcut)
                    output["cutflow"]["nEvts %s: pass MT cut" % ", ".join([lepton, jmult])] += ak.sum(MTHigh)

                        # fill hists for each systematic
                    for rewt_sys in self.reweight_systematics_to_run:
                        if to_debug: print(f"\t\tsysname: {rewt_sys}")

                        if rewt_sys == "nosys":
                            wts = (evt_weights.weight() * btag_weights["central"] * lep_SFs["central"])[cut][valid_perms][MTHigh]
                        elif rewt_sys.startswith("btag"):
                            wts = (evt_weights.weight() * btag_weights[rewt_sys.split("btag_")[-1]] * lep_SFs["central"])[cut][valid_perms][MTHigh]
                        elif rewt_sys.startswith("Lep"):
                            if rewt_sys.split("_")[-1] in lep_SFs.keys():
                                wts = (evt_weights.weight() * btag_weights["central"] * lep_SFs[rewt_sys.split("_")[-1]])[cut][valid_perms][MTHigh]
                            else:
                                print(f"{rewt_sys.split('_')[-1]} not found in {lepton} SF dict. Skipping")
                                continue
                        else:
                            if rewt_sys not in evt_weights.variations:
                                print(f"{rewt_sys} not option in event weights. Skipping")
                                continue
                            wts = (evt_weights.weight(rewt_sys) * btag_weights["central"] * lep_SFs["central"])[cut][valid_perms][MTHigh]

                        sysname = rewt_sys

                        output = self.fill_hists(acc=output, sys=sysname, jetmult=jmult, leptype=lepton,
                            perm=best_perms[valid_perms][MTHigh], jets=jets[valid_perms][MTHigh], leptons=leptons[valid_perms][MTHigh], MTvals=MT[valid_perms][MTHigh], evt_wts=wts)

        return output

    def fill_hists(self, acc, sys, jetmult, leptype, perm, jets, leptons, MTvals, evt_wts):
            ## apply alpha correction for 3Jets
        if (jetmult == "3Jets") and ("Alpha" in self.corrections):
            alpha_corr = self.corrections["Alpha"](172.5/perm["THad"].mass)
            perm["THad"] = perm["THad"].multiply(alpha_corr) # correct thad
            perm["TTbar"] = ak.flatten(perm["THad"]+perm["TLep"]) # correct ttbar

        thad_ctstar, tlep_ctstar = make_vars.ctstar(perm["THad"], perm["TLep"], flatten=True)

        pt_sorted_jets = jets[ak.argsort(jets.pt, ascending=False)]

        acc["mtt_vs_tlep_ctstar_abs"].fill(dataset=self.sample_name, sys=sys,  jmult=jetmult, leptype=leptype,
            mtt=ak.flatten(perm["TTbar"].mass), ctstar_abs=np.abs(tlep_ctstar), weight=evt_wts)

        acc["mtt"].fill(     dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype, mtt=ak.flatten(perm["TTbar"].mass), weight=evt_wts)
        acc["mthad"].fill(   dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype, mtop=ak.flatten(perm["THad"].mass), weight=evt_wts)
        acc["mWHad"].fill(   dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype, wmass=ak.flatten(perm["WHad"].mass), weight=evt_wts)
        acc["mWLep"].fill(   dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype, wmass=ak.flatten(perm["WLep"].mass), weight=evt_wts)
        acc["pt_thad"].fill( dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype, pt=ak.flatten(perm["THad"].pt), weight=evt_wts)
        acc["pt_tlep"].fill( dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype, pt=ak.flatten(perm["TLep"].pt), weight=evt_wts)
        acc["pt_tt"].fill(   dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype, pt=ak.flatten(perm["TTbar"].pt), weight=evt_wts)
        acc["eta_thad"].fill(dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype, eta=ak.flatten(perm["THad"].eta), weight=evt_wts)
        acc["eta_tlep"].fill(dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype, eta=ak.flatten(perm["TLep"].eta), weight=evt_wts)
        acc["eta_tt"].fill(  dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype, eta=ak.flatten(perm["TTbar"].eta), weight=evt_wts)

        acc["tlep_ctstar"].fill(    dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype, ctstar=tlep_ctstar, weight=evt_wts)
        acc["tlep_ctstar_abs"].fill(dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype, ctstar_abs=np.abs(tlep_ctstar), weight=evt_wts)

        acc["full_disc"].fill(dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype, prob=ak.flatten(perm["Prob"]), weight=evt_wts)
        acc["mass_disc"].fill(dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype, massdisc=ak.flatten(perm["MassDiscr"]), weight=evt_wts)
        acc["ns_disc"].fill(  dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype, nsdisc=ak.flatten(perm["NuDiscr"]), weight=evt_wts)
        acc["ns_dist"].fill(  dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype, nu_dist=ak.flatten(np.sqrt(perm["Nu"].chi2)), weight=evt_wts)

        acc["MET_pt"].fill( dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype, pt=ak.flatten(perm["MET"].pt), weight=evt_wts)
        acc["MET_phi"].fill(dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype, phi=ak.flatten(perm["MET"].phi), weight=evt_wts)

        acc["Jets_pt"].fill(   dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype, pt=ak.flatten(jets.pt), weight=ak.flatten((ak.ones_like(jets.pt)*evt_wts)))
        acc["Jets_eta"].fill(  dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype, eta=ak.flatten(jets.eta), weight=ak.flatten((ak.ones_like(jets.eta)*evt_wts)))
        acc["Jets_phi"].fill(  dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype, phi=ak.flatten(jets.phi), weight=ak.flatten((ak.ones_like(jets.phi)*evt_wts)))
        acc["Jets_njets"].fill(dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype, njets=ak.num(jets), weight=evt_wts)
        acc["Jets_phi_vs_eta"].fill(dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype,
            phi_2d=ak.flatten(jets.phi), eta_2d=ak.flatten(jets.eta), weight=ak.flatten((ak.ones_like(jets.phi)*evt_wts)))

        acc["Jets_LeadJet_pt"].fill( dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype, pt=pt_sorted_jets.pt[:, 0], weight=evt_wts)
        acc["Jets_LeadJet_eta"].fill(dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype, eta=pt_sorted_jets.eta[:, 0], weight=evt_wts)
        acc[f"Jets_{btaggers[0]}_bDisc"].fill(dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype, 
            bdisc=ak.flatten(jets[IDJet.btag_tagger_to_disc_name[btaggers[0]]]), weight=ak.flatten((ak.ones_like(jets.pt)*evt_wts)))

        acc["Lep_pt"].fill(    dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype, pt=ak.flatten(leptons.pt), weight=evt_wts)
        acc["Lep_eta"].fill(   dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype, eta=ak.flatten(leptons.eta), weight=evt_wts)
        acc["Lep_phi"].fill(   dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype, phi=ak.flatten(leptons.phi), weight=evt_wts)
        acc["Lep_iso"].fill(   dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype, 
            iso=ak.flatten(leptons["pfRelIso04_all"]) if leptype == "Muon" else ak.flatten(leptons["pfRelIso03_all"]), weight=evt_wts)
        acc["Lep_phi_vs_eta"].fill(dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype, 
            phi_2d=ak.flatten(leptons.phi), eta_2d=ak.flatten(leptons.eta), weight=evt_wts)

        acc["MT"].fill(dataset=self.sample_name, sys=sys, jmult=jetmult, leptype=leptype, mt=ak.flatten(MTvals), weight=evt_wts)

        return acc        

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
    processor_instance=htt_btag_sb_regions(),
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
