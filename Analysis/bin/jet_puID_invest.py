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
import python.GenParticleSelector as genpsel
import python.TTGenMatcher as ttmatcher
import python.TTPermutator as ttpermutator
from python.Permutations import compare_matched_best_perms
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
isInt_ = isSignal_ and ("Int" in samplename)
isData_ = samplename.startswith("data")
isSE_Data_ = samplename.startswith("data_SingleElectron")
isSM_Data_ = samplename.startswith("data_SingleMuon")
if isData_:
    lumiMask_path = os.path.join(proj_dir, "inputs", "data", base_jobid, "LumiMasks", f"{args.year}_GoldenJson_{base_jobid}.txt")

Nominal_ttJets = ["ttJetsSL", "ttJetsHad", "ttJetsDiLep"]
isNominal_ttJets_ = samplename in Nominal_ttJets
isTTShift_ = isTTbar_ and (len(samplename.split("_")) > 1)

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
    #"BTagSF" : False,
    "BTagSF" : not isData_,
    "JetCor" : jet_corrections,
    "Alpha" : alpha_corrections,
    "NNLO_Rewt" : {"Var" : cfg_pars["corrections"]["nnlo"]["var"], "Correction" : nnlo_reweighting[cfg_pars["corrections"]["nnlo"]["var"]]},
    "EWK_Rewt" : {"Correction" : ewk_reweighting, "wt" : cfg_pars["corrections"]["ewk"]["wt"]},
}

    ## parameters for b-tagging
jet_pars = cfg_pars["Jets"]
btaggers = [jet_pars["btagger"]]
#set_trace()

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

    # sideband regions are determined by dividing deepcsv medium wp values by 3 for each year
btag_regions = {}

# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class htt_btag_sb_regions(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.sys_axis = hist.Cat("sys", "Systematic")
        self.jetmult_axis = hist.Cat("jmult", "nJets")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")
        self.btag_axis = hist.Cat("btag", "BTag Region")
        self.jID_axis = hist.Cat("jetID", "jet ID type")

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
        #self.ctstar_axis = hist.Bin("ctstar", "cos($\\theta^{*}$)", np.around(np.linspace(-1., 1., 201), decimals=2))
        #self.ctstar_abs_axis = hist.Bin("ctstar_abs", "|cos($\\theta^{*}$)|", np.around(np.linspace(0., 1., 201), decimals=3))

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
        #self.jetID_types = ["tightID+PU"]
        self.jetID_types = ["tightID", "tightID+PU"]

            ## make dict of cutflow for each systematic variation
        for jID in self.jetID_types:
            histo_dict[f"cutflow_{jID}"] = processor.defaultdict_accumulator(int)
    
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
        self.regions = {jID:deepcopy(base_regions) for jID in self.jetID_types}

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
            for jetID in self.regions.keys():
                for lepton in self.regions[jetID].keys():
                    for btagregion in self.regions[jetID][lepton].keys():
                        for jmult in self.regions[jetID][lepton][btagregion].keys():
                            self.regions[jetID][lepton][btagregion][jmult].update({"semilep"})


    @property
    def accumulator(self):
        return self._accumulator


    def make_jet_hists(self):
        histo_dict = {}
        histo_dict["Jets_pt"]        = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.jID_axis, self.pt_axis)
        histo_dict["Jets_eta"]       = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.jID_axis, self.eta_axis)
        histo_dict["Jets_phi"]       = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.jID_axis, self.phi_axis)
        histo_dict["Jets_njets"]     = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.jID_axis, self.njets_axis)
        histo_dict["Jets_phi_vs_eta"]= hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.jID_axis, self.phi_2d_axis, self.eta_2d_axis)
        histo_dict["Jets_LeadJet_pt"]    = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.jID_axis, self.pt_axis)
        histo_dict["Jets_LeadJet_eta"]   = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.jID_axis, self.eta_axis)
        histo_dict[f"Jets_{btaggers[0]}_bDisc"] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.jID_axis, self.btag_discr_axis)

        return histo_dict

    def make_lep_hists(self):
        histo_dict = {}
        histo_dict["Lep_pt"]    = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.jID_axis, self.pt_axis)
        histo_dict["Lep_eta"]   = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.jID_axis, self.eta_axis)
        histo_dict["Lep_iso"]   = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.jID_axis, self.lepIso_axis)
        histo_dict["Lep_phi"]   = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.jID_axis, self.phi_axis)
        histo_dict["Lep_phi_vs_eta"] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.jID_axis, self.phi_2d_axis, self.eta_2d_axis)

        return histo_dict

    def make_selection_hists(self):
        histo_dict = {}
        histo_dict["mtt"]      = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.jID_axis, self.mtt_axis)
        histo_dict["mthad"]    = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.jID_axis, self.mtop_axis)
        histo_dict["mWHad"]    = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.jID_axis, self.wmass_axis)
        histo_dict["mWLep"]    = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.jID_axis, self.wmass_axis)
        histo_dict["pt_thad"]  = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.jID_axis, self.pt_axis)
        histo_dict["pt_tlep"]  = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.jID_axis, self.pt_axis)
        histo_dict["pt_tt"]    = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.jID_axis, self.pt_axis)
        histo_dict["eta_thad"] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.jID_axis, self.eta_axis)
        histo_dict["eta_tlep"] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.jID_axis, self.eta_axis)
        histo_dict["eta_tt"]   = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.jID_axis, self.eta_axis)

        histo_dict["tlep_ctstar"]     = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.jID_axis, self.ctstar_axis)
        histo_dict["tlep_ctstar_abs"] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.jID_axis, self.ctstar_abs_axis)

        histo_dict["full_disc"] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.jID_axis, self.probDisc_axis)
        histo_dict["mass_disc"] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.jID_axis, self.massDisc_axis)
        histo_dict["ns_disc"]   = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.jID_axis, self.nsDisc_axis)
        histo_dict["ns_dist"]   = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.jID_axis, self.nu_dist_axis)

        histo_dict["MT"] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.jID_axis, self.mt_axis)

        histo_dict["MET_pt"] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.jID_axis, self.pt_axis)
        histo_dict["MET_phi"]= hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.jID_axis, self.phi_axis)

        histo_dict["mtt_vs_tlep_ctstar_abs"] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.jID_axis, self.mtt_axis, self.ctstar_abs_axis)

        return histo_dict


    def select_tightID_jets(self, jets, muons, electrons, cutflow=None):
        #set_trace()
            ## pt and eta cuts
        pass_pt_eta_cuts = IDJet.make_pt_eta_cuts(jets)
        if cutflow is not None: cutflow["jets pass pT and eta cuts"] += ak.sum(pass_pt_eta_cuts)
    
            ## tightID
        # check https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVUL for definitions
        jet_ID = (jets.isTight)  # pass tight ID 
        if cutflow is not None: cutflow["jets pass ID"] += ak.sum(jet_ID)
    
            ## remove jets that don"t pass ID and pt/eta cuts
        jets = jets[(jet_ID & pass_pt_eta_cuts)]
    
            ## clean jets wrt veto+loose+tight el and mu
        jets_ak = ak.with_name(jets[["pt", "eta", "phi", "mass"]],"PtEtaPhiMLorentzVector")
        muons_ak = ak.with_name(muons[(muons["TIGHTMU"] | muons["LOOSEMU"] | muons["VETOMU"]) == True][["pt", "eta", "phi", "mass"]],"PtEtaPhiMLorentzVector")
        electrons_ak = ak.with_name(electrons[(electrons["TIGHTEL"] | electrons["LOOSEEL"] | electrons["VETOEL"]) == True][["pt", "eta", "phi", "mass"]],"PtEtaPhiMLorentzVector")
    
                # all mu
        jets_akc, muons_akc = ak.unzip(ak.cartesian([jets_ak, muons_ak], nested=True))
        clean_j_Mu_mask = ak.all((jets_akc.delta_r(muons_akc) >= 0.4), axis=2)
                # all el
        jets_akc, electrons_akc = ak.unzip(ak.cartesian([jets_ak, electrons_ak], nested=True))
        clean_j_El_mask = ak.all((jets_akc.delta_r(electrons_akc) >= 0.4), axis=2)
        clean_jets_mask = (clean_j_Mu_mask & clean_j_El_mask)
        if ak.sum(clean_jets_mask) > 0: jets = jets[clean_jets_mask]
    
            ## leading jet pt cut
        leadpt_cut = IDJet.make_leadjet_pt_cut(jets)
        if cutflow is not None: cutflow["jets pass lead jet pT cut"] += ak.sum(leadpt_cut)
    
            ## 3 or more jets
        njet_restriction = 3
        njets_cuts = (ak.num(jets) >= njet_restriction)
        if cutflow is not None: cutflow["nEvts with %s+ clean jets passing ID and kin selection" % njet_restriction] += ak.sum(njets_cuts)
    
        passing_jets = (leadpt_cut & njets_cuts)
        if cutflow is not None: cutflow["nEvts with jets passing selection"] += ak.sum(passing_jets)
    
        if cutflow is not None:
            return jets, passing_jets, cutflow
        else:
            return jets, passing_jets


    def select_tightPU_jets(self, jets, muons, electrons, cutflow=None):
        #set_trace()
            ## pt and eta cuts
        pass_pt_eta_cuts = IDJet.make_pt_eta_cuts(jets)
        if cutflow is not None: cutflow["jets pass pT and eta cuts"] += ak.sum(pass_pt_eta_cuts)
    
            ## tightID
        # check https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVUL for definitions
        jet_ID = (jets.isTight) # pass tight ID
        if cutflow is not None: cutflow["jets pass ID"] += ak.sum(jet_ID)
    
            ## remove jets that don"t pass ID and pt/eta cuts
        jets = jets[(jet_ID & pass_pt_eta_cuts)]
    
            ## clean jets wrt veto+loose+tight el and mu
        jets_ak = ak.with_name(jets[["pt", "eta", "phi", "mass"]],"PtEtaPhiMLorentzVector")
        muons_ak = ak.with_name(muons[(muons["TIGHTMU"] | muons["LOOSEMU"] | muons["VETOMU"]) == True][["pt", "eta", "phi", "mass"]],"PtEtaPhiMLorentzVector")
        electrons_ak = ak.with_name(electrons[(electrons["TIGHTEL"] | electrons["LOOSEEL"] | electrons["VETOEL"]) == True][["pt", "eta", "phi", "mass"]],"PtEtaPhiMLorentzVector")
    
                # all mu
        jets_akc, muons_akc = ak.unzip(ak.cartesian([jets_ak, muons_ak], nested=True))
        clean_j_Mu_mask = ak.all((jets_akc.delta_r(muons_akc) >= 0.4), axis=2)
                # all el
        jets_akc, electrons_akc = ak.unzip(ak.cartesian([jets_ak, electrons_ak], nested=True))
        clean_j_El_mask = ak.all((jets_akc.delta_r(electrons_akc) >= 0.4), axis=2)
        clean_jets_mask = (clean_j_Mu_mask & clean_j_El_mask)
        if ak.sum(clean_jets_mask) > 0: jets = jets[clean_jets_mask]

        #set_trace()
            # apply tight PU cut to jets with Pt < 50
        lowPt_PU_jets_mask = (jets.pt < 50.) & (jets.puId == 7)
        hiPt_jets_mask = jets.pt >= 50.
        jets = jets[(lowPt_PU_jets_mask | hiPt_jets_mask)]
    
            ## leading jet pt cut
        leadpt_cut = IDJet.make_leadjet_pt_cut(jets)
        if cutflow is not None: cutflow["jets pass lead jet pT cut"] += ak.sum(leadpt_cut)
    
            ## 3 or more jets
        njet_restriction = 3
        njets_cuts = (ak.num(jets) >= njet_restriction)
        if cutflow is not None: cutflow["nEvts with %s+ clean jets passing ID and kin selection" % njet_restriction] += ak.sum(njets_cuts)
    
        passing_jets = (leadpt_cut & njets_cuts)
        if cutflow is not None: cutflow["nEvts with jets passing selection"] += ak.sum(passing_jets)
    
        if cutflow is not None:
            return jets, passing_jets, cutflow
        else:
            return jets, passing_jets




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
        selection = {jID: PackedSelection() for jID in self.jetID_types}

            # get all passing leptons
        lep_and_filter_pass = objsel.select_leptons(events, year=args.year, hem_15_16=apply_hem)
        {selection[jetID].add("lep_and_filter_pass", lep_and_filter_pass) for jetID in selection.keys()} # add passing leptons requirement to all systematics

            ## build corrected jets and MET
        events["CorrectedJets"], events["CorrectedMET"] = IDJet.process_jets(events, args.year, self.corrections["JetCor"])
        jets_to_use = events["CorrectedJets"]
        if (args.year == "2018") and (apply_hem):
            jets_to_use = objsel.remove_HEM_objs(obj=jets_to_use, isData=events.run if events.metadata["dataset"].startswith("data_Single") else None)
        #set_trace()

        if isData_:
            runs = events.run
            lumis = events.luminosityBlock
            Golden_Json_LumiMask = lumi_tools.LumiMask(lumiMask_path)
            LumiMask = Golden_Json_LumiMask.__call__(runs, lumis) ## returns array of valid events
            {selection[jetID].add("lumimask", LumiMask) for jetID in selection.keys()}
                  ## object selection and add different selections
            if isSM_Data_:
                        ## muons
                {selection[jetID].add("tight_MU", ak.sum(events["Muon"]["TIGHTMU"], axis=1) == 1) for jetID in selection.keys()} # one muon passing TIGHT criteria
            if isSE_Data_:
                        ## electrons
                {selection[jetID].add("tight_EL", ak.sum(events["Electron"]["TIGHTEL"], axis=1) == 1) for jetID in selection.keys()} # one electron passing TIGHT criteria

        else:
            ## add different selections
                    ## muons
            tight_mu_sel, loose_mu_sel = ak.sum(events["Muon"]["TIGHTMU"], axis=1) == 1, ak.sum(events["Muon"]["LOOSEMU"], axis=1) == 1
            {selection[jetID].add("tight_MU", tight_mu_sel) for jetID in selection.keys()} # one muon passing TIGHT criteria
                    ## electrons
            tight_el_sel, loose_el_sel = ak.sum(events["Electron"]["TIGHTEL"], axis=1) == 1, ak.sum(events["Electron"]["LOOSEEL"], axis=1) == 1
            {selection[jetID].add("tight_EL", tight_el_sel) for jetID in selection.keys()} # one electron passing TIGHT criteria

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
                {selection[jetID].add("semilep", ak.num(events["SL"]) > 0) for jetID in selection.keys()}
            else:
                {selection[jetID].add("semilep", np.zeros(len(events), dtype=bool)) for jetID in selection.keys()}

        #if to_debug: set_trace()
            # run over jet ID types that require changes to event objects (jets+MET)
        for jetID in self.jetID_types:
            print(jetID)
            output[f"cutflow_{jetID}"]["lep_and_filter_pass"] += ak.sum(lep_and_filter_pass)

                # jet selection
            selected_jets, passing_jets, output[f"cutflow_{jetID}"]  = self.select_tightID_jets(jets=jets_to_use, muons=events["Muon"], electrons=events["Electron"], cutflow=output[f"cutflow_{jetID}"]) if jetID == "tightID" \
                else self.select_tightPU_jets(jets=jets_to_use, muons=events["Muon"], electrons=events["Electron"], cutflow=output[f"cutflow_{jetID}"])
            output[f"cutflow_{jetID}"]["nEvts passing jet and lepton obj selection"] += ak.sum(passing_jets & lep_and_filter_pass)
            selection[jetID].add("passing_jets", passing_jets)
            selection[jetID].add("jets_3",  ak.num(selected_jets) == 3)
            selection[jetID].add("jets_4p",  ak.num(selected_jets) > 3) # only for getting btag weights
            selection[jetID].add(f"{btaggers[0]}_pass", ak.sum(selected_jets[btag_wp], axis=1) >= 2)

            # sort jets by btag value
            selected_jets = selected_jets[ak.argsort(selected_jets[IDJet.btag_tagger_to_disc_name[btaggers[0]]], ascending=False)]

            #set_trace()
                ## apply btagging SFs to MC
            if (not isData_) and (corrections["BTagSF"] == True):
                btag_weights = {key : np.ones(len(events)) for key in self.corrections["BTag_Constructors"][btaggers[0]]["3Jets"].schema_.keys()}

                threeJets_cut = selection[jetID].require(lep_and_filter_pass=True, passing_jets=True, jets_3=True)
                btagger_3j_wts = self.corrections["BTag_Constructors"][btaggers[0]]["3Jets"].get_scale_factor(jets=selected_jets[threeJets_cut], passing_cut=btag_wp)
                fourplusJets_cut = selection[jetID].require(lep_and_filter_pass=True, passing_jets=True, jets_4p=True)
                btagger_4pj_wts = self.corrections["BTag_Constructors"][btaggers[0]]["4PJets"].get_scale_factor(jets=selected_jets[fourplusJets_cut], passing_cut=btag_wp)

                    # fll dict of btag weights
                for wt_name in btagger_3j_wts.keys():
                    btag_weights[wt_name][threeJets_cut] = ak.prod(btagger_3j_wts[wt_name], axis=1)
                    btag_weights[wt_name][fourplusJets_cut] = ak.prod(btagger_4pj_wts[wt_name], axis=1)
    
            elif (not isData_) and (corrections["BTagSF"] == False):
                btag_weights = {"central" : np.ones(len(events))}
                print("BTag SFs not applied to MC")
                #raise ValueError("BTag SFs not applied to MC")

            #set_trace()
            ## fill hists for each region
            for lepton in self.regions[jetID].keys():
                evt_weights = mu_evt_weights if lepton == "Muon" else el_evt_weights
                if not isData_:
                    lep_SFs = muSFs_dict if lepton == "Muon" else elSFs_dict
                for btagregion in self.regions[jetID][lepton].keys():
                    for jmult in self.regions[jetID][lepton][btagregion].keys():
                        #set_trace()
                        cut = selection[jetID].all(*self.regions[jetID][lepton][btagregion][jmult])
                        output[f"cutflow_{jetID}"]["nEvts %s" % ", ".join([lepton, btagregion, jmult])] += cut.sum()

                        if to_debug: print(lepton, btagregion, jmult)
                        if cut.sum() > 0:
                            ltype = "MU" if lepton == "Muon" else "EL"
                            if f"loose_or_tight_{ltype}" in self.regions[jetID][lepton][btagregion][jmult]:
                                leptons = events[lepton][cut][((events[lepton][cut][f"TIGHT{ltype}"] == True) | (events[lepton][cut][f"LOOSE{ltype}"] == True))]
                            elif f"tight_{ltype}" in self.regions[jetID][lepton][btagregion][jmult]:
                                leptons = events[lepton][cut][(events[lepton][cut][f"TIGHT{ltype}"] == True)]
                            elif f"loose_{ltype}" in self.regions[jetID][lepton][btagregion][jmult]:
                                leptons = events[lepton][cut][(events[lepton][cut][f"LOOSE{ltype}"] == True)]
                            else:
                                raise ValueError("Not sure what lepton type to choose for event")

                            #set_trace()
                                # get jets and MET
                            jets, met = selected_jets[cut], events["CorrectedMET"][cut]

                                # find best permutations
                            best_perms = ttpermutator.find_best_permutations(jets=jets, leptons=leptons, MET=met, btagWP=btag_wp, btag_req=True if btagregion == "btagPass" else False)
                            valid_perms = ak.num(best_perms["TTbar"].pt) > 0
                            output[f"cutflow_{jetID}"]["nEvts %s: valid perms" % ", ".join([lepton, btagregion, jmult])] += ak.sum(valid_perms)

                            bp_status = np.zeros(cut.size, dtype=int) # 0 == "" (no gen matching), 1 == "right", 2 == "matchable", 3 == "unmatchable", 4 == "sl_tau", 5 == "noslep"
                                # get matched permutation (semilep ttbar only)
                            if isTTbar_:
                                semilep_evts = selection[jetID].require(semilep=True)
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
                            output[f"cutflow_{jetID}"]["nEvts %s: pass MT cut" % ", ".join([lepton, btagregion, jmult])] += ak.sum(MTHigh)


                                # fill hists 
                            wts = evt_weights.weight()[cut][valid_perms][MTHigh] if isData_ else (evt_weights.weight() * btag_weights["central"] * lep_SFs["central"])[cut][valid_perms][MTHigh]
                            sysname = events.metadata["dataset"].split("_")[-1] if isTTShift_ else "nosys"
                                # fill hists for interference samples
                            if isInt_:
                                #set_trace()
                                    # fill hists for positive weights
                                pos_evts = np.where(wts > 0)
                                self.sample_name = "%s_pos" % events.metadata["dataset"]
                                output = self.fill_hists(acc=output, jetmult=jmult, leptype=lepton, jetID_type=jetID, permarray=bp_status[cut][valid_perms][MTHigh][pos_evts],
                                    perm=best_perms[valid_perms][MTHigh][pos_evts], jets=jets[valid_perms][MTHigh][pos_evts], leptons=leptons[valid_perms][MTHigh][pos_evts], MTvals=MT[valid_perms][MTHigh][pos_evts], evt_wts=wts[pos_evts])
                                    # fill hists for negative weights
                                neg_evts = np.where(wts < 0)
                                self.sample_name = "%s_neg" % events.metadata["dataset"]
                                output = self.fill_hists(acc=output, jetmult=jmult, leptype=lepton, jetID_type=jetID, permarray=bp_status[cut][valid_perms][MTHigh][neg_evts],
                                    perm=best_perms[valid_perms][MTHigh][neg_evts], jets=jets[valid_perms][MTHigh][neg_evts], leptons=leptons[valid_perms][MTHigh][neg_evts], MTvals=MT[valid_perms][MTHigh][neg_evts], evt_wts=wts[neg_evts])
                            else:
                                output = self.fill_hists(acc=output, jetmult=jmult, leptype=lepton, jetID_type=jetID, permarray=bp_status[cut][valid_perms][MTHigh],
                                    perm=best_perms[valid_perms][MTHigh], jets=jets[valid_perms][MTHigh], leptons=leptons[valid_perms][MTHigh], MTvals=MT[valid_perms][MTHigh], evt_wts=wts)

        return output

    def fill_hists(self, acc, jetmult, leptype, jetID_type, permarray, perm, jets, leptons, MTvals, evt_wts):
            ## apply alpha correction for 3Jets
        if (jetmult == "3Jets") and ("Alpha" in self.corrections):
            alpha_corr = self.corrections["Alpha"](172.5/perm["THad"].mass)
            perm["THad"] = perm["THad"].multiply(alpha_corr) # correct thad
            perm["TTbar"] = ak.flatten(perm["THad"]+perm["TLep"]) # correct ttbar

        thad_ctstar, tlep_ctstar = make_vars.ctstar(perm["THad"], perm["TLep"], flatten=True)

        pt_sorted_jets = jets[ak.argsort(jets.pt, ascending=False)]

        #set_trace()
        for permval in np.unique(permarray).tolist():
            perm_inds = np.where(permarray == permval)
            dataset_name = "%s_%s" % (self.sample_name, perm_cats[permval]) if permval != 0 else self.sample_name

            acc["mtt_vs_tlep_ctstar_abs"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, jetID=jetID_type,
                mtt=ak.flatten(perm["TTbar"].mass)[perm_inds], ctstar_abs=np.abs(tlep_ctstar[perm_inds]), weight=evt_wts[perm_inds])

                # only fill 2D signal hist for systematics to save space and time
            #if sys != "nosys": continue
            acc["mtt"].fill(     dataset=dataset_name, jmult=jetmult, leptype=leptype, jetID=jetID_type, mtt=ak.flatten(perm["TTbar"].mass)[perm_inds], weight=evt_wts[perm_inds])
            acc["mthad"].fill(   dataset=dataset_name, jmult=jetmult, leptype=leptype, jetID=jetID_type, mtop=ak.flatten(perm["THad"].mass)[perm_inds], weight=evt_wts[perm_inds])
            acc["mWHad"].fill(   dataset=dataset_name, jmult=jetmult, leptype=leptype, jetID=jetID_type, wmass=ak.flatten(perm["WHad"].mass)[perm_inds], weight=evt_wts[perm_inds])
            acc["mWLep"].fill(   dataset=dataset_name, jmult=jetmult, leptype=leptype, jetID=jetID_type, wmass=ak.flatten(perm["WLep"].mass)[perm_inds], weight=evt_wts[perm_inds])
            acc["pt_thad"].fill( dataset=dataset_name, jmult=jetmult, leptype=leptype, jetID=jetID_type, pt=ak.flatten(perm["THad"].pt)[perm_inds], weight=evt_wts[perm_inds])
            acc["pt_tlep"].fill( dataset=dataset_name, jmult=jetmult, leptype=leptype, jetID=jetID_type, pt=ak.flatten(perm["TLep"].pt)[perm_inds], weight=evt_wts[perm_inds])
            acc["pt_tt"].fill(   dataset=dataset_name, jmult=jetmult, leptype=leptype, jetID=jetID_type, pt=ak.flatten(perm["TTbar"].pt)[perm_inds], weight=evt_wts[perm_inds])
            acc["eta_thad"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, jetID=jetID_type, eta=ak.flatten(perm["THad"].eta)[perm_inds], weight=evt_wts[perm_inds])
            acc["eta_tlep"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, jetID=jetID_type, eta=ak.flatten(perm["TLep"].eta)[perm_inds], weight=evt_wts[perm_inds])
            acc["eta_tt"].fill(  dataset=dataset_name, jmult=jetmult, leptype=leptype, jetID=jetID_type, eta=ak.flatten(perm["TTbar"].eta)[perm_inds], weight=evt_wts[perm_inds])

            acc["tlep_ctstar"].fill(    dataset=dataset_name, jmult=jetmult, leptype=leptype, jetID=jetID_type, ctstar=tlep_ctstar[perm_inds], weight=evt_wts[perm_inds])
            acc["tlep_ctstar_abs"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, jetID=jetID_type, ctstar_abs=np.abs(tlep_ctstar[perm_inds]), weight=evt_wts[perm_inds])

            acc["full_disc"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, jetID=jetID_type, prob=ak.flatten(perm["Prob"])[perm_inds], weight=evt_wts[perm_inds])
            acc["mass_disc"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, jetID=jetID_type, massdisc=ak.flatten(perm["MassDiscr"])[perm_inds], weight=evt_wts[perm_inds])
            acc["ns_disc"].fill(  dataset=dataset_name, jmult=jetmult, leptype=leptype, jetID=jetID_type, nsdisc=ak.flatten(perm["NuDiscr"])[perm_inds], weight=evt_wts[perm_inds])
            acc["ns_dist"].fill(  dataset=dataset_name, jmult=jetmult, leptype=leptype, jetID=jetID_type, nu_dist=ak.flatten(np.sqrt(perm["Nu"].chi2))[perm_inds], weight=evt_wts[perm_inds])

            acc["MET_pt"].fill( dataset=dataset_name, jmult=jetmult, leptype=leptype, jetID=jetID_type, pt=ak.flatten(perm["MET"].pt)[perm_inds], weight=evt_wts[perm_inds])
            acc["MET_phi"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, jetID=jetID_type, phi=ak.flatten(perm["MET"].phi)[perm_inds], weight=evt_wts[perm_inds])

            acc["Jets_pt"].fill(   dataset=dataset_name, jmult=jetmult, leptype=leptype, jetID=jetID_type, pt=ak.flatten(jets.pt[perm_inds]), weight=ak.flatten((ak.ones_like(jets.pt)*evt_wts)[perm_inds]))
            acc["Jets_eta"].fill(  dataset=dataset_name, jmult=jetmult, leptype=leptype, jetID=jetID_type, eta=ak.flatten(jets.eta[perm_inds]), weight=ak.flatten((ak.ones_like(jets.eta)*evt_wts)[perm_inds]))
            acc["Jets_phi"].fill(  dataset=dataset_name, jmult=jetmult, leptype=leptype, jetID=jetID_type, phi=ak.flatten(jets.phi[perm_inds]), weight=ak.flatten((ak.ones_like(jets.phi)*evt_wts)[perm_inds]))
            acc["Jets_njets"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, jetID=jetID_type, njets=ak.num(jets)[perm_inds], weight=evt_wts[perm_inds])
            acc["Jets_phi_vs_eta"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, jetID=jetID_type,
                phi_2d=ak.flatten(jets.phi[perm_inds]), eta_2d=ak.flatten(jets.eta[perm_inds]), weight=ak.flatten((ak.ones_like(jets.phi)*evt_wts)[perm_inds]))

            acc["Jets_LeadJet_pt"].fill( dataset=dataset_name, jmult=jetmult, leptype=leptype, jetID=jetID_type, pt=pt_sorted_jets[perm_inds].pt[:, 0], weight=evt_wts[perm_inds])
            acc["Jets_LeadJet_eta"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, jetID=jetID_type, eta=pt_sorted_jets[perm_inds].eta[:, 0], weight=evt_wts[perm_inds])
            acc[f"Jets_{btaggers[0]}_bDisc"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, jetID=jetID_type,
                bdisc=ak.flatten(jets[IDJet.btag_tagger_to_disc_name[btaggers[0]]][perm_inds]), weight=ak.flatten((ak.ones_like(jets.pt)*evt_wts)[perm_inds]))

            acc["Lep_pt"].fill(    dataset=dataset_name, jmult=jetmult, leptype=leptype, jetID=jetID_type, pt=ak.flatten(leptons.pt[perm_inds]), weight=evt_wts[perm_inds])
            acc["Lep_eta"].fill(   dataset=dataset_name, jmult=jetmult, leptype=leptype, jetID=jetID_type, eta=ak.flatten(leptons.eta[perm_inds]), weight=evt_wts[perm_inds])
            acc["Lep_phi"].fill(   dataset=dataset_name, jmult=jetmult, leptype=leptype, jetID=jetID_type, phi=ak.flatten(leptons.phi[perm_inds]), weight=evt_wts[perm_inds])
            acc["Lep_iso"].fill(   dataset=dataset_name, jmult=jetmult, leptype=leptype, jetID=jetID_type,
                iso=ak.flatten(leptons["pfRelIso04_all"][perm_inds]) if leptype == "Muon" else ak.flatten(leptons["pfRelIso03_all"][perm_inds]), weight=evt_wts[perm_inds])
            acc["Lep_phi_vs_eta"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, jetID=jetID_type,
                phi_2d=ak.flatten(leptons.phi[perm_inds]), eta_2d=ak.flatten(leptons.eta[perm_inds]), weight=evt_wts[perm_inds])

            acc["MT"].fill(dataset=dataset_name, jmult=jetmult, leptype=leptype, jetID=jetID_type, mt=ak.flatten(MTvals)[perm_inds], weight=evt_wts[perm_inds])

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
