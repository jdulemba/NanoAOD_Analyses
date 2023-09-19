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
import python.TTPermutator as ttpermutator
import Utilities.systematics as systematics
import python.IDJet as IDJet
from copy import deepcopy
import uproot
from coffea.analysis_tools import Weights

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
isToponium_ = samplename.startswith("Toponium")
if not isToponium_: raise ValueError(f"Only toponium samples are to be run, not {samplename}.")

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
### find files used for reweighting to systematic variations
if args.year == "2016APV":
    year_to_use = "16pre"
elif args.year == "2016":
    year_to_use = "16post"
elif args.year == "2017":
    year_to_use = "17"
elif args.year == "2018":
    year_to_use = "18"
else:
    raise ValueError(f"{args.year} not valid")

eos_dir = f"root://eosuser.cern.ch//eos/user/j/jdulemba/NanoAOD_Analyses/toponium_reweighting/ul{year_to_use}/"
decay = "EtaTTo1L1Nu2J_meta343_mtop172p5" if samplename == "ToponiumSL" else "EtaTTo2L2Nu_meta343_mtop172p5"
topwts_fnames_dict = {fileset[samplename][idx] : os.path.join(eos_dir, decay, fileset[samplename][idx].split("/")[-1]) for idx in range(len(fileset[samplename]))}
####

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
            topwts_fnames_dict[f"/tmp/{rfile.split('/')[-1]}"] = topwts_fnames_dict[rfile]
            del topwts_fnames_dict[rfile]

            # copy toponium reweighting root files to local condor node if running on condor
        for fset_file, wt_file in topwts_fnames_dict.items():
            print(f"Attempt to copy {wt_file} to /tmp")
            cp_success = False
            wt_outfname = f"/tmp/Wt_{wt_file.split('/')[-1]}"
            try:
                output = check_output(["xrdcp", "-f", f"{wt_file}", wt_outfname], timeout=1200, stderr=STDOUT)
                print(output.decode("utf-8"))
                cp_success = True
                #set_trace()
            except:
                cp_success = False
                continue
            if not cp_success:
                raise ValueError(f"{wt_file} not copied")
            topwts_fnames_dict[fset_file] = wt_outfname


topwts_rfiles_dict = {key : uproot.open(val, array_cache="500 MB")["weights"] for key, val in topwts_fnames_dict.items()} # open toponium reweighting files

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
base_jobid = os.environ["base_jobid"]

#set_trace()
cfg_pars = prettyjson.loads(open(os.path.join(proj_dir, "cfg_files", f"cfg_pars_{jobid}.json")).read())
## load corrections for event weights
corrections = {
    "Pileup" : load(os.path.join(proj_dir, "Corrections", base_jobid, cfg_pars["corrections"]["pu"]))[args.year],
    "Prefire" : True,
    "LeptonSF" : load(os.path.join(proj_dir, "Corrections", base_jobid, cfg_pars["corrections"]["lepton"]))[args.year],
    "BTagSF" : True,
    "JetCor" : load(os.path.join(proj_dir, "Corrections", base_jobid, cfg_pars["corrections"]["jetmet"][cfg_pars["corrections"]["jetmet"]["to_use"]]))[args.year],
    "Alpha" : load(os.path.join(proj_dir, "Corrections", base_jobid, cfg_pars["corrections"]["alpha"]))[args.year]["E"]["All_2D"],
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

    ## parameters for MT
MTcut = jet_pars["MT"]

    ## cuts specific to toponium mass window
nominal_eta_mass = 343
def eta_mass_min(eta_mass):
    return eta_mass - 6
def eta_mass_max(eta_mass):
    return eta_mass + 6

#set_trace()
eta_variations = {
    "nosys"             : [eta_mass_min(nominal_eta_mass),     eta_mass_max(nominal_eta_mass),     f"{decay.split('_')[0]}_meta343_mtop172p5"],
    "BindingEnergyUp"   : [eta_mass_min(nominal_eta_mass + 1), eta_mass_max(nominal_eta_mass + 1), f"{decay.split('_')[0]}_meta344_mtop172p5"],
    "BindingEnergyDown" : [eta_mass_min(nominal_eta_mass - 1), eta_mass_max(nominal_eta_mass - 1), f"{decay.split('_')[0]}_meta342_mtop172p5"],
    "TopMassUp"         : [eta_mass_min(nominal_eta_mass + 2), eta_mass_max(nominal_eta_mass + 2), f"{decay.split('_')[0]}_meta345_mtop173p5"],
    "TopMassDown"       : [eta_mass_min(nominal_eta_mass - 2), eta_mass_max(nominal_eta_mass - 2), f"{decay.split('_')[0]}_meta341_mtop171p5"],
}

# get systematics to run
#reweight_systematics_to_run = ["nosys"]
event_systematics_to_run = eta_variations.keys()
print("Running toponium systematics:", *eta_variations.keys(), sep=", ")
#set_trace()

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

            ## for toponium reweighting
        self.weights_rfiles_dict = topwts_rfiles_dict
        #self.weights_rfiles_dict = topwts_fnames_dict
            ##

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
        #self.reweight_systematics_to_run = reweight_systematics_to_run
        self.event_systematics_to_run = event_systematics_to_run

            ## make dict of cutflow for each systematic variation
        for sys in self.event_systematics_to_run:
            histo_dict[f"cutflow_{sys}"] = processor.defaultdict_accumulator(int)
    
        self._accumulator = processor.dict_accumulator(histo_dict)

            ## define regions depending on data or MC, systematic or nosys
        base_regions = {
            "Muon" : {
                "btagPass" : {
                    "3Jets"  : {"lep_and_filter_pass", "passing_jets", "jets_3" , "tight_MU", f"{btaggers[0]}_pass", "Eta_Mass_Window"},
                    "4PJets" : {"lep_and_filter_pass", "passing_jets", "jets_4p", "tight_MU", f"{btaggers[0]}_pass", "Eta_Mass_Window"},
                },
            },
            "Electron" : {
                "btagPass" : {
                    "3Jets"  : {"lep_and_filter_pass", "passing_jets", "jets_3" , "tight_EL", f"{btaggers[0]}_pass", "Eta_Mass_Window"},
                    "4PJets" : {"lep_and_filter_pass", "passing_jets", "jets_4p", "tight_EL", f"{btaggers[0]}_pass", "Eta_Mass_Window"},
                },
            },
        }
        self.regions = {sys:deepcopy(base_regions) for sys in self.event_systematics_to_run}
            # add sideband regions for nosys
        if "nosys" in self.event_systematics_to_run:
            self.regions["nosys"]["Muon"].update({
                key : {
                    "3Jets"  : {"lep_and_filter_pass", "passing_jets", "jets_3" , "tight_MU", f"{btaggers[0]}_{key}", "Eta_Mass_Window"},
                    "4PJets" : {"lep_and_filter_pass", "passing_jets", "jets_4p", "tight_MU", f"{btaggers[0]}_{key}", "Eta_Mass_Window"},
                } for key in btag_regions.keys()
            })
            self.regions["nosys"]["Electron"].update({
                key : {
                    "3Jets"  : {"lep_and_filter_pass", "passing_jets", "jets_3" , "tight_EL", f"{btaggers[0]}_{key}", "Eta_Mass_Window"},
                    "4PJets" : {"lep_and_filter_pass", "passing_jets", "jets_4p", "tight_EL", f"{btaggers[0]}_{key}", "Eta_Mass_Window"},
                } for key in btag_regions.keys()
            })


    @property
    def accumulator(self):
        return self._accumulator


    def make_jet_hists(self):
        histo_dict = {}
        histo_dict["Jets_pt"]    = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.pt_axis)
        histo_dict["Jets_eta"]   = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.eta_axis)
        histo_dict["Jets_phi"]   = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.phi_axis)
        histo_dict["Jets_njets"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.njets_axis)
        histo_dict["Jets_LeadJet_pt"]    = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.pt_axis)
        histo_dict["Jets_LeadJet_eta"]   = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.eta_axis)
        histo_dict["Jets_phi_vs_eta"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.phi_2d_axis, self.eta_2d_axis)
        histo_dict[f"Jets_{btaggers[0]}_bDisc"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.btag_discr_axis)

        return histo_dict

    def make_lep_hists(self):
        histo_dict = {}
        histo_dict["Lep_pt"]    = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.pt_axis)
        histo_dict["Lep_eta"]   = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.eta_axis)
        histo_dict["Lep_iso"]   = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepIso_axis)
        histo_dict["Lep_phi"]   = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.phi_axis)
        histo_dict["Lep_phi_vs_eta"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.phi_2d_axis, self.eta_2d_axis)
        ##histo_dict["Lep_etaSC"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.eta_axis)
        ##histo_dict["Lep_energy"]= hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.energy_axis)

        return histo_dict

    def make_selection_hists(self):
        histo_dict = {}
        histo_dict["mtt"]      = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.mtt_axis)
        histo_dict["mthad"]    = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.mtop_axis)
        histo_dict["mWHad"]    = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.wmass_axis)
        histo_dict["mWLep"]    = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.wmass_axis)
        histo_dict["pt_thad"]  = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.pt_axis)
        histo_dict["pt_tlep"]  = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.pt_axis)
        histo_dict["pt_tt"]    = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.pt_axis)
        histo_dict["eta_thad"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.eta_axis)
        histo_dict["eta_tlep"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.eta_axis)
        histo_dict["eta_tt"]   = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.eta_axis)

        histo_dict["tlep_ctstar"]     = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.ctstar_axis)
        histo_dict["tlep_ctstar_abs"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.ctstar_abs_axis)

        histo_dict["full_disc"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.probDisc_axis)
        histo_dict["mass_disc"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.massDisc_axis)
        histo_dict["ns_disc"]   = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.nsDisc_axis)
        histo_dict["ns_dist"]   = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.nu_dist_axis)

        histo_dict["MT"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.mt_axis)

        histo_dict["MET_pt"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.pt_axis)
        histo_dict["MET_phi"]= hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.phi_axis)

        histo_dict["mtt_vs_tlep_ctstar_abs"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.mtt_axis, self.ctstar_abs_axis)

        return histo_dict


    def process(self, events):
        output = self.accumulator.identity()

        self.sample_name = events.metadata["dataset"]

        #set_trace()
            ## make event weights
        evt_weights = Weights(len(events), storeIndividual=True)
        evt_weights.add("genweight", ak.copy(events["genWeight"])) ## Generator Weights
            ## Prefire Corrections
        if (args.year != "2018") and (self.corrections["Prefire"] == True) and ("L1PreFiringWeight" in events.fields):
            evt_weights.add("Prefire", ak.copy(events["L1PreFiringWeight"]["Nom"]))

            ## Pileup Reweighting
        if "Pileup" in self.corrections.keys():
            pu_wts_dict = {
                key : ak.copy(self.corrections["Pileup"][f"{events.metadata['dataset']}_{key}"]["central"](events["Pileup"]["nTrueInt"]) ) for key in eta_variations.keys()
            }

            ## initialize selections
        selection = {evt_sys: PackedSelection() for evt_sys in self.event_systematics_to_run}

        #if to_debug: set_trace()
            ## make LHE mass cute for eta_t at +/- 6 GeV
        lhepart = events["LHEPart"][events["LHEPart"]["status"] == 1]
        mWWbb = (lhepart[:,0] + lhepart[:,1] + lhepart[:,2] + lhepart[:,3] + lhepart[:,4] + lhepart[:,5]).mass
        {selection[sys].add("Eta_Mass_Window", (mWWbb >= eta_variations[sys][0]) & (mWWbb <= eta_variations[sys][1])) for sys in selection.keys()}
            ## get event weights for the eta_tt variatons
        weightsfile = self.weights_rfiles_dict[events.metadata["filename"]]
        weights_akarray = weightsfile.arrays([val.name for val in weightsfile.values()], entry_start=events.metadata["entrystart"], entry_stop=events.metadata["entrystop"])


            # get all passing leptons
        lep_and_filter_pass = objsel.select_leptons(events, year=args.year, hem_15_16=apply_hem)
        output["cutflow_nosys"]["lep_and_filter_pass"] += ak.sum(lep_and_filter_pass)
        {selection[sys].add("lep_and_filter_pass", lep_and_filter_pass) for sys in selection.keys()} # add passing leptons requirement to all systematics

            ## build corrected jets and MET
        events["CorrectedJets"], events["CorrectedMET"] = IDJet.process_jets(events, args.year, self.corrections["JetCor"])

        ## add different selections
                ## muons
        tight_mu_sel, loose_mu_sel = ak.sum(events["Muon"]["TIGHTMU"], axis=1) == 1, ak.sum(events["Muon"]["LOOSEMU"], axis=1) == 1
        {selection[sys].add("tight_MU", tight_mu_sel) for sys in selection.keys()} # one muon passing TIGHT criteria
                ## electrons
        tight_el_sel, loose_el_sel = ak.sum(events["Electron"]["TIGHTEL"], axis=1) == 1, ak.sum(events["Electron"]["LOOSEEL"], axis=1) == 1
        {selection[sys].add("tight_EL", tight_el_sel) for sys in selection.keys()} # one electron passing TIGHT criteria

        ### apply lepton SFs to MC (only applicable to tight leptons)
        if "LeptonSF" in self.corrections.keys():
            #set_trace()
            tight_muons = events["Muon"][tight_mu_sel][(events["Muon"][tight_mu_sel]["TIGHTMU"] == True)]
            muSFs_dict =  MCWeights.get_lepton_sf(sf_dict=self.corrections["LeptonSF"]["Muons"],
                pt=ak.flatten(tight_muons["pt"]), eta=ak.flatten(tight_muons["eta"]), tight_lep_mask=tight_mu_sel, leptype="Muons")
    
            tight_electrons = events["Electron"][tight_el_sel][(events["Electron"][tight_el_sel]["TIGHTEL"] == True)]
            elSFs_dict = MCWeights.get_lepton_sf(sf_dict=self.corrections["LeptonSF"]["Electrons"],
                pt=ak.flatten(tight_electrons["pt"]), eta=ak.flatten(tight_electrons["etaSC"]), tight_lep_mask=tight_el_sel, leptype="Electrons")


        ## get jet and MET selection - nosys selection used for systematic since eta_mass_window cuts don't affect them
            # jet selection
        passing_jets = objsel.jets_selection(events, year=args.year, cutflow=output["cutflow_nosys"], shift="nosys", hem_15_16=apply_hem, met_factory=self.corrections["JetCor"]["MC"]["METFactory"])
        output["cutflow_nosys"]["nEvts passing jet and lepton obj selection"] += ak.sum(passing_jets & lep_and_filter_pass)
        {selection[sys].add("passing_jets", passing_jets) for sys in selection.keys()}
        {selection[sys].add("jets_3", ak.num(events["SelectedJets"]) == 3) for sys in selection.keys()}
        {selection[sys].add("jets_4p", ak.num(events["SelectedJets"]) > 3) for sys in selection.keys()} # only for getting btag weights
        {selection[sys].add(f"{btaggers[0]}_pass", ak.sum(events["SelectedJets"][btag_wp], axis=1) >= 2) for sys in selection.keys()}

        # sort jets by btag value
        events["SelectedJets"] = events["SelectedJets"][ak.argsort(events["SelectedJets"][IDJet.btag_tagger_to_disc_name[btaggers[0]]], ascending=False)]
            # btag sidebands
        btagger_sorted = events["SelectedJets"][ak.argsort(events["SelectedJets"][IDJet.btag_tagger_to_disc_name[btaggers[0]]], ascending=False)][IDJet.btag_tagger_to_disc_name[btaggers[0]]]
        for btag_reg, (low_bound, up_bound) in btag_regions.items():
            {selection[sys].add(f"{btaggers[0]}_{btag_reg}", (ak.max(btagger_sorted, axis=1) > low_bound) & (ak.max(btagger_sorted, axis=1) <= up_bound)) for sys in selection.keys()}

            ## apply btagging SFs to MC
        if (self.corrections["BTagSF"] == True):
            btag_weights = {key : np.ones(len(events)) for key in self.corrections["BTag_Constructors"][btaggers[0]]["3Jets"].schema_.keys()}

            threeJets_cut = selection["nosys"].require(lep_and_filter_pass=True, passing_jets=True, jets_3=True)
            btagger_3j_wts = self.corrections["BTag_Constructors"][btaggers[0]]["3Jets"].get_scale_factor(jets=events["SelectedJets"][threeJets_cut], passing_cut=btag_wp)
            fourplusJets_cut = selection["nosys"].require(lep_and_filter_pass=True, passing_jets=True, jets_4p=True)
            btagger_4pj_wts = self.corrections["BTag_Constructors"][btaggers[0]]["4PJets"].get_scale_factor(jets=events["SelectedJets"][fourplusJets_cut], passing_cut=btag_wp)

                # fll dict of btag weights
            for wt_name in btagger_3j_wts.keys():
                btag_weights[wt_name][threeJets_cut] = ak.prod(btagger_3j_wts[wt_name], axis=1)
                btag_weights[wt_name][fourplusJets_cut] = ak.prod(btagger_4pj_wts[wt_name], axis=1)
    
        else:
            btag_weights = {"central" : np.ones(len(events))}
            print("BTag SFs not applied to MC")

            # run over systematics for eta_tt mass window
        for evt_sys in self.event_systematics_to_run:
            #set_trace()
            ## fill hists for each region
            for lepton in self.regions[evt_sys].keys():
                lep_SFs = muSFs_dict if lepton == "Muon" else elSFs_dict
                for btagregion in self.regions[evt_sys][lepton].keys():
                    for jmult in self.regions[evt_sys][lepton][btagregion].keys():
                            ## loop over different event selections based on mass windows for eta_t
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
                            #if to_debug: set_trace()
                            wts = np.copy(evt_weights.weight() * btag_weights["central"] *lep_SFs["central"] * pu_wts_dict[evt_sys] * weights_akarray[eta_variations[evt_sys][2]])[cut][valid_perms][MTHigh]

                                # fill hists
                            output = self.fill_hists(acc=output, sys=evt_sys, jetmult=jmult, leptype=lepton, btagregion=btagregion, permarray=bp_status[cut][valid_perms][MTHigh],
                                perm=best_perms[valid_perms][MTHigh], jets=jets[valid_perms][MTHigh], leptons=leptons[valid_perms][MTHigh], MTvals=MT[valid_perms][MTHigh], evt_wts=wts)

        return output

    def fill_hists(self, acc, sys, jetmult, leptype, btagregion, permarray, perm, jets, leptons, MTvals, evt_wts):
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

            acc["mtt_vs_tlep_ctstar_abs"].fill(dataset=dataset_name, sys=sys,  jmult=jetmult, leptype=leptype, btag=btagregion,
                mtt=ak.flatten(perm["TTbar"].mass)[perm_inds], ctstar_abs=np.abs(tlep_ctstar[perm_inds]), weight=evt_wts[perm_inds])

                # only fill 2D signal hist for systematics to save space and time
            #if sys != "nosys": continue

            acc["mtt"].fill(     dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, mtt=ak.flatten(perm["TTbar"].mass)[perm_inds], weight=evt_wts[perm_inds])
            acc["mthad"].fill(   dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, mtop=ak.flatten(perm["THad"].mass)[perm_inds], weight=evt_wts[perm_inds])
            acc["mWHad"].fill(   dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, wmass=ak.flatten(perm["WHad"].mass)[perm_inds], weight=evt_wts[perm_inds])
            acc["mWLep"].fill(   dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, wmass=ak.flatten(perm["WLep"].mass)[perm_inds], weight=evt_wts[perm_inds])
            acc["pt_thad"].fill( dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, pt=ak.flatten(perm["THad"].pt)[perm_inds], weight=evt_wts[perm_inds])
            acc["pt_tlep"].fill( dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, pt=ak.flatten(perm["TLep"].pt)[perm_inds], weight=evt_wts[perm_inds])
            acc["pt_tt"].fill(   dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, pt=ak.flatten(perm["TTbar"].pt)[perm_inds], weight=evt_wts[perm_inds])
            acc["eta_thad"].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, eta=ak.flatten(perm["THad"].eta)[perm_inds], weight=evt_wts[perm_inds])
            acc["eta_tlep"].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, eta=ak.flatten(perm["TLep"].eta)[perm_inds], weight=evt_wts[perm_inds])
            acc["eta_tt"].fill(  dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, eta=ak.flatten(perm["TTbar"].eta)[perm_inds], weight=evt_wts[perm_inds])

            acc["tlep_ctstar"].fill(    dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, ctstar=tlep_ctstar[perm_inds], weight=evt_wts[perm_inds])
            acc["tlep_ctstar_abs"].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, ctstar_abs=np.abs(tlep_ctstar[perm_inds]), weight=evt_wts[perm_inds])

            acc["full_disc"].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, prob=ak.flatten(perm["Prob"])[perm_inds], weight=evt_wts[perm_inds])
            acc["mass_disc"].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, massdisc=ak.flatten(perm["MassDiscr"])[perm_inds], weight=evt_wts[perm_inds])
            acc["ns_disc"].fill(  dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, nsdisc=ak.flatten(perm["NuDiscr"])[perm_inds], weight=evt_wts[perm_inds])
            acc["ns_dist"].fill(  dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, nu_dist=ak.flatten(np.sqrt(perm["Nu"].chi2))[perm_inds], weight=evt_wts[perm_inds])

            acc["MET_pt"].fill( dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, pt=ak.flatten(perm["MET"].pt)[perm_inds], weight=evt_wts[perm_inds])
            acc["MET_phi"].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, phi=ak.flatten(perm["MET"].phi)[perm_inds], weight=evt_wts[perm_inds])

            acc["Jets_pt"].fill(   dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, pt=ak.flatten(jets.pt[perm_inds]), weight=ak.flatten((ak.ones_like(jets.pt)*evt_wts)[perm_inds]))
            acc["Jets_eta"].fill(  dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, eta=ak.flatten(jets.eta[perm_inds]), weight=ak.flatten((ak.ones_like(jets.eta)*evt_wts)[perm_inds]))
            acc["Jets_phi"].fill(  dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, phi=ak.flatten(jets.phi[perm_inds]), weight=ak.flatten((ak.ones_like(jets.phi)*evt_wts)[perm_inds]))
            acc["Jets_njets"].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, njets=ak.num(jets)[perm_inds], weight=evt_wts[perm_inds])
            acc["Jets_phi_vs_eta"].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion,
                phi_2d=ak.flatten(jets.phi[perm_inds]), eta_2d=ak.flatten(jets.eta[perm_inds]), weight=ak.flatten((ak.ones_like(jets.phi)*evt_wts)[perm_inds]))

            acc["Jets_LeadJet_pt"].fill( dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, pt=pt_sorted_jets[perm_inds].pt[:, 0], weight=evt_wts[perm_inds])
            acc["Jets_LeadJet_eta"].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, eta=pt_sorted_jets[perm_inds].eta[:, 0], weight=evt_wts[perm_inds])
            acc[f"Jets_{btaggers[0]}_bDisc"].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion,
                bdisc=ak.flatten(jets[IDJet.btag_tagger_to_disc_name[btaggers[0]]][perm_inds]), weight=ak.flatten((ak.ones_like(jets.pt)*evt_wts)[perm_inds]))

            acc["Lep_pt"].fill(    dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, pt=ak.flatten(leptons.pt[perm_inds]), weight=evt_wts[perm_inds])
            acc["Lep_eta"].fill(   dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, eta=ak.flatten(leptons.eta[perm_inds]), weight=evt_wts[perm_inds])
            acc["Lep_phi"].fill(   dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, phi=ak.flatten(leptons.phi[perm_inds]), weight=evt_wts[perm_inds])
            acc["Lep_iso"].fill(   dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion,
                iso=ak.flatten(leptons["pfRelIso04_all"][perm_inds]) if leptype == "Muon" else ak.flatten(leptons["pfRelIso03_all"][perm_inds]), weight=evt_wts[perm_inds])
            acc["Lep_phi_vs_eta"].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion,
                phi_2d=ak.flatten(leptons.phi[perm_inds]), eta_2d=ak.flatten(leptons.eta[perm_inds]), weight=evt_wts[perm_inds])

            acc["MT"].fill(dataset=dataset_name, sys=sys, jmult=jetmult, leptype=leptype, btag=btagregion, mt=ak.flatten(MTvals)[perm_inds], weight=evt_wts[perm_inds])

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
        os.system(f"rm {' '.join(topwts_fnames_dict.values())}")

toc = time.time()
print("Total time: %.1f" % (toc - tic))
