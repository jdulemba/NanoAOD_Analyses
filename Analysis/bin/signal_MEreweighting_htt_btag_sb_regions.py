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
import Utilities.make_variables as make_vars
import python.TTPermutator as ttpermutator
import Utilities.systematics as systematics
import python.IDJet as IDJet
from copy import deepcopy
import Utilities.final_analysis_binning as final_binning

import uproot
from heavy_higgs_reweighting.util import Scenario

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
isSignal_ = (samplename.startswith("AtoTT") or samplename.startswith("HtoTT"))
isInt_ = isSignal_ and ("Int" in samplename)
if not isSignal_:
    raise ValueError(f"Only signal datasets can be used as input, {samplename} not allowed")

# convert input string of options dictionary to actual dictionary
odict = (args.opts).replace("\'", "\"")
opts_dict = prettyjson.loads(odict)

    ## set config options passed through argparse
import ast
to_debug = ast.literal_eval(opts_dict.get("debug", "False"))
apply_hem = ast.literal_eval(opts_dict.get("apply_hem", "True"))

## init tt probs for likelihoods
ttpermutator.year_to_run(year=args.year)


######### init signal ME reweighter
def convert_scenarioTOsamplename(scenario):
    decay = "SL" if scenario.channel == "lj" else "DiLep"
    signame = f"{scenario.parity}toTTJets{decay}_M{scenario.mass}_W{str(scenario.width).replace('.', 'p')}_{scenario.part.capitalize()}"
    return signame


#set_trace()
if args.year == "2016APV":
    year_to_use = "2016pre"
elif args.year == "2016":
    year_to_use = "2016post"
else:
    year_to_use = args.year


possible_masses = [365, 380, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000]
possible_widths = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 8.0, 10.0, 13.0, 15.0, 18.0, 21.0, 25.0]
possible_masses = [str(mass) for mass in possible_masses]
possible_widths = [str(width) for width in possible_widths]

masses_to_run = opts_dict.get("allowed_masses", "All")
widths_to_run = opts_dict.get("allowed_widths", "All")

if masses_to_run == "All":
    allowed_masses = possible_masses
else:
    allowed_masses = masses_to_run.split(":")
    allowed_masses = [mass for mass in allowed_masses if mass in possible_masses]

if widths_to_run == "All":
    allowed_widths = possible_widths
else:
    allowed_widths = widths_to_run.split(":")
    allowed_widths = [width for width in allowed_widths if width in possible_widths]

#set_trace()
allowed_masses = ["M"+mass for mass in allowed_masses]
allowed_widths = ["W"+width.replace(".", "p") for width in allowed_widths]
print(f"Making target signal points for {allowed_masses} and {allowed_widths}")

eos_dir = f"root://eosuser.cern.ch//eos/user/l/ljeppe/HeavyHiggs/ME_reweighting_weights/UL{year_to_use}_lowwidth/" if (("Res" in samplename) and (float(samplename.split("_")[2][1:].replace("p", ".")) < 10.)) \
        else f"root://eosuser.cern.ch//eos/user/l/ljeppe/HeavyHiggs/ME_reweighting_weights/UL{year_to_use}/"
sigwts_fnames_dict = {fileset[samplename][idx] : os.path.join(eos_dir, "__".join(fileset[samplename][idx].split("/")[10:]) if "W9p0" in samplename else "__".join(fileset[samplename][idx].split("/")[7:])) for idx in range(len(fileset[samplename]))}
##########

        # copy fileset root files to local condor node if running on condor
if "isCondor" in opts_dict.keys():
    if ast.literal_eval(opts_dict["isCondor"]):
        import subprocess
        from subprocess import check_output, STDOUT
        sites_to_try = ["root://xrootd-cms.infn.it/", "root://cmsxrootd.fnal.gov/"]
        #sites_to_try = []
        n_retries = len(sites_to_try) + 1
        for idx, rfile in enumerate(fileset[samplename]):
            cp_success = False
            for cp_attempt in range(n_retries):
                if cp_success: continue
                cp_rfile = rfile if cp_attempt == 0 else "/".join([sites_to_try[cp_attempt-1], rfile.split("//")[-1]]) # replace whatever redirector is used to regional Bari one
                print(f"Attempt {cp_attempt+1} to copy {cp_rfile} to /tmp")
                try:
                    output = check_output(["xrdcp", "-f", f"{cp_rfile}", "/tmp"], timeout=1200, stderr=STDOUT)
                    print(output.decode("utf-8"))
                    cp_success = True
                    #if to_debug: set_trace()
                except:
                    cp_success = False
                    continue
            if not cp_success:
                raise ValueError(f"{cp_rfile} not found")
            fileset[samplename][idx] = f"/tmp/{rfile.split('/')[-1]}"
            sigwts_fnames_dict[f"/tmp/{rfile.split('/')[-1]}"] = sigwts_fnames_dict[rfile]
            del sigwts_fnames_dict[rfile]

            # copy ME reweighting root files to local condor node if running on condor
        for fset_file, MEwt_file in sigwts_fnames_dict.items():
            print(f"Attempt to copy {MEwt_file} to /tmp")
            cp_success = False
            try:
                output = check_output(["xrdcp", "-f", f"{MEwt_file}", "/tmp"], timeout=1200, stderr=STDOUT)
                print(output.decode("utf-8"))
                cp_success = True
                #set_trace()
            except:
                cp_success = False
                continue
            if not cp_success:
                raise ValueError(f"{MEwt_file} not copied")
            sigwts_fnames_dict[fset_file] = f"/tmp/{MEwt_file.split('/')[-1]}"


sigwts_rfiles_dict = {key : uproot.open(val, array_cache="500 MB")["weights"] for key, val in sigwts_fnames_dict.items()} # open ME reweighting files

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

# get systematics to run
evt_sys_to_run = opts_dict.get("evt_sys", "NONE").upper()
rewt_systs_to_run = opts_dict.get("rewt_sys", "NONE")
only_sys = ast.literal_eval(opts_dict.get("only_sys", "False"))
    
import fnmatch
if only_sys: # don't run 'nosys'
    if (evt_sys_to_run == "NONE") and (rewt_systs_to_run == "NONE"):
        raise ValueError("At least one systematic must be specified in order to run only on systematics!")

    event_systematics_to_run = ["nosys"] if (rewt_systs_to_run != "NONE") else []
    reweight_systematics_to_run = []

else:
    event_systematics_to_run = ["nosys"]
    reweight_systematics_to_run = ["nosys"]

#set_trace()
event_systematics_to_run += [systematics.event_sys_opts[args.year][name] for name in systematics.event_sys_opts[args.year].keys() if fnmatch.fnmatch(name, evt_sys_to_run)]
for rewt_sys_to_run in rewt_systs_to_run.split(","):
    reweight_systematics_to_run += [systematics.signal_reweight_opts[args.year][name] for name in systematics.signal_reweight_opts[args.year].keys() if fnmatch.fnmatch(name, rewt_sys_to_run)]

    ## check that systematics only related to ttbar events aren't used for non-ttbar events
reweight_systematics_to_run = [sys for sys in reweight_systematics_to_run if sys not in systematics.ttJets_sys.values()]
    
print("Running with event systematics:", *sorted(set(event_systematics_to_run).difference(set(["nosys"]))), sep=", ") if "nosys" in event_systematics_to_run else print("Running with event systematics:", *sorted(event_systematics_to_run), sep=", ")
print("\t\tand reweight systematics:", *sorted(reweight_systematics_to_run), sep=", ")
#set_trace()

# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class MyAnalyzer(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.sys_axis = hist.Cat("sys", "Systematic")
        self.jetmult_axis = hist.Cat("jmult", "nJets")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")
        self.mtt_axis = hist.Bin("mtt", "m($t\overline{t}$) [GeV]", np.around(final_binning.mtt_binning, decimals=0))
        self.ctstar_abs_axis = hist.Bin("ctstar_abs", "|cos($\\theta^{*}$)|", np.around(final_binning.ctstar_abs_binning, decimals=2))

            ## make dictionary of hists
        histo_dict = {}
                ## make selection plots
        selection_hists = self.make_selection_hists()
        histo_dict.update(selection_hists)        

            ## for signal ME reweighting
        self.weights_rfiles_dict = sigwts_rfiles_dict
            ##

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
                "3Jets"  : {"lep_and_filter_pass", "passing_jets", "jets_3" , "tight_MU", f"{btaggers[0]}_pass"},
                "4PJets" : {"lep_and_filter_pass", "passing_jets", "jets_4p", "tight_MU", f"{btaggers[0]}_pass"},
            },
            "Electron" : {
                "3Jets"  : {"lep_and_filter_pass", "passing_jets", "jets_3" , "tight_EL", f"{btaggers[0]}_pass"},
                "4PJets" : {"lep_and_filter_pass", "passing_jets", "jets_4p", "tight_EL", f"{btaggers[0]}_pass"},
            },
        }
        self.regions = {sys:deepcopy(base_regions) for sys in self.event_systematics_to_run}


    @property
    def accumulator(self):
        return self._accumulator


    def make_selection_hists(self):
        histo_dict = {}
        histo_dict["mtt_vs_tlep_ctstar_abs"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.mtt_axis, self.ctstar_abs_axis)

        return histo_dict


    def process(self, events):
        output = self.accumulator.identity()

            ## get ME event weights
        #if to_debug: set_trace()
        weightsfile = self.weights_rfiles_dict[events.metadata["filename"]]
        valid_signal_arrays = [val.name for val in weightsfile.values() if ((convert_scenarioTOsamplename(Scenario.fromstr(val.name)).split("_")[1] in allowed_masses) and (convert_scenarioTOsamplename(Scenario.fromstr(val.name)).split("_")[2] in allowed_widths))]
        ME_weights_akarray = weightsfile.arrays(valid_signal_arrays, entry_start=events.metadata["entrystart"], entry_stop=events.metadata["entrystop"])

        mu_evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections, isTTbar=False, isSignal=True)
        el_evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections, isTTbar=False, isSignal=True)

            ## initialize selections
        selection = {evt_sys: PackedSelection() for evt_sys in self.event_systematics_to_run}

            # get all passing leptons
        lep_and_filter_pass = objsel.select_leptons(events, year=args.year, hem_15_16=apply_hem)
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
        if "LeptonSF" in corrections.keys():
            #set_trace()
            tight_muons = events["Muon"][tight_mu_sel][(events["Muon"][tight_mu_sel]["TIGHTMU"] == True)]
            muSFs_dict =  MCWeights.get_lepton_sf(sf_dict=self.corrections["LeptonSF"]["Muons"],
                pt=ak.flatten(tight_muons["pt"]), eta=ak.flatten(tight_muons["eta"]), tight_lep_mask=tight_mu_sel, leptype="Muons")
    
            tight_electrons = events["Electron"][tight_el_sel][(events["Electron"][tight_el_sel]["TIGHTEL"] == True)]
            elSFs_dict = MCWeights.get_lepton_sf(sf_dict=self.corrections["LeptonSF"]["Electrons"],
                pt=ak.flatten(tight_electrons["pt"]), eta=ak.flatten(tight_electrons["etaSC"]), tight_lep_mask=tight_el_sel, leptype="Electrons")

        #if to_debug: set_trace()
            # run over systematics that require changes to event objects (jets+MET)
        for evt_sys in self.event_systematics_to_run:
            output[f"cutflow_{evt_sys}"]["lep_and_filter_pass"] += ak.sum(lep_and_filter_pass)
                # jet selection
            passing_jets = objsel.jets_selection(events, year=args.year, cutflow=output[f"cutflow_{evt_sys}"], shift=evt_sys, hem_15_16=apply_hem, met_factory=self.corrections["JetCor"]["MC"]["METFactory"])
            output[f"cutflow_{evt_sys}"]["nEvts passing jet and lepton obj selection"] += ak.sum(passing_jets & lep_and_filter_pass)
            selection[evt_sys].add("passing_jets", passing_jets)
            selection[evt_sys].add("jets_3",  ak.num(events["SelectedJets"]) == 3)
            selection[evt_sys].add("jets_4p",  ak.num(events["SelectedJets"]) > 3) # only for getting btag weights
            selection[evt_sys].add(f"{btaggers[0]}_pass", ak.sum(events["SelectedJets"][btag_wp], axis=1) >= 2)

                # sort jets by btag value
            events["SelectedJets"] = events["SelectedJets"][ak.argsort(events["SelectedJets"][IDJet.btag_tagger_to_disc_name[btaggers[0]]], ascending=False)]

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
    
            else:
                raise ValueError("BTag SFs not applied to MC")

            #set_trace()
            ## fill hists for each region
            for lepton in self.regions[evt_sys].keys():
                evt_weights = mu_evt_weights if lepton == "Muon" else el_evt_weights
                lep_SFs = muSFs_dict if lepton == "Muon" else elSFs_dict
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
                        best_perms = ttpermutator.find_best_permutations(jets=jets, leptons=leptons, MET=met, btagWP=btag_wp, btag_req=True)
                        valid_perms = ak.num(best_perms["TTbar"].pt) > 0
                        output[f"cutflow_{evt_sys}"]["nEvts %s: valid perms" % ", ".join([lepton, jmult])] += ak.sum(valid_perms)

                            ## create MT regions
                        MT = make_vars.MT(leptons, met)
                        MTHigh = ak.flatten(MT[valid_perms] >= MTcut)
                        output[f"cutflow_{evt_sys}"]["nEvts %s: pass MT cut" % ", ".join([lepton, jmult])] += ak.sum(MTHigh)

                        #set_trace()
                            # loop over each target signal mass+width value for each systematic
                        for idx, signame in enumerate(ME_weights_akarray.fields):
                            sigwt = ME_weights_akarray[signame]
                            signame = convert_scenarioTOsamplename(Scenario.fromstr(signame))
                                # make sure target signal is one of the allowed points
                            #set_trace()
                            if not ((signame.split("_")[1] in allowed_masses) and (signame.split("_")[2] in allowed_widths)): continue
                            if to_debug: print(f"\ttarget signal: {signame}")
                                # fill hists for each systematic
                            if to_debug: print(f"\t\tevt sys: {evt_sys}")
                            if evt_sys == "nosys":
                                for rewt_sys in self.reweight_systematics_to_run:
                                    if to_debug: print(f"\t\tsysname: {rewt_sys}")

                                    #set_trace()
                                    if rewt_sys == "nosys":
                                        wts = (evt_weights.weight() * btag_weights["central"] * lep_SFs["central"] * sigwt)[cut][valid_perms][MTHigh]
                                    elif rewt_sys.startswith("btag"):
                                        wts = (evt_weights.weight() * btag_weights[rewt_sys.split("btag_")[-1]] * lep_SFs["central"] * sigwt)[cut][valid_perms][MTHigh]
                                    elif rewt_sys.startswith("Lep"):
                                        #set_trace()
                                        if rewt_sys.split("_")[-1] in lep_SFs.keys():
                                            wts = (evt_weights.weight() * btag_weights["central"] * lep_SFs[rewt_sys.split("_")[-1]] * sigwt)[cut][valid_perms][MTHigh]
                                        else:
                                            print(f"{rewt_sys.split('_')[-1]} not found in {lepton} SF dict. Skipping")
                                            continue
                                    else:
                                        if rewt_sys not in evt_weights.variations:
                                            print(f"{rewt_sys} not option in event weights. Skipping")
                                            continue
                                        wts = (evt_weights.weight(rewt_sys) * btag_weights["central"] * lep_SFs["central"] * sigwt)[cut][valid_perms][MTHigh]

                                        # fill hists for interference samples
                                    if isInt_:
                                        #set_trace()
                                            # fill hists for positive weights
                                        pos_evts = np.where(wts > 0)
                                        self.sample_name = f"{signame}_pos"
                                        output = self.fill_hists(acc=output, sys=rewt_sys, jetmult=jmult, leptype=lepton, perm=best_perms[valid_perms][MTHigh][pos_evts], evt_wts=wts[pos_evts])
                                            # fill hists for negative weights
                                        neg_evts = np.where(wts < 0)
                                        self.sample_name = f"{signame}_neg"
                                        output = self.fill_hists(acc=output, sys=rewt_sys, jetmult=jmult, leptype=lepton, perm=best_perms[valid_perms][MTHigh][neg_evts], evt_wts=wts[neg_evts])
                                    else:
                                        self.sample_name = signame
                                        output = self.fill_hists(acc=output, sys=rewt_sys, jetmult=jmult, leptype=lepton, perm=best_perms[valid_perms][MTHigh], evt_wts=wts)

                            else:
                                if to_debug: print(f"\t\tsysname: {evt_sys}")
                                #if to_debug: set_trace()
                                wts = (evt_weights.weight() * btag_weights["central"] *lep_SFs["central"] * sigwt)[cut][valid_perms][MTHigh]

                                        # fill hists for interference samples
                                if isInt_:
                                    #set_trace()
                                        # fill hists for positive weights
                                    pos_evts = np.where(wts > 0)
                                    self.sample_name = f"{signame}_pos"
                                    output = self.fill_hists(acc=output, sys=evt_sys, jetmult=jmult, leptype=lepton, perm=best_perms[valid_perms][MTHigh][pos_evts], evt_wts=wts[pos_evts])
                                        # fill hists for negative weights
                                    neg_evts = np.where(wts < 0)
                                    self.sample_name = f"{signame}_neg"
                                    output = self.fill_hists(acc=output, sys=evt_sys, jetmult=jmult, leptype=lepton, perm=best_perms[valid_perms][MTHigh][neg_evts], evt_wts=wts[neg_evts])
                                else:
                                    self.sample_name = signame
                                    output = self.fill_hists(acc=output, sys=evt_sys, jetmult=jmult, leptype=lepton, perm=best_perms[valid_perms][MTHigh], evt_wts=wts)

        return output

    def fill_hists(self, acc, sys, jetmult, leptype, perm, evt_wts):
            ## apply alpha correction for 3Jets
        if (jetmult == "3Jets") and ("Alpha" in self.corrections):
            alpha_corr = self.corrections["Alpha"](172.5/perm["THad"].mass)
            perm["THad"] = perm["THad"].multiply(alpha_corr) # correct thad
            perm["TTbar"] = ak.flatten(perm["THad"]+perm["TLep"]) # correct ttbar

        thad_ctstar, tlep_ctstar = make_vars.ctstar(perm["THad"], perm["TLep"], flatten=True)

        acc["mtt_vs_tlep_ctstar_abs"].fill(dataset=self.sample_name, sys=sys,  jmult=jetmult, leptype=leptype,
            mtt=ak.flatten(perm["TTbar"].mass, axis=None), ctstar_abs=np.abs(tlep_ctstar), weight=evt_wts)

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
        os.system(f"rm {' '.join(sigwts_fnames_dict.values())}")

toc = time.time()
print("Total time: %.1f" % (toc - tic))
