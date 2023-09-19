#!/usr/bin/env python

import time
tic = time.time()

from coffea import hist, processor
from coffea.analysis_tools import PackedSelection
import awkward as ak
from coffea.nanoevents.methods import vector
ak.behavior.update(vector.behavior)

from pdb import set_trace
from coffea.util import save
import os
import numpy as np
import Utilities.prettyjson as prettyjson
import Utilities.make_variables as make_vars
import uproot

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
event_systematics_to_run = eta_variations.keys()
print("Running toponium systematics:", *eta_variations.keys(), sep=", ")
#set_trace()

# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class htt_btag_sb_regions(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.sys_axis = hist.Cat("sys", "Systematic")

        self.pt_axis = hist.Bin("pt", "p_{T} [GeV]", np.around(np.linspace(0., 1000., 101), decimals=0))
        self.eta_axis = hist.Bin("eta", r"$\eta$", np.around(np.linspace(-5., 5., 101), decimals=1))
        self.mtop_axis = hist.Bin("mtop", "m(top) [GeV]", np.around(np.linspace(0., 1200., 1201), decimals=0))
        self.mtt_axis = hist.Bin("mtt", "m($t\overline{t}$) [GeV]", np.around(np.arange(200., 2001., 1), decimals=0))
        #self.mtt_axis = hist.Bin("mtt", "m($t\overline{t}$) [GeV]", np.around(np.linspace(200., 2000., 361), decimals=0))
        self.ctstar_axis = hist.Bin("ctstar", "cos($\\theta^{*}$)", np.around(np.linspace(-1., 1., 41), decimals=2))
        self.ctstar_abs_axis = hist.Bin("ctstar_abs", "|cos($\\theta^{*}$)|", np.around(np.linspace(0., 1., 21), decimals=2))

            ## make dictionary of hists
        histo_dict = {}
                ## make selection plots
        selection_hists = self.make_selection_hists()
        histo_dict.update(selection_hists)

        self._accumulator = processor.dict_accumulator(histo_dict)

            ## for toponium reweighting
        self.weights_rfiles_dict = topwts_rfiles_dict
            ##

        self.sample_name = ""
        self.event_systematics_to_run = event_systematics_to_run

    @property
    def accumulator(self):
        return self._accumulator


    def make_selection_hists(self):
        histo_dict = {}
        histo_dict["mtt"]      = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.mtt_axis)
        histo_dict["mtop"]    = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.mtop_axis)
        histo_dict["pt_top"]   = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.pt_axis)
        histo_dict["pt_tt"]    = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.pt_axis)
        histo_dict["eta_top"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.eta_axis)
        histo_dict["eta_tt"]   = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.eta_axis)

        histo_dict["top_ctstar"]     = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.ctstar_axis)
        histo_dict["top_ctstar_abs"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.ctstar_abs_axis)

        histo_dict["mtt_vs_top_ctstar_abs"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.mtt_axis, self.ctstar_abs_axis)

        return histo_dict


    def process(self, events):
        output = self.accumulator.identity()

        self.sample_name = events.metadata["dataset"]

            ## initialize selections
        selection = {evt_sys: PackedSelection() for evt_sys in self.event_systematics_to_run}
            # get gen level objects and weights
        genweight = np.copy(events["genWeight"])

        genparts = events["GenPart"]
        all_gen_tops = genparts[(genparts.hasFlags(["isLastCopy"])) & (genparts.pdgId == 6)]
        all_gen_tbars = genparts[(genparts.hasFlags(["isLastCopy"])) & (genparts.pdgId == -6)]
        {selection[sys].add("GenTops_Exist", ak.num(all_gen_tops) > 0) for sys in selection.keys()}
        {selection[sys].add("GenTbars_Exist", ak.num(all_gen_tbars) > 0) for sys in selection.keys()}


            ## make LHE mass cute for eta_t at +/- 6 GeV
        lhepart = events["LHEPart"][events["LHEPart"]["status"] == 1]
        mWWbb = (lhepart[:,0] + lhepart[:,1] + lhepart[:,2] + lhepart[:,3] + lhepart[:,4] + lhepart[:,5]).mass
        {selection[sys].add("Eta_Mass_Window", (mWWbb >= eta_variations[sys][0]) & (mWWbb <= eta_variations[sys][1])) for sys in selection.keys()}
            ## get event weights for the eta_tt variatons
        weightsfile = self.weights_rfiles_dict[events.metadata["filename"]]
        weights_akarray = weightsfile.arrays([val.name for val in weightsfile.values()], entry_start=events.metadata["entrystart"], entry_stop=events.metadata["entrystop"])

        for evt_sys in self.event_systematics_to_run:
            cut = selection[evt_sys].all(*["Eta_Mass_Window", "GenTops_Exist", "GenTbars_Exist"])
            if to_debug: print(f"\t{evt_sys}")
            gen_tops, gen_tbars = all_gen_tops[cut], all_gen_tbars[cut]
            gen_ttbars = gen_tops+gen_tbars
                # calculate ctstar
            top_ctstar, tbar_ctstar = make_vars.ctstar(gen_tops, gen_tbars, flatten=True)

            evt_wt = np.copy(genweight * weights_akarray[eta_variations[evt_sys][2]])[cut]

            output["mtt"].fill(dataset=self.sample_name, sys=evt_sys, mtt=ak.flatten(gen_ttbars.mass, axis=None), weight=evt_wt)
            output["mtop"].fill(dataset=self.sample_name, sys=evt_sys, mtop=ak.flatten(gen_tops.mass, axis=None), weight=evt_wt)
            output["pt_top"].fill(dataset=self.sample_name, sys=evt_sys, pt=ak.flatten(gen_tops.pt, axis=None), weight=evt_wt)
            output["pt_tt"].fill(dataset=self.sample_name, sys=evt_sys, pt=ak.flatten(gen_ttbars.pt, axis=None), weight=evt_wt)
            output["eta_top"].fill(dataset=self.sample_name, sys=evt_sys, eta=ak.flatten(gen_tops.eta, axis=None), weight=evt_wt)
            output["eta_tt"].fill(dataset=self.sample_name, sys=evt_sys, eta=ak.flatten(gen_ttbars.eta, axis=None), weight=evt_wt)
            output["top_ctstar"].fill(dataset=self.sample_name, sys=evt_sys, ctstar=top_ctstar, weight=evt_wt)
            output["top_ctstar_abs"].fill(dataset=self.sample_name, sys=evt_sys, ctstar_abs=np.abs(top_ctstar), weight=evt_wt)
            output["mtt_vs_top_ctstar_abs"].fill(dataset=self.sample_name, sys=evt_sys, mtt=ak.flatten(gen_ttbars.mass, axis=None), ctstar_abs=np.abs(top_ctstar), weight=evt_wt)

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
