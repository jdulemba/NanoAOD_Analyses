#!/usr/bin/env python

"""
This analyzer compares the effect of applying the pileup corrections to simulation at gen-level
"""

import time
tic = time.time()

from coffea import hist, processor
import awkward as ak
from coffea.nanoevents.methods import vector
ak.behavior.update(vector.behavior)

from coffea.analysis_tools import Weights
from pdb import set_trace
from coffea.util import save, load
import os
import numpy as np
import Utilities.prettyjson as prettyjson

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("fset", type=str, help="Fileset dictionary (in string form) to be used for the processor")
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="Specify which year to run over")
parser.add_argument("outfname", type=str, help="Specify output filename, including directory and file extension")
parser.add_argument("opts", type=str, help="Fileset dictionary (in string form) to be used for the processor")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
base_jobid = os.environ["base_jobid"]

# convert input string of fileset dictionary to actual dictionary
fdict = (args.fset).replace("\'", "\"")
fileset = prettyjson.loads(fdict)

if len(fileset.keys()) > 1:
    raise ValueError("Only one topology run at a time in order to determine which corrections and systematics to run")
samplename = list(fileset.keys())[0]

# get dataset classification, used for corrections/systematics
isSignal_ = (samplename.startswith("AtoTT") or samplename.startswith("HtoTT"))
isData_ = samplename.startswith("data")
if (isData_ or isSignal_):
    raise ValueError("This should only be run with SM background simulation!")

# convert input string of options dictionary to actual dictionary
odict = (args.opts).replace("\'", "\"")
opts_dict = prettyjson.loads(odict)

    ## set config options passed through argparse
import ast
to_debug = ast.literal_eval(opts_dict.get("debug", "False"))


        # copy fileset root files to local condor node if running on condor
if "isCondor" in opts_dict.keys():
    if ast.literal_eval(opts_dict["isCondor"]):
        import subprocess
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
                    output = check_output(["xrdcp", "-f", f"{cp_rfile}", "/tmp"], timeout=600, stderr=STDOUT)
                    #output = check_output(["xrdcp", "-f", f"{cp_rfile}", "/tmp"], timeout=300, stderr=STDOUT)
                    cp_success = True
                except:
                    cp_success = False
                    continue
            if not cp_success:
                raise ValueError(f"{cp_rfile} not found")
            fileset[samplename][idx] = f"/tmp/{rfile.split('/')[-1]}"

## load corrections for event weights
cfg_pars = prettyjson.loads(open(os.path.join(proj_dir, "cfg_files", f"cfg_pars_{jobid}.json")).read())
corrections = {
    "Pileup" : load(os.path.join(proj_dir, "Corrections", base_jobid, cfg_pars["corrections"]["pu"]))[args.year],
}


# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class MyAnalyzer(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.sys_axis = hist.Cat("sys", "Systematic")

        self.nTI_axis = hist.Bin("nTI", "nTI", np.around(np.linspace(0., 100., 101), decimals=1))
        self.nPV_axis = hist.Bin("nPV", "nPV", np.around(np.linspace(0., 100., 101), decimals=1))
        self.rho_axis = hist.Bin("rho", "rho", np.around(np.linspace(0., 100., 1001), decimals=1))

            ## make dictionary of hists
        histo_dict = {}
        histo_dict["nTrueInt"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.nTI_axis)
        histo_dict["nPV"]      = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.nPV_axis)
        histo_dict["rho"]      = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.rho_axis)

        self.sample_name = ""
        self.corrections = corrections

        self._accumulator = processor.dict_accumulator(histo_dict)


    @property
    def accumulator(self):
        return self._accumulator


    def process(self, events):
        output = self.accumulator.identity()

        self.sample_name = events.metadata["dataset"]

        ## make event weights
        weights = Weights(len(events), storeIndividual=True)
            ## add gen weights
        weights.add("genweight", np.copy(events["genWeight"]))
            ## add pileup weights
        weights.add("Pileup",
            np.copy(corrections["Pileup"][events.metadata["dataset"]]["central"](events["Pileup"]["nTrueInt"])),
            np.copy(corrections["Pileup"][events.metadata["dataset"]]["up"](events["Pileup"]["nTrueInt"])),
            np.copy(corrections["Pileup"][events.metadata["dataset"]]["down"](events["Pileup"]["nTrueInt"]))
        )

        for sys in ["nosys", "Pileup", "PileupUp", "PileupDown"]:
            if sys == "nosys":
                evt_wt = np.copy(weights.partial_weight(["genweight"]))
            elif sys == "Pileup":
                evt_wt = np.copy(weights.weight())
            else:
                evt_wt = np.copy(weights.weight(sys))
            output["nTrueInt"].fill(dataset=self.sample_name, sys=sys, nTI=ak.flatten(events["Pileup"]["nTrueInt"], axis=None), weight=evt_wt)
            output["nPV"].fill(dataset=self.sample_name, sys=sys, nPV=ak.flatten(events["PV"]["npvs"], axis=None), weight=evt_wt)
            output["rho"].fill(dataset=self.sample_name, sys=sys, rho=ak.flatten(events["fixedGridRhoFastjetAll"], axis=None), weight=evt_wt)


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

toc = time.time()
print("Total time: %.1f" % (toc - tic))
