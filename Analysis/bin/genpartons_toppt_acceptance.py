#!/usr/bin/env python

import time
tic = time.time()

from coffea import hist, processor
from coffea.nanoevents import NanoAODSchema
import awkward as ak
from coffea.nanoevents.methods import vector
ak.behavior.update(vector.behavior)
from coffea.util import save
from pdb import set_trace
import os
import python.MCWeights as MCWeights
import numpy as np
import Utilities.prettyjson as prettyjson
import python.GenParticleSelector as genpsel
import python.IDJet as IDJet
from coffea.analysis_tools import PackedSelection

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

if samplename != "ttJetsSL": raise ValueError("Only l+jets sample allowed")
isTTbar_ = ["ttJetsSL", "ttJetsHad", "ttJetsDiLep"]

# convert input string of options dictionary to actual dictionary
odict = (args.opts).replace("\'", "\"")
opts_dict = prettyjson.loads(odict)

    ## set config options passed through argparse
import ast
to_debug = ast.literal_eval(opts_dict.get("debug", "False"))

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
                    cp_success = True
                except:
                    cp_success = False
                    continue
            if not cp_success:
                raise ValueError(f"{cp_rfile} not found")
            fileset[samplename][idx] = f"/tmp/{rfile.split('/')[-1]}"


## load corrections for event weights
corrections = {
    "Prefire" : False,
}

reweight_systematics_to_run = ["nosys", "TopPtUp", "TopPtDown", "RENORMUp", "RENORMDown", "FACTORUp", "FACTORDown"]

# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class Analyzer(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.sys_axis = hist.Cat("sys", "Systematic")
        self.jetmult_axis = hist.Cat("jmult", "nJets")
        self.accept_axis = hist.Bin("acceptance", "", np.around(np.linspace(0., 2., 3), decimals=0))
        #self.accept_axis = hist.Bin("acceptance", "", np.around(np.linspace(0., 1., 2), decimals=0))

            ## make dictionary of hists
        histo_dict = {}
                ## make jet hists
        gen_hists = self.genpartons_hists()
        histo_dict.update(gen_hists)

        self._accumulator = processor.dict_accumulator(histo_dict)
        self.sample_name = ""
        self.corrections = corrections
        self.reweight_systematics_to_run = reweight_systematics_to_run

        self.lep_pt_min = 30
        self.lep_eta_max = 2.4

        self.regions = {
            "4Jets" : {
                "bhad_pass", "blep_pass", "wja_pass", "wjb_pass", "lepton_pass",
            },
            #"Electron" : {
            #    "lep_and_filter_pass", "passing_jets", "jets_3" , "tight_EL", "btag_pass", "semilep"
            #},
        }

    
    @property
    def accumulator(self):
        return self._accumulator


    def genpartons_hists(self):
        histo_dict = {}
        histo_dict["Acceptance"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.accept_axis)

        return histo_dict

    def process(self, events):
        output = self.accumulator.identity()

        self.sample_name = events.metadata["dataset"]

        #set_trace()
            ## make event weights
                # data or MC distinction made internally
        evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections, isTTbar=isTTbar_)

        genpsel.select(events, mode="NORMAL") # adds SL, DL, and Had objects to events

            # calculate top pT reweighting SFs
        toppt_weights = ak.flatten(MCWeights.get_TopPt_weights(events), axis=None)
        evt_weights.add("TopPt", np.ones(len(events)), np.copy(toppt_weights))

        selection = PackedSelection()

        all_bhad = events["SL"]["BHad"]
        passing_bhad_mask = IDJet.make_pt_eta_cuts(all_bhad)
        selection.add("bhad_pass", ak.flatten(passing_bhad_mask))

        all_blep = events["SL"]["BLep"]
        passing_blep_mask = IDJet.make_pt_eta_cuts(all_blep)
        selection.add("blep_pass", ak.flatten(passing_blep_mask))

        all_wja = events["SL"]["WJa"]
        passing_wja_mask = IDJet.make_pt_eta_cuts(all_wja)
        selection.add("wja_pass", ak.flatten(passing_wja_mask))

        all_wjb = events["SL"]["WJb"]
        passing_wjb_mask = IDJet.make_pt_eta_cuts(all_wjb)
        selection.add("wjb_pass", ak.flatten(passing_wjb_mask))

        all_leptons = events["SL"]["Lepton"]
        passing_lep_mask = (all_leptons["pt"] >= self.lep_pt_min) & (np.abs(all_leptons["eta"]) <= self.lep_eta_max)
        selection.add("lepton_pass", ak.flatten(passing_lep_mask))

        for njets in self.regions.keys():
            cut = selection.all(*self.regions[njets])

            #set_trace()        
            for rewt_sys in self.reweight_systematics_to_run:
                wts = evt_weights.weight() if rewt_sys == "nosys" else evt_weights.weight(rewt_sys)
                #output["Acceptance"].fill(dataset=self.sample_name, sys=rewt_sys, jmult=njets, acceptance=cut*1)
                output["Acceptance"].fill(dataset=self.sample_name, sys=rewt_sys, jmult=njets, acceptance=cut*1, weight=wts)
            #set_trace()        

        return output

    def postprocess(self, accumulator):
        return accumulator


proc_executor = processor.iterative_executor if to_debug else processor.futures_executor
proc_exec_args = {"schema": NanoAODSchema} if to_debug else {"schema": NanoAODSchema, "workers": 8}
output = processor.run_uproot_job(
    fileset,
    treename="Events",
    processor_instance=Analyzer(),
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
