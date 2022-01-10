#!/usr/bin/env python

import time
tic = time.time()

import awkward as ak
from coffea import hist, processor
from pdb import set_trace
import os
from coffea.util import save
import Utilities.prettyjson as prettyjson
import numpy as np
import coffea.lumi_tools.lumi_tools as lumi_tools

proj_dir = os.environ["PROJECT_DIR"]
base_jobid = os.environ["base_jobid"]

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

# convert input string of options dictionary to actual dictionary
odict = (args.opts).replace("\'", "\"")
opts_dict = prettyjson.loads(odict)

    ## set config options passed through argparse
import ast
to_debug = ast.literal_eval(opts_dict.get("debug", "False"))

    ## specify ttJets samples
Nominal_ttJets = ["ttJets_PS"] if ((args.year == "2016") and (base_jobid == "NanoAODv6")) else ["ttJetsSL", "ttJetsHad", "ttJetsDiLep"]
isTTbar_ = samplename.startswith("ttJets")
isTTSL_ = samplename.startswith("ttJetsSL")
isNominal_ttJets_ = samplename in Nominal_ttJets
if not isNominal_ttJets_:
    raise ValueError("Only nominal ttbar events should be run")

ewk_reweighting = load(os.path.join(proj_dir, "Corrections", base_jobid, "EWK_to_QCD_Ratios.coffea"))


# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class EWK_Analyzer(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.pu_nTrueInt_axis = hist.Bin("pu_nTrueInt", "nTrueInt", 100, 0, 100)
        self.pu_nPU_axis = hist.Bin("pu_nPU", "nPU", 100, 0, 100)

            ## make dictionary of hists
        histo_dict = {}
        histo_dict["PU_nTrueInt"] = hist.Hist("PU_nTrueInt", self.dataset_axis, self.pu_nTrueInt_axis)
        histo_dict["PU_nPU"] = hist.Hist("PU_nPU", self.dataset_axis, self.pu_nPU_axis)

        self._accumulator = processor.dict_accumulator(histo_dict)
        self.sample_name = ""
    
    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        output = self.accumulator.identity()

        event_nums = events.event
        self.sample_name = events.metadata["dataset"]

        lheparts = events["LHEPart"]
        incoming_lhe = lheparts[lheparts.status == -1]

        gg_events = ak.all(incoming_lhe.pdgId == 21, axis=1) # events whose initial state partons are both gluons
        qqu_events = ak.all(abs(incoming_lhe.pdgId) < 7, axis=1) & ak.all(np.mod(incoming_lhe.pdgId, 2) == 0, axis=1) # events whose intial state partons are both up-type quarks
        qqd_events = ak.all(abs(incoming_lhe.pdgId) < 7, axis=1) & ak.all(np.mod(incoming_lhe.pdgId, 2) == 1, axis=1) # events whose intial state partons are both down-type quarks
        qg_events = (ak.sum(incoming_lhe.pdgId == 21, axis=1) == 1) & (ak.sum(abs(incoming_lhe.pdgId) < 7, axis=1) == 1) # events that have one gluon and one quark as inital state partons

            # classify qg events into qqu or qqd
        outgoing_lhe = lheparts[lheparts.status == 1]
        qg_qqu_events = np.mod(abs(outgoing_lhe[qg_events].pdgId[:, 0]), 2) == 0 # qg events whose extra outgoing parton is up-type quark
        qg_qqd_events = np.mod(abs(outgoing_lhe[qg_events].pdgId[:, 0]), 2) == 1 # qg events whose extra outgoing parton is down-type quark

        
        set_trace()
        #if to_debug: set_trace()
        genWeights = events["genWeight"]

        output[self.sample_name]["sumGenWeights"] += sum(genWeights)

        output["PU_nTrueInt"].fill(dataset=self.sample_name, pu_nTrueInt=events["Pileup"]["nTrueInt"], weight=genWeights)
        output["PU_nPU"].fill(dataset=self.sample_name, pu_nPU=events["Pileup"]["nPU"], weight=genWeights)

        output[self.sample_name]["nEvents"] += ak.size(event_nums)

        if "LHEPdfWeight" in events.fields:
            LHEpdfWeights = events["LHEPdfWeight"]
                ## get sum of each pdf weight over all events
            sumLHEpdfWeights = ak.sum(LHEpdfWeights, axis=0)
            output[self.sample_name]["sumLHEpdfWeights"] += ak.to_numpy(sumLHEpdfWeights)

        if "PSWeight" in events.fields:
            psweights = events["PSWeight"]
                ## get sum of each weight over all events
            sumPSweights = ak.sum(psweights, axis=0)
            output[self.sample_name]["sumPSWeights"] += ak.to_numpy(sumPSweights)

        if "LHEScaleWeight" in events.fields:
            lheweights = events["LHEScaleWeight"]
                ## get sum of each weight over all events
            sumLHEscaleWeights = ak.sum(lheweights, axis=0)
            output[self.sample_name]["sumLHEscaleWeights"] += ak.to_numpy(sumLHEscaleWeights)


        return output

    def postprocess(self, accumulator):
        return accumulator

proc_executor = processor.iterative_executor if to_debug else processor.futures_executor
proc_exec_args = {"schema": processor.NanoAODSchema} if to_debug else {"schema": processor.NanoAODSchema, "workers": 8}
output = processor.run_uproot_job(
    fileset,
    treename="Events",
    processor_instance=EWK_Analyzer(),
    executor=proc_executor,
    executor_args=proc_exec_args,
    chunksize=100000,
)

save(output, args.outfname)
print(f"{args.outfname} has been written")

toc = time.time()
print("Total time: %.1f" % (toc - tic))
