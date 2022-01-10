#!/usr/bin/env python

import time
tic = time.time()

from coffea import hist, processor
import awkward as ak
from coffea.nanoevents.methods import vector
ak.behavior.update(vector.behavior)
from pdb import set_trace
import os
from coffea.util import save
import numpy as np
import Utilities.make_variables as make_vars
import Utilities.prettyjson as prettyjson

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

# convert input string of options dictionary to actual dictionary
odict = (args.opts).replace("\'", "\"")
opts_dict = prettyjson.loads(odict)

    ## set config options passed through argparse
import ast
to_debug = ast.literal_eval(opts_dict.get("debug", "False"))

if len(fileset.keys()) > 1:
    raise ValueError("Only one topology run at a time in order to determine which corrections and systematics to run")
samplename = list(fileset.keys())[0]

   ## specify ttJets samples
nominal_ttJets_ = ["ttJetsSL", "ttJetsHad", "ttJetsDiLep"]
isNominalTTJets_ = samplename in nominal_ttJets_
if not isNominalTTJets_:
    raise ValueError("Should only be run on ttbar datasets")

# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class Analyzer(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.mtt_axis = hist.Bin("mtt", "mtt [GeV]", np.around(np.linspace(200., 4000., 381), decimals=0))
        self.ctstar_axis = hist.Bin("ctstar", "cos theta*", np.around(np.linspace(-1., 1., 41), decimals=2))

            ## make dictionary of hists
        histo_dict = {}
        histo_dict["mtt_vs_top_ctstar"] = hist.Hist("mtt_vs_top_ctstar", self.dataset_axis, self.mtt_axis, self.ctstar_axis)

        self._accumulator = processor.dict_accumulator(histo_dict)
    
    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        output = self.accumulator.identity()

        genWeights = events["genWeight"]

        genparts = events["GenPart"]
        gen_tops = genparts[(genparts.hasFlags(["isLastCopy"])) & (genparts.pdgId == 6)]
        gen_tbars = genparts[(genparts.hasFlags(["isLastCopy"])) & (genparts.pdgId == -6)]

        top_ctstar, tbar_ctstar = make_vars.ctstar(gen_tops, gen_tbars)
        mtt = (gen_tops+gen_tbars).mass
            # fill hists
        output["mtt_vs_top_ctstar"].fill(dataset=samplename, mtt=ak.flatten(mtt, axis=None), ctstar=ak.flatten(top_ctstar, axis=None), weight=genWeights)

        return output

    def postprocess(self, accumulator):
        return accumulator


proc_executor = processor.iterative_executor if to_debug else processor.futures_executor
proc_exec_args = {"schema": processor.NanoAODSchema} if to_debug else {"schema": processor.NanoAODSchema, "workers": 8}
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

toc = time.time()
print("Total time: %.1f" % (toc - tic))
