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

isTTbar_ = ["ttJetsSL", "ttJetsHad", "ttJetsDiLep"]
isSignal_ = (samplename.startswith("AtoTT") or samplename.startswith("HtoTT"))
isInt_ = isSignal_ and ("Int" in samplename)
if (samplename not in isTTbar_) and (not isSignal_):
    raise ValueError("This should only be run on SM ttbar or signal events!")

# convert input string of options dictionary to actual dictionary
odict = (args.opts).replace("\'", "\"")
opts_dict = prettyjson.loads(odict)

    ## set config options passed through argparse
import ast
to_debug = ast.literal_eval(opts_dict.get("debug", "False"))

## load corrections for event weights
corrections = {
    "Prefire" : False,
}


# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class Analyzer(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.object_axis = hist.Cat("objtype", "Gen Object")
        self.ttdecaymode_axis = hist.Cat("ttdecay", "tt Decay Mode")
        self.wdecaymode_axis = hist.Cat("wdecay", "W Decay Mode")
        self.pt_axis = hist.Bin("pt", "p_{T} [GeV]", 200, 0, 1000)
        self.eta_axis = hist.Bin("eta", r"$\eta$", 60, -3, 3)
        self.phi_axis = hist.Bin("phi", r"$\phi$", 160, -4, 4)
        self.mass_axis = hist.Bin("mass", "mass [GeV]", 3000, 0, 3000)
        self.energy_axis = hist.Bin("energy", "E [GeV]", 200, 0, 1000)

            ## make dictionary of hists
        histo_dict = {}
                ## make jet hists
        gen_hists = self.genpartons_hists()
        histo_dict.update(gen_hists)

        self._accumulator = processor.dict_accumulator(histo_dict)
        self.sample_name = ""
        self.corrections = corrections

    
    @property
    def accumulator(self):
        return self._accumulator


    def genpartons_hists(self):
        histo_dict = {}
        histo_dict["pt"]    = hist.Hist("Events", self.dataset_axis, self.object_axis, self.ttdecaymode_axis, self.pt_axis)
        histo_dict["eta"]   = hist.Hist("Events", self.dataset_axis, self.object_axis, self.ttdecaymode_axis, self.eta_axis)
        histo_dict["phi"]   = hist.Hist("Events", self.dataset_axis, self.object_axis, self.ttdecaymode_axis, self.phi_axis)
        histo_dict["mass"]  = hist.Hist("Events", self.dataset_axis, self.object_axis, self.ttdecaymode_axis, self.mass_axis)
        histo_dict["energy"]= hist.Hist("Events", self.dataset_axis, self.object_axis, self.ttdecaymode_axis, self.energy_axis)

        return histo_dict



    def process(self, events):
        output = self.accumulator.identity()

        self.sample_name = events.metadata["dataset"]

        #set_trace()
            ## make event weights
                # data or MC distinction made internally
        evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections)

        genpsel.select(events, mode="NORMAL") # adds SL, DL, and Had objects to events

        for ttdecay in ["SL", "DL", "Had"]:
            #set_trace()
            if to_debug: print(ttdecay)
            ttdecay_mask = ak.num(events[ttdecay]) > 0
            if ak.sum(ttdecay_mask) == 0: continue

            wts = evt_weights.weight()[ttdecay_mask]
            for gen_obj in events[ttdecay].fields:
                if to_debug: print(gen_obj)
                if "Int" in self.sample_name:
                    pos_evts = np.where(wts > 0)
                    output = self.fill_genp_hists(accumulator=output, dname="%s_pos" % self.sample_name, genp_type=gen_obj, ttdecaymode=ttdecay, obj=events[ttdecay][gen_obj][ttdecay_mask][pos_evts], evt_weights=wts[pos_evts])
                    neg_evts = np.where(wts < 0)
                    output = self.fill_genp_hists(accumulator=output, dname="%s_neg" % self.sample_name, genp_type=gen_obj, ttdecaymode=ttdecay, obj=events[ttdecay][gen_obj][ttdecay_mask][neg_evts], evt_weights=wts[neg_evts])
                else:
                    output = self.fill_genp_hists(accumulator=output, dname=self.sample_name, genp_type=gen_obj, ttdecaymode=ttdecay, obj=events[ttdecay][gen_obj][ttdecay_mask], evt_weights=wts)
                

        return output

    def fill_genp_hists(self, accumulator, dname,  genp_type, ttdecaymode, obj, evt_weights):
        #set_trace()
        accumulator["pt"].fill(    dataset=dname, objtype=genp_type, ttdecay=ttdecaymode, pt=ak.flatten(obj.pt, axis=None), weight=evt_weights)
        accumulator["eta"].fill(   dataset=dname, objtype=genp_type, ttdecay=ttdecaymode, eta=ak.flatten(obj.eta, axis=None), weight=evt_weights)
        accumulator["phi"].fill(   dataset=dname, objtype=genp_type, ttdecay=ttdecaymode, phi=ak.flatten(obj.phi, axis=None), weight=evt_weights)
        accumulator["mass"].fill(  dataset=dname, objtype=genp_type, ttdecay=ttdecaymode, mass=ak.flatten(obj.mass, axis=None), weight=evt_weights)
        accumulator["energy"].fill(dataset=dname, objtype=genp_type, ttdecay=ttdecaymode, energy=ak.flatten(obj.energy, axis=None), weight=evt_weights)

        return accumulator        


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
    #chunksize=10000 if to_debug else 100000,
    chunksize=100000,
)

save(output, args.outfname)
print(f"{args.outfname} has been written")

toc = time.time()
print("Total time: %.1f" % (toc - tic))
