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
#import python.GenParticleSelector as genpsel
import Utilities.make_variables as make_vars

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

isTTbar_ = ["ttJets_PS", "ttJets"] if ((args.year == "2016") and (base_jobid == "NanoAODv6")) else ["ttJetsSL", "ttJetsHad", "ttJetsDiLep"]
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
        self.flag_axis = hist.Cat("flag", "statusFlag")
        self.pt_axis = hist.Bin("pt", "p_{T} [GeV]", 200, 0, 1000)
        self.eta_axis = hist.Bin("eta", r"$\eta$", 60, -3, 3)
        self.phi_axis = hist.Bin("phi", r"$\phi$", 160, -4, 4)
        self.mass_axis = hist.Bin("mass", "mass [GeV]", 3000, 0, 3000)
        self.energy_axis = hist.Bin("energy", "E [GeV]", 600, 0, 3000)
        self.ctstar_axis = hist.Bin("ctstar", "cos($\\theta^{*}$)", 200, -1., 1.)
        self.ctstar_abs_axis = hist.Bin("ctstar_abs", "|cos($\\theta^{*}$)|", 200, 0., 1.)

            ## make dictionary of hists
        histo_dict = {}
                ## make jet hists
        gen_hists = self.gentops_comparison_hists()
        histo_dict.update(gen_hists)

        self._accumulator = processor.dict_accumulator(histo_dict)
        self.sample_name = ""
        self.corrections = corrections

        self.FLAGS = [
            #"isPrompt",
            #"isDecayedLeptonHadron",
            #"isTauDecayProduct",
            #"isPromptTauDecayProduct",
            #"isDirectTauDecayProduct",
            #"isDirectPromptTauDecayProduct",
            #"isDirectHadronDecayProduct",
            #"isHardProcess",
            #"fromHardProcess",
            #"isHardProcessTauDecayProduct",
            #"isDirectHardProcessTauDecayProduct",
            #"fromHardProcessBeforeFSR",
            "isFirstCopy",
            "isLastCopy",
            #"isLastCopyBeforeFSR",
        ]

    
    @property
    def accumulator(self):
        return self._accumulator


    def gentops_comparison_hists(self):
        histo_dict = {}
        histo_dict["pt"]    = hist.Hist("Events", self.dataset_axis, self.object_axis, self.flag_axis, self.pt_axis)
        histo_dict["eta"]   = hist.Hist("Events", self.dataset_axis, self.object_axis, self.flag_axis, self.eta_axis)
        histo_dict["phi"]   = hist.Hist("Events", self.dataset_axis, self.object_axis, self.flag_axis, self.phi_axis)
        histo_dict["mass"]  = hist.Hist("Events", self.dataset_axis, self.object_axis, self.flag_axis, self.mass_axis)
        histo_dict["energy"]= hist.Hist("Events", self.dataset_axis, self.object_axis, self.flag_axis, self.energy_axis)
        histo_dict["ctstar"]= hist.Hist("Events", self.dataset_axis, self.object_axis, self.flag_axis, self.ctstar_axis)
        histo_dict["ctstar_abs"]= hist.Hist("Events", self.dataset_axis, self.object_axis, self.flag_axis, self.ctstar_abs_axis)

        return histo_dict



    def process(self, events):
        output = self.accumulator.identity()

        self.sample_name = events.metadata["dataset"]

            ## make event weights
                # data or MC distinction made internally
        evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections)

        all_gen_tops = events["GenPart"][events["GenPart"].pdgId == 6]
        all_gen_tbars = events["GenPart"][events["GenPart"].pdgId == -6]

            # find all stages of decay/hadronization process of tops (and tbars) for which there are only at most 1 top passing per event (but not all events are empty)
        tops_dict = {flag : all_gen_tops[all_gen_tops.hasFlags([flag])] for flag in self.FLAGS if (ak.all(ak.num(all_gen_tops[all_gen_tops.hasFlags([flag])]) < 2) & ak.any(ak.num(all_gen_tops[all_gen_tops.hasFlags([flag])]) > 0))}
        tbars_dict = {flag : all_gen_tbars[all_gen_tbars.hasFlags([flag])] for flag in self.FLAGS if (ak.all(ak.num(all_gen_tbars[all_gen_tbars.hasFlags([flag])]) < 2) & ak.any(ak.num(all_gen_tbars[all_gen_tbars.hasFlags([flag])]) > 0))}

        for flag in tops_dict.keys():
            if not ak.all(ak.num(tops_dict[flag]) == ak.num(tbars_dict[flag])):
                print(f"Gen tops and tbars don't have same number of entries for {flag}")
                continue

            evt_mask = ak.num(tops_dict[flag]) > 0
            wts = evt_weights.weight()[evt_mask]

            gen_tops = tops_dict[flag][evt_mask]
            gen_tbars = tbars_dict[flag][evt_mask]
            gen_ttbars = gen_tops+gen_tbars

            #set_trace()
            top_ctstar, tbar_ctstar = make_vars.ctstar(gen_tops, gen_tbars)
            if "Int" in self.sample_name:
                pos_evts = np.where(wts > 0)
                    # tops
                output = self.fill_genp_hists(accumulator=output, dname="%s_pos" % self.sample_name, genp_type="Top", flag=flag, obj=gen_tops[pos_evts], evt_weights=wts[pos_evts])
                    # tbars
                output = self.fill_genp_hists(accumulator=output, dname="%s_pos" % self.sample_name, genp_type="Tbar", flag=flag, obj=gen_tbars[pos_evts], evt_weights=wts[pos_evts])
                    # ttbars
                output = self.fill_genp_hists(accumulator=output, dname="%s_pos" % self.sample_name, genp_type="TTbar", flag=flag, obj=gen_ttbars[pos_evts], evt_weights=wts[pos_evts])

                neg_evts = np.where(wts < 0)
                    # tops
                output = self.fill_genp_hists(accumulator=output, dname="%s_neg" % self.sample_name, genp_type="Top", flag=flag, obj=gen_tops[neg_evts], evt_weights=wts[neg_evts])
                    # tbars
                output = self.fill_genp_hists(accumulator=output, dname="%s_neg" % self.sample_name, genp_type="Tbar", flag=flag, obj=gen_tbars[neg_evts], evt_weights=wts[neg_evts])
                    # ttbars
                output = self.fill_genp_hists(accumulator=output, dname="%s_neg" % self.sample_name, genp_type="TTbar", flag=flag, obj=gen_ttbars[neg_evts], evt_weights=wts[neg_evts])

                        # fill ctstar hists
                output["ctstar"].fill(dataset="%s_pos" % self.sample_name, objtype="Top", flag=flag, ctstar=ak.flatten(top_ctstar[pos_evts], axis=None), weight=wts[pos_evts])
                output["ctstar"].fill(dataset="%s_pos" % self.sample_name, objtype="Tbar", flag=flag, ctstar=ak.flatten(tbar_ctstar[pos_evts], axis=None), weight=wts[pos_evts])
                output["ctstar"].fill(dataset="%s_neg" % self.sample_name, objtype="Top", flag=flag, ctstar=ak.flatten(top_ctstar[neg_evts], axis=None), weight=wts[neg_evts])
                output["ctstar"].fill(dataset="%s_neg" % self.sample_name, objtype="Tbar", flag=flag, ctstar=ak.flatten(tbar_ctstar[neg_evts], axis=None), weight=wts[neg_evts])
                output["ctstar_abs"].fill(dataset="%s_pos" % self.sample_name, objtype="Top", flag=flag, ctstar_abs=np.abs(ak.flatten(top_ctstar[pos_evts], axis=None)), weight=wts[pos_evts])
                output["ctstar_abs"].fill(dataset="%s_pos" % self.sample_name, objtype="Tbar", flag=flag, ctstar_abs=np.abs(ak.flatten(tbar_ctstar[pos_evts], axis=None)), weight=wts[pos_evts])
                output["ctstar_abs"].fill(dataset="%s_neg" % self.sample_name, objtype="Top", flag=flag, ctstar_abs=np.abs(ak.flatten(top_ctstar[neg_evts], axis=None)), weight=wts[neg_evts])
                output["ctstar_abs"].fill(dataset="%s_neg" % self.sample_name, objtype="Tbar", flag=flag, ctstar_abs=np.abs(ak.flatten(tbar_ctstar[neg_evts], axis=None)), weight=wts[neg_evts])
            else:
                    # tops
                output = self.fill_genp_hists(accumulator=output, dname=self.sample_name, genp_type="Top", flag=flag, obj=gen_tops, evt_weights=wts)
                    # tbars
                output = self.fill_genp_hists(accumulator=output, dname=self.sample_name, genp_type="Tbar", flag=flag, obj=gen_tbars, evt_weights=wts)
                    # ttbars
                output = self.fill_genp_hists(accumulator=output, dname=self.sample_name, genp_type="TTbar", flag=flag, obj=gen_ttbars, evt_weights=wts)
                
                        # fill ctstar hists
                output["ctstar"].fill(dataset=self.sample_name, objtype="Top", flag=flag, ctstar=ak.flatten(top_ctstar, axis=None), weight=wts)
                output["ctstar"].fill(dataset=self.sample_name, objtype="Tbar", flag=flag, ctstar=ak.flatten(tbar_ctstar, axis=None), weight=wts)
                output["ctstar_abs"].fill(dataset=self.sample_name, objtype="Top", flag=flag, ctstar_abs=np.abs(ak.flatten(top_ctstar, axis=None)), weight=wts)
                output["ctstar_abs"].fill(dataset=self.sample_name, objtype="Tbar", flag=flag, ctstar_abs=np.abs(ak.flatten(tbar_ctstar, axis=None)), weight=wts)

        return output

    def fill_genp_hists(self, accumulator, dname,  genp_type, flag, obj, evt_weights):
        #set_trace()
        accumulator["pt"].fill(    dataset=dname, objtype=genp_type, flag=flag, pt=ak.flatten(obj.pt, axis=None), weight=evt_weights)
        accumulator["eta"].fill(   dataset=dname, objtype=genp_type, flag=flag, eta=ak.flatten(obj.eta, axis=None), weight=evt_weights)
        accumulator["phi"].fill(   dataset=dname, objtype=genp_type, flag=flag, phi=ak.flatten(obj.phi, axis=None), weight=evt_weights)
        accumulator["mass"].fill(  dataset=dname, objtype=genp_type, flag=flag, mass=ak.flatten(obj.mass, axis=None), weight=evt_weights)
        accumulator["energy"].fill(dataset=dname, objtype=genp_type, flag=flag, energy=ak.flatten(obj.energy, axis=None), weight=evt_weights)

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
    chunksize=100000,
)

save(output, args.outfname)
print(f"{args.outfname} has been written")

toc = time.time()
print("Total time: %.1f" % (toc - tic))
