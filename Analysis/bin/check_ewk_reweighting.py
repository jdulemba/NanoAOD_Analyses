#!/usr/bin/env python

import time
tic = time.time()

from coffea import hist, processor
from coffea.analysis_tools import PackedSelection
import awkward as ak
from coffea.nanoevents.methods import vector
ak.behavior.update(vector.behavior)

from pdb import set_trace
import os
from coffea.util import load, save
import numpy as np
import Utilities.make_variables as make_vars
import Utilities.prettyjson as prettyjson
import python.MCWeights as MCWeights

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

   ## specify ttJets samples
nominal_ttJets_ = ["ttJetsSL", "ttJetsHad", "ttJetsDiLep"]
isNominalTTJets_ = samplename in nominal_ttJets_
if not isNominalTTJets_:
    raise ValueError("Should only be run on ttbar datasets")

# convert input string of options dictionary to actual dictionary
odict = (args.opts).replace("\'", "\"")
opts_dict = prettyjson.loads(odict)

    ## set config options passed through argparse
import ast
to_debug = ast.literal_eval(opts_dict.get("debug", "False"))

cfg_pars = prettyjson.loads(open(os.path.join(proj_dir, "cfg_files", f"cfg_pars_{jobid}.json")).read())
ewk_reweighting = load(os.path.join(proj_dir, "Corrections", base_jobid, cfg_pars["corrections"]["ewk"]))

# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class Analyzer(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.reweighting_axis = hist.Cat("rewt", "Reweighting Type")
        self.pt_axis = hist.Bin("pt", "p_{T} [GeV]",  np.around(np.linspace(0., 1000., 201), decimals=0))
        self.eta_axis = hist.Bin("eta", r"$\eta$", np.around(np.linspace(-5., 5., 201), decimals=2))
        self.phi_axis = hist.Bin("phi", r"$\phi$", np.around(np.linspace(-4., 4., 161), decimals=2))
        self.mtop_axis = hist.Bin("mtop", "m(top) [GeV]", np.around(np.linspace(0., 300., 301), decimals=0))
        self.mtt_axis = hist.Bin("mtt", "m($t\overline{t}$) [GeV]", np.around(np.linspace(200., 4000., 761), decimals=0))
        self.ctstar_axis = hist.Bin("ctstar", "cos($\\theta^{*}$)", np.around(np.linspace(-1., 1., 201), decimals=2))
        self.ctstar_abs_axis = hist.Bin("ctstar_abs", "|cos($\\theta^{*}$)|", np.around(np.linspace(0., 1., 201), decimals=3))

            ## make dictionary of hists
        histo_dict = {}
            # top
        histo_dict["pt_top"] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.pt_axis)
        histo_dict["eta_top"] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.eta_axis)
        histo_dict["phi_top"] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.phi_axis)
        histo_dict["mass_top"] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.mtop_axis)
        histo_dict["ctstar_top"] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.ctstar_axis)
        histo_dict["ctstar_abs_top"] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.ctstar_abs_axis)
        histo_dict["cpTP_top"] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.ctstar_axis)

            # tbar
        histo_dict["pt_tbar"] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.pt_axis)
        histo_dict["eta_tbar"] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.eta_axis)
        histo_dict["phi_tbar"] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.phi_axis)
        histo_dict["mass_tbar"] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.mtop_axis)
        histo_dict["ctstar_tbar"] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.ctstar_axis)
        histo_dict["ctstar_abs_tbar"] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.ctstar_abs_axis)
        histo_dict["cpTP_tbar"] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.ctstar_axis)

            # ttbar
        histo_dict["pt_tt"] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.pt_axis)
        histo_dict["eta_tt"] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.eta_axis)
        histo_dict["phi_tt"] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.phi_axis)
        histo_dict["mtt"] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.mtt_axis)

        histo_dict["mtt_vs_top_ctstar"] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.mtt_axis, self.ctstar_axis)
        histo_dict["mtt_vs_top_ctstar_abs"] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.mtt_axis, self.ctstar_abs_axis)

        self.reweighting = ewk_reweighting

        self._accumulator = processor.dict_accumulator(histo_dict)
        self.sample_name = ""
    
    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        output = self.accumulator.identity()

        self.sample_name = events.metadata["dataset"]

        ewk_wts_dict = MCWeights.get_ewk_weights(ewk_reweighting, events)

        genparts = events["GenPart"]
        gen_tops = genparts[(genparts.hasFlags(["isLastCopy"])) & (genparts.pdgId == 6)]
        gen_tbars = genparts[(genparts.hasFlags(["isLastCopy"])) & (genparts.pdgId == -6)]
        gen_ttbars = gen_tops+gen_tbars

            # calculate ctstar
        top_ctstar, tbar_ctstar = make_vars.ctstar(gen_tops, gen_tbars)
        top_ctstar, tbar_ctstar = ak.flatten(top_ctstar, axis=None), ak.flatten(tbar_ctstar, axis=None)

            # calculate cpTP
        top_cpTP, tbar_cpTP = make_vars.cpTP(gen_tops, gen_tbars)
        top_cpTP, tbar_cpTP = ak.flatten(top_cpTP, axis=None), ak.flatten(tbar_cpTP, axis=None)

        for rewt_type in ["Nominal", "yt_0", "yt_0.5", "yt_1", "yt_1.5", "yt_2"]:
            evt_wt = np.ones(len(events)) if rewt_type == "Nominal" else ewk_wts_dict[rewt_type]

                    # top
            output["pt_top"].fill(dataset=self.sample_name, rewt=rewt_type, pt=ak.flatten(gen_tops.pt, axis=None), weight=evt_wt)
            output["eta_top"].fill(dataset=self.sample_name, rewt=rewt_type, eta=ak.flatten(gen_tops.eta, axis=None), weight=evt_wt)
            output["phi_top"].fill(dataset=self.sample_name, rewt=rewt_type, phi=ak.flatten(gen_tops.phi, axis=None), weight=evt_wt)
            output["mass_top"].fill(dataset=self.sample_name, rewt=rewt_type, mtop=ak.flatten(gen_tops.mass, axis=None), weight=evt_wt)
            output["ctstar_top"].fill(dataset=self.sample_name, rewt=rewt_type, ctstar=top_ctstar, weight=evt_wt)
            output["ctstar_abs_top"].fill(dataset=self.sample_name, rewt=rewt_type, ctstar_abs=np.abs(top_ctstar), weight=evt_wt)
            output["cpTP_top"].fill(dataset=self.sample_name, rewt=rewt_type, ctstar=top_cpTP, weight=evt_wt)

                    # tbar
            output["pt_tbar"].fill(dataset=self.sample_name, rewt=rewt_type, pt=ak.flatten(gen_tbars.pt, axis=None), weight=evt_wt)
            output["eta_tbar"].fill(dataset=self.sample_name, rewt=rewt_type, eta=ak.flatten(gen_tbars.eta, axis=None), weight=evt_wt)
            output["phi_tbar"].fill(dataset=self.sample_name, rewt=rewt_type, phi=ak.flatten(gen_tbars.phi, axis=None), weight=evt_wt)
            output["mass_tbar"].fill(dataset=self.sample_name, rewt=rewt_type, mtop=ak.flatten(gen_tbars.mass, axis=None), weight=evt_wt)
            output["ctstar_tbar"].fill(dataset=self.sample_name, rewt=rewt_type, ctstar=tbar_ctstar, weight=evt_wt)
            output["ctstar_abs_tbar"].fill(dataset=self.sample_name, rewt=rewt_type, ctstar_abs=np.abs(tbar_ctstar), weight=evt_wt)
            output["cpTP_tbar"].fill(dataset=self.sample_name, rewt=rewt_type, ctstar=tbar_cpTP, weight=evt_wt)

                    # ttbar system
            output["pt_tt"].fill(dataset=self.sample_name, rewt=rewt_type, pt=ak.flatten(gen_ttbars.pt, axis=None), weight=evt_wt)
            output["eta_tt"].fill(dataset=self.sample_name, rewt=rewt_type, eta=ak.flatten(gen_ttbars.eta, axis=None), weight=evt_wt)
            output["phi_tt"].fill(dataset=self.sample_name, rewt=rewt_type, phi=ak.flatten(gen_ttbars.phi, axis=None), weight=evt_wt)
            output["mtt"].fill(dataset=self.sample_name, rewt=rewt_type, mtt=ak.flatten(gen_ttbars.mass, axis=None), weight=evt_wt)
            output["mtt_vs_top_ctstar"].fill(dataset=self.sample_name, rewt=rewt_type, mtt=ak.flatten(gen_ttbars.mass, axis=None), ctstar=top_ctstar, weight=evt_wt)
            output["mtt_vs_top_ctstar_abs"].fill(dataset=self.sample_name, rewt=rewt_type, mtt=ak.flatten(gen_ttbars.mass, axis=None), ctstar_abs=np.abs(top_ctstar), weight=evt_wt)

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
print("%s has been written" % args.outfname)

toc = time.time()
print("Total time: %.1f" % (toc - tic))