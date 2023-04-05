#!/usr/bin/env python

import time
tic = time.time()

from coffea import hist, processor
import awkward as ak
from coffea.nanoevents.methods import vector
ak.behavior.update(vector.behavior)

from pdb import set_trace
from coffea.util import save, load
import os
import python.MCWeights as MCWeights
import numpy as np
import Utilities.prettyjson as prettyjson
import Utilities.make_variables as make_vars
import python.GenParticleSelector as genpsel

import uproot
from Utilities.heavy_higgs_reweighting_util import Scenario

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

#set_trace()
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


eos_dir = f"root://eosuser.cern.ch//eos/user/l/ljeppe/HeavyHiggs/ME_reweighting_weights/UL{year_to_use}/"
sigwts_fnames_dict = {fileset[samplename][idx] : os.path.join(eos_dir, "__".join(fileset[samplename][idx].split("/")[7:])) for idx in range(len(fileset[samplename]))}
try:
    sigwts_rfiles_dict = {key : uproot.open(val) for key, val in sigwts_fnames_dict.items()}
except:
    sigwts_rfiles_dict= {}
    for fname, rewt_fname in sigwts_fnames_dict.items():
        try:
            sigwts_rfiles_dict.update({fname : uproot.open(rewt_fname)})
        except:
            print(f"Unable to find reweighting file {rewt_fname}")
            fileset[samplename].remove(fname)
##########

#set_trace()
corrections = {
    "Prefire" : False,
}
# get systematics to run
reweight_systematics_to_run = ["nosys"]
#reweight_systematics_to_run = ["nosys", "AH_FACTORUp", "AH_FACTORDown", "AH_RENORMUp", "AH_RENORMDown"]
print("\n\nRunning with systematics:", *sorted(reweight_systematics_to_run))
#set_trace()


# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class MyAnalyzer(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.sys_axis = hist.Cat("sys", "Systematic")

        self.pt_axis = hist.Bin("pt", "p_{T} [GeV]", np.around(np.linspace(0., 1000., 101), decimals=0))
        self.eta_axis = hist.Bin("eta", r"$\eta$", np.around(np.linspace(-5., 5., 101), decimals=1))
        self.mtop_axis = hist.Bin("mtop", "m(top) [GeV]", np.around(np.linspace(0., 1200., 1201), decimals=0))
        self.mtt_axis = hist.Bin("mtt", "m($t\overline{t}$) [GeV]", np.around(np.linspace(200., 2000., 361), decimals=0))
        self.ctstar_axis = hist.Bin("ctstar", "cos($\\theta^{*}$)", np.around(np.linspace(-1., 1., 41), decimals=2))
        self.ctstar_abs_axis = hist.Bin("ctstar_abs", "|cos($\\theta^{*}$)|", np.around(np.linspace(0., 1., 21), decimals=2))

            ## make dictionary of hists
        histo_dict = {}
                ## make selection plots
        selection_hists = self.make_selection_hists()
        histo_dict.update(selection_hists)        

            ## for signal ME reweighting
        self.weights_rfiles_dict = sigwts_rfiles_dict
            ##

        self.corrections = corrections
        self.reweight_systematics_to_run = reweight_systematics_to_run

        self._accumulator = processor.dict_accumulator(histo_dict)

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

            ## make event weights
        #set_trace()
        weightsfile = self.weights_rfiles_dict[events.metadata["filename"]]
        ME_weights_akarray = weightsfile["weights"].arrays(entry_start=events.metadata["entrystart"], entry_stop=events.metadata["entrystop"])

        ## make event weights
        evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections, isTTbar=False, isSignal=True)

        genpsel.select(events, mode="NORMAL")

        #set_trace()
        for ttdecay in ["SL", "DL"]:
            if not ak.any(events[ttdecay]): continue
                # calculate ctstar
            top_ctstar, tbar_ctstar = make_vars.ctstar(events[ttdecay]["Top"], events[ttdecay]["Tbar"], flatten=True)

            #set_trace()
            for signame in ME_weights_akarray.fields:
                sigwt = ME_weights_akarray[signame]
                signame = convert_scenarioTOsamplename(Scenario.fromstr(signame))
                if to_debug: print(f"\ttarget signal: {signame}")
                for rewt_sys in self.reweight_systematics_to_run:
                    if to_debug: print(f"\t\t{rewt_sys}")

                    evt_wt = evt_weights.weight() * sigwt if rewt_sys == "nosys" else evt_weights.weight(rewt_sys) * sigwt
                    if isInt_:
                        pos_evts = evt_wt > 0
                        output["mtt"].fill(dataset=f"{signame}_pos", sys=rewt_sys, mtt=ak.flatten(events[ttdecay]["TTbar"].mass, axis=None)[pos_evts], weight=evt_wt[pos_evts])
                        output["mtop"].fill(dataset=f"{signame}_pos", sys=rewt_sys, mtop=ak.flatten(events[ttdecay]["Top"].mass, axis=None)[pos_evts], weight=evt_wt[pos_evts])
                        output["pt_top"].fill(dataset=f"{signame}_pos", sys=rewt_sys, pt=ak.flatten(events[ttdecay]["Top"].pt, axis=None)[pos_evts], weight=evt_wt[pos_evts])
                        output["pt_tt"].fill(dataset=f"{signame}_pos", sys=rewt_sys, pt=ak.flatten(events[ttdecay]["TTbar"].pt, axis=None)[pos_evts], weight=evt_wt[pos_evts])
                        output["eta_top"].fill(dataset=f"{signame}_pos", sys=rewt_sys, eta=ak.flatten(events[ttdecay]["Top"].eta, axis=None)[pos_evts], weight=evt_wt[pos_evts])
                        output["eta_tt"].fill(dataset=f"{signame}_pos", sys=rewt_sys, eta=ak.flatten(events[ttdecay]["TTbar"].eta, axis=None)[pos_evts], weight=evt_wt[pos_evts])
                        output["top_ctstar"].fill(dataset=f"{signame}_pos", sys=rewt_sys, ctstar=top_ctstar[pos_evts], weight=evt_wt[pos_evts])
                        output["top_ctstar_abs"].fill(dataset=f"{signame}_pos", sys=rewt_sys, ctstar_abs=np.abs(top_ctstar)[pos_evts], weight=evt_wt[pos_evts])
                        output["mtt_vs_top_ctstar_abs"].fill(dataset=f"{signame}_pos", sys=rewt_sys, mtt=ak.flatten(events[ttdecay]["TTbar"].mass, axis=None)[pos_evts], ctstar_abs=np.abs(top_ctstar)[pos_evts], weight=evt_wt[pos_evts])

                        neg_evts = evt_wt < 0
                        output["mtt"].fill(dataset=f"{signame}_neg", sys=rewt_sys, mtt=ak.flatten(events[ttdecay]["TTbar"].mass, axis=None)[neg_evts], weight=evt_wt[neg_evts])
                        output["mtop"].fill(dataset=f"{signame}_neg", sys=rewt_sys, mtop=ak.flatten(events[ttdecay]["Top"].mass, axis=None)[neg_evts], weight=evt_wt[neg_evts])
                        output["pt_top"].fill(dataset=f"{signame}_neg", sys=rewt_sys, pt=ak.flatten(events[ttdecay]["Top"].pt, axis=None)[neg_evts], weight=evt_wt[neg_evts])
                        output["pt_tt"].fill(dataset=f"{signame}_neg", sys=rewt_sys, pt=ak.flatten(events[ttdecay]["TTbar"].pt, axis=None)[neg_evts], weight=evt_wt[neg_evts])
                        output["eta_top"].fill(dataset=f"{signame}_neg", sys=rewt_sys, eta=ak.flatten(events[ttdecay]["Top"].eta, axis=None)[neg_evts], weight=evt_wt[neg_evts])
                        output["eta_tt"].fill(dataset=f"{signame}_neg", sys=rewt_sys, eta=ak.flatten(events[ttdecay]["TTbar"].eta, axis=None)[neg_evts], weight=evt_wt[neg_evts])
                        output["top_ctstar"].fill(dataset=f"{signame}_neg", sys=rewt_sys, ctstar=top_ctstar[neg_evts], weight=evt_wt[neg_evts])
                        output["top_ctstar_abs"].fill(dataset=f"{signame}_neg", sys=rewt_sys, ctstar_abs=np.abs(top_ctstar)[neg_evts], weight=evt_wt[neg_evts])
                        output["mtt_vs_top_ctstar_abs"].fill(dataset=f"{signame}_neg", sys=rewt_sys, mtt=ak.flatten(events[ttdecay]["TTbar"].mass, axis=None)[neg_evts], ctstar_abs=np.abs(top_ctstar)[neg_evts], weight=evt_wt[neg_evts])
                    else:
                        output["mtt"].fill(dataset=signame, sys=rewt_sys, mtt=ak.flatten(events[ttdecay]["TTbar"].mass, axis=None), weight=evt_wt)
                        output["mtop"].fill(dataset=signame, sys=rewt_sys, mtop=ak.flatten(events[ttdecay]["Top"].mass, axis=None), weight=evt_wt)
                        output["pt_top"].fill(dataset=signame, sys=rewt_sys, pt=ak.flatten(events[ttdecay]["Top"].pt, axis=None), weight=evt_wt)
                        output["pt_tt"].fill(dataset=signame, sys=rewt_sys, pt=ak.flatten(events[ttdecay]["TTbar"].pt, axis=None), weight=evt_wt)
                        output["eta_top"].fill(dataset=signame, sys=rewt_sys, eta=ak.flatten(events[ttdecay]["Top"].eta, axis=None), weight=evt_wt)
                        output["eta_tt"].fill(dataset=signame, sys=rewt_sys, eta=ak.flatten(events[ttdecay]["TTbar"].eta, axis=None), weight=evt_wt)
                        output["top_ctstar"].fill(dataset=signame, sys=rewt_sys, ctstar=top_ctstar, weight=evt_wt)
                        output["top_ctstar_abs"].fill(dataset=signame, sys=rewt_sys, ctstar_abs=np.abs(top_ctstar), weight=evt_wt)
                        output["mtt_vs_top_ctstar_abs"].fill(dataset=signame, sys=rewt_sys, mtt=ak.flatten(events[ttdecay]["TTbar"].mass, axis=None), ctstar_abs=np.abs(top_ctstar), weight=evt_wt)

        return output


    def postprocess(self, accumulator):
        return accumulator


proc_executor = processor.iterative_executor if to_debug else processor.futures_executor
proc_exec_args = {"schema": processor.NanoAODSchema} if to_debug else {"schema": processor.NanoAODSchema, "workers": 8}
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

toc = time.time()
print("Total time: %.1f" % (toc - tic))
