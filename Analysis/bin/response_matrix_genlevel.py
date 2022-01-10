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
import Utilities.final_analysis_binning as final_binning
import python.GenParticleSelector as genpsel

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

isTTSL_ = samplename == "ttJetsSL"
isSignal_ = (samplename.startswith("AtoTT") or samplename.startswith("HtoTT"))
isInt_ = isSignal_ and ("Int" in samplename)

if not (isTTSL_ or isSignal_):
    raise ValueError("This analyzer should only be run with tt l+jets or signal events")

# convert input string of options dictionary to actual dictionary
odict = (args.opts).replace("\'", "\"")
opts_dict = prettyjson.loads(odict)

    ## set config options passed through argparse
import ast
to_debug = ast.literal_eval(opts_dict.get("debug", "False"))

corrections = {
    "Prefire" : False,
}

# get systematics to run
reweight_systematics_to_run = ["nosys", "ISRUp", "ISRDown", "FSRUp", "FSRDown", "FACTORUp", "FACTORDown", "RENORMUp", "RENORMDown"] if isTTSL_ else ["nosys", "AH_FACTORUp", "AH_FACTORDown", "AH_RENORMUp", "AH_RENORMDown"]
print("\n\nRunning with systematics:", *sorted(reweight_systematics_to_run))


# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class response_matrix_genlevel(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.sys_axis = hist.Cat("sys", "Systematic")

        self.gen_st_mtt_ctstar_axis = hist.Bin("gen_st_mtt_ctstar_idx", "Gen st mtt ctstar",
            np.around(np.arange((final_binning.st_binning.size-1) * (final_binning.ctstar_abs_binning.size-1) * (final_binning.mtt_binning.size-1) + 1), decimals=0))

            ## make dictionary of hists
        histo_dict = {}
                ## make selection plots
        selection_hists = self.make_selection_hists()
        histo_dict.update(selection_hists)        

        self.sample_name = ""
        self.corrections = corrections
        self.reweight_systematics_to_run = reweight_systematics_to_run

        self._accumulator = processor.dict_accumulator(histo_dict)


    @property
    def accumulator(self):
        return self._accumulator



    def make_selection_hists(self):
        histo_dict = {}
        histo_dict["Gen_mtt_x_tlep_ctstar_abs_x_st_inds"] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.gen_st_mtt_ctstar_axis)

        return histo_dict


    def process(self, events):
        output = self.accumulator.identity()

        self.sample_name = events.metadata["dataset"]

        ## make event weights
        evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections, isTTbar=isTTSL_, isSignal=isSignal_)

        genpsel.select(events, mode="NORMAL")

        #set_trace()
        #    # get kinematic cuts on charged lepton and both b quarks
        #kin_cuts =  ak.flatten((events["SL"]["Lepton"]["pt"] >= 25.) & (np.abs(events["SL"]["Lepton"]["eta"]) <= 2.5) & \
        #            (events["SL"]["BHad"]["pt"] >= 25.) & (np.abs(events["SL"]["BHad"]["eta"]) <= 2.5) & \
        #            (events["SL"]["BLep"]["pt"] >= 25.) & (np.abs(events["SL"]["BLep"]["eta"]) <= 2.5), axis=1)

            # calculate ctstar
        gen_thad = events["SL"]["THad"]
        gen_tlep = events["SL"]["TLep"]
        thad_ctstar, tlep_ctstar = make_vars.ctstar(gen_thad, gen_tlep)
        thad_ctstar, tlep_ctstar = ak.flatten(thad_ctstar, axis=None), ak.flatten(tlep_ctstar, axis=None)
        tlep_ctstar_abs = np.abs(tlep_ctstar)
            # other variables
        mtt = ak.flatten(events["SL"]["TTbar"].mass, axis=None)
        ST = ak.flatten(gen_thad.pt + gen_tlep.pt, axis=None)

        # inds and masks
            # mtt
        mtt_ind_vals = np.array([np.argmax(mtt[idx] < final_binning.mtt_binning) for idx in range(len(mtt))])
        mtt_inds_mask = (mtt > final_binning.mtt_binning[0]) & (mtt <= final_binning.mtt_binning[-1])
            # ctstar
        ctstar_ind_vals = np.array([np.argmax(tlep_ctstar_abs[idx] < final_binning.ctstar_abs_binning) for idx in range(len(tlep_ctstar_abs))])
        ctstar_inds_mask = (tlep_ctstar_abs > final_binning.ctstar_abs_binning[0]) & (tlep_ctstar_abs <= final_binning.ctstar_abs_binning[-1])
            # ST
        st_ind_vals = np.array([np.argmax(ST[idx] < final_binning.st_binning) for idx in range(len(ST))])
        st_inds_mask = (ST > final_binning.st_binning[0]) & (ST <= final_binning.st_binning[-1])
            # ST x mtt x ctstar
        st_mtt_ctstar_ind_vals = (st_ind_vals - 1) * (final_binning.mtt_binning.size - 1) * (final_binning.ctstar_abs_binning.size - 1) + (ctstar_ind_vals - 1) * (final_binning.mtt_binning.size - 1) + mtt_ind_vals
                # with kinematic cut
        st_mtt_ctstar_inds_mask = mtt_inds_mask & ctstar_inds_mask & st_inds_mask #& kin_cuts

        #set_trace()
        for rewt_sys in self.reweight_systematics_to_run:
            evt_wt = evt_weights.weight() if rewt_sys == "nosys" else evt_weights.weight(rewt_sys)
            if isInt_:
                pos_evts = evt_wt > 0
                neg_evts = evt_wt < 0
                output["Gen_mtt_x_tlep_ctstar_abs_x_st_inds"].fill(dataset=f"{self.sample_name}_pos", sys=rewt_sys,
                    gen_st_mtt_ctstar_idx=st_mtt_ctstar_ind_vals[st_mtt_ctstar_inds_mask & pos_evts]-0.5, weight=evt_wt[st_mtt_ctstar_inds_mask & pos_evts])
                output["Gen_mtt_x_tlep_ctstar_abs_x_st_inds"].fill(dataset=f"{self.sample_name}_neg", sys=rewt_sys,
                    gen_st_mtt_ctstar_idx=st_mtt_ctstar_ind_vals[st_mtt_ctstar_inds_mask & neg_evts]-0.5, weight=evt_wt[st_mtt_ctstar_inds_mask & neg_evts])
            else:
                output["Gen_mtt_x_tlep_ctstar_abs_x_st_inds"].fill(dataset=self.sample_name, sys=rewt_sys,
                    gen_st_mtt_ctstar_idx=st_mtt_ctstar_ind_vals[st_mtt_ctstar_inds_mask]-0.5, weight=evt_wt[st_mtt_ctstar_inds_mask])

        return output


    def postprocess(self, accumulator):
        return accumulator


proc_executor = processor.iterative_executor if to_debug else processor.futures_executor
proc_exec_args = {"schema": processor.NanoAODSchema} if to_debug else {"schema": processor.NanoAODSchema, "workers": 8}
output = processor.run_uproot_job(
    fileset,
    treename="Events",
    processor_instance=response_matrix_genlevel(),
    executor=proc_executor,
    executor_args=proc_exec_args,
    chunksize=100000,
)

save(output, args.outfname)
print(f"{args.outfname} has been written")

toc = time.time()
print("Total time: %.1f" % (toc - tic))
