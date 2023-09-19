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
parser.add_argument("year", choices=["2016preVFP", "2016postVFP", "2017", "2018"] if base_jobid == "ULnanoAODv9" else ["2016APV", "2016", "2017", "2018"], help="Specify which year to run over")
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

isData_ = samplename.startswith("data")
if isData_:
    lumiMask_path = os.path.join(proj_dir, "inputs", "data", base_jobid, "LumiMasks", f"{args.year}_GoldenJson_{base_jobid}.txt")


        # copy fileset root files to local condor node if running on condor
if "isCondor" in opts_dict.keys():
    if ast.literal_eval(opts_dict["isCondor"]):
        import subprocess
        sites_to_try = ["root://xrootd-cms.infn.it/", "root://cmsxrootd.fnal.gov/"]
        n_retries = len(sites_to_try) + 1
        for idx, rfile in enumerate(fileset[samplename]):
            cp_success = False
            for cp_attempt in range(n_retries):
                if cp_success: continue
                cp_rfile = rfile if cp_attempt == 0 else "/".join([sites_to_try[cp_attempt-1], rfile.split("//")[-1]]) # replace whatever redirector is used to regional Bari one
                print(f"Attempt {cp_attempt+1} to copy {cp_rfile} to /tmp")
                try:
                    subprocess.run(["xrdcp", "-f", f"{cp_rfile}", "/tmp"], timeout=300)
                    cp_success = True
                except:
                    cp_success = False
                    continue
            if not cp_success:
                raise ValueError(f"{cp_rfile} not found")
            fileset[samplename][idx] = f"/tmp/{rfile.split('/')[-1]}"

if samplename.startswith("Toponium"):
        ## cuts specific to toponium mass window
    nominal_eta_mass = 343
    def eta_mass_min(eta_mass):
        return eta_mass - 6
    def eta_mass_max(eta_mass):
        return eta_mass + 6
    
    #set_trace()
    toponium_variations = {
        "nosys"             : [eta_mass_min(nominal_eta_mass),     eta_mass_max(nominal_eta_mass)    ],# f"{decay.split('_')[0]}_meta343_mtop172p5"],
        "BindingEnergyUp"   : [eta_mass_min(nominal_eta_mass + 1), eta_mass_max(nominal_eta_mass + 1)],# f"{decay.split('_')[0]}_meta344_mtop172p5"],
        "BindingEnergyDown" : [eta_mass_min(nominal_eta_mass - 1), eta_mass_max(nominal_eta_mass - 1)],# f"{decay.split('_')[0]}_meta342_mtop172p5"],
        "TopMassUp"         : [eta_mass_min(nominal_eta_mass + 2), eta_mass_max(nominal_eta_mass + 2)],# f"{decay.split('_')[0]}_meta345_mtop173p5"],
        "TopMassDown"       : [eta_mass_min(nominal_eta_mass - 2), eta_mass_max(nominal_eta_mass - 2)],# f"{decay.split('_')[0]}_meta341_mtop171p5"],
    }

# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class Meta_Analyzer(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.pu_nTrueInt_axis = hist.Bin("pu_nTrueInt", "nTrueInt", 100, 0, 100)
        self.pu_nPU_axis = hist.Bin("pu_nPU", "nPU", 100, 0, 100)

            ## make dictionary of hists
        histo_dict = {}
        histo_dict["PU_nTrueInt"] = hist.Hist("PU_nTrueInt", self.dataset_axis, self.pu_nTrueInt_axis)
        histo_dict["PU_nPU"] = hist.Hist("PU_nPU", self.dataset_axis, self.pu_nPU_axis)

            ## construct dictionary of dictionaries to hold meta info for each sample
        for sample in fileset.keys():
            if "Int" in sample:
                histo_dict[f"{sample}_pos"] = processor.defaultdict_accumulator(int)
                histo_dict[f"{sample}_pos_runs_to_lumis"] = processor.value_accumulator(list)
                histo_dict[f"{sample}_neg"] = processor.defaultdict_accumulator(int)
                histo_dict[f"{sample}_neg_runs_to_lumis"] = processor.value_accumulator(list)
            else:
                histo_dict[sample] = processor.defaultdict_accumulator(int)
                histo_dict[f"{sample}_runs_to_lumis"] = processor.value_accumulator(list)

        self._accumulator = processor.dict_accumulator(histo_dict)
        self.sample_name = ""
    
    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        output = self.accumulator.identity()

        #set_trace()
        event_nums = events.event
        self.sample_name = events.metadata["dataset"]

        if isData_:
            runs = events.run
            lumis = events.luminosityBlock
            Golden_Json_LumiMask = lumi_tools.LumiMask(lumiMask_path)
            LumiMask = Golden_Json_LumiMask.__call__(runs, lumis) ## returns array of valid events

            if ak.size(event_nums[LumiMask]) > 0:
                valid_runs_lumis = np.unique(np.stack((ak.to_numpy(runs[LumiMask]), ak.to_numpy(lumis[LumiMask])), axis=1), axis=0) ## make 2D array of uniqe valid [[run, lumi], [run, lumi]...] pairs
                    # make dictionary of valid runs: sorted list of unique lumisections for each valid run
                lumi_map = {str(valid_run):sorted(list(set(valid_runs_lumis[:, 1][valid_runs_lumis[:, 0] == valid_run]))) for valid_run in list(set(valid_runs_lumis[:, 0]))}

                output[f"{self.sample_name}_runs_to_lumis"].add(list(lumi_map.items()))

        else:
            genWeights = events["genWeight"]

                # split interference into pos/neg weighted events
            if "Int" in self.sample_name:
                pos_evts = ak.where(genWeights > 0)
                neg_evts = ak.where(genWeights < 0)

                output[f"{self.sample_name}_pos"]["nEvents"] += ak.size(event_nums[pos_evts])
                output[f"{self.sample_name}_neg"]["nEvents"] += ak.size(event_nums[neg_evts])
                output[f"{self.sample_name}_pos"]["sumGenWeights"] += sum(genWeights[pos_evts])
                output[f"{self.sample_name}_neg"]["sumGenWeights"] += sum(genWeights[neg_evts])
                output[f"{self.sample_name}_pos"]["sumSquaredGenWeights"] += sum(np.square(genWeights[pos_evts]))
                output[f"{self.sample_name}_neg"]["sumSquaredGenWeights"] += sum(np.square(genWeights[neg_evts]))

                output["PU_nTrueInt"].fill(dataset=f"{self.sample_name}_pos", pu_nTrueInt=events["Pileup"]["nTrueInt"][pos_evts], weight=genWeights[pos_evts])
                output["PU_nPU"].fill(dataset=f"{self.sample_name}_pos", pu_nPU=events["Pileup"]["nPU"][pos_evts], weight=genWeights[pos_evts])
                output["PU_nTrueInt"].fill(dataset=f"{self.sample_name}_neg", pu_nTrueInt=events["Pileup"]["nTrueInt"][neg_evts], weight=genWeights[neg_evts])
                output["PU_nPU"].fill(dataset=f"{self.sample_name}_neg", pu_nPU=events["Pileup"]["nPU"][neg_evts], weight=genWeights[neg_evts])

                if "LHEPdfWeight" in events.fields:
                    LHEpdfWeights = events["LHEPdfWeight"]
                        ## get sum of each pdf weight over all events
                    sumLHEpdfWeights_pos = ak.sum(LHEpdfWeights[pos_evts]*genWeights[pos_evts], axis=0)
                    output[f"{self.sample_name}_pos"]["sumLHEpdfWeights"] += ak.to_numpy(sumLHEpdfWeights_pos)
                    sumLHEpdfWeights_neg = ak.sum(LHEpdfWeights[neg_evts]*genWeights[neg_evts], axis=0)
                    output[f"{self.sample_name}_neg"]["sumLHEpdfWeights"] += ak.to_numpy(sumLHEpdfWeights_neg)

                if "PSWeight" in events.fields:
                    psweights = events["PSWeight"]
                        ## get sum of each weight over all events
                    sumPSweights_pos = ak.sum(psweights[pos_evts]*genWeights[pos_evts], axis=0)
                    output[f"{self.sample_name}_pos"]["sumPSWeights"] += ak.to_numpy(sumPSweights_pos)
                    sumPSweights_neg = ak.sum(psweights[neg_evts]*genWeights[neg_evts], axis=0)
                    output[f"{self.sample_name}_neg"]["sumPSWeights"] += ak.to_numpy(sumPSweights_neg)

                if "LHEScaleWeight" in events.fields:
                    lheweights = events["LHEScaleWeight"]
                        ## get sum of each weight over all events
                    sumLHEscaleWeights_pos = ak.sum(lheweights[pos_evts]*genWeights[pos_evts], axis=0)
                    output[f"{self.sample_name}_pos"]["sumLHEscaleWeights"] += ak.to_numpy(sumLHEscaleWeights_pos)
                    sumLHEscaleWeights_neg = ak.sum(lheweights[neg_evts]*genWeights[neg_evts], axis=0)
                    output[f"{self.sample_name}_neg"]["sumLHEscaleWeights"] += ak.to_numpy(sumLHEscaleWeights_neg)

            else:
                if "Toponium" in self.sample_name: ## compute meta info for toponium events where the mass window cuts are different
                        ## make LHE mass cute for eta_t at +/- 6 GeV
                    lhepart = events["LHEPart"][events["LHEPart"]["status"] == 1]
                    mWWbb = (lhepart[:,0] + lhepart[:,1] + lhepart[:,2] + lhepart[:,3] + lhepart[:,4] + lhepart[:,5]).mass
                    for sys in toponium_variations.keys():
                        mass_cut = (mWWbb >= toponium_variations[sys][0]) & (mWWbb <= toponium_variations[sys][1])

                        output[self.sample_name][f"nEvents_{sys}"] += ak.size(event_nums[mass_cut])
                        output[self.sample_name][f"sumGenWeights_{sys}"] += sum(genWeights[mass_cut])
                        output[self.sample_name][f"sumSquaredGenWeights_{sys}"] += sum(np.square(genWeights[mass_cut]))

                        output["PU_nTrueInt"].fill(dataset=f"{self.sample_name}_{sys}", pu_nTrueInt=events["Pileup"]["nTrueInt"][mass_cut], weight=genWeights[mass_cut])
                        output["PU_nPU"].fill(dataset=f"{self.sample_name}_{sys}", pu_nPU=events["Pileup"]["nPU"][mass_cut], weight=genWeights[mass_cut])

                        if "LHEPdfWeight" in events.fields:
                            LHEpdfWeights = events["LHEPdfWeight"]
                                ## get sum of each pdf weight over all events
                            sumLHEpdfWeights = ak.sum((LHEpdfWeights*genWeights)[mass_cut], axis=0)
                            output[self.sample_name][f"sumLHEpdfWeights_{sys}"] += ak.to_numpy(sumLHEpdfWeights)

                        if "PSWeight" in events.fields:
                            psweights = events["PSWeight"]
                                ## get sum of each weight over all events
                            sumPSweights = ak.sum((psweights*genWeights)[mass_cut], axis=0)
                            output[self.sample_name][f"sumPSWeights_{sys}"] += ak.to_numpy(sumPSweights)

                        if "LHEScaleWeight" in events.fields:
                            lheweights = events["LHEScaleWeight"]
                                ## get sum of each weight over all events
                            sumLHEscaleWeights = ak.sum((lheweights*genWeights)[mass_cut], axis=0)
                            output[self.sample_name][f"sumLHEscaleWeights_{sys}"] += ak.to_numpy(sumLHEscaleWeights)

                output[self.sample_name]["nEvents"] += ak.size(event_nums)
                output[self.sample_name]["sumGenWeights"] += sum(genWeights)
                output[self.sample_name]["sumSquaredGenWeights"] += sum(np.square(genWeights))

                output["PU_nTrueInt"].fill(dataset=self.sample_name, pu_nTrueInt=events["Pileup"]["nTrueInt"], weight=genWeights)
                output["PU_nPU"].fill(dataset=self.sample_name, pu_nPU=events["Pileup"]["nPU"], weight=genWeights)

                if "LHEPdfWeight" in events.fields:
                    LHEpdfWeights = events["LHEPdfWeight"]
                        ## get sum of each pdf weight over all events
                    sumLHEpdfWeights = ak.sum(LHEpdfWeights*genWeights, axis=0)
                    output[self.sample_name]["sumLHEpdfWeights"] += ak.to_numpy(sumLHEpdfWeights)

                if "PSWeight" in events.fields:
                    psweights = events["PSWeight"]
                        ## get sum of each weight over all events
                    sumPSweights = ak.sum(psweights*genWeights, axis=0)
                    output[self.sample_name]["sumPSWeights"] += ak.to_numpy(sumPSweights)

                if "LHEScaleWeight" in events.fields:
                    lheweights = events["LHEScaleWeight"]
                        ## get sum of each weight over all events
                    sumLHEscaleWeights = ak.sum(lheweights*genWeights, axis=0)
                    output[self.sample_name]["sumLHEscaleWeights"] += ak.to_numpy(sumLHEscaleWeights)


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
    processor_instance=Meta_Analyzer(),
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
