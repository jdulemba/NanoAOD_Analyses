#!/usr/bin/env python

import time
tic = time.time()

from coffea import hist, processor
from pdb import set_trace
import os
from coffea.util import save, load
import Utilities.prettyjson as prettyjson
import numpy as np
import python.ObjectSelection as objsel
import python.MCWeights as MCWeights
import fnmatch
from coffea.analysis_tools import PackedSelection
import awkward as ak
from coffea.nanoevents.methods import vector
ak.behavior.update(vector.behavior)
import python.IDJet as IDJet

proj_dir = os.environ["PROJECT_DIR"]
base_jobid = os.environ["base_jobid"]
jobid = os.environ["jobid"]
analyzer = "htt_flav_effs"

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

allowed_samples = ["ttJets*", "WJets*", "singlet*", "ZJets*"]
for fname in fileset.keys():
    if not any([fnmatch.fnmatch(fname, sample) for sample in allowed_samples]):
        raise IOError(f"Sample {fname} not valid for finding btag efficiencies. Only {allowed_samples} are allowed")

samplename = list(fileset.keys())[0]
isTTbar_ = samplename.startswith("ttJets")

## load corrections for event weights
cfg_pars = prettyjson.loads(open(os.path.join(proj_dir, "cfg_files", f"cfg_pars_{jobid}.json")).read())
pu_correction = load(os.path.join(proj_dir, "Corrections", base_jobid, cfg_pars["corrections"]["pu"]))[args.year]
lepSF_correction = load(os.path.join(proj_dir, "Corrections", base_jobid, cfg_pars["corrections"]["lepton"]))[args.year]
jet_corrections = load(os.path.join(proj_dir, "Corrections", base_jobid, cfg_pars["corrections"]["jetmet"][cfg_pars["corrections"]["jetmet"]["to_use"]]))[args.year]
nnlo_reweighting = load(os.path.join(proj_dir, "Corrections", base_jobid, cfg_pars["corrections"]["nnlo"]["filename"]))
ewk_reweighting = load(os.path.join(proj_dir, "Corrections", base_jobid, cfg_pars["corrections"]["ewk"]["file"]))
corrections = {
    "Pileup" : pu_correction,
    "Prefire" : True,
    "LeptonSF" : lepSF_correction,
    "JetCor" : jet_corrections,
    "BTagSF" : False,
    "NNLO_Rewt" : {"Var" : cfg_pars["corrections"]["nnlo"]["var"], "Correction" : nnlo_reweighting[cfg_pars["corrections"]["nnlo"]["var"]]},
    "EWK_Rewt" : {"Correction" : ewk_reweighting, "wt" : cfg_pars["corrections"]["ewk"]["wt"]},
}

jet_pars = cfg_pars["Jets"]
wps_to_use = list(set([jet_pars["permutations"]["tightb"], jet_pars["permutations"]["looseb"]]))
if not( len(wps_to_use) == 1):
    raise IOError("Only 1 unique btag working point supported now")


# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class Htt_Flav_Effs(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.btagger_axis = hist.Cat("btagger", "B-Tag WP")
        self.jetmult_axis = hist.Cat("jmult", "nJets")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")
        self.hflav_axis = hist.Cat("hFlav", "Hadron Flavour")
        self.pt_axis = hist.Bin("pt", "p_{T} [GeV]", np.around(np.linspace(0., 1000., 201), decimals=0))
        self.eta_axis = hist.Bin("eta", r"$\eta$", np.around(np.linspace(-3., 3., 121), decimals=2))

            ## make dictionary of hists
        histo_dict = {}
                ## make jet hists
        jet_hists = self.make_jet_hists()
        histo_dict.update(jet_hists)

        histo_dict["cutflow"] = processor.defaultdict_accumulator(int)

        self._accumulator = processor.dict_accumulator(histo_dict)
        self.sample_name = ""
        self.corrections = corrections

        self.regions = {
            "Muon" : {
                "3Jets"  : {
                    "lep_and_filter_pass", "passing_jets", "jets_3", "loose_or_tight_MU"
                },
            },
            "Electron" : {
                "3Jets"  : {
                    "lep_and_filter_pass", "passing_jets", "jets_3", "loose_or_tight_EL"
                },
            },
        }
        if isTTbar_:
            ## add 4+ jets categories for ttbar events
            self.regions["Muon"].update({
                "4PJets" : {"lep_and_filter_pass", "passing_jets", "jets_4+", "loose_or_tight_MU"}
            })
            self.regions["Electron"].update({
                "4PJets" : {"lep_and_filter_pass", "passing_jets", "jets_4+", "loose_or_tight_EL"}
            })
    
    @property
    def accumulator(self):
        return self._accumulator


    def make_jet_hists(self):
        histo_dict = {}
        histo_dict["Jets_pt_eta_all"] = hist.Hist("Events", self.btagger_axis, self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.hflav_axis, self.pt_axis, self.eta_axis)
        histo_dict["Jets_pt_eta_pass"]= hist.Hist("Events", self.btagger_axis, self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.hflav_axis, self.pt_axis, self.eta_axis)

        return histo_dict



    def process(self, events):
        output = self.accumulator.identity()

        self.sample_name = events.metadata["dataset"]

            ## make event weights
                # data or MC distinction made internally
        mu_evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections)
        el_evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections)

            ## initialize selections and regions
        selection = PackedSelection()

            # get all passing leptons
        lep_and_filter_pass = objsel.select_leptons(events, year=args.year, cutflow=output["cutflow"])
        selection.add("lep_and_filter_pass", lep_and_filter_pass)
                ## muons
        selection.add("tight_MU", ak.sum(events["Muon"]["TIGHTMU"], axis=1) == 1) # one muon passing TIGHT criteria
        selection.add("loose_or_tight_MU", ak.sum(events["Muon"]["LOOSEMU"] | events["Muon"]["TIGHTMU"], axis=1) == 1) # one muon passing LOOSE or TIGHT criteria
                ## electrons
        selection.add("tight_EL", ak.sum(events["Electron"]["TIGHTEL"], axis=1) == 1) # one muon passing TIGHT criteria
        selection.add("loose_or_tight_EL", ak.sum(events["Electron"]["LOOSEEL"] | events["Electron"]["TIGHTEL"], axis=1) == 1) # one muon passing LOOSE or TIGHT criteria

            ## build corrected jets and MET
        events["Jet"], events["MET"] = IDJet.process_jets(events, args.year, self.corrections["JetCor"])

            # jet selection
        passing_jets = objsel.jets_selection(events, year=args.year, cutflow=output["cutflow"])
        selection.add("passing_jets", passing_jets)
        output["cutflow"]["nEvts passing jet and lepton obj selection"] += ak.sum(passing_jets & lep_and_filter_pass)
        selection.add("jets_3", ak.num(events["SelectedJets"]) == 3)

        ## apply lepton SFs to MC (only applicable to tight leptons)
        if "LeptonSF" in self.corrections.keys():
            tight_mu_cut = selection.require(tight_MU=True) # find events passing muon object selection with one tight muon
            tight_muons = events["Muon"][tight_mu_cut][(events["Muon"][tight_mu_cut]["TIGHTMU"] == True)]
            muSFs_dict =  MCWeights.get_lepton_sf(lepton="Muons", corrections=self.corrections["LeptonSF"],
                pt=ak.flatten(tight_muons["pt"]), eta=ak.flatten(tight_muons["eta"]))
            for source in muSFs_dict.keys():
                tmp_wts = np.ones(len(events))
                tmp_wts[tight_mu_cut] = muSFs_dict[source]["Central"]
                mu_evt_weights.add(source, np.copy(tmp_wts))

            tight_el_cut = selection.require(tight_EL=True) # find events passing electron object selection with one tight electron
            tight_electrons = events["Electron"][tight_el_cut][(events["Electron"][tight_el_cut]["TIGHTEL"] == True)]
            elSFs_dict = MCWeights.get_lepton_sf(lepton="Electrons", corrections=self.corrections["LeptonSF"],
                pt=ak.flatten(tight_electrons["pt"]), eta=ak.flatten(tight_electrons["etaSC"]))
            #set_trace()
            for source in elSFs_dict.keys():
                tmp_wts = np.ones(len(events))
                tmp_wts[tight_el_cut] = elSFs_dict[source]["Central"]
                el_evt_weights.add(source, np.copy(tmp_wts))

        if isTTbar_:
            ## add 4+ jets categories for ttbar events
            selection.add("jets_4+", ak.num(events["SelectedJets"]) > 3)

            if "NNLO_Rewt" in self.corrections.keys():
                    # find gen level particles for ttbar system
                nnlo_wts = MCWeights.get_nnlo_weights(self.corrections["NNLO_Rewt"], events)
                mu_evt_weights.add("NNLOqcd",
                    np.copy(nnlo_wts),
                )
                el_evt_weights.add("NNLOqcd",
                    np.copy(nnlo_wts),
                )

            if "EWK_Rewt" in self.corrections.keys():
                #set_trace()
                    ## NLO EW weights
                if self.corrections["EWK_Rewt"]["wt"] == "Otto":
                    ewk_wts_dict = MCWeights.get_Otto_ewk_weights(self.corrections["EWK_Rewt"]["Correction"], events)
                        # add Yukawa coupling variation
                    mu_evt_weights.add("Yukawa",  # really just varying value of Yt
                        np.copy(ewk_wts_dict["Rebinned_KFactor_1.0"]),
                        np.copy(ewk_wts_dict["Rebinned_KFactor_1.1"]),
                        np.copy(ewk_wts_dict["Rebinned_KFactor_0.9"]),
                    )
                    el_evt_weights.add("Yukawa",  # really just varying value of Yt
                        np.copy(ewk_wts_dict["Rebinned_KFactor_1.0"]),
                        np.copy(ewk_wts_dict["Rebinned_KFactor_1.1"]),
                        np.copy(ewk_wts_dict["Rebinned_KFactor_0.9"]),
                    )


        btag_wps = [wp for wp in events["SelectedJets"].fields if wps_to_use[0] in wp]
        #set_trace()
        ## fill hists for each region
        for btag_wp in btag_wps:
            for lepton in self.regions.keys():
                evt_weights = mu_evt_weights if lepton == "Muon" else el_evt_weights
                for jmult in self.regions[lepton].keys():
                    cut = selection.all(*self.regions[lepton][jmult])

                    evt_weights_to_use = evt_weights.weight()
                    jets = events["SelectedJets"][cut]
                    #if to_debug: set_trace()
                        ## hadronFlavour definitions found here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools
                    bjets = jets[(jets["hadronFlavour"] == 5)]
                    cjets = jets[(jets["hadronFlavour"] == 4)]
                    ljets = jets[(jets["hadronFlavour"] == 0)]

                    output = self.fill_jet_hists_all( acc=output, btag_wp=btag_wp, jetmult=jmult, leptype=lepton, hadFlav="bjet", obj=bjets, evt_weights=evt_weights_to_use[cut])
                    output = self.fill_jet_hists_pass(acc=output, btag_wp=btag_wp, jetmult=jmult, leptype=lepton, hadFlav="bjet", obj=bjets[(bjets[btag_wp])], evt_weights=evt_weights_to_use[cut])
                    output = self.fill_jet_hists_all( acc=output, btag_wp=btag_wp, jetmult=jmult, leptype=lepton, hadFlav="cjet", obj=cjets, evt_weights=evt_weights_to_use[cut])
                    output = self.fill_jet_hists_pass(acc=output, btag_wp=btag_wp, jetmult=jmult, leptype=lepton, hadFlav="cjet", obj=cjets[(cjets[btag_wp])], evt_weights=evt_weights_to_use[cut])
                    output = self.fill_jet_hists_all( acc=output, btag_wp=btag_wp, jetmult=jmult, leptype=lepton, hadFlav="ljet", obj=ljets, evt_weights=evt_weights_to_use[cut])
                    output = self.fill_jet_hists_pass(acc=output, btag_wp=btag_wp, jetmult=jmult, leptype=lepton, hadFlav="ljet", obj=ljets[(ljets[btag_wp])], evt_weights=evt_weights_to_use[cut])

        return output


    def fill_jet_hists_all(self, acc, btag_wp, jetmult, leptype, hadFlav, obj, evt_weights):
        #set_trace()
        acc["Jets_pt_eta_all"].fill(btagger=btag_wp, dataset=self.sample_name, jmult=jetmult, leptype=leptype, hFlav=hadFlav, pt=ak.flatten(obj["pt"]), eta=ak.flatten(obj["eta"]), weight=ak.flatten((ak.ones_like(obj["pt"])*evt_weights)))
        return acc        

    def fill_jet_hists_pass(self, acc, btag_wp, jetmult, leptype, hadFlav, obj, evt_weights):
        #set_trace()
        acc["Jets_pt_eta_pass"].fill(btagger=btag_wp, dataset=self.sample_name, jmult=jetmult, leptype=leptype, hFlav=hadFlav, pt=ak.flatten(obj["pt"]), eta=ak.flatten(obj["eta"]), weight=ak.flatten((ak.ones_like(obj["pt"])*evt_weights)))
        return acc        


    def postprocess(self, accumulator):
        return accumulator


proc_executor = processor.iterative_executor if to_debug else processor.futures_executor
proc_exec_args = {"schema": processor.NanoAODSchema} if to_debug else {"schema": processor.NanoAODSchema, "workers": 8}
output = processor.run_uproot_job(
    fileset,
    treename="Events",
    processor_instance=Htt_Flav_Effs(),
    executor=proc_executor,
    executor_args=proc_exec_args,
    chunksize=100000,
)

save(output, args.outfname)
print(f"{args.outfname} has been written")

toc = time.time()
print("Total time: %.1f" % (toc - tic))
