#!/usr/bin/env python

import time
tic = time.time()

from coffea import hist, processor
from coffea.nanoevents import NanoAODSchema
from coffea.analysis_tools import PackedSelection
import awkward as ak
from coffea.nanoevents.methods import vector
ak.behavior.update(vector.behavior)

from coffea.util import save, load
from pdb import set_trace
import os, sys
import python.ObjectSelection as objsel
import Utilities.plot_tools as plt_tools
import python.MCWeights as MCWeights
import Utilities.make_variables as make_vars
import numpy as np
import Utilities.prettyjson as prettyjson
import coffea.lumi_tools.lumi_tools as lumi_tools
import python.GenParticleSelector as genpsel
import python.TTGenMatcher as ttmatcher
import python.TTPermutator as ttpermutator
from python.Permutations import compare_matched_best_perms
import python.IDJet as IDJet

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
base_jobid = os.environ["base_jobid"]
analyzer = "ttbar_evt_cats"

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

    ## specify ttJets samples
if samplename != "ttJetsSL":
    raise ValueError("Should only be run with ttJetsSL")

# convert input string of options dictionary to actual dictionary
odict = (args.opts).replace("\'", "\"")
opts_dict = prettyjson.loads(odict)

    ## set config options passed through argparse
import ast
to_debug = ast.literal_eval(opts_dict.get("debug", "False"))
apply_hem = ast.literal_eval(opts_dict.get("apply_hem", "True"))

## init tt probs for likelihoods
ttpermutator.year_to_run(year=args.year)

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
                    #output = check_output(["xrdcp", "-f", f"{cp_rfile}", "/tmp"], timeout=300, stderr=STDOUT)
                    cp_success = True
                except:
                    cp_success = False
                    continue
            if not cp_success:
                raise ValueError(f"{cp_rfile} not found")
            fileset[samplename][idx] = f"/tmp/{rfile.split('/')[-1]}"


cfg_pars = prettyjson.loads(open(os.path.join(proj_dir, "cfg_files", f"cfg_pars_{jobid}.json")).read())
## load corrections for event weights
jet_corrections = load(os.path.join(proj_dir, "Corrections", base_jobid, cfg_pars["corrections"]["jetmet"][cfg_pars["corrections"]["jetmet"]["to_use"]]))[args.year]
corrections = {
    "Prefire" : False,
    "BTagSF" : False,
    "JetCor" : jet_corrections,
}

    ## parameters for b-tagging
jet_pars = cfg_pars["Jets"]
btaggers = [jet_pars["btagger"]]

wps_to_use = list(set([jet_pars["permutations"]["tightb"],jet_pars["permutations"]["looseb"]]))
if not( len(wps_to_use) == 1):
    raise IOError("Only 1 unique btag working point supported now")
btag_wp = btaggers[0]+wps_to_use[0]

MTcut = jet_pars["MT"]


# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class ttbar_evt_cats(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")
        self.jetmult_axis = hist.Cat("jmult", "nJets")
        self.object_axis = hist.Cat("objtype", "Gen Object")
        self.nperms_axis = hist.Bin("nperms", "n_{perms}", np.around(np.linspace(0., 1., 2), decimals=0))
        self.ntotal_correct_axis = hist.Bin("nID", "nID", np.around(np.linspace(0., 2., 3), decimals=0))
        self.permcat_axis = hist.Bin("permcat", "permcat", np.around(np.linspace(0., 5., 6), decimals=0))

            ## make dictionary of hists
        histo_dict = {}
                ## make jet hists
        gen_hists = self.make_hists()
        histo_dict.update(gen_hists)

        histo_dict["cutflow"] = processor.defaultdict_accumulator(int)

        self._accumulator = processor.dict_accumulator(histo_dict)
        self.sample_name = ""
        self.corrections = corrections

        self.regions = {
            "Muon" : {
                "Tight" : {
                    "3Jets"  : {"lep_and_filter_pass", "passing_jets", "jets_3" , "tight_MU", f"{btaggers[0]}_pass", "semilep"},
                    "4PJets"  : {"lep_and_filter_pass", "passing_jets", "jets_4p" , "tight_MU", f"{btaggers[0]}_pass", "semilep"},
                },
            },
            "Electron" : {
                "Tight" : {
                    "3Jets"  : {"lep_and_filter_pass", "passing_jets", "jets_3" , "tight_EL", f"{btaggers[0]}_pass", "semilep"},
                    "4PJets"  : {"lep_and_filter_pass", "passing_jets", "jets_4p" , "tight_EL", f"{btaggers[0]}_pass", "semilep"},
                },
            },
        }
    
    @property
    def accumulator(self):
        return self._accumulator


    def make_hists(self):
        histo_dict = {}
        histo_dict["nBestPerms_Total"]      = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.nperms_axis)
        histo_dict["nBestPerms_Solutions"]  = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.nperms_axis)
        histo_dict["Object_Identifying"]  = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.object_axis, self.ntotal_correct_axis)
        histo_dict["ValidBestPerms"]      = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.permcat_axis)

        return histo_dict



    def process(self, events):
        np.random.seed(10) # sets seed so values from random distributions are reproducible (JER corrections)
        output = self.accumulator.identity()

        self.sample_name = events.metadata["dataset"]

            ## make event weights
                # data or MC distinction made internally
        mu_evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections)
        el_evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections)

            ## initialize selections and regions
        selection = PackedSelection()

            # get all passing leptons
        lep_and_filter_pass = objsel.select_leptons(events, year=args.year, hem_15_16=apply_hem)
        output["cutflow"]["lep_and_filter_pass"] += ak.sum(lep_and_filter_pass)
        selection.add("lep_and_filter_pass", lep_and_filter_pass)

            ## build corrected jets and MET
        events["Jet"], events["MET"] = IDJet.process_jets(events, args.year, self.corrections["JetCor"])

        ## add different selections
                ## muons
        tight_mu_sel = ak.sum(events["Muon"]["TIGHTMU"], axis=1) == 1
        selection.add("tight_MU", tight_mu_sel) # one muon passing TIGHT criteria
                ## electrons
        tight_el_sel = ak.sum(events["Electron"]["TIGHTEL"], axis=1) == 1
        selection.add("tight_EL", tight_el_sel) # one electron passing TIGHT criteria

        ## find gen level particles for ttbar system
        genpsel.select(events, mode="NORMAL")
        selection.add("semilep", ak.num(events["SL"]) > 0)

            # jet selection
        passing_jets = objsel.jets_selection(events, year=args.year, cutflow=output["cutflow"], hem_15_16=apply_hem)
        output["cutflow"]["nEvts passing jet and lepton obj selection"] += ak.sum(passing_jets & lep_and_filter_pass)
        selection.add("passing_jets", passing_jets)
        selection.add("jets_3",  ak.num(events["SelectedJets"]) == 3)
        selection.add("jets_4p",  ak.num(events["SelectedJets"]) > 3) # only for getting btag weights
        selection.add(f"{btaggers[0]}_pass", ak.sum(events["SelectedJets"][btag_wp], axis=1) >= 2)

            # sort jets by btag value
        events["SelectedJets"] = events["SelectedJets"][ak.argsort(events["SelectedJets"][IDJet.btag_tagger_to_disc_name[btaggers[0]]], ascending=False)]

        ## fill hists for each region
        for lepton in self.regions.keys():
            evt_weights = mu_evt_weights if lepton == "Muon" else el_evt_weights
            for leptype in self.regions[lepton].keys():
                for jmult in self.regions[lepton][leptype].keys():

                    if to_debug: print(lepton, leptype, "btagPass", jmult)
                    cut = selection.all(*self.regions[lepton][leptype][jmult])
                    output["cutflow"]["nEvts %s" % ", ".join([lepton, leptype, "btagPass", jmult])] += cut.sum()
                    if cut.sum() == 0: continue

                    ltype = "MU" if lepton == "Muon" else "EL"
                    if f"loose_or_tight_{ltype}" in self.regions[lepton][leptype][jmult]:
                        leptons = events[lepton][cut][((events[lepton][cut][f"TIGHT{ltype}"] == True) | (events[lepton][cut][f"LOOSE{ltype}"] == True))]
                    elif f"tight_{ltype}" in self.regions[lepton][leptype][jmult]:
                        leptons = events[lepton][cut][(events[lepton][cut][f"TIGHT{ltype}"] == True)]
                    elif f"loose_{ltype}" in self.regions[lepton][leptype][jmult]:
                        leptons = events[lepton][cut][(events[lepton][cut][f"LOOSE{ltype}"] == True)]
                    else:
                        raise ValueError("Not sure what lepton type to choose for event")

                        # get jets and MET
                    jets, met = events["SelectedJets"][cut], events["SelectedMET"][cut]


                        # find best permutations of jets
                    best_perms = ttpermutator.find_best_permutations(jets=jets, leptons=leptons, MET=met, btagWP=btag_wp, btag_req=True)
                    nBestPerms_Total = ak.size(ak.num(best_perms["TTbar"].pt)) # how many events are in signal region
                    output["nBestPerms_Total"].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, nperms=np.zeros(nBestPerms_Total))
                    nBestPerms_Solutions = ak.sum(ak.num(best_perms["TTbar"].pt)) # how many events have a solution from the likelihood
                    output["nBestPerms_Solutions"].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, nperms=np.zeros(nBestPerms_Solutions))

                    valid_perms = ak.num(best_perms["TTbar"].pt) > 0
                    output["cutflow"]["nEvts %s: valid perms" % ", ".join([lepton, leptype, "btagPass", jmult])] += ak.sum(valid_perms)

                    bp_status = np.zeros(cut.size, dtype=int)
                        # find matched permutations
                    mp = ttmatcher.best_match(gen_hyp=events["SL"][cut], jets=jets, leptons=leptons, met=met)

                            ## bp perm categories
                    perm_cat_array = compare_matched_best_perms(mp, best_perms, njets=jmult)
                    bp_status[cut] = perm_cat_array
                    if ak.any(ak.num(events["SL"]["Lepton"].pdgId) != 1): raise ValueError("Number of leptons is incorrect for classifying tau+jets events")
                    sl_tau_evts = ak.where(np.abs(events["SL"]["Lepton"].pdgId) == 15)[0]
                    bp_status[sl_tau_evts] = 4

                    #set_trace()
                    output["ValidBestPerms"].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, permcat=bp_status[cut][valid_perms])

                        # compare best and matched perm objects in valid best perms
                    mp_blep_idx = ak.fill_none(ak.pad_none(mp[valid_perms]["BLep"].jetIdx, 1), -999, axis=1)
                    mp_bhad_idx = ak.fill_none(ak.pad_none(mp[valid_perms]["BHad"].jetIdx, 1), -999, axis=1)
                    mp_wja_idx = ak.fill_none(ak.pad_none(mp[valid_perms]["WJa"].jetIdx, 1), -999, axis=1)
                    mp_wjb_idx = ak.fill_none(ak.pad_none(mp[valid_perms]["WJb"].jetIdx, 1), -999, axis=1)
                    mp_lep_pt = ak.fill_none(ak.pad_none(mp[valid_perms]["Lepton"].pt, 1), -999, axis=1)
                
                    bp_blep_idx = ak.fill_none(ak.pad_none(best_perms[valid_perms]["BLep"].jetIdx, 1), -999, axis=1)
                    bp_bhad_idx = ak.fill_none(ak.pad_none(best_perms[valid_perms]["BHad"].jetIdx, 1), -999, axis=1)
                    bp_wja_idx = ak.fill_none(ak.pad_none(best_perms[valid_perms]["WJa"].jetIdx, 1), -999, axis=1)
                    bp_wjb_idx = ak.fill_none(ak.pad_none(best_perms[valid_perms]["WJb"].jetIdx, 1), -999, axis=1)
                    bp_lep_pt = ak.fill_none(ak.pad_none(best_perms[valid_perms]["Lepton"].pt, 1), -999, axis=1)
                
                        # index comparisons
                    same_blep = (mp_blep_idx == bp_blep_idx) & (mp_blep_idx >= 0)
                    same_bhad = (mp_bhad_idx == bp_bhad_idx) & (mp_bhad_idx >= 0)
                    same_wja = (mp_wja_idx == bp_wja_idx) & (mp_wja_idx >= 0)
                    same_wjb = (mp_wjb_idx == bp_wjb_idx) & (mp_wjb_idx >= 0)

                    isTLepCorrect = ak.flatten(same_blep & (bp_lep_pt == mp_lep_pt))

                    if jmult == "3Jets":
                        valid_evts = (ak.num(mp[valid_perms]["TTbar"].pt) > 0) & (ak.flatten(mp[valid_perms]["unique_matches"] >= 3))
                
                            # merged events
                        merged_evts = valid_evts & ak.flatten(mp[valid_perms]["Merged_Event"])
                        correct_merged = merged_evts & ak.flatten(mp[valid_perms]["Merged_BHadWJa"] | mp[valid_perms]["Merged_BHadWJb"] | mp[valid_perms]["Merged_WJets"])
                        wrong_merged = merged_evts & ak.flatten(~(mp[valid_perms]["Merged_BHadWJa"] | mp[valid_perms]["Merged_BHadWJb"] | mp[valid_perms]["Merged_WJets"]))
                
                            # lost events
                        lost_evts = valid_evts & ak.flatten(mp[valid_perms]["Lost_Event"])
                        correct_lost = lost_evts & ak.flatten(mp[valid_perms]["Lost_WJa"] | mp[valid_perms]["Lost_WJb"])
                        wrong_lost = lost_evts & ak.flatten(~(mp[valid_perms]["Lost_WJa"] | mp[valid_perms]["Lost_WJb"]))

                        right_thad_matching = same_bhad & (((bp_wja_idx == mp_wja_idx) | (bp_wja_idx == mp_wjb_idx)) & (bp_wja_idx >= 0))
                        isTHadCorrect = ak.flatten((correct_lost & right_thad_matching) | (correct_merged & right_thad_matching)) # matched perm is correct event type and right object matching

                    else:
                        #set_trace()
                        isWHadCorrect = ((bp_wja_idx == mp_wja_idx) & (bp_wjb_idx == mp_wjb_idx)) | ((bp_wja_idx == mp_wjb_idx) & (bp_wjb_idx == mp_wja_idx))
                        output["Object_Identifying"].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, objtype="WHad", nID=ak.flatten(isWHadCorrect)*1)
                        isTHadCorrect = ak.flatten(same_bhad & isWHadCorrect)

                    isTTbarCorrect = isTLepCorrect & isTHadCorrect

                    #set_trace()
                    output["Object_Identifying"].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, objtype="BLep", nID=ak.flatten(same_blep)*1)
                    output["Object_Identifying"].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, objtype="BHad", nID=ak.flatten(same_bhad)*1)
                    output["Object_Identifying"].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, objtype="Bs", nID=ak.flatten(same_bhad & same_blep)*1)
                    #output["Object_Identifying"].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, objtype="WJa", nID=ak.to_numpy(ak.flatten(same_wja))*1)
                    #output["Object_Identifying"].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, objtype="WJb", nID=ak.to_numpy(ak.flatten(same_wjb))*1)
                    output["Object_Identifying"].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, objtype="TLep", nID=isTLepCorrect*1)
                    output["Object_Identifying"].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, objtype="THad", nID=isTHadCorrect*1)
                    output["Object_Identifying"].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, objtype="TTbar", nID=isTTbarCorrect*1)


        return output



    def postprocess(self, accumulator):
        return accumulator

proc_executor = processor.iterative_executor if to_debug else processor.futures_executor
proc_exec_args = {"schema": NanoAODSchema} if to_debug else {"schema": NanoAODSchema, "workers": 8}
output = processor.run_uproot_job(
    fileset,
    treename="Events",
    processor_instance=ttbar_evt_cats(),
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
