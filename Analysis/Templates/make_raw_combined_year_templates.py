#! /bin/env python

import time
tic = time.time()

from coffea.util import load, save
from pdb import set_trace
import os
import numpy as np
import fnmatch
import Utilities.Plotter as Plotter
import coffea.processor as processor    
    
base_jobid = os.environ["base_jobid"]
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--only_bkg", action="store_true", help="Make background templates only.")
parser.add_argument("--only_sig", action="store_true", help="Make signal templates only.")
args = parser.parse_args()


def combine_bkg_templates(fnames_dict):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """
    histo_dict_3j = processor.dict_accumulator({"Muon" : {}, "Electron" : {}})
    histo_dict_4pj = processor.dict_accumulator({"Muon" : {}, "Electron" : {}})

        # get Muon, 3Jets hists for 2016APV
    mu_3j_16APV_dict = fnames_dict["2016APV"]["3Jets"]["Muon"]
    for tname, mu_3j_16APV_template in mu_3j_16APV_dict.items():
        proc = tname.split("_")[0] if not "data_obs" in tname else "data_obs"
        sys = sorted(filter(None, tname.split(f"{proc}_")))[0]

        print(sys, proc)
        empty_hist = mu_3j_16APV_template.copy()
        empty_hist.clear()

            # muon, 3 jets
        mu_3j_16_template = fnames_dict["2016"]["3Jets"]["Muon"][tname] if tname in fnames_dict["2016"]["3Jets"]["Muon"].keys() else empty_hist
        mu_3j_17_template = fnames_dict["2017"]["3Jets"]["Muon"][tname] if tname in fnames_dict["2017"]["3Jets"]["Muon"].keys() else empty_hist
        mu_3j_18_template = fnames_dict["2018"]["3Jets"]["Muon"][tname] if tname in fnames_dict["2018"]["3Jets"]["Muon"].keys() else empty_hist
        mu_3j_template = mu_3j_16APV_template.add(mu_3j_16_template).add(mu_3j_17_template).add(mu_3j_18_template)
        histo_dict_3j["Muon"][tname] = mu_3j_template.copy()

            # electron, 3 jets
        el_3j_16APV_template = fnames_dict["2016APV"]["3Jets"]["Electron"][tname] if tname in fnames_dict["2016APV"]["3Jets"]["Electron"].keys() else empty_hist
        el_3j_16_template = fnames_dict["2016"]["3Jets"]["Electron"][tname] if tname in fnames_dict["2016"]["3Jets"]["Electron"].keys() else empty_hist
        el_3j_17_template = fnames_dict["2017"]["3Jets"]["Electron"][tname] if tname in fnames_dict["2017"]["3Jets"]["Electron"].keys() else empty_hist
        el_3j_18_template = fnames_dict["2018"]["3Jets"]["Electron"][tname] if tname in fnames_dict["2018"]["3Jets"]["Electron"].keys() else empty_hist
        el_3j_template = el_3j_16APV_template.add(el_3j_16_template).add(el_3j_17_template).add(el_3j_18_template)
        histo_dict_3j["Electron"][tname] = el_3j_template.copy()

            # muon, 4+ jets
        mu_4pj_16APV_template = fnames_dict["2016APV"]["4PJets"]["Muon"][tname] if tname in fnames_dict["2016APV"]["4PJets"]["Muon"].keys() else empty_hist
        mu_4pj_16_template = fnames_dict["2016"]["4PJets"]["Muon"][tname] if tname in fnames_dict["2016"]["4PJets"]["Muon"].keys() else empty_hist
        mu_4pj_17_template = fnames_dict["2017"]["4PJets"]["Muon"][tname] if tname in fnames_dict["2017"]["4PJets"]["Muon"].keys() else empty_hist
        mu_4pj_18_template = fnames_dict["2018"]["4PJets"]["Muon"][tname] if tname in fnames_dict["2018"]["4PJets"]["Muon"].keys() else empty_hist
        mu_4pj_template = mu_4pj_16APV_template.add(mu_4pj_16_template).add(mu_4pj_17_template).add(mu_4pj_18_template)
        histo_dict_4pj["Muon"][tname] = mu_4pj_template.copy()

            # electron, 4+ jets
        el_4pj_16APV_template = fnames_dict["2016APV"]["4PJets"]["Electron"][tname] if tname in fnames_dict["2016APV"]["4PJets"]["Electron"].keys() else empty_hist
        el_4pj_16_template = fnames_dict["2016"]["4PJets"]["Electron"][tname] if tname in fnames_dict["2016"]["4PJets"]["Electron"].keys() else empty_hist
        el_4pj_17_template = fnames_dict["2017"]["4PJets"]["Electron"][tname] if tname in fnames_dict["2017"]["4PJets"]["Electron"].keys() else empty_hist
        el_4pj_18_template = fnames_dict["2018"]["4PJets"]["Electron"][tname] if tname in fnames_dict["2018"]["4PJets"]["Electron"].keys() else empty_hist
        el_4pj_template = el_4pj_16APV_template.add(el_4pj_16_template).add(el_4pj_17_template).add(el_4pj_18_template)
        histo_dict_4pj["Electron"][tname] = el_4pj_template.copy()

    coffea_out_3j = os.path.join(outdir, f"test_CombinedRaw_templates_lj_3Jets_bkg_{jobid}.coffea")
    save(histo_dict_3j, coffea_out_3j)
    print(f"{coffea_out_3j} written")

    coffea_out_4pj = os.path.join(outdir, f"test_CombinedRaw_templates_lj_4PJets_bkg_{jobid}.coffea")
    save(histo_dict_4pj, coffea_out_4pj)
    print(f"{coffea_out_4pj} written")


def combine_sig_templates(fnames_dict):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """
    histo_dict_3j = processor.dict_accumulator({"Muon" : {}, "Electron" :{}})
    histo_dict_4pj = processor.dict_accumulator({"Muon" : {}, "Electron" :{}})

        # get Muon, 3Jets hists for 2016APV
    mu_3j_16APV_dict = fnames_dict["2016APV"]["3Jets"]["Muon"]
    for tname, mu_3j_16APV_template in mu_3j_16APV_dict.items():
        if "Res" in tname:
            signal = "_".join([tname.split("_Res")[0], "Res"])
        else:
            signal = "_".join([tname.split("_neg")[0], "neg"]) if "neg" in tname else "_".join([tname.split("_pos")[0], "pos"])
        sys = sorted(filter(None, tname.split(f"{signal}_")))[0]
        print(sys, signal)                

        empty_hist = mu_3j_16APV_template.copy()
        empty_hist.clear()

            # muon, 3 jets
        mu_3j_16_template = fnames_dict["2016"]["3Jets"]["Muon"][tname] if tname in fnames_dict["2016"]["3Jets"]["Muon"].keys() else empty_hist
        mu_3j_17_template = fnames_dict["2017"]["3Jets"]["Muon"][tname] if tname in fnames_dict["2017"]["3Jets"]["Muon"].keys() else empty_hist
        mu_3j_18_template = fnames_dict["2018"]["3Jets"]["Muon"][tname] if tname in fnames_dict["2018"]["3Jets"]["Muon"].keys() else empty_hist
        mu_3j_template = mu_3j_16APV_template.add(mu_3j_16_template).add(mu_3j_17_template).add(mu_3j_18_template)
        histo_dict_3j["Muon"][tname] = mu_3j_template.copy()

            # electron, 3 jets
        el_3j_16APV_template = fnames_dict["2016APV"]["3Jets"]["Electron"][tname] if tname in fnames_dict["2016APV"]["3Jets"]["Electron"].keys() else empty_hist
        el_3j_16_template = fnames_dict["2016"]["3Jets"]["Electron"][tname] if tname in fnames_dict["2016"]["3Jets"]["Electron"].keys() else empty_hist
        el_3j_17_template = fnames_dict["2017"]["3Jets"]["Electron"][tname] if tname in fnames_dict["2017"]["3Jets"]["Electron"].keys() else empty_hist
        el_3j_18_template = fnames_dict["2018"]["3Jets"]["Electron"][tname] if tname in fnames_dict["2018"]["3Jets"]["Electron"].keys() else empty_hist
        el_3j_template = el_3j_16APV_template.add(el_3j_16_template).add(el_3j_17_template).add(el_3j_18_template)
        histo_dict_3j["Electron"][tname] = el_3j_template.copy()

            # muon, 4+ jets
        mu_4pj_16APV_template = fnames_dict["2016APV"]["4PJets"]["Muon"][tname] if tname in fnames_dict["2016APV"]["4PJets"]["Muon"].keys() else empty_hist
        mu_4pj_16_template = fnames_dict["2016"]["4PJets"]["Muon"][tname] if tname in fnames_dict["2016"]["4PJets"]["Muon"].keys() else empty_hist
        mu_4pj_17_template = fnames_dict["2017"]["4PJets"]["Muon"][tname] if tname in fnames_dict["2017"]["4PJets"]["Muon"].keys() else empty_hist
        mu_4pj_18_template = fnames_dict["2018"]["4PJets"]["Muon"][tname] if tname in fnames_dict["2018"]["4PJets"]["Muon"].keys() else empty_hist
        mu_4pj_template = mu_4pj_16APV_template.add(mu_4pj_16_template).add(mu_4pj_17_template).add(mu_4pj_18_template)
        histo_dict_4pj["Muon"][tname] = mu_4pj_template.copy()

            # electron, 4+ jets
        el_4pj_16APV_template = fnames_dict["2016APV"]["4PJets"]["Electron"][tname] if tname in fnames_dict["2016APV"]["4PJets"]["Electron"].keys() else empty_hist
        el_4pj_16_template = fnames_dict["2016"]["4PJets"]["Electron"][tname] if tname in fnames_dict["2016"]["4PJets"]["Electron"].keys() else empty_hist
        el_4pj_17_template = fnames_dict["2017"]["4PJets"]["Electron"][tname] if tname in fnames_dict["2017"]["4PJets"]["Electron"].keys() else empty_hist
        el_4pj_18_template = fnames_dict["2018"]["4PJets"]["Electron"][tname] if tname in fnames_dict["2018"]["4PJets"]["Electron"].keys() else empty_hist
        el_4pj_template = el_4pj_16APV_template.add(el_4pj_16_template).add(el_4pj_17_template).add(el_4pj_18_template)
        histo_dict_4pj["Electron"][tname] = el_4pj_template.copy()

    coffea_out_3j = os.path.join(outdir, f"test_CombinedRaw_templates_lj_3Jets_sig_{jobid}.coffea")
    save(histo_dict_3j, coffea_out_3j)
    print(f"{coffea_out_3j} written")

    coffea_out_4pj = os.path.join(outdir, f"test_CombinedRaw_templates_lj_4PJets_sig_{jobid}.coffea")
    save(histo_dict_4pj, coffea_out_4pj)
    print(f"{coffea_out_4pj} written")



if __name__ == "__main__":
    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]

    analyzer = "htt_btag_sb_regions"
    base_bkg_template_name = f"test_raw_templates_lj_NJETS_bkg_YEAR_{jobid}"
    base_sig_template_name = f"test_raw_templates_lj_NJETS_sig_YEAR_{jobid}"

    years_to_run = ["2016APV", "2016", "2017", "2018"]

    #set_trace()
    bkg_fnames_dict = {year : {} for year in years_to_run}
    sig_fnames_dict = {year : {} for year in years_to_run}

    outdir = os.path.join(proj_dir, "results", jobid, f"Templates_{analyzer}")
    if not os.path.isdir(outdir): os.makedirs(outdir)

    for year in years_to_run:
        input_dir = os.path.join(proj_dir, "results", f"{year}_{jobid}", f"Templates_{analyzer}")
        if os.path.isdir(input_dir):
                # define variables to get histogram for background
            bkg_fnames_dict[year].update({"3Jets" : load(os.path.join(input_dir, "%s.coffea" % base_bkg_template_name.replace("NJETS", "3Jets").replace("YEAR", year)))})
            bkg_fnames_dict[year].update({"4PJets" : load(os.path.join(input_dir, "%s.coffea" % base_bkg_template_name.replace("NJETS", "4PJets").replace("YEAR", year)))})
                # define variables to get histogram for signal
            sig_fnames_dict[year].update({"3Jets" : load(os.path.join(input_dir, "%s.coffea" % base_sig_template_name.replace("NJETS", "3Jets").replace("YEAR", year)))})
            sig_fnames_dict[year].update({"4PJets" : load(os.path.join(input_dir, "%s.coffea" % base_sig_template_name.replace("NJETS", "4PJets").replace("YEAR", year)))})
        else: print("No files found.")

    if not args.only_sig:
        try:
            print("Creating flattened background templates")
            combine_bkg_templates(bkg_fnames_dict)
        except:
            print("Could not write background templates to file")

    if not args.only_bkg:
        try:
            print("Creating flattened signal templates")
            combine_sig_templates(sig_fnames_dict)
        except:
            print("Could not write signal templates to file")

    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
