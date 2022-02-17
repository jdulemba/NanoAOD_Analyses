#! /bin/env python

import time
tic = time.time()

from coffea.util import load, save
from pdb import set_trace
import os
import numpy as np
import fnmatch
import coffea.processor as processor
import Utilities.systematics as systematics
import bkg_systs_parameters as bkg_syspar   
import signal_systs_parameters as sig_syspar
import Utilities.prettyjson as prettyjson

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("--only_bkg", action="store_true", help="Make background templates only.")
parser.add_argument("--only_sig", action="store_true", help="Make signal templates only.")
parser.add_argument("--kfactors", action="store_true", help="Apply signal k-factors to signal")
parser.add_argument("--scale_mtop3gev", action="store_true", help="Scale 3GeV mtop variations by 1/6")
args = parser.parse_args()


def final_bkg_templates(fnames):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """
    histo_dict = processor.dict_accumulator({njets : {"Muon" : {}, "Electron" :{}} for njets in njets_to_run})

    raw_hdict = load([fname for fname in fnames if sys_treatment_dict["raw"] in os.path.basename(fname)][0])
    smooth_hdict = load([fname for fname in fnames if sys_treatment_dict["smooth"] in os.path.basename(fname)][0])
    flat_hdict = load([fname for fname in fnames if sys_treatment_dict["flat"] in os.path.basename(fname)][0])
    symm_hdict = load([fname for fname in fnames if sys_treatment_dict["symm"] in os.path.basename(fname)][0])
    to_symm = prettyjson.loads(open(os.path.join(input_dir, f"templates_to_symmetrize_lj_bkg_mtopscaled_{args.year}_{jobid}.json" if args.scale_mtop3gev else f"templates_to_symmetrize_lj_bkg_{args.year}_{jobid}.json")).read())
    for jmult in raw_hdict.keys():
        for lep in raw_hdict[jmult].keys():
            for tname in raw_hdict[jmult][lep].keys():
                proc = tname.split("_")[0] if not "data_obs" in tname else "data_obs"
                sys = sorted(filter(None, tname.split(f"{proc}_")))[0]
                print(jmult, lep, sys, proc)

                if sys == "nosys":
                    treatment = "raw"
                else:
                    # choose which hist to use based on what_to_do dict in systs_parameters.py
                    sysname = systematics.sys_to_name[args.year][sys] if sys == "EWcorrUp" else "_".join(systematics.sys_to_name[args.year][sys].split("_")[:-1]) # convert sys to name used in analysis
                    treatment = bkg_syspar.what_to_do[sysname][proc][f"{lep.lower()[0]}{jmult[0]}"]

                    # check if symmetrization should be applied
                    if f"{proc}_{sysname}" in to_symm[jmult][lep].keys(): treatment = "symmetrized"

                print(f"\t{treatment}")
                if treatment == "raw":
                    hist_to_use = raw_hdict[jmult][lep][tname].copy()
                elif treatment == "flat":
                    hist_to_use = flat_hdict[jmult][lep][tname].copy()
                elif treatment == "smooth":
                    hist_to_use = smooth_hdict[jmult][lep][tname].copy()
                elif treatment == "symmetrized":
                    hist_to_use = symm_hdict[jmult][lep][tname].copy()
                else:
                    raise ValueError(f"Treatment {treatment} not supported.")

                    ## save template histos to coffea dict
                histo_dict[jmult][lep][tname] = [hist_to_use.copy(), treatment]

    #set_trace()
    coffea_out = os.path.join(outdir, f"final_templates_lj_bkg_mtopscaled_{args.year}_{jobid}.coffea" if args.scale_mtop3gev else f"final_templates_lj_bkg_{args.year}_{jobid}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")


def final_sig_templates(fnames):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """
    histo_dict = processor.dict_accumulator({njets : {"Muon" : {}, "Electron" :{}} for njets in njets_to_run})

    raw_hdict = load([fname for fname in fnames if sys_treatment_dict["raw"] in os.path.basename(fname)][0])
    smooth_hdict = load([fname for fname in fnames if sys_treatment_dict["smooth"] in os.path.basename(fname)][0])
    flat_hdict = load([fname for fname in fnames if sys_treatment_dict["flat"] in os.path.basename(fname)][0])
    for jmult in raw_hdict.keys():
        for lep in raw_hdict[jmult].keys():
            for tname in raw_hdict[jmult][lep].keys():
                if "Int" in tname:
                    proc = tname.split("pos_")[0]+"pos" if "pos" in tname else tname.split("neg_")[0]+"neg"
                else:
                    proc = tname.split("Res_")[0]+"Res"
                sys = sorted(filter(None, tname.split(f"{proc}_")))[0]
                print(lep, jmult, sys, proc)

                if args.kfactors:
                    treatment = "raw"
                else:
                    if sys == "nosys":
                        treatment = "raw"
                    else:
                        # choose which hist to use based on what_to_do dict in systs_parameters.py
                        sysname = systematics.sys_to_name[args.year][sys] # convert sys to name used in analysis
                        sys_to_check = [key for key in sig_syspar.what_to_do.keys() if fnmatch.fnmatch(sysname, f"{key}*")][0] # find sysname which corresponds to name in systs_parameters.py
                        treatment = sig_syspar.what_to_do[sys_to_check]#[f"{lep.lower()[0]}{jmult[0]}"]
                        #treatment = sig_syspar.what_to_do[sys_to_check][proc][f"{lep.lower()[0]}{jmult[0]}"]

                print(f"\t{treatment}")
                if treatment == "raw":
                    hist_to_use = raw_hdict[jmult][lep][tname].copy()
                if treatment == "flat":
                    hist_to_use = flat_hdict[jmult][lep][tname].copy()
                if treatment == "smooth":
                    hist_to_use = smooth_hdict[jmult][lep][tname].copy()

                    ## save template histos to coffea dict
                histo_dict[jmult][lep][tname] = [hist_to_use.copy(), treatment]

    coffea_out = os.path.join(outdir, f"final_templates_lj_sig_kfactors_{args.year}_{jobid}.coffea" if args.kfactors else f"final_templates_lj_sig_{args.year}_{jobid}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")


if __name__ == "__main__":
    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]

    njets_to_run = ["3Jets", "4PJets"]

    analyzer = "htt_btag_sb_regions"
    base_bkg_template_name = f"*_templates_lj_bkg_mtopscaled_{args.year}_{jobid}.coffea" if args.scale_mtop3gev else f"*_templates_lj_bkg_{args.year}_{jobid}.coffea"
    base_sig_template_name = f"*_templates_lj_sig_kfactors_{args.year}_{jobid}.coffea" if args.kfactors else f"*_templates_lj_sig_{args.year}_{jobid}.coffea"

    sys_treatment_dict = {
        "raw" : "raw", 
        "flat" : "flattened",
        "smooth" : "smoothed",
        "symm" : "symmetrized",
    }

    input_dir = os.path.join(proj_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}")
    outdir = os.path.join(proj_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}", "FINAL")
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    if not args.only_sig:
            # find background files
        bkg_fnames = fnmatch.filter(os.listdir(input_dir), base_bkg_template_name)
        bkg_fnames = [os.path.join(input_dir, fname) for fname in bkg_fnames]
        print("Creating final background templates")
        final_bkg_templates(bkg_fnames)

    if not args.only_bkg:
            # find signal files
        sig_fnames = fnmatch.filter(os.listdir(input_dir), base_sig_template_name)
        sig_fnames = [os.path.join(input_dir, fname) for fname in sig_fnames]
        print("Creating final signal templates")
        final_sig_templates(sig_fnames)

    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
