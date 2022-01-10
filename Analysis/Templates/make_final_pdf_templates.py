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
import Utilities.systematics as systematics
import bkg_systs_parameters as bkg_syspar   
import signal_systs_parameters as sig_syspar

base_jobid = os.environ["base_jobid"]
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("--njets", default="all", nargs="?", choices=["3", "4+", "all"], help="Specify which jet multiplicity to use.")
args = parser.parse_args()

njets_to_run = []
if (args.njets == "3") or (args.njets == "all"):
    njets_to_run += ["3Jets"]
if (args.njets == "4+") or (args.njets == "all"):
    njets_to_run += ["4PJets"]


def final_bkg_templates(bkg_dict):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """
    if "3Jets" in njets_to_run:
        histo_dict_3j = processor.dict_accumulator({"Muon" : {}, "Electron" :{}})
    if "4PJets" in njets_to_run:
        histo_dict_4pj = processor.dict_accumulator({"Muon" : {}, "Electron" :{}})

    for jmult in njets_to_run:
        fnames = bkg_dict[jmult]
        raw_hdict = load([fname for fname in fnames if sys_treatment_dict["raw"] in os.path.basename(fname)][0])
        smooth_hdict = load([fname for fname in fnames if sys_treatment_dict["smooth"] in os.path.basename(fname)][0])
        flat_hdict = load([fname for fname in fnames if sys_treatment_dict["flat"] in os.path.basename(fname)][0])

        for lep in raw_hdict.keys():
            for tname in raw_hdict[lep].keys():
                proc = tname.split("_")[0] if not "data_obs" in tname else "data_obs"
                sys = sorted(filter(None, tname.split(f"{proc}_")))[0]
                print(lep, jmult, sys, proc)

                if sys == "nosys":
                    hist_to_use = raw_hdict[lep][tname].copy()
                    treatment = "raw"

                else:
                    ## choose which hist to use based on what_to_do dict in systs_parameters.py
                    #sysname = systematics.sys_to_name[args.year][sys] # convert sys to name used in analysis
                    #sys_to_check = [key for key in bkg_syspar.what_to_do.keys() if fnmatch.fnmatch(sysname, f"{key}*")][0] # find sysname which corresponds to name in systs_parameters.py

                    #treatment = bkg_syspar.what_to_do[sys_to_check][proc][f"{lep.lower()[0]}{jmult[0]}"]
                    treatment = "raw"
                    if treatment == "raw":
                        hist_to_use = raw_hdict[lep][tname].copy()
                    if treatment == "flat":
                        hist_to_use = flat_hdict[lep][tname].copy()
                    if treatment == "smooth":
                        hist_to_use = smooth_hdict[lep][tname].copy()

                    ## save template histos to coffea dict
                if jmult == "3Jets":
                    histo_dict_3j[lep][tname] = [hist_to_use.copy(), treatment]
                if jmult == "4PJets":
                    histo_dict_4pj[lep][tname] = [hist_to_use.copy(), treatment]

    #set_trace()
    if "3Jets" in njets_to_run:
        coffea_out_3j = os.path.join(outdir, f"final_pdf_templates_lj_3Jets_bkg_{args.year}_{jobid}.coffea")
        save(histo_dict_3j, coffea_out_3j)
        print(f"{coffea_out_3j} written")
    if "4PJets" in njets_to_run:
        coffea_out_4pj = os.path.join(outdir, f"final_pdf_templates_lj_4PJets_bkg_{args.year}_{jobid}.coffea")
        save(histo_dict_4pj, coffea_out_4pj)
        print(f"{coffea_out_4pj} written")



if __name__ == "__main__":
    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]

    analyzer = "htt_pdfUncs"
    base_bkg_template_name = f"*_pdf_templates_lj_NJETS_bkg_{args.year}_{jobid}"

    sys_treatment_dict = {
        "raw" : "raw", 
        "flat" : "flattened",
        "smooth" : "smoothed",
    }

    input_dir = os.path.join(proj_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}")
    if os.path.isdir(input_dir):
        # define variables to get histogram for background  
            # find files for 3 jets
        bkg_3j_fnmatch = "%s.coffea" % base_bkg_template_name.replace("NJETS", "3Jets")
        bkg_3j_fnames = fnmatch.filter(os.listdir(input_dir), bkg_3j_fnmatch)
        bkg_3j_fnames = [os.path.join(input_dir, fname) for fname in bkg_3j_fnames]
            # find files for 4+ jets
        bkg_4pj_fnmatch = "%s.coffea" % base_bkg_template_name.replace("NJETS", "4PJets")
        bkg_4pj_fnames = fnmatch.filter(os.listdir(input_dir), bkg_4pj_fnmatch)
        bkg_4pj_fnames = [os.path.join(input_dir, fname) for fname in bkg_4pj_fnames]
        bkg_dict = {"3Jets" : bkg_3j_fnames, "4PJets" : bkg_4pj_fnames}
    else: print("No background file found.")

    outdir = os.path.join(proj_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}", "FINAL")
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    print("Creating final background templates")
    final_bkg_templates(bkg_dict)

    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
