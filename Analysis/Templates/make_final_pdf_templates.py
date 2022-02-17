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
#import Utilities.systematics as systematics
#import bkg_systs_parameters as bkg_syspar   

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
args = parser.parse_args()


def final_bkg_templates(fnames):
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
                proc = tname.split("_")[0] if not "data_obs" in tname else "data_obs"
                sys = sorted(filter(None, tname.split(f"{proc}_")))[0]
                print(lep, jmult, sys, proc)

                if sys == "nosys":
                    treatment = "raw"
                else:
                    ## choose which hist to use based on what_to_do dict in systs_parameters.py
                    #sysname = systematics.sys_to_name[args.year][sys] # convert sys to name used in analysis
                    #sys_to_check = [key for key in bkg_syspar.what_to_do.keys() if fnmatch.fnmatch(sysname, f"{key}*")][0] # find sysname which corresponds to name in systs_parameters.py
                    treatment = "raw"

                if treatment == "raw":
                    hist_to_use = raw_hdict[jmult][lep][tname].copy()
                elif treatment == "flat":
                    hist_to_use = flat_hdict[jmult][lep][tname].copy()
                elif treatment == "smooth":
                    hist_to_use = smooth_hdict[jmult][lep][tname].copy()
                else:
                    raise ValueError(f"Treatment {treatment} not supported.")

                    ## save template histos to coffea dict
                histo_dict[jmult][lep][tname] = [hist_to_use.copy(), treatment]

    coffea_out = os.path.join(outdir, f"final_pdf_templates_lj_bkg_{args.year}_{jobid}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")


if __name__ == "__main__":
    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]

    njets_to_run = ["3Jets", "4PJets"]

    analyzer = "htt_pdfUncs"
    base_bkg_template_name = f"*_pdf_templates_lj_bkg_{args.year}_{jobid}.coffea"

    sys_treatment_dict = {
        "raw" : "raw", 
        "flat" : "flattened",
        "smooth" : "smoothed",
    }

    input_dir = os.path.join(proj_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}")
    outdir = os.path.join(proj_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}", "FINAL")
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

        # find background files
    bkg_fnames = fnmatch.filter(os.listdir(input_dir), base_bkg_template_name)
    bkg_fnames = [os.path.join(input_dir, fname) for fname in bkg_fnames]
    print("Creating final background templates")
    final_bkg_templates(bkg_fnames)

    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
