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
    
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
args = parser.parse_args()


def flatten_bkg_templates(fname):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """
    histo_dict = processor.dict_accumulator({njets : {"Muon" : {}, "Electron" :{}} for njets in njets_to_run})
    hdict = load(fname)
    for jmult in hdict.keys():
        for lep in hdict[jmult].keys():
            for tname, orig_template in hdict[jmult][lep].items():
                proc = tname.split("_")[0] if not "data_obs" in tname else "data_obs"
                sys = sorted(filter(None, tname.split(f"{proc}_")))[0]
                print(lep, jmult, sys, proc)

                    # perform flattening
                flattened_histo = hdict[jmult][lep][f"{proc}_nosys"].copy() if sys == "nosys" else Plotter.flatten(nosys=hdict[jmult][lep][f"{proc}_nosys"].copy(), systematic=orig_template.copy())
                
                    ## save template histos to coffea dict
                histo_dict[jmult][lep][tname] = flattened_histo.copy()

    coffea_out = os.path.join(input_dir, f"flattened_pdf_templates_lj_bkg_{args.year}_{jobid}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")


if __name__ == "__main__":
    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]

    njets_to_run = ["3Jets", "4PJets"]

    analyzer = "htt_pdfUncs"
    input_dir = os.path.join(proj_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}")
    bkg_fname = os.path.join(input_dir, f"raw_pdf_templates_lj_bkg_{args.year}_{jobid}.coffea")
    if not os.path.isfile(bkg_fname): raise ValueError("No background file found.")
    print("Creating flattened background templates")
    flatten_bkg_templates(bkg_fname)

    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
