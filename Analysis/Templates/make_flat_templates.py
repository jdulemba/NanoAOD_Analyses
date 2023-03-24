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
parser.add_argument("templates_to_run", type=str, help="Choose which type of templates to run, multiple options can be input as ':' separated strings.")
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
                #if sys == "nosys": continue
                print(lep, jmult, sys, proc)

                    # perform flattening
                flattened_histo = hdict[jmult][lep][f"{proc}_nosys"].copy() if sys == "nosys" else Plotter.flatten(nosys=hdict[jmult][lep][f"{proc}_nosys"].copy(), systematic=orig_template.copy())
                
                    ## save template histos to coffea dict
                histo_dict[jmult][lep][tname] = flattened_histo.copy()

    #set_trace()
    coffea_out = os.path.join(input_dir, f"flattened_templates_lj_bkg_{args.year}_{jobid}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")


def flatten_MEsig_templates(fname):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """
    histo_dict = processor.dict_accumulator({njets : {"Muon" : {}, "Electron" :{}} for njets in njets_to_run})
    hdict = load(fname)
    for jmult in hdict.keys():
        for lep in hdict[jmult].keys():
            for tname, orig_template in hdict[jmult][lep].items():
                if "Res" in tname:
                    signal = "_".join([tname.split("_Res")[0], "Res"])
                else:
                    signal = "_".join([tname.split("_neg")[0], "neg"]) if "neg" in tname else "_".join([tname.split("_pos")[0], "pos"])
                sys = sorted(filter(None, tname.split(f"{signal}_")))[0]
                print(lep, jmult, sys, signal)                

                    # perform flattening
                flattened_histo = hdict[jmult][lep][f"{signal}_nosys"].copy() if sys == "nosys" else Plotter.flatten(nosys=hdict[jmult][lep][f"{signal}_nosys"].copy(), systematic=orig_template.copy())
                
                    ## save template histos to coffea dict
                histo_dict[jmult][lep][tname] = flattened_histo.copy()

    #set_trace()
    coffea_out = os.path.join(input_dir, f"flattened_templates_lj_MEsig_{args.year}_{jobid}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")


def flatten_pdf_templates(fname):
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
    allowed_template_options = ["bkg", "MEreweight_sig", "PDF"]
    templates_to_run = (args.templates_to_run).split(":")
    templates_to_run = [template for template in templates_to_run if template in allowed_template_options]

    jobid = os.environ["jobid"]
    eos_dir = os.environ["eos_dir"]

    njets_to_run = ["3Jets", "4PJets"]

    analyzer = "htt_btag_sb_regions"

    input_dir = os.path.join(eos_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}")
    if not os.path.isdir(input_dir): raise ValueError("No directory found.")

    if "bkg" in templates_to_run:
                # define variables to get histogram for background
        base_bkg_template_name = f"raw_templates_lj_bkg_{args.year}_{jobid}.coffea"
        bkg_fname = os.path.join(input_dir, base_bkg_template_name)
        if not os.path.isfile(bkg_fname): raise ValueError("No background file found.")
        print("Creating flattened background templates")
        flatten_bkg_templates(bkg_fname)

    if "MEreweight_sig" in templates_to_run:
        base_sig_template_name = f"raw_templates_lj_MEsig_{args.year}_{jobid}.coffea"
        sig_fname = os.path.join(input_dir, base_sig_template_name)
        if not os.path.isfile(sig_fname): raise ValueError("No signal file found.")
        print("Creating ME reweighting signal templates")
        flatten_MEsig_templates(sig_fname)

    if "PDF" in templates_to_run:
        analyzer = "htt_pdfUncs"

        input_dir = os.path.join(proj_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}")
        bkg_fname = os.path.join(input_dir, f"raw_pdf_templates_lj_bkg_{args.year}_{jobid}.coffea")
        if not os.path.isfile(bkg_fname): raise ValueError("No background file found.")
        print("Creating flattened background templates")
        flatten_pdf_templates(bkg_fname)


    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
