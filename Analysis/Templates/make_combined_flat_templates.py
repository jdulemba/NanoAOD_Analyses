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
parser.add_argument("templates_to_run", type=str, help="Choose which type of templates to run, multiple options can be input as ':' separated strings.")
parser.add_argument("channels_to_combine", type=str, help="Choose how to combine systematic templates (by lepton or eras), multiple options can be input as ':' separated strings.")
parser.add_argument("--nomSMTTxsec", action="store_true", help="Apply nominal SM cross sections to top mass and LHE scale weights")
args = parser.parse_args()


def flatten_year_and_lepton_templates(fname, process):
    histo_dict = processor.dict_accumulator({njets : {} for njets in njets_to_run})
    hdict = load(fname)
    for jmult in hdict.keys():
        for tname, orig_template in hdict[jmult].items():
            if process == "bkg":
                proc = tname.split("_")[0] if not "data_obs" in tname else "data_obs"
            else:
                if "Res" in tname:
                    proc = "_".join([tname.split("_Res")[0], "Res"])
                else:
                    proc = "_".join([tname.split("_neg")[0], "neg"]) if "neg" in tname else "_".join([tname.split("_pos")[0], "pos"])

            sys = sorted(filter(None, tname.split(f"{proc}_")))[0]
            print(jmult, sys, proc)

                # perform flattening
            flattened_histo = hdict[jmult][f"{proc}_nosys"].copy() if sys == "nosys" else Plotter.flatten(nosys=hdict[jmult][f"{proc}_nosys"].copy(), systematic=orig_template.copy(), ratio=True)
            
                ## save template histos to coffea dict
            histo_dict[jmult][tname] = flattened_histo.copy()

    #set_trace()
    coffea_out = os.path.join(input_dir, f"flattened_combined_year_and_lepton_templates_lj_{process}_nomSMTTxsec_{jobid}.coffea" if args.nomSMTTxsec else f"flattened_combined_year_and_lepton_templates_lj_{process}_{jobid}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")


def flatten_lepton_templates(fname, process):
    histo_dict = processor.dict_accumulator({njets : {} for njets in njets_to_run})
    hdict = load(fname)
    for jmult in hdict.keys():
        for tname, orig_template in hdict[jmult].items():
            if process == "bkg":
                proc = tname.split("_")[0] if not "data_obs" in tname else "data_obs"
            else:
                if "Res" in tname:
                    proc = "_".join([tname.split("_Res")[0], "Res"])
                else:
                    proc = "_".join([tname.split("_neg")[0], "neg"]) if "neg" in tname else "_".join([tname.split("_pos")[0], "pos"])

            sys = sorted(filter(None, tname.split(f"{proc}_")))[0]
            print(year, jmult, sys, proc)

                # perform flattening
            flattened_histo = hdict[jmult][f"{proc}_nosys"].copy() if sys == "nosys" else Plotter.flatten(nosys=hdict[jmult][f"{proc}_nosys"].copy(), systematic=orig_template.copy(), ratio=True)
            
                ## save template histos to coffea dict
            histo_dict[jmult][tname] = flattened_histo.copy()

    #set_trace()
    coffea_out = os.path.join(input_dir, f"flattened_combined_lep_templates_lj_{process}_{year}_{jobid}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")



if __name__ == "__main__":
    allowed_template_options = ["bkg"]
    templates_to_run = (args.templates_to_run).split(":")
    templates_to_run = [template for template in templates_to_run if template in allowed_template_options]

    allowed_combination_options = ["lepton", "era_lepton"]
    combinations_to_run = (args.channels_to_combine).split(":")
    combinations_to_run = [combination for combination in combinations_to_run if combination in allowed_combination_options]

    jobid = os.environ["jobid"]
    eos_dir = os.environ["eos_dir"]

    njets_to_run = ["3Jets", "4PJets"]
    analyzer = "htt_btag_sb_regions"

    for combination in combinations_to_run:
        if combination == "lepton":
            for year in ["2016APV", "2016", "2017", "2018"]:
                input_dir = os.path.join(eos_dir, "results", f"{year}_{jobid}", f"Templates_{analyzer}")
                if not os.path.isdir(input_dir): raise ValueError("No background file found.")

                if "bkg" in templates_to_run:
                    base_bkg_template_name = f"raw_combined_lep_templates_lj_bkg_{year}_{jobid}.coffea"
                    bkg_fname = os.path.join(input_dir, base_bkg_template_name)
                    if not os.path.isfile(bkg_fname): raise ValueError(f"{bkg_fname} not found.")
                    print(f"Flattening e+mu channels in {year} for background templates")
                    flatten_lepton_templates(bkg_fname, "bkg")

        if combination == "era_lepton":
            input_dir = os.path.join(eos_dir, "results", jobid, f"Templates_{analyzer}")
            if not os.path.isdir(input_dir): raise ValueError("No background file found.")

            if "bkg" in templates_to_run:
                bkg_fname = os.path.join(input_dir, f"raw_combined_year_and_lepton_templates_lj_bkg_nomSMTTxsec_{jobid}.coffea" if args.nomSMTTxsec \
                    else f"raw_combined_year_and_lepton_templates_lj_bkg_{jobid}.coffea")
                if not os.path.isfile(bkg_fname): raise ValueError(f"{bkg_fname} not found.")
                print("Flattening combined channels across years and leptons for background templates")
                flatten_year_and_lepton_templates(bkg_fname, "bkg")


    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
