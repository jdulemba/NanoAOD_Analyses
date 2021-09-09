#! /bin/env python

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
parser.add_argument("year", choices=["2016", "2017", "2018"] if base_jobid == "NanoAODv6" else ["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("--njets", default="all", nargs="?", choices=["3", "4+", "all"], help="Specify which jet multiplicity to use.")
args = parser.parse_args()

njets_to_run = []
if (args.njets == "3") or (args.njets == "all"):
    njets_to_run += ["3Jets"]
if (args.njets == "4+") or (args.njets == "all"):
    njets_to_run += ["4PJets"]




def flatten_bkg_templates(fnames_to_run):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """
    if "3Jets" in njets_to_run:
        histo_dict_3j = processor.dict_accumulator({"Muon" : {}, "Electron" :{}})
    if "4PJets" in njets_to_run:
        histo_dict_4pj = processor.dict_accumulator({"Muon" : {}, "Electron" :{}})

    #set_trace()
    for bkg_file in fnames_to_run:
        hdict = load(bkg_file)
        jmult = "3Jets" if "3Jets" in os.path.basename(bkg_file) else "4PJets"
        for lep in hdict.keys():
            for tname, orig_template in hdict[lep].items():
                #set_trace()
                
                proc = tname.split("_")[0] if not "data_obs" in tname else "data_obs"
                sys = sorted(filter(None, tname.split(f"{proc}_")))[0]
                #if sys == "nosys": continue
                print(lep, jmult, sys, proc)

                    # perform flattening
                flattened_histo = hdict[lep][f"{proc}_nosys"].copy() if sys == "nosys" else Plotter.flatten(nosys=hdict[lep][f"{proc}_nosys"].copy(), systematic=orig_template.copy())
                
                    ## save template histos to coffea dict
                if jmult == "3Jets":
                    histo_dict_3j[lep][tname] = flattened_histo.copy()
                if jmult == "4PJets":
                    histo_dict_4pj[lep][tname] = flattened_histo.copy()

    #set_trace()
    if "3Jets" in njets_to_run:
        coffea_out_3j = os.path.join(input_dir, f"test_flattened_templates_lj_3Jets_bkg_{args.year}_{jobid}.coffea")
        save(histo_dict_3j, coffea_out_3j)
        print(f"{coffea_out_3j} written")
    if "4PJets" in njets_to_run:
        coffea_out_4pj = os.path.join(input_dir, f"test_flattened_templates_lj_4PJets_bkg_{args.year}_{jobid}.coffea")
        save(histo_dict_4pj, coffea_out_4pj)
        print(f"{coffea_out_4pj} written")



if __name__ == "__main__":
    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]

    analyzer = "htt_btag_sb_regions"
    base_template_name = f"test_raw_templates_lj_NJETS_bkg_{args.year}_{jobid}"

        # get matching pattern based on args.njets
    njets_regex = "*" if len(njets_to_run) > 1 else njets_to_run[0]

    input_dir = os.path.join(proj_dir, "results", "%s_%s" % (args.year, jobid), f"Templates_{analyzer}")
    if os.path.isdir(input_dir):
            # define variables to get histogram for background    
        bkg_fnmatch = "%s.coffea" % base_template_name.replace("NJETS", njets_regex)
        bkg_fnames = fnmatch.filter(os.listdir(input_dir), bkg_fnmatch)
        bkg_fnames = [os.path.join(input_dir, fname) for fname in bkg_fnames]
    else: print("No background file found.")

    try:
        print("Creating flattened background templates")
        flatten_bkg_templates(bkg_fnames)
        
    except:
        print("Could not write background templates to file")
