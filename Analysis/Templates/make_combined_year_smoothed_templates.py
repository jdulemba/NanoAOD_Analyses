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
import Utilities.final_analysis_binning as final_binning
    
base_jobid = os.environ["base_jobid"]
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--njets", default="all", nargs="?", choices=["3", "4+", "all"], help="Specify which jet multiplicity to use.")
parser.add_argument("--only_bkg", action="store_true", help="Make background templates only.")
parser.add_argument("--only_sig", action="store_true", help="Make signal templates only.")
args = parser.parse_args()

njets_to_run = []
if (args.njets == "3") or (args.njets == "all"):
    njets_to_run += ["3Jets"]
if (args.njets == "4+") or (args.njets == "all"):
    njets_to_run += ["4PJets"]



def smooth_bkg_templates(fnames_to_run):
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

                    # perform smoothing
                smoothed_histo = hdict[lep][f"{proc}_nosys"].copy() if sys == "nosys" else Plotter.smoothing_mttbins(nosys=hdict[lep][f"{proc}_nosys"], systematic=orig_template, mtt_centers=mtt_centers, nbinsx=len(linearize_binning[0])-1, nbinsy=len(linearize_binning[1])-1)
                
                    ## save template histos to coffea dict
                if jmult == "3Jets":
                    histo_dict_3j[lep][tname] = smoothed_histo.copy()
                if jmult == "4PJets":
                    histo_dict_4pj[lep][tname] = smoothed_histo.copy()

    #set_trace()
    if "3Jets" in njets_to_run:
        coffea_out_3j = os.path.join(input_dir, f"test_CombinedSmoothed_templates_lj_3Jets_bkg_{jobid}.coffea")
        save(histo_dict_3j, coffea_out_3j)
        print(f"{coffea_out_3j} written")
    if "4PJets" in njets_to_run:
        coffea_out_4pj = os.path.join(input_dir, f"test_CombinedSmoothed_templates_lj_4PJets_bkg_{jobid}.coffea")
        save(histo_dict_4pj, coffea_out_4pj)
        print(f"{coffea_out_4pj} written")


def smooth_sig_templates(fnames_to_run):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """
    if "3Jets" in njets_to_run:
        histo_dict_3j = processor.dict_accumulator({"Muon" : {}, "Electron" :{}})
    if "4PJets" in njets_to_run:
        histo_dict_4pj = processor.dict_accumulator({"Muon" : {}, "Electron" :{}})

    #set_trace()
    for sig_file in fnames_to_run:
        hdict = load(sig_file)
        jmult = "3Jets" if "3Jets" in os.path.basename(sig_file) else "4PJets"
        for lep in hdict.keys():
            for tname, orig_template in hdict[lep].items():
                #set_trace()

                if "Res" in tname:
                    signal = "_".join([tname.split("_Res")[0], "Res"])
                else:
                    signal = "_".join([tname.split("_neg")[0], "neg"]) if "neg" in tname else "_".join([tname.split("_pos")[0], "pos"])
                sys = sorted(filter(None, tname.split(f"{signal}_")))[0]
                print(lep, jmult, sys, signal)

                    # perform smoothing
                smoothed_histo = hdict[lep][f"{signal}_nosys"].copy() if sys == "nosys" else Plotter.smoothing_mttbins(nosys=hdict[lep][f"{signal}_nosys"], systematic=orig_template, mtt_centers=mtt_centers, nbinsx=len(linearize_binning[0])-1, nbinsy=len(linearize_binning[1])-1)
                #set_trace()
                    ## save template histos to coffea dict
                if jmult == "3Jets":
                    histo_dict_3j[lep][tname] = smoothed_histo.copy()
                if jmult == "4PJets":
                    histo_dict_4pj[lep][tname] = smoothed_histo.copy()

    #set_trace()
    if "3Jets" in njets_to_run:
        coffea_out_3j = os.path.join(input_dir, f"test_CombinedSmoothed_templates_lj_3Jets_sig_{jobid}.coffea")
        save(histo_dict_3j, coffea_out_3j)
        print(f"{coffea_out_3j} written")
    if "4PJets" in njets_to_run:
        coffea_out_4pj = os.path.join(input_dir, f"test_CombinedSmoothed_templates_lj_4PJets_sig_{jobid}.coffea")
        save(histo_dict_4pj, coffea_out_4pj)
        print(f"{coffea_out_4pj} written")


if __name__ == "__main__":

    instructions = """
    You must follow the following instructions in order to be able to run this script within this framework!
    Source a python3 environment with 'scl enable rh-python38 bash'
    Create a virtual environment by running 'python -m venv my_env' and 'source my_env/bin/activate' within 'Analysis' directory.
    If you are creating a new environment from scratch you must 'pip install coffea' and 'pip install statsmodels' and maybe 'pip install numpy==1.20'
    Add python paths by 'source venv.sh' and then run this script
    """
    print(instructions)
    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]

    analyzer = "htt_btag_sb_regions"
    base_bkg_template_name = f"test_CombinedRaw_templates_lj_NJETS_bkg_{jobid}"
    base_sig_template_name = f"test_CombinedRaw_templates_lj_NJETS_sig_{jobid}"

        # get matching pattern based on args.njets
    njets_regex = "*" if len(njets_to_run) > 1 else njets_to_run[0]

    input_dir = os.path.join(proj_dir, "results", jobid, f"Templates_{analyzer}")
    if os.path.isdir(input_dir):
            # define variables to get histogram for background    
        bkg_fnmatch = "%s.coffea" % base_bkg_template_name.replace("NJETS", njets_regex)
        bkg_fnames = fnmatch.filter(os.listdir(input_dir), bkg_fnmatch)
        bkg_fnames = [os.path.join(input_dir, fname) for fname in bkg_fnames]
            # define variables to get histogram for background    
        sig_fnmatch = "%s.coffea" % base_sig_template_name.replace("NJETS", njets_regex)
        sig_fnames = fnmatch.filter(os.listdir(input_dir), sig_fnmatch)
        sig_fnames = [os.path.join(input_dir, fname) for fname in sig_fnames]
    else: print("No files found.")

    linearize_binning = (
        final_binning.mtt_binning,
        final_binning.ctstar_abs_binning
    )
    #linearize_binning = (
    #    #np.array([300.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 625.0, 650.0, 675.0, 700.0, 730.0, 760.0, 800.0, 850.0, 900.0, 1000.0, 1200.0]), # orig
    #    np.array([360.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 625.0, 650.0, 675.0,
    #        700.0, 730.0, 760.0, 800.0, 850.0, 900.0, 950., 1000.0, 1050., 1100.0, 1150., 1200., 1300., 1500., 2000.]),
    #    #np.array([300.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 625.0, 650.0, 675.0,
    #    #    700.0, 730.0, 760.0, 800.0, 850.0, 900.0, 950., 1000.0, 1050., 1100.0, 1150., 1200., 1300., 1500., 2000.]),
    #    #np.array([300.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 620.0, 650., 700.0, 750.0, 800.0, 850.0, 900.0, 2000.0]),
    #    np.array([0.0, 0.4, 0.6, 0.75, 0.9, 1.0])
    #)    
    mtt_centers =  np.array([(linearize_binning[0][idx]+linearize_binning[0][idx+1])/2 for idx in range(len(linearize_binning[0])-1)])


    if not args.only_sig:
        print("Creating smoothed background templates")
        smooth_bkg_templates(bkg_fnames)

    if not args.only_bkg:
        print("Creating smoothed signal templates")
        smooth_sig_templates(sig_fnames)

    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
