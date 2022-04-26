#! /bin/env python

import time
tic = time.time()

from coffea.util import load, save
from pdb import set_trace
import os
import numpy as np
import Utilities.Plotter as Plotter
import coffea.processor as processor    
import Utilities.final_analysis_binning as final_binning
    
base_jobid = os.environ["base_jobid"]
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("templates_to_run", type=str, help="Choose which type of templates to run, multiple options can be input as ':' separated strings.")
parser.add_argument("--kfactors", action="store_true", help="Apply signal k-factors to signal")
parser.add_argument("--scale_mtop3gev", action="store_true", help="Scale 3GeV mtop variations by 1/6")
args = parser.parse_args()


def smooth_bkg_templates(fname):
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

                    # perform smoothing
                smoothed_histo = hdict[jmult][lep][f"{proc}_nosys"].copy() if sys == "nosys" \
                    else Plotter.smoothing_mttbins(nosys=hdict[jmult][lep][f"{proc}_nosys"], systematic=orig_template, mtt_centers=mtt_centers, nbinsx=len(linearize_binning[0])-1, nbinsy=len(linearize_binning[1])-1)
                
                    ## save template histos to coffea dict
                histo_dict[jmult][lep][tname] = smoothed_histo.copy()

    coffea_out = os.path.join(input_dir, f"smoothed_templates_lj_bkg_mtopscaled_{args.year}_{jobid}.coffea" if args.scale_mtop3gev else f"smoothed_templates_lj_bkg_{args.year}_{jobid}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")


def smooth_sig_templates(fname):
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

                    # perform smoothing
                smoothed_histo = hdict[jmult][lep][f"{signal}_nosys"].copy() if sys == "nosys" \
                    else Plotter.smoothing_mttbins(nosys=hdict[jmult][lep][f"{signal}_nosys"], systematic=orig_template, mtt_centers=mtt_centers, nbinsx=len(linearize_binning[0])-1, nbinsy=len(linearize_binning[1])-1)
                    ## save template histos to coffea dict
                histo_dict[jmult][lep][tname] = smoothed_histo.copy()

    coffea_out = os.path.join(input_dir, f"smoothed_templates_lj_sig_kfactors_{args.year}_{jobid}.coffea" if args.kfactors else f"smoothed_templates_lj_sig_{args.year}_{jobid}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")


def smooth_MEsig_templates(fname):
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

                    # perform smoothing
                smoothed_histo = hdict[jmult][lep][f"{signal}_nosys"].copy() if sys == "nosys" \
                    else Plotter.smoothing_mttbins(nosys=hdict[jmult][lep][f"{signal}_nosys"], systematic=orig_template, mtt_centers=mtt_centers, nbinsx=len(linearize_binning[0])-1, nbinsy=len(linearize_binning[1])-1)
                    ## save template histos to coffea dict
                histo_dict[jmult][lep][tname] = smoothed_histo.copy()

    coffea_out = os.path.join(input_dir, f"smoothed_templates_lj_MEsig_kfactors_{args.year}_{jobid}.coffea" if args.kfactors else f"smoothed_templates_lj_MEsig_{args.year}_{jobid}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")


if __name__ == "__main__":

    instructions = """
    You must follow the following instructions in order to be able to run this script within this framework!
    Source a python3 environment with 'scl enable rh-python38 bash'
    Create a virtual environment by running 'python -m venv my_env' and 'source my_env/bin/activate' within 'Analysis' directory.
    If you are creating a new environment from scratch you must 'pip install coffea' and 'pip install statsmodels' and maybe 'pip install numpy==1.20'
    Add python paths by 'source venv.sh' and then run this script
    """
    print(instructions)

    allowed_template_options = ["bkg", "sig", "MEreweight_sig"]
    templates_to_run = (args.templates_to_run).split(":")
    templates_to_run = [template for template in templates_to_run if template in allowed_template_options]

    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]

    njets_to_run = ["3Jets", "4PJets"]

    analyzer = "htt_btag_sb_regions"

        # get matching pattern based on args.njets
    njets_regex = "*" if len(njets_to_run) > 1 else njets_to_run[0]

    input_dir = os.path.join(proj_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}")

    linearize_binning = (
        final_binning.mtt_binning,
        final_binning.ctstar_abs_binning
    )
    mtt_centers =  np.array([(linearize_binning[0][idx]+linearize_binning[0][idx+1])/2 for idx in range(len(linearize_binning[0])-1)])

    if "bkg" in templates_to_run:
        base_bkg_template_name = f"raw_templates_lj_bkg_mtopscaled_{args.year}_{jobid}.coffea" if args.scale_mtop3gev else f"raw_templates_lj_bkg_{args.year}_{jobid}.coffea"
        bkg_fname = os.path.join(input_dir, base_bkg_template_name)
        if not os.path.isfile(bkg_fname): raise ValueError("No background file found.")
        print("Creating smoothed background templates")
        smooth_bkg_templates(bkg_fname)

    if "sig" in templates_to_run:
        base_sig_template_name = f"raw_templates_lj_sig_kfactors_{args.year}_{jobid}.coffea" if args.kfactors else f"raw_templates_lj_sig_{args.year}_{jobid}.coffea"
        sig_fname = os.path.join(input_dir, base_sig_template_name)
        if not os.path.isfile(sig_fname): raise ValueError("No signal file found.")
        print("Creating smoothed signal templates")
        smooth_sig_templates(sig_fname)

    if "MEreweight_sig" in templates_to_run:
        base_sig_template_name = f"raw_templates_lj_MEsig_kfactors_{args.year}_{jobid}.coffea" if args.kfactors else f"raw_templates_lj_MEsig_{args.year}_{jobid}.coffea"
        sig_fname = os.path.join(input_dir, base_sig_template_name)
        if not os.path.isfile(sig_fname): raise ValueError("No signal file found.")
        print("Creating ME reweighting signal templates")
        smooth_MEsig_templates(sig_fname)

    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
