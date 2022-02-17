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
    
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
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

    coffea_out = os.path.join(input_dir, f"smoothed_pdf_templates_lj_bkg_{args.year}_{jobid}.coffea")
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
    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]

    njets_to_run = ["3Jets", "4PJets"]

    linearize_binning = (
        final_binning.mtt_binning,
        final_binning.ctstar_abs_binning
    )
    mtt_centers =  np.array([(linearize_binning[0][idx]+linearize_binning[0][idx+1])/2 for idx in range(len(linearize_binning[0])-1)])

    analyzer = "htt_pdfUncs"
    base_bkg_template_name = f"raw_pdf_templates_lj_bkg_{args.year}_{jobid}.coffea"
    input_dir = os.path.join(proj_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}")
    bkg_fname = os.path.join(input_dir, base_bkg_template_name)
    if not os.path.isfile(bkg_fname): raise ValueError("No background file found.")
    print("Creating smoothed background templates")
    smooth_bkg_templates(bkg_fname)

    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
