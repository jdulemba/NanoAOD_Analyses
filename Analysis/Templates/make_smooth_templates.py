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
import Utilities.systematics as systematics
    
base_jobid = os.environ["base_jobid"]
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("templates_to_run", type=str, help="Choose which type of templates to run, multiple options can be input as ':' separated strings.")
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
                #print(lep, jmult, sys, proc)

                if sys == "nosys":
                    smoothed_histo = hdict[jmult][lep][f"{proc}_nosys"].copy()
                else:
                    #if "sub" in sys: set_trace()
                    #set_trace()
                    sysname = "_".join(systematics.sys_to_name[args.year][sys].split("_")[:-1]) # convert sys to name used in analysis
                    if sysname in bkg_smooth_pars.treatment[args.year]["Indiv"].keys():
                        p_end = bkg_smooth_pars.treatment[args.year]["Indiv"][sysname][f"{lep[0]}{jmult[0]}"][proc]
                    else:
                        #set_trace()
                        if sysname in bkg_smooth_pars.treatment["Combined_Era_Lep"][jmult].keys():
                            p_end = bkg_smooth_pars.treatment["Combined_Era_Lep"][jmult][sysname] 
                        elif sysname in bkg_smooth_pars.treatment[args.year]["Combined_Lep"][jmult].keys():
                            p_end = bkg_smooth_pars.treatment[args.year]["Combined_Lep"][jmult][sysname][proc] if isinstance(bkg_smooth_pars.treatment[args.year]["Combined_Lep"][jmult][sysname], dict)\
                                else bkg_smooth_pars.treatment[args.year]["Combined_Lep"][jmult][sysname]
                        else:
                            continue

                    print(lep, jmult, sys, proc)

                        # perform smoothing
                    smoothed_histo = Plotter.new_smoothing(nosys=hdict[jmult][lep][f"{proc}_nosys"].copy(), systematic=orig_template.copy(), mtt_bins=linearize_binning[0],
                        nbinsx=len(linearize_binning[0])-1, nbinsy=len(linearize_binning[1])-1, p_end=p_end)

                    ## save template histos to coffea dict
                histo_dict[jmult][lep][tname] = smoothed_histo.copy()

    coffea_out = os.path.join(input_dir, f"smoothed_templates_lj_bkg_{args.year}_{jobid}.coffea")
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

    coffea_out = os.path.join(input_dir, f"smoothed_templates_lj_MEsig_{args.year}_{jobid}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")


def smooth_pdf_templates(fname):
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

                smoothed_histo = hdict[jmult][lep][f"{proc}_nosys"].copy() if sys == "nosys" else \
                    Plotter.new_smoothing(nosys=hdict[jmult][lep][f"{proc}_nosys"].copy(), systematic=orig_template.copy(), mtt_bins=linearize_binning[0], nbinsx=len(linearize_binning[0])-1, nbinsy=len(linearize_binning[1])-1, p_end=0.000001)
                    #Plotter.new_smoothing(nosys=hdict[jmult][lep][f"{proc}_nosys"].copy(), systematic=orig_template.copy(), mtt_bins=linearize_binning[0], nbinsx=len(linearize_binning[0])-1, nbinsy=len(linearize_binning[1])-1, p_end=0.001)

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

    allowed_template_options = ["bkg", "MEreweight_sig", "PDF"]
    templates_to_run = (args.templates_to_run).split(":")
    templates_to_run = [template for template in templates_to_run if template in allowed_template_options]

    jobid = os.environ["jobid"]
    eos_dir = os.environ["eos_dir"]

    if jobid == "Summer20UL_DeepJet":
        import bkg_smooth_parameters_Summer20UL_DeepJet as bkg_smooth_pars
    else:
        import bkg_smooth_parameters as bkg_smooth_pars

    njets_to_run = ["3Jets", "4PJets"]

    analyzer = "htt_btag_sb_regions"

    input_dir = os.path.join(eos_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}")

    linearize_binning = (
        final_binning.mtt_binning,
        final_binning.ctstar_abs_binning
    )
        # replace last mtt bin with 1700 instead of 2000
    #set_trace()
    linearize_binning[0][-2] = linearize_binning[0][-3] + 100.
    linearize_binning[0][-1] = linearize_binning[0][-2] + 100.
    #linearize_binning[0][-1] = 1700.

    mtt_centers =  np.array([(linearize_binning[0][idx]+linearize_binning[0][idx+1])/2 for idx in range(len(linearize_binning[0])-1)])

    if "bkg" in templates_to_run:
        base_bkg_template_name = f"raw_templates_lj_bkg_{args.year}_{jobid}.coffea"
        bkg_fname = os.path.join(input_dir, base_bkg_template_name)
        if not os.path.isfile(bkg_fname): raise ValueError("No background file found.")
        print("Creating smoothed background templates")
        smooth_bkg_templates(bkg_fname)

    if "MEreweight_sig" in templates_to_run:
        base_sig_template_name = f"raw_templates_lj_MEsig_{args.year}_{jobid}.coffea"
        sig_fname = os.path.join(input_dir, base_sig_template_name)
        if not os.path.isfile(sig_fname): raise ValueError("No signal file found.")
        print("Creating ME reweighting signal templates")
        smooth_MEsig_templates(sig_fname)

    if "PDF" in templates_to_run:
        analyzer = "htt_pdfUncs"

        input_dir = os.path.join(eos_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}")
        bkg_fname = os.path.join(input_dir, f"raw_pdf_templates_lj_bkg_{args.year}_{jobid}.coffea")
        if not os.path.isfile(bkg_fname): raise ValueError("No background file found.")
        print("Creating smoothed background templates")
        smooth_pdf_templates(bkg_fname)

    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
