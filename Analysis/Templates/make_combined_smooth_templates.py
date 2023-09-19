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
    
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("templates_to_run", type=str, help="Choose which type of templates to run, multiple options can be input as ':' separated strings.")
parser.add_argument("channels_to_combine", type=str, help="Choose how to combine systematic templates (by lepton or eras), multiple options can be input as ':' separated strings.")
parser.add_argument("--nomSMTTxsec", action="store_true", help="Apply nominal SM cross sections to top mass and LHE scale weights")
args = parser.parse_args()


def smooth_year_and_lepton_templates(fname, process):
    histo_dict = processor.dict_accumulator({njets : {} for njets in njets_to_run})
    hdict = load(fname)
    for jmult in hdict.keys():
        for tname, orig_template in hdict[jmult].items():
            if (process == "bkg") or (process == "PDF"):
                proc = tname.split("_")[0] if not "data_obs" in tname else "data_obs"
            else:
                if "Res" in tname:
                    proc = "_".join([tname.split("_Res")[0], "Res"])
                else:
                    proc = "_".join([tname.split("_neg")[0], "neg"]) if "neg" in tname else "_".join([tname.split("_pos")[0], "pos"])

            sys = sorted(filter(None, tname.split(f"{proc}_")))[0]
            if sys == "nosys":
                smoothed_histo = hdict[jmult][f"{proc}_nosys"].copy()
            else:
                #print(sys)
                #set_trace()
                if process == "PDF": p_end = 0.000001
                else:
                    if (sys == "CR1") or (sys == "CR2") or (sys == "erdON"):
                        #set_trace()
                        sysname = "_".join(systematics.combine_template_sys_to_name["2017"][sys+"Up"].split("_")).replace("Down", "").replace("Up", "")
                    else:
                        if sys not in systematics.combine_template_sys_to_name["2017"].keys(): continue
                        sysname = "_".join(systematics.combine_template_sys_to_name["2017"][sys].split("_")).replace("Down", "").replace("Up", "") # convert to name 
                    if sysname in bkg_smooth_pars.treatment["Combined_Era_Lep"][jmult].keys():
                        p_end = bkg_smooth_pars.treatment["Combined_Era_Lep"][jmult][sysname][proc] if isinstance(bkg_smooth_pars.treatment["Combined_Era_Lep"][jmult][sysname], dict) \
                            else bkg_smooth_pars.treatment["Combined_Era_Lep"][jmult][sysname]
                    else:
                        continue

                print(jmult, sys, proc)
                    # perform smoothing
                smoothed_histo = Plotter.new_smoothing(nosys=hdict[jmult][f"{proc}_nosys"].copy(), systematic=orig_template.copy(), mtt_bins=linearize_binning[0],
                    nbinsx=len(linearize_binning[0])-1, nbinsy=len(linearize_binning[1])-1, ratio=True, p_end=p_end, cutoff=5.0)# if process == "PDF" else 3.0)
            
                ## save template histos to coffea dict
            histo_dict[jmult][tname] = smoothed_histo.copy()

    #set_trace()
    coffea_out = os.path.join(input_dir, f"smoothed_combined_year_and_lepton_templates_lj_{process}_nomSMTTxsec_{jobid}.coffea" if args.nomSMTTxsec else f"smoothed_combined_year_and_lepton_templates_lj_{process}_{jobid}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")


def smooth_lepton_templates(fname, process):
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
            if sys == "nosys":
                smoothed_histo = hdict[jmult][f"{proc}_nosys"].copy()
            else:
                #set_trace()
                if sys not in systematics.combine_template_sys_to_name[year].keys():
                    #set_trace()
                    continue
                sysname = "_".join(systematics.combine_template_sys_to_name[year][sys].split("_")).replace("Down", "").replace("Up", "") # convert to name 
                if sysname in bkg_smooth_pars.treatment[year]["Combined_Lep"][jmult].keys():
                    p_end = bkg_smooth_pars.treatment[year]["Combined_Lep"][jmult][sysname][proc] if isinstance(bkg_smooth_pars.treatment[year]["Combined_Lep"][jmult][sysname], dict) \
                        else bkg_smooth_pars.treatment[year]["Combined_Lep"][jmult][sysname]
                else:
                    continue

                print(year, jmult, sys, proc)
                    # perform smoothing
                smoothed_histo = Plotter.new_smoothing(nosys=hdict[jmult][f"{proc}_nosys"].copy(), systematic=orig_template.copy(), mtt_bins=linearize_binning[0],
                    nbinsx=len(linearize_binning[0])-1, nbinsy=len(linearize_binning[1])-1, ratio=True, p_end=p_end)
            
                ## save template histos to coffea dict
            histo_dict[jmult][tname] = smoothed_histo.copy()

    #set_trace()
    coffea_out = os.path.join(input_dir, f"smoothed_combined_lep_templates_lj_{process}_{year}_{jobid}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")



if __name__ == "__main__":
    allowed_template_options = ["bkg", "PDF"]
    templates_to_run = (args.templates_to_run).split(":")
    templates_to_run = [template for template in templates_to_run if template in allowed_template_options]

    allowed_combination_options = ["lepton", "era_lepton"]
    combinations_to_run = (args.channels_to_combine).split(":")
    combinations_to_run = [combination for combination in combinations_to_run if combination in allowed_combination_options]

    jobid = os.environ["jobid"]
    eos_dir = os.environ["eos_dir"]

    if jobid == "Summer20UL_DeepJet":
        import bkg_smooth_parameters_Summer20UL_DeepJet as bkg_smooth_pars
    else:
        import bkg_smooth_parameters as bkg_smooth_pars

    njets_to_run = ["3Jets", "4PJets"]
    analyzer = "htt_btag_sb_regions"

    linearize_binning = (
        final_binning.mtt_binning,
        final_binning.ctstar_abs_binning
    )
        # replace last mtt bin with 1700 instead of 2000
    #linearize_binning[0][-2] = linearize_binning[0][-3] + 50.
    #linearize_binning[0][-1] = linearize_binning[0][-2] + 50.
    linearize_binning[0][-2] = linearize_binning[0][-3] + 100.
    linearize_binning[0][-1] = linearize_binning[0][-2] + 100.

    mtt_centers =  np.array([(linearize_binning[0][idx]+linearize_binning[0][idx+1])/2 for idx in range(len(linearize_binning[0])-1)])


    for combination in combinations_to_run:
        if combination == "lepton":
            for year in ["2016APV", "2016", "2017", "2018"]:
                input_dir = os.path.join(eos_dir, "results", f"{year}_{jobid}", f"Templates_{analyzer}")
                if not os.path.isdir(input_dir): raise ValueError("No background file found.")

                if "bkg" in templates_to_run:
                    base_bkg_template_name = f"raw_combined_lep_templates_lj_bkg_{year}_{jobid}.coffea"
                    bkg_fname = os.path.join(input_dir, base_bkg_template_name)
                    if not os.path.isfile(bkg_fname): raise ValueError(f"{bkg_fname} not found.")
                    print(f"Smoothening e+mu channels in {year} for background templates")
                    smooth_lepton_templates(bkg_fname, "bkg")

        if combination == "era_lepton":
            if "PDF" in templates_to_run:
                analyzer = "htt_pdfUncs"
                input_dir = os.path.join(eos_dir, "results", jobid, f"Templates_{analyzer}")
                if not os.path.isdir(input_dir): raise ValueError("No directory found.")

                pdf_fname = os.path.join(input_dir, f"raw_combined_year_and_lepton_templates_lj_PDF_{jobid}.coffea")
                if not os.path.isfile(pdf_fname): raise ValueError(f"{pdf_fname} not found.")
                print("Smoothening combined channels across years and leptons for PDF templates")
                smooth_year_and_lepton_templates(pdf_fname, "PDF")

            else:
                input_dir = os.path.join(eos_dir, "results", jobid, f"Templates_{analyzer}")
                if not os.path.isdir(input_dir): raise ValueError("No directory found.")

                if "bkg" in templates_to_run:
                    bkg_fname = os.path.join(input_dir, f"raw_combined_year_and_lepton_templates_lj_bkg_nomSMTTxsec_{jobid}.coffea" if args.nomSMTTxsec \
                        else f"raw_combined_year_and_lepton_templates_lj_bkg_{jobid}.coffea")
                    if not os.path.isfile(bkg_fname): raise ValueError(f"{bkg_fname} not found.")
                    print("Smoothening combined channels across years and leptons for background templates")
                    smooth_year_and_lepton_templates(bkg_fname, "bkg")


    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
