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
                #set_trace()
                smoothed_histo = unrolled_to_1d_hist(orig_template.copy())
            else:
                if process == "PDF": p_end = 0.000001
                else:
                    if (sys == "CR1") or (sys == "CR2") or (sys == "erdON"):
                        #set_trace()
                        sysname = "_".join(systematics.combine_template_sys_to_name["2017"][sys+"Up"].split("_")).replace("Down", "").replace("Up", "")
                    else:
                        if sys not in systematics.combine_template_sys_to_name["2017"].keys(): continue
                        sysname = "_".join(systematics.combine_template_sys_to_name["2017"][sys].split("_")).replace("Down", "").replace("Up", "") # convert to name
                    if sysname in bkg_smooth_pars.treatment["Combined_Era_Lep"][jmult].keys():
                        p_end = bkg_smooth_pars.treatment["Combined_Era_Lep"][jmult][sysname]
                    else:
                        continue

                print(jmult, sys, proc)
                    # perform smoothing
                smoothed_histo = Plotter.new_smoothing(nosys=unrolled_to_1d_hist(hdict[jmult][f"{proc}_nosys"].copy()), systematic=unrolled_to_1d_hist(orig_template.copy()),
                     mtt_bins=linearize_binning[0], nbinsx=len(linearize_binning[0])-1, nbinsy=len(linearize_binning[1])-1, ratio=True, p_end=p_end, cutoff=5.0)

                ## save template histos to coffea dict
            histo_dict[jmult][tname] = smoothed_histo.copy()

    #set_trace()
    coffea_out = os.path.join(input_dir, f"smoothed_1D_combined_year_and_lepton_templates_lj_{process}_{jobid}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")


def unrolled_to_1d_hist(histo):
    """
    This function take an unrolled coffea histogram object, makes it 1D in terms of mtt, and then repeats it so the dimensionality is the same as the input.
    """
    #set_trace()
    unrolled_sumw, unrolled_sumw2 = histo.values(sumw2=True)[()]
        # convert 1D unrolled array back to 2D so the ctstar bins can be integrated
    sumw_2d, sumw2_2d = np.reshape(unrolled_sumw, (linearize_binning[1].size-1, linearize_binning[0].size-1)), np.reshape(unrolled_sumw2, (linearize_binning[1].size-1, linearize_binning[0].size-1))
        # integrate ctstar bins
    sumw_1d, sumw2_1d = np.sum(sumw_2d, axis=0), np.sum(sumw2_2d, axis=0)
        # repeat mtt bins Nctstar bins times
    sumw_1d_tiled, sumw2_1d_tiled = np.tile(sumw_1d, linearize_binning[1].size-1), np.tile(sumw2_1d, linearize_binning[1].size-1)
        # make hist out of these values
    tiled_1d_hist = histo.copy()#hdict[jmult][f"{proc}_nosys"].copy()
    tiled_1d_hist.values()[()][:] = sumw_1d_tiled
    tiled_1d_hist.values(sumw2=True)[()][1][:] = sumw2_1d_tiled

    return tiled_1d_hist


if __name__ == "__main__":
    allowed_template_options = ["bkg"]
    templates_to_run = (args.templates_to_run).split(":")
    templates_to_run = [template for template in templates_to_run if template in allowed_template_options]

    allowed_combination_options = ["era_lepton"]#, "lepton"]
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
    linearize_binning[0][-2] = linearize_binning[0][-3] + 100.
    linearize_binning[0][-1] = linearize_binning[0][-2] + 100.

    mtt_centers =  np.array([(linearize_binning[0][idx]+linearize_binning[0][idx+1])/2 for idx in range(len(linearize_binning[0])-1)])

    for combination in combinations_to_run:
        if combination == "era_lepton":
            input_dir = os.path.join(eos_dir, "results", jobid, f"Templates_{analyzer}")
            if not os.path.isdir(input_dir): raise ValueError("No directory found.")

            if "bkg" in templates_to_run:
                bkg_fname = os.path.join(input_dir, f"raw_combined_year_and_lepton_templates_lj_bkg_{jobid}.coffea")
                if not os.path.isfile(bkg_fname): raise ValueError(f"{bkg_fname} not found.")
                print("Smoothening combined channels across years and leptons for background templates")
                smooth_year_and_lepton_templates(bkg_fname, "bkg")

    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
