#! /bin/env python

import time
tic = time.time()

from coffea.util import load, save
from pdb import set_trace
import os
import numpy as np
import fnmatch
import coffea.processor as processor
import Utilities.systematics as systematics
import Utilities.prettyjson as prettyjson

import argparse
class ParseKwargs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, dict())
        for value in values:
            key, value = value.split("=")
            getattr(namespace, self.dest)[key] = value

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("templates_to_run", type=str, help="Choose which type of templates to run, multiple options can be input as ':' separated strings.")
parser.add_argument("--MEopts", nargs="*", action=ParseKwargs, help="Options to specify which mass and width points to produce for ME reweighted signal.")
parser.add_argument("--nomSMTTxsec", action="store_true", help="Apply nominal SM cross sections to top mass and LHE scale weights")
args = parser.parse_args()

version = "V33"
#version = "V31"
#version = "V30"
#version = "V29"
#version = "V27"

def final_bkg_templates(hdicts):
    histo_dict = processor.dict_accumulator({njets : {"Muon" : {}, "Electron" :{}} for njets in njets_to_run})
    raw_hdict = hdicts["Indiv_Raw"]
    for jmult in sorted(raw_hdict.keys()):
        for lep in sorted(raw_hdict[jmult].keys()):
            #set_trace()
            for tname in sorted(raw_hdict[jmult][lep].keys()):
                proc = tname.split("_")[0] if not "data_obs" in tname else "data_obs"
                sys = sorted(filter(None, tname.split(f"{proc}_")))[0]

                if (sys == "CR1") or (sys == "CR2") or (sys == "erdON"):
                    #continue
                    treatment = bkg_syspar.treatment[args.year][sys][proc][f"{lep.lower()[0]}{jmult[0]}"]

                    print(f"{args.year}, {jmult}, {lep}, {sys}, {proc}:\t{treatment}")
                    if treatment not in hdicts.keys():
                        raise ValueError(f"{treatment} is not a supported treatment!")

                    #set_trace()
                    for variation in ["Up", "Down"]:
                        if "Combined_Era_Lep" in treatment:
                            if "Raw" in treatment:
                                set_trace()
                                ratio_hvals = (hdicts[treatment][jmult][tname].copy()).values()[()]/(hdicts[treatment][jmult][f"{proc}_nosys"].copy()).values()[()]
                            else:
                                ratio_hvals = (hdicts[treatment][jmult][f"{tname}{variation}"].copy()).values()[()]
                            hist_to_use = hdicts["Indiv_Raw"][jmult][lep][f"{proc}_nosys"].copy()
                            sys_vals = hist_to_use.values()[()] * ratio_hvals
                            hist_to_use.values()[()][:] = sys_vals

                        if "Combined_Lep" in treatment:
                            set_trace()
                            if "Raw" in treatment:
                                ratio_hvals = (hdicts[treatment][jmult][tname].copy()).values()[()]/(hdicts[treatment][jmult][f"{proc}_nosys"].copy()).values()[()]
                            else:
                                ratio_hvals = (hdicts[treatment][jmult][tname].copy()).values()[()]
                            hist_to_use = hdicts["Indiv_Raw"][jmult][lep][f"{proc}_nosys"].copy()
                            sys_vals = hist_to_use.values()[()] * ratio_hvals
                            hist_to_use.values()[()][:] = sys_vals

                        if "Indiv" in treatment:
                            set_trace()
                            hist_to_use = hdicts[treatment][jmult][lep][tname].copy()

                            ## save template histos to coffea dict
                        histo_dict[jmult][lep][f"{tname}{variation}"] = [hist_to_use.copy(), treatment]

                else:
                    if sys == "nosys":
                        treatment = "Indiv_Raw"
                    # choose which hist to use based on what_to_do dict in systs_parameters.py'
                    else:
                        sysname = "_".join(systematics.sys_to_name[args.year][sys].split("_")[:-1]) # convert sys to name used in analysis
                        if sysname in bkg_syspar.treatment[args.year].keys():
                            treatment = bkg_syspar.treatment[args.year][sysname][proc][f"{lep.lower()[0]}{jmult[0]}"]
                        else: continue

                    print(f"{args.year}, {jmult}, {lep}, {sys}, {proc}:\t{treatment}")
                    if treatment not in hdicts.keys():
                        raise ValueError(f"{treatment} is not a supported treatment!")

                    if "Combined_Era_Lep" in treatment:
                        if "Raw" in treatment:
                            ratio_hvals = (hdicts[treatment][jmult][tname].copy()).values()[()]/(hdicts[treatment][jmult][f"{proc}_nosys"].copy()).values()[()]
                        else:
                            ratio_hvals = (hdicts[treatment][jmult][tname].copy()).values()[()]

                            ## add 0.5 to all ratio values for ME scale variations for single top s channel processes by hand! for some reason the weights are centered around 0.5 instead of 1.
                        if ((sys.startswith("ST_RENORM")) or (sys.startswith("ST_FACTOR"))) and (proc == "TB"):
                            ratio_hvals += 0.5

                        hist_to_use = hdicts["Indiv_Raw"][jmult][lep][f"{proc}_nosys"].copy()
                        sys_vals = hist_to_use.values()[()] * ratio_hvals
                        hist_to_use.values()[()][:] = sys_vals

                    if "Combined_Lep" in treatment:
                        if "Raw" in treatment:
                            ratio_hvals = (hdicts[treatment][jmult][tname].copy()).values()[()]/(hdicts[treatment][jmult][f"{proc}_nosys"].copy()).values()[()]
                        else:
                            ratio_hvals = (hdicts[treatment][jmult][tname].copy()).values()[()]
                        hist_to_use = hdicts["Indiv_Raw"][jmult][lep][f"{proc}_nosys"].copy()
                        sys_vals = hist_to_use.values()[()] * ratio_hvals
                        hist_to_use.values()[()][:] = sys_vals

                    if "Indiv" in treatment:
                        hist_to_use = hdicts[treatment][jmult][lep][tname].copy()

                        ## save template histos to coffea dict
                    histo_dict[jmult][lep][tname] = [hist_to_use.copy(), treatment]

    coffea_out = os.path.join(outdir, f"final_templates_lj_bkg_nomSMTTxsec_{args.year}_{jobid}_{version}.coffea" if args.nomSMTTxsec else f"final_templates_lj_bkg_{args.year}_{jobid}_{version}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")


def final_pdf_templates(hdicts):
    """
    Function that chooses final templates.
    """
    histo_dict = processor.dict_accumulator({njets : {"Muon" : {}, "Electron" :{}} for njets in njets_to_run})

    #raw_hdict = hdicts["Indiv_Raw"]
    smooth_hdict = hdicts["Combined_Era_Lep_Smooth"]
    for jmult in smooth_hdict.keys():
        for tname in smooth_hdict[jmult].keys():
            proc = tname.split("_")[0] if not "data_obs" in tname else "data_obs"
            sys = sorted(filter(None, tname.split(f"{proc}_")))[0]

            #set_trace()
            if sys == "nosys": treatment = "Indiv_Nom"
            elif "alphaS" in sys: treatment = "Combined_Era_Lep_Raw"
            else: treatment = "Combined_Era_Lep_Smooth"

            for lep in ["Muon", "Electron"]:
                print(lep, jmult, sys, proc, treatment)

                if "Combined_Era_Lep" in treatment:
                    #set_trace()
                    if "Raw" in treatment:
                        ratio_hvals = hdicts[treatment][jmult][tname].values()[()]/hdicts[treatment][jmult][f"{proc}_nosys"].values()[()]
                    else:
                        ratio_hvals = hdicts[treatment][jmult][tname].values()[()]
                    hist_to_use = hdicts["Indiv_Nom"][jmult][lep][f"{proc}_nosys"].copy()
                    sys_vals = hist_to_use.values()[()] * ratio_hvals
                    hist_to_use.values()[()][:] = sys_vals

                if "Indiv" in treatment:
                    #set_trace()
                    if "Nom" in treatment:
                        ratio_hvals = np.ones(hdicts[treatment][jmult][lep][f"{proc}_nosys"].values()[()].size)
                    if "Orig" in treatment: # find original ratios and then rescale to actual ttbar
                        ratio_hvals = hdicts[treatment][jmult][lep][tname].values()[()]/hdicts[treatment][jmult][lep][f"{proc}_nosys"].values()[()]

                    hist_to_use = hdicts["Indiv_Nom"][jmult][lep][f"{proc}_nosys"].copy()
                    sys_vals = hist_to_use.values()[()] * ratio_hvals
                    hist_to_use.values()[()][:] = sys_vals

                    ## save template histos to coffea dict
                histo_dict[jmult][lep][tname] = [hist_to_use.copy(), treatment]

    coffea_out = os.path.join(outdir, f"final_pdf_templates_lj_bkg_{args.year}_{jobid}_{version}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")



def make_input_output_dir(analyzer):
    input_dir = os.path.join(eos_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}")
    output_dir = os.path.join(eos_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}", "FINAL")
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    return input_dir, output_dir


if __name__ == "__main__":
    allowed_template_options = ["bkg", "PDF"]
    templates_to_run = [template for template in (args.templates_to_run).split(":") if template in allowed_template_options]
    templates_to_not_run = [template for template in (args.templates_to_run).split(":") if template not in allowed_template_options]
    if templates_to_not_run:
        print(f"{templates_to_not_run} are not valid options for making templates, will be skipped")

    jobid = os.environ["jobid"]
    eos_dir = os.environ["eos_dir"]

    if jobid == "Summer20UL_DeepJet":
        import bkg_systs_parameters_Summer20UL_DeepJet as bkg_syspar
    else:
        import bkg_systs_parameters as bkg_syspar

    njets_to_run = sorted(["3Jets", "4PJets"])
    leps_to_run = sorted(["Muon", "Electron"])

    if "bkg" in templates_to_run:
        analyzer = "htt_btag_sb_regions"
        input_dir, outdir = make_input_output_dir(analyzer)
        if not os.path.join(input_dir): raise ValueError("No input directory found.")

        comb_era_lep_indir = os.path.join(eos_dir, "results", jobid, f"Templates_{analyzer}")
        if not os.path.join(comb_era_lep_indir): raise ValueError("No directory found for combined era and lepton channels.")

        bkg_dict = {
            "Combined_Era_Lep_Raw" : load(os.path.join(comb_era_lep_indir, f"raw_combined_year_and_lepton_templates_lj_bkg_nomSMTTxsec_{jobid}.coffea" if args.nomSMTTxsec else f"raw_combined_year_and_lepton_templates_lj_bkg_{jobid}.coffea")),
            "Combined_Era_Lep_Smooth" : load(os.path.join(comb_era_lep_indir, f"smoothed_combined_year_and_lepton_templates_lj_bkg_nomSMTTxsec_{jobid}.coffea" if args.nomSMTTxsec else f"smoothed_combined_year_and_lepton_templates_lj_bkg_{jobid}.coffea")),
            "Combined_Era_Lep_Symm_Smooth" : load(os.path.join(comb_era_lep_indir, f"symmetrized_smoothed_combined_year_and_lepton_templates_lj_bkg_nomSMTTxsec_{jobid}.coffea" if args.nomSMTTxsec else f"symmetrized_smoothed_combined_year_and_lepton_templates_lj_bkg_{jobid}.coffea")),
            "Combined_Era_Lep_Symm_Flat" : load(os.path.join(comb_era_lep_indir, f"symmetrized_flattened_combined_year_and_lepton_templates_lj_bkg_nomSMTTxsec_{jobid}.coffea" if args.nomSMTTxsec else f"symmetrized_flattened_combined_year_and_lepton_templates_lj_bkg_{jobid}.coffea")),
            "Combined_Lep_Raw" : load(os.path.join(input_dir, f"raw_combined_lep_templates_lj_bkg_{args.year}_{jobid}.coffea")),
            "Combined_Lep_Smooth" : load(os.path.join(input_dir, f"smoothed_combined_lep_templates_lj_bkg_{args.year}_{jobid}.coffea")),
            "Combined_Lep_Symm_Smooth" : load(os.path.join(input_dir, f"symmetrized_smoothed_combined_lep_templates_lj_bkg_{args.year}_{jobid}.coffea")),
            "Combined_Lep_Symm_Flat" : load(os.path.join(input_dir, f"symmetrized_flattened_combined_lep_templates_lj_bkg_{args.year}_{jobid}.coffea")),
            "Indiv_Raw" : load(os.path.join(input_dir, f"raw_templates_lj_bkg_{args.year}_{jobid}.coffea")),
            "Indiv_Symm_Smooth" : load(os.path.join(input_dir, f"symmetrized_smoothed_templates_lj_bkg_{args.year}_{jobid}.coffea")),
            "Indiv_Symm_Flat" : load(os.path.join(input_dir, f"symmetrized_flattened_templates_lj_bkg_{args.year}_{jobid}.coffea")),
            "Combined_Era_Lep_1D_Symm_Smooth" : load(os.path.join(comb_era_lep_indir, f"symmetrized_smoothed_1D_combined_year_and_lepton_templates_lj_bkg_{jobid}.coffea")),
        }
        print("Creating final background templates")
        final_bkg_templates(bkg_dict)

    if "PDF" in templates_to_run:
        analyzer = "htt_pdfUncs"
        input_dir, outdir = make_input_output_dir(analyzer)

        #set_trace()
        comb_era_lep_indir = os.path.join(eos_dir, "results", jobid, f"Templates_{analyzer}")
        if not os.path.join(comb_era_lep_indir): raise ValueError("No directory found for combined era and lepton channels.")
            # find pdf files
        pdf_dict = {
            "Indiv_Nom" : load(os.path.join(eos_dir, "results", f"{args.year}_{jobid}", f"Templates_htt_btag_sb_regions", f"raw_templates_lj_bkg_{args.year}_{jobid}.coffea")), # make sure that exact same TT_nosys templates are used for normalization as nominal ones
            "Indiv_Orig" : load(os.path.join(input_dir, f"raw_pdf_templates_lj_bkg_{args.year}_{jobid}.coffea")), # original templates
            "Combined_Era_Lep_Smooth" : load(os.path.join(comb_era_lep_indir, f"smoothed_combined_year_and_lepton_templates_lj_PDF_{jobid}.coffea")),
            "Combined_Era_Lep_Raw" : load(os.path.join(comb_era_lep_indir, f"raw_combined_year_and_lepton_templates_lj_PDF_{jobid}.coffea"))
        }
        print("Creating final PDF templates")
        final_pdf_templates(pdf_dict)


    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
