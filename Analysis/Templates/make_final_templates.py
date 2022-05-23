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
import bkg_systs_parameters as bkg_syspar   
import signal_systs_parameters as sig_syspar

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("templates_to_run", type=str, help="Choose which type of templates to run, multiple options can be input as ':' separated strings.")
args = parser.parse_args()

def final_bkg_templates(fnames):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """
    histo_dict = processor.dict_accumulator({njets : {"Muon" : {}, "Electron" :{}} for njets in njets_to_run})

    raw_hdict = load([fname for fname in fnames if sys_treatment_dict["raw"] in os.path.basename(fname)][0])
    smooth_hdict = load([fname for fname in fnames if sys_treatment_dict["smooth"] in os.path.basename(fname)][0])
    flat_hdict = load([fname for fname in fnames if sys_treatment_dict["flat"] in os.path.basename(fname)][0])
    symm_hdict = load([fname for fname in fnames if sys_treatment_dict["symm"] in os.path.basename(fname)][0])
    for jmult in raw_hdict.keys():
        for lep in raw_hdict[jmult].keys():
            for tname in raw_hdict[jmult][lep].keys():
                proc = tname.split("_")[0] if not "data_obs" in tname else "data_obs"
                sys = sorted(filter(None, tname.split(f"{proc}_")))[0]
                print(jmult, lep, sys, proc)

                #treatment = "raw"
                if sys == "nosys":
                    treatment = "raw"
                else:
                    # choose which hist to use based on what_to_do dict in systs_parameters.py'
                    sysname = "_".join(systematics.sys_to_name[args.year][sys].split("_")[:-1]) # convert sys to name used in analysis
                    #sysname = systematics.sys_to_name[args.year][sys] if sys == "EWcorrUp" else "_".join(systematics.sys_to_name[args.year][sys].split("_")[:-1]) # convert sys to name used in analysis
                    treatment = bkg_syspar.treatment[args.year][sysname][proc][f"{lep.lower()[0]}{jmult[0]}"]

                print(f"\t{treatment}")
                if treatment == "raw":
                    hist_to_use = raw_hdict[jmult][lep][tname].copy()
                elif treatment == "flat":
                    hist_to_use = flat_hdict[jmult][lep][tname].copy()
                elif treatment == "smooth":
                    hist_to_use = smooth_hdict[jmult][lep][tname].copy()
                elif treatment == "symm":
                    hist_to_use = symm_hdict[jmult][lep][tname].copy()
                else:
                    raise ValueError(f"Treatment {treatment} not supported.")

                    ## save template histos to coffea dict
                histo_dict[jmult][lep][tname] = [hist_to_use.copy(), treatment]

    #set_trace()
    coffea_out = os.path.join(outdir, f"final_templates_lj_bkg_{args.year}_{jobid}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")


def final_sig_templates(fnames):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """
    histo_dict = processor.dict_accumulator({njets : {"Muon" : {}, "Electron" :{}} for njets in njets_to_run})

    raw_hdict = load([fname for fname in fnames if sys_treatment_dict["raw"] in os.path.basename(fname)][0])
    smooth_hdict = load([fname for fname in fnames if sys_treatment_dict["smooth"] in os.path.basename(fname)][0])
    flat_hdict = load([fname for fname in fnames if sys_treatment_dict["flat"] in os.path.basename(fname)][0])
    for jmult in raw_hdict.keys():
        for lep in raw_hdict[jmult].keys():
            for tname in raw_hdict[jmult][lep].keys():
                if "Int" in tname:
                    proc = tname.split("pos_")[0]+"pos" if "pos" in tname else tname.split("neg_")[0]+"neg"
                else:
                    proc = tname.split("Res_")[0]+"Res"
                sys = sorted(filter(None, tname.split(f"{proc}_")))[0]
                print(lep, jmult, sys, proc)

                treatment = "raw"
                #if sys == "nosys":
                #    treatment = "raw"
                #else:
                #    # choose which hist to use based on what_to_do dict in systs_parameters.py
                #    sysname = systematics.sys_to_name[args.year][sys] # convert sys to name used in analysis
                #    sys_to_check = [key for key in sig_syspar.what_to_do.keys() if fnmatch.fnmatch(sysname, f"{key}*")][0] # find sysname which corresponds to name in systs_parameters.py
                #    treatment = sig_syspar.what_to_do[sys_to_check]#[f"{lep.lower()[0]}{jmult[0]}"]
                #    #treatment = sig_syspar.what_to_do[sys_to_check][proc][f"{lep.lower()[0]}{jmult[0]}"]

                #print(f"\t{treatment}")
                if treatment == "raw":
                    hist_to_use = raw_hdict[jmult][lep][tname].copy()
                #if treatment == "flat":
                #    hist_to_use = flat_hdict[jmult][lep][tname].copy()
                #if treatment == "smooth":
                #    hist_to_use = smooth_hdict[jmult][lep][tname].copy()

                    ## save template histos to coffea dict
                histo_dict[jmult][lep][tname] = [hist_to_use.copy(), treatment]

    coffea_out = os.path.join(outdir, f"final_templates_lj_sig_{args.year}_{jobid}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")


def final_MEsig_templates(fnames):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """
    histo_dict = processor.dict_accumulator({njets : {"Muon" : {}, "Electron" :{}} for njets in njets_to_run})

    raw_hdict = load([fname for fname in fnames if sys_treatment_dict["raw"] in os.path.basename(fname)][0])
    #smooth_hdict = load([fname for fname in fnames if sys_treatment_dict["smooth"] in os.path.basename(fname)][0])
    #flat_hdict = load([fname for fname in fnames if sys_treatment_dict["flat"] in os.path.basename(fname)][0])
    for jmult in raw_hdict.keys():
        for lep in raw_hdict[jmult].keys():
            for tname in raw_hdict[jmult][lep].keys():
                if "Int" in tname:
                    proc = tname.split("pos_")[0]+"pos" if "pos" in tname else tname.split("neg_")[0]+"neg"
                else:
                    proc = tname.split("Res_")[0]+"Res"
                sys = sorted(filter(None, tname.split(f"{proc}_")))[0]
                print(lep, jmult, sys, proc)

                treatment = "raw"
                #if sys == "nosys":
                #    treatment = "raw"
                #else:
                #    # choose which hist to use based on what_to_do dict in systs_parameters.py
                #    sysname = systematics.sys_to_name[args.year][sys] # convert sys to name used in analysis
                #    sys_to_check = [key for key in sig_syspar.what_to_do.keys() if fnmatch.fnmatch(sysname, f"{key}*")][0] # find sysname which corresponds to name in systs_parameters.py
                #    treatment = sig_syspar.what_to_do[sys_to_check]#[f"{lep.lower()[0]}{jmult[0]}"]
                #    #treatment = sig_syspar.what_to_do[sys_to_check][proc][f"{lep.lower()[0]}{jmult[0]}"]

                #print(f"\t{treatment}")
                if treatment == "raw":
                    hist_to_use = raw_hdict[jmult][lep][tname].copy()
                #if treatment == "flat":
                #    hist_to_use = flat_hdict[jmult][lep][tname].copy()
                #if treatment == "smooth":
                #    hist_to_use = smooth_hdict[jmult][lep][tname].copy()

                    ## save template histos to coffea dict
                histo_dict[jmult][lep][tname] = [hist_to_use.copy(), treatment]

    coffea_out = os.path.join(outdir, f"final_templates_lj_MEsig_{args.year}_{jobid}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")


def final_pdf_templates(fnames):
    """
    Function that chooses final templates.
    """
    histo_dict = processor.dict_accumulator({njets : {"Muon" : {}, "Electron" :{}} for njets in njets_to_run})

    raw_hdict = load([fname for fname in fnames if sys_treatment_dict["raw"] in os.path.basename(fname)][0])
    #smooth_hdict = load([fname for fname in fnames if sys_treatment_dict["smooth"] in os.path.basename(fname)][0])
    #flat_hdict = load([fname for fname in fnames if sys_treatment_dict["flat"] in os.path.basename(fname)][0])
    for jmult in raw_hdict.keys():
        for lep in raw_hdict[jmult].keys():
            for tname in raw_hdict[jmult][lep].keys():
                proc = tname.split("_")[0] if not "data_obs" in tname else "data_obs"
                sys = sorted(filter(None, tname.split(f"{proc}_")))[0]
                print(lep, jmult, sys, proc)

                treatment = "raw"
                #if sys == "nosys":
                #    treatment = "raw"
                #else:
                #    ## choose which hist to use based on what_to_do dict in systs_parameters.py
                #    #sysname = systematics.sys_to_name[args.year][sys] # convert sys to name used in analysis
                #    #sys_to_check = [key for key in bkg_syspar.what_to_do.keys() if fnmatch.fnmatch(sysname, f"{key}*")][0] # find sysname which corresponds to name in systs_parameters.py
                #    treatment = "raw"

                #print(f"\t{treatment}")
                if treatment == "raw":
                    hist_to_use = raw_hdict[jmult][lep][tname].copy()
                #elif treatment == "flat":
                #    hist_to_use = flat_hdict[jmult][lep][tname].copy()
                #elif treatment == "smooth":
                #    hist_to_use = smooth_hdict[jmult][lep][tname].copy()
                else:
                    raise ValueError(f"Treatment {treatment} not supported.")

                    ## save template histos to coffea dict
                histo_dict[jmult][lep][tname] = [hist_to_use.copy(), treatment]

    coffea_out = os.path.join(outdir, f"final_pdf_templates_lj_bkg_{args.year}_{jobid}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")



def make_input_output_dir(analyzer):
    input_dir = os.path.join(proj_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}")
    output_dir = os.path.join(proj_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}", "FINAL")
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    return input_dir, output_dir


if __name__ == "__main__":
    allowed_template_options = ["bkg", "sig", "MEreweight_sig", "PDF"]
    templates_to_run = [template for template in (args.templates_to_run).split(":") if template in allowed_template_options]
    templates_to_not_run = [template for template in (args.templates_to_run).split(":") if template not in allowed_template_options]
    if templates_to_not_run:
        print(f"{templates_to_not_run} are not valid options for making templates, will be skipped")

    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]

    njets_to_run = ["3Jets", "4PJets"]

    sys_treatment_dict = {
        "raw" : "raw", 
        "flat" : "flattened",
        "smooth" : "smoothed",
        "symm" : "symmetrized",
    }

    if "bkg" in templates_to_run:
        analyzer = "htt_btag_sb_regions"
        input_dir, outdir = make_input_output_dir(analyzer)

            # find background files
        base_bkg_template_name = f"*_templates_lj_bkg_{args.year}_{jobid}.coffea"
        bkg_fnames = fnmatch.filter(os.listdir(input_dir), base_bkg_template_name)
        bkg_fnames = [os.path.join(input_dir, fname) for fname in bkg_fnames]
        print("Creating final background templates")
        final_bkg_templates(bkg_fnames)

    if "sig" in templates_to_run:
        analyzer = "htt_btag_sb_regions"
        input_dir, outdir = make_input_output_dir(analyzer)

            # find signal files
        base_sig_template_name = f"*_templates_lj_sig_{args.year}_{jobid}.coffea"
        sig_fnames = fnmatch.filter(os.listdir(input_dir), base_sig_template_name)
        sig_fnames = [os.path.join(input_dir, fname) for fname in sig_fnames]
        print("Creating final signal templates")
        final_sig_templates(sig_fnames)

    if "MEreweight_sig" in templates_to_run:
        analyzer = "htt_btag_sb_regions"
        input_dir, outdir = make_input_output_dir(analyzer)

            # find signal files
        base_sig_template_name = f"*_templates_lj_MEsig_{args.year}_{jobid}.coffea"
        sig_fnames = fnmatch.filter(os.listdir(input_dir), base_sig_template_name)
        sig_fnames = [os.path.join(input_dir, fname) for fname in sig_fnames]
        print("Creating ME reweighting signal templates")
        final_MEsig_templates(sig_fnames)

    if "PDF" in templates_to_run:
        analyzer = "htt_pdfUncs"
        input_dir, outdir = make_input_output_dir(analyzer)

            # find background files
        base_pdf_template_name = f"*_pdf_templates_lj_bkg_{args.year}_{jobid}.coffea"
        pdf_fnames = fnmatch.filter(os.listdir(input_dir), base_pdf_template_name)
        pdf_fnames = [os.path.join(input_dir, fname) for fname in pdf_fnames]
        print("Creating final PDF templates")
        final_pdf_templates(pdf_fnames)


    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
