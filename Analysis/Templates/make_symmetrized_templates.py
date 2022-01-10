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
import Utilities.systematics as systematics

base_jobid = os.environ["base_jobid"]
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("treatment", choices=["flat", "smooth"], help="Specify which process to use.")
parser.add_argument("--njets", default="all", nargs="?", choices=["3", "4+", "all"], help="Specify which jet multiplicity to use.")
parser.add_argument("--only_bkg", action="store_true", help="Make background templates only.")
parser.add_argument("--only_sig", action="store_true", help="Make signal templates only.")
args = parser.parse_args()

njets_to_run = []
if (args.njets == "3") or (args.njets == "all"):
    njets_to_run += ["3Jets"]
if (args.njets == "4+") or (args.njets == "all"):
    njets_to_run += ["4PJets"]



def symmetrize_bkg_templates(fnames_to_run):
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
        systs = sorted(set(["_".join(key.split("_")[1:]) for key in sorted(hdict["Muon"].keys()) if not ("data_obs" in key or len(key.split("_")) == 1)]))
        if "nosys" in systs: systs.remove("nosys")
        systypes = sorted(set([baseSys(systematics.sys_to_name["2016APV"][sys]) for sys in systs]))
        for sys in systypes:
            up_sysname = [key for key, val in systematics.sys_to_name["2016APV"].items() if val == f"{sys}_UP"][0]
            dw_sysname = [key for key, val in systematics.sys_to_name["2016APV"].items() if val == f"{sys}_DW"][0]
            procs_sys = sorted(set([key.split(f"_{up_sysname}")[0] for key in sorted(hdict["Muon"].keys()) if up_sysname in key] + [key.split(f"_{dw_sysname}")[0] for key in sorted(hdict["Muon"].keys()) if dw_sysname in key]))
            for proc in procs_sys:
                for lep in hdict.keys():
                    print(lep, jmult, sys, proc)
                    up_histo = hdict[lep][f"{proc}_{up_sysname}"]
                    dw_histo = hdict[lep][f"{proc}_{dw_sysname}"]
                    nosys_histo = hdict[lep][f"{proc}_nosys"]

                    nosys_vals = nosys_histo.values()[()]
                    up_rel_vals = (up_histo.values()[()]-nosys_vals)/nosys_vals
                    dw_rel_vals = (dw_histo.values()[()]-nosys_vals)/nosys_vals
                    avg_rel_vals = np.sqrt((np.square(up_rel_vals)+np.square(dw_rel_vals))/2)

                    #if sys == "MTOP3GEV": set_trace()
                    symmetrized_up_vals = nosys_vals*(1 + avg_rel_vals*np.sign(up_rel_vals))
                    symmetrized_dw_vals = nosys_vals*(1 - avg_rel_vals*np.sign(up_rel_vals))
                    #symmetrized_dw_vals = nosys_vals*(1 + avg_rel_vals*np.sign(dw_rel_vals))

                    symm_up_histo = Plotter.np_array_TO_hist(sumw=symmetrized_up_vals, sumw2=np.zeros(symmetrized_up_vals.size), hist_template=nosys_histo)
                    symm_dw_histo = Plotter.np_array_TO_hist(sumw=symmetrized_dw_vals, sumw2=np.zeros(symmetrized_dw_vals.size), hist_template=nosys_histo)
                
                        ## save template histos to coffea dict
                    if jmult == "3Jets":
                        histo_dict_3j[lep][f"{proc}_{up_sysname}"] = symm_up_histo
                        histo_dict_3j[lep][f"{proc}_{dw_sysname}"] = symm_dw_histo
                        if f"{proc}_nosys" not in histo_dict_3j[lep].keys():
                            histo_dict_3j[lep][f"{proc}_nosys"] = nosys_histo
                        if "data_obs_nosys" not in histo_dict_3j[lep].keys():
                            histo_dict_3j[lep]["data_obs_nosys"] = hdict[lep]["data_obs_nosys"]
                    if jmult == "4PJets":
                        histo_dict_4pj[lep][f"{proc}_{up_sysname}"] = symm_up_histo
                        histo_dict_4pj[lep][f"{proc}_{dw_sysname}"] = symm_dw_histo
                        if f"{proc}_nosys" not in histo_dict_4pj[lep].keys():
                            histo_dict_4pj[lep][f"{proc}_nosys"] = nosys_histo
                        if "data_obs_nosys" not in histo_dict_4pj[lep].keys():
                            histo_dict_4pj[lep]["data_obs_nosys"] = hdict[lep]["data_obs_nosys"]

    #set_trace()
    if "3Jets" in njets_to_run:
        coffea_out_3j = os.path.join(input_dir, f"test_SymmetrizedSmoothed_templates_lj_3Jets_bkg_{args.year}_{jobid}.coffea" if args.treatment == "smooth" else f"test_SymmetrizedFlattened_templates_lj_3Jets_bkg_{args.year}_{jobid}.coffea")
        save(histo_dict_3j, coffea_out_3j)
        print(f"{coffea_out_3j} written")
    if "4PJets" in njets_to_run:
        coffea_out_4pj = os.path.join(input_dir, f"test_SymmetrizedSmoothed_templates_lj_4PJets_bkg_{args.year}_{jobid}.coffea" if args.treatment == "smooth" else f"test_SymmetrizedFlattened_templates_lj_4PJets_bkg_{args.year}_{jobid}.coffea")
        save(histo_dict_4pj, coffea_out_4pj)
        print(f"{coffea_out_4pj} written")


def symmetrize_sig_templates(fnames_to_run):
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
        systs = sorted(set([key.split("Res_")[1] for key in sorted(hdict["Muon"].keys()) if "Res" in key]))
        if "nosys" in systs: systs.remove("nosys")
        systypes = sorted(set([baseSys(systematics.sys_to_name["2016APV"][sys]) for sys in systs]))
        for sys in systypes:
            up_sysname = [key for key, val in systematics.sys_to_name["2016APV"].items() if val == f"{sys}_UP"][0]
            dw_sysname = [key for key, val in systematics.sys_to_name["2016APV"].items() if val == f"{sys}_DW"][0]
            procs_sys = sorted(set([key.split(f"_{up_sysname}")[0] for key in sorted(hdict["Muon"].keys()) if up_sysname in key] + [key.split(f"_{dw_sysname}")[0] for key in sorted(hdict["Muon"].keys()) if dw_sysname in key]))
            for proc in procs_sys:
                for lep in hdict.keys():
                    print(lep, jmult, sys, proc)
                    up_histo = hdict[lep][f"{proc}_{up_sysname}"]
                    dw_histo = hdict[lep][f"{proc}_{dw_sysname}"]
                    nosys_histo = hdict[lep][f"{proc}_nosys"]

                    nosys_vals = nosys_histo.values()[()]
                    up_rel_vals = (up_histo.values()[()]-nosys_vals)/nosys_vals
                    dw_rel_vals = (dw_histo.values()[()]-nosys_vals)/nosys_vals
                    avg_rel_vals = np.sqrt((np.square(up_rel_vals)+np.square(dw_rel_vals))/2)

                    symmetrized_up_vals = nosys_vals*(1 + avg_rel_vals)
                    symmetrized_dw_vals = nosys_vals*(1 - avg_rel_vals)

                    symm_up_histo = Plotter.np_array_TO_hist(sumw=symmetrized_up_vals, sumw2=np.zeros(symmetrized_up_vals.size), hist_template=nosys_histo)
                    symm_dw_histo = Plotter.np_array_TO_hist(sumw=symmetrized_dw_vals, sumw2=np.zeros(symmetrized_dw_vals.size), hist_template=nosys_histo)
                
                        ## save template histos to coffea dict
                    if jmult == "3Jets":
                        histo_dict_3j[lep][f"{proc}_{up_sysname}"] = symm_up_histo
                        histo_dict_3j[lep][f"{proc}_{dw_sysname}"] = symm_dw_histo
                        if f"{proc}_nosys" not in histo_dict_3j[lep].keys():
                            histo_dict_3j[lep][f"{proc}_nosys"] = nosys_histo
                    if jmult == "4PJets":
                        histo_dict_4pj[lep][f"{proc}_{up_sysname}"] = symm_up_histo
                        histo_dict_4pj[lep][f"{proc}_{dw_sysname}"] = symm_dw_histo
                        if f"{proc}_nosys" not in histo_dict_4pj[lep].keys():
                            histo_dict_4pj[lep][f"{proc}_nosys"] = nosys_histo

    if "3Jets" in njets_to_run:
        coffea_out_3j = os.path.join(input_dir, f"test_SymmetrizedSmoothed_templates_lj_3Jets_sig_{args.year}_{jobid}.coffea" if args.treatment == "smooth" else f"test_SymmetrizedFlattened_templates_lj_3Jets_sig_{args.year}_{jobid}.coffea")
        save(histo_dict_3j, coffea_out_3j)
        print(f"{coffea_out_3j} written")
    if "4PJets" in njets_to_run:
        coffea_out_4pj = os.path.join(input_dir, f"test_SymmetrizedSmoothed_templates_lj_4PJets_sig_{args.year}_{jobid}.coffea" if args.treatment == "smooth" else f"test_SymmetrizedFlattened_templates_lj_4PJets_sig_{args.year}_{jobid}.coffea")
        save(histo_dict_4pj, coffea_out_4pj)
        print(f"{coffea_out_4pj} written")


if __name__ == "__main__":

    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]

    analyzer = "htt_btag_sb_regions"
    base_bkg_template_name = f"test_smoothed_templates_lj_NJETS_bkg_{args.year}_{jobid}" if args.treatment == "smooth" else f"test_flattened_templates_lj_NJETS_bkg_{args.year}_{jobid}"
    base_sig_template_name = f"test_smoothed_templates_lj_NJETS_sig_{args.year}_{jobid}" if args.treatment == "smooth" else f"test_flattened_templates_lj_NJETS_sig_{args.year}_{jobid}"
    #base_bkg_template_name = f"test_CombinedSmoothed_templates_lj_NJETS_bkg_{jobid}" if args.treatment == "smooth" else f"test_CombinedFlattened_templates_lj_NJETS_bkg_{jobid}"
    #base_sig_template_name = f"test_CombinedSmoothed_templates_lj_NJETS_sig_{jobid}" if args.treatment == "smooth" else f"test_CombinedFlattened_templates_lj_NJETS_sig_{jobid}"

        # get matching pattern based on args.njets
    njets_regex = "*" if len(njets_to_run) > 1 else njets_to_run[0]

    input_dir = os.path.join(proj_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}")
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

    baseSys = lambda sys : "_".join(sys.split("_")[:-1])

    if not args.only_sig:
        try:
            print("Creating symmetrized background templates")
            symmetrize_bkg_templates(bkg_fnames)
            
        except:
            print("Could not write background templates to file")

    if not args.only_bkg:
        try:
            print("Creating symmetrized signal templates")
            symmetrize_sig_templates(sig_fnames)
            
        except:
            print("Could not write signal templates to file")

    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
