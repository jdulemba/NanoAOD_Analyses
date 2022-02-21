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
import Utilities.prettyjson as prettyjson

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("--only_bkg", action="store_true", help="Make background templates only.")
parser.add_argument("--only_sig", action="store_true", help="Make signal templates only.")
parser.add_argument("--kfactors", action="store_true", help="Apply signal k-factors to signal")
parser.add_argument("--scale_mtop3gev", action="store_true", help="Scale 3GeV mtop variations by 1/6")
args = parser.parse_args()


def check_templates(fname, isSignal=False):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """
    hdict = load(fname)
    out_dict = {njets : {"Muon" : {}, "Electron" :{}} for njets in hdict.keys()}
    histo_dict = processor.dict_accumulator({njets : {"Muon" : {}, "Electron" :{}} for njets in hdict.keys()})
    for jmult in hdict.keys():
        for lep, histo in hdict[jmult].items():
            for sys in systypes:
                up_sysname = systematics.sys_groups[args.year][sys][0]
                dw_sysname = systematics.sys_groups[args.year][sys][1]
                procs_sys = sorted(set([key.split(f"_{up_sysname}")[0] for key in histo.keys() if up_sysname in key])) if not dw_sysname \
                    else sorted(set([key.split(f"_{up_sysname}")[0] for key in histo.keys() if up_sysname in key] + [key.split(f"_{dw_sysname}")[0] for key in histo.keys() if dw_sysname in key]))

                if not procs_sys: continue

                #set_trace()
                for proc in procs_sys:
                    if not (f"{proc}_{up_sysname}" in histo.keys() and f"{proc}_{dw_sysname}" in histo.keys()): continue
                    print(jmult, lep, sys, proc)
                        # get yield vals
                    nosys_vals = histo[f"{proc}_nosys"].values()[()]
                    up_yield_vals = histo[f"{proc}_{up_sysname}"].values()[()]
                    dw_yield_vals = histo[f"{proc}_{dw_sysname}"].values()[()]

                    # create hist with symmetrized yields
                    relative_symm_yields = np.sqrt(up_yield_vals/dw_yield_vals) - 1. # find relative deviation from nominal values for symmetrized dist
                    symmetrized_yields_up = (relative_symm_yields + 1.) * nosys_vals
                    symmetrized_yields_dw = (-1.*relative_symm_yields + 1.) * nosys_vals
                    symmetrized_histo_up = histo[f"{proc}_nosys"].copy()
                    symmetrized_histo_dw = histo[f"{proc}_nosys"].copy()
                        ## fill bins
                    for xbin in range(symmetrized_yields_up.size):
                        symmetrized_histo_up.values()[()][xbin] = symmetrized_yields_up[xbin]
                        symmetrized_histo_dw.values()[()][xbin] = symmetrized_yields_dw[xbin]

                    ## save templates that need to be symmetrized to dict
                    histo_dict[jmult][lep][f"{proc}_{up_sysname}"] = symmetrized_histo_up.copy()
                    histo_dict[jmult][lep][f"{proc}_{dw_sysname}"] = symmetrized_histo_dw.copy()

                        # find relative deviation from nominal vals
                    up_rel_vals = up_yield_vals/nosys_vals - 1.
                    dw_rel_vals = dw_yield_vals/nosys_vals - 1.
                    to_symmetrize = np.sum(up_rel_vals * dw_rel_vals >= 0, axis=0) > 2
                    if to_symmetrize:
                        print("\tsymmetrize")
                        out_dict[jmult][lep][f"{proc}_{sys}"] = "symmetrize"

    #set_trace()
    if isSignal:
        coffea_out = os.path.join(input_dir, f"symmetrized_templates_lj_sig_kfactors_{args.year}_{jobid}.coffea" if args.kfactors else f"symmetrized_templates_lj_sig_{args.year}_{jobid}.coffea")
        out_3pj = os.path.join(input_dir, f"templates_to_symmetrize_lj_sig_kfactors_{args.year}_{jobid}.json" if args.kfactors else f"templates_to_symmetrize_lj_sig_{args.year}_{jobid}.json")
    else:
        coffea_out = os.path.join(input_dir, f"symmetrized_templates_lj_bkg_mtopscaled_{args.year}_{jobid}.coffea" if args.scale_mtop3gev else f"symmetrized_templates_lj_bkg_{args.year}_{jobid}.coffea")
        out_3pj = os.path.join(input_dir, f"templates_to_symmetrize_lj_bkg_mtopscaled_{args.year}_{jobid}.json" if args.scale_mtop3gev else f"templates_to_symmetrize_lj_bkg_{args.year}_{jobid}.json")
    with open(out_3pj, "w") as out:
        out.write(prettyjson.dumps(out_dict))
    print(f"{out_3pj} written")

    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")


if __name__ == "__main__":
    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]

    analyzer = "htt_btag_sb_regions"
    base_bkg_template_name = f"smoothed_templates_lj_bkg_mtopscaled_{args.year}_{jobid}.coffea" if args.scale_mtop3gev else f"smoothed_templates_lj_bkg_{args.year}_{jobid}.coffea"
    base_sig_template_name = f"smoothed_templates_lj_sig_kfactors_{args.year}_{jobid}.coffea" if args.kfactors else f"smoothed_templates_lj_sig_{args.year}_{jobid}.coffea"

    njets_to_run = ["3Jets", "4PJets"]

    systypes = systematics.sys_groups[args.year].keys()

        # get matching pattern based on args.njets
    njets_regex = "*" if len(njets_to_run) > 1 else njets_to_run[0]

    input_dir = os.path.join(proj_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}")
    if not args.only_sig:
        try:
            bkg_fname = os.path.join(input_dir, base_bkg_template_name)
            if not os.path.isfile(bkg_fname): raise ValueError("No background file found.")
            print("Checking smoothed background templates")
            check_templates(bkg_fname)
        except:
            print("Could not write background templates to file")

    if not args.only_bkg:
        try:
            sig_fname = os.path.join(input_dir, base_sig_template_name)
            if not os.path.isfile(sig_fname): raise ValueError("No background file found.")
            print("Checking smoothed signal templates")
            check_templates(sig_fname, isSignal=True)
        except:
            print("Could not write signal templates to file")

    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
