#! /bin/env python

import time
tic = time.time()

from pdb import set_trace
import os
import Utilities.systematics as systematics
import numpy as np
import uproot

import argparse
class ParseKwargs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, dict())
        for value in values:
            key, value = value.split("=")
            getattr(namespace, self.dest)[key] = value
    
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("channels_to_combine", type=str, help="Choose how to combine systematic templates (by lepton or eras), multiple options can be input as ':' separated strings.")
parser.add_argument("--MEopts", nargs="*", action=ParseKwargs, help="Options to specify which mass and width points to produce for ME reweighted signal.")
args = parser.parse_args()


def average_lepton_templates(fname, year):
    rfile = uproot.open(fname)

    for jmult in njets_to_run:
        upfout.mkdir(jmult)
        for sys in systypes:
            print(year, jmult, sys)
            #set_trace()
            up_sysname = systematics.sys_groups[year][sys][0]
            dw_sysname = systematics.sys_groups[year][sys][1]
                # filter out masses and widths for procs
            procs_sys = list(filter(lambda proc: (any([mass for mass in allowed_masses if mass in proc]) and any([width for width in allowed_widths if width in proc])),
                sorted(set([key.split(f"_{up_sysname}")[0] for key in rfile[f"/{jmult}/Muon"].keys() if up_sysname in key])) if not dw_sysname \
                else sorted(set([key.split(f"_{up_sysname}")[0] for key in rfile[f"/{jmult}/Muon"].keys() if up_sysname in key] + [key.split(f"_{dw_sysname}")[0] for key in rfile[f"/{jmult}/Muon"].keys() if dw_sysname in key]))
            ))

            #if not procs_sys: continue
            for proc in procs_sys:
                average_combined_lep_up_ratio = None
                average_combined_lep_dw_ratio = None
                idx = 0
                for lep in leps_to_run:
                    #print(year, jmult, sys, proc, lep)
                    up_ratio_vals = rfile[f"{jmult}/{lep}"][f"{proc}_{up_sysname}"].values()/rfile[f"{jmult}/{lep}"][f"{proc}_nosys"].values()
                    dw_ratio_vals = rfile[f"{jmult}/{lep}"][f"{proc}_{dw_sysname}"].values()/rfile[f"{jmult}/{lep}"][f"{proc}_nosys"].values()
                    if idx == 0:
                        average_combined_lep_up_ratio = up_ratio_vals
                        average_combined_lep_dw_ratio = dw_ratio_vals
                    else:
                        average_combined_lep_up_ratio += up_ratio_vals
                        average_combined_lep_dw_ratio += dw_ratio_vals
                    idx += 1

                #set_trace()
                    # actually average the ratios now
                average_combined_lep_up_ratio /= idx
                average_combined_lep_dw_ratio /= idx

                # convert values to hist
                    # up histo
                avg_combined_lep_up_histo = rfile[f"{jmult}/{lep}"][f"{proc}_{up_sysname}"].to_hist().copy()
                avg_combined_lep_up_histo.values()[:], avg_combined_lep_up_histo.variances()[:] = average_combined_lep_up_ratio, np.zeros(average_combined_lep_up_ratio.size)
                    # dw histo
                avg_combined_lep_dw_histo = rfile[f"{jmult}/{lep}"][f"{proc}_{dw_sysname}"].to_hist().copy()
                avg_combined_lep_dw_histo.values()[:], avg_combined_lep_dw_histo.variances()[:] = average_combined_lep_dw_ratio, np.zeros(average_combined_lep_dw_ratio.size)
                    ## save templates to root file
                upfout[jmult][f"{proc}_{up_sysname}"] = avg_combined_lep_up_histo
                upfout[jmult][f"{proc}_{dw_sysname}"] = avg_combined_lep_dw_histo

    upfout.close()
    print(f"{rname} written")



def average_year_and_lepton_templates(fnames_dict):
    #set_trace()
    varied_sysnames = {sys : [systematics.sys_groups["2018"][sys][0], systematics.sys_groups["2018"][sys][1]] for sys in systypes}
    procs_sys = None
    years_hdict = {year : {jmult : {lep : {} for lep in leps_to_run} for jmult in njets_to_run} for year in fnames_dict.keys()}

    #set_trace()
    for idx in tqdm(range(len(list(fnames_dict.keys())))):
        year = list(fnames_dict.keys())[idx]
        print(f"Loading {year}")
        rfile = uproot.open(fnames_dict[year])
        if idx == 0:
            up_sysname = sorted(varied_sysnames.values())[0][0]
            dw_sysname = sorted(varied_sysnames.values())[0][1]
                # filter out masses and widths for procs
            procs_sys = list(filter(lambda proc: (any([mass for mass in allowed_masses if mass in proc]) and any([width for width in allowed_widths if width in proc])),
                sorted(set([key.split(f"_{up_sysname}")[0] for key in rfile["/3Jets/Muon"].keys() if up_sysname in key])) if not dw_sysname \
                else sorted(set([key.split(f"_{up_sysname}")[0] for key in rfile["/3Jets/Muon"].keys() if up_sysname in key] + [key.split(f"_{dw_sysname}")[0] for key in rfile["/3Jets/Muon"].keys() if dw_sysname in key]))
            ))

        for jmult in njets_to_run:
            for lep in leps_to_run:
                for proc in procs_sys:
                    years_hdict[year][jmult][lep][f"{proc}_nosys"] = rfile[f"{jmult}/{lep}/{proc}_nosys"].to_hist().copy()
                    for (up_name, dw_name) in varied_sysnames.values():
                        years_hdict[year][jmult][lep][f"{proc}_{up_name}"] = rfile[f"{jmult}/{lep}/{proc}_{up_name}"].to_hist().copy()
                        years_hdict[year][jmult][lep][f"{proc}_{dw_name}"] = rfile[f"{jmult}/{lep}/{proc}_{dw_name}"].to_hist().copy()

    #set_trace()
    for jmult in njets_to_run:
        upfout.mkdir(jmult)
        for sys, (up_sysname, dw_sysname) in varied_sysnames.items():
            for proc in procs_sys:
                average_combined_year_lep_up_ratio = None
                average_combined_year_lep_dw_ratio = None
                idx = 0
                for year in years_hdict.keys():
                    for lep in leps_to_run:
                        print(jmult, sys, proc, year, lep)
                        up_ratio_vals = years_hdict[year][jmult][lep][f"{proc}_{up_sysname}"].values()/years_hdict[year][jmult][lep][f"{proc}_nosys"].values()
                        dw_ratio_vals = years_hdict[year][jmult][lep][f"{proc}_{dw_sysname}"].values()/years_hdict[year][jmult][lep][f"{proc}_nosys"].values()
                        if idx == 0:
                            average_combined_year_lep_up_ratio = up_ratio_vals
                            average_combined_year_lep_dw_ratio = dw_ratio_vals
                        else:
                            average_combined_year_lep_up_ratio += up_ratio_vals
                            average_combined_year_lep_dw_ratio += dw_ratio_vals
                        idx += 1

                #set_trace()
                    # actually average the ratios now
                average_combined_year_lep_up_ratio /= idx
                average_combined_year_lep_dw_ratio /= idx

                # convert values to hist
                    # up histo
                avg_combined_year_lep_up_histo = years_hdict[year][jmult][lep][f"{proc}_{up_sysname}"].copy()
                avg_combined_year_lep_up_histo.values()[:], avg_combined_year_lep_up_histo.variances()[:] = average_combined_year_lep_up_ratio, np.zeros(average_combined_year_lep_up_ratio.size)
                    # dw histo
                avg_combined_year_lep_dw_histo = years_hdict[year][jmult][lep][f"{proc}_{dw_sysname}"].copy()
                avg_combined_year_lep_dw_histo.values()[:], avg_combined_year_lep_dw_histo.variances()[:] = average_combined_year_lep_dw_ratio, np.zeros(average_combined_year_lep_dw_ratio.size)
                    ## save templates to root file
                upfout[jmult][f"{proc}_{up_sysname}"] = avg_combined_year_lep_up_histo
                upfout[jmult][f"{proc}_{dw_sysname}"] = avg_combined_year_lep_dw_histo

    upfout.close()
    print(f"{rname} written")



if __name__ == "__main__":
    allowed_combination_options = ["lepton", "era_lepton"]
    combinations_to_run = (args.channels_to_combine).split(":")
    combinations_to_run = [combination for combination in combinations_to_run if combination in allowed_combination_options]

    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]
    eos_dir = os.environ["eos_dir"]

    njets_to_run = sorted(["3Jets", "4PJets"])
    leps_to_run = sorted(["Muon", "Electron"])
    analyzer = "htt_btag_sb_regions"

    possible_masses = [365, 380, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000]
    possible_widths = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 8.0, 10.0, 13.0, 15.0, 18.0, 21.0, 25.0]
    possible_masses = [str(mass) for mass in possible_masses]
    possible_widths = [str(width) for width in possible_widths]

    masses_to_run = "All" if args.MEopts is None else args.MEopts.get("allowed_masses", "All")
    widths_to_run = "All" if args.MEopts is None else args.MEopts.get("allowed_widths", "All")

    if masses_to_run == "All":
        allowed_masses = possible_masses
    else:
        allowed_masses = masses_to_run.split(":")
        allowed_masses = [mass for mass in allowed_masses if mass in possible_masses]

    if widths_to_run == "All":
        allowed_widths = possible_widths
    else:
        allowed_widths = widths_to_run.split(":")
        allowed_widths = [width for width in allowed_widths if width in possible_widths]

    allowed_masses = ["M"+mass for mass in allowed_masses]
    allowed_widths = ["W"+width.replace(".", "p") for width in allowed_widths]
    print(f"Making target signal points for {allowed_masses} and {allowed_widths}")

    for combination in combinations_to_run:
        if combination == "lepton":
            for year in ["2016APV", "2016", "2017", "2018"]:
                input_dir = os.path.join(eos_dir, "results", f"{year}_{jobid}", f"Templates_{analyzer}")
                if not os.path.isdir(input_dir): raise ValueError("No directory found.")

                systypes = [
                    "BTAG_BC_JES", "BTAG_BC_PILEUP", "BTAG_BC_STATISTIC", "BTAG_BC_TYPE3",
                    "BTAG_BC_BFRAGMENTATION", "BTAG_BC_BTEMPCORR", "BTAG_BC_CB", "BTAG_BC_CFRAGMENTATION", "BTAG_BC_CJETS",
                    "BTAG_BC_DMUX", "BTAG_BC_GLUONSPLITTING", "BTAG_BC_JETAWAY", "BTAG_BC_KSL", "BTAG_BC_L2C",
                    "BTAG_BC_LTOTHERS", "BTAG_BC_MUDR", "BTAG_BC_MUPT", "BTAG_BC_PTREL",
                    "BTAG_BC_CORR", "BTAG_BC_UNCORR", "BTAG_L_CORR", "BTAG_L_UNCORR",
                    "JES_Absolute", f"JES_Absolute_{year}", "JES_BBEC1", f"JES_BBEC1_{year}",
                    "JES_FlavorQCD", "JES_FlavorQCDOnlyLightJets",
                    "JES_FlavorPureBottom", "JES_FlavorPureCharm", "JES_FlavorPureGluon", "JES_FlavorPureQuark",
                    "JES_FlavorPureBottomOnlyBottomJets", "JES_FlavorPureCharmOnlyCharmJets", "JES_FlavorPureGluonOnlyGluonJets", "JES_FlavorPureQuarkOnlyQuarkJets",
                    "JES_RelativeBal", f"JES_RelativeSample_{year}", "MET", "JER"
                ]
                #"PILEUP"]
                if year != "2018": systypes.append("PREFIRE")

                if masses_to_run == "All":
                    rname = os.path.join(input_dir, f"average_combined_lepton_ratios_lj_MEsig_{year}_{jobid}.root" if (widths_to_run == "All") else f"average_combined_lepton_ratios_lj_MEsig_{''.join(allowed_widths).lower()}_{year}_{jobid}.root")
                elif widths_to_run == "All":
                    rname = os.path.join(input_dir, f"average_combined_lepton_ratios_lj_MEsig_{year}_{jobid}.root" if (masses_to_run == "All") else f"average_combined_lepton_ratios_lj_MEsig_{''.join(allowed_masses).lower()}_{year}_{jobid}.root")
                else:
                    rname = os.path.join(input_dir, f"average_combined_lepton_ratios_lj_MEsig_{''.join(allowed_masses).lower()}_{''.join(allowed_widths).lower()}_{year}_{jobid}.root")

                print(f"creating {rname}")
                upfout = uproot.recreate(rname, compression=uproot.ZLIB(4)) if os.path.isfile(rname) else uproot.create(rname)

                sig_fname = os.path.join(eos_dir, "results", f"{year}_{jobid}", f"Templates_{analyzer}", f"raw_templates_lj_MEsig_{year}_{jobid}_TOT.root")
                if not os.path.isfile(sig_fname): raise ValueError(f"{sig_fname} not found")
                print(f"Combining e+mu channels in {year} for ME reweighting signal templates")
                average_lepton_templates("root://eosuser.cern.ch/" + sig_fname, year)


        if combination == "era_lepton":
            from tqdm import tqdm
            years_to_combine = ["2016APV", "2016", "2017", "2018"]

            output_dir = os.path.join(eos_dir, "results", jobid, f"Templates_{analyzer}")
            if not os.path.isdir(output_dir): os.makedirs(output_dir)

            systypes = ["AH_ISR", "AH_FSR", "AH_RENORM", "AH_FACTOR", "AH_MTOP"]

            if masses_to_run == "All":
                rname = os.path.join(output_dir, f"average_combined_year_and_lepton_ratios_lj_MEsig_{jobid}.root" if (widths_to_run == "All") else f"average_combined_year_and_lepton_ratios_lj_MEsig_{''.join(allowed_widths).lower()}_{jobid}.root")
            elif widths_to_run == "All":
                rname = os.path.join(output_dir, f"average_combined_year_and_lepton_ratios_lj_MEsig_{jobid}.root" if (masses_to_run == "All") else f"average_combined_year_and_lepton_ratios_lj_MEsig_{''.join(allowed_masses).lower()}_{jobid}.root")
            else:
                rname = os.path.join(output_dir, f"average_combined_year_and_lepton_ratios_lj_MEsig_{''.join(allowed_masses).lower()}_{''.join(allowed_widths).lower()}_{jobid}.root")
        
            print(f"creating {rname}")
            upfout = uproot.recreate(rname, compression=uproot.ZLIB(4)) if os.path.isfile(rname) else uproot.create(rname)

            #set_trace()
            sig_fnames_dict = {}
            for year in years_to_combine:
                sig_fname = os.path.join(eos_dir, "results", f"{year}_{jobid}", f"Templates_{analyzer}", f"raw_templates_lj_MEsig_{year}_{jobid}_TOT.root")
                if not os.path.isfile(sig_fname): raise ValueError(f"{sig_fname} not found")
                sig_fnames_dict[year] = "root://eosuser.cern.ch/" + sig_fname


            print(f"Combining channels across years for ME reweighting signal templates")
            average_year_and_lepton_templates(sig_fnames_dict)


    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
