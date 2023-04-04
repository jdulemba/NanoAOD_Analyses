#! /bin/env python

import time
tic = time.time()

from pdb import set_trace
import os
import fnmatch
import numpy as np
import Utilities.systematics as systematics
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
parser.add_argument("--MEopts", nargs="*", action=ParseKwargs, help="Options to specify which mass and width points to produce for ME reweighted signal.")
args = parser.parse_args()


def final_MEsig_templates(hdicts, year, allowed_mass):
    raw_hdict = hdicts["Indiv"]
    for jmult in njets_to_run:
        for lep in leps_to_run:
            orig_lepdir = "muNJETS" if lep == "Muon" else "eNJETS"
            lepdir = orig_lepdir.replace("NJETS", jmult.lower())

            if year == "2016APV": year_to_use = "2016pre"
            elif year == "2016": year_to_use = "2016post"
            else: year_to_use = year
            dirname = f"{lepdir}_{year_to_use}"
            upfout.mkdir(dirname)

            sig_to_use = sorted(set([name.split(";")[0] for name in raw_hdict[jmult][lep].keys() if any([fnmatch.fnmatch(name, f"*{allowed_mass}*{width}*") for width in allowed_widths])]))
            for tname in sig_to_use:
                if "Int" in tname:
                    proc = tname.split("pos_")[0]+"pos" if "pos" in tname else tname.split("neg_")[0]+"neg"
                    boson, mass, width, pI, wt = proc.split("_")
                else:
                    proc = tname.split("Res_")[0]+"Res"
                    boson, mass, width, pI = proc.split("_")

                if ((mass != allowed_mass) or (width not in allowed_widths)): continue
                #set_trace()
                sub_name = "_".join([boson[0], mass.lower(), width.lower(), wt]) if pI == "Int" else "_".join([boson[0], mass.lower(), width.lower(), pI.lower()])
                sys = sorted(filter(None, tname.split(f"{proc}_")))[0]

                #if sys == "mtop171p5": set_trace()
                sysname = "nosys" if sys == "nosys" else "_".join(systematics.sys_to_name[year][sys].split("_")[:-1]) # convert sys to name used in analysis
                if sysname in comb_era_lep_systypes:
                    treatment = "Combined_Era_Lep"
                    ratio_hvals = hdicts["Combined_Era_Lep"][jmult][tname].values()
                    hist_to_use = hdicts["Indiv"][jmult][lep][f"{proc}_nosys"].to_hist().copy()
                    sys_vals = hist_to_use.values() * ratio_hvals
                    hist_to_use.values()[:] = sys_vals

                elif sysname in comb_lep_systypes:
                    treatment = "Combined_Lep"
                    ratio_hvals = hdicts["Combined_Lep"][jmult][tname].values()
                    hist_to_use = hdicts["Indiv"][jmult][lep][f"{proc}_nosys"].to_hist().copy()
                    sys_vals = hist_to_use.values() * ratio_hvals
                    hist_to_use.values()[:] = sys_vals

                else:
                    treatment = "Indiv"
                    hist_to_use = hdicts["Indiv"][jmult][lep][tname].to_hist().copy()

                print(f"{year}, {jmult}, {lep}, {sys}, {proc}:\t{treatment}")

                    # get proper name for output hist
                outhname = sub_name if sys == "nosys" else "_".join(list(filter(None, [sub_name, systematics.combine_template_sys_to_name[year][sys]])))
                if "LEP" in outhname: outhname = outhname.replace("LEP", lep[0].lower())
                if "_tot" in outhname: outhname = outhname.replace("_tot", "")
                    # replace '2016APV' with '2016pre' and '2016' with '2016post' for files sent to Afiq
                if (year == "2016APV") and ("2016APV" in outhname):
                    #set_trace()
                    outhname = outhname.replace("2016APV", "2016pre")
                if (year == "2016") and ("2016" in outhname):
                    outhname = outhname.replace("2016", "2016post")

                if np.any(np.isnan(hist_to_use.values())):
                    print(f"\tSome bins in template contain nan values!")
                    set_trace()

                    ## save templates to root file
                upfout[dirname][outhname] = hist_to_use



if __name__ == "__main__":

    possible_masses = [365, 380, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000]
    possible_widths = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 8.0, 10.0, 13.0, 15.0, 18.0, 21.0, 25.0]
    possible_masses = [str(mass) for mass in possible_masses]
    possible_widths = [str(width) for width in possible_widths]

    #set_trace()
    masses_to_run = "All" if not args.MEopts else args.MEopts.get("allowed_masses", "All")
    widths_to_run = "All" if not args.MEopts else args.MEopts.get("allowed_widths", "All")

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
    #set_trace()

    jobid = os.environ["jobid"]
    eos_dir = os.environ["eos_dir"]

    njets_to_run = sorted(["3Jets", "4PJets"])
    leps_to_run = sorted(["Muon", "Electron"])

    analyzer = "htt_btag_sb_regions"

    comb_era_lep_indir = os.path.join(eos_dir, "results", jobid, f"Templates_{analyzer}")
    if not os.path.join(comb_era_lep_indir): raise ValueError("No directory found for combined era and lepton channels.")

    comb_era_lep_systypes = ["AH_ISR", "AH_FSR", "AH_RENORM", "AH_FACTOR", "AH_MTOP"]

    years_to_run = ["2016APV", "2016", "2017", "2018"]

    #set_trace()
    version = "V31"
    #version = "V30"
    #version = "V29"
    #version = "V27"
    #version = "V26"
    #version = "V25"

    if version != "V31":
        set_trace()
    output_dir = f"root://cmseos.fnal.gov//store/user/jdulemba/Htt_Templates/Summer20UL_DeepJet_NLOshapeTT_{version.lower()}"

    for mass in allowed_masses:
        print(f"Mass point: {mass}")

            ## make name of output root file
        rname = f"templates_lj_sig_{mass.lower()}.root" if (widths_to_run == "All") else f"templates_lj_sig_{mass.lower()}_{''.join(allowed_widths).lower()}.root"
        upfout = uproot.recreate(rname, compression=uproot.ZLIB(4)) if os.path.isfile(rname) else uproot.create(rname)

        for year in years_to_run:
                # define systematic groups for both types of combining templates
            comb_lep_systypes = [
                "BTAG_BC_JES", "BTAG_BC_PILEUP", "BTAG_BC_STATISTIC", "BTAG_BC_TYPE3",
                "BTAG_BC_CORR", "BTAG_BC_UNCORR", "BTAG_L_CORR", "BTAG_L_UNCORR",
                "BTAG_BC_BFRAGMENTATION", "BTAG_BC_BTEMPCORR", "BTAG_BC_CB", "BTAG_BC_CFRAGMENTATION", "BTAG_BC_CJETS",
                "BTAG_BC_DMUX", "BTAG_BC_GLUONSPLITTING", "BTAG_BC_JETAWAY", "BTAG_BC_KSL", "BTAG_BC_L2C",
                "BTAG_BC_LTOTHERS", "BTAG_BC_MUDR", "BTAG_BC_MUPT", "BTAG_BC_PTREL",
                "JES_Absolute", f"JES_Absolute_{year}", "JES_BBEC1", f"JES_BBEC1_{year}",
                "JES_FlavorQCD", "JES_FlavorQCDOnlyLightJets",
                "JES_FlavorPureBottom", "JES_FlavorPureCharm", "JES_FlavorPureGluon", "JES_FlavorPureQuark",
                "JES_FlavorPureBottomOnlyBottomJets", "JES_FlavorPureCharmOnlyCharmJets", "JES_FlavorPureGluonOnlyGluonJets", "JES_FlavorPureQuarkOnlyQuarkJets",
                "JES_RelativeBal", f"JES_RelativeSample_{year}", "MET", "JER"
            ]
            if year != "2018": comb_lep_systypes.append("PREFIRE")

            input_dir = os.path.join(eos_dir, "results", f"{year}_{jobid}", f"Templates_{analyzer}")
            #set_trace()
                # find signal files
            MEsig_dict = {
                "Combined_Era_Lep" : uproot.open("root://eosuser.cern.ch/" + os.path.join(comb_era_lep_indir, f"average_combined_year_and_lepton_ratios_lj_MEsig_{jobid}.root")),
                "Combined_Lep" : uproot.open("root://eosuser.cern.ch/" + os.path.join(input_dir, f"average_combined_lepton_ratios_lj_MEsig_{year}_{jobid}.root")),
                "Indiv" : uproot.open("root://eosuser.cern.ch/" + os.path.join(input_dir, f"raw_templates_lj_MEsig_{year}_{jobid}_TOT.root")),
            }
            final_MEsig_templates(MEsig_dict, year, mass)

        upfout.close()
        print(f"{rname} written")

            ## copy output file to input dir and then delete local copy
        exec_str = f"xrdcp {rname} {output_dir} && rm {rname}"
        try:
            print(f"Executing {exec_str}")
            os.system(exec_str)
        except:
            raise ValueError(f"Could not copy {rname} to {output_dir}")

    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
