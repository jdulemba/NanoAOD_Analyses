import time
tic = time.time()

from coffea.util import load, save
from pdb import set_trace
import os
import numpy as np
import coffea.processor as processor
import Utilities.systematics as systematics

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("channels_to_combine", type=str, help="Choose how to combine systematic templates (by lepton or eras), multiple options can be input as ':' separated strings.")
parser.add_argument("sys_treatment", choices=["smoothed", "flattened"], help="Choose which systematic treatment to symmetrize.")
parser.add_argument("--nomSMTTxsec", action="store_true", help="Apply nominal SM cross sections to top mass and LHE scale weights")
args = parser.parse_args()


def symm_year_and_lepton_templates(fname, process):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """
    hdict = load(fname)
    histo_dict = processor.dict_accumulator({njets : {} for njets in hdict.keys()})
    #set_trace()
    for jmult in hdict.keys():
        for sys in systypes:
            up_sysname = systematics.sys_groups["2017"][sys][0]
            dw_sysname = systematics.sys_groups["2017"][sys][1] if len(systematics.sys_groups["2017"][sys]) > 1 else None
            procs_sys = sorted(set([key.split(f"_{up_sysname}")[0] for key in hdict[jmult].keys() if up_sysname in key])) if not dw_sysname \
                else sorted(set([key.split(f"_{up_sysname}")[0] for key in hdict[jmult].keys() if up_sysname in key] + [key.split(f"_{dw_sysname}")[0] for key in hdict[jmult].keys() if dw_sysname in key]))

            if not procs_sys: continue

            for proc in procs_sys:
                    # skip this process/systematic if any of the variations or nominal name aren't found
                if (up_sysname is not None) & (dw_sysname is not None):
                    exists = ( (f"{proc}_{up_sysname}" in hdict[jmult].keys() and f"{proc}_{dw_sysname}" in hdict[jmult].keys()) and f"{proc}_nosys" in hdict[jmult].keys() )
                    if not exists: continue
                    print(jmult, sys, proc)

                        # get ratios
                    symmetrized_histo_up = hdict[jmult][f"{proc}_{up_sysname}"].copy()
                    symmetrized_histo_dw = hdict[jmult][f"{proc}_{dw_sysname}"].copy()
                    up_ratio_vals = symmetrized_histo_up.values()[()]
                    dw_ratio_vals = symmetrized_histo_dw.values()[()]

                        ## hard code offset of single top uR and uF uncs for TB
                    if ((sys.startswith("ST_RENORM")) or (sys.startswith("ST_FACTOR"))) and (proc == "TB"):
                        #set_trace()
                        up_ratio_vals += 0.5
                        dw_ratio_vals += 0.5

                        # find relative deviation from nominal vals
                    up_rel_vals = up_ratio_vals - 1.
                    dw_rel_vals = dw_ratio_vals - 1.

                    # create hist with symmetrized yields
                    relative_symm_yields = np.sqrt(up_ratio_vals/dw_ratio_vals) - 1. # find relative deviation from nominal values for symmetrized dist
                    symmetrized_ratios_up = np.nan_to_num((relative_symm_yields + 1.))
                    symmetrized_ratios_dw = np.nan_to_num((-1.*relative_symm_yields + 1.))

                        # check where product of relative bins is positive, i.e. the deviations are one-sided
                    rel_prods = up_rel_vals * dw_rel_vals >= 0
                        # symmetrize all bins
                    one_sided_bins = np.arange(symmetrized_ratios_up.size) if (np.sum(rel_prods, axis=0) > 10) else np.where(rel_prods)
                    #if (one_sided_bins)[0].size > 0: set_trace()
                    symmetrized_histo_up.values()[()][one_sided_bins] = symmetrized_ratios_up[one_sided_bins]
                    symmetrized_histo_dw.values()[()][one_sided_bins] = symmetrized_ratios_dw[one_sided_bins]
                    if (one_sided_bins)[0].size > 0: print("\tsymmetrize")

                if (up_sysname is not None) & (dw_sysname is None):
                    exists = f"{proc}_{up_sysname}" in hdict[jmult].keys() and f"{proc}_nosys" in hdict[jmult].keys()
                    if not exists: continue
                    print(jmult, sys, proc)

                        # get ratios
                    symmetrized_histo_up = hdict[jmult][f"{proc}_{up_sysname}"].copy()
                    symmetrized_histo_dw = symmetrized_histo_up.copy()
                    up_ratio_vals = symmetrized_histo_up.values()[()]

                        # find relative deviation from nominal vals
                    up_rel_vals = up_ratio_vals - 1.
                    dw_rel_vals = -1. * up_rel_vals

                    symmetrized_ratios_up = up_rel_vals + 1.
                    symmetrized_ratios_dw = dw_rel_vals + 1.
                    symmetrized_histo_up.values()[()][:] = symmetrized_ratios_up
                    symmetrized_histo_dw.values()[()][:] = symmetrized_ratios_dw
                    print("\tsymmetrize")
                    #set_trace()
                    dw_sysname = f"{up_sysname}Down"
                    up_sysname = f"{up_sysname}Up"

                ## save templates that need to be symmetrized to dict
                histo_dict[jmult][f"{proc}_{up_sysname}"] = symmetrized_histo_up.copy()
                histo_dict[jmult][f"{proc}_{dw_sysname}"] = symmetrized_histo_dw.copy()

    #set_trace()
    coffea_out = os.path.join(input_dir, f"symmetrized_{args.sys_treatment}_combined_year_and_lepton_templates_lj_{process}_nomSMTTxsec_{jobid}.coffea" if args.nomSMTTxsec \
            else f"symmetrized_{args.sys_treatment}_combined_year_and_lepton_templates_lj_{process}_{jobid}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")



def symm_lepton_templates(fname, process):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """
    hdict = load(fname)
    histo_dict = processor.dict_accumulator({njets : {} for njets in hdict.keys()})
    #set_trace()
    for jmult in hdict.keys():
        for sys in systypes:
            up_sysname = systematics.sys_groups[year][sys][0]
            dw_sysname = systematics.sys_groups[year][sys][1] if len(systematics.sys_groups[year][sys]) > 1 else None
            procs_sys = sorted(set([key.split(f"_{up_sysname}")[0] for key in hdict[jmult].keys() if up_sysname in key])) if not dw_sysname \
                else sorted(set([key.split(f"_{up_sysname}")[0] for key in hdict[jmult].keys() if up_sysname in key] + [key.split(f"_{dw_sysname}")[0] for key in hdict[jmult].keys() if dw_sysname in key]))

            if not procs_sys: continue

            #set_trace()
            for proc in procs_sys:
                    # skip this process/systematic if any of the variations or nominal name aren't found
                if not ( (f"{proc}_{up_sysname}" in hdict[jmult].keys() and f"{proc}_{dw_sysname}" in hdict[jmult].keys()) and f"{proc}_nosys" in hdict[jmult].keys() ): continue
                print(year, jmult, sys, proc)

                if not dw_sysname: set_trace()
                    # only symmetrize using down variation
                if sys == "SHAPE":
                    #set_trace()
                    symmetrized_histo_up = hdict[jmult][f"{proc}_{dw_sysname}"].copy()
                    symmetrized_histo_dw = hdict[jmult][f"{proc}_{dw_sysname}"].copy()
                    dw_ratio_vals = symmetrized_histo_dw.values()[()]
                        # find relative deviation from nominal vals
                    dw_rel_vals = dw_ratio_vals - 1.
                    up_rel_vals = -1 * dw_rel_vals

                    symmetrized_ratios_up = up_rel_vals + 1.
                    symmetrized_ratios_dw = dw_rel_vals + 1.
                    symmetrized_histo_up.values()[()][:] = symmetrized_ratios_up
                    symmetrized_histo_dw.values()[()][:] = symmetrized_ratios_dw

                else:
                        # get ratios
                    symmetrized_histo_up = hdict[jmult][f"{proc}_{up_sysname}"].copy()
                    symmetrized_histo_dw = hdict[jmult][f"{proc}_{dw_sysname}"].copy()
                    up_ratio_vals = symmetrized_histo_up.values()[()]
                    dw_ratio_vals = symmetrized_histo_dw.values()[()]

                        # find relative deviation from nominal vals
                    up_rel_vals = up_ratio_vals - 1.
                    dw_rel_vals = dw_ratio_vals - 1.

                    # create hist with symmetrized yields
                    relative_symm_yields = np.sqrt(up_ratio_vals/dw_ratio_vals) - 1. # find relative deviation from nominal values for symmetrized dist
                    symmetrized_ratios_up = np.nan_to_num((relative_symm_yields + 1.))
                    symmetrized_ratios_dw = np.nan_to_num((-1.*relative_symm_yields + 1.))

                        # check where product of relative bins is positive, i.e. the deviations are one-sided
                    rel_prods = up_rel_vals * dw_rel_vals >= 0
                        # symmetrize all bins
                    one_sided_bins = np.arange(symmetrized_ratios_up.size) if (np.sum(rel_prods, axis=0) > 10) else np.where(rel_prods)
                    #if (one_sided_bins)[0].size > 0: set_trace()
                    symmetrized_histo_up.values()[()][one_sided_bins] = symmetrized_ratios_up[one_sided_bins]
                    symmetrized_histo_dw.values()[()][one_sided_bins] = symmetrized_ratios_dw[one_sided_bins]
                    if (one_sided_bins)[0].size > 0: print("\tsymmetrize")

                ## save templates that need to be symmetrized to dict
                histo_dict[jmult][f"{proc}_{up_sysname}"] = symmetrized_histo_up.copy()
                histo_dict[jmult][f"{proc}_{dw_sysname}"] = symmetrized_histo_dw.copy()

    #set_trace()
    coffea_out = os.path.join(input_dir, f"symmetrized_{args.sys_treatment}_combined_lep_templates_lj_{process}_{year}_{jobid}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")


def symm_indiv_templates(fname, process):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """
    hdict = load(fname)
    histo_dict = processor.dict_accumulator({njets : {"Muon" : {}, "Electron" :{}} for njets in hdict.keys()})
    for jmult in hdict.keys():
        for lep, histo in hdict[jmult].items():
            for sys in systypes:
                up_sysname = systematics.sys_groups[year][sys][0]
                dw_sysname = systematics.sys_groups[year][sys][1] if len(systematics.sys_groups[year][sys]) > 1 else None
                procs_sys = sorted(set([key.split(f"_{up_sysname}")[0] for key in histo.keys() if up_sysname in key])) if not dw_sysname \
                    else sorted(set([key.split(f"_{up_sysname}")[0] for key in histo.keys() if up_sysname in key] + [key.split(f"_{dw_sysname}")[0] for key in histo.keys() if dw_sysname in key]))

                if not procs_sys: continue

                #set_trace()
                for proc in procs_sys:
                        # skip this process/systematic if any of the variations or nominal name aren't found
                    #if not dw_sysname: set_trace()
                    if not ( (f"{proc}_{up_sysname}" in histo.keys() and f"{proc}_{dw_sysname}" in histo.keys()) and f"{proc}_nosys" in histo.keys() ): continue
                    print(year, jmult, lep, sys, proc)

                        # get yield vals
                    nosys_vals = histo[f"{proc}_nosys"].values()[()]
                    symmetrized_histo_up = histo[f"{proc}_{up_sysname}"].copy()
                    symmetrized_histo_dw = histo[f"{proc}_{dw_sysname}"].copy()
                    up_yield_vals = symmetrized_histo_up.values()[()]
                    dw_yield_vals = symmetrized_histo_dw.values()[()]

                        # find relative deviation from nominal vals
                    up_rel_vals = up_yield_vals/nosys_vals - 1.
                    dw_rel_vals = dw_yield_vals/nosys_vals - 1.

                    # create hist with symmetrized yields
                    relative_symm_yields = np.sqrt(up_yield_vals/dw_yield_vals) - 1. # find relative deviation from nominal values for symmetrized dist
                    symmetrized_yields_up = np.nan_to_num((relative_symm_yields + 1.) * nosys_vals)
                    symmetrized_yields_dw = np.nan_to_num((-1.*relative_symm_yields + 1.) * nosys_vals)

                        # check where product of relative bins is positive, i.e. the deviations are one-sided
                    rel_prods = up_rel_vals * dw_rel_vals >= 0
                        # symmetrize all bins
                    one_sided_bins = np.arange(symmetrized_yields_up.size) if (np.sum(rel_prods, axis=0) > 10) else np.where(rel_prods)
                    #if (one_sided_bins)[0].size > 0: set_trace()
                    symmetrized_histo_up.values()[()][one_sided_bins] = symmetrized_yields_up[one_sided_bins]
                    symmetrized_histo_dw.values()[()][one_sided_bins] = symmetrized_yields_dw[one_sided_bins]
                    if (one_sided_bins)[0].size > 0: print("\tsymmetrize")

                    ## save templates that need to be symmetrized to dict
                    histo_dict[jmult][lep][f"{proc}_{up_sysname}"] = symmetrized_histo_up.copy()
                    histo_dict[jmult][lep][f"{proc}_{dw_sysname}"] = symmetrized_histo_dw.copy()

    #set_trace()
    coffea_out = os.path.join(input_dir, f"symmetrized_{args.sys_treatment}_templates_lj_{process}_{year}_{jobid}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")



def symm_1D_year_and_lepton_templates(fname, process):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """
    hdict = load(fname)
    histo_dict = processor.dict_accumulator({njets : {} for njets in hdict.keys()})
    #set_trace()
    for jmult in hdict.keys():
        for sys in systypes:
            up_sysname = systematics.sys_groups["2017"][sys][0]
            dw_sysname = systematics.sys_groups["2017"][sys][1] if len(systematics.sys_groups["2017"][sys]) > 1 else None
            procs_sys = sorted(set([key.split(f"_{up_sysname}")[0] for key in hdict[jmult].keys() if up_sysname in key])) if not dw_sysname \
                else sorted(set([key.split(f"_{up_sysname}")[0] for key in hdict[jmult].keys() if up_sysname in key] + [key.split(f"_{dw_sysname}")[0] for key in hdict[jmult].keys() if dw_sysname in key]))

            if not procs_sys: continue

            #set_trace()
            for proc in procs_sys:
                    # skip this process/systematic if any of the variations or nominal name aren't found
                if (up_sysname is not None) & (dw_sysname is not None):
                    exists = ( (f"{proc}_{up_sysname}" in hdict[jmult].keys() and f"{proc}_{dw_sysname}" in hdict[jmult].keys()) and f"{proc}_nosys" in hdict[jmult].keys() )
                    if not exists: continue
                    print(jmult, sys, proc)

                        # get ratios
                    symmetrized_histo_up = hdict[jmult][f"{proc}_{up_sysname}"].copy()
                    symmetrized_histo_dw = hdict[jmult][f"{proc}_{dw_sysname}"].copy()
                    up_ratio_vals = symmetrized_histo_up.values()[()]
                    dw_ratio_vals = symmetrized_histo_dw.values()[()]

                        # find relative deviation from nominal vals
                    up_rel_vals = up_ratio_vals - 1.
                    dw_rel_vals = dw_ratio_vals - 1.

                    # create hist with symmetrized yields
                    relative_symm_yields = np.sqrt(up_ratio_vals/dw_ratio_vals) - 1. # find relative deviation from nominal values for symmetrized dist
                    symmetrized_ratios_up = np.nan_to_num((relative_symm_yields + 1.))
                    symmetrized_ratios_dw = np.nan_to_num((-1.*relative_symm_yields + 1.))

                        # check where product of relative bins is positive, i.e. the deviations are one-sided
                    rel_prods = up_rel_vals * dw_rel_vals >= 0
                        # symmetrize all bins
                    one_sided_bins = np.arange(symmetrized_ratios_up.size) if (np.sum(rel_prods, axis=0) > 10) else np.where(rel_prods)
                    #if (one_sided_bins)[0].size > 0: set_trace()
                    symmetrized_histo_up.values()[()][one_sided_bins] = symmetrized_ratios_up[one_sided_bins]
                    symmetrized_histo_dw.values()[()][one_sided_bins] = symmetrized_ratios_dw[one_sided_bins]
                    if (one_sided_bins)[0].size > 0: print("\tsymmetrize")

                if (up_sysname is not None) & (dw_sysname is None):
                    exists = f"{proc}_{up_sysname}" in hdict[jmult].keys() and f"{proc}_nosys" in hdict[jmult].keys()
                    if not exists: continue
                    print(jmult, sys, proc)

                        # get ratios
                    symmetrized_histo_up = hdict[jmult][f"{proc}_{up_sysname}"].copy()
                    symmetrized_histo_dw = symmetrized_histo_up.copy()
                    up_ratio_vals = symmetrized_histo_up.values()[()]

                        # find relative deviation from nominal vals
                    up_rel_vals = up_ratio_vals - 1.
                    dw_rel_vals = -1. * up_rel_vals

                    symmetrized_ratios_up = up_rel_vals + 1.
                    symmetrized_ratios_dw = dw_rel_vals + 1.
                    symmetrized_histo_up.values()[()][:] = symmetrized_ratios_up
                    symmetrized_histo_dw.values()[()][:] = symmetrized_ratios_dw
                    print("\tsymmetrize")
                    #set_trace()
                    dw_sysname = f"{up_sysname}Down"
                    up_sysname = f"{up_sysname}Up"

                ## save templates that need to be symmetrized to dict
                histo_dict[jmult][f"{proc}_{up_sysname}"] = symmetrized_histo_up.copy()
                histo_dict[jmult][f"{proc}_{dw_sysname}"] = symmetrized_histo_dw.copy()

    #set_trace()
    coffea_out = os.path.join(input_dir, f"symmetrized_{args.sys_treatment}_1D_combined_year_and_lepton_templates_lj_{process}_{jobid}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")



if __name__ == "__main__":
    allowed_combination_options = ["lepton", "era_lepton", "indiv"]
    combinations_to_run = (args.channels_to_combine).split(":")
    combinations_to_run = [combination for combination in combinations_to_run if combination in allowed_combination_options]

    jobid = os.environ["jobid"]
    eos_dir = os.environ["eos_dir"]

    analyzer = "htt_btag_sb_regions"

    for combination in combinations_to_run:
        if combination == "era_lepton":
            systypes = ["EWK_scheme", "EWK_yukawa", "ISR", "FSR", "FACTOR", "RENORM", "HDAMP", "UE", "MTOP3GEV", "MTOP1GEV", "CR1", "CR2", "erdON", "dQCD", "ST_ISR", "ST_FSR", "ST_RENORM", "ST_FACTOR"]

            #set_trace()
            input_dir = os.path.join(eos_dir, "results", jobid, f"Templates_{analyzer}")
            if not os.path.isdir(input_dir): raise ValueError("No directory found.")

            bkg_fname = os.path.join(input_dir, f"{args.sys_treatment}_combined_year_and_lepton_templates_lj_bkg_nomSMTTxsec_{jobid}.coffea" if args.nomSMTTxsec else f"{args.sys_treatment}_combined_year_and_lepton_templates_lj_bkg_{jobid}.coffea")
            if not os.path.isfile(bkg_fname): raise ValueError(f"{bkg_fname} not found.")
            print("Symmetrizing combined channels across years and leptons for background templates")
            symm_year_and_lepton_templates(bkg_fname, "bkg")

            if args.sys_treatment == "smoothed":
                bkg_1d_fname = os.path.join(input_dir, f"{args.sys_treatment}_1D_combined_year_and_lepton_templates_lj_bkg_{jobid}.coffea")
                print("Symmetrizing combined channels across years and leptons for background templates")
                symm_1D_year_and_lepton_templates(bkg_1d_fname, "bkg")


        if combination == "lepton":
            for year in ["2016APV", "2016", "2017", "2018"]:
                input_dir = os.path.join(eos_dir, "results", f"{year}_{jobid}", f"Templates_{analyzer}")
                if not os.path.isdir(input_dir): raise ValueError("No directory found.")

                systypes = [
                    "BTAG_BC_JES", "BTAG_BC_PILEUP", "BTAG_BC_STATISTIC", "BTAG_BC_TYPE3",
                    "BTAG_BC_CORR", "BTAG_BC_UNCORR", "BTAG_L_CORR", "BTAG_L_UNCORR",
                    "BTAG_BC_BFRAGMENTATION", "BTAG_BC_BTEMPCORR", "BTAG_BC_CB", "BTAG_BC_CFRAGMENTATION", "BTAG_BC_CJETS",
                    "BTAG_BC_DMUX", "BTAG_BC_GLUONSPLITTING", "BTAG_BC_JETAWAY", "BTAG_BC_KSL", "BTAG_BC_L2C",
                    "BTAG_BC_LTOTHERS", "BTAG_BC_MUDR", "BTAG_BC_MUPT", "BTAG_BC_PTREL",
                    "JES_Absolute", f"JES_Absolute_{year}", "JES_BBEC1", f"JES_BBEC1_{year}",
                    "JES_FlavorQCD", "JES_FlavorQCDOnlyLightJets",
                    "JES_FlavorPureBottom", "JES_FlavorPureCharm", "JES_FlavorPureGluon", "JES_FlavorPureQuark",
                    "JES_FlavorPureBottomOnlyBottomJets", "JES_FlavorPureCharmOnlyCharmJets", "JES_FlavorPureGluonOnlyGluonJets", "JES_FlavorPureQuarkOnlyQuarkJets",
                    "JES_RelativeBal", f"JES_RelativeSample_{year}", "MET", "JER", "SHAPE",
                ]                
                if year != "2018": systypes.append("PREFIRE")

                bkg_fname = os.path.join(input_dir, f"{args.sys_treatment}_combined_lep_templates_lj_bkg_{year}_{jobid}.coffea")
                if not os.path.isfile(bkg_fname): raise ValueError(f"{bkg_fname} not found.")
                print(f"Symmetrizing e+mu channels in {year} for background templates")
                symm_lepton_templates(bkg_fname, "bkg")

        if combination == "indiv":
            for year in ["2016APV", "2016", "2017", "2018"]:
                input_dir = os.path.join(eos_dir, "results", f"{year}_{jobid}", f"Templates_{analyzer}")
                if not os.path.isdir(input_dir): raise ValueError("No directory found.")
                
                systypes = systematics.sys_groups[year].keys()

                bkg_fname = os.path.join(input_dir, f"{args.sys_treatment}_templates_lj_bkg_{year}_{jobid}.coffea")
                if not os.path.isfile(bkg_fname): raise ValueError(f"{bkg_fname} not found.")
                print(f"Symmetrizing all channels in {year} for background templates")
                symm_indiv_templates(bkg_fname, "bkg")

    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
