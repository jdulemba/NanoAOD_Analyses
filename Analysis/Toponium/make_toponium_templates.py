#! /bin/env python
import time
tic = time.time()

from coffea.util import load, save
from pdb import set_trace
import os
from coffea import hist
import numpy as np
import fnmatch
import Utilities.Plotter as Plotter
import coffea.processor as processor    
import Utilities.final_analysis_binning as final_binning
from Utilities import systematics

import argparse
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("template_type", choices=["Raw", "Combined_Era_Lepton", "Combined_Lepton"], help="What type of template do you want to make?")
parser.add_argument("--additional_treatment", choices=["smooth", "flatten", "smooth_symmetrize", "flatten_symmetrize"], help="Specify which additional treatment to apply to the template types.")
parser.add_argument("--no_root", action="store_true", help="Suppress writing toponium templates to a root file (default is to write).")
args = parser.parse_args()

version = "V2"
#version = "V1"

def make_output_dir(analyzer):
    outdir = os.path.join(eos_dir, "results", jobid, analyzer, "Templates", version)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    return outdir
  

def make_raw_templates():
    tot_hdict = {year : {} for year in years_to_run}
    
    if args.additional_treatment == "smooth": smoothed_output = processor.dict_accumulator({year : {njets : {"Muon" : {}, "Electron" :{}} for njets in njets_to_run} for year in years_to_run})
    if args.additional_treatment == "flatten": flattened_output = processor.dict_accumulator({year : {njets : {"Muon" : {}, "Electron" :{}} for njets in njets_to_run} for year in years_to_run})

    # loop over each year to lumi scale systematics
    for year in years_to_run:
        input_dir = os.path.join(eos_dir, "results", jobid, "Toponium_HBSR", year)
        if os.path.isdir(input_dir):
            fnames = fnmatch.filter(os.listdir(input_dir), "*TOT.coffea")
            fnames = [os.path.join(input_dir, fname) for fname in fnames]
            if len(fnames) > 1: raise ValueError("Multiple input files found")
        else: raise ValueError("No input file found.")

        toponium_fname = fnames[0]
        tot_hdict[year] = make_year_raw_templates(toponium_fname, year)

            # smooth templates
        if args.additional_treatment == "smooth":
            #set_trace()
            for jmult in tot_hdict[year].keys():
                for lep in tot_hdict[year][jmult].keys():
                    for hname in tot_hdict[year][jmult][lep].keys():
                        smoothed_output[year][jmult][lep][hname] = smooth_templates(nosys=tot_hdict[year][jmult][lep]["Toponium_nosys"].copy(), variation=tot_hdict[year][jmult][lep][hname].copy())

            # flatten templates
        if args.additional_treatment == "flatten":
            #set_trace()
            for jmult in tot_hdict[year].keys():
                for lep in tot_hdict[year][jmult].keys():
                    for hname in tot_hdict[year][jmult][lep].keys():
                        flattened_output[year][jmult][lep][hname] = Plotter.flatten(nosys=tot_hdict[year][jmult][lep]["Toponium_nosys"].copy(), systematic=tot_hdict[year][jmult][lep][hname].copy())
    
    coffea_out = os.path.join(outdir, f"raw_templates_lj_toponium_{jobid}_{version}.coffea")
    save(tot_hdict, coffea_out)
    print(f"{coffea_out} written")

        ## save smooth templates
    if args.additional_treatment == "smooth":
        smooth_outname = os.path.join(outdir, f"smoothed_indiv_templates_lj_toponium_{jobid}_{version}.coffea")
        save(smoothed_output, smooth_outname)
        print(f"{smooth_outname} written")

        ## save flattened templates
    if args.additional_treatment == "flatten":
        flat_outname = os.path.join(outdir, f"flattened_indiv_templates_lj_toponium_{jobid}_{version}.coffea")
        save(flattened_output, flat_outname)
        print(f"{flat_outname} written")


def make_year_raw_templates(input_fname, year):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """
    #set_trace()
    hdict = load(input_fname)

        # get correct hist and rebin
    hname_to_use = "mtt_vs_tlep_ctstar_abs"
    if hname_to_use not in hdict.keys():
        raise ValueError(f"{hname_to_use} not found in file")
    xrebinning, yrebinning = linearize_binning
    histo = hdict[hname_to_use][Plotter.nonsignal_samples][:, :, :, :, "btagPass"].integrate("btag") # process, sys, jmult, leptype, btag, lepcat
    
    xaxis_name = histo.dense_axes()[0].name
    yaxis_name = histo.dense_axes()[1].name
        ## rebin x axis
    if isinstance(xrebinning, np.ndarray):
        new_xbins = hist.Bin(xaxis_name, xaxis_name, xrebinning)
    elif isinstance(xrebinning, float) or isinstance(xrebinning, int):
        new_xbins = xrebinning
    histo = histo.rebin(xaxis_name, new_xbins)
        ## rebin y axis
    if isinstance(yrebinning, np.ndarray):
        new_ybins = hist.Bin(yaxis_name, yaxis_name, yrebinning)
    elif isinstance(yrebinning, float) or isinstance(yrebinning, int):
        new_ybins = yrebinning
    rebin_histo = histo.rebin(yaxis_name, new_ybins)

        ## make groups based on process
    process = hist.Cat("process", "Process", sorting="placement")
    process_cat = "dataset"
    process_groups = {"Toponium" : ["ToponiumSL", "ToponiumDiLep"]}

    histo_dict = processor.dict_accumulator({njets : {"Muon" : {}, "Electron" :{}} for njets in njets_to_run})

    for lep in ["Muon", "Electron"]:
        histo = rebin_histo.copy()

        ## make groups and normalize based on processes
            # scale toponium events by mass window cut for lumi correction
        lumi_correction = lumi_corr_dict[year][f"{lep}s"]
        lumi_correction.update({"_".join(key) : lumi_correction[f"{key[0]}_nosys"] for key in sorted(set([(key[0], key[1]) for key in histo.values().keys()])) if key[1] not in toponium_lumi_systs})
            # scale {(ToponiumSL, nosys) : ToponiumSL_nosys}
        histo.scale({key : lumi_correction["_".join(key)] for key in sorted(set([(key[0], key[1]) for key in histo.values().keys()]))}, axis=("dataset", "sys"))
        histo = histo.group(process_cat, process, process_groups)[:, :, :, lep, :, :].integrate("leptype")

        systs = sorted(set([key[1] for key in histo.values().keys()]))

            # loop over each jet multiplicity
        for jmult in njets_to_run:
            linearized_histo = Plotter.linearize_hist(histo[:, :, jmult].integrate("jmult"))
                # loop over each systematic
            for sys in systs:
                sys_histo = (linearized_histo[:, sys].integrate("sys")).copy()

                    ## write nominal and systematic variations for each topology to file
                for proc in sorted(set([key[0] for key in sys_histo.values().keys()])):
                    if proc != "Toponium": raise ValueError(f"Only valid proc is Toponium! Not {proc}")
                    if not sys_histo[proc].values().keys():
                        print(f"Systematic {sys} for {lep} {jmult} {proc} not found, skipping")
                        continue

                    print(year, lep, jmult, sys, proc)
                    template_histo = (sys_histo[proc].integrate("process")).copy()
                    sumw, sumw2 = template_histo.values(sumw2=True)[()]
                    rel_err = np.sqrt(sumw2)/np.abs(sumw)
                    rel_err_mask = rel_err > 10
                    if np.any(rel_err_mask):
                        set_trace()

                        ## save template histos to coffea dict
                    histo_dict[jmult][lep][f"{proc}_nosys" if sys == "nosys" else f"{proc}_{systematics.sys_to_name[year][sys]}"] = template_histo.copy()

    return histo_dict


def combine_lepton_templates(input_dict):
    #set_trace()
        ## combine templates across lepton flavour channels
    combined_lep_acc = processor.dict_accumulator({year : {njets : {} for njets in njets_to_run} for year in input_dict.keys()})
    if args.additional_treatment == "smooth": smoothed_output = processor.dict_accumulator({year : {njets : {} for njets in njets_to_run} for year in input_dict.keys()})
    if args.additional_treatment == "flatten": flattened_output = processor.dict_accumulator({year : {njets : {} for njets in njets_to_run} for year in input_dict.keys()})

    for year in input_dict.keys():
        for jmult in njets_to_run:
                # find common names across lepton channels
            mu_tnames = sorted(input_dict[year][jmult]["Muon"].keys())
            el_tnames = sorted(input_dict[year][jmult]["Electron"].keys())
            tnames = sorted(set(mu_tnames) & set(el_tnames))
            for tname in tnames:
                print(year, jmult, tname)

                combined_lep_template = input_dict[year][jmult]["Muon"][tname].copy()
                combined_lep_template = combined_lep_template.add(input_dict[year][jmult]["Electron"][tname].copy())

                    ## save template histos to coffea dict
                combined_lep_acc[year][jmult][tname] = combined_lep_template.copy()

                # smooth templates
            if args.additional_treatment == "smooth":
                for tname in tnames:
                    smoothed_output[year][jmult][tname] = smooth_templates(nosys=combined_lep_acc[year][jmult]["Toponium_nosys"].copy(), variation=combined_lep_acc[year][jmult][tname].copy())

                # flatten templates
            if args.additional_treatment == "flatten":
                for tname in tnames:
                    flattened_output[year][jmult][tname] = Plotter.flatten(nosys=combined_lep_acc[year][jmult]["Toponium_nosys"].copy(), systematic=combined_lep_acc[year][jmult][tname].copy())

    coffea_out = os.path.join(outdir, f"raw_combined_lepton_templates_lj_toponium_{jobid}_{version}.coffea")
    save(combined_lep_acc, coffea_out)
    print(f"{coffea_out} written")

        ## save smooth templates
    if args.additional_treatment == "smooth":
        smooth_outname = os.path.join(outdir, f"smoothed_combined_lepton_templates_lj_toponium_{jobid}_{version}.coffea")
        save(smoothed_output, smooth_outname)
        print(f"{smooth_outname} written")

        ## save flattened templates
    if args.additional_treatment == "flatten":
        flat_outname = os.path.join(outdir, f"flattened_combined_lepton_templates_lj_toponium_{jobid}_{version}.coffea")
        save(flattened_output, flat_outname)
        print(f"{flat_outname} written")



def combine_era_lepton_templates(input_dict):
        ## combine templates across years and lepton flavour channels
    combined_acc = processor.dict_accumulator({njets : {} for njets in njets_to_run})
    if args.additional_treatment == "smooth": smoothed_output = processor.dict_accumulator({njets : {} for njets in njets_to_run})
    if args.additional_treatment == "flatten": flattened_output = processor.dict_accumulator({njets : {} for njets in njets_to_run})

    for jmult in njets_to_run:
            # find common names across lepton channels
        mu_tnames_dict = {year : sorted(input_dict[year][jmult]["Muon"].keys()) for year in input_dict.keys() }
        el_tnames_dict = {year : sorted(input_dict[year][jmult]["Electron"].keys()) for year in input_dict.keys() }
        tnames = sorted( set.intersection(*map(set, sorted(mu_tnames_dict.values()))) & set.intersection(*map(set, sorted(el_tnames_dict.values()))) )
        for tname in tnames:
            print(jmult, tname)

            combined_year_template = None
            idx = 0
            for year in input_dict.keys():
                for lep in sorted(["Muon", "Electron"]):
                    combined_year_template = input_dict[year][jmult][lep][tname].copy() if idx == 0 else combined_year_template.add(input_dict[year][jmult][lep][tname].copy())
                    idx += 1

                ## save template histos to coffea dict
            combined_acc[jmult][tname] = combined_year_template.copy()

            # smooth templates
        if args.additional_treatment == "smooth":
            for tname in tnames:
                smoothed_output[jmult][tname] = smooth_templates(nosys=combined_acc[jmult]["Toponium_nosys"].copy(), variation=combined_acc[jmult][tname].copy())

            # flatten templates
        if args.additional_treatment == "flatten":
            for tname in tnames:
                flattened_output[jmult][tname] = Plotter.flatten(nosys=combined_acc[jmult]["Toponium_nosys"].copy(), systematic=combined_acc[jmult][tname].copy())

    coffea_out = os.path.join(outdir, f"raw_combined_era_lepton_templates_lj_toponium_{jobid}_{version}.coffea")
    save(combined_acc, coffea_out)
    print(f"{coffea_out} written")

        ## save smooth templates
    if args.additional_treatment == "smooth":
        smooth_outname = os.path.join(outdir, f"smoothed_combined_era_lepton_templates_lj_toponium_{jobid}_{version}.coffea")
        save(smoothed_output, smooth_outname)
        print(f"{smooth_outname} written")

        ## save flattened templates
    if args.additional_treatment == "flatten":
        flat_outname = os.path.join(outdir, f"flattened_combined_era_lepton_templates_lj_toponium_{jobid}_{version}.coffea")
        save(flattened_output, flat_outname)
        print(f"{flat_outname} written")



def smooth_templates(nosys, variation, p_end=0.1):
    """
    Smooth templates
    """
        # perform smoothing
    smoothed_histo = Plotter.new_smoothing(nosys=nosys.copy(), systematic=variation.copy(), mtt_bins=smooth_binning[0],
        nbinsx=len(smooth_binning[0])-1, nbinsy=len(smooth_binning[1])-1, p_end=p_end)

    return smoothed_histo



def write_rootfile(dist_type):
    import uproot
        ## save templates to root file for each distribution for each year
    if dist_type == "Raw":
        input_tfile = os.path.join(outdir, f"raw_templates_lj_toponium_{jobid}_{version}.coffea")
        if not os.path.isfile(input_tfile): raise ValueError(f"{input_tfile} not found")
        tdict = load(input_tfile)
    else:
        raise ValueError("Only 'Raw' allowed as distribution type at this point")

    output_rname = input_tfile.replace(".coffea", ".root")
    rfile = uproot.recreate(output_rname, compression=uproot.ZLIB(4)) if os.path.isfile(output_rname) else uproot.create(output_rname)

    for year in tdict.keys():
        if year == "2016APV": year_to_use = "2016pre"
        elif year == "2016": year_to_use = "2016post"
        else: year_to_use = year

        for jmult in tdict[year].keys():
            for lep in tdict[year][jmult].keys():
                orig_lepdir = "muNJETS" if lep == "Muon" else "eNJETS"
                lepdir = orig_lepdir.replace("NJETS", jmult.lower())

                dirname = f"{lepdir}_{year_to_use}"
                rfile.mkdir(dirname)

                for hname, template in tdict[year][jmult][lep].items():
                    if "nosys" in hname:
                        rfile[dirname]["EtaT"] = (template.copy()).to_hist()
                    elif "BindingEnergy" in hname:
                        rfile[dirname][hname.replace("Toponium", "EtaT").replace("BindingEnergy", "Eb")] = (template.copy()).to_hist()
                    elif "TopMass" in hname:
                        rfile[dirname][hname.replace("Toponium", "EtaT").replace("TopMass", "tmass_EtaT")] = (template.copy()).to_hist()
                    else:
                        continue
                        print(hname)
                        #set_trace()
                        rfile[dirname][f"EtaT_{systematics.combine_template_sys_to_name[year][hname.replace('Toponium_', '')]}"] = (template.copy()).to_hist()

    rfile.close()
    print(f"{output_rname} written")




if __name__ == "__main__":
    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]
    base_jobid = os.environ["base_jobid"]
    eos_dir = os.environ["eos_dir"]

    njets_to_run = sorted(["3Jets", "4PJets"])
    years_to_run = ["2016APV", "2016", "2017", "2018"]

        ## initialize lumi scaling files
    lumi_name = "MC_LumiWeights"
    lumi_corr_dict = load(os.path.join(proj_dir, "Corrections", base_jobid, f"{lumi_name}.coffea"))

    linearize_binning = (
        final_binning.mtt_binning,
        final_binning.ctstar_abs_binning
    ) 

    if args.additional_treatment == "smooth":
            # replace last mtt bin with 1700 instead of 2000
        smooth_binning = (np.copy(linearize_binning[0]), np.copy(linearize_binning[1]))
        smooth_binning[0][-2] = smooth_binning[0][-3] + 100.
        smooth_binning[0][-1] = smooth_binning[0][-2] + 100.

        mtt_centers =  np.array([(smooth_binning[0][idx]+smooth_binning[0][idx+1])/2 for idx in range(len(smooth_binning[0])-1)])

    toponium_lumi_systs = ["nosys", "BindingEnergyUp", "BindingEnergyDown", "TopMassUp", "TopMassDown"]

    outdir = make_output_dir("Toponium_HBSR")

    if args.template_type == "Raw":
        make_raw_templates()
        if not args.no_root: write_rootfile("Raw")

    if args.template_type == "Combined_Lepton":
            # make sure raw template file exists
        raw_tfile = os.path.join(outdir, f"raw_templates_lj_toponium_{jobid}_{version}.coffea")
        if not os.path.isfile(raw_tfile):
            make_raw_templates()
        combine_lepton_templates(load(raw_tfile))

    if args.template_type == "Combined_Era_Lepton":
            # make sure raw template file exists
        raw_tfile = os.path.join(outdir, f"raw_templates_lj_toponium_{jobid}_{version}.coffea")
        if not os.path.isfile(raw_tfile):
            make_raw_templates()
        combine_era_lepton_templates(load(raw_tfile))




    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
