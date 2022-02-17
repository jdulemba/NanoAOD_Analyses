#! /bin/env python

import time
tic = time.time()

from coffea.util import load, save
from pdb import set_trace
import os
import Utilities.plot_tools as plt_tools
from coffea import hist
import numpy as np
import fnmatch
import Utilities.Plotter as Plotter
import coffea.processor as processor    
import Utilities.final_analysis_binning as final_binning

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
args = parser.parse_args()


def get_bkg_templates():
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """
    ## variables that only need to be defined/evaluated once
    hdict = plt_tools.add_coffea_files(bkg_fnames) if len(bkg_fnames) > 1 else load(bkg_fnames[0])
    histo_dict = processor.dict_accumulator({njets : {"Muon" : {}, "Electron" :{}} for njets in njets_to_run})

        # get correct hist and rebin
    hname_to_use = "mtt_vs_tlep_ctstar_abs"
    if hname_to_use not in hdict.keys():
        raise ValueError(f"{hname_to_use} not found in file")
    xrebinning, yrebinning = linearize_binning
    histo = hdict[hname_to_use][Plotter.nonsignal_samples] # process, sys, jmult, leptype, btag, lepcat
    
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
    
        ## scale ttJets events, split by reconstruction type, by normal ttJets lumi correction
    ttJets_permcats = ["*right", "*matchable", "*unmatchable", "*sl_tau", "*other"]
    names = [dataset for dataset in sorted(set([key[0] for key in histo.values().keys()]))] # get dataset names in hists
    ttJets_cats = [name for name in names if any([fnmatch.fnmatch(name, cat) for cat in ttJets_permcats])] # gets ttJets(_PS)_other, ...

        ## make groups based on process
    process = hist.Cat("process", "Process", sorting="placement")
    process_cat = "dataset"

    for lep in ["Muon", "Electron"]:
        #set_trace()    
        ## make groups based on process
        process_groups = plt_tools.make_dataset_groups(lep, args.year, samples=names, gdict="templates")
        
        lumi_correction = lumi_corr_dict[args.year][f"{lep}s"]
                # scale ttJets events, split by reconstruction type, by normal ttJets lumi correction
        if len(ttJets_cats) > 0:
            for tt_cat in ttJets_cats:
                ttJets_lumi_topo = "_".join(tt_cat.split("_")[:-2]) if "sl_tau" in tt_cat else "_".join(tt_cat.split("_")[:-1]) # gets ttJets[SL, Had, DiLep] or ttJets_PS
                ttJets_eff_lumi = lumi_correction[ttJets_lumi_topo]
                lumi_correction.update({tt_cat: ttJets_eff_lumi})
    
        histo = rebin_histo.copy()
        histo.scale(lumi_correction, axis="dataset")
        histo = histo.group(process_cat, process, process_groups)[:, :, :, lep, :, :].integrate("leptype")

            # remove 0th and last 2 pdf replicas 0th is base set compatible with 1, last two sets are variations in alpha_S
        histo = histo.remove(["pdf_0"], "sys")
        systs = sorted(set([key[1] for key in histo.values().keys()]))
        systs.insert(0, systs.pop(systs.index("nosys"))) # move "nosys" to the front

            # loop over each jet multiplicity
        for jmult in njets_to_run:
            sig_histo = Plotter.linearize_hist(histo[:, :, jmult].integrate("jmult"))

                # loop over each systematic
            for sys in systs:
                sys_histo = sig_histo[:, sys].integrate("sys")

                    ## write nominal and systematic variations for each topology to file
                for proc in sorted(set([key[0] for key in sys_histo.values().keys()])):
                    if not sys_histo[proc].values().keys():
                        print(f"Systematic {sys} for {lep} {jmult} {proc} not found, skipping")
                        continue
                    if sys == "pdf_101": sys = "alphaSDown"
                    if sys == "pdf_102" : sys = "alphaSUp"

                    print(args.year, lep, jmult, sys, proc)
                    template_histo = sys_histo[proc].integrate("process")

                        ## save template histos to coffea dict
                    histo_dict[jmult][lep][f"{proc}_{sys}"] = template_histo.copy()

    coffea_out = os.path.join(outdir, f"raw_pdf_templates_lj_bkg_{args.year}_{jobid}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")


if __name__ == "__main__":
    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]
    base_jobid = os.environ["base_jobid"]

    njets_to_run = ["3Jets", "4PJets"]
    
        ## initialize lumi scaling files 
    lumi_corr_dict = load(os.path.join(proj_dir, "Corrections", base_jobid, "MC_LumiWeights.coffea"))

        # define variables to get histogram for background    
    bkg_analyzer = "htt_pdfUncs"
    bkg_input_dir = os.path.join(proj_dir, "results", f"{args.year}_{jobid}", bkg_analyzer)
    if os.path.isdir(bkg_input_dir):
        bkg_fnames = fnmatch.filter(os.listdir(bkg_input_dir), "*TOT.coffea")
        bkg_fnames = [os.path.join(bkg_input_dir, fname) for fname in bkg_fnames]
    else: print("No background file found.")

    outdir = os.path.join(proj_dir, "results", f"{args.year}_{jobid}", f"Templates_{bkg_analyzer}")
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
  
    linearize_binning = (
        final_binning.mtt_binning,
        final_binning.ctstar_abs_binning
    ) 

    print("Creating background templates")
    get_bkg_templates()

    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
