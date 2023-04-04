#! /bin/env python

import time
tic = time.time()

from coffea.util import load
from pdb import set_trace
import os
import Utilities.plot_tools as plt_tools
from coffea import hist
import numpy as np
import Utilities.Plotter as Plotter
import Utilities.final_analysis_binning as final_binning
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
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("--MEopts", nargs="*", action=ParseKwargs, help="Options to specify which mass and width points to produce for ME reweighted signal.")
args = parser.parse_args()


def get_MEsig_templates():
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """
    #set_trace()
    fname_dir = f"{eos_dir}/results/{args.year}_{jobid}/signal_ME_evtReweighting/RecoLevel"
    print("Check which input file is being used and output filename")

    #    # mtop variations
    #mtop_fnames = {
    #    "2016APV" : "BATCH_MEsig_m365TOm1000_relw0p5TO25p0_MTopUncs_17March2023_2016APV_Summer20UL_DeepJet.coffea",
    #    "2016"    : "BATCH_MEsig_m365TOm1000_relw0p5TO25p0_MTopUncs_17March2023_2016_Summer20UL_DeepJet.coffea",
    #    "2017"    : "BATCH_MEsig_m365TOm1000_relw0p5TO25p0_MTopUncs_17March2023_2017_Summer20UL_DeepJet.coffea",
    #    "2018"    : "BATCH_MEsig_m365TOm1000_relw0p5TO25p0_MTopUncs_16March2023_2018_Summer20UL_DeepJet.coffea"
    #}
    #    # btag subsources
    #btagSub_fnames = {
    #    "2016APV" : "BATCH_MEsig_m365TOm1000_relw0p5TO25p0_BtagType3Subsources_10February2023_2016APV_Summer20UL_DeepJet.coffea",
    #    "2016" : "BATCH_MEsig_m365TOm1000_relw0p5TO25p0_BtagType3Subsources_11February2023_2016_Summer20UL_DeepJet.coffea",
    #    "2017" : "BATCH_MEsig_m365TOm1000_relw0p5TO25p0_BtagType3Subsources_10February2023_2017_Summer20UL_DeepJet.coffea",
    #    "2018" : "BATCH_MEsig_m365TOm1000_relw0p5TO25p0_BtagType3Subsources_31January2023_2018_Summer20UL_DeepJet.coffea"
    #}
    #    
    #fnames_to_run = [btagSub_fnames[args.year]]
    ##fnames_to_run = [mtop_fnames[args.year]]

    #set_trace()
    ##fnames_to_run = [fname for fname in os.listdir(fname_dir) if (fname.startswith("BATCH") and fname.endswith("TOT.coffea"))]
    fnames_to_run = [fname for fname in os.listdir(fname_dir) if (fname.startswith("m") and fname.endswith("TOT.coffea"))]

    if not fnames_to_run: raise ValueError(f"No file found in {fname_dir}")
    if len(fnames_to_run) > 1: raise ValueError(f"{len(fnames_to_run)} found to run. There should only be 1")

    if masses_to_run == "All":
        rname = os.path.join(outdir, f"raw_templates_lj_MEsig_{args.year}_{jobid}.root" if (widths_to_run == "All") else f"raw_templates_lj_MEsig_{''.join(allowed_widths).lower()}_{args.year}_{jobid}.root")
    elif widths_to_run == "All":
        rname = os.path.join(outdir, f"raw_templates_lj_MEsig_{args.year}_{jobid}.root" if (masses_to_run == "All") else f"raw_templates_lj_MEsig_{''.join(allowed_masses).lower()}_{args.year}_{jobid}.root")
    else:
        rname = os.path.join(outdir, f"raw_templates_lj_MEsig_{''.join(allowed_masses).lower()}_{''.join(allowed_widths).lower()}_{args.year}_{jobid}.root")

    #set_trace()
    print(f"creating {rname}")
    upfout = uproot.recreate(rname, compression=uproot.ZLIB(4)) if os.path.isfile(rname) else uproot.create(rname)

    for fname in fnames_to_run:
        print(f"Loading {fname}")
        hdict = load(os.path.join(fname_dir, fname))
        print(f"\tLoaded {fname}")

            # get correct hist and rebin
        if hname_to_use not in hdict.keys():
            raise ValueError(f"{hname_to_use} not found in file")
        xrebinning, yrebinning = linearize_binning
        histo = hdict[hname_to_use] # process, sys, jmult, leptype, btag, lepcat

        #set_trace()    
        xaxis_name = histo.dense_axes()[0].name
        yaxis_name = histo.dense_axes()[1].name

            ## rebin axes
        if not (np.array_equal(histo.axis(xaxis_name).edges(), xrebinning) & np.array_equal(histo.axis(yaxis_name).edges(), yrebinning)):
            print("rebinning histograms")
                ## rebin x axis
            if isinstance(xrebinning, np.ndarray):
                new_xbins = hist.Bin(xaxis_name, xaxis_name, xrebinning)
            elif isinstance(xrebinning, float) or isinstance(xrebinning, int):
                new_xbins = xrebinning

                ## rebin y axis
            if isinstance(yrebinning, np.ndarray):
                new_ybins = hist.Bin(yaxis_name, yaxis_name, yrebinning)
            elif isinstance(yrebinning, float) or isinstance(yrebinning, int):
                new_ybins = yrebinning
            rebin_histo = histo.rebin(xaxis_name, new_xbins).rebin(yaxis_name, new_ybins)

        else:
            rebin_histo = histo.copy()

        names = sorted(set([id.name for id in rebin_histo.axis("dataset").identifiers()]))
        names = [name for name in names for mass in allowed_masses for width in allowed_widths if f"{mass}_{width}" in name] # filter names based on chosen masses and widths
        process_groups = plt_tools.make_dataset_groups("Muon", args.year, samples=names, bkgdict="templates", sigdict="MEreweight_combined")
        print("Regrouping hists")
        rebin_histo = rebin_histo.group("dataset", hist.Cat("process", "Process", sorting="placement"), process_groups)
        print("\tRegrouped hists")

        signals = sorted(set([key[0] for key in rebin_histo.values().keys()]))
        #set_trace()

        systs = sorted(set([key[1] for key in rebin_histo.values().keys()]))

            # write signal dists to temp file
        for lep in leps_to_run:
            histo = rebin_histo[:, :, :, lep].integrate("leptype")

            for jmult in njets_to_run:
                dirname = f"{jmult}/{lep}"
                upfout.mkdir(dirname)
                for signal in signals:
                        # only run over allowed masses and width points for now
                    for sys in systs:
                        template_histo = Plotter.linearize_hist(histo[signal, sys, jmult].integrate("jmult").integrate("process").integrate("sys"))
                        if not template_histo.values().keys():
                            print(f"No template found for {args.year}, {lep}, {jmult}, {signal}, {sys}")
                            #set_trace()
                            continue
                        print(args.year, lep, jmult, signal, sys)

                            ## save template histos to coffea dict
                        upfout[dirname][f"{signal}_{sys}"] = template_histo.to_hist()

    #set_trace()
    upfout.close()
    print(f"{rname} written")



def make_output_dir(analyzer):
    outdir = os.path.join(eos_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}")
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    return outdir
  

if __name__ == "__main__":

    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]
    base_jobid = os.environ["base_jobid"]
    eos_dir = os.environ["eos_dir"]

    leps_to_run = sorted(["Muon", "Electron"])
    njets_to_run = sorted(["3Jets", "4PJets"])

    linearize_binning = (
        final_binning.mtt_binning,
        final_binning.ctstar_abs_binning
    ) 

    #set_trace()
    analyzer = "htt_btag_sb_regions"
    outdir = make_output_dir(analyzer)

    ## variables that only need to be defined/evaluated once
    hname_to_use = "mtt_vs_tlep_ctstar_abs"

    possible_masses = [365, 380, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000]
    possible_widths = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 8.0, 10.0, 13.0, 15.0, 18.0, 21.0, 25.0]
    possible_masses = [str(mass) for mass in possible_masses]
    possible_widths = [str(width) for width in possible_widths]

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
    
    #set_trace()
    allowed_masses = ["M"+mass for mass in allowed_masses]
    allowed_widths = ["W"+width.replace(".", "p") for width in allowed_widths]
    print(f"Making target signal points for {allowed_masses} and {allowed_widths}")

    print("Creating ME reweighting signal templates")
    get_MEsig_templates()


    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
