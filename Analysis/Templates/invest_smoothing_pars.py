#! /bin/env python

"""
This script compares the output from the LOWESS smoothing algorithm using different input parameters for "frac" and "it"
Created by Joseph Dulemba on 1 September 2021
"""
# matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend("agg")
from matplotlib import rcParams
rcParams["font.size"] = 20
rcParams["savefig.format"] = "png"
rcParams["savefig.bbox"] = "tight"
from cycler import cycler
rcParams["axes.prop_cycle"] = cycler(color=["k", "#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33"])

from coffea.util import load, save
from pdb import set_trace
import os
import Utilities.plot_tools as plt_tools
import Utilities.prettyjson as prettyjson
from coffea import hist
import numpy as np
import fnmatch
import Utilities.Plotter as Plotter
import coffea.processor as processor    
import Utilities.systematics as systematics
import Utilities.final_analysis_binning as final_binning
    
base_jobid = os.environ["base_jobid"]
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016", "2017", "2018"] if base_jobid == "NanoAODv6" else ["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("--njets", default="all", nargs="?", choices=["3", "4+", "all"], help="Specify which jet multiplicity to use.")
args = parser.parse_args()

njets_to_run = []
if (args.njets == "3") or (args.njets == "all"):
    njets_to_run += ["3Jets"]
if (args.njets == "4+") or (args.njets == "all"):
    njets_to_run += ["4PJets"]



def find_chi2(h_fitted, h_unc):
    chi2 = np.sum(np.square(h_fitted.values()[()]-h_unc.values()[()]))/np.sum(h_unc.values(sumw2=True)[()][1])
    return chi2



def get_bkg_templates(fnames_to_run):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """

    #set_trace()
    for bkg_file in fnames_to_run:
        hdict = load(bkg_file)
        jmult = "3Jets" if "3Jets" in os.path.basename(bkg_file) else "4PJets"
        for lep in hdict.keys():
            for tname, orig_template in hdict[lep].items():
                
                proc = tname.split("_")[0] if not "data_obs" in tname else "data_obs"
                sys = sorted(filter(None, tname.split(f"{proc}_")))[0]

                #if not ((sys == "ueDOWN") and (proc == "ttJets")): continue
                if sys == "nosys": continue
                print(lep, jmult, sys, proc)

                nominal_hist = hdict[lep][f"{proc}_nosys"].copy()

                x_lims = (0, nominal_hist.dense_axes()[0].centers().size)

                    # perform smoothing
                smoothed_histos_list = [ (Plotter.smoothing_mttbins(nosys=nominal_hist, systematic=orig_template, mtt_centers=mtt_centers, nbinsx=len(linearize_binning[0])-1, nbinsy=len(linearize_binning[1])-1, **{"frac" : frac_val/10.}), frac_val/10.)
                    for frac_val in np.arange(2, 7, 2)
                ]
                #smoothed_histos_chi2 = {frac_val :  find_chi2(h_fitted=smooth_histo, h_unc=orig_template) for smooth_histo, frac_val in smoothed_histos_list}
                    # perform flattening
                flattened_histo = Plotter.flatten(nosys=nominal_hist, systematic=orig_template)
                #flat_chi2 = find_chi2(h_fitted=flattened_histo, h_unc=orig_template)

                # plot relative deviation
                fig, ax = plt.subplots()
                fig.subplots_adjust(hspace=.07)

                    # plot original dist
                orig_masked_vals, orig_masked_bins = Plotter.get_ratio_arrays(num_vals=orig_template.values()[()]-nominal_hist.values()[()], denom_vals=nominal_hist.values()[()], input_bins=nominal_hist.dense_axes()[0].edges())
                ax.fill_between(orig_masked_bins, orig_masked_vals, facecolor="k", step="post", alpha=0.5, label="Unsmoothed")

                    # plot smoothed versions
                for smooth_histo, frac_val in smoothed_histos_list:
                    smooth_masked_vals, smooth_masked_bins = Plotter.get_ratio_arrays(num_vals=smooth_histo.values()[()]-nominal_hist.values()[()], denom_vals=nominal_hist.values()[()], input_bins=nominal_hist.dense_axes()[0].edges())
                    ax.step(smooth_masked_bins, smooth_masked_vals, where="post", **{"linestyle":"-", "label":f"Frac={frac_val}", "linewidth":2})

                # plot flattened val
                flat_masked_vals, flat_masked_bins = Plotter.get_ratio_arrays(num_vals=flattened_histo.values()[()]-nominal_hist.values()[()], denom_vals=nominal_hist.values()[()], input_bins=nominal_hist.dense_axes()[0].edges())
                ax.step(flat_masked_bins, flat_masked_vals, where="post", **{"linestyle":"-", "label":"Flat", "linewidth":2})

                ax.legend(loc="upper right", title=f"{sys}, {proc}")
                ax.axhline(0, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
                ax.autoscale()
                ax.set_xlim(x_lims)
                ax.set_xlabel("$m_{t\\bar{t}}$ $\otimes$ |cos($\\theta^{*}_{t_{l}}$)|")
                ax.set_ylabel("Rel. Deviaton from Nominal")

                    # add lepton/jet multiplicity label
                ax.text(
                    0.02, 0.94, f"{leptypes[lep]}, {jet_mults[jmult]}",
                    fontsize=rcParams["font.size"]*0.9, horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                )
                    ## draw vertical lines for distinguishing different ctstar bins
                vlines = [x_lims[1]*ybin/5 for ybin in range(1, 5)]
                for vline in vlines:
                    ax.axvline(vline, color="k", linestyle="--")
                hep.cms.label(ax=ax, data=False, paper=False, year=args.year, lumi=round(data_lumi_year[f"{lep}s"]/1000., 1))

                #set_trace()
                pltdir = os.path.join(outdir, lep, jmult, sys)
                if not os.path.isdir(pltdir):
                    os.makedirs(pltdir)

                #figname = os.path.join(pltdir, "_".join([jmult, lep, sys, proc, "BinWidths_Comp"]))
                #figname = os.path.join(pltdir, "_".join([jmult, lep, sys, proc, "SmoothValues_Comp"]))
                #figname = os.path.join(pltdir, "_".join([jmult, lep, sys, proc, "MttBinWidths_SmoothValues_Comp"]))
                figname = os.path.join(pltdir, "_".join([jmult, lep, sys, proc, "SmoothedFlatVals_Comp"]))
                fig.savefig(figname)
                print(f"{figname} written")
                plt.close()


if __name__ == "__main__":
    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]

    analyzer = "htt_btag_sb_regions"
    base_template_name = f"test_raw_templates_lj_NJETS_bkg_{args.year}_{jobid}"
    #base_template_name = "raw_templates_lj_NJETS_SIG_%s_%s" % (jobid, args.year)

        # get matching pattern based on args.njets
    njets_regex = "*" if len(njets_to_run) > 1 else njets_to_run[0]

    input_dir = os.path.join(proj_dir, "results", "%s_%s" % (args.year, jobid), f"Templates_{analyzer}")
    if os.path.isdir(input_dir):
            # define variables to get histogram for background    
        bkg_fnmatch = "%s.coffea" % base_template_name.replace("NJETS", njets_regex)
        bkg_fnames = fnmatch.filter(os.listdir(input_dir), bkg_fnmatch)
        bkg_fnames = [os.path.join(input_dir, fname) for fname in bkg_fnames]
    else: print("No background file found.")

    data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, "inputs", "%s_lumis_data.json" % base_jobid)).read())[args.year]
    outdir = os.path.join(proj_dir, "plots", "%s_%s" % (args.year, jobid), "Templates_%s" % analyzer, "Smoothing_Invest")
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    
    jet_mults = {
        "3Jets" : "3 jets",
        "4PJets" : "4+ jets"
    }
    
    leptypes = {
        "Muon" : "$\\mu$",
        "Electron" : "$e$",
    }

    linearize_binning = (
        final_binning.mtt_binning,
        final_binning.ctstar_abs_binning
    )
    #linearize_binning = (
    #    #np.array([300.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 625.0, 650.0, 675.0, 700.0, 730.0, 760.0, 800.0, 850.0, 900.0, 1000.0, 1200.0]), # orig
    #    np.array([300.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 625.0, 650.0, 675.0,
    #        700.0, 730.0, 760.0, 800.0, 850.0, 900.0, 950., 1000.0, 1050., 1100.0, 1150., 1200., 1300., 1500., 2000.]),
    #    #np.array([300.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 625.0, 650.0, 675.0,
    #    #    700.0, 730.0, 760.0, 800.0, 850.0, 900.0, 950., 1000.0, 1050., 1100.0, 1150., 1200., 1300., 1500., 2000.]),
    #    #np.array([300.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 620.0, 650., 700.0, 750.0, 800.0, 850.0, 900.0, 2000.0]),
    #    np.array([0.0, 0.4, 0.6, 0.75, 0.9, 1.0])
    #)
    mtt_centers =  [(linearize_binning[0][idx]+linearize_binning[0][idx+1])/2 for idx in range(len(linearize_binning[0])-1)]

    try:
        print("Creating background templates")
        get_bkg_templates(bkg_fnames)
        
    except:
        print("Could not write background templates to file")
