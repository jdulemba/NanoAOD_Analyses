#! /bin/env python

from coffea.util import load#, save
from pdb import set_trace
import os
import Utilities.prettyjson as prettyjson
import numpy as np
import fnmatch
import Utilities.systematics as systematics
import Utilities.Plotter as Plotter
from scipy.stats import chisquare

# matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend("agg")
from matplotlib import rcParams
rcParams["font.size"] = 20
rcParams["savefig.format"] = "png"
rcParams["savefig.bbox"] = "tight"

base_jobid = os.environ["base_jobid"]    
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016", "2017", "2018"] if base_jobid == "NanoAODv6" else ["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("lepton", choices=["Electron", "Muon"], help="Choose which lepton to make plots for")
parser.add_argument("--njets", default="all", nargs="?", choices=["3", "4+", "all"], help="Specify which jet multiplicity to use.")
args = parser.parse_args()

njets_to_run = []
if (args.njets == "3") or (args.njets == "all"):
    njets_to_run += ["3Jets"]
if (args.njets == "4+") or (args.njets == "all"):
    njets_to_run += ["4PJets"]


def plot_bkg_templates(fnames_to_run):
    """
    Runs LOWESS smoothing algorithm ntoys times and finds 1 and 2 sigma bands for interpolation
    """

    for bkg_file in fnames_to_run:
        hdict = load(bkg_file)
        jmult = "3Jets" if "3Jets" in os.path.basename(bkg_file) else "4PJets"
        for tname, orig_template in hdict[args.lepton].items():

            proc = tname.split("_")[0] if not "data_obs" in tname else "data_obs"
            sys = sorted(filter(None, tname.split(f"{proc}_")))[0]

            if proc == "BKG": continue
            #if sys not in ["hdampUP", "hdampDOWN", "mtop1665", "mtop1695", "mtop1715", "mtop1735", "mtop1755", "mtop1785", "ueUP", "ueDOWN"]: continue
            if sys == "nosys": continue
            print(args.lepton, jmult, sys, proc)

            nosys_hist = hdict[args.lepton][f"{proc}_nosys"].copy()
            orig_smooth_hist = Plotter.smoothing_mttbins(nosys=nosys_hist, systematic=orig_template, mtt_centers=mtt_centers, nbinsx=nbinsx, nbinsy=nbinsy)

            x_lims = (0, nosys_hist.dense_axes()[0].centers().size)

                # get vals and errors of systematic variation
            sys_histo_vals, sys_histo_sumw2 = orig_template.values(sumw2=True)[()]
            sys_histo_errs = np.sqrt(sys_histo_sumw2)

                # make toys based on Gaussian distribution of mu=bin_val, sigma=bin_error
            toy_arrays = np.zeros((nbins, ntoys))
            for idx in range(nbins):
                toy_arrays[idx] = np.random.normal(sys_histo_vals[idx], sys_histo_errs[idx], size=ntoys)

                # get smoothed relative deviation distributions from toys
            smoothed_rel_dev_arrays = np.zeros((ntoys, nbins))
            chi2_pvals = np.zeros((ntoys, 2))
            for idx in range(ntoys):
                smoothed_array = Plotter.smoothing_mttbins(nosys=nosys_hist, systematic=(toy_arrays.T)[idx], mtt_centers=mtt_centers, nbinsx=nbinsx, nbinsy=nbinsy)
                chi2_pval = chisquare(f_obs=smoothed_array, f_exp=orig_smooth_hist.values()[()]) # convert to expected yields so inputs are greater than 5
                chi2_pvals[idx] = np.array([chi2_pval.statistic, chi2_pval.pvalue])
                smoothed_rel_dev_arrays[idx] = (smoothed_array - nosys_hist.values()[()])/nosys_hist.values()[()]

                ## find 68% and 95% intervals
            plus_one_sigma_smooth_vals, minus_one_sigma_smooth_vals = np.zeros(nbins), np.zeros(nbins)
            plus_two_sigma_smooth_vals, minus_two_sigma_smooth_vals = np.zeros(nbins), np.zeros(nbins)
            for bin in range(nbins):
                plus_one_sigma_smooth_vals[bin] = np.sort(smoothed_rel_dev_arrays[:,bin])[plus_one_sigma_ind]
                minus_one_sigma_smooth_vals[bin] = np.sort(smoothed_rel_dev_arrays[:,bin])[minus_one_sigma_ind]
                plus_two_sigma_smooth_vals[bin] = np.sort(smoothed_rel_dev_arrays[:,bin])[plus_two_sigma_ind]
                minus_two_sigma_smooth_vals[bin] = np.sort(smoothed_rel_dev_arrays[:,bin])[minus_two_sigma_ind]


            # plot relative deviation
            fig, ax = plt.subplots()
            fig.subplots_adjust(hspace=.07)

                # original relative deviations
            orig_masked_vals, orig_masked_bins = Plotter.get_ratio_arrays(num_vals=orig_template.values()[()]-nosys_hist.values()[()], denom_vals=nosys_hist.values()[()], input_bins=nosys_hist.dense_axes()[0].edges())
            ax.step(orig_masked_bins, orig_masked_vals, where="post", **{"color":"k", "linestyle":"-", "label":"Original"})
                # original smoothing relative deviations
            orig_smoothed_masked_vals, orig_smoothed_masked_bins = Plotter.get_ratio_arrays(num_vals=orig_smooth_hist.values()[()]-nosys_hist.values()[()], denom_vals=nosys_hist.values()[()], input_bins=nosys_hist.dense_axes()[0].edges())
            ax.step(orig_smoothed_masked_bins, orig_smoothed_masked_vals, where="post", **{"color":"r", "linestyle":"-", "label":"Original Smoothing"})
                # plot 68 and 95% intervals for yields 
            ax.fill_between(nosys_hist.dense_axes()[0].edges(), np.r_[minus_one_sigma_smooth_vals, minus_one_sigma_smooth_vals[-1]], np.r_[plus_one_sigma_smooth_vals, plus_one_sigma_smooth_vals[-1]],
                where=np.r_[plus_one_sigma_smooth_vals, plus_one_sigma_smooth_vals[-1]] > np.r_[minus_one_sigma_smooth_vals, minus_one_sigma_smooth_vals[-1]], step="post", **{"label":"68%", "facecolor":"#00cc00", "alpha":0.5})
            ax.fill_between(nosys_hist.dense_axes()[0].edges(), np.r_[minus_two_sigma_smooth_vals, minus_two_sigma_smooth_vals[-1]], np.r_[plus_two_sigma_smooth_vals, plus_two_sigma_smooth_vals[-1]],
                where=np.r_[plus_two_sigma_smooth_vals, plus_two_sigma_smooth_vals[-1]] > np.r_[minus_two_sigma_smooth_vals, minus_two_sigma_smooth_vals[-1]], step="post", **{"label":"95%", "facecolor":"#ffcc00", "alpha":0.5})

            ax.legend(loc="upper right", title=f"{sys}, {proc}")
            ax.axhline(0, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
            ax.autoscale()
            ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.15)
            ax.set_xlim(x_lims)
            ax.set_xlabel("$m_{t\\bar{t}}$ $\otimes$ |cos($\\theta^{*}_{t_{l}}$)|")
            ax.set_ylabel("Rel. Deviaton from Nominal")

                # add lepton/jet multiplicity label
            ax.text(
                0.02, 0.94, f"{leptypes[args.lepton]}, {jet_mults[jmult]}",
                fontsize=rcParams["font.size"]*0.9, horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
            )
                ## draw vertical lines for distinguishing different ctstar bins
            vlines = [x_lims[1]*ybin/5 for ybin in range(1, 5)]
            for vline in vlines:
                ax.axvline(vline, color="k", linestyle="--")
            hep.cms.label(ax=ax, data=False, paper=False, year=args.year, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))

            #set_trace()
            pltdir = os.path.join(outdir, args.lepton, jmult, sys)
            if not os.path.isdir(pltdir):
                os.makedirs(pltdir)

            figname = os.path.join(pltdir, "_".join([jmult, args.lepton, sys, proc, "SmoothingConfidenceIntervals"]))
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

    outdir = os.path.join(proj_dir, "plots", "%s_%s" % (args.year, jobid), "Templates_%s" % analyzer, "Smoothing_Checks")
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    baseSys = lambda sys : "_".join(sys.split("_")[:-1])

        # define variables used for plotting
    data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, "inputs", "%s_lumis_data.json" % base_jobid)).read())[args.year]

    jet_mults = {
        "3Jets" : "3 jets",
        "4PJets" : "4+ jets"
    }
    
    leptypes = {
        "Muon" : "$\\mu$",
        "Electron" : "$e$",
    }

        # define variables used for smoothing
    linearize_binning = (
        #np.array([300.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 625.0, 650.0, 675.0, 700.0, 730.0, 760.0, 800.0, 850.0, 900.0, 1000.0, 1200.0]), # orig
        np.array([300.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 625.0, 650.0, 675.0,
            700.0, 730.0, 760.0, 800.0, 850.0, 900.0, 950., 1000.0, 1050., 1100.0, 1150., 1200., 1300., 1500., 2000.]),
        #np.array([300.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 620.0, 650., 700.0, 750.0, 800.0, 850.0, 900.0, 2000.0]),
        np.array([0.0, 0.4, 0.6, 0.75, 0.9, 1.0])
    )
    mtt_centers =  [(linearize_binning[0][idx]+linearize_binning[0][idx+1])/2 for idx in range(len(linearize_binning[0])-1)]

    nbinsx, nbinsy = len(linearize_binning[0])-1, len(linearize_binning[1])-1
    nbins = nbinsx*nbinsy
    #ntoys = 10
    ntoys = 1000
    plus_one_sigma_ind, minus_one_sigma_ind = int((ntoys/2)+(ntoys*0.34))-1, int((ntoys/2)-(ntoys*0.34))-1
    plus_two_sigma_ind, minus_two_sigma_ind = int((ntoys/2)+(ntoys*0.475))-1, int((ntoys/2)-(ntoys*0.475))-1

    try:
        print(f"Running smoothing algorith using {ntoys} for background templates")
        plot_bkg_templates(bkg_fnames)

    except:
        print("Could not perform smoothing")
