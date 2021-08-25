#!/usr/bin/env python

# matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend("agg")
from matplotlib import rcParams
rcParams["font.size"] = 20
rcParams["savefig.format"] = "png"
rcParams["savefig.bbox"] = "tight"

import Utilities.Plotter as Plotter
from coffea.util import load, save
from pdb import set_trace
import os
import numpy as np
from coffea.lookup_tools.dense_lookup import dense_lookup
from rootpy.plotting import Hist2D, Hist
import rootpy
rootpy.log["/"].setLevel(rootpy.log.INFO)

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
base_jobid = os.environ["base_jobid"]

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--save_ratios", action="store_true", help="Save NNLO/tune ratios for all years")
args = parser.parse_args()

fname = os.path.join(proj_dir, "NNLO_files", "NNLO_to_Tune_Ratios_%s.coffea" % base_jobid)
ratios = load(fname)

var_opts = {
    "thad_pt" : {"xtitle" : "$p_{T}$($t_{h}$) [GeV]", "ytitle" : "NNLO/Tune"},
    "mtt_vs_thad_ctstar" : {"xtitle" : "m($t\\bar{t}$) [GeV]", "ytitle" : "cos($\\theta^{*}_{t_{h}}$)"},
}
year_opts = {
    "2016" : {"leg_title" : "2016 T4" if base_jobid == "NanoAODv6" else "2016 post-VFP", "col" : "#377eb8"}, ## blue
    "2017" : {"leg_title" : "2017", "col" : "#e41a1c"}, ## red
    "2018" : {"leg_title" : "2018", "col" : "#4daf4a"}, ## green
}
if base_jobid != "NanoAODv6":
    year_opts["2016APV"] = {"leg_title" : "2016 pre-VFP", "col" : "#984ea3"} # purple

outdir = os.path.join(proj_dir, "plots", base_jobid, "make_nnlo_reweighting_vars", "Interp")
if not os.path.isdir(outdir): os.makedirs(outdir)

years_to_run = ["2016", "2017", "2018"] if base_jobid == "NanoAODv6" else ["2016APV", "2016", "2017", "2018"]

for var in var_opts.keys():
    if var == "thad_pt":
        fig, ax = plt.subplots()
        fig.subplots_adjust(hspace=.07)
        ax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
        ax.text(
            0.02, 0.90, "$t\\bart \\rightarrow e/\mu + jets$\nparton level",
            fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
        )
    
    #set_trace()
    for year in years_to_run:
        histo = ratios[year][var] # get histogram

        if histo._dimension == 1:
            orig_bins, orig_vals = histo._axes, histo._values
                # plot original distribution
            Plotter.plot_1D(orig_vals, orig_bins, xlabel=var_opts[var]["xtitle"], ylabel=var_opts[var]["ytitle"], color=year_opts[year]["col"], ax=ax, label="%s Original" % year_opts[year]["leg_title"])

                # get interpolated values from ROOT Interpolate
            orig_ratio_hist = Hist(orig_bins, name="", title="")
            for xbin in range(1, orig_bins.size):
                orig_ratio_hist[xbin] = orig_vals[xbin-1]

            output_bins = np.arange(min(orig_bins), max(orig_bins)+10, 10)
            interped_array = np.zeros(output_bins.size-1)
            for xbin in range(output_bins.size-1):
                interped_array[xbin] = orig_ratio_hist.Interpolate(output_bins[xbin])

            Plotter.plot_1D(interped_array, output_bins, xlabel=var_opts[var]["xtitle"], ylabel=var_opts[var]["ytitle"], color=year_opts[year]["col"], ax=ax, label="%s Interp" % (year_opts[year]["leg_title"]), linestyle="--")

            if args.save_ratios:
                lookup = dense_lookup(*(interped_array, output_bins))
                ratios[year].update({"%s_Interp" % var : lookup})

        else:
            (orig_xbins, orig_ybins), orig_vals = histo._axes, histo._values
            orig_ratio_hist = Hist2D(orig_xbins, orig_ybins, name="mtt_ctstar", title="mtt_ctstar")
            for ybin in range(1, orig_ybins.size):
                for xbin in range(1, orig_xbins.size):
                    orig_ratio_hist[xbin, ybin] = orig_vals[xbin-1][ybin-1]

            output_xbins = np.arange(min(orig_xbins), max(orig_xbins)+10, 10)
            output_ybins = np.arange(min(orig_ybins), max(orig_ybins)+0.05, 0.05)
            interped_array = np.zeros((output_xbins.size-1, output_ybins.size-1))
            for ybin in range(output_ybins.size-1):
                for xbin in range(output_xbins.size-1):
                    interped_array[xbin, ybin] = orig_ratio_hist.Interpolate(output_xbins[xbin], output_ybins[ybin])
            
            if args.save_ratios:
                lookup = dense_lookup(*(interped_array, (output_xbins, output_ybins)))
                ratios[year].update({"%s_Interp" % var : lookup})

                # plot original distribution
            fig, ax = plt.subplots()
            fig.subplots_adjust(hspace=.07)
            Plotter.plot_2d_norm(0, xbins=orig_xbins, ybins=orig_ybins, values=orig_vals,
                xlimits=(min(orig_xbins), max(orig_xbins)), ylimits=(min(orig_ybins), max(orig_ybins)),
                xlabel=var_opts[var]["xtitle"], ylabel=var_opts[var]["ytitle"], ax=ax, **{"cmap_label" : "NNLO/Tune"},
            )
            ax.text(
                0.02, 0.90, "$t\\bart \\rightarrow e/\mu + jets$\nparton level",
                fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
            )
            ax.text(0.98, 0.95, "Original", fontsize=rcParams["font.size"], horizontalalignment="right", verticalalignment="bottom", transform=ax.transAxes)
            hep.cms.label(ax=ax, year=year)

            figname = os.path.join(outdir, "%s_%s_ttJets_Orig_NNLO_to_Tune_%s" % (year, jobid, var))
            fig.savefig(figname)
            print(f"{figname} written")

                # plot interpolated distribution
            fig_int, ax_int = plt.subplots()
            fig_int.subplots_adjust(hspace=.07)
            Plotter.plot_2d_norm(0, xbins=output_xbins, ybins=output_ybins, values=interped_array,
                xlimits=(min(orig_xbins), max(orig_xbins)), ylimits=(min(orig_ybins), max(orig_ybins)),
                xlabel=var_opts[var]["xtitle"], ylabel=var_opts[var]["ytitle"], ax=ax_int, **{"cmap_label" : "NNLO/Tune"},
            )
            ax_int.text(
                0.02, 0.90, "$t\\bart \\rightarrow e/\mu + jets$\nparton level",
                fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax_int.transAxes
            )
            ax_int.text(0.98, 0.95, "Interpolated", fontsize=rcParams["font.size"], horizontalalignment="right", verticalalignment="bottom", transform=ax_int.transAxes)
            hep.cms.label(ax=ax_int, year=year)

            figname_int = os.path.join(outdir, "%s_%s_ttJets_Interp_NNLO_to_Tune_%s" % (year, jobid, var))
            fig_int.savefig(figname_int)
            print(f"{figname_int} written")
            plt.close()


    if var == "thad_pt":
            # format axes
        ax.autoscale(axis="x", tight=True)
        ax.set_ylim(0.6, ax.get_ylim()[1])
        ax.legend(loc="upper right")
        hep.cms.label(ax=ax, year="")

        figname = os.path.join(outdir, "%s_ttJets_Interpolated_NNLO_to_Tune_%s_Comp" % (jobid, var))
        fig.savefig(figname)
        print(f"{figname} written")
        plt.close()

if args.save_ratios:
    orig_interp_ratios_fname = os.path.join(proj_dir, "NNLO_files", "NNLO_to_Tune_Orig_Interp_Ratios_%s.coffea" % base_jobid)
    save(ratios, orig_interp_ratios_fname)
    print(f"\n{orig_interp_ratios_fname} written")
