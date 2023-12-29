#! /bin/env python

import time
tic = time.time()

# matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend("agg")
from matplotlib import rcParams
rcParams["font.size"] = 20
rcParams["savefig.format"] = "pdf"
rcParams["savefig.bbox"] = "tight"

from coffea.util import load
from pdb import set_trace
import os
import numpy as np
import Utilities.Plotter as Plotter
import coffea.processor as processor    
import Utilities.final_analysis_binning as final_binning
import Utilities.prettyjson as prettyjson
import Utilities.common_features as cfeatures
from Utilities import systematics

import argparse
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("template_type", choices=["Indiv", "Combined_Era_Lepton", "Combined_Lepton"], help="What type of template do you want to make?")
parser.add_argument("--year", choices=["2016APV", "2016", "2017", "2018"], help="Specify which year to use if combination isn't across year and lepton channels.")
args = parser.parse_args()


if (args.template_type == "Combined_Lepton") and (not args.year):
    raise ValueError("Year must be specified when the template_type is only across leptons")

if (args.template_type == "Indiv") and (not args.year):
    raise ValueError("Year must be specified when the individual channel templates are used")


version = "V2"
#version = "V1"



## make plots of systematic uncs
def make_combined_era_lepton_plots(hdict):
    year_to_use = cfeatures.year_labels["Total"]
    lumi_to_use = (data_lumi["TOT"]["Electrons"] + data_lumi["TOT"]["Muons"])/2000.


        # loop over jmults
    for jmult in hdict.keys():
        plt_outdir = os.path.join(plot_outdir, jobid, "Toponium_HBSR", version, "Combined_Era_Lepton", jmult)
        if not os.path.isdir(plt_outdir):
            os.makedirs(plt_outdir)    

        nosys = hdict[jmult]["Toponium_nosys"].copy()

        systs = sorted(set([key.replace("Toponium_", "").replace("_Down", "").replace("_DW", "").replace("Down", "").replace("_down", "").replace("_Up", "").replace("_UP", "").replace("Up", "").replace("_up", "") for key in hdict[jmult].keys()]))
        for sys in systs:
            if sys not in systematics.combined_era_lepton: continue
            fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True, figsize=(15.0, 10.0))
            fig.subplots_adjust(hspace=.07)

            #set_trace()
                ## plot yields
            up_sysname = [key for key in hdict[jmult].keys() if (sys in key) and ("UP" in key.upper())][0]
            dw_sysname = [key for key in hdict[jmult].keys() if (sys in key) and ("DW" in key.upper() or "DOWN" in key.upper())][0]
            up_masked_vals, up_masked_bins = Plotter.get_ratio_arrays(num_vals=hdict[jmult][up_sysname].values()[()], denom_vals=np.ones_like(hdict[jmult][up_sysname].values()[()]), input_bins=nosys.dense_axes()[0].edges())
            ax.step(up_masked_bins, up_masked_vals, where="post", **{"color" : "r", "label" : "Up"})
            ax.step(up_masked_bins, np.r_[nosys.values()[()], nosys.values()[()][-1]], where="post", **{"color" : "k", "label" : "Nominal"})
            dw_masked_vals, dw_masked_bins = Plotter.get_ratio_arrays(num_vals=hdict[jmult][dw_sysname].values()[()], denom_vals=np.ones_like(hdict[jmult][dw_sysname].values()[()]), input_bins=nosys.dense_axes()[0].edges())
            ax.step(dw_masked_bins, dw_masked_vals, where="post", **{"color" : "b", "label" : "Down"})

                # plot ratios
            up_ratio_masked_vals, up_ratio_masked_bins = Plotter.get_ratio_arrays(num_vals=hdict[jmult][up_sysname].values()[()], denom_vals=nosys.values()[()], input_bins=nosys.dense_axes()[0].edges())
            dw_ratio_masked_vals, dw_ratio_masked_bins = Plotter.get_ratio_arrays(num_vals=hdict[jmult][dw_sysname].values()[()], denom_vals=nosys.values()[()], input_bins=nosys.dense_axes()[0].edges())
            rax.step(up_ratio_masked_bins, up_ratio_masked_vals, where="post", **{"color" : "r"})
            rax.step(dw_ratio_masked_bins, dw_ratio_masked_vals, where="post", **{"color" : "b"})
            #rax.fill_between(up_masked_bins, up_masked_vals, y2=1., facecolor=up_style["color"], label=up_style["label"], step="post", alpha=0.5) if use_fill_between\
            #    else rax.step(up_masked_bins, up_masked_vals, where="post", **up_style)

            ax.legend(loc="upper right", title="$\\eta_{t}$", ncol=1)
            ax.autoscale()
            ax.set_ylabel("Events")
            ax.set_ylim(0, ax.get_ylim()[1]*1.3)
            ax.set_xlabel(None)
            ax.set_xlim(x_lims)

            rax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
            rax.autoscale()
            rax.set_ylabel("Ratio to Nominal")
            rax.set_ylim(max(rax.get_ylim()[0], 0.5), min(rax.get_ylim()[1], 1.5))
            rax.set_xlim(x_lims)
            rax.set_xlabel("$m_{t\\bar{t}}$ [GeV]")

                # add lepton/jet multiplicity label
            ax.text(
                0.02, 0.86, f"{sys}\n{cfeatures.channel_labels[f'Lepton_{jmult}']}",
                fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
            )
                ## draw vertical lines for distinguishing different ctstar bins
            vlines = [x_lims[1]*ybin/5 for ybin in range(1, 5)]
            [ax.axvline(vline, color="k", linestyle="--") for vline in vlines]
            [rax.axvline(vline, color="k", linestyle="--") for vline in vlines]

            for idx, label in enumerate(ctstar_binlabels):
                ax.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                    xytext=(0, 25), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)

            ax.set_xticks(mtt_bin_inds_to_plot)
            ax.set_xticklabels(mtt_bins_to_plot)
            rax.set_xticks(mtt_bin_inds_to_plot)
            rax.set_xticklabels(mtt_bins_to_plot)
            hep.cms.label(ax=ax, data=False, label="Preliminary", year=year_to_use, lumi=int(round(lumi_to_use, 0)))

            figname = os.path.join(plt_outdir, "_".join(["Toponium", jmult, "Combined_Era_Lepton", sys,"Comp"]))
            fig.savefig(figname)
            print(f"{figname} written")
            plt.close(fig)



def make_indiv_plots(fname):
    hdict = load(fname)[args.year]
    year_to_use = cfeatures.year_labels[args.year]

        # loop over jmults
    for jmult in hdict.keys():
            # loop over leptons
        for lep in hdict[jmult].keys():
            plt_outdir = os.path.join(plot_outdir, jobid, "Toponium_HBSR", version, "Indiv", args.year, lep, jmult)
            if not os.path.isdir(plt_outdir):
                os.makedirs(plt_outdir)    

            lumi_to_use = (data_lumi[args.year][f"{lep}s"])/1000.
    
            nosys = hdict[jmult][lep]["Toponium_nosys"].copy()

            #set_trace()
            systs = sorted(set([key.replace("Toponium_", "").replace("_Down", "").replace("_DW", "").replace("Down", "").replace("_down", "").replace("_Up", "").replace("_UP", "").replace("Up", "").replace("_up", "") for key in hdict[jmult][lep].keys()]))
            for sys in systs:
                if sys == "nosys": continue
                #if sys == "BTAG_BC_MUPT": set_trace()
                #if sys not in systematics.combined_lepton[args.year]: continue
                fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True, figsize=(15.0, 10.0))
                fig.subplots_adjust(hspace=.07)

                #set_trace()
                    ## plot yields
                up_sysname = [key for key in hdict[jmult][lep].keys() if (sys in key) and (("_UP" in key) or ("_Up" in key) or ("Up" in key) or ("_up" in key) )][0]
                #up_sysname = [key for key in hdict[jmult][lep].keys() if (sys in key) and ("UP" in key.upper())][0]
                dw_sysname = [key for key in hdict[jmult][lep].keys() if (sys in key) and (("_Down" in key) or ("_DW" in key) or ("Down" in key) or ("_down" in key) )][0]
                #dw_sysname = [key for key in hdict[jmult][lep].keys() if (sys in key) and ("DW" in key.upper() or "DOWN" in key.upper())][0]
                up_masked_vals, up_masked_bins = Plotter.get_ratio_arrays(num_vals=hdict[jmult][lep][up_sysname].values()[()], denom_vals=np.ones_like(hdict[jmult][lep][up_sysname].values()[()]), input_bins=nosys.dense_axes()[0].edges())
                dw_masked_vals, dw_masked_bins = Plotter.get_ratio_arrays(num_vals=hdict[jmult][lep][dw_sysname].values()[()], denom_vals=np.ones_like(hdict[jmult][lep][dw_sysname].values()[()]), input_bins=nosys.dense_axes()[0].edges())
                ax.step(up_masked_bins, up_masked_vals, where="post", **{"color" : "r", "label" : "Up"})
                ax.step(up_masked_bins, np.r_[nosys.values()[()], nosys.values()[()][-1]], where="post", **{"color" : "k", "label" : "Nominal"})
                ax.step(dw_masked_bins, dw_masked_vals, where="post", **{"color" : "b", "label" : "Down"})

                    # plot ratios
                up_ratio_masked_vals, up_ratio_masked_bins = Plotter.get_ratio_arrays(num_vals=hdict[jmult][lep][up_sysname].values()[()], denom_vals=nosys.values()[()], input_bins=nosys.dense_axes()[0].edges())
                dw_ratio_masked_vals, dw_ratio_masked_bins = Plotter.get_ratio_arrays(num_vals=hdict[jmult][lep][dw_sysname].values()[()], denom_vals=nosys.values()[()], input_bins=nosys.dense_axes()[0].edges())
                rax.step(up_ratio_masked_bins, up_ratio_masked_vals, where="post", **{"color" : "r"})
                rax.step(dw_ratio_masked_bins, dw_ratio_masked_vals, where="post", **{"color" : "b"})
                #rax.fill_between(up_masked_bins, up_masked_vals, y2=1., facecolor=up_style["color"], label=up_style["label"], step="post", alpha=0.5) if use_fill_between\
                #    else rax.step(up_masked_bins, up_masked_vals, where="post", **up_style)

                ax.legend(loc="upper right", title="$\\eta_{t}$", ncol=1)
                ax.autoscale()
                ax.set_ylabel("Events")
                ax.set_ylim(0, ax.get_ylim()[1]*1.3)
                ax.set_xlabel(None)
                ax.set_xlim(x_lims)
    
                rax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
                rax.autoscale()
                rax.set_ylabel("Ratio to Nominal")
                rax.set_ylim(max(rax.get_ylim()[0], 0.5), min(rax.get_ylim()[1], 1.5))
                rax.set_xlim(x_lims)
                rax.set_xlabel("$m_{t\\bar{t}}$ [GeV]")
    
                    # add lepton/jet multiplicity label
                ax.text(
                    0.02, 0.86, f"{sys}\n{cfeatures.channel_labels[f'{lep}_{jmult}']}",
                    fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                )
                    ## draw vertical lines for distinguishing different ctstar bins
                vlines = [x_lims[1]*ybin/5 for ybin in range(1, 5)]
                [ax.axvline(vline, color="k", linestyle="--") for vline in vlines]
                [rax.axvline(vline, color="k", linestyle="--") for vline in vlines]
    
                for idx, label in enumerate(ctstar_binlabels):
                    ax.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                        xytext=(0, 25), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)
    
                ax.set_xticks(mtt_bin_inds_to_plot)
                ax.set_xticklabels(mtt_bins_to_plot)
                rax.set_xticks(mtt_bin_inds_to_plot)
                rax.set_xticklabels(mtt_bins_to_plot)
                hep.cms.label(ax=ax, data=False, label="Preliminary", year=year_to_use, lumi=round(lumi_to_use, 1))
    
                figname = os.path.join(plt_outdir, "_".join(["Toponium", lep, jmult, "Indiv", args.year, sys, "Comp"]))
                fig.savefig(figname)
                print(f"{figname} written")
                plt.close(fig)




if __name__ == "__main__":
    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]
    base_jobid = os.environ["base_jobid"]
    eos_dir = os.environ["eos_dir"]
    plot_outdir = os.environ["plots_dir"]

    linearize_binning = (
        final_binning.mtt_binning,
        final_binning.ctstar_abs_binning
    ) 

    ctstar_binlabels = [r"%s $\in [%s, %s)$" % (cfeatures.variable_names_to_labels["tlep_ctstar_abs"], linearize_binning[1][bin], linearize_binning[1][bin+1]) for bin in range(len(linearize_binning[1])-1)]
    ctstar_binlabels[-1] = r"%s $\in [%s, %s]$" % (cfeatures.variable_names_to_labels["tlep_ctstar_abs"], linearize_binning[1][-2], linearize_binning[1][-1])
    ctstar_bin_locs = np.linspace((len(linearize_binning[0])-1)/2, (len(linearize_binning[0])-1)*(len(linearize_binning[1])-1) - (len(linearize_binning[0])-1)/2, len(linearize_binning[1])-1)
    mtt_vals_to_plot = np.array([400, 600, 1000])
    mtt_tiled_labels = np.tile(linearize_binning[0][:-1], linearize_binning[1].size-1)
    mtt_bin_inds_to_plot = np.where(np.in1d(mtt_tiled_labels, mtt_vals_to_plot))[0]
    mtt_bins_to_plot = np.tile(mtt_vals_to_plot, linearize_binning[1].size-1)

    x_lims = (0, (linearize_binning[1].size - 1)* (linearize_binning[0].size - 1))

    data_lumi = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())

    input_dir = os.path.join(eos_dir, "results", jobid, "Toponium_HBSR", "Templates", version)
    if args.template_type == "Indiv":
        cfile = os.path.join(input_dir, f"raw_templates_lj_toponium_{jobid}_{version}.coffea")
        if not os.path.isfile(cfile):
            raise ValueError(f"{cfile} not found.")
        make_indiv_plots(cfile)

    ##if args.template_type == "Combined_Lepton":
    #        # make sure raw template file exists
    #    raw_tfile = os.path.join(outdir, f"raw_templates_lj_toponium_{jobid}_{version}.coffea")
    #    if not os.path.isfile(raw_tfile):
    #        make_raw_templates()
    #    combine_lepton_templates(load(raw_tfile))

    if args.template_type == "Combined_Era_Lepton":
        cfile = os.path.join(input_dir, f"raw_combined_era_lepton_templates_lj_toponium_{jobid}_{version}.coffea")
        if not os.path.isfile(cfile):
            raise ValueError(f"{cfile} not found.")
        make_combined_era_lepton_plots(load(cfile))

    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
