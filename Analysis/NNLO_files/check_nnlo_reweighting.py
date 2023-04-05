#!/usr/bin/env python

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

import Utilities.Plotter as Plotter
from coffea.hist import plot
from coffea import hist
from coffea.util import load
from pdb import set_trace
import os
import Utilities.prettyjson as prettyjson
import Utilities.plot_tools as plt_tools
import numpy as np
import Utilities.final_analysis_binning as final_binning
import Utilities.common_features as cfeatures

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
base_jobid = os.environ["base_jobid"]
analyzer = "check_nnlo_reweighting"
eos_dir = os.environ["eos_dir"]
plots_dir = os.environ["plots_dir"]

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="Specify which year to run over")
args = parser.parse_args()

f_ext = "TOT.coffea"
input_dir = os.path.join(eos_dir, "results", f"{args.year}_{base_jobid}", analyzer)
fnames = sorted(["%s/%s" % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)])
hdict = plt_tools.add_coffea_files(fnames) if len(fnames) > 1 else load(fnames[0])
    
outdir = os.path.join(plots_dir, base_jobid, analyzer, args.year)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

#set_trace()    
rewt_style_dict = {
    "Nominal" : {"label" : "Powheg", "color" : "k", "linestyle" : "-"},
    "mtt_vs_top_ctstar" : {"label" : "Powheg x $W$($m_{t\\bar{t}}$, cos($\\theta^{*}_{t}$))", "color" : "#377eb8", "linestyle" : "-"}, ## blue
    "top_pt" : {"label" : "Powheg x $W$($p_{T}$($t$))", "color" : "#e42a2c", "linestyle" : "-"}, ## red
}
#rewt_style_dict = {
#    "Nominal" : ("Nominal", "k"),
#    #"thad_pt" : ("$W_{Orig}$($p_{T}$($t_{h}$))", "#e42a2c"), ## red
#    #"thad_pt_Interp" : ("$W_{Int}$($p_{T}$($t_{h}$))", "#4daf4a"), ## green
#    "mtt_vs_thad_ctstar" : ("$W_{Orig}$($m_{t\\bar{t}}$, cos($\\theta^{*}_{t_{h}}$))", "#377eb8"), ## blue
#    #"mtt_vs_thad_ctstar_Interp" : ("$W_{Int}$($m_{t\\bar{t}}$, cos($\\theta^{*}_{t_{h}}$))", "#ff7f00"), ## orange
#}

variables = {
        # top
    "pt_top" : (cfeatures.variable_names_to_labels["pt_top"], 2, (0., 500.)),
    "eta_top" : (cfeatures.variable_names_to_labels["eta_top"], 2, (-4., 4.)),
    "phi_top" : (cfeatures.variable_names_to_labels["phi_top"], 1, (-3.2, 3.2)),
    "mass_top" : (cfeatures.variable_names_to_labels["mtop"], 1, (150., 200.)),
    "ctstar_top" : (cfeatures.variable_names_to_labels["top_ctstar"], 2, (-1., 1.)),
    "ctstar_abs_top" : (cfeatures.variable_names_to_labels["top_ctstar_abs"], 1, (0., 1.)),
        # tbar
    "pt_tbar" : (cfeatures.variable_names_to_labels["pt_tbar"], 2, (0., 500.)),
    "eta_tbar" : (cfeatures.variable_names_to_labels["eta_tbar"], 2, (-4., 4.)),
    "phi_tbar" : (cfeatures.variable_names_to_labels["phi_tbar"], 1, (-3.2, 3.2)),
    "mass_tbar" : (cfeatures.variable_names_to_labels["mtbar"], 1, (150., 200.)),
    "ctstar_tbar" : (cfeatures.variable_names_to_labels["tbar_ctstar"], 2, (-1., 1.)),
    "ctstar_abs_tbar" : (cfeatures.variable_names_to_labels["tbar_ctstar_abs"], 1, (0., 1.)),
        # ttbar
    "pt_tt" : (cfeatures.variable_names_to_labels["pt_tt"], 2, (0., 500.)),
    "eta_tt" : (cfeatures.variable_names_to_labels["eta_tt"], 2, (-4., 4.)),
    "phi_tt" : (cfeatures.variable_names_to_labels["phi_tt"], 1, (-3.2, 3.2)),
    #"mtt" : (cfeatures.variable_names_to_labels["mtt"], final_binning.mtt_binning, (final_binning.mtt_binning[0], final_binning.mtt_binning[-1])),
    "mtt" : (cfeatures.variable_names_to_labels["mtt"], 10, (200., 2000.)),
    "mtt_vs_top_ctstar" : (cfeatures.variable_names_to_labels["mtt"], 1, None),
    "mtt_vs_top_ctstar_abs" : (cfeatures.variable_names_to_labels["mtt"], 1, None),
}

mtt_vs_top_ctstar_abs_binning = (
    final_binning.mtt_binning,
    final_binning.ctstar_abs_binning
)

mtt_vs_top_ctstar_binning = (
    final_binning.mtt_binning,
    np.array([-1.0, -0.6, -0.2, 0.2, 0.6, 1.0])
)

# mtt x ctstar labels
ctstar_binlabels = [r"%s $\in [%s, %s)$" % (cfeatures.variable_names_to_labels["top_ctstar"], mtt_vs_top_ctstar_binning[1][bin], mtt_vs_top_ctstar_binning[1][bin+1]) for bin in range(len(mtt_vs_top_ctstar_binning[1])-1)]
ctstar_binlabels[-1] = r"%s $\in [%s, %s]$" % (cfeatures.variable_names_to_labels["top_ctstar"], mtt_vs_top_ctstar_binning[1][-2], mtt_vs_top_ctstar_binning[1][-1])
ctstar_bin_locs = np.linspace((len(mtt_vs_top_ctstar_binning[0])-1)/2, (len(mtt_vs_top_ctstar_binning[0])-1)*(len(mtt_vs_top_ctstar_binning[1])-1) - (len(mtt_vs_top_ctstar_binning[0])-1)/2, len(mtt_vs_top_ctstar_binning[1])-1)
mtt_ct_vals_to_plot = np.array([400, 600, 1000])
mtt_ct_tiled_labels = np.tile(mtt_vs_top_ctstar_binning[0][:-1], mtt_vs_top_ctstar_binning[1].size-1)
mtt_ct_bin_inds_to_plot = np.where(np.in1d(mtt_ct_tiled_labels, mtt_ct_vals_to_plot))[0]
mtt_ct_bins_to_plot = np.tile(mtt_ct_vals_to_plot, mtt_vs_top_ctstar_binning[1].size-1)

# mtt x ctstar abs labels
ctstar_abs_binlabels = [r"%s $\in [%s, %s)$" % (cfeatures.variable_names_to_labels["top_ctstar_abs"], mtt_vs_top_ctstar_abs_binning[1][bin], mtt_vs_top_ctstar_abs_binning[1][bin+1]) for bin in range(len(mtt_vs_top_ctstar_abs_binning[1])-1)]
ctstar_abs_binlabels[-1] = r"%s $\in [%s, %s]$" % (cfeatures.variable_names_to_labels["top_ctstar_abs"], mtt_vs_top_ctstar_abs_binning[1][-2], mtt_vs_top_ctstar_abs_binning[1][-1])
ctstar_abs_bin_locs = np.linspace((len(mtt_vs_top_ctstar_abs_binning[0])-1)/2, (len(mtt_vs_top_ctstar_abs_binning[0])-1)*(len(mtt_vs_top_ctstar_abs_binning[1])-1) - (len(mtt_vs_top_ctstar_abs_binning[0])-1)/2, len(mtt_vs_top_ctstar_abs_binning[1])-1)
mtt_ctabs_vals_to_plot = np.array([400, 600, 1000])
mtt_ctabs_tiled_labels = np.tile(mtt_vs_top_ctstar_abs_binning[0][:-1], mtt_vs_top_ctstar_abs_binning[1].size-1)
mtt_ctabs_bin_inds_to_plot = np.where(np.in1d(mtt_ctabs_tiled_labels, mtt_ctabs_vals_to_plot))[0]
mtt_ctabs_bins_to_plot = np.tile(mtt_ctabs_vals_to_plot, mtt_vs_top_ctstar_abs_binning[1].size-1)

#set_trace()
    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())[args.year]
lumi_to_use = (data_lumi_year["Muons"]+data_lumi_year["Electrons"])/2000.

lumi_correction_dict = load(os.path.join(proj_dir, "Corrections", jobid, "MC_LumiWeights.coffea"))[args.year]
avg_lumi_scale_dict = {tt : (lumi_correction_dict["Muons"][tt]+lumi_correction_dict["Electrons"][tt])/2 for tt in ["ttJetsSL", "ttJetsHad", "ttJetsDiLep"]}

process_axis = hist.Cat("process", "Process")
reweighting_axis = hist.Cat("rewt", "Reweighting Type")

for hname in variables.keys():
    if hname not in hdict.keys():
        raise ValueError(f"{hname} not found in file")
    histo = hdict[hname].copy()
    histo.scale(avg_lumi_scale_dict, axis="dataset")
    histo = histo.integrate("dataset")

    is2d = histo.dense_dim() == 2
    axes_to_sum = (histo.dense_axes()[0].name,) if not is2d else (histo.dense_axes()[0].name, histo.dense_axes()[1].name)

    #set_trace()
    if is2d:
        if hname == "mtt_vs_top_ctstar":
            xrebinning, yrebinning = mtt_vs_top_ctstar_binning
        if hname == "mtt_vs_top_ctstar_abs":
            xrebinning, yrebinning = mtt_vs_top_ctstar_abs_binning

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
        histo = histo.rebin(yaxis_name, new_ybins)

            # vars only for linearized mtt ctstar hist
        mtt_binning, ctstar_binning = histo.axis(axes_to_sum[0]).edges(), histo.axis(axes_to_sum[1]).edges()
        nbins = (len(mtt_binning)-1)*(len(ctstar_binning)-1)
        vlines = [nbins*ybin/(len(ctstar_binning)-1) for ybin in range(1, len(ctstar_binning)-1)]

            # easier to rename sparse axis than change linearize()
        tmp_histo = histo.copy()
        tmp_histo = tmp_histo.group(histo.sparse_axes()[0].name, process_axis, {key[0]:key[0] for key in histo.values().keys()})
        hline = Plotter.linearize_hist(tmp_histo)
            # revert sparse axis name to original
        hline = hline.group(hline.sparse_axes()[0].name, reweighting_axis, {key[0]:key[0] for key in histo.values().keys()})
        histo = hline

    else:
        xtitle, rebinning, x_lims = variables[hname]
        #set_trace()
        if isinstance(rebinning, np.ndarray):
            new_xbins = hist.Bin(axes_to_sum[0], axes_to_sum[0], rebinning)
        elif isinstance(rebinning, float) or isinstance(rebinning, int):
            new_xbins = rebinning
        histo = histo.rebin(*axes_to_sum, new_xbins)

    #set_trace()
        # yield and ratio plots
    fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True, figsize=(15.0, 10.0)) if is2d else plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
    #fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
    fig.subplots_adjust(hspace=.07)

        # just ratio plots
    fig_ratio, ax_ratio = plt.subplots(figsize=(15.0, 10.0)) if is2d else plt.subplots()
    #fig_ratio, ax_ratio = plt.subplots()
    fig_ratio.subplots_adjust(hspace=.07)

    for rewt_key in sorted(histo.values().keys()):
        hep.plot.histplot(histo[rewt_key].integrate("rewt").values()[()], histo.dense_axes()[0].edges(), ax=ax, histtype="step", **rewt_style_dict[rewt_key[0]]) # nosys template
            ## plot ratios
        if rewt_key == ("Nominal",): continue
                # orig
        ratio_masked_vals, ratio_bins = Plotter.get_ratio_arrays(num_vals=histo.values()[rewt_key], denom_vals=histo.values()[("Nominal",)], input_bins=histo.dense_axes()[0].edges())
        rax.step(ratio_bins, ratio_masked_vals, where="post", **rewt_style_dict[rewt_key[0]])

                # ratio
        ax_ratio.step(ratio_bins, ratio_masked_vals, where="post", **rewt_style_dict[rewt_key[0]])

    #set_trace()
    # plot yields
        # format axes
    ax.autoscale()
    ax.set_xlim((0, nbins) if is2d else x_lims)
    ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.2)
    ax.set_xlabel(None)
    ax.set_ylabel("Events")
    
    rax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
    rax.set_ylabel("Ratio to Powheg")
    #rax.set_ylabel("W(var)/Nominal")
    rax.set_xlabel(xtitle)
    
    ax_ratio.autoscale()
    ax_ratio.set_xlim((0, nbins) if is2d else x_lims)
    ax_ratio.set_ylim(ax_ratio.get_ylim()[0], ax_ratio.get_ylim()[1])
    #ax_ratio.set_ylim(ax_ratio.get_ylim()[0], ax_ratio.get_ylim()[1]*1.05)
    ax_ratio.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
    ax_ratio.set_ylabel("W(var)/Nominal")
    ax_ratio.set_xlabel(xtitle)
    
    ## set plotting styles
        ## set legend and corresponding colors
    #set_trace()
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles,labels, loc="upper right")
    ax.text(
        0.02, 0.86,
        "$t\\bar{t}$\nparton level",
        fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
    )

    handles, labels = ax_ratio.get_legend_handles_labels()
    ax_ratio.legend(handles,labels, loc="upper right")
    ax_ratio.text(
        0.02, 0.90,
        #0.02, 0.86,
        "$t\\bar{t}$\nparton level",
        fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax_ratio.transAxes
    )
    if is2d:
        # plot vertical lines for mtt vs ctstar
        for vline in vlines:
                # orig
            ax.axvline(vline, color="k", linestyle="--")
            rax.axvline(vline, color="k", linestyle="--")
            ax_ratio.axvline(vline, color="k", linestyle="--")

        ## ctstar abs labels
        if hname == "mtt_vs_top_ctstar_abs":
            for idx, label in enumerate(ctstar_abs_binlabels):
                ax.annotate(label, xy=(ctstar_abs_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                    xytext=(0, 250), textcoords="offset points", va="top", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)
                ax_ratio.annotate(label, xy=(ctstar_abs_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                    xytext=(0, 30), textcoords="offset points", va="top", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)
            #set_trace()
            rax.set_xticks(mtt_ctabs_bin_inds_to_plot)
            rax.set_xticklabels(mtt_ctabs_bins_to_plot)
            ax_ratio.set_xticks(mtt_ctabs_bin_inds_to_plot)
            ax_ratio.set_xticklabels(mtt_ctabs_bins_to_plot)

        ## ctstar labels
        if hname == "mtt_vs_top_ctstar":
            for idx, label in enumerate(ctstar_binlabels):
                ax.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                    xytext=(0, 250), textcoords="offset points", va="top", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)
                ax_ratio.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                    xytext=(0, 30), textcoords="offset points", va="top", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)
            #set_trace()
            rax.set_xticks(mtt_ct_bin_inds_to_plot)
            rax.set_xticklabels(mtt_ct_bins_to_plot)
            ax_ratio.set_xticks(mtt_ct_bin_inds_to_plot)
            ax_ratio.set_xticklabels(mtt_ct_bins_to_plot)


    hep.cms.label(ax=ax, data=False, label="Preliminary", year=args.year, lumi=round(lumi_to_use, 1))    

    #figname = os.path.join(outdir, f"TTbar_RewtComp_{hname}_ANbinning")
    figname = os.path.join(outdir, f"TTbar_RewtComp_{hname}")
    fig.savefig(figname)
    plt.close(fig)
    #set_trace()
    print(f"{figname} written")

    hep.cms.label(ax=ax_ratio, data=False, label="Preliminary", year=args.year, lumi=round(lumi_to_use, 1))    
    #ratio_figname = os.path.join(outdir, f"TTbar_RewtComp_{hname}_ratio_ANbinning")
    ratio_figname = os.path.join(outdir, f"TTbar_RewtComp_{hname}_ratio")
    fig_ratio.savefig(ratio_figname)
    plt.close(fig_ratio)
    #set_trace()
    print(f"{ratio_figname} written")

toc = time.time()
print("Total time: %.1f" % (toc - tic))
