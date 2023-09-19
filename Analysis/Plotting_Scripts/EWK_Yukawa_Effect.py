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
analyzer = "EWK_Yukawa_Effect"
plot_outdir = os.environ["plots_dir"]
eos_dir = os.environ["eos_dir"]

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2018"], help="Specify which year to run over")
#parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="Specify which year to run over")
args = parser.parse_args()

f_ext = "TOT.coffea"
input_dir = os.path.join(eos_dir, "results", f"{args.year}_{base_jobid}", analyzer)
fnames = sorted(["%s/%s" % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)])
hdict = plt_tools.add_coffea_files(fnames) if len(fnames) > 1 else load(fnames[0])
    
outdir = os.path.join(plot_outdir, base_jobid, analyzer, args.year)
if not os.path.isdir(outdir):
    os.makedirs(outdir)
    
rewt_style_dict = {
    "nosys" : {"label" : "$Y_t = 1.0$", "color" : "k", "linestyle" : "-"},
    "YukawaUp" : {"label" : "$Y_t = 1.1$", "color" : "#ff7f00", "linestyle" : "-"}, ## orange
    "YukawaDown" : {"label" : "$Y_t = 0.9$", "color" : "#984ea3", "linestyle" : "-"}, ## purple
    "NewYukawaUp" : {"label" : "$Y_t = 1.11$", "color" : "r", "linestyle" : "-"}, ## red
    "NewYukawaDown" : {"label" : "$Y_t = 0.88$", "color" : "b", "linestyle" : "-"}, ## blue

}

variables = {
        # top
    "pt_top" : ("$p_{T}$($t$) [GeV]", 2, (0., 500.)),
    "eta_top" : ("$\eta$($t$)", 2, (-4., 4.)),
    "phi_top" : ("$\phi$($t$)", 1, (-3.2, 3.2)),
    "mass_top" : ("m($t$) [GeV]", 1, (150., 200.)),
    "ctstar_top" : ("cos($\\theta^{*}_{t}$)", 2, (-1., 1.)),
    "ctstar_abs_top" : ("|cos($\\theta^{*}_{t}$)|", 1, (0., 1.)),
    "cpTP_top" : ("cpTP($t$)", 2, (-1., 1.)),
        # tbar
    "pt_tbar" : ("$p_{T}$($\\bar{t}$) [GeV]", 2, (0., 500.)),
    "eta_tbar" : ("$\eta$($\\bar{t}$)", 2, (-4., 4.)),
    "phi_tbar" : ("$\phi$($\\bar{t}$)", 1, (-3.2, 3.2)),
    "mass_tbar" : ("m($\\bar{t}$) [GeV]", 1, (150., 200.)),
    "ctstar_tbar" : ("cos($\\theta^{*}_{\\bar{t}}$)", 2, (-1., 1.)),
    "ctstar_abs_tbar" : ("|cos($\\theta^{*}_{\\bar{t}}$)|", 1, (0., 1.)),
    "cpTP_tbar" : ("cpTP($\\bar{t}$)", 2, (-1., 1.)),
        # ttbar
    "pt_tt" : ("$p_{T}$($t\\bar{t}$) [GeV]", 2, (0., 500.)),
    "eta_tt" : ("$\\eta$($t\\bar{t}$)", 2, (-4., 4.)),
    "phi_tt" : ("$\phi$($t\\bar{t}$)", 1, (-3.2, 3.2)),
    "mtt" : ("$m_{t\\bar{t}}$ [GeV]", 10, (350., 4000.)),
    "mtt_vs_top_ctstar" : (cfeatures.variable_names_to_labels["mtt"], cfeatures.variable_names_to_labels["ctstar"], 1, None),
    "mtt_vs_top_ctstar_abs" : (cfeatures.variable_names_to_labels["mtt"], cfeatures.variable_names_to_labels["ctstar_abs"], 1, None),
}

mtt_vals_to_plot = np.array([400, 600, 1000])
mtt_vs_top_ctstar_abs_binning = (
    final_binning.mtt_binning,
    final_binning.ctstar_abs_binning
)

ctstar_abs_binlabels = [r"%s $\in [%s, %s)$" % (cfeatures.variable_names_to_labels["ctstar_abs"], mtt_vs_top_ctstar_abs_binning[1][bin], mtt_vs_top_ctstar_abs_binning[1][bin+1]) for bin in range(len(mtt_vs_top_ctstar_abs_binning[1])-1)]
ctstar_abs_binlabels[-1] = r"%s $\in [%s, %s]$" % (cfeatures.variable_names_to_labels["ctstar_abs"], mtt_vs_top_ctstar_abs_binning[1][-2], mtt_vs_top_ctstar_abs_binning[1][-1])
ctstar_abs_bin_locs = np.linspace((len(mtt_vs_top_ctstar_abs_binning[0])-1)/2, (len(mtt_vs_top_ctstar_abs_binning[0])-1)*(len(mtt_vs_top_ctstar_abs_binning[1])-1) - (len(mtt_vs_top_ctstar_abs_binning[0])-1)/2, len(mtt_vs_top_ctstar_abs_binning[1])-1)
mtt_abs_tiled_labels = np.tile(mtt_vs_top_ctstar_abs_binning[0][:-1], mtt_vs_top_ctstar_abs_binning[1].size-1)
mtt_abs_bin_inds_to_plot = np.where(np.in1d(mtt_abs_tiled_labels, mtt_vals_to_plot))[0]
mtt_abs_bins_to_plot = np.tile(mtt_vals_to_plot, mtt_vs_top_ctstar_abs_binning[1].size-1)

mtt_vs_top_ctstar_binning = (
    final_binning.mtt_binning,
    np.array([-1.0, -0.9, -0.75, -0.6, -0.4, 0.0, 0.4, 0.6, 0.75, 0.9, 1.0])
)

ctstar_binlabels = [r"%s $\in [%s, %s)$" % (cfeatures.variable_names_to_labels["ctstar"], mtt_vs_top_ctstar_binning[1][bin], mtt_vs_top_ctstar_binning[1][bin+1]) for bin in range(len(mtt_vs_top_ctstar_binning[1])-1)]
ctstar_binlabels[-1] = r"%s $\in [%s, %s]$" % (cfeatures.variable_names_to_labels["ctstar"], mtt_vs_top_ctstar_binning[1][-2], mtt_vs_top_ctstar_binning[1][-1])
ctstar_bin_locs = np.linspace((len(mtt_vs_top_ctstar_binning[0])-1)/2, (len(mtt_vs_top_ctstar_binning[0])-1)*(len(mtt_vs_top_ctstar_binning[1])-1) - (len(mtt_vs_top_ctstar_binning[0])-1)/2, len(mtt_vs_top_ctstar_binning[1])-1)
mtt_tiled_labels = np.tile(mtt_vs_top_ctstar_binning[0][:-1], mtt_vs_top_ctstar_binning[1].size-1)
mtt_bin_inds_to_plot = np.where(np.in1d(mtt_tiled_labels, mtt_vals_to_plot))[0]
mtt_bins_to_plot = np.tile(mtt_vals_to_plot, mtt_vs_top_ctstar_binning[1].size-1)

    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())[args.year]
lumi_to_use = (data_lumi_year["Muons"]+data_lumi_year["Electrons"])/2000.

lumi_correction_dict = load(os.path.join(proj_dir, "Corrections", jobid, "MC_LumiWeights.coffea"))[args.year]
avg_lumi_scale_dict = {tt : (lumi_correction_dict["Muons"][tt]+lumi_correction_dict["Electrons"][tt])/2 for tt in ["ttJetsSL", "ttJetsHad", "ttJetsDiLep"]}

process_axis = hist.Cat("process", "Process")
reweighting_axis = hist.Cat("rewt", "Reweighting Type")

#set_trace()
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
        #set_trace()
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

        xtitle, ytitle, rebinning, x_lims = variables[hname]

    else:
        xtitle, rebinning, x_lims = variables[hname]
    if rebinning != 1:
        histo = histo.rebin(*axes_to_sum, rebinning)

    #set_trace()
        # get yields
    nosys_vals = histo.values()[("nosys"),]
    #YukawaUp_vals = histo.values()[("YukawaUp",)]
    #YukawaDown_vals = histo.values()[("YukawaDown",)]
    NewYukawaUp_vals = histo.values()[("NewYukawaUp",)]
    NewYukawaDown_vals = histo.values()[("NewYukawaDown",)]

    #set_trace()
    fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True, figsize=(15.0, 10.0)) if is2d else plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
    fig.subplots_adjust(hspace=.07)

        # plot yields
    hep.plot.histplot(nosys_vals, histo.dense_axes()[0].edges(), ax=ax, histtype="step", **rewt_style_dict["nosys"])
    #hep.plot.histplot(YukawaUp_vals, histo.dense_axes()[0].edges(), ax=ax, histtype="step", **rewt_style_dict["YukawaUp"])
    #hep.plot.histplot(YukawaDown_vals, histo.dense_axes()[0].edges(), ax=ax, histtype="step", **rewt_style_dict["YukawaDown"])
    hep.plot.histplot(NewYukawaUp_vals, histo.dense_axes()[0].edges(), ax=ax, histtype="step", **rewt_style_dict["NewYukawaUp"])
    hep.plot.histplot(NewYukawaDown_vals, histo.dense_axes()[0].edges(), ax=ax, histtype="step", **rewt_style_dict["NewYukawaDown"])

        # plot relative deviations
    #up_ratio_masked_vals, up_ratio_masked_bins = Plotter.get_ratio_arrays(num_vals=YukawaUp_vals, denom_vals=nosys_vals, input_bins=histo.dense_axes()[0].edges())
    #rax.step(up_ratio_masked_bins, up_ratio_masked_vals, where='post', **{"color" : rewt_style_dict["YukawaUp"]["color"], "linestyle" : "-"})
    #dw_ratio_masked_vals, dw_ratio_masked_bins = Plotter.get_ratio_arrays(num_vals=YukawaDown_vals, denom_vals=nosys_vals, input_bins=histo.dense_axes()[0].edges())
    #rax.step(dw_ratio_masked_bins, dw_ratio_masked_vals, where='post', **{"color" : rewt_style_dict["YukawaDown"]["color"], "linestyle" : "-"})
    up_ratio_masked_vals, up_ratio_masked_bins = Plotter.get_ratio_arrays(num_vals=NewYukawaUp_vals, denom_vals=nosys_vals, input_bins=histo.dense_axes()[0].edges())
    rax.step(up_ratio_masked_bins, up_ratio_masked_vals, where='post', **{"color" : rewt_style_dict["NewYukawaUp"]["color"], "linestyle" : "-"})
    dw_ratio_masked_vals, dw_ratio_masked_bins = Plotter.get_ratio_arrays(num_vals=NewYukawaDown_vals, denom_vals=nosys_vals, input_bins=histo.dense_axes()[0].edges())
    rax.step(dw_ratio_masked_bins, dw_ratio_masked_vals, where='post', **{"color" : rewt_style_dict["NewYukawaDown"]["color"], "linestyle" : "-"})

    #set_trace()
    # plot yields
        # format axes
    ax.set_yscale("log")
    ax.autoscale()
    ax.set_xlim((0, nbins) if is2d else x_lims)
    ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.2)
    ax.set_xlabel(None)
    ax.set_ylabel("Events")

    rax.autoscale()
    rax.axhline(1.0, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
    rax.set_ylabel("Ratio to $Y_t = 1.0$")
    #rax.set_ylabel("Ratio to Nominal")
    rax.set_xlabel(xtitle)
    rax.set_xlim((0, nbins) if is2d else x_lims)
    if hname == "mtt": rax.set_ylim(-0.3, 0.3)
    
    ## set plotting styles
        ## set legend and corresponding colors
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles,labels, loc="upper right")

    if is2d:
        # plot vertical lines for mtt vs ctstar
        [ax.axvline(vline, color="k", linestyle="--") for vline in vlines]
        [rax.axvline(vline, color="k", linestyle="--") for vline in vlines]

        y_binlabels = ctstar_abs_binlabels if hname == "mtt_vs_top_ctstar_abs" else ctstar_binlabels
        y_bin_locs = ctstar_abs_bin_locs if hname == "mtt_vs_top_ctstar_abs" else ctstar_bin_locs
        x_inds_to_plot = mtt_abs_bin_inds_to_plot if hname == "mtt_vs_top_ctstar_abs" else mtt_bin_inds_to_plot
        x_bins_to_plot = mtt_abs_bins_to_plot if hname == "mtt_vs_top_ctstar_abs" else mtt_bins_to_plot

        ## plot y labels
        for idx, label in enumerate(y_binlabels):
            ax.annotate(label, xy=(y_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                xytext=(0, 30), textcoords="offset points", va="top", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)

        rax.set_xticks(x_inds_to_plot)
        rax.set_xticklabels(x_bins_to_plot)


    ax.text(
        0.05, 0.85, "$t\\bar{t}$\nparton level",
        fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
    )
    hep.cms.label(ax=ax, data=False, label="Preliminary", year=cfeatures.year_labels[args.year], lumi=round(lumi_to_use, 1))    

    #set_trace()
    figname = os.path.join(outdir, f"{args.year}_TTbar_EWKYukawa_Effect_{hname}")
    fig.savefig(figname)
    #set_trace()
    print(f"{figname} written")
    plt.close(fig)

toc = time.time()
print("Total time: %.1f" % (toc - tic))
