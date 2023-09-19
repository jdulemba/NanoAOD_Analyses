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
import Utilities.common_features as cfeatures
import Utilities.final_analysis_binning as final_binning

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2018"], help="Specify which year to run over")
parser.add_argument("--comp_sys", action="store_true", help="Make plots comparing all toponium uncs")
#parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="Specify which year to run over")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
base_jobid = os.environ["base_jobid"]
analyzer = "Toponium_GenLevel"
plot_outdir = os.environ["plots_dir"]
eos_dir = os.environ["eos_dir"]

f_ext = "TOT.coffea"
input_dir = os.path.join(eos_dir, "results", base_jobid, analyzer, args.year)
fnames = sorted(["%s/%s" % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)])
hdict = plt_tools.add_coffea_files(fnames) if len(fnames) > 1 else load(fnames[0])
    
outdir = os.path.join(plot_outdir, base_jobid, analyzer, args.year)#, args.comp)
if not os.path.isdir(outdir):
    os.makedirs(outdir)
    

sys_style_dict = {
    "nosys"             : {"color" : "k", "label" : "Nominal"},
    #"BindingEnergyUp"   : {"color" : "#e41a1c", "label" : "Binding Energy Up"}, # red
    #"BindingEnergyDown" : {"color" : "#377eb8", "label" : "Binding Energy Down"}, #blue
    #"TopMassUp"         : {"color" : "#ff7f00", "label" : "Top Mass Up"}, # orange
    #"TopMassDown"       : {"color" : "#4daf4a", "label" : "Top Mass Down"}, # green
}
if args.comp_sys:
    sys_style_dict.update({
        "BindingEnergyUp"   : {"color" : "#e41a1c", "label" : "Binding Energy Up"}, # red
        "BindingEnergyDown" : {"color" : "#377eb8", "label" : "Binding Energy Down"}, #blue
        "TopMassUp"         : {"color" : "#ff7f00", "label" : "Top Mass Up"}, # orange
        "TopMassDown"       : {"color" : "#4daf4a", "label" : "Top Mass Down"}, # green
    })

var_names = cfeatures.variable_names_to_labels
variables = {
        # top
    "pt_top" : (var_names["pt_top"], 2, (0., 500.)),
    "eta_top" : (var_names["eta_top"], 2, (-4., 4.)),
    "mtop" : (var_names["mtop"], 1, (150., 200.)),
    "top_ctstar" : (var_names["top_ctstar"], 2, (-1., 1.)),
    "top_ctstar_abs" : (var_names["top_ctstar_abs"], 1, (0., 1.)),
        # ttbar
    "pt_tt" : (var_names["pt_tt"], 2, (0., 500.)),
    "eta_tt" : (var_names["eta_tt"], 2, (-4., 4.)),
    "mtt" : (var_names["mtt"], 1, (320., 360.)),
    ##"mtt" : (var_names["mtt"], 1, (300., 2000.)),
    "mtt_vs_top_ctstar_abs" : (var_names["mtt"], var_names["ctstar_abs"], 1, None),
    ##"mttANbinning" : (var_names["mtt"], final_binning.mtt_binning, (final_binning.mtt_binning[0], final_binning.mtt_binning[-1])),
}

mtt_vals_to_plot = np.array([400, 600, 1000])
mtt_vs_top_ctstar_abs_binning = (
    final_binning.mtt_binning,
    final_binning.ctstar_abs_binning
)

ctstar_abs_binlabels = [r"%s $\in [%s, %s)$" % (var_names["ctstar_abs"], mtt_vs_top_ctstar_abs_binning[1][bin], mtt_vs_top_ctstar_abs_binning[1][bin+1]) for bin in range(len(mtt_vs_top_ctstar_abs_binning[1])-1)]
ctstar_abs_binlabels[-1] = r"%s $\in [%s, %s]$" % (var_names["ctstar_abs"], mtt_vs_top_ctstar_abs_binning[1][-2], mtt_vs_top_ctstar_abs_binning[1][-1])
ctstar_abs_bin_locs = np.linspace((len(mtt_vs_top_ctstar_abs_binning[0])-1)/2, (len(mtt_vs_top_ctstar_abs_binning[0])-1)*(len(mtt_vs_top_ctstar_abs_binning[1])-1) - (len(mtt_vs_top_ctstar_abs_binning[0])-1)/2, len(mtt_vs_top_ctstar_abs_binning[1])-1)
mtt_abs_tiled_labels = np.tile(mtt_vs_top_ctstar_abs_binning[0][:-1], mtt_vs_top_ctstar_abs_binning[1].size-1)
mtt_abs_bin_inds_to_plot = np.where(np.in1d(mtt_abs_tiled_labels, mtt_vals_to_plot))[0]
mtt_abs_bins_to_plot = np.tile(mtt_vals_to_plot, mtt_vs_top_ctstar_abs_binning[1].size-1)


    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())[args.year]
lumi_to_use = (data_lumi_year["Muons"]+data_lumi_year["Electrons"])/2000.

#set_trace()
# scale {(ToponiumSL, nosys) : ToponiumSL_nosys}
lumi_correction = load(os.path.join(proj_dir, "Corrections", jobid, "MC_LumiWeights.coffea"))[args.year]
avg_lumi_corr = {key : (lumi_correction["Muons"]["_".join(key)] + lumi_correction["Electrons"]["_".join(key)])/2 for key in sorted(set([(key[0], key[1]) for key in hdict["mtt"].values().keys()]))}
    ## make groups based on process
process_axis = hist.Cat("process", "Process", sorting="placement")
process_cat = "dataset"
#etaT_samples = ["ToponiumSL", "ToponiumDiLep"]
#process_groups = {"Toponium" : ["ToponiumSL", "ToponiumDiLep"]}

#set_trace()
for hname in variables.keys():
    #if hname not in hdict.keys():
    #    raise ValueError(f"{hname} not found in file")
    histo = hdict["mtt"].copy() if hname == "mttANbinning" else hdict[hname].copy()
    histo.scale(avg_lumi_corr, axis=("dataset", "sys"))
    histo = (histo.integrate("dataset")).copy()

    is2d = histo.dense_dim() == 2
    axes_to_sum = (histo.dense_axes()[0].name,) if not is2d else (histo.dense_axes()[0].name, histo.dense_axes()[1].name)

    #set_trace()
    if is2d:
        xtitle, ytitle, rebinning, x_lims = variables[hname]
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

        #set_trace()
            # easier to rename sparse axis than change linearize()
        tmp_histo = histo.copy()
        tmp_histo = tmp_histo.group(histo.sparse_axes()[0].name, process_axis, {key[0]:key[0] for key in histo.values().keys()})
        histo = Plotter.linearize_hist(tmp_histo).copy()

    else:
        xtitle, rebinning, x_lims = variables[hname]
        #set_trace()
        if isinstance(rebinning, np.ndarray):
            new_xbins = hist.Bin(axes_to_sum[0], axes_to_sum[0], rebinning)
        elif isinstance(rebinning, float) or isinstance(rebinning, int):
            new_xbins = rebinning
        histo = histo.rebin(*axes_to_sum, new_xbins)

    #set_trace()
    if args.comp_sys:
        fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True, figsize=(15.0, 10.0)) if is2d else plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
    else:
        fig, ax = plt.subplots(figsize=(15.0, 10.0)) if is2d else plt.subplots()
        rax = None
    fig.subplots_adjust(hspace=.07)

    for sysname in sys_style_dict.keys():
        hep.plot.histplot(histo.values()[(sysname,)], histo.dense_axes()[0].edges(), ax=ax, histtype="step", **sys_style_dict[sysname])
        if args.comp_sys:
            if sysname == "nosys": continue
                # plot ratio to nominal powheg
            ratio_masked_vals, ratio_masked_bins = Plotter.get_ratio_arrays(num_vals=histo.values()[(sysname,)], denom_vals=histo.values()[("nosys",)], input_bins=histo.dense_axes()[0].edges())
            rax.step(ratio_masked_bins, ratio_masked_vals, where='post', **sys_style_dict[sysname])

    ## set plotting styles
        ## set legend and corresponding colors
    if rax: ax.legend(loc="upper right", ncol=2)
        # format axes
    #ax.set_yscale("log")
    ax.autoscale()
    ax.set_xlim((0, nbins) if is2d else x_lims)
    ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.3 if rax else ax.get_ylim()[1]*1.1)
    ax.set_xlabel(None if rax else xtitle)
    ax.set_ylabel("Events")

    if rax:
        rax.autoscale()
        rax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
        rax.set_ylabel("Ratio to Nominal")
        rax.set_xlabel(xtitle)
        rax.set_xlim((0, nbins) if is2d else x_lims)
    

    if is2d:
        [ax.axvline(vline, color="k", linestyle="--") for vline in vlines]
        if rax: [rax.axvline(vline, color="k", linestyle="--") for vline in vlines]

        y_binlabels = ctstar_abs_binlabels
        y_bin_locs = ctstar_abs_bin_locs
        x_inds_to_plot = mtt_abs_bin_inds_to_plot
        x_bins_to_plot = mtt_abs_bins_to_plot

        ## plot y labels
        for idx, label in enumerate(y_binlabels):
            ax.annotate(label, xy=(y_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                xytext=(0, 250), textcoords="offset points", va="top", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)

        if rax:
            rax.set_xticks(x_inds_to_plot)
            rax.set_xticklabels(x_bins_to_plot)
        else:
            ax.set_xticks(x_inds_to_plot)
            ax.set_xticklabels(x_bins_to_plot)

    ax.text(
        0.02, 0.86 if rax else 0.90, "$\\eta_{t}$\nparton level",
        fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
    )
    hep.cms.label(ax=ax, data=False, label="Preliminary", year=cfeatures.year_labels[args.year], lumi=round(lumi_to_use, 1))    

    figname = os.path.join(outdir, f"{args.year}_{analyzer}_{hname}_SysComp" if args.comp_sys else f"{args.year}_{analyzer}_{hname}")
    fig.savefig(figname)
    #set_trace()
    print(f"{figname} written")
    plt.close(fig)


toc = time.time()
print("Total time: %.1f" % (toc - tic))
