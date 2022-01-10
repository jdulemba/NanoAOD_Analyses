#!/usr/bin/env python

# matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend("agg")
from matplotlib import rcParams
rcParams["font.size"] = 20
rcParams["savefig.format"] = "png"
rcParams["savefig.bbox"] = "tight"

import Utilities.Plotter as Plotter
from coffea import hist
from coffea.util import load
from pdb import set_trace
import os
import Utilities.prettyjson as prettyjson
import Utilities.plot_tools as plt_tools
import numpy as np

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
base_jobid = os.environ["base_jobid"]
analyzer = "check_ewk_uncs"

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="Specify which year to run over")
parser.add_argument("plot", choices=["Otto2Afiq"], help="Specify which distribution to compare yt reweighting values to")
parser.add_argument("comp", choices=["Nominal", "yt1"], help="Specify which distribution to compare yt reweighting values to")
args = parser.parse_args()

f_ext = "TOT.coffea"
input_dir = os.path.join(proj_dir, "results", f"{args.year}_{base_jobid}", analyzer)
fnames = sorted(["%s/%s" % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)])
hdict = plt_tools.add_coffea_files(fnames) if len(fnames) > 1 else load(fnames[0])
    
outdir = os.path.join(proj_dir, "plots", base_jobid, analyzer, args.year)
if not os.path.isdir(outdir):
    os.makedirs(outdir)
    
rewt_style_dict = {
    "Nominal" : {"label" : "SM $t\\bar{t}$", "color" : "k", "linestyle" : "-"},
    "yt_0" : {"label" : "$y^{t}/y^{t}_{SM}$ = 0", "color" : "#e41a1c", "linestyle" : "-"}, ## red
    "yt_0.5" : {"label" : "$y^{t}/y^{t}_{SM}$ = 0.5", "color" : "#4daf4a", "linestyle" : "-"}, ## green
    "yt_1" : {"label" : "$y^{t}/y^{t}_{SM}$ = 1", "color" : "#377eb8", "linestyle" : "-"}, ## blue
    "yt_1.5" : {"label" : "$y^{t}/y^{t}_{SM}$ = 1.5", "color" : "#ff7f00", "linestyle" : "-"}, ## orange
    "yt_2" : {"label" : "$y^{t}/y^{t}_{SM}$ = 2", "color" : "#984ea3", "linestyle" : "-"}, ## purple
    #"0" : {"label" : "$y^{t}/y^{t}_{SM}$ = 0", "color" : "#e41a1c", "linestyle" : "-"}, ## red
    "0.5" : {"label" : "$y^{t}/y^{t}_{SM}$ = 0.5", "color" : "#4daf4a", "linestyle" : "-"}, ## green
    "0.9" : {"label" : "$y^{t}/y^{t}_{SM}$ = 0.9", "color" : "#377eb8", "linestyle" : "-"}, ## blue
    "1.0" : {"label" : "$y^{t}/y^{t}_{SM}$ = 1", "color" : "k", "linestyle" : "-"}, ## black
    "1.1" : {"label" : "$y^{t}/y^{t}_{SM}$ = 1.1", "color" : "#984ea3", "linestyle" : "-"}, ## purple
    "1.5" : {"label" : "$y^{t}/y^{t}_{SM}$ = 1.5", "color" : "#ff7f00", "linestyle" : "-"}, ## orange
    #"2.0" : {"label" : "$y^{t}/y^{t}_{SM}$ = 2", "color" : "#984ea3", "linestyle" : "-"}, ## purple
    #"3.0" : {"label" : "$y^{t}/y^{t}_{SM}$ = 3", "color" : "#984ea3", "linestyle" : "-"}, ## purple
    #"4.0" : {"label" : "$y^{t}/y^{t}_{SM}$ = 4", "color" : "#984ea3", "linestyle" : "-"}, ## purple
    "Rebinned_0.5" : {"label" : "$y^{t}/y^{t}_{SM}$ = 0.5", "color" : "#4daf4a", "linestyle" : "--"}, ## green
    "Rebinned_0.9" : {"label" : "$y^{t}/y^{t}_{SM}$ = 0.9", "color" : "#377eb8", "linestyle" : "--"}, ## blue
    "Rebinned_1.0" : {"label" : "$y^{t}/y^{t}_{SM}$ = 1", "color" : "k", "linestyle" : "--"}, ## black
    "Rebinned_1.1" : {"label" : "$y^{t}/y^{t}_{SM}$ = 1.1", "color" : "#984ea3", "linestyle" : "--"}, ## purple
    "Rebinned_1.5" : {"label" : "$y^{t}/y^{t}_{SM}$ = 1.5", "color" : "#ff7f00", "linestyle" : "--"}, ## orange
}
compare_otto_to_afiq_dict = {
    "Nominal" : {"label" : "SM $t\\bar{t}$", "color" : "k", "linestyle" : "-"},
    "yt_0.5" : {"label" : "$y^{t}/y^{t}_{SM}$ = 0.5", "color" : "#4daf4a", "linestyle" : "-"}, ## green
    "yt_1" : {"label" : "$y^{t}/y^{t}_{SM}$ = 1", "color" : "#377eb8", "linestyle" : "-"}, ## blue
    "yt_1.5" : {"label" : "$y^{t}/y^{t}_{SM}$ = 1.5", "color" : "#e41a1c", "linestyle" : "-"}, ## red
    "0.5" : {"label" : "Otto 0.5", "color" : "#4daf4a", "linestyle" : "--"}, ## green
    "1.0" : {"label" : "Otto 1", "color" : "#377eb8", "linestyle" : "--"}, ## blue
    "1.5" : {"label" : "Otto 1.5", "color" : "#e41a1c", "linestyle" : "--"}, ## red
    "Rebinned_0.5" : {"label" : "Otto coarse 0.5", "color" : "#4daf4a", "linestyle" : "dotted"}, ## green
    "Rebinned_1.0" : {"label" : "Otto coarse 1", "color" : "#377eb8", "linestyle" : "dotted"}, ## blue
    "Rebinned_1.5" : {"label" : "Otto coarse 1.5", "color" : "#e41a1c", "linestyle" : "dotted"}, ## red
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
    "mtt" : ("$m_{t\\bar{t}}$ [GeV]", 10, (200., 4000.)),
    "mtt_vs_top_ctstar" : ("$m_{t\\bar{t}}$ $\otimes$ cos($\\theta^{*}_{t}$)", 1, None),
    "mtt_vs_top_ctstar_abs" : ("$m_{t\\bar{t}}$ $\otimes$ |cos($\\theta^{*}_{t}$)|", 1, None),
}

mtt_vs_top_ctstar_abs_binning = (
    np.array([360.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 625.0, 650.0, 675.0,
        700.0, 730.0, 760.0, 800.0, 850.0, 900.0, 950., 1000.0, 1050., 1100.0, 1150., 1200., 1300., 1500., 2000.]),
    np.array([0.0, 0.4, 0.6, 0.75, 0.9, 1.0])
)

mtt_vs_top_ctstar_binning = (
    np.array([360.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 625.0, 650.0, 675.0,
        700.0, 730.0, 760.0, 800.0, 850.0, 900.0, 950., 1000.0, 1050., 1100.0, 1150., 1200., 1300., 1500., 2000.]),
    np.array([-1.0, -0.9, -0.75, -0.6, -0.4, 0.0, 0.4, 0.6, 0.75, 0.9, 1.0])
)

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

    xtitle, rebinning, x_lims = variables[hname]
    if rebinning != 1:
        histo = histo.rebin(*axes_to_sum, rebinning)

        # orig
    fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
    fig.subplots_adjust(hspace=.07)

    #set_trace()
    if args.plot == "Otto2Afiq":
        for rewt_key in compare_otto_to_afiq_dict.keys():
                # plot yields    
            #if (rewt_key != ("yt_1",)) or (rewt_key == ("Nominal",)): continue
            hep.plot.histplot(histo[(rewt_key),].integrate("rewt").values()[()], histo.dense_axes()[0].edges(), ax=ax, histtype="step", **compare_otto_to_afiq_dict[rewt_key]) 

                ## plot ratios
            #if rewt_key != ("yt_1",): continue
            #if rewt_key == "Nominal": continue
            if rewt_key == args.comp: continue
                    # orig
            ratio_masked_vals, ratio_bins = Plotter.get_ratio_arrays(num_vals=histo.values()[(rewt_key),], denom_vals=histo.values()[(args.comp,)], input_bins=histo.dense_axes()[0].edges())
            rax.step(ratio_bins, ratio_masked_vals, where="post", **compare_otto_to_afiq_dict[rewt_key])

    #set_trace()
    # plot yields
        # format axes
    ax.autoscale()
    ax.set_xlim((0, nbins) if is2d else x_lims)
    ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.2)
    ax.set_xlabel(None)
    ax.set_ylabel("Events")
    
    rax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
    rax.set_ylabel("Yukawa Weight/$t\\bar{t}_{SM}$" if args.comp == "Nominal" else "Yukawa Weight/$y^{t}$ = 1")
    rax.set_xlabel(xtitle)
    
    ## set plotting styles
        ## set legend and corresponding colors
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles,labels, loc="upper right")

    if is2d:
        # plot vertical lines for mtt vs ctstar
        for vline in vlines:
                # orig
            ax.axvline(vline, color="k", linestyle="--")
            rax.axvline(vline, color="k", linestyle="--")

    hep.cms.label(ax=ax, data=False, year=args.year, lumi=round(lumi_to_use, 1))    

    #set_trace()
    figname = os.path.join(outdir, f"{args.year}_TTbar_{args.plot}_EWKRewtComp_to%s_{hname}" % (args.comp).upper())
    fig.savefig(figname)
    #set_trace()
    print(f"{figname} written")
    plt.close()
