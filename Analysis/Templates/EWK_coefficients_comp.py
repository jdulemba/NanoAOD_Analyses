import time
tic = time.time()

import os
import uproot
import argparse
import numpy as np
from collections import defaultdict
from pdb import set_trace

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

parser = argparse.ArgumentParser()
parser.add_argument("channel", choices=["lj", "ll"], help="Choose which channel to make plots for")
parser.add_argument("comp", choices=["channel", "year"], help="Choose which comparison to make plots for")
args = parser.parse_args()

input_dir = "/eos/user/j/jdulemba/NanoAOD_Analyses/results/Summer20UL_DeepJet/Combine_Templates/dyscales_fix_231020"

ll_rfile = uproot.open(os.path.join(input_dir, "ll", "bkg_ll_3D-33_rate_mtuX_pca_ewk.root"))
lj_rfile = uproot.open(os.path.join(input_dir, "lj", "templates_lj_bkg_rate_mtuX_pca_ewk.root"))

lj_chan_labels = {
    "mu3jets" : ("$\\mu$/3 jets", "r"),
    "mu4pjets" : ("$\\mu$/4+ jets", "g"),
    "e3jets" : ("$e$/3 jets", "k"),
    "e4pjets" : ("$e$/4+ jets", "b"),
}

ll_chan_labels = {
    "em" : ("$e/\\mu$", "r"),
    "ee" : ("$ee$", "b"),
    "mm" : ("$\\mu \\mu$", "k"),
}
#years = ["2018"]
years = ["2016pre", "2016post", "2017", "2018"]
years_labels = {
    "2016pre" : "g",
    "2016post" : "b",
    "2017" : "r",
    "2018" : "k",
}

lj_edges = np.arange(111, dtype=float)
ll_edges = np.arange(181, dtype=float)
#set_trace()

plot_outdir = os.environ["plots_dir"]
jobid = os.environ["jobid"]
outdir = os.path.join(plot_outdir, jobid, "EWK_Coefficients_Comp")
if not os.path.isdir(outdir):
    os.makedirs(outdir)


### make plots comparing channels per year
if args.comp == "channel":
    for year in years:
        pltdir = os.path.join(outdir, year)
        if not os.path.isdir(pltdir):
            os.makedirs(pltdir)
    
                ## lj plots    
        if args.channel == "lj":
            b0_fig, b0_ax = plt.subplots(figsize=(15.0, 10.0))
            b0_fig.subplots_adjust(hspace=.07)
            b1_fig, b1_ax = plt.subplots(figsize=(15.0, 10.0))
            b1_fig.subplots_adjust(hspace=.07)
            b2_fig, b2_ax = plt.subplots(figsize=(15.0, 10.0))
            b2_fig.subplots_adjust(hspace=.07)
        
            for channel, (label, color) in lj_chan_labels.items():
                cat = channel + "_" + year
                
                TT = np.copy(lj_rfile[cat]["TT"].values())
                b0_pos, b0_neg = np.copy(lj_rfile[cat]["EWK_TT_const_pos"].values()), np.copy(lj_rfile[cat]["EWK_TT_const_neg"].values())
                b1_pos, b1_neg = np.copy(lj_rfile[cat]["EWK_TT_lin_pos"].values()), np.copy(lj_rfile[cat]["EWK_TT_lin_neg"].values())
                b2_pos, b2_neg = np.copy(lj_rfile[cat]["EWK_TT_quad_pos"].values()), np.copy(lj_rfile[cat]["EWK_TT_quad_neg"].values())
        
                #set_trace()
            
                Plotter.plot_1D((b0_pos-b0_neg)/TT, lj_edges, xlimits=(lj_edges[0], lj_edges[-1]),
                    color=color, ax=b0_ax, label=label, ylabel="$b_{0}$/$TT_{NNLO}$")
                Plotter.plot_1D((b1_pos-b1_neg)/TT, lj_edges, xlimits=(lj_edges[0], lj_edges[-1]),
                    color=color, ax=b1_ax, label=label, ylabel="$b_{1}$/$TT_{NNLO}$")
                Plotter.plot_1D((b2_pos-b2_neg)/TT, lj_edges, xlimits=(lj_edges[0], lj_edges[-1]),
                    color=color, ax=b2_ax, label=label, ylabel="$b_{2}$/$TT_{NNLO}$")
            #set_trace()
            b0_ax.legend(loc="upper right", title="$b_{0}$ Coefficient", ncol=2)
            b0_ax.axhline(0, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
            b0_ax.autoscale()
            b0_ax.set_ylim(b0_ax.get_ylim()[0], b0_ax.get_ylim()[1]*1.15)
            b0_ax.set_xlim((lj_edges[0], lj_edges[-1]))
            hep.cms.label(ax=b0_ax, data=False, year=year)
            b0_figname = os.path.join(pltdir, "_".join(["lj", "b0", year]))
            b0_fig.savefig(b0_figname)
            print(f"{b0_figname} written")
            plt.close(b0_fig)
        
            b1_ax.legend(loc="upper right", title="$b_{1}$ Coefficient", ncol=2)
            b1_ax.axhline(0, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
            b1_ax.autoscale()
            b1_ax.set_ylim(b1_ax.get_ylim()[0], b1_ax.get_ylim()[1]*1.15)
            b1_ax.set_xlim((lj_edges[0], lj_edges[-1]))
            hep.cms.label(ax=b1_ax, data=False, year=year)
            b1_figname = os.path.join(pltdir, "_".join(["lj", "b1", year]))
            b1_fig.savefig(b1_figname)
            print(f"{b1_figname} written")
            plt.close(b1_fig)
        
            b2_ax.legend(loc="upper right", title="$b_{2}$ Coefficient", ncol=2)
            b2_ax.axhline(0, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
            b2_ax.autoscale()
            b2_ax.set_ylim(b2_ax.get_ylim()[0], b2_ax.get_ylim()[1]*1.15)
            b2_ax.set_xlim((lj_edges[0], lj_edges[-1]))
            hep.cms.label(ax=b2_ax, data=False, year=year)
            b2_figname = os.path.join(pltdir, "_".join(["lj", "b2", year]))
            b2_fig.savefig(b2_figname)
            print(f"{b2_figname} written")
            plt.close(b2_fig)
    
    
            ## ll plots    
        if args.channel == "ll":
            b0_fig, b0_ax = plt.subplots(figsize=(15.0, 10.0))
            b0_fig.subplots_adjust(hspace=.07)
            b1_fig, b1_ax = plt.subplots(figsize=(15.0, 10.0))
            b1_fig.subplots_adjust(hspace=.07)
            b2_fig, b2_ax = plt.subplots(figsize=(15.0, 10.0))
            b2_fig.subplots_adjust(hspace=.07)
    
            for channel, (label, color) in ll_chan_labels.items():
                cat = channel + "_" + year
                
                TT = np.copy(ll_rfile[cat]["TT"].values())
                b0_pos, b0_neg = np.copy(ll_rfile[cat]["EWK_TT_const_pos"].values()), np.copy(ll_rfile[cat]["EWK_TT_const_neg"].values())
                b1_pos, b1_neg = np.copy(ll_rfile[cat]["EWK_TT_lin_pos"].values()), np.copy(ll_rfile[cat]["EWK_TT_lin_neg"].values())
                b2_pos, b2_neg = np.copy(ll_rfile[cat]["EWK_TT_quad_pos"].values()), np.copy(ll_rfile[cat]["EWK_TT_quad_neg"].values())
    
                #set_trace()
            
                Plotter.plot_1D((b0_pos-b0_neg)/TT, ll_edges, xlimits=(ll_edges[0], ll_edges[-1]),
                    color=color, ax=b0_ax, label=label, ylabel="$b_{0}$/$TT_{NNLO}$")
                Plotter.plot_1D((b1_pos-b1_neg)/TT, ll_edges, xlimits=(ll_edges[0], ll_edges[-1]),
                    color=color, ax=b1_ax, label=label, ylabel="$b_{1}$/$TT_{NNLO}$")
                Plotter.plot_1D((b2_pos-b2_neg)/TT, ll_edges, xlimits=(ll_edges[0], ll_edges[-1]),
                    color=color, ax=b2_ax, label=label, ylabel="$b_{2}$/$TT_{NNLO}$")
    
            #set_trace()
            b0_ax.legend(loc="upper right", title="$b_{0}$ Coefficient", ncol=2)
            b0_ax.axhline(0, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
            b0_ax.autoscale()
            b0_ax.set_ylim(b0_ax.get_ylim()[0], b0_ax.get_ylim()[1]*1.15)
            b0_ax.set_xlim((ll_edges[0], ll_edges[-1]))
            hep.cms.label(ax=b0_ax, data=False, year=year)
            b0_figname = os.path.join(pltdir, "_".join(["ll", "b0", year]))
            b0_fig.savefig(b0_figname)
            print(f"{b0_figname} written")
            plt.close(b0_fig)
    
            b1_ax.legend(loc="upper right", title="$b_{1}$ Coefficient", ncol=2)
            b1_ax.axhline(0, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
            b1_ax.autoscale()
            b1_ax.set_ylim(b1_ax.get_ylim()[0], b1_ax.get_ylim()[1]*1.15)
            b1_ax.set_xlim((ll_edges[0], ll_edges[-1]))
            hep.cms.label(ax=b1_ax, data=False, year=year)
            b1_figname = os.path.join(pltdir, "_".join(["ll", "b1", year]))
            b1_fig.savefig(b1_figname)
            print(f"{b1_figname} written")
            plt.close(b1_fig)
    
            b2_ax.legend(loc="upper right", title="$b_{2}$ Coefficient", ncol=2)
            b2_ax.axhline(0, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
            b2_ax.autoscale()
            b2_ax.set_ylim(b2_ax.get_ylim()[0], b2_ax.get_ylim()[1]*1.15)
            b2_ax.set_xlim((ll_edges[0], ll_edges[-1]))
            hep.cms.label(ax=b2_ax, data=False, year=year)
            b2_figname = os.path.join(pltdir, "_".join(["ll", "b2", year]))
            b2_fig.savefig(b2_figname)
            print(f"{b2_figname} written")
            plt.close(b2_fig)


    ## make plots comparing years per channel
if args.comp == "year":
            ## lj plots    
    channel_labels = lj_chan_labels if args.channel == "lj" else ll_chan_labels
    for channel, (chan_label, _) in channel_labels.items():
    
    #if args.channel == "lj":
    #    for channel, (chan_label, _) in lj_chan_labels.items():
            pltdir = os.path.join(outdir, channel)
            if not os.path.isdir(pltdir):
                os.makedirs(pltdir)
                
            b0_fig, b0_ax = plt.subplots(figsize=(15.0, 10.0))
            b0_fig.subplots_adjust(hspace=.07)
            b1_fig, b1_ax = plt.subplots(figsize=(15.0, 10.0))
            b1_fig.subplots_adjust(hspace=.07)
            b2_fig, b2_ax = plt.subplots(figsize=(15.0, 10.0))
            b2_fig.subplots_adjust(hspace=.07)
        
            for year, color in years_labels.items():
                cat = channel + "_" + year
    
                if args.channel == "lj":
                    TT = np.copy(lj_rfile[cat]["TT"].values())
                    b0_pos, b0_neg = np.copy(lj_rfile[cat]["EWK_TT_const_pos"].values()), np.copy(lj_rfile[cat]["EWK_TT_const_neg"].values())
                    b1_pos, b1_neg = np.copy(lj_rfile[cat]["EWK_TT_lin_pos"].values()), np.copy(lj_rfile[cat]["EWK_TT_lin_neg"].values())
                    b2_pos, b2_neg = np.copy(lj_rfile[cat]["EWK_TT_quad_pos"].values()), np.copy(lj_rfile[cat]["EWK_TT_quad_neg"].values())
                if args.channel == "ll":
                    TT = np.copy(ll_rfile[cat]["TT"].values())
                    b0_pos, b0_neg = np.copy(ll_rfile[cat]["EWK_TT_const_pos"].values()), np.copy(ll_rfile[cat]["EWK_TT_const_neg"].values())
                    b1_pos, b1_neg = np.copy(ll_rfile[cat]["EWK_TT_lin_pos"].values()), np.copy(ll_rfile[cat]["EWK_TT_lin_neg"].values())
                    b2_pos, b2_neg = np.copy(ll_rfile[cat]["EWK_TT_quad_pos"].values()), np.copy(ll_rfile[cat]["EWK_TT_quad_neg"].values())
        
                #set_trace()
            
                Plotter.plot_1D((b0_pos-b0_neg)/TT, lj_edges if args.channel == "lj" else ll_edges,
                    xlimits=(lj_edges[0], lj_edges[-1]) if args.channel == "lj" else (ll_edges[0], ll_edges[-1]), color=color, ax=b0_ax, label=year, ylabel="$b_{0}$/$TT_{NNLO}$")
                Plotter.plot_1D((b1_pos-b1_neg)/TT, lj_edges if args.channel == "lj" else ll_edges,
                    xlimits=(lj_edges[0], lj_edges[-1]) if args.channel == "lj" else (ll_edges[0], ll_edges[-1]), color=color, ax=b1_ax, label=year, ylabel="$b_{1}$/$TT_{NNLO}$")
                Plotter.plot_1D((b2_pos-b2_neg)/TT, lj_edges if args.channel == "lj" else ll_edges,
                    xlimits=(lj_edges[0], lj_edges[-1]) if args.channel == "lj" else (ll_edges[0], ll_edges[-1]), color=color, ax=b2_ax, label=year, ylabel="$b_{2}$/$TT_{NNLO}$")
            #set_trace()
            b0_ax.legend(loc="upper right", title="$b_{0}$ Coefficient", ncol=2)
            b0_ax.axhline(0, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
            b0_ax.autoscale()
            b0_ax.set_ylim(b0_ax.get_ylim()[0], b0_ax.get_ylim()[1]*1.15)
            b0_ax.set_xlim((lj_edges[0], lj_edges[-1]) if args.channel == "lj" else (ll_edges[0], ll_edges[-1]))
            b0_ax.text(0.02, 0.92, chan_label, fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=b0_ax.transAxes)
            hep.cms.label(ax=b0_ax, data=False, year=year)
            b0_figname = os.path.join(pltdir, "_".join([channel, "b0", "YearComp"]))
            b0_fig.savefig(b0_figname)
            print(f"{b0_figname} written")
            plt.close(b0_fig)
        
            b1_ax.legend(loc="upper right", title="$b_{1}$ Coefficient", ncol=2)
            b1_ax.axhline(0, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
            b1_ax.autoscale()
            b1_ax.set_ylim(b1_ax.get_ylim()[0], b1_ax.get_ylim()[1]*1.15)
            b1_ax.set_xlim((lj_edges[0], lj_edges[-1]) if args.channel == "lj" else (ll_edges[0], ll_edges[-1]))
            b1_ax.text(0.02, 0.92, chan_label, fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=b1_ax.transAxes)
            hep.cms.label(ax=b1_ax, data=False, year=year)
            b1_figname = os.path.join(pltdir, "_".join([channel, "b1", "YearComp"]))
            b1_fig.savefig(b1_figname)
            print(f"{b1_figname} written")
            plt.close(b1_fig)
        
            b2_ax.legend(loc="upper right", title="$b_{2}$ Coefficient", ncol=2)
            b2_ax.axhline(0, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
            b2_ax.autoscale()
            b2_ax.set_ylim(b2_ax.get_ylim()[0], b2_ax.get_ylim()[1]*1.15)
            b2_ax.set_xlim((lj_edges[0], lj_edges[-1]) if args.channel == "lj" else (ll_edges[0], ll_edges[-1]))
            b2_ax.text(0.02, 0.92, chan_label, fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=b2_ax.transAxes)
            hep.cms.label(ax=b2_ax, data=False, year=year)
            b2_figname = os.path.join(pltdir, "_".join([channel, "b2", "YearComp"]))
            b2_fig.savefig(b2_figname)
            print(f"{b2_figname} written")
            plt.close(b2_fig)
    


toc = time.time()
print("Total time: %.1f" % (toc - tic))
