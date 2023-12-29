#! /bin/env python

import time
tic = time.time()

# matplotlib
import matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.style.CMS)
plt.switch_backend("agg")
from matplotlib import rcParams
rcParams["font.size"] = 20
rcParams["savefig.format"] = "pdf"
rcParams["savefig.bbox"] = "tight"

from pdb import set_trace
import os
import Utilities.plot_tools as plt_tools
import Utilities.prettyjson as prettyjson
import hist
import numpy as np
#import Utilities.Plotter as Plotter
import Utilities.final_analysis_binning as final_binning
import Utilities.common_features as cfeatures
import uproot

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
args = parser.parse_args()

jobid = os.environ["jobid"]
base_jobid = os.environ["base_jobid"]
proj_dir = os.environ["PROJECT_DIR"]
plot_outdir = os.environ["plots_dir"]

basedir = "root://cmseos.fnal.gov//store/user/jdulemba/HeavyHiggsFitter/Thesis_constrainedYt_noEtaT_v1"
rname = f"IMPACTSv35_lj_data_A_m400_relw5p0_sig1.root"
rfile = uproot.open(os.path.join(basedir, rname))

year_label = cfeatures.year_labels[args.year]

outdir = os.path.join(plot_outdir, jobid, "Limits_htt_btag_sb_regions", os.path.basename(basedir).split(f"{jobid}_")[-1], "prefit", "new", year_label.replace("VFP",""))
if not os.path.isdir(outdir):
    os.makedirs(outdir)
    
linearize_binning = (
    final_binning.mtt_binning,
    final_binning.ctstar_abs_binning
)

#set_trace()
ctstar_binlabels = [r"%s $\in [%s, %s)$" % (cfeatures.variable_names_to_labels["tlep_ctstar_abs"], linearize_binning[1][bin], linearize_binning[1][bin+1]) for bin in range(len(linearize_binning[1])-1)]
ctstar_binlabels[-1] = r"%s $\in [%s, %s]$" % (cfeatures.variable_names_to_labels["tlep_ctstar_abs"], linearize_binning[1][-2], linearize_binning[1][-1])
ctstar_bin_locs = np.linspace((len(linearize_binning[0])-1)/2, (len(linearize_binning[0])-1)*(len(linearize_binning[1])-1) - (len(linearize_binning[0])-1)/2, len(linearize_binning[1])-1)
mtt_vals_to_plot = np.array([400, 600, 1000])
mtt_tiled_labels = np.tile(linearize_binning[0][:-1], linearize_binning[1].size-1)
mtt_bin_inds_to_plot = np.where(np.in1d(mtt_tiled_labels, mtt_vals_to_plot))[0]
mtt_bins_to_plot = np.tile(mtt_vals_to_plot, linearize_binning[1].size-1)
xtitle, ytitle = cfeatures.variable_names_to_labels["mtt"], cfeatures.variable_names_to_labels["tlep_ctstar_abs"]

## get data lumi and scale MC by lumi
data_lumi_file = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())
data_lumi_year = data_lumi_file[args.year]

mc_opts = {
    "mcorder" : ["ttJets", "singlet", "EWQCD"],
    "maskData" : False,
    "overflow" : "none",
    "ncol" : 2,
}

vlines = [(len(linearize_binning[0])-1)*ybin for ybin in range(1, len(linearize_binning[1])-1)]
    # hardcoded groups
procs_groups_dict = {
    "data" : ["data"],
    "ttJets" : ["TT", "TT_ytconst", "TT_ytlin", "TT_ytquad"],
    "singlet" : ["TB", "TW", "TQ"],
    "EWQCD" : ["EWQCD"],
}

procs_list = sorted({x for v in procs_groups_dict.values() for x in v})

#for lepton in ["Muon", "Electron"]:
#    for jmult in ["3Jets", "4PJets"]:
for lepton in ["Muon"]:
    for jmult in ["3Jets"]:
        print(f"Making plots for {lepton} {jmult}\n")

        channel = f"mu{(jmult).lower()}" if lepton == "Muon" else f"e{(jmult).lower()}"
        channel_label = f"{cfeatures.channel_labels[f'{lepton}_{jmult}']}"
        dirname = f"{channel}_{year_label.replace('VFP','')}_prefit"
        lumi_to_use = round(data_lumi_year[f"{lepton}s"]/1000., 1)

        #set_trace()        
        #histo = None
        histo_dict = {}

        new_histo = hist.Hist(
            hist.axis.StrCategory([], growth=True, name="dataset", label="dataset"),
            hist.axis.Variable(bin_edges, name="x", label="mtt_ctstar"),
            "Weight"
            #storage="weight",
            #name="Counts"
        )
        for proc in procs_list:
            vals = np.copy(rfile[dirname][proc].values()) if proc == "data" else np.copy(rfile[dirname][f"{channel}_{year_label.replace('VFP','')}_{proc}"].values())
            variances = np.copy(rfile[dirname][proc].variances()) if proc == "data" else np.copy(rfile[dirname][f"{channel}_{year_label.replace('VFP','')}_{proc}"].variances())
        
            if not histo_dict:
                bin_edges = np.arange(len(vals)+1, dtype=int)
            new_histo.fill(dataset=proc, x=np.zeros(vals.size))
            new_histo[{"dataset" : proc}].values()[:], new_histo[{"dataset" : proc}].variances()[:] = vals, variances if proc == "data" else np.zeros_like(variances)
            #new_histo.view().value[:], new_histo.view().variance[:] = vals, variances if proc == "data" else np.zeros_like(variances)
            #new_histo.fill(dataset=proc, x=np.zeros(0))
            histo_dict[proc] =  hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([vals, variances], axis=-1)) if proc == "data" else\
                hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([vals, np.zeros_like(variances)], axis=-1))
                #hist2 = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([vals, variances], axis=-1))
                hist1 = hist.Hist(hist.axis.Variable(bin_edges), hist.axis.StrCategory([proc], name="proc"), data=np.stack([vals, variances], axis=-1))
                hist2 = hist.Hist(hist.axis.Variable(bin_edges), hist.axis.StrCategory([proc], name="proc"), data=np.stack([vals, np.zeros_like(variances)], axis=-1))
                #hist1 = hist.Hist(hist.axis.Variable(bin_edges), hist.axis.StrCategory([proc], name="proc"), data=np.stack([vals, np.zeros_like(variances)], axis=-1))
            #new_histo =  hist.Hist(
            #    hist.axis.StrCategory([], growth=True, name="dataset", label="dataset"),
            #    hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([vals, variances], axis=-1)) if proc == "data" else\
            #    hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([vals, np.zeros_like(variances)], axis=-1))
            #)
        
        set_trace()
        mc_histo = hist.Stack.from_dict({"ttJets" : histo_dict["TT"]+histo_dict["TT_ytconst"]+histo_dict["TT_ytlin"]+histo_dict["TT_ytquad"], "singlet" : histo_dict["TB"]+histo_dict["TW"]+histo_dict["TQ"], "EWQCD" : histo_dict["EWQCD"]})
        data_histo = hist.Stack.from_dict({"data" : histo_dict["data"]})
        #histo = hist.Stack.from_dict({"data" : hist2, "EWQCD" : hist1})
        nbins = len(bin_edges)-1
        
        #set_trace()
        fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True, figsize=(15.0, 10.0))
        fig.subplots_adjust(hspace=.07)

        mc_histo.plot(histtype="fill", stack=True, ax=ax)
        #histo[{"x": sum, "y": sum}].plot1d(overlay="samp", histtype="fill", stack=True)        
        #### plot non-signal MC prediction
        #plot.plot1d(histo[Plotter.mc_samples], overlay=histo[Plotter.mc_samples].axes()[0].name,
        #    ax=ax, clear=False, stack=True, line_opts=None, fill_opts=Plotter.stack_fill_opts,
        #    error_opts=None, #stack_error_opts,
        #    order=mc_opts["mcorder"], overflow="none", overlay_overflow="none",
        #)
        pred_error = np.copy(np.sqrt(rfile[dirname]["model"].variances()))
        pred_yield = histo[Plotter.mc_samples].integrate("process").values()[()]
        ax.fill_between(
            x=histo.dense_axes()[0].edges(),
            y1=np.r_[pred_yield-pred_error, (pred_yield-pred_error)[-1]],
            y2=np.r_[pred_yield, pred_yield[-1]],
            **{"step": "post", "label": "Total Unc.", "hatch": "///",
                "facecolor": "none", "edgecolor": (0, 0, 0, 0.5), "linewidth": 0,
            }
        )
        ax.fill_between(
            x=histo.dense_axes()[0].edges(),
            y1=np.r_[pred_yield+pred_error, (pred_yield+pred_error)[-1]],
            y2=np.r_[pred_yield, pred_yield[-1]],
            **{"step": "post", "hatch": "///",
                "facecolor": "none", "edgecolor": (0, 0, 0, 0.5), "linewidth": 0,
            }
        )
        
        #### plot data
        plot.plot1d(histo[Plotter.data_samples], overlay=histo[Plotter.data_samples].axes()[0].name,
            ax=ax, clear=False, error_opts=Plotter.hstyles["data_err_opts"], overflow="none", overlay_overflow="none",
        )
        
        ax.autoscale()
        ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.3)
        #ax.set_ylim(0, ax.get_ylim()[1]*1.3)
        ax.set_xlabel(None)
        ax.set_xlim((0, nbins))
        
            ## set legend and corresponding colors
        handles, labels = ax.get_legend_handles_labels()
        #set_trace()
        for idx, sample in enumerate(labels):
            if sample == "data" or sample == "Observed": continue
            if isinstance(handles[idx], matplotlib.lines.Line2D): continue
            facecolor, legname = plt_tools.get_styles(sample, Plotter.hstyles)
            handles[idx].set_facecolor(facecolor)
            labels[idx] = legname
        # call ax.legend() with the new values
        ax.legend(handles,labels, loc="upper right", title=mc_opts["legend_title"], ncol=mc_opts["ncol"]) if "legend_title" in mc_opts.keys() else ax.legend(handles,labels, loc="upper right", ncol=mc_opts["ncol"])
        
            ## plot data/MC ratio
        data_yield, data_variance = histo[Plotter.data_samples].integrate("process").values(sumw2=True)[()]
        #set_trace()
        rax.errorbar(
            x=histo.dense_axes()[0].centers(),
            y=data_yield/pred_yield,
            yerr=(data_yield+np.sqrt(data_variance))/pred_yield - data_yield/pred_yield,
            **{"markersize": 10., "elinewidth": 1, "fmt" : ".k"}
        )
            ## plot MC/MC ratio 
        rax.fill_between(
            x=histo.dense_axes()[0].edges(),
            y1=np.r_[(pred_yield+pred_error)/pred_yield, ((pred_yield+pred_error)/pred_yield)[-1]],
            y2=np.r_[(pred_yield-pred_error)/pred_yield, ((pred_yield-pred_error)/pred_yield)[-1]],
            **{"step": "post", "facecolor": (0, 0, 0, 0.3), "linewidth": 0}
        )
        rax.axhline(y=1., color="k", linestyle="-") # plot ratio at 1
        
            ## set axes labels and titles
        rax.set_ylabel("$\dfrac{data}{Pred.}$")
        rax.set_ylim(0.8, 1.2)
        #rax.set_ylim(0.5, 1.5)
        rax.set_xlim((0, nbins))
        rax.set_xlabel(xtitle)
        
                ## draw vertical lines for distinguishing different ctstar bins
        [ax.axvline(vline, color="k", linestyle="--") for vline in vlines]
        [rax.axvline(vline, color="k", linestyle="--") for vline in vlines]
        
            # plot unrolled x and y labels for each bin
        for idx, label in enumerate(ctstar_binlabels):
            ax.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                xytext=(0, 250), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)
        
        rax.set_xticks(mtt_bin_inds_to_plot)
        rax.set_xticklabels(mtt_bins_to_plot)
        
        
            # add lepton/jet multiplicity label
        ax.text(
            0.02, 0.86, f"{channel_label}\nprefit",
            fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
        )
        hep.cms.label(ax=ax, data=True, label="Preliminary", year=year_label, lumi=lumi_to_use)
        
        #set_trace()
        figname = os.path.join(outdir, dirname)
        fig.savefig(figname)
        print(f"{figname} written")
        plt.close(fig)


toc = time.time()
print("Total time: %.1f" % (toc - tic))
