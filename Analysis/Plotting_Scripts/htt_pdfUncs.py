import time
tic = time.time()

# matplotlib
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend("agg")
from matplotlib import rcParams
rcParams["font.size"] = 20
rcParams["savefig.format"] = "pdf"
rcParams["savefig.bbox"] = "tight"

from coffea.hist import plot
from coffea.util import load
from pdb import set_trace
import os
import Utilities.plot_tools as plt_tools
import Utilities.prettyjson as prettyjson
import re
from coffea import hist
import numpy as np
import fnmatch
import Utilities.Plotter as Plotter
from sklearn.decomposition import PCA
import itertools
import Utilities.final_analysis_binning as final_binning
import Utilities.btag_sideband_regions as btag_sidebands
import Utilities.common_features as cfeatures

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("lepton", choices=["Electron", "Muon"], help="Choose which lepton to make plots for")
parser.add_argument("--FINAL", action="store_true", help="Make plots using final PDF variations")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
base_jobid = os.environ["base_jobid"]
jobid = os.environ["jobid"]
analyzer = "htt_pdfUncs"
eos_dir = os.environ["eos_dir"]
plots_dir = os.environ["plots_dir"]

input_dir = os.path.join(eos_dir, "results", f"{args.year}_{jobid}", analyzer)
f_ext = "TOT.coffea"
fnames = sorted(["%s/%s" % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)])
hdict = plt_tools.add_coffea_files(fnames) if len(fnames) > 1 else load(fnames[0])

outdir = os.path.join(plots_dir, jobid, analyzer, args.year)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

#set_trace()
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


variables = {
    "mtt_vs_tlep_ctstar_abs" : ("$m_{t\\bar{t}}$", "|cos($\\theta^{*}_{t_{l}}$)|", linearize_binning[0], linearize_binning[1], (linearize_binning[0][0], linearize_binning[0][-1]),
        (linearize_binning[1][0], linearize_binning[1][-1]), True, None, (ctstar_binlabels, ctstar_bin_locs)),
}


    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())[args.year][f"{args.lepton}s"]
year_to_use = cfeatures.year_labels[args.year]
lumi_correction = load(os.path.join(proj_dir, "Corrections", base_jobid, "MC_LumiWeights.coffea"))[args.year][f"{args.lepton}s"]
        # scale ttJets events, split by reconstruction type, by normal ttJets lumi correction
ttJets_permcats = ["*right", "*matchable", "*unmatchable", "*sl_tau", "*other"]
names = [dataset for dataset in sorted(set([key[0] for key in hdict["mtt_vs_tlep_ctstar_abs"].values().keys()]))] # get dataset names in hists
ttJets_cats = [name for name in names if any([fnmatch.fnmatch(name, cat) for cat in ttJets_permcats])] # gets ttJets(_PS)_other, ...
if len(ttJets_cats) > 0:
    for tt_cat in ttJets_cats:
        ttJets_lumi_topo = "_".join(tt_cat.split("_")[:-2]) if "sl_tau" in tt_cat else "_".join(tt_cat.split("_")[:-1]) # gets ttJets[SL, Had, DiLep] or ttJets_PS
        ttJets_eff_lumi = lumi_correction[ttJets_lumi_topo]
        lumi_correction.update({tt_cat: ttJets_eff_lumi})

## make groups based on process
process = hist.Cat("process", "Process", sorting="placement")
process_cat = "dataset"
process_groups = {"TT" : names}

    ## make plots
    # scale and group hists by process
for hname in hdict.keys():
    if "cutflow" in hname: continue
    #set_trace()
    hdict[hname].scale(lumi_correction, axis="dataset") # scale hists
    hdict[hname] = hdict[hname].group(process_cat, process, process_groups) # group by process
    hdict[hname] = hdict[hname][:, :, :, args.lepton].integrate("leptype").integrate("process") # only pick out specified lepton

for hname in variables.keys():
    if hname not in hdict.keys():
        raise ValueError(f"{hname} not found in file")
    #set_trace()
    histo = hdict[hname].copy() # process, sys, jmult

    xtitle, ytitle, xrebinning, yrebinning, x_lims, y_lims, withData, plot_xlabels, plot_ylabels = variables[hname]

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

    ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
    for jmult in sorted(set([key[1] for key in histo.values().keys()])):
        pltdir = os.path.join(outdir, args.lepton, jmult)
        if not os.path.isdir(pltdir):
            os.makedirs(pltdir)

        nosys_hline = Plotter.linearize_hist(histo["nosys", jmult].integrate("jmult").integrate("sys").copy())
        x_lims = (0, nosys_hline.dense_axes()[0].centers().size)

        ## init dict of for highest variations
        pdf_vars_dict = {}

        ### plot alphaS variations
        fig, ax = plt.subplots(figsize=(15.0, 10.0))
        fig.subplots_adjust(hspace=.07)

        alphaSup = Plotter.linearize_hist(histo["alphaSUp", jmult].integrate("jmult").integrate("sys"))
        alphaSdw = Plotter.linearize_hist(histo["alphaSDown", jmult].integrate("jmult").integrate("sys"))
        alphaSUp_masked_vals, alphaSUp_masked_bins = Plotter.get_ratio_arrays(num_vals=alphaSup.values()[()], denom_vals=nosys_hline.values()[()], input_bins=nosys_hline.dense_axes()[0].edges())
        alphaSDown_masked_vals, alphaSDown_masked_bins = Plotter.get_ratio_arrays(num_vals=alphaSdw.values()[()], denom_vals=nosys_hline.values()[()], input_bins=nosys_hline.dense_axes()[0].edges())
        ax.fill_between(alphaSUp_masked_bins, alphaSUp_masked_vals, y2=1., facecolor="r", step="post", alpha=0.5, label="Up")
        ax.fill_between(alphaSDown_masked_bins, alphaSDown_masked_vals, y2=1., facecolor="b", step="post", alpha=0.5, label="Down")
        ax.legend(loc="upper right", title="$\\alpha_{S}$")
        ax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
        ax.autoscale()
        ax.set_xlim(x_lims)
        ax.set_xlabel(xtitle)
        ax.set_ylabel("Ratio to Nominal")

            # draw vertical lines separating ctstar bins
        bin_sep_lines = [histo["nosys", jmult].integrate("jmult").integrate("sys").values()[()].shape[0]*ybin for ybin in range(1, histo["nosys", jmult].integrate("jmult").integrate("sys").values()[()].shape[1])]
        [ax.axvline(binline, color="k", linestyle="--") for binline in bin_sep_lines]

        [ax.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
            xytext=(0, 400), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0) for idx, label in enumerate(ctstar_binlabels)]
            #xytext=(0, 10), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0) for idx, label in enumerate(ctstar_binlabels)]

        #set_trace()
        ax.set_xticks(mtt_bin_inds_to_plot)
        ax.set_xticklabels(mtt_bins_to_plot)

            # add lepton/jet multiplicity label
        ax.text(
            0.02, 0.92, cfeatures.channel_labels[f'{args.lepton}_{jmult}'],
            fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
        )
        hep.cms.label(ax=ax, label="Preliminary", data=False, year=year_to_use, lumi=round(data_lumi_year/1000., 1))

        #set_trace()
        figname = os.path.join(pltdir, "_".join([jmult, args.lepton, hname, "alphaS"]))
        fig.savefig(figname)
        print(f"{figname} written")
        plt.close(fig)
        #set_trace()


        for pdf_sys in sorted(set([key[0] for key in histo.values().keys()])):
            if (pdf_sys == "nosys") or (pdf_sys == "pdf_0"): continue
            if "alphaS" in pdf_sys: continue # plot alphaS up/down together
            print(f"{args.year}, {args.lepton}, {jmult}, {pdf_sys}")
            #print(*[jmult, hname, pdf_sys], sep=", ")
            hslice = histo[pdf_sys, jmult].integrate("jmult").integrate("sys")

            # plot linearized view 
            fig, ax = plt.subplots(figsize=(15.0, 10.0))
            fig.subplots_adjust(hspace=.07)

            hline = Plotter.linearize_hist(hslice)
            hline_masked_vals, hline_masked_bins = Plotter.get_ratio_arrays(num_vals=hline.values()[()], denom_vals=nosys_hline.values()[()], input_bins=nosys_hline.dense_axes()[0].edges())
                # save sum of distance from nominal and ratio-to-nominal array
            pdf_vars_dict[np.sum(np.sqrt(np.square(hline.values()[()]-nosys_hline.values()[()])))] = [pdf_sys, hline_masked_vals]
            #set_trace()

            ax.fill_between(hline_masked_bins, hline_masked_vals, y2=1., facecolor="r", step="post", alpha=0.5, label=pdf_sys)
            ax.legend(loc="upper right")
            ax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
            ax.autoscale()
            ax.set_xlim(x_lims)
            ax.set_xlabel(xtitle)
            ax.set_ylabel("Ratio to Nominal")

                # draw vertical lines separating ctstar bins
            bin_sep_lines = [hslice.values()[()].shape[0]*ybin for ybin in range(1, hslice.values()[()].shape[1])]
            [ax.axvline(binline, color="k", linestyle="--") for binline in bin_sep_lines]

            [ax.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                xytext=(0, 400), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0) for idx, label in enumerate(ctstar_binlabels)]

            #set_trace()
            ax.set_xticks(mtt_bin_inds_to_plot)
            ax.set_xticklabels(mtt_bins_to_plot)

                # add lepton/jet multiplicity label
            ax.text(
                0.02, 0.92, cfeatures.channel_labels[f'{args.lepton}_{jmult}'],
                fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
            )
            hep.cms.label(ax=ax, label="Preliminary", data=False, year=year_to_use, lumi=round(data_lumi_year/1000., 1))

            #set_trace()
            figname = os.path.join(pltdir, "_".join([jmult, args.lepton, hname, pdf_sys]))
            fig.savefig(figname)
            print(f"{figname} written")
            plt.close(fig)
            #set_trace()

        #set_trace()
        ### make plots comparing top 5 highest impact pdf uncs + alphaS variations
        top_vars = np.sort(np.array(list(pdf_vars_dict.keys())))[::-1][:5]
        fig, ax = plt.subplots(figsize=(15.0, 10.0))
        fig.subplots_adjust(hspace=.07)

        ax.step(alphaSUp_masked_bins, alphaSUp_masked_vals, where="post", **{"linestyle": "-", "label": "$\\alpha_{S}$ Up"})
        ax.step(alphaSDown_masked_bins, alphaSDown_masked_vals, where="post", **{"linestyle": "-", "label": "$\\alpha_{S}$ Down"})

        for pdf_var in top_vars:
            pdf_sys, pdf_vals = pdf_vars_dict[pdf_var]
            ax.step(alphaSDown_masked_bins, pdf_vals, where="post", **{"linestyle": "-", "label": pdf_sys})

        ax.legend(loc="upper right", ncol=3)
        ax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
        ax.autoscale()
        #set_trace()
        ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.01)
        ax.set_xlim(x_lims)
        ax.set_xlabel(xtitle)
        ax.set_ylabel("Ratio to Nominal")

            # draw vertical lines separating ctstar bins
        [ax.axvline(binline, color="k", linestyle="--") for binline in bin_sep_lines]

        [ax.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
            xytext=(0, 50), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0) for idx, label in enumerate(ctstar_binlabels)]
            #xytext=(0, 400), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0) for idx, label in enumerate(ctstar_binlabels)]

        #set_trace()
        ax.set_xticks(mtt_bin_inds_to_plot)
        ax.set_xticklabels(mtt_bins_to_plot)

            # add lepton/jet multiplicity label
        ax.text(
            0.02, 0.92, cfeatures.channel_labels[f'{args.lepton}_{jmult}'],
            fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
        )
        hep.cms.label(ax=ax, label="Preliminary", data=False, year=year_to_use, lumi=round(data_lumi_year/1000., 1))

        #set_trace()
        figname = os.path.join(pltdir, "_".join([jmult, args.lepton, hname, "PDF_Comp"]))
        fig.savefig(figname)
        print(f"{figname} written")
        plt.close(fig)


toc = time.time()
print("Total time: %.1f" % (toc - tic))
