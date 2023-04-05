#! /bin/env python

import time
tic = time.time()

# matplotlib
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import mplhep as hep
plt.style.use(hep.style.CMS)
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
import Utilities.systematics as systematics
from copy import deepcopy
from Utilities.styles import styles as styles
import Utilities.final_analysis_binning as final_binning
import Utilities.common_features as cfeatures

base_jobid = os.environ["base_jobid"]

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2018"], help="What year is the ntuple from.")
#parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("lepton", choices=["Electron", "Muon"], help="Choose which lepton to make plots for")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
analyzer = "HBSR_TopPt_Reweighting"
eos_dir = os.environ["eos_dir"]
plots_dir = os.environ["plots_dir"]

#set_trace()
base_ext = "*TOT.coffea"
input_dir = os.path.join(eos_dir, "results", f"{args.year}_{jobid}", analyzer)
fnames = fnmatch.filter(os.listdir(input_dir), base_ext)
fnames = [os.path.join(input_dir, fname) for fname in fnames]
hdict = load(fnames[0])
#set_trace()

outdir = os.path.join(plots_dir, jobid, analyzer, args.year)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

year_to_use = cfeatures.year_labels[args.year]

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

# binning and bin labels for phi x eta
phi_eta_binning = (
    np.array([-3.2, -2.4, -1.57, -0.87, 0., 0.87, 1.57, 2.4, 3.2]), # phi
    np.array([-2.5, -1.3, -0.7, 0., 0.7, 1.3, 2.5]) # eta binning
)
eta_binlabels = ["%s $\leq$ $\\eta$ $\leq$ %s" % (phi_eta_binning[1][bin], phi_eta_binning[1][bin+1]) for bin in range(len(phi_eta_binning[1])-1)]
eta_bin_locs = np.linspace((len(phi_eta_binning[0])-1)/2, (len(phi_eta_binning[0])-1)*(len(phi_eta_binning[1])-1) - (len(phi_eta_binning[0])-1)/2, len(phi_eta_binning[1])-1)
phi_binlabels = ["%s $\leq$ $\\phi$ $\leq$ %s" % (phi_eta_binning[0][bin], phi_eta_binning[0][bin+1]) for bin in range(len(phi_eta_binning[0])-1)]*len(eta_binlabels)


variables = {
    "mtt_vs_tlep_ctstar_abs" : (cfeatures.variable_names_to_labels["mtt"], cfeatures.variable_names_to_labels["tlep_ctstar_abs"], linearize_binning[0], linearize_binning[1],
        (linearize_binning[0][0], linearize_binning[0][-1]), (linearize_binning[1][0], linearize_binning[1][-1]), False, None, (ctstar_binlabels, ctstar_bin_locs)),
        #(linearize_binning[0][0], linearize_binning[0][-1]), (linearize_binning[1][0], linearize_binning[1][-1]), True, None, (ctstar_binlabels, ctstar_bin_locs)),
    "Jets_phi_vs_eta" : (cfeatures.variable_names_to_labels["Jets_phi"], cfeatures.variable_names_to_labels["Jets_eta"],
        phi_eta_binning[0], phi_eta_binning[1], (-3.3, 3.3), (-2.6, 2.6), True, phi_binlabels, (eta_binlabels, eta_bin_locs)),
    "Lep_phi_vs_eta" : (cfeatures.variable_names_to_labels["Lep_phi"] % cfeatures.objtypes[args.lepton], cfeatures.variable_names_to_labels["Lep_eta"] % cfeatures.objtypes[args.lepton],
        phi_eta_binning[0], phi_eta_binning[1], (-3.3, 3.3), (-2.6, 2.6), True, phi_binlabels, (eta_binlabels, eta_bin_locs)),
    "Jets_njets" : (cfeatures.variable_names_to_labels["Jets_njets"], 1, (0, 15), True),
    "mtt" : (cfeatures.variable_names_to_labels["mtt"], 4, (200., 2000.), False),
    "tlep_ctstar_abs" : (cfeatures.variable_names_to_labels["tlep_ctstar_abs"], 1, (0., 1.), False),
    "tlep_ctstar_abs" : (cfeatures.variable_names_to_labels["tlep_ctstar_abs"], 1, (0., 1.), True),
    "mthad" : (cfeatures.variable_names_to_labels["mthad"], 2, (0., 300.), True),
    "mWHad" : (cfeatures.variable_names_to_labels["mWHad"], 2, (0., 300.), True),
    "mWLep" : (cfeatures.variable_names_to_labels["mWLep"], 2, (0., 300.), True),
    "pt_thad" : (cfeatures.variable_names_to_labels["pt_thad"], 2, (0., 500.), True),
    "pt_tlep" : (cfeatures.variable_names_to_labels["pt_tlep"], 2, (0., 500.), True),
    "pt_tt" : (cfeatures.variable_names_to_labels["pt_tt"], 2, (0., 500.), True),
    "eta_thad" : (cfeatures.variable_names_to_labels["eta_thad"], 2, (-4., 4.), True),
    "eta_tlep" : (cfeatures.variable_names_to_labels["eta_tlep"], 2, (-4., 4.), True),
    "eta_tt" : (cfeatures.variable_names_to_labels["eta_tt"], 2, (-4., 4.), True),
    "tlep_ctstar" : (cfeatures.variable_names_to_labels["tlep_ctstar"], 2, (-1., 1.), False),
    #"tlep_ctstar" : (cfeatures.variable_names_to_labels["tlep_ctstar"], 2, (-1., 1.), True),
    "full_disc" : (cfeatures.variable_names_to_labels["full_disc"], 2, (5, 25.), True),
    "mass_disc" : (cfeatures.variable_names_to_labels["mass_disc"], 2, (0, 20.), True),
    "ns_disc" : (cfeatures.variable_names_to_labels["ns_disc"], 2, (3., 10.), True),
    "ns_dist" : (cfeatures.variable_names_to_labels["ns_dist"], 1, (0., 150.), True),
    "Jets_pt" : (cfeatures.variable_names_to_labels["Jets_pt"], 1, (0., 300.), True),
    "Jets_eta" : (cfeatures.variable_names_to_labels["Jets_eta"], 2, (-2.6, 2.6), True),
    "Jets_phi" : (cfeatures.variable_names_to_labels["Jets_phi"], 2, (-4., 4.), True),
    "Jets_LeadJet_pt" : (cfeatures.variable_names_to_labels["Jets_LeadJet_pt"], 1, (0., 300.), True),
    "Jets_LeadJet_eta" : (cfeatures.variable_names_to_labels["Jets_LeadJet_eta"], 2, (-2.6, 2.6), True),
    "Jets_DeepCSV_bDisc" : (cfeatures.variable_names_to_labels["Jets_DeepCSV_bDisc"], 1, (-0.01, 1.), True),
    "Jets_DeepJet_bDisc" : (cfeatures.variable_names_to_labels["Jets_DeepJet_bDisc"], 1, (-0.01, 1.), True),
    "Lep_pt" : (cfeatures.variable_names_to_labels["Lep_pt"] % cfeatures.objtypes[args.lepton], 1, (0., 300.), True),
    "Lep_eta" : (cfeatures.variable_names_to_labels["Lep_eta"] % cfeatures.objtypes[args.lepton], 2, (-2.6, 2.6), True),
    "Lep_phi" : (cfeatures.variable_names_to_labels["Lep_phi"] % cfeatures.objtypes[args.lepton], 2, (-4, 4), True),
    "Lep_iso" : (cfeatures.variable_names_to_labels["Lep_iso"] % cfeatures.objtypes[args.lepton], 1, (0., 1.), True),
    "MT" : (cfeatures.variable_names_to_labels["MT"], 1, (0., 300.), True),
    "MET_pt" : (cfeatures.variable_names_to_labels["MET_pt"], 1, (0., 300.), True),
    "MET_phi" : (cfeatures.variable_names_to_labels["MET_phi"], 1, (-3.2, 3.2), True),
}


    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())[args.year]
lumi_correction = load(os.path.join(proj_dir, "Corrections", base_jobid, "MC_LumiWeights.coffea"))[args.year][f"{args.lepton}s"]
        # scale ttJets events, split by reconstruction type, by normal ttJets lumi correction
ttJets_permcats = ["*right", "*matchable", "*unmatchable", "*sl_tau", "*other"]
#set_trace()
names = [dataset for dataset in sorted(set([key[0] for key in hdict["mtt_vs_tlep_ctstar_abs"].values().keys()]))] # get dataset names in hists
ttJets_cats = [name for name in names if any([fnmatch.fnmatch(name, cat) for cat in ttJets_permcats])] # gets ttJets(_PS)_other, ...
if len(ttJets_cats) > 0:
    for tt_cat in ttJets_cats:
        ttJets_lumi_topo = "_".join(tt_cat.split("_")[:-2]) if "sl_tau" in tt_cat else "_".join(tt_cat.split("_")[:-1]) # gets ttJets[SL, Had, DiLep] or ttJets_PS
        ttJets_eff_lumi = lumi_correction[ttJets_lumi_topo]
        lumi_correction.update({tt_cat: ttJets_eff_lumi})

## make groups based on process
process_groups = plt_tools.make_dataset_groups(args.lepton, args.year, samples=names, bkgdict="templates", sigdict="MC_indiv")
tt_LHEscale_wts_name_dict = {
    "FACTORDown" : "uF_down",
    "FACTORUp"   : "uF_up",
    "RENORMDown" : "uR_down",
    "RENORMUp"   : "uR_up",
}
if sorted(set(process_groups.keys())) != ["TT"]:
    raise ValueError("Expected only ttbar events in this file")

    # scale and group hists by process
for hname in variables.keys():
    if hname not in hdict.keys():
        print(f"{hname} not found in file")
        continue
    if "cutflow" in hname: continue
    #set_trace()
    hdict[hname].scale(lumi_correction, axis="dataset") # scale hists
    for nom_tt in ["ttJetsSL", "ttJetsDiLep", "ttJetsHad"]:
        for cat in ttJets_permcats:
            # rescale LHEscale systematics correctly
            for sysname, dname in tt_LHEscale_wts_name_dict.items():
                if f"{nom_tt}_{cat[1:]}" not in ttJets_cats: continue
                #print(f"{nom_tt}_{cat[1:]}_{dname}")
                lhe_scale = lumi_correction[f"{nom_tt}_{dname}"]/lumi_correction[f"{nom_tt}_{cat[1:]}"]
                hdict[hname].scale({(f"{nom_tt}_{cat[1:]}", sysname) : lhe_scale}, axis=("dataset", "sys"))
    hdict[hname] = hdict[hname][:, :, :, args.lepton].integrate("btag").integrate("dataset").integrate("leptype")

styles_dict = {
    "nosys" : {"color":"k", "linestyle":"-", "label":"Nominal"},
    "TopPtUp" : {"color":"r", "linestyle":"-", "label":"$p_{T}$($t$) Up"},
    "TopPtDown" : {"color":"r", "linestyle":"--", "label":"$p_{T}$($t$) Down"},
    "RENORMUp" : {"color":"b", "linestyle":"-", "label":"$\\mu_{R}$ Up"},
    "RENORMDown" : {"color":"b", "linestyle":"--", "label":"$\\mu_{R}$ Down"},
    "FACTORUp" : {"color":"g", "linestyle":"-", "label":"$\\mu_{F}$ Up"},
    "FACTORDown" : {"color":"g", "linestyle":"--", "label":"$\\mu_{F}$ Down"},
}

    ## make plots
for hname in variables.keys():
    if hname not in hdict.keys():
        print(f"{hname} not found in file")
        continue
        #raise ValueError(f"{hname} not found in file")

    #set_trace()
    histo = hdict[hname]

    is2d = histo.dense_dim() == 2
    if histo.dense_dim() == 1:
        xtitle, rebinning, x_lims, withData = variables[hname]
        orig_xtitle = xtitle
        if rebinning != 1:
            xaxis_name = histo.dense_axes()[0].name
            histo = histo.rebin(xaxis_name, rebinning)

    if is2d:
        xtitle, ytitle, xrebinning, yrebinning, x_lims, y_lims, withData, plot_xlabels, plot_ylabels = variables[hname]

        xaxis_name = histo.dense_axes()[0].name
        yaxis_name = histo.dense_axes()[1].name
            ## rebin x axis
        if isinstance(xrebinning, np.ndarray):
            new_xbins = hist.Bin(xaxis_name, xaxis_name, xrebinning)
        elif isinstance(xrebinning, float) or isinstance(xrebinning, int):
            new_xbins = xrebinning
            ## rebin y axis
        if isinstance(yrebinning, np.ndarray):
            new_ybins = hist.Bin(yaxis_name, yaxis_name, yrebinning)
        elif isinstance(yrebinning, float) or isinstance(yrebinning, int):
            new_ybins = yrebinning
        histo = histo.rebin(yaxis_name, new_ybins).rebin(xaxis_name, new_xbins)

    ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
    for jmult in sorted(set([key[1] for key in histo.values().keys()])):
        print(*[jmult, hname], sep=", ")
        pltdir = os.path.join(outdir, args.lepton, jmult)
        if not os.path.isdir(pltdir):
            os.makedirs(pltdir)

        fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True, figsize=(15.0, 10.0)) if is2d else plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
        fig.subplots_adjust(hspace=.07)

        ratio_fig, ratio_ax = plt.subplots(figsize=(15.0, 10.0)) if is2d else plt.subplots()
        ratio_fig.subplots_adjust(hspace=.07)

        hslice = histo[:, jmult].integrate("jmult")

        if "disc" in hname:
            xtitle = orig_xtitle.replace("j", "3") if jmult == "3Jets" else orig_xtitle.replace("j", "4")

        if hname == "Lep_iso":
            x_lims = (0., 0.15) if (args.lepton == "Muon") else (0., 0.1)

            # plot nominal
        hep.plot.histplot(Plotter.linearize_hist(hslice["nosys"].integrate("sys")).values()[()], Plotter.linearize_hist(hslice["nosys"].integrate("sys")).dense_axes()[0].edges(), ax=ax, histtype="step", **styles_dict["nosys"]) if is2d \
            else hep.plot.histplot(hslice.values()[("nosys",)], hslice.dense_axes()[0].edges(), ax=ax, histtype="step", **styles_dict["nosys"])

            # plot systematics
        for sys in sorted(set([key[0] for key in hslice.values().keys()])):
            if sys == "nosys": continue
            hep.plot.histplot(Plotter.linearize_hist(hslice[sys].integrate("sys")).values()[()], Plotter.linearize_hist(hslice[sys].integrate("sys")).dense_axes()[0].edges(), ax=ax, histtype="step", **styles_dict[sys]) if is2d \
                else hep.plot.histplot(hslice.values()[(sys,)], hslice.dense_axes()[0].edges(), ax=ax, histtype="step", **styles_dict[sys])

            masked_vals, masked_bins = Plotter.get_ratio_arrays(num_vals=Plotter.linearize_hist(hslice[sys].integrate("sys")).values()[()], denom_vals=Plotter.linearize_hist(hslice["nosys"].integrate("sys")).values()[()],
                    input_bins=Plotter.linearize_hist(hslice["nosys"].integrate("sys")).dense_axes()[0].edges()) if is2d \
                else Plotter.get_ratio_arrays(num_vals=hslice.values()[(sys,)], denom_vals=hslice.values()[("nosys",)], input_bins=hslice.dense_axes()[0].edges())
            rax.step(masked_bins, masked_vals, where="post", **styles_dict[sys])
            ratio_ax.step(masked_bins, masked_vals, where="post", **styles_dict[sys])

        ax.set_yscale("log")
        ax.legend(loc="upper right", ncol=2)
        ax.autoscale()
        ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.5)
        #ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.15)
        ax.set_xlim((0, Plotter.linearize_hist(hslice["nosys"].integrate("sys")).dense_axes()[0].centers().size) if is2d else x_lims)
        ax.set_xlabel(None)
        ax.set_ylabel("Events")

        rax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
        rax.autoscale()
        rax.set_xlim((0, Plotter.linearize_hist(hslice["nosys"].integrate("sys")).dense_axes()[0].centers().size) if is2d else x_lims)
        #rax.set_xlim(x_lims)
        rax.set_xlabel(xtitle)
        rax.set_ylabel("Ratio to Nominal")

        ratio_ax.legend(loc="upper right", ncol=2)
        ratio_ax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
        ratio_ax.autoscale()
        ratio_ax.set_ylim(ratio_ax.get_ylim()[0], ratio_ax.get_ylim()[1]*1.15)
        ratio_ax.set_xlim((0, Plotter.linearize_hist(hslice["nosys"].integrate("sys")).dense_axes()[0].centers().size) if is2d else x_lims)
        #ratio_ax.set_xlim(x_lims)
        ratio_ax.set_xlabel(xtitle)
        ratio_ax.set_ylabel("Ratio to Nominal")

        if is2d:
            #set_trace()
                ## draw vertical lines for distinguishing different ctstar bins
            vlines = [Plotter.linearize_hist(hslice["nosys"].integrate("sys")).dense_axes()[0].centers().size * ybin/5 for ybin in range(1, 5)]
            #vlines = [x_lims[1]*ybin/5 for ybin in range(1, 5)]
            [ax.axvline(vline, color="k", linestyle="--") for vline in vlines]
            [rax.axvline(vline, color="k", linestyle="--") for vline in vlines]
            [ratio_ax.axvline(vline, color="k", linestyle="--") for vline in vlines]

                # plot unrolled x and y labels for each bin
            ## plot x labels
            if plot_xlabels is not None:
                #set_trace()
                rax.set_xticks(np.arange(len(plot_xlabels)))
                rax.set_xticklabels(plot_xlabels)
                ax.tick_params(which="minor", bottom=False, top=False)
                rax.tick_params(which="minor", bottom=False, top=False)
                plt.setp(rax.get_xticklabels(), rotation=90, ha="left", fontsize=rcParams["font.size"]*0.5)

            ## plot y labels
            if plot_ylabels is not None: # (binlabels, bin_locs)
                for idx, label in enumerate(plot_ylabels[0]):
                    ax.annotate(label, xy=(plot_ylabels[1][idx], 0), xycoords=("data", "axes fraction"),
                        xytext=(0, 250), textcoords="offset points", va="top", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)
                        #xytext=(0, 250 if plot_xlabels is None else 120), textcoords="offset points", va="top", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)
                    ratio_ax.annotate(label, xy=(plot_ylabels[1][idx], 0), xycoords=("data", "axes fraction"),
                        xytext=(0, 250), textcoords="offset points", va="top", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)

                if hname == "mtt_vs_tlep_ctstar_abs":
                    #set_trace()
                    rax.set_xticks(mtt_bin_inds_to_plot)
                    rax.set_xticklabels(mtt_bins_to_plot)
                    ratio_ax.set_xticks(mtt_bin_inds_to_plot)
                    ratio_ax.set_xticklabels(mtt_bins_to_plot)


        ##set_trace() 
        #if hname == "Jets_njets":
        #    print(jmult)
        #    yields_txt, yields_json = Plotter.get_samples_yield_and_frac(hslice, data_lumi_year[f"{args.lepton}s"]/1000., promptmc=True)
        #    frac_name = os.path.join(pltdir, f"{jmult}_{args.lepton}_yields_and_fracs.txt")
        #    plt_tools.print_table(yields_txt, filename=frac_name, print_output=True)
        #    print(f"{frac_name} written")
        #    with open(frac_name.replace(".txt", ".json"), "w") as out:
        #        out.write(prettyjson.dumps(yields_json))

            # add lepton/jet multiplicity label
        #set_trace()
        ax.text(
            0.02, 0.92, cfeatures.channel_labels[f'{args.lepton}_{jmult}'],
            fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
        )
        hep.cms.label(ax=ax, data=withData, label="Preliminary", year=year_to_use, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))

        ratio_ax.text(
            0.02, 0.92, cfeatures.channel_labels[f'{args.lepton}_{jmult}'],
            fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ratio_ax.transAxes
        )
        hep.cms.label(ax=ratio_ax, data=withData, label="Preliminary", year=year_to_use, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))

        #set_trace()
        figname = os.path.join(pltdir, "_".join([jmult, args.lepton, hname]))
        fig.savefig(figname)
        print(f"{figname} written")
        plt.close(fig)
    
        ratio_figname = os.path.join(pltdir, "_".join([jmult, args.lepton, hname, "ratio"]))
        ratio_fig.savefig(ratio_figname)
        print(f"{ratio_figname} written")
        plt.close(ratio_fig)
        #set_trace()

toc = time.time()
print("Total time: %.1f" % (toc - tic))
