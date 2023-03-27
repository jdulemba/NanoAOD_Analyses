#!/usr/bin/env python

# matplotlib
import matplotlib
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
from coffea.util import load#, save
from pdb import set_trace
import os
import fnmatch
import Utilities.prettyjson as prettyjson
import Utilities.plot_tools as plt_tools
import numpy as np
import Utilities.common_features as cfeatures
import Utilities.final_analysis_binning as final_binning

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2018"], help="Specify which year to run over")
#parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="Specify which year to run over")
parser.add_argument("lepton", choices=["Electron", "Muon"], help="Choose which lepton to make plots for")
parser.add_argument("--test", action="store_true", help="Use test file")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
base_jobid = os.environ["base_jobid"]
analyzer = "pileup_effect"
plot_outdir = os.environ["plots_dir"]
eos_dir = os.environ["eos_dir"]

if args.test:
    #fname = "ttJetsSL_2018_pileup_effect.only_pu_wts_test.coffea"
    fname = f"ttJetsSL_{args.year}_pileup_effect.test.coffea"
    #fname = "ttJetsSL_2018_pileup_effect.test.coffea"
    hdict = load(fname)
else:
    input_dir = os.path.join(eos_dir, "results", f"{args.year}_{jobid}", analyzer)
    f_ext = "TOT.coffea"
    fnames = [os.path.join(input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)]
    if len(fnames) > 1: raise ValueError("Too many input files found")
    hdict = load(fnames[0])

outdir = os.path.join(plot_outdir, f"{args.year}_{jobid}", analyzer)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

variables = {
    #"nPU" : ("Number of Pileup Interactions", 1, (0, 100), False),
    "nPV" : ("Number of Primary Vertices", 1, (0, 50), True),
    "rho" : ("Fastjet $\\rho$", 10, (0, 50), True),
    #"rho" : ("Fastjet $\\rho$", 1, (0, 100), True),
    #"pu_wt" : ("Pileup Weight", 1, (0, 2), False),
}

#if args.test:
#    variables.update({
#        "nTrueInt" : ("Number of True Interactions", 1, (0, 100), False),
#    })

year_to_use = cfeatures.year_labels[args.year]

#set_trace()
    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())[args.year]
lumi_correction = load(os.path.join(proj_dir, "Corrections", base_jobid, "MC_LumiWeights.coffea"))[args.year][f"{args.lepton}s"]
        # scale ttJets events, split by reconstruction type, by normal ttJets lumi correction
ttJets_permcats = ["*right", "*matchable", "*unmatchable", "*sl_tau", "*other"]
names = [dataset for dataset in sorted(set([key[0] for key in hdict[sorted(variables.keys())[0]].values().keys()]))] # get dataset names in hists
ttJets_cats = [name for name in names if any([fnmatch.fnmatch(name, cat) for cat in ttJets_permcats])] # gets ttJets(_PS)_other, ...
if len(ttJets_cats) > 0:
    for tt_cat in ttJets_cats:
        ttJets_lumi_topo = "_".join(tt_cat.split("_")[:-2]) if "sl_tau" in tt_cat else "_".join(tt_cat.split("_")[:-1]) # gets ttJets[SL, Had, DiLep] or ttJets_PS
        ttJets_eff_lumi = lumi_correction[ttJets_lumi_topo]
        lumi_correction.update({tt_cat: ttJets_eff_lumi})

## make groups based on process
process = hist.Cat("process", "Process", sorting="placement")
process_cat = "dataset"
process_groups = plt_tools.make_dataset_groups(args.lepton, args.year, samples=names, bkgdict="dataset")

#set_trace()
    # scale and group hists by process
for hname in hdict.keys():
    if "cutflow" in hname: continue
    #process_groups = {key : [key] for key in sorted(set([key[0] for key in hdict[hname].values().keys()]))}
    #set_trace()
    hdict[hname].scale(lumi_correction, axis="dataset") # scale hists
    hdict[hname] = hdict[hname].group(process_cat, process, process_groups) # group by process
    if args.test:
        hdict[hname] = hdict[hname].integrate("process")
    else:
        if "TOT" not in hname:
            hdict[hname] = hdict[hname][:, :, :, args.lepton].integrate("leptype") # only pick out specified lepton

#set_trace()

## plot histograms
if args.test:
    for hname in variables.keys():
        if hname not in hdict.keys():
            raise ValueError(f"{hname} not found in file")

        print(hname)    
        histo = hdict[hname].copy()
        tot_histo = hdict[f"TOT_{hname}"].copy() if f"TOT_{hname}" in hdict.keys() else None
    
        xtitle, rebinning, x_lims, withData = variables[hname]
        xaxis_name = histo.dense_axes()[0].name
        if rebinning != 1:
            histo = histo.rebin(xaxis_name, rebinning)
            if tot_histo: tot_histo = tot_histo.rebin(xaxis_name, rebinning)
    
        #set_trace()
        if tot_histo and (hname != "pu_wt"):
            pltdir = os.path.join(outdir, "TOT")
            if not os.path.isdir(pltdir):
                os.makedirs(pltdir)

                # plot expected yields and ratios for total events
            fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
            fig.subplots_adjust(hspace=.07)

            hep.plot.histplot(tot_histo["nosys"].integrate("sys").values()[()], tot_histo.axis(xaxis_name).edges(), ax=ax, histtype="step", **{"linestyle" : "-", "color" : "k", "label" : "Central"})
            hep.plot.histplot(tot_histo["PileupUp"].integrate("sys").values()[()], tot_histo.axis(xaxis_name).edges(), ax=ax, histtype="step", **{"linestyle" : "-", "color" : "r", "label" : "Up"})
            hep.plot.histplot(tot_histo["PileupDown"].integrate("sys").values()[()], tot_histo.axis(xaxis_name).edges(), ax=ax, histtype="step", **{"linestyle" : "-", "color" : "b", "label" : "Down"})
            hep.plot.histplot(tot_histo["No_Pileup"].integrate("sys").values()[()], tot_histo.axis(xaxis_name).edges(), ax=ax, histtype="step", **{"linestyle" : "-", "color" : "g", "label" : "Uncorrected"})

            ## ratios of correction/uncorrected
            cen_masked_vals, cen_masked_bins = Plotter.get_ratio_arrays(num_vals=tot_histo["nosys"].integrate("sys").values()[()],
                denom_vals=tot_histo["No_Pileup"].integrate("sys").values()[()], input_bins=tot_histo.axis(xaxis_name).edges())
            up_masked_vals, up_masked_bins = Plotter.get_ratio_arrays(num_vals=tot_histo["PileupUp"].integrate("sys").values()[()],
                denom_vals=tot_histo["No_Pileup"].integrate("sys").values()[()], input_bins=tot_histo.axis(xaxis_name).edges())
            dw_masked_vals, dw_masked_bins = Plotter.get_ratio_arrays(num_vals=tot_histo["PileupDown"].integrate("sys").values()[()],
                denom_vals=tot_histo["No_Pileup"].integrate("sys").values()[()], input_bins=tot_histo.axis(xaxis_name).edges())

            rax.step(cen_masked_bins, cen_masked_vals, where="post", **{"linestyle" : "-", "color" : "k", "label" : "Central"})
            rax.step(up_masked_bins, up_masked_vals, where="post", **{"linestyle" : "-", "color" : "r", "label" : "Up"})
            rax.step(dw_masked_bins, dw_masked_vals, where="post", **{"linestyle" : "-", "color" : "b", "label" : "Down"})

            ax.legend(ncol=2)
            ax.autoscale()
            ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.3)
            ax.set_xlim(x_lims)
            ax.set_xlabel(None)
            ax.set_ylabel("Events")

            rax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
            rax.autoscale()
            rax.set_xlim(x_lims)
            rax.set_xlabel(xtitle)
            rax.set_ylabel("Ratio to Uncorrected")
            rax.set_ylim(0.5, 1.5)

            hep.cms.label(ax=ax, data=False, year=year_to_use, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))

            figname = os.path.join(pltdir, f"New_{hname}_only_pu_wts") if "only_pu_wts" in fname else os.path.join(pltdir, f"New_{hname}")
            #figname = os.path.join(pltdir, f"{hname}_only_pu_wts") if "only_pu_wts" in fname else os.path.join(pltdir, hname)
            fig.savefig(figname)
            print(f"{figname} written")
            plt.close(fig)

        ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
        for jmult in sorted(set([key[1] for key in histo.values().keys()])):
            pltdir = os.path.join(outdir, args.lepton, jmult)
            if not os.path.isdir(pltdir):
                os.makedirs(pltdir)
            
            print(f"\t{jmult}")
    
            #set_trace()
            channel = histo[:, jmult, args.lepton].integrate("jmult").integrate("leptype")
    
                # plot normalized shapes
            fig, ax = plt.subplots()
            fig.subplots_adjust(hspace=.07)

            Plotter.plot_1D(channel["nosys"].integrate("sys").values()[()]/np.sum(channel["nosys"].integrate("sys").values()[()]), channel.axis(xaxis_name).edges(),
                xlimits=x_lims, xlabel=xtitle, ylabel="Probability Density", color="k", ax=ax, label=f"{cfeatures.channel_labels[f'{args.lepton}_{jmult}']} Central")
            Plotter.plot_1D(channel["PileupUp"].integrate("sys").values()[()]/np.sum(channel["PileupUp"].integrate("sys").values()[()]), channel.axis(xaxis_name).edges(),
                xlimits=x_lims, xlabel=xtitle, ylabel="Probability Density", color="r", ax=ax, label=f"{cfeatures.channel_labels[f'{args.lepton}_{jmult}']} Up")
            Plotter.plot_1D(channel["PileupDown"].integrate("sys").values()[()]/np.sum(channel["PileupDown"].integrate("sys").values()[()]), channel.axis(xaxis_name).edges(),
                xlimits=x_lims, xlabel=xtitle, ylabel="Probability Density", color="b", ax=ax, label=f"{cfeatures.channel_labels[f'{args.lepton}_{jmult}']} Down")

            if hname != "pu_wt":
                #set_trace()
                Plotter.plot_1D(channel["No_Pileup"].integrate("sys").values()[()]/np.sum(channel["No_Pileup"].integrate("sys").values()[()]), channel.axis(xaxis_name).edges(),
                    xlimits=x_lims, xlabel=xtitle, ylabel="Probability Density", color="g", ax=ax, label=f"{cfeatures.channel_labels[f'{args.lepton}_{jmult}']} Uncorrected")

            if tot_histo:
                Plotter.plot_1D(tot_histo["nosys"].integrate("sys").values()[()]/np.sum(tot_histo["nosys"].integrate("sys").values()[()]), tot_histo.axis(xaxis_name).edges(),
                    xlimits=x_lims, xlabel=xtitle, ylabel="Probability Density", linestyle="--", color="k", ax=ax, label="Central")
                Plotter.plot_1D(tot_histo["PileupUp"].integrate("sys").values()[()]/np.sum(tot_histo["PileupUp"].integrate("sys").values()[()]), tot_histo.axis(xaxis_name).edges(),
                    xlimits=x_lims, xlabel=xtitle, ylabel="Probability Density", linestyle="--", color="r", ax=ax, label="Up")
                Plotter.plot_1D(tot_histo["PileupDown"].integrate("sys").values()[()]/np.sum(tot_histo["PileupDown"].integrate("sys").values()[()]), tot_histo.axis(xaxis_name).edges(),
                    xlimits=x_lims, xlabel=xtitle, ylabel="Probability Density", linestyle="--", color="b", ax=ax, label="Down")
                if hname != "pu_wt":
                    Plotter.plot_1D(tot_histo["No_Pileup"].integrate("sys").values()[()]/np.sum(tot_histo["No_Pileup"].integrate("sys").values()[()]), tot_histo.axis(xaxis_name).edges(),
                        xlimits=x_lims, xlabel=xtitle, ylabel="Probability Density", linestyle="--", color="g", ax=ax, label="Uncorrected")

            ax.legend()

            figname = os.path.join(pltdir, "_".join(["New", args.lepton, jmult, hname, "NORM", "only_pu_wts"])) if "only_pu_wts" in fname else os.path.join(pltdir, "_".join(["New", args.lepton, jmult, hname, "NORM"]))
            #figname = os.path.join(pltdir, "_".join([args.lepton, jmult, hname, "NORM", "only_pu_wts"])) if "only_pu_wts" in fname else os.path.join(pltdir, "_".join([args.lepton, jmult, hname, "NORM"]))
            fig.savefig(figname)
            print(f"{figname} written")
            plt.close(fig)

            #set_trace()
            if hname != "pu_wt":
                    # plot expected yields and ratios for total events
                fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                fig.subplots_adjust(hspace=.07)

                hep.plot.histplot(channel["nosys"].integrate("sys").values()[()], channel.axis(xaxis_name).edges(), ax=ax, histtype="step",
                    **{"linestyle" : "-", "color" : "k", "label" : "Central"})
                hep.plot.histplot(channel["PileupUp"].integrate("sys").values()[()], channel.axis(xaxis_name).edges(), ax=ax, histtype="step",
                    **{"linestyle" : "-", "color" : "r", "label" : "Up"})
                hep.plot.histplot(channel["PileupDown"].integrate("sys").values()[()], channel.axis(xaxis_name).edges(), ax=ax, histtype="step",
                    **{"linestyle" : "-", "color" : "b", "label" : "Down"})
                hep.plot.histplot(channel["No_Pileup"].integrate("sys").values()[()], channel.axis(xaxis_name).edges(), ax=ax, histtype="step",
                    **{"linestyle" : "-", "color" : "g", "label" : "Uncorrected"})

                ## ratios of correction/uncorrected
                cen_masked_vals, cen_masked_bins = Plotter.get_ratio_arrays(num_vals=channel["nosys"].integrate("sys").values()[()],
                    denom_vals=channel["No_Pileup"].integrate("sys").values()[()], input_bins=channel.axis(xaxis_name).edges())
                up_masked_vals, up_masked_bins = Plotter.get_ratio_arrays(num_vals=channel["PileupUp"].integrate("sys").values()[()],
                    denom_vals=channel["No_Pileup"].integrate("sys").values()[()], input_bins=channel.axis(xaxis_name).edges())
                dw_masked_vals, dw_masked_bins = Plotter.get_ratio_arrays(num_vals=channel["PileupDown"].integrate("sys").values()[()],
                    denom_vals=channel["No_Pileup"].integrate("sys").values()[()], input_bins=channel.axis(xaxis_name).edges())

                rax.step(cen_masked_bins, cen_masked_vals, where="post", **{"linestyle" : "-", "color" : "k", "label" : "Central"})
                rax.step(up_masked_bins, up_masked_vals, where="post", **{"linestyle" : "-", "color" : "r", "label" : "Up"})
                rax.step(dw_masked_bins, dw_masked_vals, where="post", **{"linestyle" : "-", "color" : "b", "label" : "Down"})

                ax.legend(ncol=2)
                ax.autoscale()
                ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.3)
                ax.set_xlim(x_lims)
                ax.set_xlabel(None)
                ax.set_ylabel("Events")

                rax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
                rax.autoscale()
                rax.set_xlim(x_lims)
                rax.set_xlabel(xtitle)
                rax.set_ylabel("Ratio to Uncorrected")
                rax.set_ylim(0.5, 1.5)

                ax.text(
                    0.02, 0.90, f"{cfeatures.channel_labels[f'{args.lepton}_{jmult}']}",
                    fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                )
                hep.cms.label(ax=ax, data=False, year=year_to_use, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))

                figname = os.path.join(pltdir, "_".join(["New", args.lepton, jmult, hname, "yields", "only_pu_wts"])) if "only_pu_wts" in fname else os.path.join(pltdir, "_".join(["New", args.lepton, jmult, hname, "yields"]))
                #figname = os.path.join(pltdir, "_".join([args.lepton, jmult, hname, "yields", "only_pu_wts"])) if "only_pu_wts" in fname else os.path.join(pltdir, "_".join([args.lepton, jmult, hname, "yields"]))
                fig.savefig(figname)
                print(f"{figname} written")
                plt.close(fig)
    
else:
    for hname in variables.keys():
        if hname not in hdict.keys():
            raise ValueError(f"{hname} not found in file")
    
        print(hname)    
        histo = hdict[hname].copy()
    
        xtitle, rebinning, x_lims, withData = variables[hname]
        xaxis_name = histo.dense_axes()[0].name
        if rebinning != 1:
            histo = histo.rebin(xaxis_name, rebinning)
    
        ##set_trace()
        #nosys_histo = hdict[hname][:, "nosys"].integrate("sys")
        #pu_up_histo = hdict[hname][Plotter.mc_samples, "PileupUp"].integrate("process").integrate("sys")
        #pu_dw_histo = hdict[hname][Plotter.mc_samples, "PileupDown"].integrate("process").integrate("sys")
        #uncorr_histo= hdict[hname][Plotter.mc_samples, "No_Pileup"].integrate("process").integrate("sys")
        ##set_trace()
    
        #xtitle, rebinning, x_lims, withData = variables[hname]
        #if rebinning != 1:
        #    xaxis_name = nosys_histo.dense_axes()[0].name
        #    nosys_histo = nosys_histo.rebin(xaxis_name, rebinning)
        #    pu_up_histo = pu_up_histo.rebin(xaxis_name, rebinning)
        #    pu_dw_histo = pu_dw_histo.rebin(xaxis_name, rebinning)
        #    uncorr_histo = uncorr_histo.rebin(xaxis_name, rebinning)
    
        ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
        #for jmult in sorted(set([key[1] for key in nosys_histo.values().keys()])):
        for jmult in sorted(set([key[2] for key in histo.values().keys()])):
            pltdir = os.path.join(outdir, args.lepton, jmult)
            if not os.path.isdir(pltdir):
                os.makedirs(pltdir)
            
            print(f"\t{jmult}")
    
            ##set_trace()        
            #nosys_hslice = nosys_histo[:, jmult].integrate("jmult")
    
            #fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
            #fig.subplots_adjust(hspace=.07)
    
            ### plot nosys MC
            #hist.plot.plot1d(nosys_hslice[Plotter.mc_samples], overlay=nosys_hslice.axes()[0].name,
            #    ax=ax, clear=False, stack=True, line_opts=None, fill_opts=Plotter.stack_fill_opts,
            #    error_opts=Plotter.stack_error_opts,
            #    overflow="none", overlay_overflow="none",
            #)
            ### plot data
            #hist.plot.plot1d(nosys_hslice[Plotter.data_samples], overlay=nosys_hslice[Plotter.data_samples].axes()[0].name,
            #    ax=ax, clear=False, error_opts=Plotter.hstyles["data_err_opts"],
            #    overflow="none", overlay_overflow="none",
            #)
            ### plot uncorrected MC
            #hist.plot.plot1d(uncorr_histo[jmult].integrate("jmult"),
            #    ax=ax, clear=False, stack=False, fill_opts=None, error_opts=None,
            #    line_opts={"linestyle" : "-", "color" : "r"},
            #    overflow="none", overlay_overflow="none",
            #)
            #### plot pileup up variation
            ##hist.plot.plot1d(pu_up_histo[jmult].integrate("jmult"),
            ##    ax=ax, clear=False, stack=False, fill_opts=None, error_opts=None,
            ##    line_opts={"linestyle" : "-", "color" : "r"},
            ##    overflow="none", overlay_overflow="none",
            ##)
            #### plot pileup down variation
            ##hist.plot.plot1d(pu_dw_histo[jmult].integrate("jmult"),
            ##    ax=ax, clear=False, stack=False, fill_opts=None, error_opts=None,
            ##    line_opts={"linestyle" : "-", "color" : "b"},
            ##    overflow="none", overlay_overflow="none",
            ##)
    
            ##set_trace()
            #ax.autoscale()
            #ax.set_ylim(0, ax.get_ylim()[1]*1.3)
            #ax.set_xlabel(None)
            #ax.set_xlim(x_lims)
    
            #    ## set legend and corresponding colors
            #handles, labels = ax.get_legend_handles_labels()
            ##set_trace()
            #for idx, sample in enumerate(labels):
            #    if sample == "data" or sample == "Observed": continue
            #    if isinstance(handles[idx], matplotlib.lines.Line2D): continue
            #    if sample == "None": labels[idx] = "Uncorrected Sim."
            #    else:
            #        facecolor, legname = plt_tools.get_styles(sample, Plotter.hstyles)
            #        handles[idx].set_facecolor(facecolor)
            #        labels[idx] = legname
            ## call ax.legend() with the new values
            #ax.legend(handles,labels, loc="upper right", ncol=2)
    
            ### plot data/central MC ratio
            #hist.plot.plotratio(nosys_hslice[Plotter.data_samples].sum(nosys_hslice[Plotter.data_samples].axes()[0].name),
            #    nosys_hslice[Plotter.mc_samples].sum(nosys_hslice[Plotter.mc_samples].axes()[0].name),
            #    ax=rax, error_opts=Plotter.hstyles["data_err_opts"], denom_fill_opts={},
            #    guide_opts={}, unc="num", overflow="none",
            #)
    
            ### plot data/uncorrected central MC ratio
            ##set_trace()
            #hist.plot.plotratio(nosys_hslice[Plotter.data_samples].sum(nosys_hslice[Plotter.data_samples].axes()[0].name),
            #    uncorr_histo[jmult].integrate("jmult"),
            #    ax=rax, error_opts={"linestyle" : "-", "color" : "r"}, denom_fill_opts={},
            #    clear=False,
            #    #ax=rax, error_opts=Plotter.hstyles["data_err_opts"], denom_fill_opts={},
            #    guide_opts={}, unc="num", overflow="none",
            #)
    
            ##### plot data/pileup variations ratio
            ##up_masked_vals, up_masked_bins = Plotter.get_ratio_arrays(num_vals=nosys_hslice[Plotter.data_samples].sum(nosys_hslice[Plotter.data_samples].axes()[0].name).values()[()],
            ##    denom_vals=pu_up_histo[jmult].integrate("jmult").values()[()], input_bins=nosys_hslice[Plotter.data_samples].dense_axes()[0].edges())
            ##dw_masked_vals, dw_masked_bins = Plotter.get_ratio_arrays(num_vals=nosys_hslice[Plotter.data_samples].sum(nosys_hslice[Plotter.data_samples].axes()[0].name).values()[()],
            ##    denom_vals=pu_dw_histo[jmult].integrate("jmult").values()[()], input_bins=nosys_hslice[Plotter.data_samples].dense_axes()[0].edges())
    
            ##rax.step(up_masked_bins, up_masked_vals, where="post", label="Pileup Up", color="r")
            ##rax.step(dw_masked_bins, dw_masked_vals, where="post", label="Pileup Down", color="b")
            ###rax.fill_between(dw_masked_bins, dw_masked_vals, y2=up_masked_vals, facecolor="r", label="Pileup", step="post", alpha=0.5)
            ####rax.fill_between(dw_masked_bins, dw_masked_vals, y2=1., facecolor=dw_style["color"], step="post", alpha=0.5)
    
            #    ## set axes labels and titles
            #rax.legend(loc="upper right")
            #rax.set_ylabel("data/MC")
            #rax.set_ylim(0.7, 1.2)
            ##rax.set_ylim(0.5, 1.5)
            #rax.set_xlim(x_lims)
            #rax.set_xlabel(xtitle)
    
            #    # add lepton/jet multiplicity label
            #ax.text(
            #    0.02, 0.88, f"{cfeatures.channel_labels[f'{args.lepton}_{jmult}']}",
            #    fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
            #)
            #hep.cms.label(ax=ax, label="Preliminary", data=withData, year=year_to_use, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))
            #
            #figname = os.path.join(pltdir, "_".join([args.lepton, jmult, hname]))
            #fig.savefig(figname)
            #print(f"{figname} written")
            #plt.close(fig)

                # plot expected yields and ratios for total events
            fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
            fig.subplots_adjust(hspace=.07)

            ## plot total simulation
            hep.plot.histplot(histo[Plotter.mc_samples, "nosys", jmult].integrate("process").integrate("sys").integrate("jmult").values()[()],
                histo.axis(xaxis_name).edges(), ax=ax, histtype="step", **{"linestyle" : "-", "color" : "k", "label" : "Sim., Central"})
            hep.plot.histplot(histo[Plotter.mc_samples, "PileupUp", jmult].integrate("process").integrate("sys").integrate("jmult").values()[()],
                histo.axis(xaxis_name).edges(), ax=ax, histtype="step", **{"linestyle" : "-", "color" : "r", "label" : "Sim., Up"})
            hep.plot.histplot(histo[Plotter.mc_samples, "PileupDown", jmult].integrate("process").integrate("sys").integrate("jmult").values()[()],
                histo.axis(xaxis_name).edges(), ax=ax, histtype="step", **{"linestyle" : "-", "color" : "b", "label" : "Sim., Down"})
            hep.plot.histplot(histo[Plotter.mc_samples, "No_Pileup", jmult].integrate("process").integrate("sys").integrate("jmult").values()[()],
                histo.axis(xaxis_name).edges(), ax=ax, histtype="step", **{"linestyle" : "-", "color" : "g", "label" : "Sim., Uncorrected"})
            ## plot data
            hist.plot.plot1d(histo[Plotter.data_samples, :, jmult].integrate("sys").integrate("jmult"), overlay=histo.axes()[0].name,
                ax=ax, clear=False, error_opts=Plotter.hstyles["data_err_opts"], overflow="none", overlay_overflow="none",
            )

            #set_trace()
            ## ratios of correction/uncorrected
            cen_masked_vals, cen_masked_bins = Plotter.get_ratio_arrays(num_vals=histo[Plotter.mc_samples, "nosys", jmult].integrate("process").integrate("sys").integrate("jmult").values()[()],
                denom_vals=histo[Plotter.mc_samples, "No_Pileup", jmult].integrate("process").integrate("sys").integrate("jmult").values()[()], input_bins=histo.axis(xaxis_name).edges())
            up_masked_vals, up_masked_bins = Plotter.get_ratio_arrays(num_vals=histo[Plotter.mc_samples, "PileupUp", jmult].integrate("process").integrate("sys").integrate("jmult").values()[()],
                denom_vals=histo[Plotter.mc_samples, "No_Pileup", jmult].integrate("process").integrate("sys").integrate("jmult").values()[()], input_bins=histo.axis(xaxis_name).edges())
            dw_masked_vals, dw_masked_bins = Plotter.get_ratio_arrays(num_vals=histo[Plotter.mc_samples, "PileupDown", jmult].integrate("process").integrate("sys").integrate("jmult").values()[()],
                denom_vals=histo[Plotter.mc_samples, "No_Pileup", jmult].integrate("process").integrate("sys").integrate("jmult").values()[()], input_bins=histo.axis(xaxis_name).edges())

            rax.step(cen_masked_bins, cen_masked_vals, where="post", **{"linestyle" : "-", "color" : "k", "label" : "Central"})
            rax.step(up_masked_bins, up_masked_vals, where="post", **{"linestyle" : "-", "color" : "r", "label" : "Up"})
            rax.step(dw_masked_bins, dw_masked_vals, where="post", **{"linestyle" : "-", "color" : "b", "label" : "Down"})

            ax.legend(ncol=2)
            ax.autoscale()
            ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.3)
            ax.set_xlim(x_lims)
            ax.set_xlabel(None)
            ax.set_ylabel("Events")

            rax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
            rax.autoscale()
            rax.set_xlim(x_lims)
            rax.set_xlabel(xtitle)
            rax.set_ylabel("Ratio to Uncorrected")
            rax.set_ylim(0.5, 1.5)

            ax.text(
                0.02, 0.90, f"{cfeatures.channel_labels[f'{args.lepton}_{jmult}']}",
                fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
            )
            hep.cms.label(ax=ax, data=True, label="Preliminary", year=year_to_use, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))

            #figname = os.path.join(pltdir, "_".join([args.lepton, jmult, hname, "yields", "only_pu_wts"])) if "only_pu_wts" in fname else os.path.join(pltdir, "_".join([args.lepton, jmult, hname, "yields"]))
            figname = os.path.join(pltdir, "_".join(["New", args.lepton, jmult, hname, "AllSim_Data_Comp"]))
            #figname = os.path.join(pltdir, "_".join([args.lepton, jmult, hname, "AllSim_Data_Comp"]))
            fig.savefig(figname)
            print(f"{figname} written")
            plt.close(fig)
    

                # plot expected yields and ratios for individual simulation
            for proc in sorted(set([key[0] for key in histo[Plotter.mc_samples].values().keys()])):
                #if "ttJets" not in proc: continue
                print(hname, jmult, proc)

                fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                fig.subplots_adjust(hspace=.07)

                ## plot total simulation
                hep.plot.histplot(histo[proc, "nosys", jmult].integrate("process").integrate("sys").integrate("jmult").values()[()],
                    histo.axis(xaxis_name).edges(), ax=ax, histtype="step", **{"linestyle" : "-", "color" : "k", "label" : "Sim., Central"})
                hep.plot.histplot(histo[proc, "PileupUp", jmult].integrate("process").integrate("sys").integrate("jmult").values()[()],
                    histo.axis(xaxis_name).edges(), ax=ax, histtype="step", **{"linestyle" : "-", "color" : "r", "label" : "Sim., Up"})
                hep.plot.histplot(histo[proc, "PileupDown", jmult].integrate("process").integrate("sys").integrate("jmult").values()[()],
                    histo.axis(xaxis_name).edges(), ax=ax, histtype="step", **{"linestyle" : "-", "color" : "b", "label" : "Sim., Down"})
                hep.plot.histplot(histo[proc, "No_Pileup", jmult].integrate("process").integrate("sys").integrate("jmult").values()[()],
                    histo.axis(xaxis_name).edges(), ax=ax, histtype="step", **{"linestyle" : "-", "color" : "g", "label" : "Sim., Uncorrected"})
                ### plot data
                #hist.plot.plot1d(histo[Plotter.data_samples, :, jmult].integrate("sys").integrate("jmult"), overlay=histo.axes()[0].name,
                #    ax=ax, clear=False, error_opts=Plotter.hstyles["data_err_opts"], overflow="none", overlay_overflow="none",
                #)

                #set_trace()
                ## ratios of correction/uncorrected
                cen_masked_vals, cen_masked_bins = Plotter.get_ratio_arrays(num_vals=histo[proc, "nosys", jmult].integrate("process").integrate("sys").integrate("jmult").values()[()],
                    denom_vals=histo[proc, "No_Pileup", jmult].integrate("process").integrate("sys").integrate("jmult").values()[()], input_bins=histo.axis(xaxis_name).edges())
                up_masked_vals, up_masked_bins = Plotter.get_ratio_arrays(num_vals=histo[proc, "PileupUp", jmult].integrate("process").integrate("sys").integrate("jmult").values()[()],
                    denom_vals=histo[proc, "No_Pileup", jmult].integrate("process").integrate("sys").integrate("jmult").values()[()], input_bins=histo.axis(xaxis_name).edges())
                dw_masked_vals, dw_masked_bins = Plotter.get_ratio_arrays(num_vals=histo[proc, "PileupDown", jmult].integrate("process").integrate("sys").integrate("jmult").values()[()],
                    denom_vals=histo[proc, "No_Pileup", jmult].integrate("process").integrate("sys").integrate("jmult").values()[()], input_bins=histo.axis(xaxis_name).edges())

                rax.step(cen_masked_bins, cen_masked_vals, where="post", **{"linestyle" : "-", "color" : "k", "label" : "Central"})
                rax.step(up_masked_bins, up_masked_vals, where="post", **{"linestyle" : "-", "color" : "r", "label" : "Up"})
                rax.step(dw_masked_bins, dw_masked_vals, where="post", **{"linestyle" : "-", "color" : "b", "label" : "Down"})

                ax.legend(ncol=2, title=proc)
                ax.autoscale()
                ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.3)
                ax.set_xlim(x_lims)
                ax.set_xlabel(None)
                ax.set_ylabel("Events")

                rax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
                rax.autoscale()
                rax.set_xlim(x_lims)
                rax.set_xlabel(xtitle)
                rax.set_ylabel("Ratio to Uncorrected")
                rax.set_ylim(0.5, 1.5)

                ax.text(
                    0.02, 0.90, f"{cfeatures.channel_labels[f'{args.lepton}_{jmult}']}",
                    fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                )
                hep.cms.label(ax=ax, data=False, label="Preliminary", year=year_to_use, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))

                figname = os.path.join(pltdir, "_".join(["New", args.lepton, jmult, hname, proc]))
                #figname = os.path.join(pltdir, "_".join([args.lepton, jmult, hname, proc]))
                fig.savefig(figname)
                print(f"{figname} written")
                plt.close(fig)
    

