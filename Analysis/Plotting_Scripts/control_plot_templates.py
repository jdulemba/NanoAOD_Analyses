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

from coffea.hist import plot
from pdb import set_trace
import os
import Utilities.plot_tools as plt_tools
import Utilities.prettyjson as prettyjson
import re
from coffea import hist
import numpy as np
import Utilities.Plotter as Plotter
import Utilities.final_analysis_binning as final_binning
import Utilities.common_features as cfeatures
from coffea.lookup_tools.root_converters import convert_histo_root_file

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("lepton", choices=["Electron", "Muon", "Lepton"], help="Choose which lepton to make plots for")
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018", "Run2"], help="What year is the ntuple from.")
parser.add_argument("--no_est", action="store_true", help="Background estimation not performed, default is True.")
parser.add_argument("--group_tt", action="store_true", help="Group all ttbar events together")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
base_jobid = os.environ["base_jobid"]
jobid = os.environ["jobid"]
analyzer = "htt_btag_sb_regions"
eos_dir = os.environ["eos_dir"]
plots_dir = os.environ["plots_dir"]

dir_version = "BtagType3Subsources_V29"
ewqcd_str = "noEWQCDest" if args.no_est else "EWQCDest"
fname = f"run2_uncs_{ewqcd_str}_MCstat_{dir_version}.root" if args.year == "Run2" else f"splityear_uncs_{ewqcd_str}_MCstat_{dir_version}.root"
outdir = os.path.join(plots_dir, jobid, analyzer, "SysUncs", dir_version, args.year, args.lepton)

infname = os.path.join(eos_dir, "results", jobid, f"Control_Plot_Templates_{analyzer}", "Fit_Output", dir_version, fname)
rfile = convert_histo_root_file(infname)
#set_trace()

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

# binning and bin labels for phi x eta
phi_eta_binning = (
    np.array([-3.2, -2.4, -1.57, -0.87, 0., 0.87, 1.57, 2.4, 3.2]), # phi
    np.array([-2.5, -1.3, -0.7, 0., 0.7, 1.3, 2.5]) # eta binning
)
eta_binlabels = ["%s $\leq$ $\\eta$ $\leq$ %s" % (phi_eta_binning[1][bin], phi_eta_binning[1][bin+1]) for bin in range(len(phi_eta_binning[1])-1)]
eta_bin_locs = np.linspace((len(phi_eta_binning[0])-1)/2, (len(phi_eta_binning[0])-1)*(len(phi_eta_binning[1])-1) - (len(phi_eta_binning[0])-1)/2, len(phi_eta_binning[1])-1)
phi_binlabels = ["%s $\leq$ $\\phi$ $\leq$ %s" % (phi_eta_binning[0][bin], phi_eta_binning[0][bin+1]) for bin in range(len(phi_eta_binning[0])-1)]*len(eta_binlabels)


variables = {
    "Jets_njets" : (cfeatures.variable_names_to_labels["Jets_njets"], (0, 15), True),
    "mtt_vs_tlep_ctstar_abs" : (cfeatures.variable_names_to_labels["mtt"], cfeatures.variable_names_to_labels["tlep_ctstar_abs"], linearize_binning[0], linearize_binning[1],
        (linearize_binning[0][0], linearize_binning[0][-1]), (linearize_binning[1][0], linearize_binning[1][-1]), False, None, (ctstar_binlabels, ctstar_bin_locs)),
    "Jets_phi_vs_eta" : (cfeatures.variable_names_to_labels["Jets_phi"], cfeatures.variable_names_to_labels["Jets_eta"],
        phi_eta_binning[0], phi_eta_binning[1], (-3.3, 3.3), (-2.6, 2.6), True, phi_binlabels, (eta_binlabels, eta_bin_locs)),
    "Lep_phi_vs_eta" : (cfeatures.variable_names_to_labels["Lep_phi"] % cfeatures.objtypes[args.lepton], cfeatures.variable_names_to_labels["Lep_eta"] % cfeatures.objtypes[args.lepton],
        phi_eta_binning[0], phi_eta_binning[1], (-3.3, 3.3), (-2.6, 2.6), True, phi_binlabels, (eta_binlabels, eta_bin_locs)),
    "mtt" : (cfeatures.variable_names_to_labels["mtt"], (200., 2000.), False),
    "tlep_ctstar_abs" : (cfeatures.variable_names_to_labels["tlep_ctstar_abs"], (0., 1.), False),
    "mthad" : (cfeatures.variable_names_to_labels["mthad"], (0., 300.), True),
    "mWHad" : (cfeatures.variable_names_to_labels["mWHad"], (0., 300.), True),
    "mWLep" : (cfeatures.variable_names_to_labels["mWLep"], (0., 300.), True),
    "pt_thad" : (cfeatures.variable_names_to_labels["pt_thad"], (0., 500.), True),
    "pt_tlep" : (cfeatures.variable_names_to_labels["pt_tlep"], (0., 500.), True),
    "pt_tt" : (cfeatures.variable_names_to_labels["pt_tt"], (0., 500.), True),
    "eta_thad" : (cfeatures.variable_names_to_labels["eta_thad"], (-4., 4.), True),
    "eta_tlep" : (cfeatures.variable_names_to_labels["eta_tlep"], (-4., 4.), True),
    "eta_tt" : (cfeatures.variable_names_to_labels["eta_tt"], (-4., 4.), True),
    "tlep_ctstar" : (cfeatures.variable_names_to_labels["tlep_ctstar"], (-1., 1.), False),
    "full_disc" : (cfeatures.variable_names_to_labels["full_disc"], (5, 25.), True),
    "mass_disc" : (cfeatures.variable_names_to_labels["mass_disc"], (0, 20.), True),
    "ns_disc" : (cfeatures.variable_names_to_labels["ns_disc"], (3., 10.), True),
    "ns_dist" : (cfeatures.variable_names_to_labels["ns_dist"], (0., 150.), True),
    "Jets_pt" : (cfeatures.variable_names_to_labels["Jets_pt"], (0., 300.), True),
    "Jets_eta" : (cfeatures.variable_names_to_labels["Jets_eta"], (-2.6, 2.6), True),
    "Jets_phi" : (cfeatures.variable_names_to_labels["Jets_phi"], (-4., 4.), True),
    "Jets_LeadJet_pt" : (cfeatures.variable_names_to_labels["Jets_LeadJet_pt"], (0., 300.), True),
    "Jets_LeadJet_eta" : (cfeatures.variable_names_to_labels["Jets_LeadJet_eta"], (-2.6, 2.6), True),
    "Lep_pt" : (cfeatures.variable_names_to_labels["Lep_pt"] % cfeatures.objtypes[args.lepton], (0., 300.), True),
    "Lep_eta" : (cfeatures.variable_names_to_labels["Lep_eta"] % cfeatures.objtypes[args.lepton], (-2.6, 2.6), True),
    "Lep_phi" : (cfeatures.variable_names_to_labels["Lep_phi"] % cfeatures.objtypes[args.lepton], (-4, 4), True),
    "Lep_iso" : (cfeatures.variable_names_to_labels["Lep_iso"] % cfeatures.objtypes[args.lepton], (0., 1.), True),
    "MT" : (cfeatures.variable_names_to_labels["MT"], (0., 300.), True),
    "MET_pt" : (cfeatures.variable_names_to_labels["MET_pt"], (0., 300.), True),
    "MET_phi" : (cfeatures.variable_names_to_labels["MET_phi"], (-3.2, 3.2), True),
}


## get data lumi and scale MC by lumi
data_lumi_file = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())
data_lumi_year = data_lumi_file["TOT" if args.year == "Run2" else args.year]
year_label = cfeatures.year_labels["Total" if args.year == "Run2" else args.year]

if args.no_est:
    procs = ["data", "ttJets_right", "ttJets_matchable", "ttJets_unmatchable", "ttJets_sl_tau", "ttJets_other", "TB", "TW", "TQ", "EWK", "QCD"]

    if args.group_tt:
        procs_groups_dict = {
            "data" : ["data"],
            "ttJets" : ["ttJets_right", "ttJets_matchable", "ttJets_unmatchable", "ttJets_sl_tau", "ttJets_other"],
            "singlet" : ["TB", "TW", "TQ"],
            "EWK" : ["EWK"],
            "QCD" : ["QCD"],
        }
    else:
        procs_groups_dict = {
            "data" : ["data"],
            "ttJets_right" : ["ttJets_right"],
            "ttJets_matchable" : ["ttJets_matchable"],
            "ttJets_unmatchable" : ["ttJets_unmatchable"],
            "ttJets_sl_tau" : ["ttJets_sl_tau"],
            "ttJets_other" : ["ttJets_other"],
            "singlet" : ["TB", "TW", "TQ"],
            "EWK" : ["EWK"],
            "QCD" : ["QCD"],
        }
else:
    procs = ["data", "TT", "TB", "TW", "TQ", "EWQCD"]
    procs_groups_dict = {
        "data" : ["data"],
        "ttJets" : [ "TT"],
        "singlet" : ["TB", "TW", "TQ"],
        "EWQCD" : ["EWQCD"],
    }


orig_lepdir_dict = {
    "Muon" : "muNJETS",
    "Electron" : "eNJETS",
}

#set_trace()
# make plots
if args.year == "Run2":
    for hname in variables.keys():
        for jmult in ["3Jets", "4PJets"]:
            print(hname, args.year, args.lepton, jmult)
            if args.lepton == "Lepton":
                #set_trace()
                histo = None
                for lep in ["Muon", "Electron"]:
                    dirname = f"{hname}_{orig_lepdir_dict[lep].replace('NJETS', jmult.lower())}_prefit"

                    #histo = None
                    for proc in procs:
                        if proc == "data":
                            vals, bins = rfile[(f"{dirname}/{proc}", "dense_lookup")]
                            error, _ = rfile[(f"{dirname}/{proc}_error", "dense_lookup")]
                            #print(f"{proc}: {vals}")
                            if not histo:
                                histo = hist.Hist("Events", hist.Cat("dataset", "dataset"), hist.Bin("x", "x", *bins))
                                histo.fill(**{"dataset" : proc, "x": np.zeros(0), "weight" : np.zeros(0)})
                                histo.values()[(proc,)][:] = vals
                                histo.values(sumw2=True)[(proc,)][1][:] = np.square(error)
                            else:
                                tmp_histo = hist.Hist("Events", hist.Cat("dataset", "dataset"), hist.Bin("x", "x", *bins))
                                tmp_histo.fill(**{"dataset" : proc, "x": np.zeros(0), "weight" : np.zeros(0)})
                                tmp_histo.values()[(proc,)][:] = vals
                                tmp_histo.values(sumw2=True)[(proc,)][1][:] = np.square(error)
                                histo.add(tmp_histo)

                        else:
                            for year in ["2016pre", "2016post", "2017", "2018"]:
    
                                vals, bins = rfile[(f"{dirname}/{proc}_{year}", "dense_lookup")]
                                #print(f"{proc}_{year}: {vals}")
    
                                tmp_histo = hist.Hist("Events", hist.Cat("dataset", "dataset"), hist.Bin("x", "x", *bins))
                                tmp_histo.fill(**{"dataset" : proc, "x": np.zeros(0), "weight" : np.zeros(0)})
                                tmp_histo.values()[(proc,)][:] = vals

                                if histo: histo.add(tmp_histo)
            else:
                dirname = f"{hname}_{orig_lepdir_dict[args.lepton].replace('NJETS', jmult.lower())}_prefit"

                histo = None
                for proc in procs:
                    if proc == "data":
                        vals, bins = rfile[(f"{dirname}/{proc}", "dense_lookup")]
                        #print(f"{proc}: {vals}")
                        if not histo:
                            histo = hist.Hist("Events", hist.Cat("dataset", "dataset"), hist.Bin("x", "x", *bins))
                        histo.fill(**{"dataset" : proc, "x": np.zeros(0), "weight" : np.zeros(0)})
                        histo.values()[(proc,)][:] = vals
                        error, _ = rfile[(f"{dirname}/{proc}_error", "dense_lookup")]
                        histo.values(sumw2=True)[(proc,)][1][:] = np.square(error)

                    else:
                        for year in ["2016pre", "2016post", "2017", "2018"]:
    
                            vals, bins = rfile[(f"{dirname}/{proc}_{year}", "dense_lookup")]
                            #print(f"{proc}_{year}: {vals}")
    
                            tmp_histo = hist.Hist("Events", hist.Cat("dataset", "dataset"), hist.Bin("x", "x", *bins))
                            tmp_histo.fill(**{"dataset" : proc, "x": np.zeros(0), "weight" : np.zeros(0)})
                            tmp_histo.values()[(proc,)][:] = vals

                            if histo: histo.add(tmp_histo)
            #set_trace()
            histo = histo.group("dataset", hist.Cat("process", "Process", sorting="placement"), procs_groups_dict)
            
            #set_trace()
            is2d = "_vs_" in hname
            if is2d:
                xtitle, ytitle, xrebinning, yrebinning, x_lims, y_lims, withData, plot_xlabels, plot_ylabels = variables[hname]
                if (hname == "mtt_vs_tlep_ctstar_abs"):
                    withData = False
            else:
                xtitle, x_lims, withData = variables[hname]
                orig_xtitle = xtitle
    
            #set_trace()
            if args.no_est:
                mc_opts = {
                    "mcorder" : ["ttJets", "singlet", "EWK", "QCD"] if args.group_tt else ["ttJets_right", "ttJets_matchable", "ttJets_unmatchable", "ttJets_sl_tau", "ttJets_other", "singlet", "EWK", "QCD"],
                    "maskData" : not withData,
                    "overflow" : "none",
                    "ncol" : 2,
                }
            else:
                mc_opts = {
                    "mcorder" : ["ttJets", "singlet", "EWQCD"],
                    "maskData" : not withData,
                    "overflow" : "none",
                    "ncol" : 2,
                }
    
            vlines = [(len(xrebinning)-1)*ybin for ybin in range(1, len(yrebinning)-1)] if is2d else None
            nbins = (len(xrebinning)-1)*(len(yrebinning)-1) if is2d else None
    
            if hname == "Lep_iso":
                x_lims = (0., 0.1) if args.lepton == "Electron" else (0., 0.15)
                #x_lims = (0., 0.15) if args.lepton == "Muon" else (0., 0.1)
    
            if "disc" in hname:
                xtitle = orig_xtitle.replace("j", "3") if jmult == "3Jets" else orig_xtitle.replace("j", "4")
    
            #set_trace()
            fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True, figsize=(15.0, 10.0)) if is2d else plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
            fig.subplots_adjust(hspace=.07)
    
            #set_trace()
            #### plot prediction
            plot.plot1d(histo[Plotter.mc_samples], overlay=histo[Plotter.mc_samples].axes()[0].name,
                ax=ax, clear=False, stack=True, line_opts=None, fill_opts=Plotter.stack_fill_opts,
                error_opts=None, #stack_error_opts,
                order=mc_opts["mcorder"], overflow="none", overlay_overflow="none",
            )
   
            if args.lepton == "Lepton":
                #set_trace()
                pred_error = np.sqrt(sum([np.square(rfile[(f"{hname}_{orig_lepdir_dict[lep].replace('NJETS', jmult.lower())}_prefit/model_error", "dense_lookup")][0]) for lep in ["Muon", "Electron"]]))
            else:
                pred_error, _ = rfile[(f"{dirname}/model_error", "dense_lookup")]
            pred_yield = histo[Plotter.mc_samples].integrate("process").values()[()]
            ax.fill_between(
                x=histo.dense_axes()[0].edges(),
                y1=np.r_[pred_yield-pred_error, (pred_yield-pred_error)[-1]],
                y2=np.r_[pred_yield, pred_yield[-1]],
                #**{"step": "post", "label": "Total Unc.", "hatch": "///",
                **{"step": "post", "label": "Total MEscale Unc." if dir_version == "MEscaleUncs" else "Total Unc.", "hatch": "///",
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
            if withData:
                plot.plot1d(histo[Plotter.data_samples], overlay=histo[Plotter.data_samples].axes()[0].name,
                    ax=ax, clear=False, error_opts=Plotter.hstyles["data_err_opts"], overflow="none", overlay_overflow="none",
                )
    
            ax.autoscale()
            ax.set_ylim(0, ax.get_ylim()[1]*1.3)
            ax.set_xlabel(None)
            ax.set_xlim((0, nbins) if is2d else x_lims)
    
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
            #set_trace()
            
                ## plot data/MC ratio
            if not withData:
                rax.axhspan(0.5, 1.5, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
            else:
                plot.plotratio(histo[Plotter.data_samples].sum(histo[Plotter.data_samples].axes()[0].name),
                    histo[Plotter.mc_samples].sum(histo[Plotter.mc_samples].axes()[0].name),
                    ax=rax, error_opts=Plotter.hstyles["data_err_opts"],
                    denom_fill_opts=None, guide_opts={}, unc="num", overflow="none",
                )
                ##set_trace()
                sumw_denom = histo[Plotter.mc_samples].sum(histo[Plotter.mc_samples].axes()[0].name).values()[()]
                    # recompute what denom_fill_opts does
                unity = np.ones_like(sumw_denom)
                denom_unc = plot.poisson_interval(unity, np.square(pred_error) / sumw_denom**2)
                #set_trace()
                first_valid_bin, last_valid_bin = np.where(~np.isnan(denom_unc[0]))[0][0], np.where(~np.isnan(denom_unc[0]))[0][-1]+1
                low_vals = np.ma.masked_where(np.isnan(denom_unc[0][first_valid_bin:last_valid_bin]), denom_unc[0][first_valid_bin:last_valid_bin])
                hi_vals = np.ma.masked_where(np.isnan(denom_unc[1][first_valid_bin:last_valid_bin]), denom_unc[1][first_valid_bin:last_valid_bin])
                rax.fill_between(
                    histo.dense_axes()[0].edges()[first_valid_bin:last_valid_bin+1],
                    np.r_[low_vals, low_vals[-1]],
                    np.r_[hi_vals, hi_vals[-1]],
                    **{"step": "post", "facecolor": (0, 0, 0, 0.3), "linewidth": 0}
                )
    
                ## set axes labels and titles
            rax.set_ylabel("$\dfrac{data}{Pred.}$")
            rax.set_ylim(0.5, 1.5)
            rax.set_xlim((0, nbins) if is2d else x_lims)
            rax.set_xlabel(xtitle)
    
                    ## draw vertical lines for distinguishing different ctstar bins
            if vlines is not None:
                for vline in vlines:
                    ax.axvline(vline, color="k", linestyle="--")
                    if rax is not None: rax.axvline(vline, color="k", linestyle="--")
    
                # plot unrolled x and y labels for each bin
            if is2d:
                ## plot x labels
                if plot_xlabels is not None:
                    rax.set_xticks(np.arange(len(plot_xlabels)))
                    rax.set_xticklabels(plot_xlabels)
                    ax.tick_params(which="minor", bottom=False, top=False)
                    rax.tick_params(which="minor", bottom=False, top=False)
                    plt.setp(rax.get_xticklabels(), rotation=90, ha="left", fontsize=rcParams["font.size"]*0.5)
    
                if plot_ylabels is not None: # (binlabels, bin_locs)
                    for idx, label in enumerate(plot_ylabels[0]):
                        ax.annotate(label, xy=(plot_ylabels[1][idx], 0), xycoords=("data", "axes fraction"),
                            xytext=(0, 250 if plot_xlabels is None else 120), textcoords="offset points", va="top", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)
    
                    if hname == "mtt_vs_tlep_ctstar_abs":
                        #set_trace()
                        rax.set_xticks(mtt_bin_inds_to_plot)
                        rax.set_xticklabels(mtt_bins_to_plot)
    
    
                # add lepton/jet multiplicity label
            ax.text(
                0.02, 0.88, f"{cfeatures.channel_labels[f'{args.lepton}_{jmult}']}",
                fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
            )
            hep.cms.label(ax=ax, data=withData, label="Preliminary", year=year_label, lumi="138")
    
            bkg_dir = os.path.join(outdir, jmult, "btagPass") if args.no_est else os.path.join(outdir, jmult, "BKG_Est_orthog", "Sideband_Norm")
            if not os.path.isdir(bkg_dir):
                os.makedirs(bkg_dir)

            if args.no_est:
                figname = os.path.join(bkg_dir, "_".join([args.year, jmult, args.lepton, hname, "btagPass", "CombinedTT"])) if args.group_tt else os.path.join(bkg_dir, "_".join([args.year, jmult, args.lepton, hname, "btagPass"]))
            else:
                figname = os.path.join(bkg_dir, "%s_BKG_Est_orthog_Sideband_Norm" % ("_".join([args.year, jmult, args.lepton, hname])))
            fig.savefig(figname)
            print(f"{figname} written")
            plt.close(fig)
            #set_trace()

    

#set_trace()
## individual year plots
else:
    for hname in variables.keys():
        for jmult in ["3Jets", "4PJets"]:
            print(hname, args.year, args.lepton, jmult)
    
            if args.year == "2016APV": year_to_use = "2016pre"
            elif args.year == "2016": year_to_use = "2016post"
            else: year_to_use = args.year

            if args.lepton == "Lepton":
                #set_trace()
                histo = None
                for lep in ["Muon", "Electron"]:
                    dirname = f"{hname}_{orig_lepdir_dict[lep].replace('NJETS', jmult.lower())}_{year_to_use}_prefit"

                    for proc in procs:
                        vals, bins = rfile[(f"{dirname}/{proc}", "dense_lookup")]
                        error, _ = rfile[(f"{dirname}/{proc}_error", "dense_lookup")]
                        if not histo:
                            histo = hist.Hist("Events", hist.Cat("dataset", "dataset"), hist.Bin("x", "x", *bins))
                            histo.fill(**{"dataset" : proc, "x": np.zeros(0), "weight" : np.zeros(0)})
                            histo.values()[(proc,)][:] = vals
                            if proc == "data": histo.values(sumw2=True)[(proc,)][1][:] = np.square(error)
                        else:
                            tmp_histo = hist.Hist("Events", hist.Cat("dataset", "dataset"), hist.Bin("x", "x", *bins))
                            tmp_histo.fill(**{"dataset" : proc, "x": np.zeros(0), "weight" : np.zeros(0)})
                            tmp_histo.values()[(proc,)][:] = vals
                            if proc == "data": tmp_histo.values(sumw2=True)[(proc,)][1][:] = np.square(error)
                            histo.add(tmp_histo)
            else:
                dirname = f"{hname}_{orig_lepdir_dict[args.lepton].replace('NJETS', jmult.lower())}_{year_to_use}_prefit"
                histo = None
                for proc in procs:
                    vals, bins = rfile[(f"{dirname}/{proc}", "dense_lookup")]
    
                    if not histo:
                        histo = hist.Hist("Events", hist.Cat("dataset", "dataset"), hist.Bin("x", "x", *bins))
                    histo.fill(**{"dataset" : proc, "x": np.zeros(0), "weight" : np.zeros(0)})
                    histo.values()[(proc,)][:] = vals
                    if proc == "data":
                        error, _ = rfile[(f"{dirname}/{proc}_error", "dense_lookup")]
                        histo.values(sumw2=True)[(proc,)][1][:] = np.square(error)
    
            histo = histo.group("dataset", hist.Cat("process", "Process", sorting="placement"), procs_groups_dict)
            
            #set_trace()
            is2d = "_vs_" in hname
            if is2d:
                xtitle, ytitle, xrebinning, yrebinning, x_lims, y_lims, withData, plot_xlabels, plot_ylabels = variables[hname]
                if (hname == "mtt_vs_tlep_ctstar_abs"):
                    withData = False
            else:
                xtitle, x_lims, withData = variables[hname]
                orig_xtitle = xtitle
    
            #set_trace()
            if args.no_est:
                mc_opts = {
                    "mcorder" : ["ttJets", "singlet", "EWK", "QCD"] if args.group_tt else ["ttJets_right", "ttJets_matchable", "ttJets_unmatchable", "ttJets_sl_tau", "ttJets_other", "singlet", "EWK", "QCD"],
                    "maskData" : not withData,
                    "overflow" : "none",
                    "ncol" : 2,
                }
            else:
                mc_opts = {
                    "mcorder" : ["ttJets", "singlet", "EWQCD"],
                    "maskData" : not withData,
                    "overflow" : "none",
                    "ncol" : 2,
                }
    
            vlines = [(len(xrebinning)-1)*ybin for ybin in range(1, len(yrebinning)-1)] if is2d else None
            nbins = (len(xrebinning)-1)*(len(yrebinning)-1) if is2d else None
    
            if hname == "Lep_iso":
                x_lims = (0., 0.15) if args.lepton == "Muon" else (0., 0.1)
    
            if "disc" in hname:
                xtitle = orig_xtitle.replace("j", "3") if jmult == "3Jets" else orig_xtitle.replace("j", "4")
    
            #set_trace()
            fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True, figsize=(15.0, 10.0)) if is2d else plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
            fig.subplots_adjust(hspace=.07)
    
            #set_trace()
            #### plot prediction
            plot.plot1d(histo[Plotter.mc_samples], overlay=histo[Plotter.mc_samples].axes()[0].name,
                ax=ax, clear=False, stack=True, line_opts=None, fill_opts=Plotter.stack_fill_opts,
                error_opts=None, #stack_error_opts,
                order=mc_opts["mcorder"], overflow="none", overlay_overflow="none",
            )
    
            if args.lepton == "Lepton":
                #set_trace()
                pred_error = np.sqrt(sum([np.square(rfile[(f"{hname}_{orig_lepdir_dict[lep].replace('NJETS', jmult.lower())}_{year_to_use}_prefit/model_error", "dense_lookup")][0]) for lep in ["Muon", "Electron"]]))
            else:
                pred_error, _ = rfile[(f"{dirname}/model_error", "dense_lookup")]
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
            if withData:
                plot.plot1d(histo[Plotter.data_samples], overlay=histo[Plotter.data_samples].axes()[0].name,
                    ax=ax, clear=False, error_opts=Plotter.hstyles["data_err_opts"], overflow="none", overlay_overflow="none",
                )
    
            ax.autoscale()
            ax.set_ylim(0, ax.get_ylim()[1]*1.3)
            ax.set_xlabel(None)
            ax.set_xlim((0, nbins) if is2d else x_lims)
    
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
            #set_trace()
            
                ## plot data/MC ratio
            if not withData:
                rax.axhspan(0.5, 1.5, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
            else:
                plot.plotratio(histo[Plotter.data_samples].sum(histo[Plotter.data_samples].axes()[0].name),
                    histo[Plotter.mc_samples].sum(histo[Plotter.mc_samples].axes()[0].name),
                    ax=rax, error_opts=Plotter.hstyles["data_err_opts"],
                    denom_fill_opts=None, guide_opts={}, unc="num", overflow="none",
                )
                ##set_trace()
                sumw_denom = histo[Plotter.mc_samples].sum(histo[Plotter.mc_samples].axes()[0].name).values()[()]
                    # recompute what denom_fill_opts does
                unity = np.ones_like(sumw_denom)
                denom_unc = plot.poisson_interval(unity, np.square(pred_error) / sumw_denom**2)
                #set_trace()
                first_valid_bin, last_valid_bin = np.where(~np.isnan(denom_unc[0]))[0][0], np.where(~np.isnan(denom_unc[0]))[0][-1]+1
                low_vals = np.ma.masked_where(np.isnan(denom_unc[0][first_valid_bin:last_valid_bin]), denom_unc[0][first_valid_bin:last_valid_bin])
                hi_vals = np.ma.masked_where(np.isnan(denom_unc[1][first_valid_bin:last_valid_bin]), denom_unc[1][first_valid_bin:last_valid_bin])
                rax.fill_between(
                    histo.dense_axes()[0].edges()[first_valid_bin:last_valid_bin+1],
                    np.r_[low_vals, low_vals[-1]],
                    np.r_[hi_vals, hi_vals[-1]],
                    **{"step": "post", "facecolor": (0, 0, 0, 0.3), "linewidth": 0}
                )
    
                ## set axes labels and titles
            rax.set_ylabel("$\dfrac{data}{Pred.}$")
            rax.set_ylim(0.5, 1.5)
            rax.set_xlim((0, nbins) if is2d else x_lims)
            rax.set_xlabel(xtitle)
    
                    ## draw vertical lines for distinguishing different ctstar bins
            if vlines is not None:
                for vline in vlines:
                    ax.axvline(vline, color="k", linestyle="--")
                    if rax is not None: rax.axvline(vline, color="k", linestyle="--")
    
                # plot unrolled x and y labels for each bin
            if is2d:
                ## plot x labels
                if plot_xlabels is not None:
                    rax.set_xticks(np.arange(len(plot_xlabels)))
                    rax.set_xticklabels(plot_xlabels)
                    ax.tick_params(which="minor", bottom=False, top=False)
                    rax.tick_params(which="minor", bottom=False, top=False)
                    plt.setp(rax.get_xticklabels(), rotation=90, ha="left", fontsize=rcParams["font.size"]*0.5)
    
                if plot_ylabels is not None: # (binlabels, bin_locs)
                    for idx, label in enumerate(plot_ylabels[0]):
                        ax.annotate(label, xy=(plot_ylabels[1][idx], 0), xycoords=("data", "axes fraction"),
                            xytext=(0, 250 if plot_xlabels is None else 120), textcoords="offset points", va="top", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)
    
                    if hname == "mtt_vs_tlep_ctstar_abs":
                        #set_trace()
                        rax.set_xticks(mtt_bin_inds_to_plot)
                        rax.set_xticklabels(mtt_bins_to_plot)
    
    
                # add lepton/jet multiplicity label
            ax.text(
                0.02, 0.88, f"{cfeatures.channel_labels[f'{args.lepton}_{jmult}']}",
                fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
            )
            hep.cms.label(ax=ax, data=withData, label="Preliminary", year=year_label, lumi=round((data_lumi_year["Muons"]+data_lumi_year["Electrons"])/2000., 1) if args.lepton == "Lepton" else round(data_lumi_year[f"{args.lepton}s"]/1000., 1))
    
            bkg_dir = os.path.join(outdir, jmult, "btagPass") if args.no_est else os.path.join(outdir, jmult, "BKG_Est_orthog", "Sideband_Norm")
            if not os.path.isdir(bkg_dir):
                os.makedirs(bkg_dir)

            if args.no_est:
                figname = os.path.join(bkg_dir, "_".join([args.year, jmult, args.lepton, hname, "btagPass", "CombinedTT"])) if args.group_tt else os.path.join(bkg_dir, "_".join([args.year, jmult, args.lepton, hname, "btagPass"]))
            else:
                figname = os.path.join(bkg_dir, "%s_BKG_Est_orthog_Sideband_Norm" % ("_".join([args.year, jmult, args.lepton, hname])))
            fig.savefig(figname)
            print(f"{figname} written")
            plt.close(fig)
            #set_trace()

toc = time.time()
print("Total time: %.1f" % (toc - tic))
