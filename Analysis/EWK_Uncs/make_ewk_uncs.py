#!/usr/bin/env python

# matplotlib
import matplotlib as mpl
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
from pdb import set_trace
import os
import numpy as np
from coffea.lookup_tools.root_converters import convert_histo_root_file
from coffea.lookup_tools.dense_lookup import dense_lookup
from itertools import product
import Utilities.common_features as cfeatures

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("uncs", choices=["EW", "NNLO", "All"], help="Which EW uncertainty to calculate.")
parser.add_argument("--save_ratios", action="store_true", help="Save NNLO/tune ratios for all years (default is to not save them).")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
base_jobid = os.environ["base_jobid"]
base_dir = "EWK_Uncs"
plots_dir = os.environ["plots_dir"]
eos_dir = os.environ["eos_dir"]

save_dict = {}
    
fname = "HATHORscaled.root"
ewk_dists = convert_histo_root_file(os.path.join(proj_dir, base_dir, fname))

# LO QCD dist is same for all directories
lo_xsec_vals, (tot_mtt_binning, tot_ctstar_binning) = ewk_dists[("EWno/tot", "dense_lookup")]
rebinned_lo_xsec_vals, (rebinned_ctstar_binning, rebinned_mtt_binning) = ewk_dists[("EWno/tot_rebinned", "dense_lookup")]
lo_ztitle = "$\\sigma_{LO\ QCD}$ [pb]"

max_mtt_val = 2000.
        # histo plotting params
mtt_title, ctstar_title = cfeatures.variable_names_to_labels["mtt"], cfeatures.variable_names_to_labels["ctstar"]
tot_mtt_lims = (np.min(tot_mtt_binning), min(np.max(tot_mtt_binning), max_mtt_val))
tot_ctstar_lims = (np.min(tot_ctstar_binning), np.max(tot_ctstar_binning))
rebinned_mtt_lims = tot_mtt_lims
#rebinned_mtt_lims = (np.min(rebinned_mtt_binning), min(np.max(rebinned_mtt_binning), max_mtt_val))
rebinned_ctstar_lims = (np.min(rebinned_ctstar_binning), np.max(rebinned_ctstar_binning))
#set_trace()
    
if (args.uncs == "EW") or (args.uncs == "All"):
    outdir = os.path.join(plots_dir, base_dir, "EW")
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
       
    yt_vals_dict = {
        "0" : 0.0,
        "5" : 0.5,
        "9" : 0.9,
        "10" : 1.0,
        "11" : 1.1,
        "15" : 1.5,
        "20" : 2.0,
        "30" : 3.0,
        "40" : 4.0,
        "111" : 1.11,
        "088" : 0.88,
    }
        
    #nameTOval = lambda val : str(val).replace('p', '.')
    nlo_ewk_xsec_base_label = "$\\sigma_{LO\ QCD\ +\ NLO\ EW}$ [pb], $y^{t}/y^{t}_{SM}$ = YT"
    kfactor_base_label = "$K^{YT}_{NLO\ EW} \\equiv \dfrac{\\sigma_{LO\ QCD\ +\ NLO\ EW}}{\\sigma_{LO\ QCD}}$, $y^{t}/y^{t}_{SM}$ = YT"
    deltaEW_base_label = "$\\delta^{YT}_{EW} \\equiv \dfrac{\\sigma_{LO\ QCD\ +\ NLO\ EW} - \\sigma_{LO\ QCD}}{\\sigma_{LO\ QCD}}$, $y^{t}/y^{t}_{SM}$ = YT"

    yt1_nlo_ewk_xsec_base_label = "$\\sigma_{LO\ QCD\ +\ NLO\ EW}$ [pb]"
    yt1_kfactor_base_label = "$K_{NLO\ EW} \\equiv \dfrac{\\sigma_{LO\ QCD\ +\ NLO\ EW}}{\\sigma_{LO\ QCD}}$"
    yt1_deltaEW_base_label = "$\\delta_{EW} \\equiv \dfrac{\\sigma_{LO\ QCD\ +\ NLO\ EW} - \\sigma_{LO\ QCD}}{\\sigma_{LO\ QCD}}$"

            # only get mtt values < 2000 GeV
    masked_lo_vals = lo_xsec_vals[np.where(tot_mtt_binning <= max_mtt_val)[0][0]:np.where(tot_mtt_binning <= max_mtt_val)[0][-1], :]
    masked_tot_mtt_binning = tot_mtt_binning[np.where(tot_mtt_binning <= max_mtt_val)[0]]
    masked_rebinned_lo_vals = rebinned_lo_xsec_vals[:, np.where(rebinned_mtt_binning <= max_mtt_val)[0][0]:np.where(rebinned_mtt_binning <= max_mtt_val)[0][-1]]
    masked_rebinned_mtt_binning = rebinned_mtt_binning[np.where(rebinned_mtt_binning <= max_mtt_val)[0]]

            # get binning for differential plots
    max_mtt_val_bin = np.argwhere(masked_rebinned_mtt_binning == max_mtt_val)[0][0]
    mtt_binwidths = np.array([masked_rebinned_mtt_binning[i+1] - masked_rebinned_mtt_binning[i] for i in range(max_mtt_val_bin)])
    ctstar_binwidths = np.array([rebinned_ctstar_binning[i+1] - rebinned_ctstar_binning[i] for i in range(len(rebinned_ctstar_binning) - 1)])
    tiled_mtt_binwidths = np.tile(mtt_binwidths, (ctstar_binwidths.size, 1)).T
    tiled_ctstar_binwidths = np.tile(ctstar_binwidths, (mtt_binwidths.size, 1))
    #set_trace()

        ### plot original LO QCD hist
    fig, ax = plt.subplots()
    fig.subplots_adjust(hspace=.07)
    Plotter.plot_2d_norm(values=np.ma.masked_where(masked_lo_vals <= 0.0, masked_lo_vals), # mask nonzero probabilities for plotting
        xbins=masked_tot_mtt_binning, ybins=tot_ctstar_binning,
        xlimits=tot_mtt_lims, ylimits=tot_ctstar_lims, xlabel=mtt_title, ylabel=ctstar_title,
        ax=ax, **{"cmap_label" : lo_ztitle})
    ax.text(
        0.98, 0.90, "$t\\bar{t}$\nparton level",
        fontsize=rcParams["font.size"], horizontalalignment="right", verticalalignment="bottom", transform=ax.transAxes, color="w",
    )
    hep.label.exp_label(ax=ax, exp="MADGRAPH", rlabel="")
    figname = os.path.join(plots_dir, base_dir, f"LO_QCD_xsec_mtt_vs_ctstar")
    fig.savefig(figname)
    print(f"{figname} written")
    plt.close(fig)
    
        # plot rebinneded LO QCD hist
    fig_rebinned, ax_rebinned = plt.subplots()
    fig_rebinned.subplots_adjust(hspace=.07)
    Plotter.plot_2d_norm(values=np.ma.masked_where(masked_rebinned_lo_vals.T <= 0.0, masked_rebinned_lo_vals.T), # mask nonzero probabilities for plotting
        xbins=masked_rebinned_mtt_binning, ybins=rebinned_ctstar_binning,
        xlimits=rebinned_mtt_lims, ylimits=rebinned_ctstar_lims, xlabel=mtt_title, ylabel=ctstar_title,
        ax=ax_rebinned, **{"cmap_label" : lo_ztitle})
    ax_rebinned.text(
        0.98, 0.90, "$t\\bar{t}$\nparton level",
        fontsize=rcParams["font.size"], horizontalalignment="right", verticalalignment="bottom", transform=ax_rebinned.transAxes, color="w",
    )
    hep.label.exp_label(ax=ax_rebinned, exp="MADGRAPH", rlabel="")
    figname_rebinned = os.path.join(plots_dir, base_dir, f"LO_QCD_xsec_mtt_vs_ctstar_rebinned")
    fig_rebinned.savefig(figname_rebinned)
    print(f"{figname_rebinned} written")
    plt.close(fig_rebinned)
    
        # plot rebinneded differential LO QCD hist (divide by bin width)
    fig_rebinned_diff, ax_rebinned_diff = plt.subplots()
    fig_rebinned_diff.subplots_adjust(hspace=.07)
            # make values differential
    diff_norm_vals = np.where(masked_rebinned_lo_vals.T < 1e-10, .0, masked_rebinned_lo_vals.T)
    diff_norm_vals = diff_norm_vals/(tiled_ctstar_binwidths*tiled_mtt_binwidths)
    Plotter.plot_2d_norm(values=np.ma.masked_where(diff_norm_vals <= 0.0, diff_norm_vals), # mask nonzero probabilities for plotting
        xbins=masked_rebinned_mtt_binning, ybins=rebinned_ctstar_binning,
        xlimits=rebinned_mtt_lims, ylimits=rebinned_ctstar_lims, xlabel=mtt_title, ylabel=ctstar_title,
        ax=ax_rebinned_diff, **{"cmap_label" : lo_ztitle})
    ax_rebinned_diff.text(
        0.98, 0.90, "$t\\bar{t}$\nparton level",
        fontsize=rcParams["font.size"], horizontalalignment="right", verticalalignment="bottom", transform=ax_rebinned_diff.transAxes, color="w",
    )
    hep.label.exp_label(ax=ax_rebinned_diff, exp="MADGRAPH", rlabel="")
    figname_rebinned_diff = os.path.join(plots_dir, base_dir, f"LO_QCD_xsec_mtt_vs_ctstar_rebinned_Differential")
    fig_rebinned_diff.savefig(figname_rebinned_diff)
    print(f"{figname_rebinned_diff} written")
    plt.close(fig_rebinned)
    
    #set_trace()
    for yt_name, yt_val in yt_vals_dict.items():
            # original NLO EW dists
        if (yt_name == "111") or (yt_name == "088"): # extract value throught interpolation of yt=1 and yt=1.1
            #set_trace()
                # original
            _, (mtt_binning, ctstar_binning) = ewk_dists[("EWYt0/tot", "dense_lookup")]
            yt0_nlo = dense_lookup(*ewk_dists[("EWYt0/tot", "dense_lookup")])
            yt10_nlo = dense_lookup(*ewk_dists[("EWYt10/tot", "dense_lookup")])
            yt11_nlo = dense_lookup(*ewk_dists[("EWYt11/tot", "dense_lookup")])
            nlo_xsec_vals = np.zeros(yt0_nlo._values.shape)
            fit_xcenters = np.array( [(yt0_nlo._axes[0][i+1]+yt0_nlo._axes[0][i])/2 for i in range( yt0_nlo._axes[0].size - 1)] )
            fit_ycenters = np.array( [(yt0_nlo._axes[1][i+1]+yt0_nlo._axes[1][i])/2 for i in range( yt0_nlo._axes[1].size - 1)] )
            for xbin, xval in enumerate(fit_xcenters):
                for ybin, yval in enumerate(fit_ycenters):
                    nlo_xsec_vals[xbin, ybin] = np.poly1d(np.polyfit(np.array([0., 1.0, 1.1]), np.array([yt0_nlo(xval, yval), yt10_nlo(xval, yval), yt11_nlo(xval, yval)]), 2))(yt_val)

                # rebinned
            _, (rebinned_ctstar_binning, rebinned_mtt_binning) = ewk_dists[("EWYt0/tot_rebinned", "dense_lookup")]
            rebinned_yt0_nlo = dense_lookup(*ewk_dists[("EWYt0/tot_rebinned", "dense_lookup")])
            rebinned_yt10_nlo = dense_lookup(*ewk_dists[("EWYt10/tot_rebinned", "dense_lookup")])
            rebinned_yt11_nlo = dense_lookup(*ewk_dists[("EWYt11/tot_rebinned", "dense_lookup")])
            rebinned_nlo_xsec_vals = np.zeros(rebinned_yt0_nlo._values.shape)
            rebinned_fit_xcenters = np.array( [(rebinned_yt0_nlo._axes[0][i+1]+rebinned_yt0_nlo._axes[0][i])/2 for i in range( rebinned_yt0_nlo._axes[0].size - 1)] )
            rebinned_fit_ycenters = np.array( [(rebinned_yt0_nlo._axes[1][i+1]+rebinned_yt0_nlo._axes[1][i])/2 for i in range( rebinned_yt0_nlo._axes[1].size - 1)] )
            for xbin, xval in enumerate(rebinned_fit_xcenters):
                for ybin, yval in enumerate(rebinned_fit_ycenters):
                    rebinned_nlo_xsec_vals[xbin, ybin] = np.poly1d(np.polyfit(np.array([0., 1.0, 1.1]), np.array([rebinned_yt0_nlo(xval, yval), rebinned_yt10_nlo(xval, yval), rebinned_yt11_nlo(xval, yval)]), 2))(yt_val)

        else:
            nlo_xsec_vals, (mtt_binning, ctstar_binning) = ewk_dists[(f"EWYt{yt_name}/tot", "dense_lookup")]
            rebinned_nlo_xsec_vals, (rebinned_ctstar_binning, rebinned_mtt_binning) = ewk_dists[(f"EWYt{yt_name}/tot_rebinned", "dense_lookup")]
        #set_trace()
        NLO_to_LO_KFactor = np.where(nlo_xsec_vals == 0., 1, nlo_xsec_vals/lo_xsec_vals)
        DeltaEW_vals = np.where(lo_xsec_vals == 0., 1, (nlo_xsec_vals - lo_xsec_vals)/lo_xsec_vals)
    
            # rebinned NLO EW dists
        rebinned_NLO_to_LO_KFactor = np.where(rebinned_nlo_xsec_vals == 0., 1, rebinned_nlo_xsec_vals/rebinned_lo_xsec_vals)
        rebinned_DeltaEW_vals = np.where(rebinned_lo_xsec_vals == 0., 1, (rebinned_nlo_xsec_vals - rebinned_lo_xsec_vals)/rebinned_lo_xsec_vals)
        
        if (args.uncs == "All") or args.save_ratios:
        #if (args.uncs == "All") and args.save_ratios:
            #set_trace()
            save_dict[f"EW_KFactor_{yt_val}"] = dense_lookup(NLO_to_LO_KFactor, (mtt_binning, ctstar_binning))
            save_dict[f"DeltaEW_{yt_val}"] = dense_lookup(DeltaEW_vals, (mtt_binning, ctstar_binning))
            save_dict[f"Rebinned_EW_KFactor_{yt_val}"] = dense_lookup(rebinned_NLO_to_LO_KFactor.T, (rebinned_mtt_binning, rebinned_ctstar_binning))
            save_dict[f"Rebinned_DeltaEW_{yt_val}"] = dense_lookup(rebinned_DeltaEW_vals.T, (rebinned_mtt_binning, rebinned_ctstar_binning))
    
        ############### Plot hists ################
                # histo plotting params
        if yt_val == 1.0:
            nlo_ztitle = yt1_nlo_ewk_xsec_base_label
            kfactor_ztitle = yt1_kfactor_base_label
            deltaEW_ztitle = yt1_deltaEW_base_label
        else:
            nlo_ztitle = nlo_ewk_xsec_base_label.replace("YT", f"{yt_val}")
            kfactor_ztitle = kfactor_base_label.replace("YT", f"{yt_val}")
            deltaEW_ztitle = deltaEW_base_label.replace("YT", f"{yt_val}")

                # only get mtt values < 2000 GeV
        masked_DeltaEW_vals = DeltaEW_vals[np.where(mtt_binning <= max_mtt_val)[0][0]:np.where(mtt_binning <= max_mtt_val)[0][-1], :]
        masked_mtt_binning = mtt_binning[np.where(mtt_binning <= max_mtt_val)[0]]
        masked_rebinned_DeltaEW_vals = rebinned_DeltaEW_vals[:, np.where(rebinned_mtt_binning <= max_mtt_val)[0][0]:np.where(rebinned_mtt_binning <= max_mtt_val)[0][-1]]
        masked_rebinned_mtt_binning = rebinned_mtt_binning[np.where(rebinned_mtt_binning <= max_mtt_val)[0]]
        masked_nlo_vals = nlo_xsec_vals[np.where(mtt_binning <= max_mtt_val)[0][0]:np.where(mtt_binning <= max_mtt_val)[0][-1], :]
        masked_rebinned_nlo_vals = rebinned_nlo_xsec_vals[:, np.where(rebinned_mtt_binning <= max_mtt_val)[0][0]:np.where(rebinned_mtt_binning <= max_mtt_val)[0][-1]]
        masked_NLO_to_LO_KFactor = NLO_to_LO_KFactor[np.where(mtt_binning <= max_mtt_val)[0][0]:np.where(mtt_binning <= max_mtt_val)[0][-1], :]
        masked_rebinned_NLO_to_LO_KFactor = rebinned_NLO_to_LO_KFactor[:, np.where(rebinned_mtt_binning <= max_mtt_val)[0][0]:np.where(rebinned_mtt_binning <= max_mtt_val)[0][-1]]
        #set_trace()    

            ### plot original DeltaEW hist
        fig, ax = plt.subplots()
        fig.subplots_adjust(hspace=.07)
        Plotter.plot_2d_norm(values=np.ma.masked_where(masked_DeltaEW_vals == 1.0, masked_DeltaEW_vals),
            xbins=masked_mtt_binning, ybins=ctstar_binning,
            xlimits=tot_mtt_lims, ylimits=tot_ctstar_lims, xlabel=mtt_title, ylabel=ctstar_title,
            ax=ax, **{"cmap_label" : deltaEW_ztitle})
        ax.text(
            0.98, 0.90, "$t\\bar{t}$\nparton level",
            fontsize=rcParams["font.size"], horizontalalignment="right", verticalalignment="bottom", transform=ax.transAxes, color="w",
        )
        hep.label.exp_label(ax=ax, rlabel="")
        figname = os.path.join(outdir, f"DeltaEW_xsec_{yt_name}_mtt_vs_ctstar")
        fig.savefig(figname)
        print(f"{figname} written")
        plt.close(fig)
    
            ### plot rebinned DeltaEW hist
        fig_rebinned, ax_rebinned = plt.subplots()
        fig_rebinned.subplots_adjust(hspace=.07)
        Plotter.plot_2d_norm(values=np.ma.masked_where(masked_rebinned_DeltaEW_vals.T == 1.0, masked_rebinned_DeltaEW_vals.T),
            xbins=masked_rebinned_mtt_binning, ybins=rebinned_ctstar_binning,
            xlimits=rebinned_mtt_lims, ylimits=rebinned_ctstar_lims, xlabel=mtt_title, ylabel=ctstar_title,
            ax=ax_rebinned, **{"cmap_label" : deltaEW_ztitle})
        ax_rebinned.text(
            0.98, 0.90, "$t\\bar{t}$\nparton level",
            fontsize=rcParams["font.size"], horizontalalignment="right", verticalalignment="bottom", transform=ax_rebinned.transAxes, color="w",
        )
        hep.label.exp_label(ax=ax_rebinned, rlabel="")
        figname_rebinned = os.path.join(outdir, f"DeltaEW_xsec_{yt_name}_mtt_vs_ctstar_rebinned")
        fig_rebinned.savefig(figname_rebinned)
        print(f"{figname_rebinned} written")
        plt.close(fig_rebinned)
    
    
            ### plot original nlo ewk hist
        fig, ax = plt.subplots()
        fig.subplots_adjust(hspace=.07)
        Plotter.plot_2d_norm(values=np.ma.masked_where(masked_nlo_vals <= 0.0, masked_nlo_vals), # mask nonzero probabilities for plotting
            xbins=masked_mtt_binning, ybins=ctstar_binning,
            xlimits=tot_mtt_lims, ylimits=tot_ctstar_lims, xlabel=mtt_title, ylabel=ctstar_title,
            ax=ax, **{"cmap_label" : nlo_ztitle})
        ax.text(
            0.98, 0.90, "$t\\bar{t}$\nparton level",
            fontsize=rcParams["font.size"], horizontalalignment="right", verticalalignment="bottom", transform=ax.transAxes, color="w",
        )
        hep.label.exp_label(ax=ax, exp="HATHOR", rlabel="")
        figname = os.path.join(outdir, f"NLO_EW_xsec_{yt_name}_mtt_vs_ctstar")
        fig.savefig(figname)
        print(f"{figname} written")
        plt.close(fig)

    
            ### plot rebinned nlo ewk hist
        fig_rebinned, ax_rebinned = plt.subplots()
        fig_rebinned.subplots_adjust(hspace=.07)
        Plotter.plot_2d_norm(values=np.ma.masked_where(masked_rebinned_nlo_vals.T <= 0.0, masked_rebinned_nlo_vals.T), # mask nonzero probabilities for plotting
            xbins=masked_rebinned_mtt_binning, ybins=rebinned_ctstar_binning,
            xlimits=rebinned_mtt_lims, ylimits=rebinned_ctstar_lims, xlabel=mtt_title, ylabel=ctstar_title,
            ax=ax_rebinned, **{"cmap_label" : nlo_ztitle})
        ax_rebinned.text(
            0.98, 0.90, "$t\\bar{t}$\nparton level",
            fontsize=rcParams["font.size"], horizontalalignment="right", verticalalignment="bottom", transform=ax_rebinned.transAxes, color="w",
        )
        hep.label.exp_label(ax=ax_rebinned, exp="HATHOR", rlabel="")
        figname_rebinned = os.path.join(outdir, f"NLO_EW_xsec_{yt_name}_mtt_vs_ctstar_rebinned")
        fig_rebinned.savefig(figname_rebinned)
        print(f"{figname_rebinned} written")
        plt.close(fig_rebinned)
   
        #set_trace() 
            ### plot rebinned differential nlo ewk hist (divide by bin width)
        fig_rebinned_diff, ax_rebinned_diff = plt.subplots()
        fig_rebinned_diff.subplots_adjust(hspace=.07)
                # make values differential
        diff_norm_vals = np.where(masked_rebinned_nlo_vals.T < 1e-10, .0, masked_rebinned_nlo_vals.T)
        diff_norm_vals = diff_norm_vals/(tiled_ctstar_binwidths*tiled_mtt_binwidths)
        Plotter.plot_2d_norm(values=np.ma.masked_where(diff_norm_vals <= 0.0, diff_norm_vals), # mask nonzero probabilities for plotting
            xbins=masked_rebinned_mtt_binning, ybins=rebinned_ctstar_binning,
            xlimits=rebinned_mtt_lims, ylimits=rebinned_ctstar_lims, xlabel=mtt_title, ylabel=ctstar_title,
            ax=ax_rebinned_diff, **{"cmap_label" : nlo_ztitle})
        ax_rebinned_diff.text(
            0.98, 0.90, "$t\\bar{t}$\nparton level",
            fontsize=rcParams["font.size"], horizontalalignment="right", verticalalignment="bottom", transform=ax_rebinned_diff.transAxes, color="w",
        )
        hep.label.exp_label(ax=ax_rebinned_diff, exp="HATHOR", rlabel="")
        figname_rebinned_diff = os.path.join(outdir, f"NLO_EW_xsec_{yt_name}_mtt_vs_ctstar_rebinned_Differential")
        fig_rebinned_diff.savefig(figname_rebinned_diff)
        print(f"{figname_rebinned_diff} written")
        plt.close(fig_rebinned_diff)
    
    
    
            ### plot original kfactor hist
        fig, ax = plt.subplots()
        fig.subplots_adjust(hspace=.07)
        Plotter.plot_2d_norm(values=np.ma.masked_where(masked_NLO_to_LO_KFactor == 1.0, masked_NLO_to_LO_KFactor),
            xbins=masked_mtt_binning, ybins=ctstar_binning,
            xlimits=tot_mtt_lims, ylimits=tot_ctstar_lims, xlabel=mtt_title, ylabel=ctstar_title,
            ax=ax, **{"cmap_label" : kfactor_ztitle})
        ax.text(
            0.98, 0.90, "$t\\bar{t}$\nparton level",
            fontsize=rcParams["font.size"], horizontalalignment="right", verticalalignment="bottom", transform=ax.transAxes, color="w",
        )
        hep.label.exp_label(ax=ax, exp="", rlabel="")
        figname = os.path.join(outdir, f"EW_KFactor_{yt_name}_mtt_vs_ctstar")
        fig.savefig(figname)
        print(f"{figname} written")
        plt.close(fig)
    

            # plot rebinned kfactor hist
        fig_rebinned, ax_rebinned = plt.subplots()
        fig_rebinned.subplots_adjust(hspace=.07)
        Plotter.plot_2d_norm(values=np.ma.masked_where(masked_rebinned_NLO_to_LO_KFactor.T == 1.0, masked_rebinned_NLO_to_LO_KFactor.T),
            xbins=masked_rebinned_mtt_binning, ybins=rebinned_ctstar_binning,
            xlimits=rebinned_mtt_lims, ylimits=rebinned_ctstar_lims, xlabel=mtt_title, ylabel=ctstar_title,
            ax=ax_rebinned, **{"cmap_label" : kfactor_ztitle})
        ax_rebinned.text(
            0.98, 0.90, "$t\\bar{t}$\nparton level",
            fontsize=rcParams["font.size"], horizontalalignment="right", verticalalignment="bottom", transform=ax_rebinned.transAxes, color="w",
        )
        hep.label.exp_label(ax=ax_rebinned, exp="", rlabel="")
        figname_rebinned = os.path.join(outdir, f"EW_KFactor_{yt_name}_mtt_vs_ctstar_rebinned")
        fig_rebinned.savefig(figname_rebinned)
        print(f"{figname_rebinned} written")
        plt.close(fig_rebinned)



if (args.uncs == "NNLO") or (args.uncs == "All"):
    outdir = os.path.join(plots_dir, base_dir, "NNLO")
    #outdir = os.path.join(proj_dir, "plots", base_dir, "NNLO")
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
       
    nnlo_ztitle = "$\\sigma_{NNLO\ QCD}$"
    deltaQCD_ztitle = "$\\delta_{QCD} \\equiv \dfrac{\\sigma_{NNLO\ QCD} - \\sigma_{LO\ QCD}}{\\sigma_{NNLO\ QCD}}$"
    #set_trace()
        ## get values from NNLO root file
    nnlo_fname = "MATRIX_ttmVStheta.root" # "xsec_central" dist has only statistical uncs
    nnlo_dists = convert_histo_root_file(os.path.join(proj_dir, "NNLO_files", nnlo_fname))
    nnlo_xsec_vals, (nnlo_ctstar_binning, nnlo_mtt_binning) = nnlo_dists[("xsec_central", "dense_lookup")]

        # get DeltaQCD distribution
    DeltaQCD_vals = np.where(rebinned_lo_xsec_vals == 0., 1., (nnlo_xsec_vals - rebinned_lo_xsec_vals)/nnlo_xsec_vals)

    if (args.uncs == "All") or args.save_ratios:
    #if (args.uncs == "All") and args.save_ratios:
        #set_trace()
        save_dict[f"DeltaQCD"] = dense_lookup(DeltaQCD_vals.T, (nnlo_mtt_binning, nnlo_ctstar_binning))

    #set_trace()
            # only get mtt values < 2000 GeV
    DeltaQCD_vals = DeltaQCD_vals[:, np.where(nnlo_mtt_binning <= max_mtt_val)[0][0]:np.where(nnlo_mtt_binning <= max_mtt_val)[0][-1]]
    nnlo_xsec_vals = nnlo_xsec_vals[:, np.where(nnlo_mtt_binning <= max_mtt_val)[0][0]:np.where(nnlo_mtt_binning <= max_mtt_val)[0][-1]]
    nnlo_mtt_binning = nnlo_mtt_binning[np.where(nnlo_mtt_binning <= max_mtt_val)[0]]

        ### plot NNLO QCD hist
    fig, ax = plt.subplots()
    fig.subplots_adjust(hspace=.07)
    Plotter.plot_2d_norm(values=np.ma.masked_where(nnlo_xsec_vals.T == 0.0, nnlo_xsec_vals.T),
        xbins=nnlo_mtt_binning, ybins=nnlo_ctstar_binning,
        xlimits=rebinned_mtt_lims, ylimits=rebinned_ctstar_lims, xlabel=mtt_title, ylabel=ctstar_title,
        ax=ax, **{"cmap_label" : nnlo_ztitle})
    hep.label.exp_label(ax=ax, exp="MATRIX")
    figname = os.path.join(outdir, f"MATRIX_NNLO_QCD_TTbar_xsec_central_mtt_vs_ctstar")
    fig.savefig(figname)
    print(f"{figname} written")
    plt.close(fig)


        ### plot DeltaQCD hist
    fig, ax = plt.subplots()
    fig.subplots_adjust(hspace=.07)
    Plotter.plot_2d_norm(values=np.ma.masked_where(DeltaQCD_vals.T == 1.0, DeltaQCD_vals.T),
        xbins=nnlo_mtt_binning, ybins=nnlo_ctstar_binning,
        xlimits=rebinned_mtt_lims, ylimits=rebinned_ctstar_lims, xlabel=mtt_title, ylabel=ctstar_title,
        ax=ax, **{"cmap_label" : deltaQCD_ztitle})
    hep.label.exp_label(ax=ax, exp="")
    figname = os.path.join(outdir, f"DeltaQCDEW_xsec_mtt_vs_ctstar")
    fig.savefig(figname)
    print(f"{figname} written")
    plt.close(fig)
    

#set_trace()
if (args.uncs == "All") or args.save_ratios:
    #set_trace()
    from coffea.util import save
    ratios_fname = os.path.join(proj_dir, "EWK_Uncs", f"EWK_Corrections.coffea")
    save(save_dict, ratios_fname)
    print(f"\n{ratios_fname} written")
