#!/usr/bin/env python

# matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend("agg")
from matplotlib import rcParams
rcParams["font.size"] = 20
rcParams["savefig.format"] = "png"
rcParams["savefig.bbox"] = "tight"

import Utilities.Plotter as Plotter
#from coffea.hist import plot
from coffea import hist
from pdb import set_trace
import os
import numpy as np
from coffea.lookup_tools.root_converters import convert_histo_root_file
from coffea.lookup_tools.dense_lookup import dense_lookup
from itertools import product

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("uncs", choices=["EW", "NNLO", "All"], help="Which EW uncertainty to calculate.")
parser.add_argument("--save_ratios", action="store_true", help="Save NNLO/tune ratios for all years (default is to not save them).")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
base_jobid = os.environ["base_jobid"]
base_dir = "EWK_Uncs"


save_dict = {}
    
fname = "HATHORscaled.root"
ewk_dists = convert_histo_root_file(os.path.join(proj_dir, "EWK_TT_HATHOR", fname))

# LO QCD dist is same for all directories
lo_xsec_vals, (tot_mtt_binning, tot_ctstar_binning) = ewk_dists[("EWno/tot", "dense_lookup")]
rebinned_lo_xsec_vals, (rebinned_ctstar_binning, rebinned_mtt_binning) = ewk_dists[("EWno/tot_rebinned", "dense_lookup")]
lo_ztitle = "$\\sigma_{LO\ QCD}$"

        # histo plotting params
mtt_title, ctstar_title = "$m_{t\\bar{t}}$", "cos($\\theta^{*}$)"
tot_mtt_lims = (np.min(tot_mtt_binning), np.max(tot_mtt_binning))
tot_ctstar_lims = (np.min(tot_ctstar_binning), np.max(tot_ctstar_binning))
rebinned_mtt_lims = (np.min(rebinned_mtt_binning), np.max(rebinned_mtt_binning))
rebinned_ctstar_lims = (np.min(rebinned_ctstar_binning), np.max(rebinned_ctstar_binning))
    
if (args.uncs == "EW") or (args.uncs == "All"):
    outdir = os.path.join(proj_dir, "plots", base_dir, "EW")
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
        "40" : 4.0
    }
        
    #nameTOval = lambda val : str(val).replace('p', '.')
    nlo_ewk_xsec_base_label = "$\\sigma_{LO\ QCD\ +\ NLO\ EW}$, $y^{t}/y^{t}_{SM}$ = YT"
    kfactor_base_label = "$K^{YT}_{NLO\ EW} \\equiv \dfrac{\\sigma_{LO\ QCD\ +\ NLO\ EW}}{\\sigma_{LO\ QCD}}$, $y^{t}/y^{t}_{SM}$ = YT"
    deltaEW_base_label = "$\\delta^{YT}_{EW} \\equiv \dfrac{\\sigma_{LO\ QCD\ +\ NLO\ EW} - \\sigma_{LO\ QCD}}{\\sigma_{LO\ QCD}}$, $y^{t}/y^{t}_{SM}$ = YT"
    
    
        ### plot original LO QCD hist
    fig, ax = plt.subplots()
    fig.subplots_adjust(hspace=.07)
    Plotter.plot_2d_norm(values=np.ma.masked_where(lo_xsec_vals <= 0.0, lo_xsec_vals), # mask nonzero probabilities for plotting
        xbins=tot_mtt_binning, ybins=tot_ctstar_binning,
        xlimits=tot_mtt_lims, ylimits=tot_ctstar_lims, xlabel=mtt_title, ylabel=ctstar_title,
        ax=ax, **{"cmap_label" : lo_ztitle})
    hep.label.exp_label(ax=ax, exp="MADGRAPH")
    figname = os.path.join(proj_dir, "plots", base_dir, f"LO_QCD_xsec_mtt_vs_ctstar")
    fig.savefig(figname)
    print(f"{figname} written")
    plt.close()    
    
        # plot rebinneded LO QCD hist
    fig_rebinned, ax_rebinned = plt.subplots()
    fig_rebinned.subplots_adjust(hspace=.07)
    Plotter.plot_2d_norm(values=np.ma.masked_where(rebinned_lo_xsec_vals.T <= 0.0, rebinned_lo_xsec_vals.T), # mask nonzero probabilities for plotting
        xbins=rebinned_mtt_binning, ybins=rebinned_ctstar_binning,
        xlimits=rebinned_mtt_lims, ylimits=rebinned_ctstar_lims, xlabel=mtt_title, ylabel=ctstar_title,
        ax=ax_rebinned, **{"cmap_label" : lo_ztitle})
    hep.label.exp_label(ax=ax_rebinned, exp="MADGRAPH")
    figname_rebinned = os.path.join(proj_dir, "plots", base_dir, f"LO_QCD_xsec_mtt_vs_ctstar_rebinned")
    fig_rebinned.savefig(figname_rebinned)
    print(f"{figname_rebinned} written")
    plt.close()    
    
    #set_trace()
    for yt_name, yt_val in yt_vals_dict.items():
            # original NLO EW dists
        nlo_xsec_vals, (mtt_binning, ctstar_binning) = ewk_dists[(f"EWYt{yt_name}/tot", "dense_lookup")]
        NLO_to_LO_KFactor = np.where(nlo_xsec_vals == 0., 1, nlo_xsec_vals/lo_xsec_vals)
        DeltaEW_vals = np.where(lo_xsec_vals == 0., 1, (nlo_xsec_vals - lo_xsec_vals)/lo_xsec_vals)
    
            # rebinned NLO EW dists
        rebinned_nlo_xsec_vals, (rebinned_ctstar_binning, rebinned_mtt_binning) = ewk_dists[(f"EWYt{yt_name}/tot_rebinned", "dense_lookup")]
        rebinned_NLO_to_LO_KFactor = np.where(rebinned_nlo_xsec_vals == 0., 1, rebinned_nlo_xsec_vals/rebinned_lo_xsec_vals)
        rebinned_DeltaEW_vals = np.where(rebinned_lo_xsec_vals == 0., 1, (rebinned_nlo_xsec_vals - rebinned_lo_xsec_vals)/rebinned_lo_xsec_vals)
        
        if (args.uncs == "All") and args.save_ratios:
            #set_trace()
            save_dict[f"EW_KFactor_{yt_val}"] = dense_lookup(NLO_to_LO_KFactor, (mtt_binning, ctstar_binning))
            save_dict[f"DeltaEW_{yt_val}"] = dense_lookup(DeltaEW_vals, (mtt_binning, ctstar_binning))
            save_dict[f"Rebinned_EW_KFactor_{yt_val}"] = dense_lookup(rebinned_NLO_to_LO_KFactor.T, (rebinned_mtt_binning, rebinned_ctstar_binning))
            save_dict[f"Rebinned_DeltaEW_{yt_val}"] = dense_lookup(rebinned_DeltaEW_vals.T, (rebinned_mtt_binning, rebinned_ctstar_binning))
    
        ############### Plot hists ################
                # histo plotting params
        nlo_ztitle = nlo_ewk_xsec_base_label.replace("YT", f"{yt_val}")
        kfactor_ztitle = kfactor_base_label.replace("YT", f"{yt_val}")
        deltaEW_ztitle = deltaEW_base_label.replace("YT", f"{yt_val}")
    
            ### plot original DeltaEW hist
        fig, ax = plt.subplots()
        fig.subplots_adjust(hspace=.07)
        #Plotter.plot_2d_norm(values=DeltaEW_vals,
        Plotter.plot_2d_norm(values=np.ma.masked_where(DeltaEW_vals == 1.0, DeltaEW_vals),
            xbins=mtt_binning, ybins=ctstar_binning,
            xlimits=tot_mtt_lims, ylimits=tot_ctstar_lims, xlabel=mtt_title, ylabel=ctstar_title,
            ax=ax, **{"cmap_label" : deltaEW_ztitle})
        hep.label.exp_label(ax=ax, exp="HATHOR")
        figname = os.path.join(outdir, f"DeltaEW_xsec_{yt_name}_mtt_vs_ctstar")
        fig.savefig(figname)
        print(f"{figname} written")
        plt.close()    
    
            ### plot rebinneded DeltaEW hist
        fig_rebinned, ax_rebinned = plt.subplots()
        fig_rebinned.subplots_adjust(hspace=.07)
        Plotter.plot_2d_norm(values=np.ma.masked_where(rebinned_DeltaEW_vals.T == 1.0, rebinned_DeltaEW_vals.T),
            xbins=rebinned_mtt_binning, ybins=rebinned_ctstar_binning,
            xlimits=rebinned_mtt_lims, ylimits=rebinned_ctstar_lims, xlabel=mtt_title, ylabel=ctstar_title,
            ax=ax_rebinned, **{"cmap_label" : deltaEW_ztitle})
        hep.label.exp_label(ax=ax_rebinned, exp="HATHOR")
        figname_rebinned = os.path.join(outdir, f"DeltaEW_xsec_{yt_name}_mtt_vs_ctstar_rebinned")
        fig_rebinned.savefig(figname_rebinned)
        print(f"{figname_rebinned} written")
        plt.close()
    
    
            ### plot original nlo ewk hist
        fig, ax = plt.subplots()
        fig.subplots_adjust(hspace=.07)
        Plotter.plot_2d_norm(values=np.ma.masked_where(nlo_xsec_vals <= 0.0, nlo_xsec_vals), # mask nonzero probabilities for plotting
            xbins=mtt_binning, ybins=ctstar_binning,
            xlimits=tot_mtt_lims, ylimits=tot_ctstar_lims, xlabel=mtt_title, ylabel=ctstar_title,
            ax=ax, **{"cmap_label" : nlo_ztitle})
        hep.label.exp_label(ax=ax, exp="HATHOR")
        figname = os.path.join(outdir, f"NLO_EW_xsec_{yt_name}_mtt_vs_ctstar")
        fig.savefig(figname)
        print(f"{figname} written")
        plt.close()    

    
            ### plot rebinned nlo ewk hist
        fig_rebinned, ax_rebinned = plt.subplots()
        fig_rebinned.subplots_adjust(hspace=.07)
        Plotter.plot_2d_norm(values=np.ma.masked_where(rebinned_nlo_xsec_vals.T <= 0.0, rebinned_nlo_xsec_vals.T), # mask nonzero probabilities for plotting
            xbins=rebinned_mtt_binning, ybins=rebinned_ctstar_binning,
            xlimits=rebinned_mtt_lims, ylimits=rebinned_ctstar_lims, xlabel=mtt_title, ylabel=ctstar_title,
            ax=ax_rebinned, **{"cmap_label" : nlo_ztitle})
        hep.label.exp_label(ax=ax_rebinned, exp="HATHOR")
        figname_rebinned = os.path.join(outdir, f"NLO_EW_xsec_{yt_name}_mtt_vs_ctstar_rebinned")
        fig_rebinned.savefig(figname_rebinned)
        print(f"{figname_rebinned} written")
        plt.close()
    
    
            ### plot original kfactor hist
        fig, ax = plt.subplots()
        fig.subplots_adjust(hspace=.07)
        #Plotter.plot_2d_norm(values=NLO_to_LO_KFactor,
        Plotter.plot_2d_norm(values=np.ma.masked_where(NLO_to_LO_KFactor == 1.0, NLO_to_LO_KFactor),
            xbins=mtt_binning, ybins=ctstar_binning,
            xlimits=tot_mtt_lims, ylimits=tot_ctstar_lims, xlabel=mtt_title, ylabel=ctstar_title,
            ax=ax, **{"cmap_label" : kfactor_ztitle})
        hep.label.exp_label(ax=ax, exp="HATHOR")
        figname = os.path.join(outdir, f"EW_KFactor_{yt_name}_mtt_vs_ctstar")
        fig.savefig(figname)
        print(f"{figname} written")
        plt.close()    
    

            # plot rebinned kfactor hist
        fig_rebinned, ax_rebinned = plt.subplots()
        fig_rebinned.subplots_adjust(hspace=.07)
        Plotter.plot_2d_norm(values=np.ma.masked_where(rebinned_NLO_to_LO_KFactor.T == 1.0, rebinned_NLO_to_LO_KFactor.T),
            xbins=rebinned_mtt_binning, ybins=rebinned_ctstar_binning,
            xlimits=rebinned_mtt_lims, ylimits=rebinned_ctstar_lims, xlabel=mtt_title, ylabel=ctstar_title,
            ax=ax_rebinned, **{"cmap_label" : kfactor_ztitle})
        hep.label.exp_label(ax=ax_rebinned, exp="HATHOR")
        figname_rebinned = os.path.join(outdir, f"EW_KFactor_{yt_name}_mtt_vs_ctstar_rebinned")
        fig_rebinned.savefig(figname_rebinned)
        print(f"{figname_rebinned} written")
        plt.close()



if (args.uncs == "NNLO") or (args.uncs == "All"):
    outdir = os.path.join(proj_dir, "plots", base_dir, "NNLO")
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

    if (args.uncs == "All") and args.save_ratios:
        #set_trace()
        save_dict[f"DeltaQCD"] = dense_lookup(DeltaQCD_vals.T, (nnlo_mtt_binning, nnlo_ctstar_binning))

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
    plt.close()    


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
    plt.close()    
    



if args.save_ratios:
    #set_trace()
    from coffea.util import save
    ratios_fname = os.path.join(proj_dir, "EWK_Uncs", f"EWK_Corrections.coffea")
    save(save_dict, ratios_fname)
    print(f"\n{ratios_fname} written")
