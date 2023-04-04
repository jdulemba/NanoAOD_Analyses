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
import Utilities.final_analysis_binning as final_binning
import Utilities.common_features as cfeatures

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="Specify which year to run over")
parser.add_argument("--saveTemplates", action="store_true", help="Save nominal templates for use in fit.")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
base_jobid = os.environ["base_jobid"]
analyzer = "NNLO_reweighting_recolevel"
eos_dir = os.environ["eos_dir"]
plots_dir = os.environ["plots_dir"]

f_ext = "TOT.coffea"
input_dir = os.path.join(eos_dir, "results", jobid, analyzer, args.year)
fnames = sorted(["%s/%s" % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)])
hdict = plt_tools.add_coffea_files(fnames) if len(fnames) > 1 else load(fnames[0])

if args.saveTemplates:
    import coffea.processor as processor
    histo_dict = processor.dict_accumulator({njets : {"Muon" : {}, "Electron" :{}} for njets in ["3Jets", "4PJets"]})

outdir = os.path.join(plots_dir, jobid, analyzer, args.year)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

#set_trace()    
rewt_style_dict = {
    "nosys" : {"label" : "nominal", "color" : "k", "linestyle" : "-"},
    "Mtt_vs_top_ctstarUp" : {"label" : "nominal x $W$($m_{t\\bar{t}}$, cos($\\theta^{*}_{t}$))", "color" : "#377eb8", "linestyle" : "-"}, ## blue
    "TopPtUp" : {"label" : "nominal x $W$($p_{T}$($t$))", "color" : "#e42a2c", "linestyle" : "-"}, ## red
}

year_to_use = cfeatures.year_labels[args.year]

    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())[args.year]
lumi_correction = load(os.path.join(proj_dir, "Corrections", base_jobid, "MC_LumiWeights.coffea"))[args.year]

linearize_binning = (
    final_binning.mtt_binning,
    final_binning.ctstar_abs_binning
)

# binning and bin labels for mtt x costheta
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

for lep in sorted(["Electron", "Muon"]):
    
    variables = {
        "mtt_vs_tlep_ctstar_abs" : (cfeatures.variable_names_to_labels["mtt"], cfeatures.variable_names_to_labels["tlep_ctstar_abs"], linearize_binning[0], linearize_binning[1],
            (linearize_binning[0][0], linearize_binning[0][-1]), (linearize_binning[1][0], linearize_binning[1][-1]), None, (ctstar_binlabels, ctstar_bin_locs)),
        "Jets_phi_vs_eta" : (cfeatures.variable_names_to_labels["Jets_phi"], cfeatures.variable_names_to_labels["Jets_eta"],
            phi_eta_binning[0], phi_eta_binning[1], (-3.3, 3.3), (-2.6, 2.6), phi_binlabels, (eta_binlabels, eta_bin_locs)),
        "Lep_phi_vs_eta" : (cfeatures.variable_names_to_labels["Lep_phi"] % cfeatures.objtypes[lep], cfeatures.variable_names_to_labels["Lep_eta"] % cfeatures.objtypes[lep],
            phi_eta_binning[0], phi_eta_binning[1], (-3.3, 3.3), (-2.6, 2.6), phi_binlabels, (eta_binlabels, eta_bin_locs)),
        "Jets_njets" : (cfeatures.variable_names_to_labels["Jets_njets"], 1, (0, 15)),
        "mtt" : (cfeatures.variable_names_to_labels["mtt"], 4, (200., 2000.)),
        "tlep_ctstar_abs" : (cfeatures.variable_names_to_labels["tlep_ctstar_abs"], 1, (0., 1.)),
        "mthad" : (cfeatures.variable_names_to_labels["mthad"], 2, (0., 300.)),
        "mWHad" : (cfeatures.variable_names_to_labels["mWHad"], 2, (0., 300.)),
        "mWLep" : (cfeatures.variable_names_to_labels["mWLep"], 2, (0., 300.)),
        "pt_thad" : (cfeatures.variable_names_to_labels["pt_thad"], 2, (0., 500.)),
        "pt_tlep" : (cfeatures.variable_names_to_labels["pt_tlep"], 2, (0., 500.)),
        "pt_tt" : (cfeatures.variable_names_to_labels["pt_tt"], 2, (0., 500.)),
        "eta_thad" : (cfeatures.variable_names_to_labels["eta_thad"], 2, (-4., 4.)),
        "eta_tlep" : (cfeatures.variable_names_to_labels["eta_tlep"], 2, (-4., 4.)),
        "eta_tt" : (cfeatures.variable_names_to_labels["eta_tt"], 2, (-4., 4.)),
        "tlep_ctstar" : (cfeatures.variable_names_to_labels["tlep_ctstar"], 2, (-1., 1.)),
        "full_disc" : (cfeatures.variable_names_to_labels["full_disc"], 2, (5, 25.)),
        "mass_disc" : (cfeatures.variable_names_to_labels["mass_disc"], 2, (0, 20.)),
        "ns_disc" : (cfeatures.variable_names_to_labels["ns_disc"], 2, (3., 10.)),
        "ns_dist" : (cfeatures.variable_names_to_labels["ns_dist"], 1, (0., 150.)),
        "Jets_pt" : (cfeatures.variable_names_to_labels["Jets_pt"], 1, (0., 300.)),
        "Jets_eta" : (cfeatures.variable_names_to_labels["Jets_eta"], 2, (-2.6, 2.6)),
        "Jets_phi" : (cfeatures.variable_names_to_labels["Jets_phi"], 2, (-4., 4.)),
        "Jets_LeadJet_pt" : (cfeatures.variable_names_to_labels["Jets_LeadJet_pt"], 1, (0., 300.)),
        "Jets_LeadJet_eta" : (cfeatures.variable_names_to_labels["Jets_LeadJet_eta"], 2, (-2.6, 2.6)),
        "Jets_DeepJet_bDisc" : (cfeatures.variable_names_to_labels["Jets_DeepJet_bDisc"], 1, (-0.01, 1.)),
        "Lep_pt" : (cfeatures.variable_names_to_labels["Lep_pt"] % cfeatures.objtypes[lep], 1, (0., 300.)),
        "Lep_eta" : (cfeatures.variable_names_to_labels["Lep_eta"] % cfeatures.objtypes[lep], 2, (-2.6, 2.6)),
        "Lep_phi" : (cfeatures.variable_names_to_labels["Lep_phi"] % cfeatures.objtypes[lep], 2, (-4, 4)),
        "Lep_iso" : (cfeatures.variable_names_to_labels["Lep_iso"] % cfeatures.objtypes[lep], 1, (0., 1.)),
        "MT" : (cfeatures.variable_names_to_labels["MT"], 1, (0., 300.)),
        "MET_pt" : (cfeatures.variable_names_to_labels["MET_pt"], 1, (0., 300.)),
        "MET_phi" : (cfeatures.variable_names_to_labels["MET_phi"], 1, (-3.2, 3.2)),
    }
    
    #set_trace()
    for hname in variables.keys():
        if hname not in hdict.keys():
            raise ValueError(f"{hname} not found in file")
        #set_trace()
        histo = hdict[hname].copy()
            ## scale and integrate out axes
        histo.scale(lumi_correction[f"{lep}s"], axis="dataset") # scale hists
        histo = histo[:, :, :, lep].integrate("leptype").integrate("dataset") # integrate out ttbar procs and lepton flavor
    
        is2d = histo.dense_dim() == 2
        axes_to_sum = (histo.dense_axes()[0].name,) if not is2d else (histo.dense_axes()[0].name, histo.dense_axes()[1].name)
    
        #set_trace()
        if is2d:
            xtitle, ytitle, xrebinning, yrebinning, x_lims, y_lims, plot_xlabels, plot_ylabels = variables[hname]
    
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
    
                # vars only for hist
            x_binning, y_binning = histo.axis(axes_to_sum[0]).edges(), histo.axis(axes_to_sum[1]).edges()
            nbins = (len(x_binning)-1)*(len(y_binning)-1)
            vlines = [nbins*ybin/(len(y_binning)-1) for ybin in range(1, len(y_binning)-1)]
    
        else:
            xtitle, rebinning, x_lims = variables[hname]
            #set_trace()
            if isinstance(rebinning, np.ndarray):
                new_xbins = hist.Bin(axes_to_sum[0], axes_to_sum[0], rebinning)
            elif isinstance(rebinning, float) or isinstance(rebinning, int):
                new_xbins = rebinning
            histo = histo.rebin(*axes_to_sum, new_xbins)
    
        #set_trace()
        for jmult in sorted(set([key[1] for key in histo.values().keys()])):
            print(f"{hname}: {lep}, {jmult}")
            pltdir = os.path.join(outdir, lep, jmult)
            if not os.path.isdir(pltdir):
                os.makedirs(pltdir)
    
            hslice = histo[:, jmult].integrate("jmult")
            if is2d:
                #set_trace()
                    # easier to rename sparse axis than change linearize()
                process_axis = hist.Cat("process", "Process")
                orig_axis = hslice.sparse_axes()[0]
                tmp_histo = hslice.copy()
                tmp_histo = tmp_histo.group(hslice.sparse_axes()[0].name, process_axis, {key[0]:key[0] for key in hslice.values().keys()})
                hline = Plotter.linearize_hist(tmp_histo)
                    # revert sparse axis name to original
                hline = hline.group(hline.sparse_axes()[0].name, orig_axis, {key[0]:key[0] for key in hslice.values().keys()})
                hslice = hline
    
                    # add template to coffea dict
                if args.saveTemplates and (hname == "mtt_vs_tlep_ctstar_abs"):
                    #set_trace()
                        ## symmetrize up and down variations, 'up' variation is defined as the actual template and down is symmetrized
                        # get yield vals
                    nom_vals = hslice["Mtt_vs_top_ctstarUp"].integrate("sys").values()[()]
                    up_histo = hslice["nosys"].integrate("sys").copy()
                    dw_histo = up_histo.copy()
                    up_vals = up_histo.values()[()]

                        # find relative deviation from nominal vals
                    up_rel_vals = up_vals/nom_vals - 1.
                    dw_rel_vals = -1. * up_rel_vals
                    symmetrized_ratios_dw = dw_rel_vals + 1.
                    dw_histo.values()[()][:] = symmetrized_ratios_dw * nom_vals

                    histo_dict[jmult][lep]["TT_NLOshapeUp"] = hslice["nosys"].integrate("sys").copy()
                    histo_dict[jmult][lep]["TT_NLOshapeDown"] = dw_histo.copy()
    
                # yield and ratio plots
            fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True, figsize=(15.0, 10.0)) if is2d else plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
            fig.subplots_adjust(hspace=.07)
    
                # just ratio plots
            fig_ratio, ax_ratio = plt.subplots(figsize=(15.0, 10.0)) if is2d else plt.subplots()
            fig_ratio.subplots_adjust(hspace=.07)
    
            #set_trace()
            for rewt_key in sorted(set([key[0] for key in hslice.values().keys()])):
                if rewt_key not in rewt_style_dict.keys(): continue
                hep.plot.histplot(hslice[rewt_key].integrate("sys").values()[()], hslice.dense_axes()[0].edges(), ax=ax, histtype="step", **rewt_style_dict[rewt_key]) # nosys template
                    ## plot ratios
                if rewt_key == "nosys": continue
                        # orig
                ratio_masked_vals, ratio_bins = Plotter.get_ratio_arrays(num_vals=hslice[rewt_key].integrate("sys").values()[()],
                    denom_vals=hslice["nosys"].integrate("sys").values()[()], input_bins=hslice.dense_axes()[0].edges())
                rax.step(ratio_bins, ratio_masked_vals, where="post", **rewt_style_dict[rewt_key])
    
                        # ratio
                ax_ratio.step(ratio_bins, ratio_masked_vals, where="post", **rewt_style_dict[rewt_key])
    
            #set_trace()
            # plot yields
                # format axes
            ax.autoscale()
            ax.set_xlim((0, nbins) if is2d else x_lims)
            ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.2)
            ax.set_xlabel(None)
            ax.set_ylabel("Events")
            
            rax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
            rax.set_ylabel("Ratio to Nominal")
            rax.set_xlabel(xtitle)
            
            ax_ratio.autoscale()
            ax_ratio.set_xlim((0, nbins) if is2d else x_lims)
            ax_ratio.set_ylim(ax_ratio.get_ylim()[0], ax_ratio.get_ylim()[1])
            #ax_ratio.set_ylim(ax_ratio.get_ylim()[0], ax_ratio.get_ylim()[1]*1.05)
            ax_ratio.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
            ax_ratio.set_ylabel("Ratio to Nominal")
            ax_ratio.set_xlabel(xtitle)
            
            ## set plotting styles
                ## set legend and corresponding colors
            #set_trace()
            handles, labels = ax.get_legend_handles_labels()
            ax.legend(handles,labels, loc="upper right")
            ax.text(
                0.02, 0.84, "$t\\bar{t}$"+f"\n{cfeatures.channel_labels[f'{lep}_{jmult}']}",
                fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
            )
    
            handles, labels = ax_ratio.get_legend_handles_labels()
            ax_ratio.legend(handles,labels, loc="upper right")
            ax_ratio.text(
                0.02, 0.88, "$t\\bar{t}$"+f"\n{cfeatures.channel_labels[f'{lep}_{jmult}']}",
                fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax_ratio.transAxes
            )
            if is2d:
                #set_trace()
                # plot vertical lines
                [ax.axvline(vline, color="k", linestyle="--") for vline in vlines]
                [rax.axvline(vline, color="k", linestyle="--") for vline in vlines]
                [ax_ratio.axvline(vline, color="k", linestyle="--") for vline in vlines]
    
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
                        ax_ratio.annotate(label, xy=(plot_ylabels[1][idx], 0), xycoords=("data", "axes fraction"),
                            xytext=(0, 30), textcoords="offset points", va="top", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)
    
                    if hname == "mtt_vs_tlep_ctstar_abs":
                        #set_trace()
                        rax.set_xticks(mtt_bin_inds_to_plot)
                        rax.set_xticklabels(mtt_bins_to_plot)
                        ax_ratio.set_xticks(mtt_bin_inds_to_plot)
                        ax_ratio.set_xticklabels(mtt_bins_to_plot)
    
            hep.cms.label(ax=ax, data=False, label="Preliminary", year=year_to_use, lumi=round(data_lumi_year[f"{lep}s"]/1000., 1))
    
            #set_trace()
            figname = os.path.join(pltdir, "_".join([jmult, lep, "TTbar", "NNLORewtComp", hname]))
            fig.savefig(figname)
            plt.close(fig)
            #set_trace()
            print(f"{figname} written")
    
            hep.cms.label(ax=ax_ratio, data=False, label="Preliminary", year=year_to_use, lumi=round(data_lumi_year[f"{lep}s"]/1000., 1))
            ratio_figname = os.path.join(pltdir, "_".join([jmult, lep, "TTbar", "NNLORewtComp", hname, "ratio"]))
            fig_ratio.savefig(ratio_figname)
            plt.close(fig_ratio)
            #set_trace()
            print(f"{ratio_figname} written")

    # save template to coffea dict
if args.saveTemplates:
    from coffea.util import save
    #set_trace()
    coffea_outname = os.path.join(eos_dir, "results", f"{args.year}_{jobid}", "Templates_htt_btag_sb_regions", f"raw_templates_TT_NLOshape_{args.year}_{jobid}.coffea")
    save(histo_dict, coffea_outname)
    print(f"\n\n{coffea_outname} written")

toc = time.time()
print("Total time: %.1f" % (toc - tic))
