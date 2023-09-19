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

from coffea.util import load
from pdb import set_trace
import os
import Utilities.prettyjson as prettyjson
from coffea import hist
import numpy as np
import fnmatch
import Utilities.Plotter as Plotter
import Utilities.final_analysis_binning as final_binning
import Utilities.common_features as cfeatures

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2018"], help="What year is the ntuple from.")
#parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("lepton", choices=["Electron", "Muon"], help="Choose which lepton to make plots for")
args = parser.parse_args()

base_jobid = os.environ["base_jobid"]
proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
eos_dir = os.environ["eos_dir"]
plots_dir = os.environ["plots_dir"]

## get file for SM ttbar dists
SMtt_fnames = fnmatch.filter(os.listdir(os.path.join(eos_dir, "results", f"{args.year}_{jobid}", "htt_btag_sb_regions")), "*BKG*TOT.coffea")
SMtt_fnames = [os.path.join(eos_dir, "results", f"{args.year}_{jobid}", "htt_btag_sb_regions", fname) for fname in SMtt_fnames]
if len(SMtt_fnames) > 1: raise ValueError("Multiple SM files found")
SMtt_hdict = load(SMtt_fnames[0])
print("Loaded SM ttbar file")

## get toponium file
EtaT_fnames = fnmatch.filter(os.listdir(os.path.join(eos_dir, "results", jobid, "Toponium_HBSR", args.year)), "*TOT.coffea")
EtaT_fnames = [os.path.join(eos_dir, "results", jobid, "Toponium_HBSR", args.year, fname) for fname in EtaT_fnames]
if len(EtaT_fnames) > 1: raise ValueError("Multiple toponium files found")
EtaT_hdict = load(EtaT_fnames[0])
print("Loaded toponium file")

outdir = os.path.join(plots_dir, jobid, "Compare_SMtt_EtaT", args.year, args.lepton)
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
        (linearize_binning[0][0], linearize_binning[0][-1]), (linearize_binning[1][0], linearize_binning[1][-1]), None, (ctstar_binlabels, ctstar_bin_locs)),
    "Jets_phi_vs_eta" : (cfeatures.variable_names_to_labels["Jets_phi"], cfeatures.variable_names_to_labels["Jets_eta"],
        phi_eta_binning[0], phi_eta_binning[1], (-3.3, 3.3), (-2.6, 2.6), phi_binlabels, (eta_binlabels, eta_bin_locs)),
    "Lep_phi_vs_eta" : (cfeatures.variable_names_to_labels["Lep_phi"] % cfeatures.objtypes[args.lepton], cfeatures.variable_names_to_labels["Lep_eta"] % cfeatures.objtypes[args.lepton],
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
    "Jets_DeepCSV_bDisc" : (cfeatures.variable_names_to_labels["Jets_DeepCSV_bDisc"], 1, (-0.01, 1.)),
    "Jets_DeepJet_bDisc" : (cfeatures.variable_names_to_labels["Jets_DeepJet_bDisc"], 1, (-0.01, 1.)),
    "Lep_pt" : (cfeatures.variable_names_to_labels["Lep_pt"] % cfeatures.objtypes[args.lepton], 1, (0., 300.)),
    "Lep_eta" : (cfeatures.variable_names_to_labels["Lep_eta"] % cfeatures.objtypes[args.lepton], 2, (-2.6, 2.6)),
    "Lep_phi" : (cfeatures.variable_names_to_labels["Lep_phi"] % cfeatures.objtypes[args.lepton], 2, (-4, 4)),
    "Lep_iso" : (cfeatures.variable_names_to_labels["Lep_iso"] % cfeatures.objtypes[args.lepton], 1, (0., 1.)),
    "MT" : (cfeatures.variable_names_to_labels["MT"], 1, (0., 300.)),
    "MET_pt" : (cfeatures.variable_names_to_labels["MET_pt"], 1, (0., 300.)),
    "MET_phi" : (cfeatures.variable_names_to_labels["MET_phi"], 1, (-3.2, 3.2)),
}


    ## get data lumi and scale MC by lumi
lumi_to_use = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())[args.year][f"{args.lepton}s"]
lumi_correction = load(os.path.join(proj_dir, "Corrections", base_jobid, "MC_LumiWeights.coffea"))[args.year][f"{args.lepton}s"]
        # scale ttJets events, split by reconstruction type, by normal ttJets lumi correction
ttJets_groups = ["ttJetsDiLep_other", "ttJetsHad_other", "ttJetsSL_right", "ttJetsSL_matchable", "ttJetsSL_unmatchable", "ttJetsSL_sl_tau"]
#set_trace()
for tt_cat in ttJets_groups:
    ttJets_lumi_topo = "_".join(tt_cat.split("_")[:-2]) if "sl_tau" in tt_cat else "_".join(tt_cat.split("_")[:-1]) # gets ttJets[SL, Had, DiLep] or ttJets_PS
    ttJets_eff_lumi = lumi_correction[ttJets_lumi_topo]
    lumi_correction.update({tt_cat: ttJets_eff_lumi})


styles_dict = {
    "ttJets" : {"color" : "#e41a1c", "label" : "$t\\bar{t}$"},
    "EtaT"   : {"color" : "#377eb8", "label" : "$\\eta_{t}$"},
}

scaled_hists = {}
    # scale and group hists by process
for hname in SMtt_hdict.keys():
    if "cutflow" in hname: continue
    if hname not in variables.keys(): continue
    #set_trace()
    tmp_SM_EtaT_hist = (SMtt_hdict[hname]["ttJets*", "nosys", :, args.lepton, "btagPass"].integrate("sys").integrate("leptype").integrate("btag")).copy()
    tmp_EtaT_hist = (EtaT_hdict[hname][:, "nosys", :, args.lepton, "btagPass"].integrate("sys").integrate("leptype").integrate("btag")).copy()
    tmp_SM_EtaT_hist.add(tmp_EtaT_hist)
    tmp_SM_EtaT_hist.scale(lumi_correction, axis="dataset")
    tmp_SM_EtaT_hist = tmp_SM_EtaT_hist.group("dataset", hist.Cat("process", "Process", sorting="placement"), {"ttJets" : ttJets_groups, "EtaT" : ["ToponiumSL", "ToponiumDiLep"]})
    scaled_hists[hname] = tmp_SM_EtaT_hist.copy()



    ## make plots
for hname, histo in scaled_hists.items():
    if histo.dense_dim() == 1:
        xtitle, rebinning, x_lims= variables[hname]
        orig_xtitle = xtitle
        if rebinning != 1:
            xaxis_name = histo.dense_axes()[0].name
            histo = histo.rebin(xaxis_name, rebinning)

        #for jmult in ["3Jets"]:
        for jmult in sorted(set([key[1] for key in histo.values().keys()])):
            print(*[jmult, hname], sep=", ")
            pltdir = os.path.join(outdir, jmult)
            if not os.path.isdir(pltdir):
                os.makedirs(pltdir)

            fig, ax = plt.subplots()
            fig.subplots_adjust(hspace=.07)

            #set_trace()
            hslice = (histo[:, jmult].integrate("jmult")).copy()
            for topo in sorted(set([key[0] for key in hslice.values().keys()])):
                masked_vals, masked_bins = Plotter.get_ratio_arrays(num_vals=hslice.values()[(topo,)], denom_vals=np.ones_like(hslice.values()[(topo,)]), input_bins=hslice.dense_axes()[0].edges())
                ax.step(masked_bins, masked_vals, where="post", **styles_dict[topo])

            if "disc" in hname:
                xtitle = orig_xtitle.replace("j", "3") if jmult == "3Jets" else orig_xtitle.replace("j", "4")

            if hname == "Lep_iso":
                x_lims = (0., 0.15) if (args.lepton == "Muon") else (0., 0.1)

            ax.legend(loc="upper right", ncol=1, fontsize=24)
            ax.autoscale()
            ax.set_ylim(1e-3, ax.get_ylim()[1]*1.3)
            #ax.set_ylim(0, ax.get_ylim()[1]*1.3)
            ax.set_yscale("log")
            ax.set_xlabel(xtitle)
            ax.set_xlim(x_lims)

                # add lepton/jet multiplicity label
            ax.text(
                0.02, 0.92, cfeatures.channel_labels[f"{args.lepton}_{jmult}"],
                fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
            )
            hep.cms.label(ax=ax, data=False, label="Preliminary", year=year_to_use, lumi=round(lumi_to_use/1000., 1))

            #set_trace()
            figname = os.path.join(pltdir, "_".join([jmult, args.lepton, hname, "Compare_SMtt_EtaT"]))
            fig.savefig(figname)
            print(f"{figname} written")
            plt.close(fig)
    
    if histo.dense_dim() == 2:
        #set_trace()
        xtitle, ytitle, xrebinning, yrebinning, x_lims, y_lims, plot_xlabels, plot_ylabels = variables[hname]

        xaxis_name = histo.dense_axes()[0].name
        yaxis_name = histo.dense_axes()[1].name
            ## rebin x axis
        if isinstance(xrebinning, np.ndarray):
            new_xbins = hist.Bin(xaxis_name, xaxis_name, xrebinning)
        elif isinstance(xrebinning, float) or isinstance(xrebinning, int):
            new_xbins = xrebinning
        #histo = histo.rebin(xaxis_name, new_xbins)
            ## rebin y axis
        if isinstance(yrebinning, np.ndarray):
            new_ybins = hist.Bin(yaxis_name, yaxis_name, yrebinning)
        elif isinstance(yrebinning, float) or isinstance(yrebinning, int):
            new_ybins = yrebinning
        histo = histo.rebin(yaxis_name, new_ybins).rebin(xaxis_name, new_xbins)

        for jmult in sorted(set([key[1] for key in histo.values().keys()])):
            print(*[jmult, hname], sep=", ")
            pltdir = os.path.join(outdir, jmult)
            if not os.path.isdir(pltdir):
                os.makedirs(pltdir)

            hslice = (histo[:, jmult].integrate("jmult")).copy()
            hline = Plotter.linearize_hist(hslice)
                # plot linearized view 
            fig, ax = plt.subplots(figsize=(15.0, 10.0))
            fig.subplots_adjust(hspace=.07)

            for topo in sorted(set([key[0] for key in hline.values().keys()])):
                masked_vals, masked_bins = Plotter.get_ratio_arrays(num_vals=hline.values()[(topo,)], denom_vals=np.ones_like(hline.values()[(topo,)]), input_bins=hline.dense_axes()[0].edges())
                ax.step(masked_bins, masked_vals, where="post", **styles_dict[topo])

            #set_trace()
            ax.legend(loc="upper right", ncol=1, fontsize=24)
            ax.autoscale()
            ax.set_ylim(1e-3, ax.get_ylim()[1]*1.5)
            #ax.set_ylim(0, ax.get_ylim()[1]*1.3)
            ax.set_yscale("log")
            ax.set_xlim(0, hline.dense_axes()[0].centers().size)
            ax.set_xlabel(xtitle)

                ## draw vertical lines for distinguishing different ctstar bins
            vlines = [hslice.values()[("EtaT",)].shape[0]*ybin for ybin in range(1, hslice.values()[("EtaT",)].shape[1])]
            [ax.axvline(vline, color="k", linestyle="--") for vline in vlines]

                # plot unrolled x and y labels for each bin
            ## plot x labels
            if plot_xlabels is not None:
                #set_trace()
                ax.set_xticks(np.arange(len(plot_xlabels)))
                ax.set_xticklabels(plot_xlabels)
                ax.tick_params(which="minor", bottom=False, top=False)
                plt.setp(ax.get_xticklabels(), rotation=90, ha="left", fontsize=rcParams["font.size"]*0.5)

            ## plot y labels
            if plot_ylabels is not None: # (binlabels, bin_locs)
                for idx, label in enumerate(plot_ylabels[0]):
                    ax.annotate(label, xy=(plot_ylabels[1][idx], 0), xycoords=("data", "axes fraction"),
                        xytext=(0, 50), textcoords="offset points", va="top", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)

                if hname == "mtt_vs_tlep_ctstar_abs":
                    #set_trace()
                    ax.set_xticks(mtt_bin_inds_to_plot)
                    ax.set_xticklabels(mtt_bins_to_plot)

                # add lepton/jet multiplicity label
            ax.text(
                0.02, 0.92, cfeatures.channel_labels[f"{args.lepton}_{jmult}"],
                fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
            )
            hep.cms.label(ax=ax, data=False, label="Preliminary", year=year_to_use, lumi=round(lumi_to_use/1000., 1))

            #set_trace()
            figname = os.path.join(pltdir, "_".join([jmult, args.lepton, hname, "Compare_SMtt_EtaT"]))
            fig.savefig(figname)
            print(f"{figname} written")
            plt.close(fig)

toc = time.time()
print("Total time: %.1f" % (toc - tic))
