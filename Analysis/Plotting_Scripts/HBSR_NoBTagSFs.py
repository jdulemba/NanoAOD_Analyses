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
import Utilities.btag_sideband_regions as btag_sidebands
import Utilities.common_features as cfeatures

base_jobid = os.environ["base_jobid"]

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("lepton", choices=["Electron", "Muon"], help="Choose which lepton to make plots for")
parser.add_argument("--nosys", action="store_true", help="Make plots without systematics and no qcd estimation")
parser.add_argument("--bkg_shapes", action="store_true", help="Make plots comparing EKW+QCD MC and data-driven shapes.")
parser.add_argument("--bkg_shapes_comp", action="store_true", help="Make plots comparing EKW+QCD MC and data-driven shapes with and without b-tagging SFs applied.")
#parser.add_argument("--bkg_est", action="store_true", help="Estimate ewk+qcd contribution")
parser.add_argument("--group_tt", action="store_true", help="Group all ttbar events together")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
analyzer = "HBSR_NoBTagSFs"
eos_dir = os.environ["eos_dir"]
plots_dir = os.environ["plots_dir"]

#set_trace()
base_ext = "*BKG*TOT.coffea"
input_dir = os.path.join(eos_dir, "results", f"{args.year}_{jobid}", analyzer)
fnames = fnmatch.filter(os.listdir(input_dir), base_ext)
fnames = [os.path.join(input_dir, fname) for fname in fnames]
hdict = load(fnames[0])


if args.bkg_shapes_comp:
    #set_trace()
    fnames = fnmatch.filter(os.listdir(os.path.join(eos_dir, "results", f"{args.year}_{jobid}", "htt_btag_sb_regions")), base_ext)
    fnames = [os.path.join(os.path.join(eos_dir, "results", f"{args.year}_{jobid}", "htt_btag_sb_regions"), fname) for fname in fnames]
    btagSFs_dict = load(fnames[0])

outdir = os.path.join(plots_dir, f"{args.year}_{jobid}", analyzer)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

jet_pars = prettyjson.loads(open(os.path.join(proj_dir, "cfg_files", f"cfg_pars_{jobid}.json")).read())["Jets"]
btagger = jet_pars["btagger"]
wps_to_use = list(set([jet_pars["permutations"]["tightb"],jet_pars["permutations"]["looseb"]]))
if not( len(wps_to_use) == 1):
    raise IOError("Only 1 unique btag working point supported now")
btag_wp = btagger+wps_to_use[0]
btag_cats = btag_sidebands.btag_cats[args.year]
btag_reg_names_dict = btag_sidebands.btag_reg_names_dict[args.year][btag_wp]

bkg_groups = {
    "EWQCD" : ["EWK", "QCD"],
    #"BKG" : ["EWK", "QCD"],
    "ttJets" : ["ttJets_right", "ttJets_matchable", "ttJets_unmatchable", "ttJets_sl_tau", "ttJets_other"],
    "singlet" : ["singlet"],
    "data" : ["data"],
}

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
    #"mtt" : (cfeatures.variable_names_to_labels["mtt"], 4, (200., 2000.), True),
    "tlep_ctstar_abs" : (cfeatures.variable_names_to_labels["tlep_ctstar_abs"], 1, (0., 1.), False),
    #"tlep_ctstar_abs" : (cfeatures.variable_names_to_labels["tlep_ctstar_abs"], 1, (0., 1.), True),
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
#names = [dataset for dataset in sorted(set([key[0] for key in hdict[sorted(variables.keys())[0]].values().keys()]))] # get dataset names in hists
ttJets_cats = [name for name in names if any([fnmatch.fnmatch(name, cat) for cat in ttJets_permcats])] # gets ttJets(_PS)_other, ...
if len(ttJets_cats) > 0:
    for tt_cat in ttJets_cats:
        ttJets_lumi_topo = "_".join(tt_cat.split("_")[:-2]) if "sl_tau" in tt_cat else "_".join(tt_cat.split("_")[:-1]) # gets ttJets[SL, Had, DiLep] or ttJets_PS
        ttJets_eff_lumi = lumi_correction[ttJets_lumi_topo]
        lumi_correction.update({tt_cat: ttJets_eff_lumi})

## make groups based on process
process = hist.Cat("process", "Process", sorting="placement")
process_cat = "dataset"
process_groups = plt_tools.make_dataset_groups(args.lepton, args.year, samples=names, bkgdict="dataset", sigdict="MC_indiv")

if args.group_tt:
    if not args.nosys: raise ValueError("This should only be run for --nosys")
    tt_categories = []
    tmp_groups = deepcopy(process_groups)
    for key, proc_list in tmp_groups.items():
        if "ttJets" not in key: continue
        [tt_categories.append(proc) for proc in proc_list]
        del process_groups[key]
    process_groups["ttJets"] = tt_categories

## make btagging sideband groups
btag_bin = hist.Cat("btag", "btag", sorting="placement")
btag_cat = "btag"
btag_groups = btag_sidebands.btag_groups[args.year]

    # scale and group hists by process
for hname in hdict.keys():
    if "cutflow" in hname: continue
    #set_trace()
    hdict[hname].scale(lumi_correction, axis="dataset") # scale hists
    hdict[hname] = hdict[hname].group(btag_cat, btag_bin, btag_groups) # group by btagging sideband regions
    hdict[hname] = hdict[hname].group(process_cat, process, process_groups) # group by process
    hdict[hname] = hdict[hname][:, :, :, :, args.lepton].integrate("leptype") # only pick out specified lepton

    if args.bkg_shapes_comp:
        btagSFs_dict[hname].scale(lumi_correction, axis="dataset") # scale hists
        btagSFs_dict[hname] = btagSFs_dict[hname].group(btag_cat, btag_bin, btag_groups) # group by btagging sideband regions
        btagSFs_dict[hname] = btagSFs_dict[hname].group(process_cat, process, process_groups) # group by process
        btagSFs_dict[hname] = btagSFs_dict[hname][:, :, :, :, args.lepton].integrate("leptype") # only pick out specified lepton


def plot_ewk_qcd_dmtop_shape(cen_sb={}, up_sb={}, dw_sb={}, **opts):
        # get opts
    xlimits = opts.get("xlims")
    xlabel = opts.get("xtitle")
    jmult = opts.get("jmult")
    pltdir = opts.get("pltdir")
    vlines = opts.get("vlines")
    overflow = opts.get("overflow", "none")

    fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True, figsize=(15.0, 10.0)) if vlines is not None else plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
    fig.subplots_adjust(hspace=.07)

    labels = []
    if cen_sb:
        cen_dmp_hist = cen_sb["histo"].integrate("process").copy()

            # get normalized arrays of data-promptMC    
        cen_dmp_sumw, cen_dmp_sumw2 = Plotter.get_qcd_shape(Plotter.data_minus_top(cen_sb["histo"], isForTemplates=False))

        cen_dmp_hist.values(overflow="all", sumw2=True)[()][0][:] = cen_dmp_sumw
        cen_dmp_hist.values(overflow="all", sumw2=True)[()][1][:] = cen_dmp_sumw2

        hist.plot1d(cen_dmp_hist, ax=ax, clear=False,
            line_opts={"color" : cen_sb["color"]},
            overflow=overflow,
            overlay_overflow=overflow,
        )
        hist.plot1d(cen_dmp_hist, ax=ax, clear=False,
            error_opts={"color" : cen_sb["color"]},
            overflow=overflow,
            overlay_overflow=overflow,
        )
        labels.append(cen_sb["label"])

    if up_sb:
        up_dmp_hist = up_sb["histo"].integrate("process").copy()

            # get normalized arrays of data-promptMC    
        up_dmp_sumw, up_dmp_sumw2 = Plotter.get_qcd_shape(Plotter.data_minus_top(up_sb["histo"], isForTemplates=False))

        up_dmp_hist.values(overflow="all", sumw2=True)[()][0][:] = up_dmp_sumw
        up_dmp_hist.values(overflow="all", sumw2=True)[()][1][:] = up_dmp_sumw2

        hist.plot1d(up_dmp_hist, ax=ax, clear=False,
            line_opts={"color" : up_sb["color"]},
            overflow=overflow,
            overlay_overflow=overflow,
        )
        hist.plot1d(up_dmp_hist, ax=ax, clear=False,
            error_opts={"color" : up_sb["color"]},
            overflow=overflow,
            overlay_overflow=overflow,
        )
        labels.append(up_sb["label"])

            # plot num/denom ratio
        hist.plotratio(
            up_dmp_hist, cen_dmp_hist, error_opts={"marker": ".", "markersize": 10., "color": up_sb["color"], "elinewidth": 1},
            unc="num", clear=False, ax=rax, guide_opts={},
            overflow=overflow,
            #overlay_overflow=overflow,
        )

    if dw_sb:
        dw_dmp_hist = dw_sb["histo"].integrate("process").copy()

            # get normalized arrays of data-promptMC    
        dw_dmp_sumw, dw_dmp_sumw2 = Plotter.get_qcd_shape(Plotter.data_minus_top(dw_sb["histo"], isForTemplates=False))

        dw_dmp_hist.values(overflow="all", sumw2=True)[()][0][:] = dw_dmp_sumw
        dw_dmp_hist.values(overflow="all", sumw2=True)[()][1][:] = dw_dmp_sumw2

        hist.plot1d(dw_dmp_hist, ax=ax, clear=False,
            line_opts={"color" : dw_sb["color"]},
            overflow=overflow,
            overlay_overflow=overflow,
        )
        hist.plot1d(dw_dmp_hist, ax=ax, clear=False,
            error_opts={"color" : dw_sb["color"]},
            overflow=overflow,
            overlay_overflow=overflow,
        )
        labels.append(dw_sb["label"])
        
            # plot num/denom ratio
        hist.plotratio(
            dw_dmp_hist, cen_dmp_hist, error_opts={"marker": ".", "markersize": 10., "color": dw_sb["color"], "elinewidth": 1},
            unc="num", clear=False, ax=rax, guide_opts={},
            overflow=overflow,
            #overlay_overflow=overflow,
        )


    #set_trace()
        ## set legend and corresponding colors
    handles, def_labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc="upper right", title="data-$MC_{st/t\\bar{t}}$")

    ax.autoscale()
    ax.set_ylabel("Probability Density")
    ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.5)
    #ax.set_ylim(0, None)
    ax.set_xlabel(None)
    ax.set_xlim(xlimits)

    rax.autoscale()
    rax.set_ylabel("Sys/Cen")
    rax.set_ylim(max(0.0, rax.get_ylim()[0]), min(1.5, rax.get_ylim()[1]))
    #rax.set_ylim(max(0.0, rax.get_ylim()[0]), min(5.0, rax.get_ylim()[1]))
    rax.set_xlim(xlimits)
    rax.set_xlabel(xlabel)

    if vlines is not None:
        #set_trace()
            # draw vertical lines separating ctstar bins
        for vline in vlines:
            ax.axvline(vline, color="k", linestyle="--")
            rax.axvline(vline, color="k", linestyle="--")

    if "ytext" in opts.keys():
        if opts["ytext"] is not None:
            ybinlabels, ybin_locs = opts["ytext"]
            #set_trace()
            for idx, label in enumerate(ybinlabels):
                ax.annotate(label, xy=(ybin_locs[idx], 0), xycoords=("data", "axes fraction"),
                    xytext=(0, 280), textcoords="offset points", va="top", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)
                    #xytext=(0, 250), textcoords="offset points", va="top", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)
    if "xtext" in opts.keys():
        if opts["xtext"] is not None:
            #set_trace()
            xbin_inds, xbin_labels = opts["xtext"]
            rax.set_xticks(xbin_inds)
            rax.set_xticklabels(xbin_labels)

    ax.text(
        0.02, 0.88, f"{cfeatures.channel_labels[f'{args.lepton}_{jmult}']}",
        fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
    )
    hep.cms.label(ax=ax, data=withData, label="Preliminary", year=year_to_use, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))

    #set_trace()
    if not os.path.isdir(pltdir):
        os.makedirs(pltdir)
    pltname = opts.get("fname")
    figname = os.path.join(pltdir, pltname)
    fig.savefig(figname)
    print(f"{figname} written")
    plt.close(fig)



def compare_dmtop_shapes_btagSFs(cen_sb_btagSF={}, up_sb_btagSF={}, dw_sb_btagSF={}, cen_sb_nobtagSF={}, up_sb_nobtagSF={}, dw_sb_nobtagSF={}, **opts):
    #set_trace()
        # get opts
    xlimits = opts.get("xlims")
    xlabel = opts.get("xtitle")
    jmult = opts.get("jmult")
    pltdir = opts.get("pltdir")
    vlines = opts.get("vlines")
    overflow = opts.get("overflow", "none")

    fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True, figsize=(15.0, 10.0)) if vlines is not None else plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
    fig.subplots_adjust(hspace=.07)

    labels = []
    if cen_sb_btagSF:
        cen_btagSF_dmt_hist = cen_sb_btagSF["histo"].integrate("process").copy()

            # get normalized arrays of data-promptMC    
        cen_btagSF_dmt_sumw, cen_btagSF_dmt_sumw2 = Plotter.get_qcd_shape(Plotter.data_minus_top(cen_sb_btagSF["histo"], isForTemplates=False))

        cen_btagSF_dmt_hist.values(overflow="all", sumw2=True)[()][0][:] = cen_btagSF_dmt_sumw
        cen_btagSF_dmt_hist.values(overflow="all", sumw2=True)[()][1][:] = cen_btagSF_dmt_sumw2

        hist.plot1d(cen_btagSF_dmt_hist, ax=ax, clear=False,
            line_opts={"color" : cen_sb_btagSF["color"]},
            overflow=overflow,
            overlay_overflow=overflow,
        )
        hist.plot1d(cen_btagSF_dmt_hist, ax=ax, clear=False,
            error_opts={"color" : cen_sb_btagSF["color"]},
            overflow=overflow,
            overlay_overflow=overflow,
        )
        labels.append(cen_sb_btagSF["label"])

    if up_sb_btagSF:
        up_btagSF_dmt_hist = up_sb_btagSF["histo"].integrate("process").copy()

            # get normalized arrays of data-promptMC    
        up_btagSF_dmt_sumw, up_btagSF_dmt_sumw2 = Plotter.get_qcd_shape(Plotter.data_minus_top(up_sb_btagSF["histo"], isForTemplates=False))

        up_btagSF_dmt_hist.values(overflow="all", sumw2=True)[()][0][:] = up_btagSF_dmt_sumw
        up_btagSF_dmt_hist.values(overflow="all", sumw2=True)[()][1][:] = up_btagSF_dmt_sumw2

        hist.plot1d(up_btagSF_dmt_hist, ax=ax, clear=False,
            line_opts={"color" : up_sb_btagSF["color"]},
            overflow=overflow,
            overlay_overflow=overflow,
        )
        hist.plot1d(up_btagSF_dmt_hist, ax=ax, clear=False,
            error_opts={"color" : up_sb_btagSF["color"]},
            overflow=overflow,
            overlay_overflow=overflow,
        )
        labels.append(up_sb_btagSF["label"])

            # plot num/denom ratio
        hist.plotratio(
            up_btagSF_dmt_hist, cen_btagSF_dmt_hist, error_opts={"marker": ".", "markersize": 10., "color": up_sb_btagSF["color"], "elinewidth": 1},
            unc="num", clear=False, ax=rax, guide_opts={},
            overflow=overflow,
        )

    if dw_sb_btagSF:
        dw_btagSF_dmt_hist = dw_sb_btagSF["histo"].integrate("process").copy()

            # get normalized arrays of data-promptMC    
        dw_btagSF_dmt_sumw, dw_btagSF_dmt_sumw2 = Plotter.get_qcd_shape(Plotter.data_minus_top(dw_sb_btagSF["histo"], isForTemplates=False))

        dw_btagSF_dmt_hist.values(overflow="all", sumw2=True)[()][0][:] = dw_btagSF_dmt_sumw
        dw_btagSF_dmt_hist.values(overflow="all", sumw2=True)[()][1][:] = dw_btagSF_dmt_sumw2

        hist.plot1d(dw_btagSF_dmt_hist, ax=ax, clear=False,
            line_opts={"color" : dw_sb_btagSF["color"]},
            overflow=overflow,
            overlay_overflow=overflow,
        )
        hist.plot1d(dw_btagSF_dmt_hist, ax=ax, clear=False,
            error_opts={"color" : dw_sb_btagSF["color"]},
            overflow=overflow,
            overlay_overflow=overflow,
        )
        labels.append(dw_sb_btagSF["label"])

            # plot num/denom ratio
        hist.plotratio(
            dw_btagSF_dmt_hist, cen_btagSF_dmt_hist, error_opts={"marker": ".", "markersize": 10., "color": dw_sb_btagSF["color"], "elinewidth": 1},
            unc="num", clear=False, ax=rax, guide_opts={},
            overflow=overflow,
        )

    if cen_sb_nobtagSF:
        cen_nobtagSF_dmt_hist = cen_sb_nobtagSF["histo"].integrate("process").copy()

            # get normalized arrays of data-promptMC    
        cen_nobtagSF_dmt_sumw, cen_nobtagSF_dmt_sumw2 = Plotter.get_qcd_shape(Plotter.data_minus_top(cen_sb_nobtagSF["histo"], isForTemplates=False))

        cen_nobtagSF_dmt_hist.values(overflow="all", sumw2=True)[()][0][:] = cen_nobtagSF_dmt_sumw
        cen_nobtagSF_dmt_hist.values(overflow="all", sumw2=True)[()][1][:] = cen_nobtagSF_dmt_sumw2

        hist.plot1d(cen_nobtagSF_dmt_hist, ax=ax, clear=False,
            line_opts={"color" : cen_sb_nobtagSF["color"]},
            overflow=overflow,
            overlay_overflow=overflow,
        )
        hist.plot1d(cen_nobtagSF_dmt_hist, ax=ax, clear=False,
            error_opts={"color" : cen_sb_nobtagSF["color"]},
            overflow=overflow,
            overlay_overflow=overflow,
        )
        labels.append(cen_sb_nobtagSF["label"])

    if up_sb_nobtagSF:
        up_nobtagSF_dmt_hist = up_sb_nobtagSF["histo"].integrate("process").copy()

            # get normalized arrays of data-promptMC    
        up_nobtagSF_dmt_sumw, up_nobtagSF_dmt_sumw2 = Plotter.get_qcd_shape(Plotter.data_minus_top(up_sb_nobtagSF["histo"], isForTemplates=False))

        up_nobtagSF_dmt_hist.values(overflow="all", sumw2=True)[()][0][:] = up_nobtagSF_dmt_sumw
        up_nobtagSF_dmt_hist.values(overflow="all", sumw2=True)[()][1][:] = up_nobtagSF_dmt_sumw2

        hist.plot1d(up_nobtagSF_dmt_hist, ax=ax, clear=False,
            line_opts={"color" : up_sb_nobtagSF["color"]},
            overflow=overflow,
            overlay_overflow=overflow,
        )
        hist.plot1d(up_nobtagSF_dmt_hist, ax=ax, clear=False,
            error_opts={"color" : up_sb_nobtagSF["color"]},
            overflow=overflow,
            overlay_overflow=overflow,
        )
        labels.append(up_sb_nobtagSF["label"])

            # plot num/denom ratio
        hist.plotratio(
            up_nobtagSF_dmt_hist, cen_nobtagSF_dmt_hist, error_opts={"marker": ".", "markersize": 10., "color": up_sb_nobtagSF["color"], "elinewidth": 1},
            unc="num", clear=False, ax=rax, guide_opts={},
            overflow=overflow,
        )

    if dw_sb_nobtagSF:
        dw_nobtagSF_dmt_hist = dw_sb_nobtagSF["histo"].integrate("process").copy()

            # get normalized arrays of data-promptMC    
        dw_nobtagSF_dmt_sumw, dw_nobtagSF_dmt_sumw2 = Plotter.get_qcd_shape(Plotter.data_minus_top(dw_sb_nobtagSF["histo"], isForTemplates=False))

        dw_nobtagSF_dmt_hist.values(overflow="all", sumw2=True)[()][0][:] = dw_nobtagSF_dmt_sumw
        dw_nobtagSF_dmt_hist.values(overflow="all", sumw2=True)[()][1][:] = dw_nobtagSF_dmt_sumw2

        hist.plot1d(dw_nobtagSF_dmt_hist, ax=ax, clear=False,
            line_opts={"color" : dw_sb_nobtagSF["color"]},
            overflow=overflow,
            overlay_overflow=overflow,
        )
        hist.plot1d(dw_nobtagSF_dmt_hist, ax=ax, clear=False,
            error_opts={"color" : dw_sb_nobtagSF["color"]},
            overflow=overflow,
            overlay_overflow=overflow,
        )
        labels.append(dw_sb_nobtagSF["label"])

            # plot num/denom ratio
        hist.plotratio(
            dw_nobtagSF_dmt_hist, cen_nobtagSF_dmt_hist, error_opts={"marker": ".", "markersize": 10., "color": dw_sb_nobtagSF["color"], "elinewidth": 1},
            unc="num", clear=False, ax=rax, guide_opts={},
            overflow=overflow,
        )

    #set_trace()
        ## set legend and corresponding colors
    handles, def_labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc="upper right", title="data-$MC_{st/t\\bar{t}}$", ncol=2)

    ax.autoscale()
    ax.set_ylabel("Probability Density")
    ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.5)
    #ax.set_ylim(0, None)
    ax.set_xlabel(None)
    ax.set_xlim(xlimits)

    #rax.set_yscale("log")
    rax.autoscale()
    rax.set_ylabel("Sys/Cen")
    #rax.set_ylim(max(1e-5, rax.get_ylim()[0]), min(1.5, rax.get_ylim()[1]))
    rax.set_ylim(max(0.0, rax.get_ylim()[0]), min(1.5, rax.get_ylim()[1]))
    rax.set_xlim(xlimits)
    rax.set_xlabel(xlabel)

    if vlines is not None:
        #set_trace()
            # draw vertical lines separating ctstar bins
        for vline in vlines:
            ax.axvline(vline, color="k", linestyle="--")
            rax.axvline(vline, color="k", linestyle="--")

    if "ytext" in opts.keys():
        if opts["ytext"] is not None:
            ybinlabels, ybin_locs = opts["ytext"]
            #set_trace()
            for idx, label in enumerate(ybinlabels):
                ax.annotate(label, xy=(ybin_locs[idx], 0), xycoords=("data", "axes fraction"),
                    xytext=(0, 280), textcoords="offset points", va="top", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)
                    #xytext=(0, 250), textcoords="offset points", va="top", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)
    if "xtext" in opts.keys():
        if opts["xtext"] is not None:
            #set_trace()
            xbin_inds, xbin_labels = opts["xtext"]
            rax.set_xticks(xbin_inds)
            rax.set_xticklabels(xbin_labels)

    ax.text(
        0.02, 0.88, f"{cfeatures.channel_labels[f'{args.lepton}_{jmult}']}",
        fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
    )
    hep.cms.label(ax=ax, data=withData, label="Preliminary", year=year_to_use, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))

    #set_trace()
    if not os.path.isdir(pltdir):
        os.makedirs(pltdir)
    pltname = opts.get("fname")
    figname = os.path.join(pltdir, pltname)
    fig.savefig(figname)
    print(f"{figname} written")
    plt.close(fig)



    ## make plots
if args.nosys:
    for hname in variables.keys():
        if hname not in hdict.keys():
            print(f"{hname} not found in file")
            continue
            #raise ValueError(f"{hname} not found in file")

        #set_trace()
        histo = hdict[hname][:, :, "nosys", :].integrate("sys") # process, sys, jmult, btag
        #set_trace()
    
        if histo.dense_dim() == 1:
            xtitle, rebinning, x_lims, withData = variables[hname]
            orig_xtitle = xtitle
            if rebinning != 1:
                xaxis_name = histo.dense_axes()[0].name
                histo = histo.rebin(xaxis_name, rebinning)
    
            ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
            #for jmult in ["3Jets"]:
            #for jmult in ["4PJets"]:
            for jmult in sorted(set([key[2] for key in histo.values().keys()])):
                #for btagregion in ["btagPass"]:
                for btagregion in sorted(set([key[1] for key in histo.values().keys()]) - set(["btagPass"])):
                #for btagregion in sorted(set([key[1] for key in histo.values().keys()])):
                    print(*[jmult, btagregion, hname], sep=", ")
                    pltdir = os.path.join(outdir, args.lepton, jmult, btagregion)
                    if not os.path.isdir(pltdir):
                        os.makedirs(pltdir)
   
                    fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                    fig.subplots_adjust(hspace=.07)
    
                    #set_trace()
                    hslice = histo[:, btagregion, jmult].integrate("jmult").integrate("btag")
                    nonsignal_hslice = hslice[Plotter.nonsignal_samples]

                    if "disc" in hname:
                        xtitle = orig_xtitle.replace("j", "3") if jmult == "3Jets" else orig_xtitle.replace("j", "4")
   
                    if hname == "Lep_iso":
                        x_lims = (0., 0.15) if (args.lepton == "Muon") else (0., 0.1)
    
                    mc_opts = {
                        "mcorder" : ["ttJets", "singlet", "EWK", "QCD"] if args.group_tt else ["ttJets_right", "ttJets_matchable", "ttJets_unmatchable", "ttJets_sl_tau", "ttJets_other", "singlet", "EWK", "QCD"],
                        "maskData" : not withData if btagregion == "btagPass" else False,
                        #"maskData" : not withData,
                        "overflow" : "under" if "DeepCSV_bDisc" in hname else "none",
                    }

                    Plotter.plot_stack1d(ax, rax, nonsignal_hslice, xlabel=xtitle, xlimits=x_lims, **mc_opts)
    
                    #set_trace() 
                    if hname == "Jets_njets":
                        print(jmult)
                        yields_txt, yields_json = Plotter.get_samples_yield_and_frac(hslice, data_lumi_year[f"{args.lepton}s"]/1000., promptmc=True)
                        frac_name = os.path.join(pltdir, f"{jmult}_{args.lepton}_{btagregion}_yields_and_fracs_CombinedTT.txt") if args.group_tt \
                            else os.path.join(pltdir, f"{jmult}_{args.lepton}_{btagregion}_yields_and_fracs.txt")
                        plt_tools.print_table(yields_txt, filename=frac_name, print_output=True)
                        print(f"{frac_name} written")
                        with open(frac_name.replace(".txt", ".json"), "w") as out:
                            out.write(prettyjson.dumps(yields_json))
    
                        # add lepton/jet multiplicity label
                    #set_trace()
                    ax.text(
                        0.02, 0.84, f"{cfeatures.channel_labels[f'{args.lepton}_{jmult}']}\n{btag_cats[btagregion]}",
                        fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                    )
                    hep.cms.label(ax=ax, data=withData, label="Preliminary", year=year_to_use, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))
    
                    #set_trace()
                    figname = os.path.join(pltdir, "_".join([jmult, args.lepton, btagregion, hname, "CombinedTT"])) if args.group_tt else os.path.join(pltdir, "_".join([jmult, args.lepton, btagregion, hname]))
                    fig.savefig(figname)
                    print(f"{figname} written")
                    plt.close(fig)
        
        if histo.dense_dim() == 2:
            xtitle, ytitle, xrebinning, yrebinning, x_lims, y_lims, withData, plot_xlabels, plot_ylabels = variables[hname]

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
   
            ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
            #for jmult in ["3Jets"]:
            ##for jmult in ["4PJets"]:
            #   for btagregion in ["btagPass"]:
            for jmult in sorted(set([key[2] for key in histo.values().keys()])):
                #for btagregion in ["btagPass"]:
                for btagregion in sorted(set([key[1] for key in histo.values().keys()]) - set(["btagPass"])):
                #for btagregion in sorted(set([key[1] for key in histo.values().keys()])):
                    print(*[jmult, btagregion, hname], sep=", ")
                    pltdir = os.path.join(outdir, args.lepton, jmult, btagregion)
                    if not os.path.isdir(pltdir):
                        os.makedirs(pltdir)
  
                    if (btagregion == "btagPass") and (hname == "mtt_vs_tlep_ctstar_abs"): 
                        withData = False
                    else:
                        withData = variables[hname][-1]

                    hslice = histo[:, btagregion, jmult].integrate("jmult").integrate("btag")
                    hslice = hslice[Plotter.nonsignal_samples]
   
                    mc_opts = {
                        "mcorder" : ["ttJets", "singlet", "EWK", "QCD"] if args.group_tt else ["ttJets_right", "ttJets_matchable", "ttJets_unmatchable", "ttJets_sl_tau", "ttJets_other", "singlet", "EWK", "QCD"],
                        "maskData" : not withData,
                    }

                    if "DeepCSV_bDisc" in hname:
                        underflow_hist = hslice[:, :, :0.00].integrate("deepcsv_bdisc", overflow="under")
                        normal_reg_hist = hslice.integrate("deepcsv_bdisc")

                            # plot 1D projection of underflow bin
                        fig_und, (ax_und, rax_und) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                        fig_und.subplots_adjust(hspace=.07)

                        Plotter.plot_stack1d(ax_und, rax_und, underflow_hist, xlabel="%s DeepCSV underflow" % xtitle, xlimits=x_lims, **mc_opts)
    
                            # add lepton/jet multiplicity label
                        ax_und.text(
                            0.02, 0.84, f"{cfeatures.channel_labels[f'{args.lepton}_{jmult}']}\n{btag_cats[btagregion]}",
                            fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax_und.transAxes
                        )
                        hep.cms.label(ax=ax_und, data=withData, label="Preliminary", year=year_to_use, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))
    
                        #set_trace()
                        und_figname = os.path.join(pltdir, "_".join([jmult, args.lepton, btagregion, hname, "DeepCSV_underflow", "proj", "CombinedTT"])) if args.group_tt \
                            else os.path.join(pltdir, "_".join([jmult, args.lepton, btagregion, hname, "DeepCSV_underflow", "proj"]))
                        fig_und.savefig(und_figname)
                        print(f"{und_figname} written")
                        plt.close(fig_und)
 
                            # plot 1D projection of normal bin region
                        fig_norm, (ax_norm, rax_norm) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True, figsize=(15.0, 10.0))
                        fig_norm.subplots_adjust(hspace=.07)

                        Plotter.plot_stack1d(ax_norm, rax_norm, normal_reg_hist, xlabel="%s DeepCSV normal region" % xtitle, xlimits=x_lims, **mc_opts)
    
                            # add lepton/jet multiplicity label
                        ax_norm.text(
                            0.02, 0.84, f"{cfeatures.channel_labels[f'{args.lepton}_{jmult}']}\n{btag_cats[btagregion]}",
                            fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax_norm.transAxes
                        )
                        hep.cms.label(ax=ax_norm, data=withData, label="Preliminary", year=year_to_use, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))
    
                        #set_trace()
                        norm_figname = os.path.join(pltdir, "_".join([jmult, args.lepton, btagregion, hname, "DeepCSV_normal", "proj", "CombinedTT"])) if args.group_tt \
                            else os.path.join(pltdir, "_".join([jmult, args.lepton, btagregion, hname, "DeepCSV_normal", "proj"]))
                        fig_norm.savefig(norm_figname)
                        print(f"{norm_figname} written")
                        plt.close(fig_norm)
 
                        continue


                        # make 1D projection along dense axes
                    for dax in range(2):
                        if dax == 0:
                            xlabel = ytitle
                            xlimits = y_lims
                            figname = os.path.join(pltdir, "_".join([jmult, args.lepton, btagregion, yaxis_name, "proj", "CombinedTT"])) if args.group_tt \
                                else os.path.join(pltdir, "_".join([jmult, args.lepton, btagregion, yaxis_name, "proj"]))
                        else:
                            xlabel = xtitle
                            xlimits = x_lims
                            figname = os.path.join(pltdir, "_".join([jmult, args.lepton, btagregion, xaxis_name, "proj", "CombinedTT"])) if args.group_tt \
                                else os.path.join(pltdir, "_".join([jmult, args.lepton, btagregion, xaxis_name, "proj"]))

                        fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                        fig.subplots_adjust(hspace=.07)

                        hproj = hslice.integrate(hslice.dense_axes()[dax].name)
    
                        Plotter.plot_stack1d(ax, rax, hproj, xlabel=xlabel, xlimits=xlimits, **mc_opts)
    
                            # add lepton/jet multiplicity label
                        ax.text(
                            0.02, 0.84, f"{cfeatures.channel_labels[f'{args.lepton}_{jmult}']}\n{btag_cats[btagregion]}",
                            fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                        )
                        hep.cms.label(ax=ax, data=withData, label="Preliminary", year=year_to_use, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))
    
                        #set_trace()
                        fig.savefig(figname)
                        print(f"{figname} written")
                        plt.close(fig)
 
                        # plot linearized view 
                    fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True, figsize=(15.0, 10.0))
                    fig.subplots_adjust(hspace=.07)

                    hline = Plotter.linearize_hist(hslice)
                    
                    Plotter.plot_stack1d(ax, rax, hline, xlabel=xtitle, xlimits=(0, len(hline.axis(hline.dense_axes()[0].name).edges())-1), **mc_opts)
    
                        # draw vertical lines separating ctstar bins
                    bin_sep_lines = [hslice.values()[("EWK",)].shape[0]*ybin for ybin in range(1, hslice.values()[("EWK",)].shape[1])]
                    for binline in bin_sep_lines:
                        ax.axvline(binline, color="k", linestyle="--")
                        if rax is not None: rax.axvline(binline, color="k", linestyle="--")

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

                        if hname == "mtt_vs_tlep_ctstar_abs":
                            #set_trace()
                            rax.set_xticks(mtt_bin_inds_to_plot)
                            rax.set_xticklabels(mtt_bins_to_plot)

                        # add lepton/jet multiplicity label
                    ax.text(
                        0.02, 0.84, f"{cfeatures.channel_labels[f'{args.lepton}_{jmult}']}\n{btag_cats[btagregion]}",
                        fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                    )
                    hep.cms.label(ax=ax, data=withData, label="Preliminary", year=year_to_use, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))
    
                    #set_trace()
                    figname = os.path.join(pltdir, "_".join([jmult, args.lepton, btagregion, hname, "CombinedTT"])) if args.group_tt \
                        else os.path.join(pltdir, "_".join([jmult, args.lepton, btagregion, hname]))
                    fig.savefig(figname)
                    print(f"{figname} written")
                    plt.close(fig)



if args.bkg_shapes:
    for hname in variables.keys():
        #if hname != "mtt_vs_tlep_ctstar_abs": continue
        if hname not in hdict.keys():
            print(f"{hname} not found in file")
            continue
            #raise ValueError(f"{hname} not found in file")
        histo = hdict[hname][:, :, "nosys", :].integrate("sys") # process, sys, jmult, btag
        #set_trace()
        
        if histo.dense_dim() == 1:
            xtitle, xrebinning, x_lims, withData = variables[hname]
            orig_xtitle = xtitle
            xaxis_name = histo.dense_axes()[0].name
                ## rebin x axis
            if isinstance(xrebinning, np.ndarray):
                new_xbins = hist.Bin(xaxis_name, xaxis_name, xrebinning)
            elif isinstance(xrebinning, float) or isinstance(xrebinning, int):
                new_xbins = xrebinning
            histo = histo.rebin(xaxis_name, new_xbins)
    
        if histo.dense_dim() == 2:
            xtitle, ytitle, xrebinning, yrebinning, x_lims, y_lims, withData, plot_xlabels, plot_ylabels = variables[hname]
            if (hname == "mtt_vs_tlep_ctstar_abs"):
                #set_trace()
                xrebinning = np.array([linearize_binning[0][0], 400., 500.0, 600.0, 800.0, 1000.0, 1500.0, linearize_binning[0][-1]])
                x_lims = (xrebinning[0], xrebinning[-1])
                ctstar_binlabels = [r"|$\cos (\theta^{*}_{t_{l}})$| $\in [%s, %s)$" % (yrebinning[bin], yrebinning[bin+1]) for bin in range(len(yrebinning)-1)]
                ctstar_binlabels[-1] = r"|$\cos (\theta^{*}_{t_{l}})$| $\in [%s, %s]$" % (yrebinning[-2], yrebinning[-1])
                tmp_ctstar_bin_locs = np.linspace((len(xrebinning)-1)/2, (len(xrebinning)-1)*(len(yrebinning)-1) - (len(xrebinning)-1)/2, len(yrebinning)-1)
                #mtt_vals_to_plot = np.array([400, 600, 1000])
                tmp_mtt_tiled_labels = np.tile(xrebinning[:-1], yrebinning.size-1)
                tmp_mtt_bin_inds_to_plot = np.where(np.in1d(tmp_mtt_tiled_labels, mtt_vals_to_plot))[0]
                tmp_mtt_bins_to_plot = np.tile(mtt_vals_to_plot, yrebinning.size-1)

            if "DeepCSV_bDisc" not in hname:
                vlines = [(len(xrebinning)-1)*ybin for ybin in range(1, len(yrebinning)-1)] if histo.dense_dim() == 2 else None
                nbins = (len(xrebinning)-1)*(len(yrebinning)-1) if histo.dense_dim() == 2 else None
    
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

        #set_trace()   
        if hname == "Lep_iso":
            x_lims = (0., 0.15) if (args.lepton == "Muon") else (0., 0.1)
    
        ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
        #for jmult in ["3Jets"]:
        #for jmult in ["4PJets"]:
        for jmult in sorted(set([key[2] for key in histo.values().keys()])):
            if "disc" in hname:
                xtitle = orig_xtitle.replace("j", "3") if jmult == "3Jets" else orig_xtitle.replace("j", "4")
   
            shape_dir = os.path.join(outdir, args.lepton, jmult, "BKG_Est_orthog", "Shapes", "DMT")
            if not os.path.isdir(shape_dir): os.makedirs(shape_dir)
            #comp_dir = os.path.join(outdir, args.lepton, jmult, "BKG_Est_orthog", "Shapes", "MC_DD_Comp")
            #if not os.path.isdir(comp_dir): os.makedirs(comp_dir)

                # get sideband and signal region hists
            cen_sb_histo = histo[:, btag_reg_names_dict["Central"]["reg"], jmult].integrate("jmult").integrate("btag")
            up_sb_histo = histo[:, btag_reg_names_dict["Up"]["reg"], jmult].integrate("jmult").integrate("btag")
            dw_sb_histo = histo[:, btag_reg_names_dict["Down"]["reg"], jmult].integrate("jmult").integrate("btag")
            #sig_histo = histo[:, btag_reg_names_dict["Signal"]["reg"], jmult].integrate("jmult").integrate("btag")

            if histo.dense_dim() == 1:
                # make data-(single top/ttbar) shape plots
                plot_ewk_qcd_dmtop_shape(
                    cen_sb={"histo" : cen_sb_histo, "label" : btag_reg_names_dict["Central"]["label"], "color" : btag_reg_names_dict["Central"]["color"]},
                    up_sb={"histo" : up_sb_histo, "label" : btag_reg_names_dict["Up"]["label"], "color" : btag_reg_names_dict["Up"]["color"]},
                    dw_sb={"histo" : dw_sb_histo, "label" : btag_reg_names_dict["Down"]["label"], "color" : btag_reg_names_dict["Down"]["color"]},
                   **{"xtitle":xtitle, "xlims":x_lims, "pltdir":shape_dir, "hname": hname, "jmult":jmult, 
                        "fname": "BKG_DD_Shapes_%s_%s_%s" % (args.lepton, jmult, hname), "overflow" : "under" if "DeepCSV_bDisc" in hname else "none"})

            if histo.dense_dim() == 2:
                if "DeepCSV_bDisc" in hname:
                        ## plot 1D projection of underflow bin
                    plot_ewk_qcd_dmtop_shape(
                        cen_sb={"histo" : cen_sb_histo[:, :, :0.00].integrate("deepcsv_bdisc", overflow="under"), "label" : btag_reg_names_dict["Central"]["label"], "color" : btag_reg_names_dict["Central"]["color"]},
                        up_sb={"histo" : up_sb_histo[:, :, :0.00].integrate("deepcsv_bdisc", overflow="under"), "label" : btag_reg_names_dict["Up"]["label"], "color" : btag_reg_names_dict["Up"]["color"]},
                        dw_sb={"histo" : dw_sb_histo[:, :, :0.00].integrate("deepcsv_bdisc", overflow="under"), "label" : btag_reg_names_dict["Down"]["label"], "color" : btag_reg_names_dict["Down"]["color"]},
                       **{"xtitle":"%s DeepCSV underflow" % xtitle, "xlims":x_lims, "pltdir":shape_dir, "jmult":jmult, "fname": "_".join(["BKG_DD_Shapes", args.lepton, jmult, hname, "DeepCSV_underflow", "proj"])})

                        ## plot 1D projection of normal bin content
                    plot_ewk_qcd_dmtop_shape(
                        cen_sb={"histo" : cen_sb_histo.integrate("deepcsv_bdisc"), "label" : btag_reg_names_dict["Central"]["label"], "color" : btag_reg_names_dict["Central"]["color"]},
                        up_sb={"histo" : up_sb_histo.integrate("deepcsv_bdisc"), "label" : btag_reg_names_dict["Up"]["label"], "color" : btag_reg_names_dict["Up"]["color"]},
                        dw_sb={"histo" : dw_sb_histo.integrate("deepcsv_bdisc"), "label" : btag_reg_names_dict["Down"]["label"], "color" : btag_reg_names_dict["Down"]["color"]},
                       **{"xtitle":"%s DeepCSV normal region" % xtitle, "xlims":x_lims, "pltdir":shape_dir, "jmult":jmult, "fname": "_".join(["BKG_DD_Shapes", args.lepton, jmult, hname, "DeepCSV_normal", "proj"])})

                    continue


                    # linearize sideband and signal region hists (hardcoded)
                cen_sb_lin = Plotter.linearize_hist(cen_sb_histo)
                up_sb_lin = Plotter.linearize_hist(up_sb_histo)
                dw_sb_lin = Plotter.linearize_hist(dw_sb_histo)
                #sig_lin = Plotter.linearize_hist(sig_histo)

                # make data-(single top/ttbar) shape plots
                plot_ewk_qcd_dmtop_shape(
                    cen_sb={"histo" : cen_sb_lin, "label" : btag_reg_names_dict["Central"]["label"], "color" : btag_reg_names_dict["Central"]["color"]},
                    up_sb={"histo" : up_sb_lin, "label" : btag_reg_names_dict["Up"]["label"], "color" : btag_reg_names_dict["Up"]["color"]},
                    dw_sb={"histo" : dw_sb_lin, "label" : btag_reg_names_dict["Down"]["label"], "color" : btag_reg_names_dict["Down"]["color"]},
                   **{"xtitle": xtitle, "xlims":(0, nbins), "pltdir":shape_dir, "hname": hname, "jmult":jmult, "vlines":vlines, "fname": f"BKG_DD_Shapes_{args.lepton}_{jmult}_{hname}",
                   #**{"xtitle":f"{xtitle} [GeV]", "xlims":(0, nbins), "pltdir":shape_dir, "hname": hname, "jmult":jmult, "vlines":vlines, "fname": f"BKG_DD_Shapes_{args.lepton}_{jmult}_{hname}",
                        "ytext":(ctstar_binlabels, tmp_ctstar_bin_locs) if hname == "mtt_vs_tlep_ctstar_abs" else None, "xtext":(tmp_mtt_bin_inds_to_plot, tmp_mtt_bins_to_plot) if hname == "mtt_vs_tlep_ctstar_abs" else None})
                   #**{"xtitle":"%s $\otimes$ %s" % (xtitle, ytitle), "xlims":(0, nbins), "pltdir":shape_dir, "hname": hname, "jmult":jmult, "vlines":vlines, "fname": "BKG_DD_Shapes_%s_%s_%s" % (args.lepton, jmult, hname)})

                # make 1D projection along dense axes
                for dax in range(2):
                    if dax == 0:
                        ## plot projection onto yaxis
                        plot_ewk_qcd_dmtop_shape(
                            cen_sb={"histo" : cen_sb_histo.integrate(cen_sb_histo.dense_axes()[dax].name), "label" : btag_reg_names_dict["Central"]["label"], "color" : btag_reg_names_dict["Central"]["color"]},
                            up_sb={"histo" : up_sb_histo.integrate(up_sb_histo.dense_axes()[dax].name), "label" : btag_reg_names_dict["Up"]["label"], "color" : btag_reg_names_dict["Up"]["color"]},
                            dw_sb={"histo" : dw_sb_histo.integrate(dw_sb_histo.dense_axes()[dax].name), "label" : btag_reg_names_dict["Down"]["label"], "color" : btag_reg_names_dict["Down"]["color"]},
                           **{"xtitle":ytitle, "xlims":y_lims, "pltdir":shape_dir, "jmult":jmult, "fname": "BKG_DD_Shapes_%s_%s_%s_proj" % (args.lepton, jmult, hname.split("_vs_")[1])})
                    else:
                        ## plot projection onto xaxis
                        plot_ewk_qcd_dmtop_shape(
                            cen_sb={"histo" : cen_sb_histo.integrate(cen_sb_histo.dense_axes()[dax].name), "label" : btag_reg_names_dict["Central"]["label"], "color" : btag_reg_names_dict["Central"]["color"]},
                            up_sb={"histo" : up_sb_histo.integrate(up_sb_histo.dense_axes()[dax].name), "label" : btag_reg_names_dict["Up"]["label"], "color" : btag_reg_names_dict["Up"]["color"]},
                            dw_sb={"histo" : dw_sb_histo.integrate(dw_sb_histo.dense_axes()[dax].name), "label" : btag_reg_names_dict["Down"]["label"], "color" : btag_reg_names_dict["Down"]["color"]},
                           **{"xtitle":xtitle, "xlims":x_lims, "pltdir":shape_dir, "jmult":jmult, "fname": "BKG_DD_Shapes_%s_%s_%s_proj" % (args.lepton, jmult, hname.split("_vs_")[0])})


if args.bkg_shapes_comp:
    for hname in variables.keys():
        #if hname != "mtt_vs_tlep_ctstar_abs": continue
        if hname not in hdict.keys():
            print(f"{hname} not found in file")
            continue
            #raise ValueError(f"{hname} not found in file")
        histo = hdict[hname][:, :, "nosys", :].integrate("sys") # process, sys, jmult, btag
        #set_trace()
        btagSFs_histo = btagSFs_dict[hname][:, :, "nosys", :].integrate("sys") # process, sys, jmult, btag
        
        if histo.dense_dim() == 1:
            xtitle, xrebinning, x_lims, withData = variables[hname]
            orig_xtitle = xtitle
            xaxis_name = histo.dense_axes()[0].name
                ## rebin x axis
            if isinstance(xrebinning, np.ndarray):
                new_xbins = hist.Bin(xaxis_name, xaxis_name, xrebinning)
            elif isinstance(xrebinning, float) or isinstance(xrebinning, int):
                new_xbins = xrebinning
            histo = histo.rebin(xaxis_name, new_xbins)
            btagSFs_histo = btagSFs_histo.rebin(xaxis_name, new_xbins)
    
        if histo.dense_dim() == 2:
            xtitle, ytitle, xrebinning, yrebinning, x_lims, y_lims, withData, plot_xlabels, plot_ylabels = variables[hname]

            if "DeepCSV_bDisc" not in hname:
                vlines = [(len(xrebinning)-1)*ybin for ybin in range(1, len(yrebinning)-1)] if histo.dense_dim() == 2 else None
                nbins = (len(xrebinning)-1)*(len(yrebinning)-1) if histo.dense_dim() == 2 else None
    
            xaxis_name = histo.dense_axes()[0].name
            yaxis_name = histo.dense_axes()[1].name
                ## rebin x axis
            if isinstance(xrebinning, np.ndarray):
                new_xbins = hist.Bin(xaxis_name, xaxis_name, xrebinning)
            elif isinstance(xrebinning, float) or isinstance(xrebinning, int):
                new_xbins = xrebinning
            histo = histo.rebin(xaxis_name, new_xbins)
            btagSFs_histo = btagSFs_histo.rebin(xaxis_name, new_xbins)
                ## rebin y axis
            if isinstance(yrebinning, np.ndarray):
                new_ybins = hist.Bin(yaxis_name, yaxis_name, yrebinning)
            elif isinstance(yrebinning, float) or isinstance(yrebinning, int):
                new_ybins = yrebinning
            histo = histo.rebin(yaxis_name, new_ybins)
            btagSFs_histo = btagSFs_histo.rebin(yaxis_name, new_ybins)

        #set_trace()   
        if hname == "Lep_iso":
            x_lims = (0., 0.15) if (args.lepton == "Muon") else (0., 0.1)
    
        ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
        #for jmult in ["3Jets"]:
        #for jmult in ["4PJets"]:
        for jmult in sorted(set([key[2] for key in histo.values().keys()])):
            if "disc" in hname:
                xtitle = orig_xtitle.replace("j", "3") if jmult == "3Jets" else orig_xtitle.replace("j", "4")
   
            shape_dir = os.path.join(outdir, args.lepton, jmult, "BKG_Est_orthog", "BTagSFS_Shapes_Comp")
            if not os.path.isdir(shape_dir): os.makedirs(shape_dir)

                # get sideband and signal region hists
            cen_sb_histo = histo[:, btag_reg_names_dict["Central"]["reg"], jmult].integrate("jmult").integrate("btag")
            up_sb_histo = histo[:, btag_reg_names_dict["Up"]["reg"], jmult].integrate("jmult").integrate("btag")
            dw_sb_histo = histo[:, btag_reg_names_dict["Down"]["reg"], jmult].integrate("jmult").integrate("btag")
                    # with btagSFs
            btagSFs_cen_sb_histo = btagSFs_histo[:, btag_reg_names_dict["Central"]["reg"], jmult].integrate("jmult").integrate("btag")
            btagSFs_up_sb_histo = btagSFs_histo[:, btag_reg_names_dict["Up"]["reg"], jmult].integrate("jmult").integrate("btag")
            btagSFs_dw_sb_histo = btagSFs_histo[:, btag_reg_names_dict["Down"]["reg"], jmult].integrate("jmult").integrate("btag")

            #set_trace()
            if histo.dense_dim() == 1:
                # make data-(single top/ttbar) shape plots
                compare_dmtop_shapes_btagSFs(
                    cen_sb_btagSF={"histo" : btagSFs_cen_sb_histo, "label" : "Cen: btag SFs", "color" : "k"},
                    cen_sb_nobtagSF={"histo" : cen_sb_histo, "label" : "Cen: no btag SFs", "color" : "g"},
                    up_sb_btagSF={"histo" : btagSFs_up_sb_histo, "label" : "Up: btag SFs", "color" : "r"},
                    up_sb_nobtagSF={"histo" : up_sb_histo, "label" : "Up: no btag SFs", "color" : "#984ea3"},
                    dw_sb_btagSF={"histo" : btagSFs_dw_sb_histo, "label" : "Down: btag SFs", "color" : "b"},
                    dw_sb_nobtagSF={"histo" : dw_sb_histo, "label" : "Down: no btag SFs", "color" : "#ff7f00"},
                   **{"xtitle":xtitle, "xlims":x_lims, "pltdir":shape_dir, "hname": hname, "jmult":jmult, 
                        "fname": "BTagSFs_Comp_BKG_DD_Shapes_%s_%s_%s" % (args.lepton, jmult, hname)})

            #set_trace()
            if histo.dense_dim() == 2:
                    # linearize sideband and signal region hists (hardcoded)
                cen_sb_lin = Plotter.linearize_hist(cen_sb_histo)
                up_sb_lin = Plotter.linearize_hist(up_sb_histo)
                dw_sb_lin = Plotter.linearize_hist(dw_sb_histo)
                    # with btagSFs
                btagSFs_cen_sb_lin = Plotter.linearize_hist(btagSFs_cen_sb_histo)
                btagSFs_up_sb_lin = Plotter.linearize_hist(btagSFs_up_sb_histo)
                btagSFs_dw_sb_lin = Plotter.linearize_hist(btagSFs_dw_sb_histo)

                # make data-(single top/ttbar) shape plots
                compare_dmtop_shapes_btagSFs(
                    #cen_sb={"histo" : cen_sb_lin, "label" : btag_reg_names_dict["Central"]["label"], "color" : btag_reg_names_dict["Central"]["color"]},
                    #up_sb={"histo" : up_sb_lin, "label" : btag_reg_names_dict["Up"]["label"], "color" : btag_reg_names_dict["Up"]["color"]},
                    #dw_sb={"histo" : dw_sb_lin, "label" : btag_reg_names_dict["Down"]["label"], "color" : btag_reg_names_dict["Down"]["color"]},
                    cen_sb_btagSF={"histo" : btagSFs_cen_sb_lin, "label" : "Cen: btag SFs", "color" : "k"},
                    cen_sb_nobtagSF={"histo" : cen_sb_lin, "label" : "Cen: no btag SFs", "color" : "g"},
                    up_sb_btagSF={"histo" : btagSFs_up_sb_lin, "label" : "Up: btag SFs", "color" : "r"},
                    up_sb_nobtagSF={"histo" : up_sb_lin, "label" : "Up: no btag SFs", "color" : "#984ea3"},
                    dw_sb_btagSF={"histo" : btagSFs_dw_sb_lin, "label" : "Down: btag SFs", "color" : "b"},
                    dw_sb_nobtagSF={"histo" : dw_sb_lin, "label" : "Down: no btag SFs", "color" : "#ff7f00"},
                   **{"xtitle": xtitle, "xlims":(0, nbins), "pltdir":shape_dir, "hname": hname, "jmult":jmult, "vlines":vlines, "fname": f"BTagSFs_Comp_BKG_DD_Shapes_{args.lepton}_{jmult}_{hname}",
                        "ytext":(ctstar_binlabels, ctstar_bin_locs) if hname == "mtt_vs_tlep_ctstar_abs" else None, "xtext":(mtt_bin_inds_to_plot, mtt_bins_to_plot) if hname == "mtt_vs_tlep_ctstar_abs" else None})
                        #"ytext":(ctstar_binlabels, tmp_ctstar_bin_locs) if hname == "mtt_vs_tlep_ctstar_abs" else None, "xtext":(tmp_mtt_bin_inds_to_plot, tmp_mtt_bins_to_plot) if hname == "mtt_vs_tlep_ctstar_abs" else None})

                ## make 1D projection along dense axes
                #for dax in range(2):
                #    if dax == 0:
                #        ## plot projection onto yaxis
                #        compare_dmtop_shapes_btagSFs(
                #            cen_sb={"histo" : cen_sb_histo.integrate(cen_sb_histo.dense_axes()[dax].name), "label" : btag_reg_names_dict["Central"]["label"], "color" : btag_reg_names_dict["Central"]["color"]},
                #            up_sb={"histo" : up_sb_histo.integrate(up_sb_histo.dense_axes()[dax].name), "label" : btag_reg_names_dict["Up"]["label"], "color" : btag_reg_names_dict["Up"]["color"]},
                #            dw_sb={"histo" : dw_sb_histo.integrate(dw_sb_histo.dense_axes()[dax].name), "label" : btag_reg_names_dict["Down"]["label"], "color" : btag_reg_names_dict["Down"]["color"]},
                #           **{"xtitle":ytitle, "xlims":y_lims, "pltdir":shape_dir, "jmult":jmult, "fname": "BTagSFs_Comp_BKG_DD_Shapes_%s_%s_%s_proj" % (args.lepton, jmult, hname.split("_vs_")[1])})
                #    else:
                #        ## plot projection onto xaxis
                #        compare_dmtop_shapes_btagSFs(
                #            cen_sb={"histo" : cen_sb_histo.integrate(cen_sb_histo.dense_axes()[dax].name), "label" : btag_reg_names_dict["Central"]["label"], "color" : btag_reg_names_dict["Central"]["color"]},
                #            up_sb={"histo" : up_sb_histo.integrate(up_sb_histo.dense_axes()[dax].name), "label" : btag_reg_names_dict["Up"]["label"], "color" : btag_reg_names_dict["Up"]["color"]},
                #            dw_sb={"histo" : dw_sb_histo.integrate(dw_sb_histo.dense_axes()[dax].name), "label" : btag_reg_names_dict["Down"]["label"], "color" : btag_reg_names_dict["Down"]["color"]},
                #           **{"xtitle":xtitle, "xlims":x_lims, "pltdir":shape_dir, "jmult":jmult, "fname": "BTagSFs_Comp_BKG_DD_Shapes_%s_%s_%s_proj" % (args.lepton, jmult, hname.split("_vs_")[0])})




toc = time.time()
print("Total time: %.1f" % (toc - tic))
