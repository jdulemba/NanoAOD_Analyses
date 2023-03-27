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
parser.add_argument("year", choices=["2018"], help="What year is the ntuple from.")
#parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("lepton", choices=["Electron", "Muon"], help="Choose which lepton to make plots for")
parser.add_argument("--nosys", action="store_true", help="Make plots without systematics and no qcd estimation")
parser.add_argument("--alpha_comp", action="store_true", help="Make plots comparing MC data with and without alpha correction applied.")
parser.add_argument("--group_tt", action="store_true", help="Group all ttbar events together")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
analyzer = "HBSR_NoAlphaCorr"
eos_dir = os.environ["eos_dir"]
plots_dir = os.environ["plots_dir"]

#set_trace()
base_ext = "*BKG*TOT.coffea"
input_dir = os.path.join(eos_dir, "results", f"{args.year}_{jobid}", analyzer)
fnames = fnmatch.filter(os.listdir(input_dir), base_ext)
fnames = [os.path.join(input_dir, fname) for fname in fnames]
hdict = load(fnames[0])
#set_trace()

if args.alpha_comp:
    #set_trace()
    fnames = fnmatch.filter(os.listdir(os.path.join(eos_dir, "results", f"{args.year}_{jobid}", "htt_btag_sb_regions")), base_ext)
    fnames = [os.path.join(os.path.join(eos_dir, "results", f"{args.year}_{jobid}", "htt_btag_sb_regions"), fname) for fname in fnames]
    alpha_dict = load(fnames[0])

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
    #"mtt_vs_tlep_ctstar_abs" : (cfeatures.variable_names_to_labels["mtt"], cfeatures.variable_names_to_labels["tlep_ctstar_abs"], linearize_binning[0], linearize_binning[1],
    #    (linearize_binning[0][0], linearize_binning[0][-1]), (linearize_binning[1][0], linearize_binning[1][-1]), False, None, (ctstar_binlabels, ctstar_bin_locs)),
    #    #(linearize_binning[0][0], linearize_binning[0][-1]), (linearize_binning[1][0], linearize_binning[1][-1]), True, None, (ctstar_binlabels, ctstar_bin_locs)),
    #"Jets_phi_vs_eta" : (cfeatures.variable_names_to_labels["Jets_phi"], cfeatures.variable_names_to_labels["Jets_eta"],
    #    phi_eta_binning[0], phi_eta_binning[1], (-3.3, 3.3), (-2.6, 2.6), True, phi_binlabels, (eta_binlabels, eta_bin_locs)),
    #"Lep_phi_vs_eta" : (cfeatures.variable_names_to_labels["Lep_phi"] % cfeatures.objtypes[args.lepton], cfeatures.variable_names_to_labels["Lep_eta"] % cfeatures.objtypes[args.lepton],
    #    phi_eta_binning[0], phi_eta_binning[1], (-3.3, 3.3), (-2.6, 2.6), True, phi_binlabels, (eta_binlabels, eta_bin_locs)),
    #"Jets_njets" : (cfeatures.variable_names_to_labels["Jets_njets"], 1, (0, 15), True),
    #"mtt" : (cfeatures.variable_names_to_labels["mtt"], 4, (200., 2000.), False),
    ##"mtt" : (cfeatures.variable_names_to_labels["mtt"], 4, (200., 2000.), True),
    #"tlep_ctstar_abs" : (cfeatures.variable_names_to_labels["tlep_ctstar_abs"], 1, (0., 1.), False),
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
    #"tlep_ctstar" : (cfeatures.variable_names_to_labels["tlep_ctstar"], 2, (-1., 1.), False),
    ##"tlep_ctstar" : (cfeatures.variable_names_to_labels["tlep_ctstar"], 2, (-1., 1.), True),
    #"full_disc" : (cfeatures.variable_names_to_labels["full_disc"], 2, (5, 25.), True),
    #"mass_disc" : (cfeatures.variable_names_to_labels["mass_disc"], 2, (0, 20.), True),
    #"ns_disc" : (cfeatures.variable_names_to_labels["ns_disc"], 2, (3., 10.), True),
    #"ns_dist" : (cfeatures.variable_names_to_labels["ns_dist"], 1, (0., 150.), True),
    #"Jets_pt" : (cfeatures.variable_names_to_labels["Jets_pt"], 1, (0., 300.), True),
    #"Jets_eta" : (cfeatures.variable_names_to_labels["Jets_eta"], 2, (-2.6, 2.6), True),
    #"Jets_phi" : (cfeatures.variable_names_to_labels["Jets_phi"], 2, (-4., 4.), True),
    #"Jets_LeadJet_pt" : (cfeatures.variable_names_to_labels["Jets_LeadJet_pt"], 1, (0., 300.), True),
    #"Jets_LeadJet_eta" : (cfeatures.variable_names_to_labels["Jets_LeadJet_eta"], 2, (-2.6, 2.6), True),
    #"Jets_DeepCSV_bDisc" : (cfeatures.variable_names_to_labels["Jets_DeepCSV_bDisc"], 1, (-0.01, 1.), True),
    #"Jets_DeepJet_bDisc" : (cfeatures.variable_names_to_labels["Jets_DeepJet_bDisc"], 1, (-0.01, 1.), True),
    #"Lep_pt" : (cfeatures.variable_names_to_labels["Lep_pt"] % cfeatures.objtypes[args.lepton], 1, (0., 300.), True),
    #"Lep_eta" : (cfeatures.variable_names_to_labels["Lep_eta"] % cfeatures.objtypes[args.lepton], 2, (-2.6, 2.6), True),
    #"Lep_phi" : (cfeatures.variable_names_to_labels["Lep_phi"] % cfeatures.objtypes[args.lepton], 2, (-4, 4), True),
    #"Lep_iso" : (cfeatures.variable_names_to_labels["Lep_iso"] % cfeatures.objtypes[args.lepton], 1, (0., 1.), True),
    #"MT" : (cfeatures.variable_names_to_labels["MT"], 1, (0., 300.), True),
    #"MET_pt" : (cfeatures.variable_names_to_labels["MET_pt"], 1, (0., 300.), True),
    #"MET_phi" : (cfeatures.variable_names_to_labels["MET_phi"], 1, (-3.2, 3.2), True),
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

    if args.alpha_comp:
        alpha_dict[hname].scale(lumi_correction, axis="dataset") # scale hists
        alpha_dict[hname] = alpha_dict[hname].group(btag_cat, btag_bin, btag_groups) # group by btagging sideband regions
        alpha_dict[hname] = alpha_dict[hname].group(process_cat, process, process_groups) # group by process
        alpha_dict[hname] = alpha_dict[hname][:, :, :, :, args.lepton].integrate("leptype") # only pick out specified lepton


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
                for btagregion in sorted(set([key[1] for key in histo.values().keys()])):
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
                for btagregion in sorted(set([key[1] for key in histo.values().keys()])):
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
                        #"maskData" : not withData,
                        "maskData" : not withData if btagregion == "btagPass" else False,
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


    ## make plots
if args.alpha_comp:
    for hname in variables.keys():
        if hname not in hdict.keys():
            print(f"{hname} not found in file")
            continue
            #raise ValueError(f"{hname} not found in file")

        #set_trace()
        histo = hdict[hname][:, :, "nosys", :].integrate("sys") # process, sys, jmult, btag
        alpha_histo = alpha_dict[hname][:, :, "nosys", :].integrate("sys") # process, sys, jmult, btag
        #set_trace()
    
        if histo.dense_dim() == 1:
            xtitle, rebinning, x_lims, withData = variables[hname]
            orig_xtitle = xtitle
            if rebinning != 1:
                xaxis_name = histo.dense_axes()[0].name
                histo = histo.rebin(xaxis_name, rebinning)
                alpha_histo = alpha_histo.rebin(xaxis_name, rebinning)
    
            ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
            for jmult in sorted(set([key[2] for key in histo.values().keys()])):
                #for btagregion in ["btagPass"]:
                for btagregion in sorted(set([key[1] for key in histo.values().keys()])):
                    print(*[jmult, btagregion, hname], sep=", ")
                    pltdir = os.path.join(outdir, args.lepton, jmult, f"Alpha_Comp_{btagregion}")
                    if not os.path.isdir(pltdir):
                        os.makedirs(pltdir)
   
                    fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                    fig.subplots_adjust(hspace=.07)
    
                    hslice = histo[Plotter.nonsignal_samples, btagregion, jmult].integrate("jmult").integrate("btag")
                    alpha_hslice = alpha_histo[Plotter.nonsignal_samples, btagregion, jmult].integrate("jmult").integrate("btag")

                    if "disc" in hname:
                        xtitle = orig_xtitle.replace("j", "3") if jmult == "3Jets" else orig_xtitle.replace("j", "4")
   
                    if hname == "Lep_iso":
                        x_lims = (0., 0.15) if (args.lepton == "Muon") else (0., 0.1)
    
                    #set_trace()

                    # plot MC
                    hep.plot.histplot(hslice[Plotter.mc_samples].integrate("process").values()[()], hslice.axis(xaxis_name).edges(), ax=ax, histtype="step", **{"linestyle" : "-", "color" : "r", "label" : "Sim.: No $\\alpha$"})
                    hep.plot.histplot(hslice[Plotter.data_samples].integrate("process").values()[()], hslice.axis(xaxis_name).edges(), ax=ax, histtype="step", **{"linestyle" : "-", "color" : "k", "label" : "Data: No $\\alpha$"})
                    hep.plot.histplot(alpha_hslice[Plotter.mc_samples].integrate("process").values()[()], hslice.axis(xaxis_name).edges(), ax=ax, histtype="step", **{"linestyle" : "-", "color" : "g", "label" : "Sim.: $\\alpha$"})
                    hep.plot.histplot(alpha_hslice[Plotter.data_samples].integrate("process").values()[()], hslice.axis(xaxis_name).edges(), ax=ax, histtype="step", **{"linestyle" : "-", "color" : "b", "label" : "Data: $\\alpha$"})
   
                    noAlpha_masked_vals, noAlpha_masked_bins = Plotter.get_ratio_arrays(num_vals=hslice[Plotter.data_samples].integrate("process").values()[()],
                        denom_vals=hslice[Plotter.mc_samples].integrate("process").values()[()], input_bins=hslice.axis(xaxis_name).edges())
                    Alpha_masked_vals, Alpha_masked_bins = Plotter.get_ratio_arrays(num_vals=alpha_hslice[Plotter.data_samples].integrate("process").values()[()],
                        denom_vals=alpha_hslice[Plotter.mc_samples].integrate("process").values()[()], input_bins=alpha_hslice.axis(xaxis_name).edges())

                    rax.step(noAlpha_masked_bins, noAlpha_masked_vals, where="post", **{"linestyle" : "-", "color" : "k", "label" : "No $\\alpha$"})
                    rax.step(Alpha_masked_bins, Alpha_masked_vals, where="post", **{"linestyle" : "-", "color" : "b", "label" : "$\\alpha$"})

                    ax.legend(ncol=2)
                    ax.autoscale()
                    ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.3)
                    ax.set_xlim(x_lims)
                    ax.set_xlabel(None)
                    ax.set_ylabel("Events")

                    rax.legend(ncol=2)
                    rax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
                    rax.autoscale()
                    rax.set_xlim(x_lims)
                    rax.set_xlabel(xtitle)
                    rax.set_ylabel("$\dfrac{data}{Pred.}$")
                    rax.set_ylim(0.5, 1.5) 
                    #set_trace() 
    
                        # add lepton/jet multiplicity label
                    #set_trace()
                    ax.text(
                        0.02, 0.88, f"{cfeatures.channel_labels[f'{args.lepton}_{jmult}']}",
                        #0.02, 0.84, f"{cfeatures.channel_labels[f'{args.lepton}_{jmult}']}\n{btag_cats[btagregion]}",
                        fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                    )
                    hep.cms.label(ax=ax, data=withData, label="Preliminary", year=year_to_use, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))
    
                    #set_trace()
                    figname = os.path.join(pltdir, "_".join([jmult, args.lepton, "Alpha_Comp", btagregion, hname]))
                    #figname = os.path.join(pltdir, "_".join([jmult, args.lepton, btagregion, hname, "CombinedTT"])) if args.group_tt else os.path.join(pltdir, "_".join([jmult, args.lepton, btagregion, hname]))
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
            alpha_histo = alpha_histo.rebin(yaxis_name, new_ybins).rebin(xaxis_name, new_xbins)
   
            ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
            #for jmult in ["3Jets"]:
            ##for jmult in ["4PJets"]:
            #   for btagregion in ["btagPass"]:
            for jmult in sorted(set([key[2] for key in histo.values().keys()])):
                #for btagregion in ["btagPass"]:
                for btagregion in sorted(set([key[1] for key in histo.values().keys()])):
                    print(*[jmult, btagregion, hname], sep=", ")
                    pltdir = os.path.join(outdir, args.lepton, jmult, f"Alpha_Comp_{btagregion}")
                    if not os.path.isdir(pltdir):
                        os.makedirs(pltdir)
  
                    if (btagregion == "btagPass") and (hname == "mtt_vs_tlep_ctstar_abs"): 
                        withData = False
                    else:
                        withData = variables[hname][-1]

                    hslice = histo[Plotter.nonsignal_samples, btagregion, jmult].integrate("jmult").integrate("btag")
                    alpha_hslice = alpha_histo[Plotter.nonsignal_samples, btagregion, jmult].integrate("jmult").integrate("btag")
   
                    mc_opts = {
                        "mcorder" : ["ttJets", "singlet", "EWK", "QCD"] if args.group_tt else ["ttJets_right", "ttJets_matchable", "ttJets_unmatchable", "ttJets_sl_tau", "ttJets_other", "singlet", "EWK", "QCD"],
                        #"maskData" : not withData,
                        "maskData" : not withData if btagregion == "btagPass" else False,
                    }

                        # plot linearized view 
                    fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True, figsize=(15.0, 10.0))
                    fig.subplots_adjust(hspace=.07)

                    hline = Plotter.linearize_hist(hslice)
                    alpha_hline = Plotter.linearize_hist(alpha_hslice)
                    
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




toc = time.time()
print("Total time: %.1f" % (toc - tic))