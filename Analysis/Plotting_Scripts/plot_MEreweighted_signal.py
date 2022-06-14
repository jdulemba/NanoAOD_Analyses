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
rcParams["savefig.format"] = "png"
rcParams["savefig.bbox"] = "tight"

from coffea.util import load
from pdb import set_trace
import os
import Utilities.prettyjson as prettyjson
import numpy as np
import Utilities.Plotter as Plotter
import Utilities.final_analysis_binning as final_binning
import uproot

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("--nosys", action="store_true", help="Make plots for nominal templates.")
parser.add_argument("--sys", action="store_true", help="Make plots for systematic variation templates.")
parser.add_argument("--allowed_masses", type=str, default="All", help="Choose which mass points to make plots for, multiple options can be input as ':' separated strings.")
parser.add_argument("--allowed_widths", type=str, default="All", help="Choose which width points to make plots for, multiple options can be input as ':' separated strings.")
args = parser.parse_args()


def plot_comparison(coffea_histo, uproot_histo, sig_type, signal_label):
        # plot all signal components together
    fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True, figsize=(15.0, 10.0))
    fig.subplots_adjust(hspace=.07)
    
        # plot original yields
    ch_sumw, ch_sumw2 = coffea_histo.values(sumw2=True)[()]
    #set_trace()
    ch_xerrs = np.array([coffea_histo.dense_axes()[0].centers()[idx] - coffea_histo.dense_axes()[0].edges()[idx] for idx in range(len(coffea_histo.dense_axes()[0].centers()))])
    ax.errorbar(coffea_histo.dense_axes()[0].centers(), ch_sumw, xerr=ch_xerrs, yerr=np.sqrt(ch_sumw2), fmt="none", **{"color" : "r", "label" : "ME reweighted"}) # ME reweighted template
    hep.plot.histplot(uproot_histo.values(), coffea_histo.dense_axes()[0].edges(), ax=ax, histtype="step", **{"color" : "b", "linestyle" : "-", "label" : "Morphed"}) # morphed template

    #if (sig_type == "Res") or (sig_type == "Int_pos"):
    #    #set_trace()
    #    if np.any(coffea_histo.values()[()] < 0):
    #        print(f"{sig_type} ME reweighted hist has {np.where(coffea_histo.values()[()] < 0)[0].size} negative bins")
    #    if np.any(uproot_histo.values() < 0):
    #        print(f"{sig_type} Morphed hist has {np.where(uproot_histo.values() < 0)[0].size} negative bins")

        # plot Morphed/MEreweight ratio
    ratio_masked_vals, ratio_masked_bins = Plotter.get_ratio_arrays(num_vals=uproot_histo.values(), denom_vals=coffea_histo.values()[()], input_bins=coffea_histo.dense_axes()[0].edges())
    rax.step(ratio_masked_bins, ratio_masked_vals, where="post", **{"color" : "k"})

    ax.legend(title=signal_label, loc="upper right", title_fontsize=24) # only show title, not line
    ax.autoscale()
    ax.set_ylabel("Events")
    ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.3)
    ax.axhline(0, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
    ax.set_xlabel(None)
    ax.set_xlim(x_lims)

    rax.autoscale()
    rax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
    rax.set_ylabel("Morphed/ME")
    #rax.set_ylabel("Morphed/ME reweighted")
    rax.set_ylim(max(0.0, rax.get_ylim()[0]), min(1.5, rax.get_ylim()[1]))
    #rax.set_ylim(max(0.0, rax.get_ylim()[0]), min(5.0, rax.get_ylim()[1]))
    rax.set_xlim(x_lims)
    rax.set_xlabel("$m_{t\\bar{t}}$ [GeV]")

        # add lepton/jet multiplicity label
    ax.text(
        0.02, 0.88, f"{objtypes['Lep'][lepton]}, {jet_mults[jmult]}",
        fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
    )
    
        ## draw vertical lines for distinguishing different ctstar bins
    vlines = [x_lims[1]*ybin/5 for ybin in range(1, 5)]
    [ax.axvline(vline, color="k", linestyle="--") for vline in vlines]
    [rax.axvline(vline, color="k", linestyle="--") for vline in vlines]
    
    for idx, label in enumerate(ctstar_binlabels):
        ax.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
            xytext=(0, 250), textcoords="offset points", va="top", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)
            #xytext=(0, 450), textcoords="offset points", va="top", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)
    
    #set_trace()
    rax.set_xticks(mtt_bin_inds_to_plot)
    rax.set_xticklabels(mtt_bins_to_plot)
    
    hep.cms.label(ax=ax, data=False, year=year_to_use, lumi=round(data_lumi_year[f"{lepton}s"]/1000., 1))
    
    #set_trace()
    sig_figname = os.path.join(signal_dir, "_".join([f"{sigpoint}_{sig_type}", jmult, lepton, "MEreweight_Morphed_Comparison"]))
    fig.savefig(sig_figname)
    print(f"{sig_figname} written")
    plt.close(fig)
        

if __name__ == "__main__":
    proj_dir = os.environ["PROJECT_DIR"]
    base_jobid = os.environ["base_jobid"]
    jobid = os.environ["jobid"]
    analyzer = "htt_btag_sb_regions"
    
    #set_trace()
    if args.nosys:
        input_dir = os.path.join(proj_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}")
        outdir = os.path.join(proj_dir, "plots", f"{args.year}_{jobid}", analyzer)
        #set_trace()
        fname = os.path.join(input_dir, f"raw_templates_lj_MEsig_{args.year}_{jobid}.coffea")
        if not os.path.isfile(fname): raise ValueError(f"{fname} not found")
        hdict = load(fname)
    
        nosys_morphed_rfiles_dict = {fname.strip(".root") : uproot.open(os.path.join(f"root://eosuser.cern.ch//eos/user/j/jdulemba/NanoAOD_Analyses/results/Morphed_Signal", fname)) for fname in os.listdir("/eos/user/j/jdulemba/NanoAOD_Analyses/results/Morphed_Signal")}
    if args.sys:
        outdir = os.path.join(proj_dir, "plots", f"{args.year}_{jobid}", f"Templates_{analyzer}", "MEreweighted_morphed_Comp_SIG")
        final_mew_rfile = uproot.open("root://cmseos.fnal.gov//store/user/lpcbtagging/UR_ntuples/heavyhiggsinputs/v13_mew_smoothed/templates_lj_sig.root")
        original_rfile = uproot.open(f"{proj_dir}/results/{args.year}_{jobid}/Templates_{analyzer}/FINAL/final_templates_lj_MEsig_{args.year}_{jobid}.root")
        #set_trace()
        morphed_signal_rfiles_dict = {fname.strip(".root") : uproot.open(os.path.join("/eos/user/j/jdulemba/NanoAOD_Analyses/results/Morphed_Signal/Uncs/Smoothed", fname)) for fname in os.listdir("/eos/user/j/jdulemba/NanoAOD_Analyses/results/Morphed_Signal/Uncs/Smoothed")}
        #morphed_signal_rfiles_dict = {fname.strip(".root") : uproot.open(os.path.join("/eos/user/j/jdulemba/NanoAOD_Analyses/results/Morphed_Signal/Uncs", fname)) for fname in os.listdir("/eos/user/j/jdulemba/NanoAOD_Analyses/results/Morphed_Signal/Uncs")}
    
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
        
    jet_mults = {
        "3Jets" : "3 jets",
        "4PJets" : "4+ jets"
    }
    
    objtypes = {
        "Jets" : "jets",
        "Lep" :  {
            "Muon" : "$\\mu$",
            "Electron" : "$e$",
        }
    }
    
    if args.year == "2016APV":
        year_to_use = "2016pre"
    elif args.year == "2016":
        year_to_use = "2016post"
    else:
        year_to_use = args.year
    
    linearize_binning = (
        final_binning.mtt_binning,
        final_binning.ctstar_abs_binning
    )
    
    ctstar_binlabels = [r"|$\cos (\theta^{*}_{t_{l}})$| $\in [%s, %s)$" % (linearize_binning[1][bin], linearize_binning[1][bin+1]) for bin in range(len(linearize_binning[1])-1)]
    ctstar_binlabels[-1] = r"|$\cos (\theta^{*}_{t_{l}})$| $\in [%s, %s]$" % (linearize_binning[1][-2], linearize_binning[1][-1])
    ctstar_bin_locs = np.linspace((len(linearize_binning[0])-1)/2, (len(linearize_binning[0])-1)*(len(linearize_binning[1])-1) - (len(linearize_binning[0])-1)/2, len(linearize_binning[1])-1)
    #mtt_bins_to_plot = np.array([400., 500., 600., 800., 1200.])
    #set_trace()
    mtt_vals_to_plot = np.array([400, 600, 1000])
    mtt_tiled_labels = np.tile(linearize_binning[0][:-1], linearize_binning[1].size-1)
    mtt_bin_inds_to_plot = np.where(np.in1d(mtt_tiled_labels, mtt_vals_to_plot))[0]
    mtt_bins_to_plot = np.tile(mtt_vals_to_plot, linearize_binning[1].size-1)
    
    
        ## get data lumi and scale MC by lumi
    data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())[args.year]
    kfactors_corr = prettyjson.loads(open(os.path.join(proj_dir, "inputs", "signal_kfactors_ulkfactor_final_220129.json")).read())
    
    #possible_masses = [365, 380, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000]
    possible_masses = [400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000]
    #possible_masses = [400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000]
    possible_widths = [2.5, 10.0, 15.0]
    #possible_widths = [2.5, 3.0, 4.0, 5.0, 8.0, 10.0, 13.0, 15.0, 18.0, 21.0, 25.0]
    #possible_widths = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 8.0, 10.0, 13.0, 15.0, 18.0, 21.0, 25.0]
    possible_masses = [str(mass) for mass in possible_masses]
    possible_widths = [str(width) for width in possible_widths]

    masses_to_run = args.allowed_masses
    widths_to_run = args.allowed_widths

    if masses_to_run == "All":
        allowed_masses = possible_masses
    else:
        allowed_masses = masses_to_run.split(":")
        allowed_masses = [mass for mass in allowed_masses if mass in possible_masses]

    if widths_to_run == "All":
        allowed_widths = possible_widths
    else:
        allowed_widths = widths_to_run.split(":")
        allowed_widths = [width for width in allowed_widths if width in possible_widths]

    if not args.sys:
        allowed_masses = ["M"+mass for mass in allowed_masses]
        allowed_widths = ["W"+width.replace(".", "p") for width in allowed_widths]
        print(f"Making target signal points for {allowed_masses} and {allowed_widths}")

    #set_trace()
    
    res_styles = {"color" : "r", "linestyle" : "-", "label" :"Res"}
    int_neg_styles = {"color" : "b", "linestyle" : "-", "label" : "Int, $w < 0$"}
    int_pos_styles = {"color" : "g", "linestyle" : "-", "label" : "Int, $w > 0$"}
    combined_styles = {"color" : "k", "linestyle" : "-", "label" : "Combined"}


    #set_trace()
    if args.nosys:
        #for jmult in ["4PJets"]:
        #    for lepton in ["Muon"]:
        for jmult in sorted(hdict.keys()):
            for lepton in sorted(hdict[jmult].keys()):
                nosys_dists = sorted([key for key in hdict[jmult][lepton].keys() for mass in allowed_masses for width in allowed_widths if ( ("nosys" in key) and (mass in key) and (width in key) )])
        
                for sigpoint in sorted(set(["_".join(key.split("_")[:3]) for key in nosys_dists])): # get boson_mass_width
                    proc, mass, width = sigpoint.split("_")
                    if proc.startswith("H"): continue
                    print(f"\t{args.year}, {jmult} {lepton} {sigpoint}")
                    signal_dir = os.path.join(outdir, lepton, jmult, "MEreweighted_SIG", sigpoint)
                    if not os.path.isdir(signal_dir):
                        os.makedirs(signal_dir)
        
                    res_kfactor = kfactors_corr[f"AtoTTJetsSL_{mass}_{width}_Res"] # SL and DiLep kfactors are the same
                    int_kfactor = kfactors_corr[f"AtoTTJetsSL_{mass}_{width}_Int"] # SL and DiLep kfactors are the same for both parts of interference
                        # get ME reweighted signal
                    me_res = hdict[jmult][lepton][f"{sigpoint}_Res_nosys"].copy()
                    me_neg = hdict[jmult][lepton][f"{sigpoint}_Int_neg_nosys"].copy()
                    me_pos = hdict[jmult][lepton][f"{sigpoint}_Int_pos_nosys"].copy()
                    #set_trace()
                        # scale components by kfactor
                    me_res.scale(res_kfactor)
                    me_neg.scale(int_kfactor)
                    me_pos.scale(int_kfactor)

                        # combine all signal into one    
                    me_total = me_res.copy()
                    me_total.add(me_neg.copy())
                    me_total.add(me_pos.copy())
        
        
                    x_lims = (0, me_res.dense_axes()[0].centers().size)
        
                    sig_label = "$%s_{%s\ GeV}^{%s\%%}$" % (proc[0], mass.strip("M"), width.strip("W").replace("p", "."))
        
                        # get morphed signal
                    morphed_rfile = nosys_morphed_rfiles_dict[f"signal_{mass}_{width}"]
                    morphed_res = morphed_rfile[f"{lepton[0:2].lower()}{jmult.lower()}_{year_to_use}/res" if lepton == "Muon" else f"{lepton[0].lower()}{jmult.lower()}_{year_to_use}/res"]
                    morphed_neg = morphed_rfile[f"{lepton[0:2].lower()}{jmult.lower()}_{year_to_use}/neg" if lepton == "Muon" else f"{lepton[0].lower()}{jmult.lower()}_{year_to_use}/neg"]
                    morphed_pos = morphed_rfile[f"{lepton[0:2].lower()}{jmult.lower()}_{year_to_use}/pos" if lepton == "Muon" else f"{lepton[0].lower()}{jmult.lower()}_{year_to_use}/pos"]
        
                    morphed_total = morphed_res.to_boost().copy()
                    morphed_total = morphed_total.__add__(morphed_neg.to_boost().copy())
                    morphed_total = morphed_total.__add__(morphed_pos.to_boost().copy())
        
                    #set_trace()
                        # compare resonant signal
                    plot_comparison(coffea_histo=me_res, uproot_histo=morphed_res, sig_type="Res", signal_label=f"{sig_label}, Res")
                        # compare positive interference signal
                    plot_comparison(coffea_histo=me_pos, uproot_histo=morphed_pos, sig_type="Int_pos", signal_label=f"{sig_label}, $w > 0$")
                        # compare negative interference signal
                    plot_comparison(coffea_histo=me_neg, uproot_histo=morphed_neg, sig_type="Int_neg", signal_label=f"{sig_label}, $w < 0$")
                        # compare combined signal
                    plot_comparison(coffea_histo=me_total, uproot_histo=morphed_total, sig_type="Total", signal_label=f"{sig_label}, Total")


    if args.sys:
        import Utilities.systematics as systematics
            # only look at M=775 and W=2.5, 10.0 for systematic variations
        allowed_masses = ["m775"]
        allowed_widths = ["w10p0"]
        #allowed_widths = ["w2p5", "w10p0"]
        print(f"Making target signal points for {allowed_masses} and {allowed_widths}")
        #set_trace()
        
        for jmult in sorted(["3Jets", "4PJets"]):
            for lepton in (["Electron", "Muon"]):
        #for jmult in sorted(["3Jets"]):
        #    for lepton in (["Muon"]):
                dirname = f"{lepton[0:2].lower()}{jmult.lower()}_{year_to_use}" if lepton == "Muon" else f"{lepton[0].lower()}{jmult.lower()}_{year_to_use}"
                signal_dists = sorted([key for key in original_rfile[dirname].keys() for mass in allowed_masses for width in allowed_widths if ( (mass in key) and (width in key) and ("Up" not in key) and ("Down" not in key) )])
        
                for sigpoint in sorted(set(["_".join(key.split("_")[:3]) for key in signal_dists])): # get boson_mass_width
                    proc, mass, width = sigpoint.split("_")
                    if proc == "H": continue
                    print(f"\t{args.year}, {jmult} {lepton} {sigpoint}")

                        # get original templates
                    orig_res_nosys = original_rfile[dirname][f"{sigpoint}_res"]
                    orig_neg_nosys = original_rfile[dirname][f"{sigpoint}_neg"]
                    orig_pos_nosys = original_rfile[dirname][f"{sigpoint}_pos"]

                    x_edges = orig_res_nosys.axes[0].edges()
        
                    sig_label = "$%s_{%s\ GeV}^{%s\%%}$" % (proc, mass.strip("m"), width.strip("w").replace("p", "."))
        
                    #set_trace()
                    if f"signal_{mass}_{width}" in morphed_signal_rfiles_dict.keys():
                        morphed_res_nosys = morphed_signal_rfiles_dict[f"signal_{mass}_{width}"][dirname]["res"]
                        morphed_neg_nosys = morphed_signal_rfiles_dict[f"signal_{mass}_{width}"][dirname]["neg"]
                        morphed_pos_nosys = morphed_signal_rfiles_dict[f"signal_{mass}_{width}"][dirname]["pos"]
                    else:
                        morphed_res_nosys = None
                        morphed_neg_nosys = None
                        morphed_pos_nosys = None
                    #set_trace()

                    up_systs = sorted([key.split(f"{sigpoint}_res_")[-1].split(";")[0] for key in original_rfile[dirname].keys() if (f"{sigpoint}_res" in key) and ("Up" in key)])
                    for up_sysname in up_systs:
                        #set_trace()

                        if ";" in up_sysname: set_trace()
                        print(f"\t\t{up_sysname.replace('Up', '')}")
                        signal_dir = os.path.join(outdir, lepton, jmult, sigpoint, up_sysname.replace("Up", ""))
                        if not os.path.isdir(signal_dir):
                            os.makedirs(signal_dir)
                        #if up_sysname == "CMS_eff_m_reco_totUp": set_trace()

                        dw_sysname = up_sysname.replace("Up", "Down")

                            # get original variations
                        orig_res_sys_up, orig_res_sys_dw = original_rfile[dirname][f"{sigpoint}_res_{up_sysname}"], original_rfile[dirname][f"{sigpoint}_res_{dw_sysname}"]
                        orig_neg_sys_up, orig_neg_sys_dw = original_rfile[dirname][f"{sigpoint}_neg_{up_sysname}"], original_rfile[dirname][f"{sigpoint}_neg_{dw_sysname}"]
                        orig_pos_sys_up, orig_pos_sys_dw = original_rfile[dirname][f"{sigpoint}_pos_{up_sysname}"], original_rfile[dirname][f"{sigpoint}_pos_{dw_sysname}"]
                        
                            # get final variations
                        try:
                            final_res_sys_up, final_res_sys_dw = final_mew_rfile[dirname][f"{sigpoint}_res_{up_sysname}"], final_mew_rfile[dirname][f"{sigpoint}_res_{dw_sysname}"]
                            final_neg_sys_up, final_neg_sys_dw = final_mew_rfile[dirname][f"{sigpoint}_neg_{up_sysname}"], final_mew_rfile[dirname][f"{sigpoint}_neg_{dw_sysname}"]
                            final_pos_sys_up, final_pos_sys_dw = final_mew_rfile[dirname][f"{sigpoint}_pos_{up_sysname}"], final_mew_rfile[dirname][f"{sigpoint}_pos_{dw_sysname}"]
                        except:
                            print(f"{up_sysname.replace('Up', '')} not found in smoothed templates file")
                            #set_trace()
                            final_res_sys_up, final_res_sys_dw = None, None
                            final_neg_sys_up, final_neg_sys_dw = None, None
                            final_pos_sys_up, final_pos_sys_dw = None, None
                            #continue

                        if f"signal_{mass}_{width}" in morphed_signal_rfiles_dict.keys():
                            #set_trace()
                            try:
                                morphed_res_sys_up, morphed_res_sys_dw = morphed_signal_rfiles_dict[f"signal_{mass}_{width}"][dirname][f"res_{up_sysname}"], morphed_signal_rfiles_dict[f"signal_{mass}_{width}"][dirname][f"res_{dw_sysname}"]
                                morphed_neg_sys_up, morphed_neg_sys_dw = morphed_signal_rfiles_dict[f"signal_{mass}_{width}"][dirname][f"neg_{up_sysname}"], morphed_signal_rfiles_dict[f"signal_{mass}_{width}"][dirname][f"neg_{dw_sysname}"]
                                morphed_pos_sys_up, morphed_pos_sys_dw = morphed_signal_rfiles_dict[f"signal_{mass}_{width}"][dirname][f"pos_{up_sysname}"], morphed_signal_rfiles_dict[f"signal_{mass}_{width}"][dirname][f"pos_{dw_sysname}"]
                            except:
                                morphed_res_sys_up, morphed_res_sys_dw = None, None
                                morphed_neg_sys_up, morphed_neg_sys_dw = None, None
                                morphed_pos_sys_up, morphed_pos_sys_dw = None, None
                        else:
                            morphed_res_sys_up, morphed_res_sys_dw = None, None
                            morphed_neg_sys_up, morphed_neg_sys_dw = None, None
                            morphed_pos_sys_up, morphed_pos_sys_dw = None, None
                        
                        histos_dict = {
                            "Res" : {
                                "Up": [
                                    (orig_res_sys_up, orig_res_nosys, {"color": "r", "linestyle": "-", "label": "Original ME Up"}, True),
                                    (final_res_sys_up, orig_res_nosys, {"color": "r", "linestyle": "-", "label": "Smoothed ME Up"}, False),
                                    (morphed_res_sys_up, morphed_res_nosys, {"color": "r", "linestyle": "--", "label": "Smoothed Morphed Up"}, False),
                                ],
                                "Down": [
                                    (orig_res_sys_dw, orig_res_nosys, {"color": "b", "linestyle": "-", "label": "Original ME Down"}, True),
                                    (final_res_sys_dw, orig_res_nosys, {"color": "b", "linestyle": "-", "label": "Smoothed ME Down"}, False),
                                    (morphed_res_sys_dw, morphed_res_nosys, {"color": "b", "linestyle": "--", "label": "Smoothed Morphed Down"}, False),
                                ],
                            },
                            "Int_neg" : {
                                "Up": [
                                    (orig_neg_sys_up, orig_neg_nosys, {"color": "r", "linestyle": "--", "label": "Original ME Up"}, True),
                                    (final_neg_sys_up, orig_neg_nosys, {"color": "r", "linestyle": "-", "label": "Smoothed ME Up"}, False),
                                    (morphed_neg_sys_up, morphed_neg_nosys, {"color": "r", "linestyle": "--", "label": "Smoothed Morphed Up"}, False),
                                ],
                                "Down": [
                                    (orig_neg_sys_dw, orig_neg_nosys, {"color": "b", "linestyle": "--", "label": "Original ME Down"}, True),
                                    (final_neg_sys_dw, orig_neg_nosys, {"color": "b", "linestyle": "-", "label": "Smoothed ME Down"}, False),
                                    (morphed_neg_sys_dw, morphed_neg_nosys, {"color": "b", "linestyle": "--", "label": "Smoothed Morphed Down"}, False),
                                ],
                            },
                            "Int_pos" : {
                                "Up": [
                                    (orig_pos_sys_up, orig_pos_nosys, {"color": "r", "linestyle": "--", "label": "Original ME Up"}, True),
                                    (final_pos_sys_up, orig_pos_nosys, {"color": "r", "linestyle": "-", "label": "Smoothed ME Up"}, False),
                                    (morphed_pos_sys_up, morphed_pos_nosys, {"color": "r", "linestyle": "--", "label": "Smoothed Morphed Up"}, False),
                                ],
                                "Down": [
                                    (orig_pos_sys_dw, orig_pos_nosys, {"color": "b", "linestyle": "--", "label": "Original ME Down"}, True),
                                    (final_pos_sys_dw, orig_pos_nosys, {"color": "b", "linestyle": "-", "label": "Smoothed ME Down"}, False),
                                    (morphed_pos_sys_dw, morphed_pos_nosys, {"color": "b", "linestyle": "--", "label": "Smoothed Morphed Down"}, False),
                                ],
                            },
                        }

                        for sig_type in histos_dict.keys():
                            if sig_type == "Res":
                                leg_title = f"{sig_label}, Res"
                            elif sig_type == "Int_neg":
                                leg_title = f"{sig_label}, $w < 0$"
                            elif sig_type == "Int_pos":
                                leg_title = f"{sig_label}, $w > 0$"
                            elif sig_type == "Total":
                                leg_title = f"{sig_label}, Total"

                            up_histos = histos_dict[sig_type]["Up"]
                            dw_histos = histos_dict[sig_type]["Down"]

    
                            fig, ax = plt.subplots(figsize=(15.0, 10.0))
                            fig.subplots_adjust(hspace=.07)
    
                                ## plot relative deviations
                            for up_histo, nominal, up_style, use_fill_between in up_histos:
                                    # there is at least one actual value
                                if up_histo:
                                    up_masked_vals, up_masked_bins = Plotter.get_ratio_arrays(num_vals=up_histo.values(), denom_vals=nominal.values(), input_bins=x_edges)
                                    ax.fill_between(up_masked_bins, up_masked_vals, y2=1., facecolor=up_style["color"], step="post", alpha=0.5) if use_fill_between else ax.step(up_masked_bins, up_masked_vals, where="post", **up_style)
    
                            for dw_histo, nominal, dw_style, use_fill_between in dw_histos:
                                    # there is at least one actual value
                                if dw_histo:
                                    dw_masked_vals, dw_masked_bins = Plotter.get_ratio_arrays(num_vals=dw_histo.values(), denom_vals=nominal.values(), input_bins=x_edges)
                                    ax.fill_between(dw_masked_bins, dw_masked_vals, y2=1., facecolor=dw_style["color"], step="post", alpha=0.5) if use_fill_between else ax.step(dw_masked_bins, dw_masked_vals, where="post", **dw_style)
    

                            ax.legend(loc="upper right", title=leg_title, ncol=int(((len(up_histos)-1) + (len(dw_histos)-1))/2))
                            ax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
                            ax.autoscale()
                            #if "CMS_eff_m_iso" in up_sysname: set_trace()
                            ##ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.3)
                            #ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.02)
                            ##ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.1)
                            ax.set_xlim(x_edges[0], x_edges[-1])
                            ax.set_xlabel("$m_{t\\bar{t}}$ [GeV]")
                            ax.set_ylabel("Ratio to Nominal")

                                # add lepton/jet multiplicity label
                            ax.text(
                                0.02, 0.88, f"{up_sysname.replace('Up', '')}\n{objtypes['Lep'][lepton]}, {jet_mults[jmult]}",
                                fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                            )
                                ## draw vertical lines for distinguishing different ctstar bins
                            vlines = [x_edges[-1]*ybin/5 for ybin in range(1, 5)]
                            [ax.axvline(vline, color="k", linestyle="--") for vline in vlines]

                            for idx, label in enumerate(ctstar_binlabels):
                                ax.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                                    xytext=(0, 25), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)
                                    #xytext=(0, 450), textcoords="offset points", va="top", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)

                            ax.set_xticks(mtt_bin_inds_to_plot)
                            ax.set_xticklabels(mtt_bins_to_plot)

                            hep.cms.label(ax=ax, data=False, year=args.year, lumi=round(data_lumi_year[f"{lepton}s"]/1000., 1))

                            #set_trace()
                            figname = os.path.join(signal_dir, "_".join([sigpoint, sig_type, jmult, lepton, up_sysname.replace("Up", ""), "SysTemplates_Comp"]))
                            #set_trace()
                            fig.savefig(figname)
                            print(f"{figname} written")
                            plt.close(fig)



    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
