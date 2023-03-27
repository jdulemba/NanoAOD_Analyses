import time
tic = time.time()

# matplotlib
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend("agg")
from matplotlib import rcParams
rcParams["font.size"] = 20
rcParams["savefig.format"] = "png"
rcParams["savefig.bbox"] = "tight"

from coffea.hist import plot
from coffea.util import load, save
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
import coffea.processor as processor

base_jobid = os.environ["base_jobid"]

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("lepton", choices=["Electron", "Muon"], help="Choose which lepton to make plots for")
parser.add_argument("level", choices=["Gen", "Reco"], help="Choose which level to make plots for.")
parser.add_argument("--comparison", action="store_true", help="Make plots comparing ME reweighting dists to original MC signal points.")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
analyzer = "signal_ME_evtReweighting"

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

btag_cats = btag_sidebands.btag_cats

lep_cats = {
    "Tight" : "tight %s" % objtypes["Lep"][args.lepton],
    "Loose" : "loose %s" % objtypes["Lep"][args.lepton],
}

linearize_binning = (
    final_binning.mtt_binning,
    final_binning.ctstar_abs_binning
)

#ctstar_binlabels = ["%s $\leq$ |cos($\\theta^{*}_{t_{l}}$)| $\leq$ %s" % (linearize_binning[1][bin], linearize_binning[1][bin+1]) for bin in range(len(linearize_binning[1])-1)]
#ctstar_bin_locs = np.linspace((len(linearize_binning[0])-1)/2, (len(linearize_binning[0])-1)*(len(linearize_binning[1])-1) - (len(linearize_binning[0])-1)/2, len(linearize_binning[1])-1)
#mtt_binlabels = ["%s $\leq$ m($t\\bar{t}$) $\leq$ %s" % (linearize_binning[0][bin], linearize_binning[0][bin+1]) for bin in range(len(linearize_binning[0])-1)]*len(ctstar_binlabels)
#
## binning and bin labels for phi x eta
#phi_eta_binning = (
#    np.array([-3.2, -2.4, -1.57, -0.87, 0., 0.87, 1.57, 2.4, 3.2]), # phi
#    np.array([-2.5, -1.3, -0.7, 0., 0.7, 1.3, 2.5]) # eta binning
#)
#eta_binlabels = ["%s $\leq$ $\\eta$ $\leq$ %s" % (phi_eta_binning[1][bin], phi_eta_binning[1][bin+1]) for bin in range(len(phi_eta_binning[1])-1)]
#eta_bin_locs = np.linspace((len(phi_eta_binning[0])-1)/2, (len(phi_eta_binning[0])-1)*(len(phi_eta_binning[1])-1) - (len(phi_eta_binning[0])-1)/2, len(phi_eta_binning[1])-1)
#phi_binlabels = ["%s $\leq$ $\\phi$ $\leq$ %s" % (phi_eta_binning[0][bin], phi_eta_binning[0][bin+1]) for bin in range(len(phi_eta_binning[0])-1)]*len(eta_binlabels)


reco_variables = {
    "mtt_vs_tlep_ctstar_abs" : ("$m_{t\\bar{t}}$", "|cos($\\theta^{*}_{t_{l}}$)|", linearize_binning[0], linearize_binning[1], (linearize_binning[0][0], linearize_binning[0][-1]),
        (linearize_binning[1][0], linearize_binning[1][-1])),
    #"Jets_phi_vs_eta" : ("$\\phi$(jets)", "$\\eta$(jets)", phi_eta_binning[0], phi_eta_binning[1], (-3.3, 3.3), (-2.6, 2.6), False, phi_binlabels, (eta_binlabels, eta_bin_locs)),
    #"Lep_phi_vs_eta" : ("$\\phi$(%s)" % objtypes["Lep"][args.lepton], "$\\eta$(%s)" % objtypes["Lep"][args.lepton], phi_eta_binning[0], phi_eta_binning[1],
    #    (-3.3, 3.3), (-2.6, 2.6), False, phi_binlabels, (eta_binlabels, eta_bin_locs)),
    #"Jets_njets" : ("$n_{jets}$", 1, (0, 15)),
    #"mtt" : ("$m_{t\\bar{t}}$ [GeV]", 4, (200., 2000.)),
    #"tlep_ctstar_abs" : ("|cos($\\theta^{*}_{t_{l}}$)|", 1, (0., 1.)),
    #"mthad" : ("$m_{t_{h}}$ [GeV]", 2, (0., 300.)),
    #"mWHad" : ("$m_{W_{h}}$ [GeV]", 2, (0., 300.)),
    #"mWLep" : ("$m_{W_{l}}$ [GeV]", 2, (0., 300.)),
    #"pt_thad" : ("$p_{T}$($t_{h}$) [GeV]", 2, (0., 500.)),
    #"pt_tlep" : ("$p_{T}$($t_{l}$) [GeV]", 2, (0., 500.)),
    #"pt_tt" : ("$p_{T}$($t\\bar{t}$) [GeV]", 2, (0., 500.)),
    #"eta_thad" : ("$\\eta$($t_{h}$)", 2, (-4., 4.)),
    #"eta_tlep" : ("$\\eta$($t_{l}$)", 2, (-4., 4.)),
    #"eta_tt" : ("$\\eta$($t\\bar{t}$)", 2, (-4., 4.)),
    #"tlep_ctstar" : ("cos($\\theta^{*}_{t_{l}}$)", 2, (-1., 1.)),
    #"full_disc" : ("$\\lambda_{j}$", 2, (5, 25.)),
    #"mass_disc" : ("$\\lambda_{Mj}$", 2, (0, 20.)),
    #"ns_disc" : ("$\\lambda_{NSj}$", 2, (3., 10.)),
    #"ns_dist" : ("$D_{\\nu, min}$", 1, (0., 150.)),
    #"Jets_pt" : ("$p_{T}$(jets) [GeV]", 1, (0., 300.)),
    #"Jets_eta" : ("$\\eta$(jets)", 2, (-2.6, 2.6)),
    #"Jets_phi" : ("$\\phi$(jets)", 2, (-4., 4.)),
    #"Jets_LeadJet_pt" : ("$p_{T}$(leading jet) [GeV]", 1, (0., 300.)),
    #"Jets_LeadJet_eta" : ("$\\eta$(leading jet)", 2, (-2.6, 2.6)),
    #"Jets_DeepCSV_bDisc" : ("DeepCSV b Disc", 1, (-0.01, 1.)),
    #"Lep_pt" : ("$p_{T}$(%s) [GeV]" % objtypes["Lep"][args.lepton], 1, (0., 300.)),
    #"Lep_eta" : ("$\\eta$(%s)" % objtypes["Lep"][args.lepton], 2, (-2.6, 2.6)),
    #"Lep_phi" : ("$\\phi$(%s)" % objtypes["Lep"][args.lepton], 2, (-4, 4)),
    #"Lep_iso" : ("pfRelIso, %s" % objtypes["Lep"][args.lepton], 1, (0., 1.)),
    #"MT" : ("$M_{T}$ [GeV]", 1, (0., 300.)),
    #"MET_pt" : ("$p_{T}$(MET) [GeV]", 1, (0., 300.)),
    #"MET_phi" : ("$\phi$(MET)", 1, (-3.2, 3.2)),
}
gen_variables = {
    "mtt_vs_top_ctstar_abs" : ("$m_{t\\bar{t}}$", "|cos($\\theta^{*}_{t}$)|", linearize_binning[0], linearize_binning[1], (linearize_binning[0][0], linearize_binning[0][-1]),
        (linearize_binning[1][0], linearize_binning[1][-1])),
    "mtt" : ("$m_{t\\bar{t}}$ [GeV]", 4, (200., 2000.)),
    "top_ctstar_abs" : ("|cos($\\theta^{*}_{t}$)|", 1, (0., 1.)),
    "top_ctstar" : ("cos($\\theta^{*}_{t}$)", 1, (-1., 1.)),
    "mtop" : ("$m_{t}$ [GeV]", 2, (150., 200.)),
    "pt_top" : ("$p_{T}^{t}$ [GeV]", 2, (0., 500.)),
    "pt_tt" : ("$p_{T}^{t\\bar{t}}$ [GeV]", 2, (0., 500.)),
    "eta_top" : ("$\\eta_{t}$", 2, (-4., 4.)),
    "eta_tt" : ("$\\eta_{t\\bar{t}}$", 2, (-4., 4.)),
}
variables_dict = reco_variables if args.level == "Reco" else gen_variables

    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())[args.year]

outdir = os.path.join(proj_dir, "plots", f"{args.year}_{jobid}", analyzer, f"{args.level}Level", "Comparison" if args.comparison else "MEreweighting")
if not os.path.isdir(outdir):
    os.makedirs(outdir)

input_dir = f"/eos/user/j/jdulemba/NanoAOD_Analyses/results/{args.year}_{jobid}/{analyzer}/{args.level}Level"
base_ext = "BATCH*TOT.coffea"
fnames = fnmatch.filter(os.listdir(input_dir), base_ext)
fnames = [os.path.join(input_dir, fname) for fname in fnames]
averaged_output_histos = plt_tools.add_coffea_files(fnames)

## make groups based on process
#set_trace()
names = [dataset for dataset in sorted(set([key[0] for key in averaged_output_histos[sorted(variables_dict.keys())[0]].values().keys()]))] # get dataset names in hists
process = hist.Cat("dataset", "dataset", sorting="placement")
process_cat = "dataset"
process_groups = plt_tools.make_dataset_groups(args.lepton, args.year, samples=names, sigdict="MC_combined" if args.comparison else "MEreweight_combined")

if args.comparison:
    #signal_coffea_regex = re.compile(r"((?:%s))" % "|".join(["AtoTTJetsSL_M800_W10p0_Res", "AtoTTJetsSL_M800_W10p0_Int*", "HtoTTJetsSL_M800_W10p0_Res", "HtoTTJetsSL_M800_W10p0_Int*",\
    #    "AtoTTJetsDiLep_M800_W10p0_Res", "AtoTTJetsDiLep_M800_W10p0_Int*", "HtoTTJetsDiLep_M800_W10p0_Res", "HtoTTJetsDiLep_M800_W10p0_Int*"]))
    if args.level == "Gen":
        lumi_correction = load(os.path.join(proj_dir, "Corrections", base_jobid, "MC_LumiWeights.coffea"))[args.year][f"{args.lepton}s"]
        mc_hnames = fnmatch.filter(os.listdir(os.path.join("/".join(input_dir.split("/")[:-1]), "MC")), "BATCH*TOT.coffea")
        mc_hnames = [os.path.join("/".join(input_dir.split("/")[:-1]), "MC", fname) for fname in mc_hnames]
        mc_hdict = plt_tools.add_coffea_files(mc_hnames)
        #set_trace()
        signal_coffea_regex = re.compile(r"((?:%s))" % "|".join(sorted([key[0] for key in mc_hdict[sorted(mc_hdict.keys())[0]].integrate("sys").values().keys()])))
            # filter out non wanted signal points
        for hname in mc_hdict.keys():
            if isinstance(mc_hdict[hname], processor.accumulator.defaultdict_accumulator): continue
            mc_hdict[hname].scale(lumi_correction, axis="dataset")
            mc_hdict[hname] = mc_hdict[hname][signal_coffea_regex]
            mc_hdict[hname] = mc_hdict[hname].group(process_cat, process, process_groups)

    if args.level == "Reco":
        #set_trace()
        signal_coffea_regex = re.compile(r"((?:%s))" % "|".join(sorted(plt_tools.create_sig_groups("MC_indiv")[args.year].keys())))
            # get signal templates (no kfactors applied)
        mc_hname = os.path.join(proj_dir, "results", f"{args.year}_{jobid}", "Templates_htt_btag_sb_regions", f"raw_templates_lj_sig_{args.year}_{jobid}.coffea")
        mc_hdict = load(mc_hname)

    # scale and group hists by process
for hname in averaged_output_histos.keys():
    if "cutflow" in hname: continue
    if isinstance(averaged_output_histos[hname], processor.accumulator.defaultdict_accumulator): continue

    if args.comparison:
        averaged_output_histos[hname] = averaged_output_histos[hname][signal_coffea_regex]
    #set_trace()
    averaged_output_histos[hname] = averaged_output_histos[hname].group(process_cat, process, process_groups)
    if args.level == "Reco":
        averaged_output_histos[hname] = averaged_output_histos[hname][:, :, :, args.lepton, "btagPass"].integrate("leptype").integrate("btag") # only pick out specified lepton (and only btag region)


    ## make plots
#set_trace()
for hname in variables_dict.keys():
    if hname not in averaged_output_histos.keys():
        raise ValueError(f"{hname} not found in file")
    histo = averaged_output_histos[hname][:, "nosys", :].integrate("sys") # dataset, sys, jmult
    #set_trace()

    if histo.dense_dim() == 1:
        xtitle, rebinning, x_lims = variables_dict[hname]
        orig_xtitle = xtitle
        if rebinning != 1:
            xaxis_name = histo.dense_axes()[0].name
            histo = histo.rebin(xaxis_name, rebinning)

    elif histo.dense_dim() == 2:
        xtitle, ytitle, xrebinning, yrebinning, x_lims, y_lims = variables_dict[hname]

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

    else: raise ValueError("3D histograms or more dimensions are not supported currently.")

    if args.level == "Gen":
        print(hname)

        #set_trace()
        all_possible_signal = sorted(set([key[0] for key in histo.values().keys()]))
            # combine pos/neg interference into Int
        signal_to_use = ["_".join(signal.split("_")[:-1]) for signal in all_possible_signal if "neg" in signal] + [signal for signal in all_possible_signal if "Res" in signal]
        for signal in sorted(signal_to_use):
            print(f"\t{signal}")
            signal_dir = os.path.join(outdir, signal)
            if not os.path.isdir(signal_dir):
                os.makedirs(signal_dir)

            sig_label = styles[signal]["name"]
            signal_mask = re.compile(r"((?:%s*))" % signal)

            sig_histo = histo[signal_mask] # get all components from signal (pos and neg for Int, just one for Res)
                # make linearized view 
            if histo.dense_dim() == 2: 
                #set_trace()
                proc_histo = sig_histo.group("dataset", hist.Cat("process", "Process", sorting="placement"), {key[0] : [key[0]] for key in sorted(sig_histo.values().keys())})
                hline = Plotter.linearize_hist(proc_histo)

                ## make plots comparing original MC signal to ME reweighting signal
            if args.comparison:
                #set_trace()
                fig_comp, ax_comp = plt.subplots()
                fig_comp.subplots_adjust(hspace=.07)

                mc_sig_histo = mc_hdict[hname][signal_mask, "nosys"].integrate("sys")
                if histo.dense_dim() == 2:
                    mc_sig_histo = mc_sig_histo.rebin(xaxis_name, new_xbins)
                    mc_sig_histo = mc_sig_histo.rebin(yaxis_name, new_ybins)
                    mc_proc_histo = mc_sig_histo.group("dataset", hist.Cat("process", "Process", sorting="placement"), {key[0] : [key[0]] for key in sorted(mc_sig_histo.values().keys())})
                    mc_hline = Plotter.linearize_hist(mc_proc_histo)

                    Plotter.plot_1D(mc_hline.integrate("process").values()[()], mc_hline.dense_axes()[0].edges(), xlimits=(0, mc_hline.dense_axes()[0].edges()[-1]), ax=ax_comp, histtype="errorbar", label="Official MC", color="#e41a1c") # red
                    Plotter.plot_1D(hline.integrate("process").values()[()], hline.dense_axes()[0].edges(), xlimits=(0, hline.dense_axes()[0].edges()[-1]), ax=ax_comp, histtype="errorbar", label="ME Reweighting", color="#377eb8") # blue
                else:
                    mc_sig_histo = mc_sig_histo.rebin(mc_sig_histo.dense_axes()[0].name, rebinning)

                    Plotter.plot_1D(mc_sig_histo.integrate(mc_sig_histo.sparse_axes()[0].name).values()[()], mc_sig_histo.dense_axes()[0].edges(), xlimits=x_lims, ax=ax_comp, histtype="errorbar", label="Official MC", color="#e41a1c") # red
                    Plotter.plot_1D(sig_histo.integrate(sig_histo.sparse_axes()[0].name).values()[()], sig_histo.dense_axes()[0].edges(), xlimits=x_lims, ax=ax_comp, histtype="errorbar", label="ME Reweighting", color="#377eb8") # blue

                    # format axes/legend
                handles, labels = ax_comp.get_legend_handles_labels()
                #for idx, sample in enumerate(labels):
                #    if "Error" in sample:
                #        labels[idx] = None
                ax_comp.legend(handles, labels, title=sig_label, loc="upper right", title_fontsize=24)
                ax_comp.autoscale()
                ax_comp.set_ylim(ax_comp.get_ylim()[0], ax_comp.get_ylim()[1]*1.15)
                ax_comp.set_xlabel(f"{xtitle} $\otimes$ {ytitle}" if histo.dense_dim() == 2 else xtitle)
                ax_comp.set_xlim((0, len(hline.axis(hline.dense_axes()[0].name).edges())-1) if histo.dense_dim() == 2 else x_lims)
                ax_comp.axhline(y=0., color="k", linestyle="-")

                if histo.dense_dim() == 2:
                        # draw vertical lines separating ctstar bins
                    bin_sep_lines = [histo.integrate("dataset").values()[()].shape[0]*ybin for ybin in range(1, histo.integrate("dataset").values()[()].shape[1])]
                    for binline in bin_sep_lines:
                        ax_comp.axvline(binline, ymin=0., ymax=0.90, color="k", linestyle="--")

                    # add lepton/jet multiplicity label
                ax_comp.text(
                    0.02, 0.94, "parton level",
                    fontsize=rcParams["font.size"]*0.90, horizontalalignment="left", verticalalignment="bottom", transform=ax_comp.transAxes
                )
                hep.cms.label(ax=ax_comp, data=False, year=args.year, lumi=round((data_lumi_year[f"Muons"]+data_lumi_year[f"Electrons"])/2000., 1))

                #set_trace()
                sig_figname = os.path.join(signal_dir, "_".join([signal, hname, "Comparison"]))
                fig_comp.savefig(sig_figname)
                plt.close(fig_comp)
                print(f"{sig_figname} written")
        
            else:
                fig, ax = plt.subplots()
                fig.subplots_adjust(hspace=.07)

                if "Int" in signal:
                        # plot neg and pos dists
                    plot.plot1d(
                        hline if histo.dense_dim() == 2 else sig_histo,
                        ax=ax, clear=False, stack=False, line_opts={"color" : ["#e41a1c", "#377eb8"]}, # red and blue
                    )
                    # plot combined contribution
                plot.plot1d(
                    hline.integrate("process") if histo.dense_dim() == 2 else sig_histo.integrate("dataset"),
                    ax=ax, clear=False, stack=False, line_opts={"color" : ["k"]},
                )

                    # format axes/legend
                handles, labels = ax.get_legend_handles_labels()
                for idx, sample in enumerate(labels):
                    if "neg" in sample:
                        labels[idx] = "w$<$0"
                    elif "pos" in sample:
                        labels[idx] = "w$>$0"
                    else:
                        labels[idx] = "Combined"
                ax.legend(handles, labels if "Int" in signal else [], title=sig_label, loc="upper right", title_fontsize=24)
                ax.autoscale()
                ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.15)
                ax.set_xlabel(f"{xtitle} $\otimes$ {ytitle}" if histo.dense_dim() == 2 else xtitle)
                ax.set_xlim((0, len(hline.axis(hline.dense_axes()[0].name).edges())-1) if histo.dense_dim() == 2 else x_lims)
                ax.axhline(y=0., color="k", linestyle="-")

                if histo.dense_dim() == 2:
                        # draw vertical lines separating ctstar bins
                    bin_sep_lines = [histo.integrate("dataset").values()[()].shape[0]*ybin for ybin in range(1, histo.integrate("dataset").values()[()].shape[1])]
                    for binline in bin_sep_lines:
                        ax.axvline(binline, ymin=0., ymax=0.90, color="k", linestyle="--")

                    # add lepton/jet multiplicity label
                ax.text(
                    0.02, 0.94, "parton level",
                    fontsize=rcParams["font.size"]*0.90, horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                )
                hep.cms.label(ax=ax, data=False, year=args.year, lumi=round((data_lumi_year[f"Muons"]+data_lumi_year[f"Electrons"])/2000., 1))

                #set_trace()
                sig_figname = os.path.join(signal_dir, "_".join([signal, hname]))
                fig.savefig(sig_figname)
                plt.close(fig)
                print(f"{sig_figname} written")
    

    if args.level == "Reco":
        ## hists should have 2 category axes (dataset, jet multiplicity ) followed by variable
        #for jmult in ["3Jets"]:
        for jmult in sorted(set([key[1] for key in histo.values().keys()])):
            print(*[jmult, hname], sep=", ")

            #set_trace()
            hslice = histo[:, jmult].integrate("jmult")

            if "disc" in hname:
                xtitle = orig_xtitle.replace("j", "3") if jmult == "3Jets" else orig_xtitle.replace("j", "4")

            if hname == "Lep_iso":
                x_lims = (0., 0.15) if (args.lepton == "Muon") else (0., 0.1)

            #set_trace()
            all_possible_signal = sorted(set([key[0] for key in hslice.values().keys()]))
                # combine pos/neg interference into Int
            signal_to_use = ["_".join(signal.split("_")[:-1]) for signal in all_possible_signal if "neg" in signal] + [signal for signal in all_possible_signal if "Res" in signal]
            for signal in sorted(signal_to_use):
                print(f"\t{signal}")
                signal_dir = os.path.join(outdir, args.lepton, jmult, signal)
                if not os.path.isdir(signal_dir):
                    os.makedirs(signal_dir)

                sig_label = styles[signal]["name"]
                signal_mask = re.compile(r"((?:%s*))" % signal)

                sig_histo = hslice[signal_mask] # get all components from signal (pos and neg for Int, just one for Res)
                    # make linearized view 
                if histo.dense_dim() == 2: 
                    proc_histo = sig_histo.group("dataset", hist.Cat("process", "Process", sorting="placement"), {key[0] : [key[0]] for key in sorted(sig_histo.values().keys())})
                    hline = Plotter.linearize_hist(proc_histo)


                    ## make plots comparing original MC signal to ME reweighting signal
                if (hname == "mtt_vs_tlep_ctstar_abs") and args.comparison:
                    #set_trace()
                    fig_comp, ax_comp = plt.subplots()
                    fig_comp.subplots_adjust(hspace=.07)

                    if "Int" in signal:
                        mc_sig_neg_histo = mc_hdict[jmult][args.lepton][f"{signal}_neg_nosys"].copy()
                        mc_sig_pos_histo = mc_hdict[jmult][args.lepton][f"{signal}_pos_nosys"].copy()
                        tmp_neg = mc_sig_neg_histo.copy()
                        mc_sig_histo = tmp_neg.add(mc_sig_pos_histo)
                    else:
                        mc_sig_histo = mc_hdict[jmult][args.lepton][f"{signal}_nosys"].copy()

                    if histo.dense_dim() == 2:
                        #Plotter.plot_1D(mc_sig_histo.values()[()], mc_sig_histo.dense_axes()[0].edges(), xlimits=(0, len(hline.axis(hline.dense_axes()[0].name).edges())-1), ax=ax_comp, histtype="errorbar", label="Official MC", color="#e41a1c") # red
                        #Plotter.plot_1D(hline.integrate("process").values()[()], hline.dense_axes()[0].edges(), xlimits=(0,hline.dense_axes()[0].edges()[-1]), ax=ax_comp, histtype="errorbar", label="ME Reweighting", color="#377eb8") # blue
                        Plotter.plot_1D(mc_sig_histo.values()[()], mc_sig_histo.dense_axes()[0].edges(), xlimits=(0, len(hline.axis(hline.dense_axes()[0].name).edges())-1), ax=ax_comp, histtype="step", label="Official MC", color="#e41a1c") # red
                        Plotter.plot_1D(hline.integrate("process").values()[()], hline.dense_axes()[0].edges(), xlimits=(0,hline.dense_axes()[0].edges()[-1]), ax=ax_comp, histtype="step", label="ME Reweighting", color="#377eb8") # blue
                    else:
                        mc_sig_histo = mc_sig_histo.rebin(mc_sig_histo.dense_axes()[0].name, rebinning)

                        Plotter.plot_1D(mc_sig_histo.integrate(mc_sig_histo.sparse_axes()[0].name).values()[()], mc_sig_histo.dense_axes()[0].edges(), xlimits=x_lims, ax=ax_comp, histtype="errorbar", label="Official MC", color="#e41a1c") # red
                        Plotter.plot_1D(sig_histo.integrate(sig_histo.sparse_axes()[0].name).values()[()], sig_histo.dense_axes()[0].edges(), xlimits=x_lims, ax=ax_comp, histtype="errorbar", label="ME Reweighting", color="#377eb8") # blue

                        # format axes/legend
                    ax_comp.legend(title=sig_label, loc="upper right", title_fontsize=24)
                    ax_comp.autoscale()
                    ax_comp.set_ylim(ax_comp.get_ylim()[0], ax_comp.get_ylim()[1]*1.15)
                    ax_comp.set_xlabel(f"{xtitle} $\otimes$ {ytitle}" if histo.dense_dim() == 2 else xtitle)
                    ax_comp.set_xlim((0, len(hline.axis(hline.dense_axes()[0].name).edges())-1) if histo.dense_dim() == 2 else x_lims)
                    ax_comp.axhline(y=0., color="k", linestyle="-")

                    if histo.dense_dim() == 2:
                            # draw vertical lines separating ctstar bins
                        bin_sep_lines = [hslice.integrate("dataset").values()[()].shape[0]*ybin for ybin in range(1, hslice.integrate("dataset").values()[()].shape[1])]
                        for binline in bin_sep_lines:
                            ax_comp.axvline(binline, ymin=0., ymax=0.90, color="k", linestyle="--")

                        # add lepton/jet multiplicity label
                    ax_comp.text(
                        0.02, 0.90, "%s, %s\n%s" % (objtypes["Lep"][args.lepton], jet_mults[jmult], btag_cats["btagPass"]),
                        fontsize=rcParams["font.size"]*0.90, horizontalalignment="left", verticalalignment="bottom", transform=ax_comp.transAxes
                    )
                    hep.cms.label(ax=ax_comp, data=False, year=args.year, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))

                    #set_trace()
                    sig_figname = os.path.join(signal_dir, "_".join([signal, jmult, args.lepton, hname, "Comparison"]))
                    fig_comp.savefig(sig_figname)
                    plt.close(fig_comp)
                    print(f"{sig_figname} written")
        
                else:
                    #set_trace()
                    fig, ax = plt.subplots()
                    fig.subplots_adjust(hspace=.07)

                    if "Int" in signal:
                            # plot neg and pos dists
                        plot.plot1d(
                            hline if histo.dense_dim() == 2 else sig_histo,
                            ax=ax, clear=False, stack=False, line_opts={"color" : ["#e41a1c", "#377eb8"]}, # red and blue
                        )
                        # plot combined contribution
                    plot.plot1d(
                        hline.integrate("process") if histo.dense_dim() == 2 else sig_histo.integrate("dataset"),
                        ax=ax, clear=False, stack=False, line_opts={"color" : ["k"]},
                    )

                        # format axes/legend
                    handles, labels = ax.get_legend_handles_labels()
                    for idx, sample in enumerate(labels):
                        if "neg" in sample:
                            labels[idx] = "w$<$0"
                        elif "pos" in sample:
                            labels[idx] = "w$>$0"
                        else:
                            labels[idx] = "Combined"
                    ax.legend(handles, labels if "Int" in signal else [], title=sig_label, loc="upper right", title_fontsize=24)
                    ax.autoscale()
                    ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.15)
                    ax.set_xlabel(f"{xtitle} $\otimes$ {ytitle}" if histo.dense_dim() == 2 else xtitle)
                    ax.set_xlim((0, len(hline.axis(hline.dense_axes()[0].name).edges())-1) if histo.dense_dim() == 2 else x_lims)
                    ax.axhline(y=0., color="k", linestyle="-")

                    if histo.dense_dim() == 2:
                            # draw vertical lines separating ctstar bins
                        bin_sep_lines = [hslice.integrate("dataset").values()[()].shape[0]*ybin for ybin in range(1, hslice.integrate("dataset").values()[()].shape[1])]
                        for binline in bin_sep_lines:
                            ax.axvline(binline, ymin=0., ymax=0.90, color="k", linestyle="--")

                        # add lepton/jet multiplicity label
                    ax.text(
                        0.02, 0.90, "%s, %s\n%s" % (objtypes["Lep"][args.lepton], jet_mults[jmult], btag_cats["btagPass"]),
                        fontsize=rcParams["font.size"]*0.90, horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                    )
                    hep.cms.label(ax=ax, data=False, year=args.year, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))

                    #set_trace()
                    sig_figname = os.path.join(signal_dir, "_".join([signal, jmult, args.lepton, hname]))
                    fig.savefig(sig_figname)
                    plt.close(fig)
                    print(f"{sig_figname} written")

toc = time.time()
print("Total time: %.1f" % (toc - tic))    
