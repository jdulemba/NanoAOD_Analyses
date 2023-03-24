#! /bin/env python

import time
tic = time.time()

# matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend("agg")
from matplotlib import rcParams
rcParams["font.size"] = 20
rcParams["savefig.format"] = "png"
rcParams["savefig.bbox"] = "tight"
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
import Utilities.common_features as cfeatures
import Utilities.final_analysis_binning as final_binning

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("lepton", choices=["Electron", "Muon"], help="Choose which lepton to make plots for")
parser.add_argument("--hem", choices=["none", "scaled", "removed"], default="not", help="Make plots when HEM correction is applied or not (default is not applied).")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
base_jobid = os.environ["base_jobid"]
analyzer = "hem_invest"
plot_outdir = os.environ["plots_dir"]
eos_dir = os.environ["eos_dir"]

year_to_use = "2018"
input_dir = os.path.join(eos_dir, "results", f"{year_to_use}_{jobid}", analyzer)
f_ext = "TOT.coffea"

if args.hem == "none":
    hem_odir = "NoTreatment"
    fname_str = "NoTreatment"
    treatment = "No Treatment"
if args.hem == "scaled":
    hem_odir = "Scaled"
    fname_str = "SCALED"
    treatment = "Objects Scaled"
if args.hem == "removed":
    hem_odir = "Removed"
    fname_str = "REMOVE"
    treatment = "Objects Removed"
    
#set_trace()
fnames = [os.path.join(input_dir, fname) for fname in os.listdir(input_dir) if (fname.endswith(f_ext) and (fname_str in fname))]
if len(fnames) > 1: raise ValueError("Too many input files found")
hdict = load(fnames[0])

outdir = os.path.join(plot_outdir, f"{year_to_use}_{jobid}", analyzer, hem_odir)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

# binning and bin labels for mtt x ctstar
linearize_binning = (
    final_binning.mtt_binning,
    final_binning.ctstar_abs_binning
)
#set_trace()
ctstar_binlabels = [r"|$\cos (\theta^{*}_{t_{\ell}})$| $\in [%s, %s)$" % (linearize_binning[1][bin], linearize_binning[1][bin+1]) for bin in range(len(linearize_binning[1])-1)]
ctstar_binlabels[-1] = r"|$\cos (\theta^{*}_{t_{\ell}})$| $\in [%s, %s]$" % (linearize_binning[1][-2], linearize_binning[1][-1])
ctstar_bin_locs = np.linspace((len(linearize_binning[0])-1)/2, (len(linearize_binning[0])-1)*(len(linearize_binning[1])-1) - (len(linearize_binning[0])-1)/2, len(linearize_binning[1])-1)
#mtt_binlabels = ["%s-%s" % (linearize_binning[0][bin], linearize_binning[0][bin+1]) for bin in range(len(linearize_binning[0])-1)]*len(ctstar_binlabels)
#mtt_bins_to_plot = np.array([400., 500., 600., 800., 1200.])
#set_trace()
mtt_vals_to_plot = np.array([400, 600, 1000])
mtt_tiled_labels = np.tile(linearize_binning[0][:-1], linearize_binning[1].size-1)
mtt_bin_inds_to_plot = np.where(np.in1d(mtt_tiled_labels, mtt_vals_to_plot))[0]
mtt_bins_to_plot = np.tile(mtt_vals_to_plot, linearize_binning[1].size-1)

# binning and bin labels for phi x eta
phi_eta_binning = (
    #np.array([-3.2, -2.4, -1.6, -0.85, 0., 0.85, 1.6, 2.4, 3.2]), # phi binning
    np.array([-3.2, -2.4, -1.57, -0.87, 0., 0.87, 1.57, 2.4, 3.2]), # phi
    np.array([-2.5, -1.3, -0.7, 0., 0.7, 1.3, 2.5]) # eta binning
)
eta_binlabels = ["%s $\leq$ $\\eta$ $\leq$ %s" % (phi_eta_binning[1][bin], phi_eta_binning[1][bin+1]) for bin in range(len(phi_eta_binning[1])-1)]
eta_bin_locs = np.linspace((len(phi_eta_binning[0])-1)/2, (len(phi_eta_binning[0])-1)*(len(phi_eta_binning[1])-1) - (len(phi_eta_binning[0])-1)/2, len(phi_eta_binning[1])-1)
phi_binlabels = ["%s $\leq$ $\\phi$ $\leq$ %s" % (phi_eta_binning[0][bin], phi_eta_binning[0][bin+1]) for bin in range(len(phi_eta_binning[0])-1)]*len(eta_binlabels)


variables = {
    "mtt_vs_tlep_ctstar_abs" : ("$m_{t\\bar{t}}$ [GeV]", "|cos($\\theta^{*}_{t_{\ell}}$)|", linearize_binning[0], linearize_binning[1], (linearize_binning[0][0], linearize_binning[0][-1]),
        (linearize_binning[1][0], linearize_binning[1][-1]), True, None, (ctstar_binlabels, ctstar_bin_locs)),
    ##"mtt_vs_tlep_ctstar_abs" : ("$m_{t\\bar{t}}$ [GeV]", "|cos($\\theta^{*}_{t_{l}}$)|", linearize_binning[0], linearize_binning[1], (linearize_binning[0][0], linearize_binning[0][-1]),
    ##    (linearize_binning[1][0], linearize_binning[1][-1]), True, None, (ctstar_binlabels, ctstar_bin_locs)),
    "Jets_phi_vs_eta" : ("$\\phi$(jets)", "$\\eta$(jets)", phi_eta_binning[0], phi_eta_binning[1], (-3.3, 3.3), (-2.6, 2.6), True, phi_binlabels, (eta_binlabels, eta_bin_locs)),
    "Lep_phi_vs_eta" : ("$\\phi$(%s)" % cfeatures.objtypes["Lep"][args.lepton], "$\\eta$(%s)" % cfeatures.objtypes["Lep"][args.lepton], phi_eta_binning[0], phi_eta_binning[1],
        (-3.3, 3.3), (-2.6, 2.6), True, phi_binlabels, (eta_binlabels, eta_bin_locs)),
    ##"Jets_njets" : ("$n_{jets}$", 1, (0, 15), True),
    ##"mtt" : ("m($t\\bar{t}$) [GeV]", 4, (200., 2000.), True),
    ##"tlep_ctstar_abs" : ("|cos($\\theta^{*}_{t_{l}}$)|", 1, (0., 1.), True),
    ##"mthad" : ("m($t_{h}$) [GeV]", 2, (0., 300.), True),
    ##"pt_thad" : ("$p_{T}$($t_{h}$) [GeV]", 2, (0., 500.), True),
    ##"pt_tlep" : ("$p_{T}$($t_{l}$) [GeV]", 2, (0., 500.), True),
    ##"pt_tt" : ("$p_{T}$($t\\bar{t}$) [GeV]", 2, (0., 500.), True),
    ##"eta_thad" : ("$\\eta$($t_{h}$)", 2, (-4., 4.), True),
    ##"eta_tlep" : ("$\\eta$($t_{l}$)", 2, (-4., 4.), True),
    ##"eta_tt" : ("$\\eta$($t\\bar{t}$)", 2, (-4., 4.), True),
    ##"tlep_ctstar" : ("cos($\\theta^{*}_{t_{l}}$)", 2, (-1., 1.), True),
    ##"full_disc" : ("$\\lambda_{C}$", 2, (5, 25.), True),
    ##"mass_disc" : ("$\\lambda_{M}$", 2, (0, 20.), True),
    ##"ns_disc" : ("$\\lambda_{NS}$", 2, (3., 10.), True),
    ##"ns_dist" : ("$D_{\\nu, min}$", 1, (0., 150.), True),
    ##"Jets_pt" : ("$p_{T}$(jets) [GeV]", 2, (0., 300.), True),
    ##"Jets_eta" : ("$\\eta$(jets)", 2, (-2.6, 2.6), True),
    ##"Jets_phi" : ("$\\phi$(jets)", 2, (-4., 4.), True),
    ##"Jets_LeadJet_pt" : ("$p_{T}$(leading jet) [GeV]", 2, (0., 300.), True),
    ##"Jets_LeadJet_eta" : ("$\\eta$(leading jet)", 2, (-2.6, 2.6), True),
    ##"Lep_pt" : ("$p_{T}$(%s) [GeV]" % objtypes["Lep"][args.lepton], 2, (0., 300.), True),
    ##"Lep_eta" : ("$\\eta$(%s)" % objtypes["Lep"][args.lepton], 2, (-2.6, 2.6), True),
    ##"Lep_phi" : ("$\\phi$(%s)" % objtypes["Lep"][args.lepton], 2, (-4, 4), True),
    ##"Lep_iso" : ("pfRelIso, %s" % objtypes["Lep"][args.lepton], 1, (0., 1.), True),
    ##"MT" : ("$M_{T}$ [GeV]", 1, (0., 300.), True),
    ##"MET_pt" : ("$p_{T}$(MET) [GeV]", 1, (0., 300.), True),
    ##"MET_phi" : ("$\phi$(MET)", 1, (-3.2, 3.2), True),
}


    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())[year_to_use]
lumi_correction = load(os.path.join(proj_dir, "Corrections", base_jobid, "MC_LumiWeights.coffea"))[year_to_use][f"{args.lepton}s"]
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
process_groups = plt_tools.make_dataset_groups(args.lepton, year_to_use, samples=names, bkgdict="dataset")

    # scale and group hists by process
for hname in hdict.keys():
    if "cutflow" in hname: continue
    #set_trace()
    hdict[hname].scale(lumi_correction, axis="dataset") # scale hists
    hdict[hname] = hdict[hname].group(process_cat, process, process_groups) # group by process
    hdict[hname] = hdict[hname][:, :, args.lepton].integrate("leptype") # only pick out specified lepton


## make nosys plots
for hname in variables.keys():
    if hname not in hdict.keys():
        raise ValueError(f"{hname} not found in file")
    histo = hdict[hname].copy()
    #set_trace()

    if histo.dense_dim() == 1:
        set_trace()
        xtitle, rebinning, x_lims, withData = variables[hname]
        if rebinning != 1:
            xaxis_name = histo.dense_axes()[0].name
            histo = histo.rebin(xaxis_name, rebinning)

        ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
        for jmult in sorted(set([key[1] for key in histo.values().keys()])):
            for lepcat in sorted(set([key[3] for key in histo.values().keys()])):
            #for lepcat in ["Tight"]:
                #for btagregion in ["btagPass"]:
                for btagregion in sorted(set([key[2] for key in histo.values().keys()])):
                    pltdir = os.path.join(outdir, args.lepton, jmult, lepcat, btagregion)
                    if not os.path.isdir(pltdir):
                        os.makedirs(pltdir)

                    print(", ".join([jmult, lepcat, btagregion, hname])) 
                    fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                    fig.subplots_adjust(hspace=.07)

                    hslice = histo[:, jmult, btagregion, lepcat].integrate("jmult").integrate("lepcat").integrate("btag")

                    if hname == "Lep_iso":
                        if args.lepton == "Muon":
                            x_lims = (0., 0.15) if lepcat == "Tight" else (0.15, 1.)
                        if args.lepton == "Electron":
                            x_lims = (0., 0.1) if lepcat == "Tight" else (0., 0.5)

                    mc_opts = {
                        #"mcorder" : ["QCD", "EWK", "singlet", "ttJets"] if not ttJets_cats else ["QCD", "EWK", "singlet", "ttJets_other", "ttJets_unmatchable", "ttJets_matchable", "ttJets_right"]
                        #"maskData" : not withData
                    }

                    Plotter.plot_stack1d(ax, rax, hslice, xlabel=xtitle, xlimits=x_lims, **mc_opts)

                    #set_trace() 
                    if hname == "Jets_njets":
                        print(jmult)
                        yields_txt, yields_json = Plotter.get_samples_yield_and_frac(hslice, data_lumi_year["%ss" % args.lepton]/1000., promptmc=True)
                        frac_name = "%s_yields_and_fracs" % "_".join([jmult, args.lepton, lepcat, btagregion])
                        plt_tools.print_table(yields_txt, filename="%s/%s.txt" % (pltdir, frac_name), print_output=True)
                        print("%s/%s.txt written" % (pltdir, frac_name))
                        with open("%s/%s.json" % (pltdir, frac_name), "w") as out:
                            out.write(prettyjson.dumps(yields_json))

                        # add lepton/jet multiplicity label
                    #set_trace()
                    ax.text(
                        0.02, 0.88, "%s, %s\n%s" % (lep_cats[lepcat], jet_mults[jmult], btag_cats[btagregion]),
                        fontsize=rcParams["font.size"]*0.75, horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                    )
                    hep.cms.label(ax=ax, data=withData, paper=False, year=year_to_use, lumi=round(data_lumi_year["%ss" % args.lepton]/1000., 1))

                    #set_trace()
                    figname = os.path.join(pltdir, "_".join([jmult, args.lepton, lepcat, btagregion, hname]))
                    fig.savefig(figname)
                    print("%s written" % figname)
                    plt.close()
    
    if histo.dense_dim() == 2:
        xtitle, ytitle, xrebinning, yrebinning, x_lims, y_lims, withData, plot_xlabels, plot_ylabels = variables[hname]

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

        ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
        for jmult in sorted(set([key[1] for key in histo.values().keys()])):
            print(f"{hname} {jmult}")
            pltdir = os.path.join(outdir, args.lepton, jmult)
            if not os.path.isdir(pltdir):
                os.makedirs(pltdir)

            withData = True
            #withData = False if hname == "mtt_vs_tlep_ctstar_abs" else variables[hname][-1]

            hslice = histo[:, jmult].integrate("jmult")

            mc_opts = {
                "maskData" : not withData
            }

                # make 1D projection along dense axes
            for dax in range(2):
                if dax == 0:
                    xlabel = ytitle
                    xlimits = y_lims
                    figname = os.path.join(pltdir, "_".join([args.lepton, jmult, yaxis_name, "proj"]))
                else:
                    xlabel = xtitle
                    #xlabel = "%s [GeV]" % xtitle
                    xlimits = x_lims
                    figname = os.path.join(pltdir, "_".join([args.lepton, jmult, xaxis_name, "proj"]))

                fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                fig.subplots_adjust(hspace=.07)

                hproj = hslice.integrate(hslice.dense_axes()[dax].name)

                Plotter.plot_stack1d(ax, rax, hproj, xlabel=xlabel, xlimits=xlimits, **{"leg_ncols" : 2})

                    # add lepton/jet multiplicity label
                ax.text(
                    0.02, 0.86, f"{cfeatures.channel_labels[f'{args.lepton}_{jmult}']}\n{treatment}",
                    fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                )
                hep.cms.label(ax=ax, label="Preliminary", data=withData, year=year_to_use, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))

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
            [ax.axvline(binline, color="k", linestyle="--") for binline in bin_sep_lines]
            if rax is not None: [rax.axvline(binline, color="k", linestyle="--") for binline in bin_sep_lines]

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
                        xytext=(0, 250 if plot_xlabels is None else 50), textcoords="offset points", va="top", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)
                        #xytext=(0, 250 if plot_xlabels is None else 120), textcoords="offset points", va="top", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)

                if hname == "mtt_vs_tlep_ctstar_abs":
                    #set_trace()
                    rax.set_xticks(mtt_bin_inds_to_plot)
                    rax.set_xticklabels(mtt_bins_to_plot)

                # add lepton/jet multiplicity label
            ax.text(
                0.02, 0.86, f"{cfeatures.channel_labels[f'{args.lepton}_{jmult}']}\n{treatment}",
                fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
            )
            hep.cms.label(ax=ax, label="Preliminary", data=withData, year=year_to_use, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))

            #set_trace()
            figname = os.path.join(pltdir, "_".join([args.lepton, jmult, hname]))
            fig.savefig(figname)
            print(f"{figname} written")
            plt.close(fig)


toc = time.time()
print("Total time: %.1f" % (toc - tic))
