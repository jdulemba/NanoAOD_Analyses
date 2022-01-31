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
from coffea import hist
import numpy as np
import Utilities.Plotter as Plotter
import Utilities.systematics as systematics
from coffea.hist import plot

base_jobid = os.environ["base_jobid"]

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("lepton", choices=["Electron", "Muon"], help="Choose which lepton to make plots for")
parser.add_argument("process", choices=["bkg", "sig"], help="Specify which process to use.")
parser.add_argument("--indiv", action="store_true", help="Plot individual sys variations")
parser.add_argument("--comp", action="store_true", help="Plot up/down sys variations with nominal mc")
parser.add_argument("--scale_mtop3gev", action="store_true", help="Scale 3GeV mtop variations by 1/6")
parser.add_argument("--kfactors", action="store_true", help="Use signal files scaled by kfactors")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
analyzer = "htt_btag_sb_regions"

input_dir = os.path.join(proj_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}")
base_outdir = os.path.join(proj_dir, "plots", f"{args.year}_{jobid}", f"Templates_{analyzer}")
if (args.process == "sig") and (args.kfactors):
    outdir = os.path.join(base_outdir, "SIG_KFACTORS", args.lepton)
else:
    outdir = os.path.join(base_outdir, args.lepton)

if not os.path.isdir(outdir):
    os.makedirs(outdir)

jet_mults = {
    "3Jets" : "3 jets",
    "4PJets" : "4+ jets"
}

leptypes = {
    "Muon" : "$\\mu$",
    "Electron" : "$e$",
}
lepdir = "mujets" if args.lepton == "Muon" else "ejets"

baseSys = lambda sys : "_".join(sys.split("_")[:-1])

templates_names = {
    "3Jets" : {
        "bkg" : {
            "Orig" : load(os.path.join(input_dir, f"raw_templates_lj_3Jets_bkg_mtopscaled_{args.year}_{jobid}.coffea" if args.scale_mtop3gev else f"raw_templates_lj_3Jets_bkg_{args.year}_{jobid}.coffea")),
            "Smooth" : load(os.path.join(input_dir, f"smoothed_templates_lj_3Jets_bkg_mtopscaled_{args.year}_{jobid}.coffea" if args.scale_mtop3gev else f"smoothed_templates_lj_3Jets_bkg_{args.year}_{jobid}.coffea")),
            "Flat" : load(os.path.join(input_dir, f"flattened_templates_lj_3Jets_bkg_mtopscaled_{args.year}_{jobid}.coffea" if args.scale_mtop3gev else f"flattened_templates_lj_3Jets_bkg_{args.year}_{jobid}.coffea")),
            },
        "sig" : {
            "Orig" : load(os.path.join(input_dir, f"raw_templates_lj_3Jets_sig_kfactors_{args.year}_{jobid}.coffea" if args.kfactors else f"raw_templates_lj_3Jets_sig_{args.year}_{jobid}.coffea")),
            "Smooth" : load(os.path.join(input_dir, f"smoothed_templates_lj_3Jets_sig_kfactors_{args.year}_{jobid}.coffea" if args.kfactors else f"smoothed_templates_lj_3Jets_sig_{args.year}_{jobid}.coffea")),
            "Flat" : load(os.path.join(input_dir, f"flattened_templates_lj_3Jets_sig_kfactors_{args.year}_{jobid}.coffea" if args.kfactors else f"flattened_templates_lj_3Jets_sig_{args.year}_{jobid}.coffea")),
            },
    },
    "4PJets" : {
        "bkg" : {
            "Orig" : load(os.path.join(input_dir, f"raw_templates_lj_4PJets_bkg_mtopscaled_{args.year}_{jobid}.coffea" if args.scale_mtop3gev else f"raw_templates_lj_4PJets_bkg_{args.year}_{jobid}.coffea")),
            "Smooth" : load(os.path.join(input_dir, f"smoothed_templates_lj_4PJets_bkg_mtopscaled_{args.year}_{jobid}.coffea" if args.scale_mtop3gev else f"smoothed_templates_lj_4PJets_bkg_{args.year}_{jobid}.coffea")),
            "Flat" : load(os.path.join(input_dir, f"flattened_templates_lj_4PJets_bkg_mtopscaled_{args.year}_{jobid}.coffea" if args.scale_mtop3gev else f"flattened_templates_lj_4PJets_bkg_{args.year}_{jobid}.coffea")),
            },
        "sig" : {
            "Orig" : load(os.path.join(input_dir, f"raw_templates_lj_4PJets_sig_kfactors_{args.year}_{jobid}.coffea" if args.kfactors else f"raw_templates_lj_4PJets_sig_{args.year}_{jobid}.coffea")),
            "Smooth" : load(os.path.join(input_dir, f"smoothed_templates_lj_4PJets_sig_kfactors_{args.year}_{jobid}.coffea" if args.kfactors else f"smoothed_templates_lj_4PJets_sig_{args.year}_{jobid}.coffea")),
            "Flat" : load(os.path.join(input_dir, f"flattened_templates_lj_4PJets_sig_kfactors_{args.year}_{jobid}.coffea" if args.kfactors else f"flattened_templates_lj_4PJets_sig_{args.year}_{jobid}.coffea")),
            },
    },
}

data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())[args.year]

nom_styles = {"color":"k", "linestyle":"-", "label":"Nominal"}
orig_styles = {"color":"b", "linestyle":"-", "label":"Original"}
smooth_styles = {"color":"r", "linestyle":"-", "label":"Smoothed"}
flat_styles = {"color":"g", "linestyle":"-", "label":"Flattened"}

#set_trace()
    ## make plots for background templates
for jmult in templates_names.keys():
    orig_dict = templates_names[jmult][args.process]["Orig"][args.lepton]
    smoothed_dict = templates_names[jmult][args.process]["Smooth"][args.lepton]
    flattened_dict = templates_names[jmult][args.process]["Flat"][args.lepton]

        # get all keys from both files to make sure they"re the same    
    orig_keys = sorted(orig_dict.keys())
    smoothed_keys = sorted(smoothed_dict.keys())
    diff = list(set(orig_keys) - set(smoothed_keys))
    if diff:
        raise ValueError(f"Input templates for smoothed and original distributions not the same for {jmult}")

    #set_trace()
    systs = sorted(set([key.split("Res_")[1] for key in orig_keys if "Res" in key])) if args.process == "sig" else sorted(set(["_".join(key.split("_")[1:]) for key in orig_keys if not ("data_obs" in key or len(key.split("_")) == 1)]))
    if "nosys" in systs: systs.remove("nosys")

    if args.indiv:
        #set_trace()
        for sys in systs:
            if sys != "EWcorrUp": continue

            pltdir = os.path.join(outdir, jmult, "Individual", sys)
            if not os.path.isdir(pltdir):
                os.makedirs(pltdir)
                # find histograms of associated systematics and their processes
            procs_sys = [key for key in orig_keys if sys in key]
            for proc_sys in procs_sys:
                proc = proc_sys.split(f"_{sys}")[0]

                print(jmult, sys, proc)
                nominal = orig_dict[f"{proc}_nosys"]
                orig_sys = orig_dict[f"{proc}_{sys}"]
                smooth_sys = smoothed_dict[f"{proc}_{sys}"]
                flat_sys = flattened_dict[f"{proc}_{sys}"]
                if (not nominal.values()) or (not orig_sys.values()): continue

                x_lims = (0, nominal.dense_axes()[0].centers().size)
    
                fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)            
                fig.subplots_adjust(hspace=.07)
    
                    ## plot normal hists
                hep.plot.histplot(nominal.values()[()], nominal.dense_axes()[0].edges(), ax=ax, histtype="step", **nom_styles) # nosys template
                hep.plot.histplot(orig_sys.values()[()], nominal.dense_axes()[0].edges(), ax=ax, histtype="step", **orig_styles) # original template
                if not np.array_equal(smooth_sys.values()[()], orig_sys.values()[()]):
                    hep.plot.histplot(smooth_sys.values()[()], nominal.dense_axes()[0].edges(), ax=ax, histtype="step", **smooth_styles) # smoothed template
                #if not np.array_equal(flat_sys.values()[()], orig_sys.values()[()]):
                #    hep.plot.histplot(flat_sys.values()[()], nominal.dense_axes()[0].edges(), ax=ax, histtype="step", **flat_styles) # flattened template
                ax.legend(loc="upper right", title=f"{sys}, {proc}")
                ax.set_ylabel("Events")
                ax.autoscale()
                ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.15)

                    # plot relative deviation
                orig_masked_vals, orig_masked_bins = Plotter.get_ratio_arrays(num_vals=orig_sys.values()[()]-nominal.values()[()], denom_vals=nominal.values()[()], input_bins=nominal.dense_axes()[0].edges())
                rax.step(orig_masked_bins, orig_masked_vals, where='post', **orig_styles)
                if not np.array_equal(smooth_sys.values()[()], orig_sys.values()[()]):
                    smooth_masked_vals, smooth_masked_bins = Plotter.get_ratio_arrays(num_vals=smooth_sys.values()[()]-nominal.values()[()], denom_vals=nominal.values()[()], input_bins=nominal.dense_axes()[0].edges())
                    rax.step(smooth_masked_bins, smooth_masked_vals, where='post', **smooth_styles)
                #if not np.array_equal(flat_sys.values()[()], orig_sys.values()[()]):
                #    flat_masked_vals, flat_masked_bins = Plotter.get_ratio_arrays(num_vals=flat_sys.values()[()]-nominal.values()[()], denom_vals=nominal.values()[()], input_bins=nominal.dense_axes()[0].edges())
                #    rax.step(flat_masked_bins, flat_masked_vals, where='post', **flat_styles)
                rax.axhline(0, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
                rax.autoscale()
                rax.set_xlim(x_lims)
                rax.set_xlabel("$m_{t\\bar{t}}$ $\otimes$ |cos($\\theta^{*}_{t_{l}}$)|")
                rax.set_ylabel("Rel. Deviaton from Nominal")
 
                    # add lepton/jet multiplicity label
                ax.text(
                    0.02, 0.94, f"{leptypes[args.lepton]}, {jet_mults[jmult]}",
                    fontsize=rcParams["font.size"]*0.9, horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                )
                    ## draw vertical lines for distinguishing different ctstar bins
                vlines = [x_lims[1]*ybin/5 for ybin in range(1, 5)]
                for vline in vlines:
                    ax.axvline(vline, color="k", linestyle="--")
                    rax.axvline(vline, color="k", linestyle="--")
                hep.cms.label(ax=ax, data=False, year=args.year, lumi=round(data_lumi_year["%ss" % args.lepton]/1000., 1))
                
                #set_trace()
                figname = os.path.join(pltdir, "_".join([jmult, args.lepton, sys, proc, "Sys_Comp"]))
                fig.savefig(figname)
                print(f"{figname} written")
                plt.close()

    
    if args.comp:
        #set_trace()
        systypes = systematics.sys_groups[args.year].keys()
        #systypes = ["deltaQCDdeltaEW"]
        for sys in systypes:
                # find histograms of associated systematics and their processes
            up_sysname = systematics.sys_groups[args.year][sys][0]
            dw_sysname = systematics.sys_groups[args.year][sys][1]
            procs_sys = sorted(set([key.split(f"_{up_sysname}")[0] for key in orig_keys if up_sysname in key])) if not dw_sysname \
                else sorted(set([key.split(f"_{up_sysname}")[0] for key in orig_keys if up_sysname in key] + [key.split(f"_{dw_sysname}")[0] for key in orig_keys if dw_sysname in key]))

            if not procs_sys: continue

            pltdir = os.path.join(outdir, jmult, "Comp", sys)
            if not os.path.isdir(pltdir):
                os.makedirs(pltdir)

            for proc in procs_sys:
                print(jmult, sys, proc)

                nominal = orig_dict[f"{proc}_nosys"]
                orig_up = orig_dict[f"{proc}_{up_sysname}"]
                orig_dw = orig_dict[f"{proc}_{dw_sysname}"] if dw_sysname is not None else None
                smooth_up = smoothed_dict[f"{proc}_{up_sysname}"]
                smooth_dw = smoothed_dict[f"{proc}_{dw_sysname}"] if dw_sysname is not None else None
                flat_up = flattened_dict[f"{proc}_{up_sysname}"]
                flat_dw = flattened_dict[f"{proc}_{dw_sysname}"] if dw_sysname is not None else None
                if (not nominal.values()) or (not orig_up.values()): continue
                up_histos = [(orig_up, {"color": "r", "linestyle": "-", "label": "Up"}, False)] if np.array_equal(smooth_up.values()[()], orig_up.values()[()]) else \
                    [
                        (orig_up, {"color": "r", "linestyle": "--", "label": "Original Up"}, True),
                        (smooth_up, {"color": "r", "linestyle": "-", "label": "Smooth Up"}, False),
                        (flat_up, {"color": "r", "linestyle": "--", "label": "Flat Up"}, False),
                    ]

                if orig_dw:
                    dw_histos = [(orig_dw, {"color": "b", "linestyle": "-", "label": "Down"}, False)] if np.array_equal(smooth_dw.values()[()], orig_dw.values()[()]) else \
                        [
                            (orig_dw, {"color": "b", "linestyle": "--", "label": "Original Down"}, True),
                            (smooth_dw, {"color": "b", "linestyle": "-", "label": "Smooth Down"}, False),
                            (flat_dw, {"color": "b", "linestyle": "--", "label": "Flat Down"}, False),
                        ]
                else: dw_histos = None

                x_lims = (0, nominal.dense_axes()[0].centers().size)
    
                fig, ax = plt.subplots()
                fig.subplots_adjust(hspace=.07)

                    ## plot relative deviations
                for up_histo, up_style, use_fill_between in up_histos:
                        # there is at least one actual value
                    if np.any(~np.isnan(up_histo.values()[()])):
                        up_masked_vals, up_masked_bins = Plotter.get_ratio_arrays(num_vals=up_histo.values()[()]-nominal.values()[()], denom_vals=nominal.values()[()], input_bins=nominal.dense_axes()[0].edges())
                        ax.fill_between(up_masked_bins, up_masked_vals, facecolor=up_style["color"], step="post", alpha=0.5) if use_fill_between else ax.step(up_masked_bins, up_masked_vals, where="post", **up_style)

                if dw_histos:
                    for dw_histo, dw_style, use_fill_between in dw_histos:
                            # there is at least one actual value
                        if np.any(~np.isnan(dw_histo.values()[()])):
                            dw_masked_vals, dw_masked_bins = Plotter.get_ratio_arrays(num_vals=dw_histo.values()[()]-nominal.values()[()], denom_vals=nominal.values()[()], input_bins=nominal.dense_axes()[0].edges())
                            ax.fill_between(dw_masked_bins, dw_masked_vals, facecolor=dw_style["color"], step="post", alpha=0.5) if use_fill_between else ax.step(dw_masked_bins, dw_masked_vals, where="post", **dw_style)

                ax.legend(loc="upper right", title=f"{sys}, {proc}")
                ax.axhline(0, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
                ax.autoscale()
                ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.15)
                ax.set_xlim(x_lims)
                ax.set_xlabel("$m_{t\\bar{t}}$ $\otimes$ |cos($\\theta^{*}_{t_{l}}$)|")
                ax.set_ylabel("Rel. Deviaton from Nominal")
                
                    # add lepton/jet multiplicity label
                ax.text(
                    0.02, 0.94, f"{leptypes[args.lepton]}, {jet_mults[jmult]}",
                    fontsize=rcParams["font.size"]*0.9, horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                )
                    ## draw vertical lines for distinguishing different ctstar bins
                vlines = [x_lims[1]*ybin/5 for ybin in range(1, 5)]
                for vline in vlines:
                    ax.axvline(vline, color="k", linestyle="--")
                hep.cms.label(ax=ax, data=False, year=args.year, lumi=round(data_lumi_year["%ss" % args.lepton]/1000., 1))
                
                #if (("MTOP" in sys.upper()) and (args.scale_mtop3gev)): set_trace()
                figname = os.path.join(pltdir, "_".join([jmult, args.lepton, sys, "scaled", proc, "SysTemplates_Comp"])) if (("MTOP" in sys.upper()) and (args.scale_mtop3gev)) else os.path.join(pltdir, "_".join([jmult, args.lepton, sys, proc, "SysTemplates_Comp"]))
                #figname = os.path.join(pltdir, "_".join([jmult, args.lepton, sys, "scaled", proc, "SysTemplates_OriginalSys"])) if (("MTOP" in sys.upper()) and (args.scale_mtop3gev)) else os.path.join(pltdir, "_".join([jmult, args.lepton, sys, proc, "SysTemplates_OriginalSys"]))
                #set_trace()
                fig.savefig(figname)
                print(f"{figname} written")
                plt.close()
