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
#parser.add_argument("plot", choices=["indiv", "comp", "all"], help="Specify which systematic plots to create.")
parser.add_argument("--scale_mtop3gev", action="store_true", help="Scale 3GeV mtop variations by 1/6")
parser.add_argument("--kfactors", action="store_true", help="Use signal files scaled by kfactors")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
analyzer = "htt_btag_sb_regions"

input_dir = os.path.join(proj_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}", "FINAL")
base_outdir = os.path.join(proj_dir, "plots", f"{args.year}_{jobid}", f"Templates_{analyzer}", "FINAL")
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

baseSys = lambda sys : "_".join(sys.split("_")[:-1])

if args.process == "bkg":
    template_hdict = load(os.path.join(input_dir, f"final_templates_lj_bkg_mtopscaled_{args.year}_{jobid}.coffea" if args.scale_mtop3gev else f"final_templates_lj_bkg_{args.year}_{jobid}.coffea"))

#if args.process == "sig":
#    templates_names = {
#        "Orig" : load(os.path.join(input_dir, f"raw_templates_lj_sig_kfactors_{args.year}_{jobid}.coffea" if args.kfactors else f"raw_templates_lj_sig_{args.year}_{jobid}.coffea")),
#        "Smooth" : load(os.path.join(input_dir, f"smoothed_templates_lj_sig_kfactors_{args.year}_{jobid}.coffea" if args.kfactors else f"smoothed_templates_lj_sig_{args.year}_{jobid}.coffea")),
#        "Flat" : load(os.path.join(input_dir, f"flattened_templates_lj_sig_kfactors_{args.year}_{jobid}.coffea" if args.kfactors else f"flattened_templates_lj_sig_{args.year}_{jobid}.coffea")),
#    }

data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())[args.year]

treatment_to_name_dict = {
    "raw" : "Original",
    "smooth" : "Smoothed",
    "flat" : "Flattened",
    "symmetrized" : "Smoothed+Symmetrized",
}

    ## make plots for background templates
for jmult in template_hdict.keys():
    for lep, hdict in template_hdict[jmult].items():
        if lep != args.lepton: continue
        #set_trace()
        systypes = systematics.sys_groups[args.year].keys()
        for sys in systypes:
                # find histograms of associated systematics and their processes
            up_sysname = systematics.sys_groups[args.year][sys][0]
            dw_sysname = systematics.sys_groups[args.year][sys][1]
            procs_sys = sorted(set([key.split(f"_{up_sysname}")[0] for key in hdict.keys() if up_sysname in key])) if not dw_sysname \
                else sorted(set([key.split(f"_{up_sysname}")[0] for key in hdict.keys() if up_sysname in key] + [key.split(f"_{dw_sysname}")[0] for key in hdict.keys() if dw_sysname in key]))

            if not procs_sys: continue

            pltdir = os.path.join(outdir, jmult, "Comp", sys)
            if not os.path.isdir(pltdir):
                os.makedirs(pltdir)

            for proc in procs_sys:
                print(lep, jmult, sys, proc)

                nominal, _ = hdict[f"{proc}_nosys"]
                x_lims = (0, nominal.dense_axes()[0].centers().size)
                final_up, treatment_up = hdict[f"{proc}_{up_sysname}"]
                final_dw, treatment_dw = hdict[f"{proc}_{dw_sysname}"]
                if (not nominal.values()) or (not final_up.values()): continue
                up_histos = [
                    (final_up, {"color": "r", "linestyle": "-", "label": "Up"}, True)
                ]
                dw_histos = [
                    (final_dw, {"color": "b", "linestyle": "-", "label": "Down"}, True)
                ]

                fig, ax = plt.subplots()
                fig.subplots_adjust(hspace=.07)

                    ## plot relative deviations
                for up_histo, up_style, use_fill_between in up_histos:
                        # there is at least one actual value
                    if np.any(~np.isnan(up_histo.values()[()])):
                        up_masked_vals, up_masked_bins = Plotter.get_ratio_arrays(num_vals=up_histo.values()[()]-nominal.values()[()], denom_vals=nominal.values()[()], input_bins=nominal.dense_axes()[0].edges())
                        ax.fill_between(up_masked_bins, up_masked_vals, step="post", alpha=0.5, **up_style) if use_fill_between else ax.step(up_masked_bins, up_masked_vals, where="post", **up_style)

                if dw_histos:
                    for dw_histo, dw_style, use_fill_between in dw_histos:
                            # there is at least one actual value
                        if np.any(~np.isnan(dw_histo.values()[()])):
                            dw_masked_vals, dw_masked_bins = Plotter.get_ratio_arrays(num_vals=dw_histo.values()[()]-nominal.values()[()], denom_vals=nominal.values()[()], input_bins=nominal.dense_axes()[0].edges())
                            ax.fill_between(dw_masked_bins, dw_masked_vals, step="post", alpha=0.5, **dw_style) if use_fill_between else ax.step(dw_masked_bins, dw_masked_vals, where="post", **dw_style)

                ax.legend(loc="upper right", title=f"{sys}, {proc}\n{treatment_to_name_dict[treatment_up]}")
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
                hep.cms.label(ax=ax, data=False, year=args.year, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))
                
                figname = os.path.join(pltdir, "_".join([jmult, args.lepton, sys, "scaled", proc, "FinalSysTemplates"])) if (("MTOP" in sys.upper()) and (args.scale_mtop3gev)) else os.path.join(pltdir, "_".join([jmult, args.lepton, sys, proc, "FinalSysTemplates"]))
                #set_trace()
                fig.savefig(figname)
                print(f"{figname} written")
                plt.close()
