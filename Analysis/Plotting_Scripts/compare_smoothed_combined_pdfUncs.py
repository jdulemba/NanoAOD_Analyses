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
from coffea.util import load
from pdb import set_trace
import os
from coffea import hist
import numpy as np
import Utilities.Plotter as Plotter
from coffea.hist import plot
from Utilities.styles import styles as hstyles
import Utilities.final_analysis_binning as final_binning
import Utilities.common_features as cfeatures
import Utilities.prettyjson as prettyjson

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("combination", choices=["lepton", "era_lepton", "indiv"], help="What type of combination has been performed.")
parser.add_argument("--year", choices=["2016APV", "2016", "2017", "2018"], help="Specify which year to use if combination is only across leptons.")
args = parser.parse_args()

if (args.combination == "lepton") and (not args.year):
    raise ValueError("Year must be specified when the combination is only across leptons")

if (args.combination == "indiv") and (not args.year):
    raise ValueError("Year must be specified when the individual channel templates are used")

base_jobid = os.environ["base_jobid"]
proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
analyzer = "htt_pdfUncs"
plot_outdir = os.environ["plots_dir"]
eos_dir = os.environ["eos_dir"]

year_to_use = cfeatures.year_labels[args.year] if args.year else cfeatures.year_labels["Total"]
data_lumi_dict = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())

#set_trace()
#if args.combination == "lepton":
#    data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())[args.year]
#    lumi_to_use = (data_lumi_year["Muons"]+data_lumi_year["Electrons"])/2000
#
#    orig_hdict = load(os.path.join(eos_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}", f"raw_combined_lep_templates_lj_bkg_{args.year}_{jobid}.coffea"))
#    smoothed_hdict = load(os.path.join(eos_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}", f"smoothed_combined_lep_templates_lj_bkg_{args.year}_{jobid}.coffea"))
#    #year_use = args.year

if args.combination == "era_lepton":
    orig_hdict = load(os.path.join(eos_dir, "results", jobid, f"Templates_{analyzer}", f"raw_combined_year_and_lepton_templates_lj_PDF_{jobid}.coffea"))
    smoothed_hdict = load(os.path.join(eos_dir, "results", jobid, f"Templates_{analyzer}", f"smoothed_combined_year_and_lepton_templates_lj_PDF_{jobid}.coffea"))
    lumi_to_use = int(round((data_lumi_dict["TOT"]["Muons"]+data_lumi_dict["TOT"]["Electrons"])/2000, 0))

if args.combination == "indiv":
    orig_hdict = load(os.path.join(eos_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}", f"raw_pdf_templates_lj_bkg_{args.year}_{jobid}.coffea"))
    smoothed_hdict = load(os.path.join(eos_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}", f"smoothed_pdf_templates_lj_bkg_{args.year}_{jobid}.coffea"))

outdir = os.path.join(plot_outdir, jobid, "Compare_Smoothed_Combined_Templates_htt_btag_sb_regions", "PDF", (args.combination).upper()) if args.combination == "era_lepton" \
        else os.path.join(plot_outdir, jobid, "Compare_Smoothed_Combined_Templates_htt_btag_sb_regions", "PDF", args.year, (args.combination).upper())
if not os.path.isdir(outdir):
    os.makedirs(outdir)

linearize_binning = (
    final_binning.mtt_binning,
    final_binning.ctstar_abs_binning
)
ctstar_binlabels = [r"%s $\in [%s, %s)$" % (cfeatures.variable_names_to_labels["tlep_ctstar_abs"], linearize_binning[1][bin], linearize_binning[1][bin+1]) for bin in range(len(linearize_binning[1])-1)]
ctstar_binlabels[-1] = r"%s $\in [%s, %s]$" % (cfeatures.variable_names_to_labels["tlep_ctstar_abs"], linearize_binning[1][-2], linearize_binning[1][-1])
ctstar_bin_locs = np.linspace((len(linearize_binning[0])-1)/2, (len(linearize_binning[0])-1)*(len(linearize_binning[1])-1) - (len(linearize_binning[0])-1)/2, len(linearize_binning[1])-1)
mtt_vals_to_plot = np.array([400, 600, 1000])
mtt_tiled_labels = np.tile(linearize_binning[0][:-1], linearize_binning[1].size-1)
mtt_bin_inds_to_plot = np.where(np.in1d(mtt_tiled_labels, mtt_vals_to_plot))[0]
mtt_bins_to_plot = np.tile(mtt_vals_to_plot, linearize_binning[1].size-1)

#set_trace()

if args.combination == "indiv":
    for jmult in sorted(orig_hdict.keys()):
        for lep in sorted(orig_hdict[jmult].keys()):
            lumi_to_use = round(data_lumi_dict[args.year][f"{lep}s"]/1000, 1)

            orig_nominal = orig_hdict[jmult][lep]["TT_nosys"].copy()
            x_lims = (0, orig_nominal.dense_axes()[0].centers().size)
    
                # make sure output directory exists
            pltdir = os.path.join(outdir, lep, jmult)
            if not os.path.isdir(pltdir):
                os.makedirs(pltdir)
    
            ### plot alphaS variations
            print(f"{args.year} {lep} {jmult} alphaS TT")
            fig, ax = plt.subplots(figsize=(15.0, 10.0))
            fig.subplots_adjust(hspace=.07)

            #set_trace()
            orig_alphaSup = orig_hdict[jmult][lep]["TT_alphaSUp"].copy()
            orig_alphaSdw = orig_hdict[jmult][lep]["TT_alphaSDown"].copy()
            smooth_alphaSup = smoothed_hdict[jmult][lep]["TT_alphaSUp"].copy()
            smooth_alphaSdw = smoothed_hdict[jmult][lep]["TT_alphaSDown"].copy()
   
            orig_asUp_masked_vals, orig_asUp_masked_bins = Plotter.get_ratio_arrays(num_vals=orig_alphaSup.values()[()], denom_vals=orig_nominal.values()[()], input_bins=orig_nominal.dense_axes()[0].edges())
            orig_asDw_masked_vals, orig_asDw_masked_bins = Plotter.get_ratio_arrays(num_vals=orig_alphaSdw.values()[()], denom_vals=orig_nominal.values()[()], input_bins=orig_nominal.dense_axes()[0].edges())
            smooth_asUp_masked_vals, smooth_asUp_masked_bins = Plotter.get_ratio_arrays(num_vals=smooth_alphaSup.values()[()], denom_vals=orig_nominal.values()[()], input_bins=orig_nominal.dense_axes()[0].edges())
            smooth_asDw_masked_vals, smooth_asDw_masked_bins = Plotter.get_ratio_arrays(num_vals=smooth_alphaSdw.values()[()], denom_vals=orig_nominal.values()[()], input_bins=orig_nominal.dense_axes()[0].edges())
            #set_trace()    

            ax.fill_between(orig_asUp_masked_bins, orig_asUp_masked_vals, y2=1., facecolor="r", step="post", alpha=0.5, label="Original Up")
            ax.step(smooth_asUp_masked_bins, smooth_asUp_masked_vals, where="post", **{"color": "r", "linestyle": "-", "label":"Smoothed Up"})
            ax.fill_between(orig_asDw_masked_bins, orig_asDw_masked_vals, y2=1., facecolor="b", step="post", alpha=0.5, label="Original Down")
            ax.step(smooth_asDw_masked_bins, smooth_asDw_masked_vals, where="post", **{"color": "b", "linestyle": "-", "label":"Smoothed Down"})
            ax.legend(loc="upper right", title="$\\alpha_{S}$", ncol=2)
            ax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
            ax.autoscale()
            ax.set_xlim(x_lims)
            ax.set_xlabel("$m_{t\\bar{t}}$")
            ax.set_ylabel("Ratio to Nominal")
    
                # draw vertical lines separating ctstar bins
            vlines = [x_lims[1]*ybin/5 for ybin in range(1, 5)]
            [ax.axvline(vline, color="k", linestyle="--") for vline in vlines]
    
            [ax.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                xytext=(0, 400), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0) for idx, label in enumerate(ctstar_binlabels)]
    
            #set_trace()
            ax.set_xticks(mtt_bin_inds_to_plot)
            ax.set_xticklabels(mtt_bins_to_plot)
                # add lepton/jet multiplicity label
            ax.text(
                0.02, 0.92, cfeatures.channel_labels[f"{lep}_{jmult}"],
                fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
            )
            hep.cms.label(ax=ax, label="Preliminary", data=False, year=year_to_use, lumi=lumi_to_use)
 
            #set_trace()
            figname = os.path.join(pltdir, "_".join([jmult, lep, "pdf_alphaS"]))
            fig.savefig(figname)
            print(f"{figname} written")
            plt.close(fig)
            #set_trace()

            ### plot other PDF variations
            for tname in orig_hdict[jmult][lep].keys():
                proc = tname.split("_")[0]
                sys = sorted(filter(None, tname.split(f"{proc}_")))[0]
                if (sys == "nosys") or (sys == "pdf_0"): continue
                if "alphaS" in sys: continue # plot alphaS up/down together
                print(f"{args.year} {lep} {jmult} {sys} {proc}")
    
                orig_hist = orig_hdict[jmult][lep][tname].copy()
                smooth_hist = smoothed_hdict[jmult][lep][tname].copy()
    
                orig_masked_vals, orig_masked_bins = Plotter.get_ratio_arrays(num_vals=orig_hist.values()[()], denom_vals=orig_nominal.values()[()], input_bins=orig_nominal.dense_axes()[0].edges())
                smooth_masked_vals, smooth_masked_bins = Plotter.get_ratio_arrays(num_vals=smooth_hist.values()[()], denom_vals=orig_nominal.values()[()], input_bins=orig_nominal.dense_axes()[0].edges())
    
                # plot
                fig, ax = plt.subplots(figsize=(15.0, 10.0))
                fig.subplots_adjust(hspace=.07)
    
                ax.fill_between(orig_masked_bins, orig_masked_vals, y2=1., facecolor="r", step="post", alpha=0.5, label="Original")
                ax.step(smooth_masked_bins, smooth_masked_vals, where="post", **{"color": "r", "linestyle": "-", "label":"Smoothed"})
                ax.legend(loc="upper right", title=sys)
                ax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
                ax.autoscale()
                ax.set_xlim(x_lims)
                ax.set_xlabel("$m_{t\\bar{t}}$")
                ax.set_ylabel("Ratio to Nominal")
    
                    # draw vertical lines separating ctstar bins
                vlines = [x_lims[1]*ybin/5 for ybin in range(1, 5)]
                [ax.axvline(vline, color="k", linestyle="--") for vline in vlines]
    
                [ax.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                    xytext=(0, 400), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0) for idx, label in enumerate(ctstar_binlabels)]
    
                #set_trace()
                ax.set_xticks(mtt_bin_inds_to_plot)
                ax.set_xticklabels(mtt_bins_to_plot)
                    # add lepton/jet multiplicity label
                ax.text(
                    0.02, 0.92, cfeatures.channel_labels[f"{lep}_{jmult}"],
                    fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                )
                hep.cms.label(ax=ax, label="Preliminary", data=False, year=year_to_use, lumi=lumi_to_use)
    
                #set_trace()
                figname = os.path.join(pltdir, "_".join([jmult, lep, sys]))
                fig.savefig(figname)
                print(f"{figname} written")
                plt.close(fig)
                #set_trace()


if args.combination == "era_lepton":
    #set_trace()
    for jmult in sorted(orig_hdict.keys()):
            orig_nominal = orig_hdict[jmult]["TT_nosys"].copy()
            x_lims = (0, orig_nominal.dense_axes()[0].centers().size)
    
                # make sure output directory exists
            pltdir = os.path.join(outdir, jmult)
            if not os.path.isdir(pltdir):
                os.makedirs(pltdir)

            ### plot alphaS variations
            print(f"Run 2 {jmult} alphaS TT")
            fig, ax = plt.subplots(figsize=(15.0, 10.0))
            fig.subplots_adjust(hspace=.07)

            #set_trace()
            orig_alphaSup = orig_hdict[jmult]["TT_alphaSUp"].copy()
            orig_alphaSdw = orig_hdict[jmult]["TT_alphaSDown"].copy()
            smooth_alphaSup = smoothed_hdict[jmult]["TT_alphaSUp"].copy()
            smooth_alphaSdw = smoothed_hdict[jmult]["TT_alphaSDown"].copy()
   
            orig_asUp_masked_vals, orig_asUp_masked_bins = Plotter.get_ratio_arrays(num_vals=orig_alphaSup.values()[()], denom_vals=orig_nominal.values()[()], input_bins=orig_nominal.dense_axes()[0].edges())
            orig_asDw_masked_vals, orig_asDw_masked_bins = Plotter.get_ratio_arrays(num_vals=orig_alphaSdw.values()[()], denom_vals=orig_nominal.values()[()], input_bins=orig_nominal.dense_axes()[0].edges())
            smooth_asUp_masked_vals = np.r_[smooth_alphaSup.values()[()], smooth_alphaSup.values()[()][-1]] 
            smooth_asDw_masked_vals = np.r_[smooth_alphaSdw.values()[()], smooth_alphaSdw.values()[()][-1]] 
            #set_trace()    

            ax.fill_between(orig_asUp_masked_bins, orig_asUp_masked_vals, y2=1., facecolor="r", step="post", alpha=0.5, label="Original Up")
            ax.step(orig_asUp_masked_bins, smooth_asUp_masked_vals, where="post", **{"color": "r", "linestyle": "-", "label":"Smoothed Up"})
            ax.fill_between(orig_asDw_masked_bins, orig_asDw_masked_vals, y2=1., facecolor="b", step="post", alpha=0.5, label="Original Down")
            ax.step(orig_asDw_masked_bins, smooth_asDw_masked_vals, where="post", **{"color": "b", "linestyle": "-", "label":"Smoothed Down"})
            ax.legend(loc="upper right", title="$\\alpha_{S}$", ncol=2)
            ax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
            ax.autoscale()
            ax.set_xlim(x_lims)
            ax.set_xlabel("$m_{t\\bar{t}}$")
            ax.set_ylabel("Ratio to Nominal")
    
                # draw vertical lines separating ctstar bins
            vlines = [x_lims[1]*ybin/5 for ybin in range(1, 5)]
            [ax.axvline(vline, color="k", linestyle="--") for vline in vlines]
    
            [ax.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                xytext=(0, 400), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0) for idx, label in enumerate(ctstar_binlabels)]
            ax.set_xticks(mtt_bin_inds_to_plot)
            ax.set_xticklabels(mtt_bins_to_plot)

                # add lepton/jet multiplicity label
            ax.text(
                0.02, 0.92, cfeatures.channel_labels[f"Lepton_{jmult}"],
                fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
            )
            hep.cms.label(ax=ax, label="Preliminary", data=False, year=year_to_use, lumi=lumi_to_use)
 
            #set_trace()
            figname = os.path.join(pltdir, "_".join([jmult, "pdf_alphaS"]))
            fig.savefig(figname)
            print(f"{figname} written")
            plt.close(fig)
            #set_trace()
    
            ### plot all other PDF variations
            for tname in orig_hdict[jmult].keys():
                proc = tname.split("_")[0]
                sys = sorted(filter(None, tname.split(f"{proc}_")))[0]
                if (sys == "nosys") or (sys == "pdf_0"): continue
                if "alphaS" in sys: continue # plot alphaS up/down together
                print(f"Run 2 {jmult} {sys} {proc}")
    
                orig_hist = orig_hdict[jmult][tname].copy()
                smooth_hist = smoothed_hdict[jmult][tname].copy()
    
                orig_masked_vals, orig_masked_bins = Plotter.get_ratio_arrays(num_vals=orig_hist.values()[()], denom_vals=orig_nominal.values()[()], input_bins=orig_nominal.dense_axes()[0].edges())
    
                # plot
                fig, ax = plt.subplots(figsize=(15.0, 10.0))
                fig.subplots_adjust(hspace=.07)

                #set_trace()    
                ax.fill_between(orig_masked_bins, orig_masked_vals, y2=1., facecolor="r", step="post", alpha=0.5, label="Original")
                ax.step(orig_masked_bins, np.r_[smooth_hist.values()[()], smooth_hist.values()[()][-1]], where="post", **{"color": "r", "linestyle": "-", "label":"Smoothed"})
                ax.legend(loc="upper right", title=sys)
                ax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
                ax.autoscale()
                ax.set_xlim(x_lims)
                ax.set_xlabel("$m_{t\\bar{t}}$")
                ax.set_ylabel("Ratio to Nominal")
    
                    # draw vertical lines separating ctstar bins
                vlines = [x_lims[1]*ybin/5 for ybin in range(1, 5)]
                [ax.axvline(vline, color="k", linestyle="--") for vline in vlines]
    
                [ax.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                    xytext=(0, 400), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0) for idx, label in enumerate(ctstar_binlabels)]
                #set_trace()
                ax.set_xticks(mtt_bin_inds_to_plot)
                ax.set_xticklabels(mtt_bins_to_plot)

                    # add lepton/jet multiplicity label
                ax.text(
                    0.02, 0.92, cfeatures.channel_labels[f"Lepton_{jmult}"],
                    fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                )
                hep.cms.label(ax=ax, label="Preliminary", data=False, year=year_to_use, lumi=lumi_to_use)
    
                #set_trace()
                figname = os.path.join(pltdir, "_".join([jmult, sys]))
                fig.savefig(figname)
                print(f"{figname} written")
                plt.close(fig)
                #set_trace()


toc = time.time()
print("Total time: %.1f" % (toc - tic))
