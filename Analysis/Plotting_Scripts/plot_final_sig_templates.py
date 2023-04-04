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
rcParams["savefig.format"] = "pdf"
rcParams["savefig.bbox"] = "tight"
from pdb import set_trace
import os
import Utilities.plot_tools as plt_tools
import Utilities.prettyjson as prettyjson
import numpy as np
import Utilities.Plotter as Plotter
import Utilities.systematics as systematics
from Utilities.styles import styles as hstyles
import Utilities.final_analysis_binning as final_binning
import uproot
import Utilities.common_features as cfeatures

import argparse
class ParseKwargs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, dict())
        for value in values:
            key, value = value.split("=")
            getattr(namespace, self.dest)[key] = value


from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("lepton", choices=["Electron", "Muon"], help="Choose which lepton to make plots for")
parser.add_argument("--opts", nargs="*", action=ParseKwargs, help="Options to specify which mass and width points to produce for ME reweighted signal.")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
base_jobid = os.environ["base_jobid"]
jobid = os.environ["jobid"]
analyzer = "htt_btag_sb_regions"
eos_dir = os.environ["eos_dir"]
plots_dir = os.environ["plots_dir"]

if args.year == "2016APV":
    year_to_use = "2016pre"
elif args.year == "2016":
    year_to_use = "2016post"
else:
    year_to_use = args.year
year_to_use = cfeatures.year_labels[args.year]

#template_version = "v20"
template_version = "v29"
template_dname = "root://cmseos.fnal.gov//store/user/lpcbtagging/UR_ntuples/heavyhiggsinputs"

possible_masses = [365, 380, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000]
#possible_widths = [2.5, 3.0, 4.0, 5.0, 8.0, 10.0, 13.0, 15.0, 18.0, 21.0, 25.0]
possible_widths = [2.5, 5.0, 10.0, 15.0, 25.0]
#possible_widths = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 8.0, 10.0, 13.0, 15.0, 18.0, 21.0, 25.0]
possible_masses = [str(mass) for mass in possible_masses]
possible_widths = [str(width) for width in possible_widths]

#set_trace()
masses_to_run = "All" if not args.opts else args.opts.get("allowed_masses", "All")
widths_to_run = "All" if not args.opts else args.opts.get("allowed_widths", "All")

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

allowed_masses = ["m"+mass for mass in allowed_masses]
allowed_widths = ["w"+width.replace(".", "p") for width in allowed_widths]
print(f"Making target signal points for {allowed_masses} and {allowed_widths}")
#set_trace()

output_dir = os.path.join(plots_dir, jobid, f"Templates_{analyzer}", (template_version).upper(),  args.year, "FINAL", "MEreweighted_LO_SIG", args.lepton)
smoothed_template_fname = os.path.join(template_dname, f"{template_version}_mew_smoothed/templates_lj_sig_{year_to_use}.root")
smoothed_rfile = uproot.open(smoothed_template_fname)
orig_template_fname = os.path.join(template_dname, f"{template_version}_mew/templates_lj_sig_{year_to_use}.root")
rfile = uproot.open(orig_template_fname)

#set_trace()
if not os.path.isdir(output_dir):
    os.makedirs(output_dir)

jet_mults = {
    "3Jets" : "3 jets",
    "4PJets" : "4+ jets"
}

leptypes = {
    "Muon" : "$\\mu$",
    "Electron" : "$e$",
}
lepdirs_dict = {
    "Muon, 3Jets" : "mu3jets",
    "Muon, 4PJets" : "mu4pjets",
    "Electron, 3Jets" : "e3jets",
    "Electron, 4PJets" : "e4pjets",
}

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


data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())[args.year][f"{args.lepton}s"]


#set_trace()
    ## make plots for background templates
for jmult in ["3Jets", "4PJets"]:
    dirname = lepdirs_dict[f"{args.lepton}, {jmult}"]+f"_{year_to_use}"
    signal_dists = sorted([key for key in rfile[dirname].keys() for mass in allowed_masses for width in allowed_widths if ( (mass in key) and (width in key) and ("Up" not in key) and ("Down" not in key) )])

    for sigpoint in sorted(set(["_".join(key.split("_")[:3]) for key in signal_dists])): # get boson_mass_width
        proc, mass, width = sigpoint.split("_")
        #if proc == "H": continue
        print(f"\t{args.year}, {jmult} {args.lepton} {sigpoint}")

        x_edges = rfile[dirname][f"{sigpoint}_res"].axes[0].edges()
        sig_label = "$%s_{%s\ GeV}^{%s\%%}$" % (proc, mass.strip("m"), width.strip("w").replace("p", "."))

            # get original templates for ME reweighted signal
        orig_res_nosys = rfile[dirname][f"{sigpoint}_res"]
        orig_neg_nosys = rfile[dirname][f"{sigpoint}_neg"]
        orig_pos_nosys = rfile[dirname][f"{sigpoint}_pos"]

        #set_trace()
            # find uncs which are in both files
        systs = sorted(set([key.split(f"{sigpoint}_res_")[-1].split(";")[0].split("Up")[0] for key in rfile[dirname].keys() if (f"{sigpoint}_res" in key) and ("Up" in key)]).intersection(
            [key.split(f"{sigpoint}_res_")[-1].split(";")[0].split("Up")[0] for key in smoothed_rfile[dirname].keys() if (f"{sigpoint}_res" in key) and ("Up" in key)]))
        for sys in systs:
            #set_trace()
            if sys == "CMS_eff_m_reco_tot":
                print(f"\t\tSkipping {sys} for now")
                continue

            if ";" in sys: set_trace()
            print(f"\t\t{sys}")
            signal_dir = os.path.join(output_dir, jmult, sigpoint, sys)
            if not os.path.isdir(signal_dir):
                os.makedirs(signal_dir)

            #set_trace()
            histos_dict = {
                "Res" : {
                    "Up": [
                        (rfile[dirname][f"{sigpoint}_res_{sys}Up"], {"color": "r", "linestyle": "-", "label": "Original Up"}, True),
                        (smoothed_rfile[dirname][f"{sigpoint}_res_{sys}Up"], {"color": "r", "linestyle": "-", "label": "Final Up"}, False),
                    ],
                    "Down": [
                        (rfile[dirname][f"{sigpoint}_res_{sys}Down"], {"color": "b", "linestyle": "-", "label": "Original Down"}, True),
                        (smoothed_rfile[dirname][f"{sigpoint}_res_{sys}Down"], {"color": "b", "linestyle": "-", "label": "Final Down"}, False),
                    ],
                },
                "Int_neg" : {
                    "Up": [
                        (rfile[dirname][f"{sigpoint}_neg_{sys}Up"], {"color": "r", "linestyle": "-", "label": "Original Up"}, True),
                        (smoothed_rfile[dirname][f"{sigpoint}_neg_{sys}Up"], {"color": "r", "linestyle": "-", "label": "Final Up"}, False),
                    ],
                    "Down": [
                        (rfile[dirname][f"{sigpoint}_neg_{sys}Down"], {"color": "b", "linestyle": "-", "label": "Original Down"}, True),
                        (smoothed_rfile[dirname][f"{sigpoint}_neg_{sys}Down"], {"color": "b", "linestyle": "-", "label": "Final Down"}, False),
                    ],
                },
                "Int_pos" : {
                    "Up": [
                        (rfile[dirname][f"{sigpoint}_pos_{sys}Up"], {"color": "r", "linestyle": "-", "label": "Original Up"}, True),
                        (smoothed_rfile[dirname][f"{sigpoint}_pos_{sys}Up"], {"color": "r", "linestyle": "-", "label": "Final Up"}, False),
                    ],
                    "Down": [
                        (rfile[dirname][f"{sigpoint}_pos_{sys}Down"], {"color": "b", "linestyle": "-", "label": "Original Down"}, True),
                        (smoothed_rfile[dirname][f"{sigpoint}_pos_{sys}Down"], {"color": "b", "linestyle": "-", "label": "Final Down"}, False),
                    ],
                },
            }

            for sig_type in histos_dict.keys():

                if sig_type == "Res":
                    leg_title = f"{sig_label}, Res"
                        # get original templates for ME reweighted signal
                    orig_nosys = orig_res_nosys.values()
                elif sig_type == "Int_neg":
                    leg_title = f"{sig_label}, $w < 0$"
                        # get original templates for ME reweighted signal
                    orig_nosys = orig_neg_nosys.values()
                elif sig_type == "Int_pos":
                    leg_title = f"{sig_label}, $w > 0$"
                        # get original templates for ME reweighted signal
                    orig_nosys = orig_pos_nosys.values()
                #elif sig_type == "Total":
                #    leg_title = f"{sig_label}, Total"

                up_histos = histos_dict[sig_type]["Up"]
                dw_histos = histos_dict[sig_type]["Down"]


                fig, ax = plt.subplots(figsize=(15.0, 10.0))
                fig.subplots_adjust(hspace=.07)

                #set_trace()
                    ## plot relative deviations
                for up_histo, up_style, use_fill_between in up_histos:
                        # there is at least one actual value
                    if up_histo:
                        up_masked_vals, up_masked_bins = Plotter.get_ratio_arrays(num_vals=up_histo.values(), denom_vals=orig_nosys, input_bins=x_edges)
                        ax.fill_between(up_masked_bins, up_masked_vals, y2=1., facecolor=up_style["color"], step="post", alpha=0.5) if use_fill_between else ax.step(up_masked_bins, up_masked_vals, where="post", **up_style)

                for dw_histo, dw_style, use_fill_between in dw_histos:
                        # there is at least one actual value
                    if dw_histo:
                        dw_masked_vals, dw_masked_bins = Plotter.get_ratio_arrays(num_vals=dw_histo.values(), denom_vals=orig_nosys, input_bins=x_edges)
                        ax.fill_between(dw_masked_bins, dw_masked_vals, y2=1., facecolor=dw_style["color"], step="post", alpha=0.5) if use_fill_between else ax.step(dw_masked_bins, dw_masked_vals, where="post", **dw_style)


                ax.legend(loc="upper right", title=leg_title, ncol=2)
                ax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
                ax.autoscale()
                ax.set_xlim(x_edges[0], x_edges[-1])
                ax.set_xlabel(cfeatures.variable_names_to_labels["mtt"])
                ax.set_ylabel("Ratio to Nominal")

                    # add lepton/jet multiplicity label
                ax.text(
                    0.02, 0.90, f"{sys}\n{cfeatures.channel_labels[f'{args.lepton}_{jmult}']}",
                    fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                )

                    ## draw vertical lines for distinguishing different ctstar bins
                vlines = [x_edges[-1]*ybin/5 for ybin in range(1, 5)]
                [ax.axvline(vline, color="k", linestyle="--") for vline in vlines]

                [ax.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                    xytext=(0, 10), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0) for idx, label in enumerate(ctstar_binlabels)]

                ax.set_xticks(mtt_bin_inds_to_plot)
                ax.set_xticklabels(mtt_bins_to_plot)

                hep.cms.label(ax=ax, label="Preliminary", data=False, year=cfeatures.year_labels[args.year], lumi=round(data_lumi_year/1000., 1))

                #set_trace()
                figname = os.path.join(signal_dir, "_".join([sigpoint, sig_type, jmult, args.lepton, sys, "FinalSysTemplates_RatioOnly"]))
                fig.savefig(figname)
                print(f"{figname} written")
                plt.close(fig)
                #set_trace()


toc = time.time()
print("Total time: %.1f" % (toc - tic))
