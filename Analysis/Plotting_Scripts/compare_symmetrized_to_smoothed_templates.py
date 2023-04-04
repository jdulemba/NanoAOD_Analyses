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
import Utilities.prettyjson as prettyjson
import numpy as np
import Utilities.Plotter as Plotter
import Utilities.systematics as systematics
from coffea.hist import plot
from Utilities.styles import styles as hstyles
import Utilities.final_analysis_binning as final_binning
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
parser.add_argument("process", choices=["bkg"], help="Specify which process to use.")
parser.add_argument("combination", choices=["era_lepton", "lepton", "indiv"], help="Specify which combination across templates to use.")
parser.add_argument("--ME_opts", nargs="*", action=ParseKwargs, help="Options to specify which mass and width points to produce for ME reweighted signal.")
args = parser.parse_args()

jobid = os.environ["jobid"]
base_jobid = os.environ["base_jobid"]
analyzer = "htt_btag_sb_regions"
plot_outdir = os.environ["plots_dir"]
eos_dir = os.environ["eos_dir"]
proj_dir = os.environ["PROJECT_DIR"]

input_dir = os.path.join(eos_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}")

outdir = os.path.join(plot_outdir, jobid, f"Compare_Symmetrized_to_Smoothed_Templates_{analyzer}", args.year, (args.process).upper(), (args.combination).upper(), args.lepton)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

year_to_use = cfeatures.year_labels[args.year]

if args.process == "bkg":
    orig_hdict = load(os.path.join(input_dir, f"raw_templates_lj_bkg_{args.year}_{jobid}.coffea"))

    if args.combination == "era_lepton":
        smooth_hdict = load(os.path.join(eos_dir, "results", jobid, f"Templates_{analyzer}", f"smoothed_combined_year_and_lepton_templates_lj_bkg_{jobid}.coffea"))
        symm_hdict = load(os.path.join(eos_dir, "results", jobid, f"Templates_{analyzer}", f"symmetrized_smoothed_combined_year_and_lepton_templates_lj_bkg_{jobid}.coffea"))
    if args.combination == "lepton":
        smooth_hdict = load(os.path.join(input_dir, f"smoothed_combined_lep_templates_lj_bkg_{args.year}_{jobid}.coffea"))
        symm_hdict = load(os.path.join(input_dir, f"symmetrized_smoothed_combined_lep_templates_lj_bkg_{args.year}_{jobid}.coffea"))
    if args.combination == "indiv":
        smooth_hdict = load(os.path.join(input_dir, f"smoothed_templates_lj_bkg_{args.year}_{jobid}.coffea"))
        symm_hdict = load(os.path.join(input_dir, f"symmetrized_smoothed_templates_lj_bkg_{args.year}_{jobid}.coffea"))
    #comb_lep_hdict = load(os.path.join(input_dir, f"raw_combined_lep_templates_lj_bkg_{args.year}_{jobid}.coffea"))
    #comb_era_lep_hdict = load(os.path.join(eos_dir, "results", jobid, f"Templates_{analyzer}", f"raw_combined_year_and_lepton_templates_lj_bkg_{jobid}.coffea"))

    #era_lep_systypes = ["EWK_scheme", "EWK_yukawa", "ISR", "FSR", "FACTOR", "RENORM", "HDAMP", "UE", "MTOP3GEV"]
    #lep_systypes = ["BTAG_BC_CORR", "BTAG_BC_UNCORR", "BTAG_L_CORR", "BTAG_L_UNCORR", "JER",
    #"JES_Absolute", f"JES_Absolute_{args.year}", "JES_BBEC1", f"JES_BBEC1_{args.year}",
    #"JES_FlavorQCD", "JES_RelativeBal", f"JES_RelativeSample_{args.year}", "MET", "PILEUP"]
    #if args.year != "2018": lep_systypes.append("PREFIRE")



data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())[args.year]

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
for jmult in sorted(orig_hdict.keys()):
    orig_dict = orig_hdict[jmult][args.lepton]

        # get all keys from both files to make sure they"re the same    
    orig_keys = sorted(orig_dict.keys())

    #set_trace()
    systypes = systematics.sys_groups[args.year].keys()
    for sys in systypes:
            # find histograms of associated systematics and their processes
        up_sysname = systematics.sys_groups[args.year][sys][0]
        dw_sysname = systematics.sys_groups[args.year][sys][1]
        procs_sys = sorted(set([key.split(f"_{up_sysname}")[0] for key in orig_keys if up_sysname in key])) if not dw_sysname \
            else sorted(set([key.split(f"_{up_sysname}")[0] for key in orig_keys if up_sysname in key] + [key.split(f"_{dw_sysname}")[0] for key in orig_keys if dw_sysname in key]))

        if not procs_sys: continue

        #if sys == "CR1": set_trace()
        if (sys == "CR1") or (sys == "CR2") or (sys == "erdON"):
            combine_sysname = systematics.combine_template_sys_to_name[args.year][f"{up_sysname}Up"].replace("Up", "")
        else:
            combine_sysname = systematics.combine_template_sys_to_name[args.year][up_sysname].replace("Up", "")
        if "LEP" in combine_sysname: combine_sysname = combine_sysname.replace("LEP", args.lepton[0].lower())
        if (args.year == "2016APV") and ("2016APV" in combine_sysname): combine_sysname = combine_sysname.replace("2016APV", "2016pre")
        if (args.year == "2016") and ("2016" in combine_sysname): combine_sysname = combine_sysname.replace("2016", "2016post")
        if "CHAN_" in combine_sysname: combine_sysname = combine_sysname.replace("CHAN_", "")

        for proc in procs_sys:
            if ((args.process == "sig") or (args.process == "MEreweight_sig")):
                if not (any([mass for mass in allowed_masses if mass in proc]) and any([width for width in allowed_widths if width in proc])): continue

            #set_trace()
            # check if same hists exist in other files
            if args.process == "MEreweight_sig":
                comb_lep_exists = (f"{proc}_{up_sysname}" in comb_lep_hdict[jmult].keys()) & (f"{proc}_{dw_sysname}" in comb_lep_hdict[jmult].keys()) if comb_lep_hdict else False
                comb_era_lep_exists = (f"{proc}_{up_sysname}" in comb_era_lep_hdict[jmult].keys()) & (f"{proc}_{dw_sysname}" in comb_era_lep_hdict[jmult].keys())
            else:
                smooth_exists = (f"{proc}_{up_sysname}" in smooth_hdict[jmult].keys()) & (f"{proc}_{dw_sysname}" in smooth_hdict[jmult].keys())
                symm_exists = (f"{proc}_{up_sysname}" in symm_hdict[jmult].keys()) & (f"{proc}_{dw_sysname}" in symm_hdict[jmult].keys())
            if not (smooth_exists and symm_exists): continue

            if not ( (f"{proc}_{up_sysname}" in orig_dict.keys() and f"{proc}_{dw_sysname}" in orig_dict.keys()) and f"{proc}_nosys" in orig_dict.keys() ): continue
            
            print(jmult, combine_sysname, proc)

            #set_trace()
                # make sure output directory exists
            pltdir = os.path.join(outdir, jmult, combine_sysname) if args.process == "bkg" else os.path.join(outdir, jmult, "_".join(proc.split("_")[:3]), combine_sysname)
            if not os.path.isdir(pltdir):
                os.makedirs(pltdir)

            up_histos, dw_histos = [], []

            #set_trace()
            orig_nominal, orig_up, orig_dw = orig_dict[f"{proc}_nosys"].copy(), orig_dict[f"{proc}_{up_sysname}"].copy(), orig_dict[f"{proc}_{dw_sysname}"].copy()
            up_histos.append((orig_nominal, orig_up, {"color": "r", "linestyle": "--", "label": f"Original {year_to_use} Up, {cfeatures.channel_labels[f'{args.lepton}_{jmult}']}"}, True))
            dw_histos.append((orig_nominal, orig_dw, {"color": "b", "linestyle": "--", "label": f"Original {year_to_use} Down, {cfeatures.channel_labels[f'{args.lepton}_{jmult}']}"}, True))

            smooth_up, smooth_dw = smooth_hdict[jmult][f"{proc}_{up_sysname}"].copy(), smooth_hdict[jmult][f"{proc}_{dw_sysname}"].copy()
            up_histos.append((None, smooth_up, {"color": "r", "linestyle": "-", "label": f"Smoothed Up"}, False))
            dw_histos.append((None, smooth_dw, {"color": "b", "linestyle": "-", "label": f"Smoothed Down"}, False))

            symm_up, symm_dw = symm_hdict[jmult][f"{proc}_{up_sysname}"].copy(), symm_hdict[jmult][f"{proc}_{dw_sysname}"].copy()
            up_histos.append((None, symm_up, {"color": "g", "linestyle": "-", "label": f"Symmetrized Up"}, False))
            dw_histos.append((None, symm_dw, {"color": "k", "linestyle": "-", "label": f"Symmetrized Down"}, False))

            x_lims = (0, orig_nominal.dense_axes()[0].centers().size)
    
            fig, ax = plt.subplots(figsize=(15.0, 10.0))
            fig.subplots_adjust(hspace=.07)

                ## plot relative deviations
            for nominal, up_histo, up_style, use_fill_between in up_histos:
                    # there is at least one actual value
                if np.any(~np.isnan(up_histo.values()[()])):
                    up_masked_vals, up_masked_bins = Plotter.get_ratio_arrays(num_vals=up_histo.values()[()], denom_vals=nominal.values()[()] if nominal else np.ones(up_histo.values()[()].size), input_bins=up_histo.dense_axes()[0].edges())
                    ax.fill_between(up_masked_bins, up_masked_vals, y2=1., facecolor=up_style["color"], step="post", alpha=0.5) if use_fill_between\
                        else ax.step(up_masked_bins, up_masked_vals, where="post", **up_style)

            for nominal, dw_histo, dw_style, use_fill_between in dw_histos:
                    # there is at least one actual value
                if np.any(~np.isnan(dw_histo.values()[()])):
                    dw_masked_vals, dw_masked_bins = Plotter.get_ratio_arrays(num_vals=dw_histo.values()[()], denom_vals=nominal.values()[()] if nominal else np.ones(dw_histo.values()[()].size), input_bins=dw_histo.dense_axes()[0].edges())
                    ax.fill_between(dw_masked_bins, dw_masked_vals, y2=1., facecolor=dw_style["color"], step="post", alpha=0.5) if use_fill_between \
                        else ax.step(dw_masked_bins, dw_masked_vals, where="post", **dw_style)

            #set_trace()
            ax.legend(loc="upper right", title=proc, ncol=2)
            #ax.legend(loc="upper right", title=proc, ncol=max(int(((len(up_histos)-1) + (len(dw_histos)-1))/2), 1))
            ax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
            ax.autoscale()
            
            ax.set_xlim(x_lims)
            ax.set_xlabel("$m_{t\\bar{t}}$ [GeV]")
            ax.set_ylabel("Ratio to Nominal")
            
                # add lepton/jet multiplicity label
            ax.text(
                0.02, 0.90, f"{combine_sysname}\n{cfeatures.channel_labels[f'{args.lepton}_{jmult}']}",
                #0.02, 0.92, combine_sysname,
                fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
            )
                ## draw vertical lines for distinguishing different ctstar bins
            vlines = [x_lims[1]*ybin/5 for ybin in range(1, 5)]
            [ax.axvline(vline, color="k", linestyle="--") for vline in vlines]

            for idx, label in enumerate(ctstar_binlabels):
                ax.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                    xytext=(0, 25), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)

            #set_trace()
            ax.set_xticks(mtt_bin_inds_to_plot)
            ax.set_xticklabels(mtt_bins_to_plot)

            hep.cms.label(ax=ax, data=False, label="Preliminary", year=year_to_use, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))
            
            #set_trace()
            figname = os.path.join(pltdir, "_".join([proc, jmult, args.lepton, combine_sysname, "Smooth_to_Symm_Comp"]))
            fig.savefig(figname)
            print(f"{figname} written")
            plt.close(fig)

toc = time.time()
print("Total time: %.1f" % (toc - tic))
