import time
tic = time.time()

# matplotlib
import matplotlib
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
import Utilities.prettyjson as prettyjson
import numpy as np
from Utilities.styles import styles as hstyles
import Utilities.final_analysis_binning as final_binning
import Utilities.common_features as cfeatures

from argparse import ArgumentParser
parser = ArgumentParser()
#parser.add_argument("year", choices=["2018"], help="What year is the ntuple from.")
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("lepton", choices=["Electron", "Muon", "All"], help="Choose which lepton to make plots for")
parser.add_argument("jmult", choices=["3", "4+", "All"], help="Choose which jet multiplicity to make plots for")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
base_jobid = os.environ["base_jobid"]
jobid = os.environ["jobid"]
plot_outdir = os.environ["plots_dir"]
eos_dir = os.environ["eos_dir"]

leps_to_run = ["Electron", "Muon"] if args.lepton == "All" else [args.lepton]
if args.jmult == "3": jmults_to_run = ["3Jets"]
elif args.jmult == "4+": jmults_to_run = ["4PJets"]
elif args.jmult == "All": jmults_to_run = ["3Jets", "4PJets"]
#set_trace()
sig_to_run = None

outdir = os.path.join(plot_outdir, jobid, "Signal_Shapes_Varied_Coupling", args.year)
if not os.path.isdir(outdir):
    os.makedirs(outdir)


linearize_binning = (
    final_binning.mtt_binning,
    final_binning.ctstar_abs_binning
)

nbins = (len(linearize_binning[0]) - 1) * (len(linearize_binning[1]) - 1)
x_lims = (0, nbins)
vlines = [x_lims[1]*ybin/5 for ybin in range(1, 5)]

#set_trace()
tlep_ctstar_binlabels = [r"%s $\in [%s, %s)$" % (cfeatures.variable_names_to_labels["tlep_ctstar_abs"], linearize_binning[1][bin], linearize_binning[1][bin+1]) for bin in range(len(linearize_binning[1])-1)]
tlep_ctstar_binlabels[-1] = r"%s $\in [%s, %s]$" % (cfeatures.variable_names_to_labels["tlep_ctstar_abs"], linearize_binning[1][-2], linearize_binning[1][-1])
tlep_ctstar_bin_locs = np.linspace((len(linearize_binning[0])-1)/2, (len(linearize_binning[0])-1)*(len(linearize_binning[1])-1) - (len(linearize_binning[0])-1)/2, len(linearize_binning[1])-1)
mtt_vals_to_plot = np.array([400, 600, 1000])
mtt_tiled_labels = np.tile(linearize_binning[0][:-1], linearize_binning[1].size-1)
mtt_bin_inds_to_plot = np.where(np.in1d(mtt_tiled_labels, mtt_vals_to_plot))[0]
mtt_bins_to_plot = np.tile(mtt_vals_to_plot, linearize_binning[1].size-1)

year_to_use = cfeatures.year_labels[args.year]

    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())[args.year]

coupling_vals_dict = {
    0.1 : {"label" : "0.1", "color" : "#ff7f00"}, # orange
    0.2 : {"label" : "0.2", "color" : "#984ea3"}, # purple
    0.5 : {"label" : "0.5", "color" : "#4daf4a"}, # green
    0.8 : {"label" : "0.8", "color" : "#377eb8"}, # blue
    0.9 : {"label" : "0.9", "color" : "#e41a1c"}, # red
    1.0 : {"label" : "1.0", "color" : "k"},
}

hdict = load(os.path.join(eos_dir, "results", f"{args.year}_{jobid}", "Templates_htt_btag_sb_regions", f"raw_templates_lj_MEsig_{args.year}_{jobid}.coffea"))
for jmult in sorted(jmults_to_run):
    for lepton in sorted(leps_to_run):
        if not sig_to_run:
            sig_to_run = sorted(set([key.split("_Res_nosys")[0] for key in hdict[jmult][lepton].keys() if "Res_nosys" in key]))

        for sig in sig_to_run:
            print(args.year, jmult, lepton, sig)
            pltdir = os.path.join(outdir, lepton, jmult, sig.split("_W")[0]) # make directory for each mass point (not mass+width)
            if not os.path.isdir(pltdir):
                os.makedirs(pltdir)

            res_histo, neg_int_histo, pos_int_histo = hdict[jmult][lepton][f"{sig}_Res_nosys"].copy(), hdict[jmult][lepton][f"{sig}_Int_neg_nosys"].copy(), hdict[jmult][lepton][f"{sig}_Int_pos_nosys"].copy()
            orig_res_vals, orig_int_vals = res_histo.values()[()], neg_int_histo.values()[()]+pos_int_histo.values()[()]

            fig, ax = plt.subplots(figsize=(15.0, 10.0))
            fig.subplots_adjust(hspace=.07)
            for g_val in coupling_vals_dict.keys():
                sig_vals = orig_res_vals * (g_val**4) + orig_int_vals * (g_val**2)
                ax.step(res_histo.dense_axes()[0].edges(), np.r_[sig_vals, sig_vals[-1]], where="post", **coupling_vals_dict[g_val])

            ax.legend(loc="upper right", title="g", ncol=2)
            ax.autoscale()
            ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.1)
            ax.set_xlim(x_lims)
            ax.set_xlabel("$m_{t\\bar{t}}$ [GeV]")
            ax.set_ylabel("Events")
            ax.axhline(0., **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})

                # draw vertical lines separating ctstar bins
            [ax.axvline(vline, color="k", linestyle="--") for vline in vlines]

            [ax.annotate(label, xy=(tlep_ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                xytext=(0, 400), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0) for idx, label in enumerate(tlep_ctstar_binlabels)]

            #set_trace()
            ax.set_xticks(mtt_bin_inds_to_plot)
            ax.set_xticklabels(mtt_bins_to_plot)

                # add lepton/jet multiplicity label

            ax.text(
                0.02, 0.90, hstyles[f"{sig}_Total"]["name"]+"\n"+f"{cfeatures.channel_labels[f'{lepton}_{jmult}']}",
                fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
            )
            hep.cms.label(ax=ax, label="Preliminary", data=False, year=year_to_use, lumi=round(data_lumi_year[f"{lepton}s"]/1000., 1))

            figname = os.path.join(pltdir, "_".join([sig, jmult, lepton, "CoupingComp"]))
            fig.savefig(figname)
            print(f"{figname} written")
            plt.close(fig)
            #set_trace()


toc = time.time()
print("Total time: %.1f" % (toc - tic))
