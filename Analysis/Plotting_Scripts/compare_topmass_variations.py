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
import Utilities.systematics as systematics
from coffea.hist import plot
from Utilities.styles import styles as hstyles
import Utilities.final_analysis_binning as final_binning
import Utilities.common_features as cfeatures

base_jobid = os.environ["base_jobid"]

from argparse import ArgumentParser
parser = ArgumentParser()
#parser.add_argument("combination", choices=["lepton", "era_lepton", "indiv"], help="What type of combination has been performed.")
#parser.add_argument("process", choices=["bkg"], help="Specify which process to use.")
#parser.add_argument("process", choices=["bkg", "MEreweight_sig", "sig"], help="Specify which process to use.")
#parser.add_argument("--year", choices=["2016APV", "2016", "2017", "2018"], help="Specify which year to use if combination is only across leptons.")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
analyzer = "htt_btag_sb_regions"
plot_outdir = os.environ["plots_dir"]
eos_dir = os.environ["eos_dir"]

outdir = os.path.join(plot_outdir, jobid, f"Compare_TopMassVariations", "BKG", "ERA_LEPTON")
if not os.path.isdir(outdir):
    os.makedirs(outdir)

year_to_use = cfeatures.year_labels["Total"]

#set_trace()
orig_hdict = load(os.path.join(eos_dir, "results", jobid, f"Templates_{analyzer}", f"raw_combined_year_and_lepton_templates_lj_bkg_{jobid}.coffea"))
smoothed_hdict = load(os.path.join(eos_dir, "results", jobid, f"Templates_{analyzer}", f"smoothed_combined_year_and_lepton_templates_lj_bkg_{jobid}.coffea"))

systypes = ["MTOP3GEV", "MTOP1GEV"]
year_use = "2017"


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
    orig_dict = orig_hdict[jmult]

        # get all keys from both files to make sure they"re the same    
    orig_keys = sorted(orig_dict.keys())

    #set_trace()
    orig_nominal = orig_dict["TT_nosys"].copy()
    orig_1gev_up, orig_1gev_dw = orig_dict["TT_mtop1735"].copy(), orig_dict["TT_mtop1715"].copy()
    orig_3gev_up, orig_3gev_dw = orig_dict["TT_mtop1755"].copy(), orig_dict["TT_mtop1695"].copy()

    x_lims = (0, orig_nominal.dense_axes()[0].centers().size)
    fig, ax = plt.subplots(figsize=(15.0, 10.0))
    fig.subplots_adjust(hspace=.07)

    if np.any(~np.isnan(orig_1gev_up.values()[()])):
        up_masked_vals, up_masked_bins = Plotter.get_ratio_arrays(num_vals=orig_1gev_up.values()[()], denom_vals=orig_nominal.values()[()], input_bins=orig_1gev_up.dense_axes()[0].edges())
        ax.step(up_masked_bins, up_masked_vals, where="post", **{"color": "r", "linestyle": "-", "label": f"Original 1GeV Up"})
    if np.any(~np.isnan(orig_1gev_dw.values()[()])):
        dw_masked_vals, dw_masked_bins = Plotter.get_ratio_arrays(num_vals=orig_1gev_dw.values()[()], denom_vals=orig_nominal.values()[()], input_bins=orig_1gev_dw.dense_axes()[0].edges())
        ax.step(dw_masked_bins, dw_masked_vals, where="post", **{"color": "b", "linestyle": "-", "label": f"Original 1GeV Down"})

        # scale 3GeV variations by 1/3 to compare to 1GeV
    if np.any(~np.isnan(orig_3gev_up.values()[()])):
        up_masked_vals, up_masked_bins = Plotter.get_ratio_arrays(num_vals=orig_3gev_up.values()[()], denom_vals=orig_nominal.values()[()], input_bins=orig_3gev_up.dense_axes()[0].edges())
        scaled_up_masked_vals = (up_masked_vals - 1)/3 + 1
        ax.step(up_masked_bins, scaled_up_masked_vals, where="post", **{"color": "k", "linestyle": "-", "label": f"1/3 Scaled 3GeV Up"})
    if np.any(~np.isnan(orig_3gev_dw.values()[()])):
        dw_masked_vals, dw_masked_bins = Plotter.get_ratio_arrays(num_vals=orig_3gev_dw.values()[()], denom_vals=orig_nominal.values()[()], input_bins=orig_3gev_dw.dense_axes()[0].edges())
        scaled_dw_masked_vals = (dw_masked_vals - 1)/3 + 1
        ax.step(dw_masked_bins, scaled_dw_masked_vals, where="post", **{"color": "g", "linestyle": "-", "label": f"1/3 Scaled 3GeV Down"})

    ax.legend(loc="upper right", title=cfeatures.channel_labels[f'Lepton_{jmult}'], ncol=2)
    ax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
    ax.autoscale()
    #set_trace()
    
    ax.set_xlim(x_lims)
    ax.set_xlabel("$m_{t\\bar{t}}$ [GeV]")
    ax.set_ylabel("Ratio to Nominal")
    ax.set_ylim(max(ax.get_ylim()[0], 0.65), min(ax.get_ylim()[1], 1.25))

        ## draw vertical lines for distinguishing different ctstar bins
    vlines = [x_lims[1]*ybin/5 for ybin in range(1, 5)]
    [ax.axvline(vline, color="k", linestyle="--") for vline in vlines]

    for idx, label in enumerate(ctstar_binlabels):
        ax.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
            xytext=(0, 25), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)

    #set_trace()
    ax.set_xticks(mtt_bin_inds_to_plot)
    ax.set_xticklabels(mtt_bins_to_plot)

    hep.cms.label(ax=ax, data=False, year=year_to_use)
    
    figname = os.path.join(outdir, f"TopMassVariations_Scaled3GeV_Comp_{jmult}")
    fig.savefig(figname)
    print(f"{figname} written")
    plt.close(fig)


toc = time.time()
print("Total time: %.1f" % (toc - tic))
