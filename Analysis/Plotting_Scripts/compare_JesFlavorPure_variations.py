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
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="Specify which year to use if combination is only across leptons.")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
analyzer = "htt_btag_sb_regions"
plot_outdir = os.environ["plots_dir"]
eos_dir = os.environ["eos_dir"]

year_to_use = cfeatures.year_labels[args.year]

orig_hdict = load(os.path.join(eos_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}", f"raw_combined_lep_templates_lj_bkg_{args.year}_{jobid}.coffea"))
smoothed_hdict = load(os.path.join(eos_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}", f"smoothed_combined_lep_templates_lj_bkg_{args.year}_{jobid}.coffea"))

sys_groups = {
    "PureBottom" : ["JES_FlavorPureBottom", "JES_FlavorPureBottomOnlyBottomJets"],
    "PureCharm" : ["JES_FlavorPureCharm", "JES_FlavorPureCharmOnlyCharmJets"],
    "PureQuark" : ["JES_FlavorPureQuark", "JES_FlavorPureQuarkOnlyQuarkJets"],
    "PureGluon" : ["JES_FlavorPureGluon", "JES_FlavorPureGluonOnlyGluonJets"],
    "FlavorQCD" : ["JES_FlavorQCD", "JES_FlavorQCDOnlyLightJets"],
}


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

outdir = os.path.join(plot_outdir, jobid, f"Compare_JesFlavorPureVariations", args.year)
if not os.path.isdir(outdir):
    os.makedirs(outdir)


for jmult in sorted(orig_hdict.keys()):
    #set_trace()
    orig_nominal = orig_hdict[jmult]["TT_nosys"].copy()
    x_lims = (0, orig_nominal.dense_axes()[0].centers().size)

    for systype, (nom_sys, jet_sys) in sys_groups.items(): # get systematic type name, nominal application name, and jet-specific name
        smoothed_nom_up = smoothed_hdict[jmult][f"TT_{nom_sys}_UP"].copy()
        smoothed_jet_up = smoothed_hdict[jmult][f"TT_{jet_sys}_UP"].copy()
        smoothed_nom_dw = smoothed_hdict[jmult][f"TT_{nom_sys}_DW"].copy()
        smoothed_jet_dw = smoothed_hdict[jmult][f"TT_{jet_sys}_DW"].copy()

        #set_trace()
        fig, ax = plt.subplots(figsize=(15.0, 10.0))
        fig.subplots_adjust(hspace=.07)

            # nominal variations
        if np.any(~np.isnan(smoothed_nom_up.values()[()])):
            ax.step(smoothed_nom_up.dense_axes()[0].edges(), np.r_[smoothed_nom_up.values()[()], smoothed_nom_up.values()[()][-1]], where="post", **{"color": "r", "linestyle": "-", "label": "Smoothed All Jets Up"})
        if np.any(~np.isnan(smoothed_nom_dw.values()[()])):
            ax.step(smoothed_nom_dw.dense_axes()[0].edges(), np.r_[smoothed_nom_dw.values()[()], smoothed_nom_dw.values()[()][-1]], where="post", **{"color": "b", "linestyle": "-", "label": "Smoothed All Jets Down"})

        #set_trace()
        jet_type = "non-Bottom" if systype == "FlavorQCD" else systype.split("Pure")[-1]
            # jet-specific variations
        if np.any(~np.isnan(smoothed_jet_up.values()[()])):
            ax.step(smoothed_jet_up.dense_axes()[0].edges(), np.r_[smoothed_jet_up.values()[()], smoothed_jet_up.values()[()][-1]], where="post", **{"color": "k", "linestyle": "-", "label": f"Smoothed {jet_type} Jets Up"})
        if np.any(~np.isnan(smoothed_jet_dw.values()[()])):
            ax.step(smoothed_jet_dw.dense_axes()[0].edges(), np.r_[smoothed_jet_dw.values()[()], smoothed_jet_dw.values()[()][-1]], where="post", **{"color": "g", "linestyle": "-", "label": f"Smoothed {jet_type} Jets Down"})

        ax.legend(loc="upper right", ncol=2)
        ax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
        ax.autoscale()
        #set_trace()
        
        ax.set_xlim(x_lims)
        ax.set_xlabel("$m_{t\\bar{t}}$ [GeV]")
        ax.set_ylabel("Ratio to Nominal")
        ax.set_ylim(max(ax.get_ylim()[0], 0.65), min(ax.get_ylim()[1], 1.25))
        
            # add lepton/jet multiplicity label
        ax.text(
            0.02, 0.92, f"TT {systype}, {cfeatures.channel_labels[f'Lepton_{jmult}']}",
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

        hep.cms.label(ax=ax, data=False, label="Preliminary", year=year_to_use)
        
        #set_trace()
        pltdir = os.path.join(outdir, jmult, systype)
        if not os.path.isdir(pltdir):
            os.makedirs(pltdir)

        figname = os.path.join(pltdir, f"TT_{jmult}_{systype}_JetApplications_Comp")
        fig.savefig(figname)
        print(f"{figname} written")
        plt.close(fig)

toc = time.time()
print("Total time: %.1f" % (toc - tic))
