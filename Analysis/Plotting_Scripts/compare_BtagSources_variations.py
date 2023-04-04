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
import numpy as np
import Utilities.Plotter as Plotter
import Utilities.final_analysis_binning as final_binning
import Utilities.common_features as cfeatures
import Utilities.prettyjson as prettyjson

proj_dir = os.environ["PROJECT_DIR"]
base_jobid = os.environ["base_jobid"]
jobid = os.environ["jobid"]
analyzer = "htt_btag_sb_regions"
plot_outdir = os.environ["plots_dir"]
eos_dir = os.environ["eos_dir"]

hdicts = {year : load(os.path.join(eos_dir, "results", f"{year}_{jobid}", f"Templates_{analyzer}", f"raw_combined_lep_templates_lj_bkg_{year}_{jobid}.coffea")) for year in ["2016APV", "2016", "2017", "2018"]}

sys_groups = {
    "CMS_eff_b_13TeV_JEC" : ["btag_bc_jes_up", "btag_bc_jes_down"],
    "CMS_eff_b_13TeV_Pileup" : ["btag_bc_pileup_up", "btag_bc_pileup_down"],
    "CMS_eff_b_13TeV_Type3" : ["btag_bc_type3_up", "btag_bc_type3_down"],
    "CMS_eff_b_13TeV_Statistic" : ["btag_bc_statistic_up", "btag_bc_statistic_down"],
    "CMS_eff_b_13TeV" : ["btag_bc_up_correlated", "btag_bc_down_correlated"],
    "CMS_eff_b_13TeV_uncorr" : ["btag_bc_up_uncorrelated", "btag_bc_down_uncorrelated"],
        ## Type3 subsources
    "CMS_eff_b_13TeV_BottomFragmentation" : ["btag_bc_bfragmentation_up", "btag_bc_bfragmentation_down"],
    "CMS_eff_b_13TeV_BottomTemplateCorrection" : ["btag_bc_btempcorr_up", "btag_bc_btempcorr_down"],
    "CMS_eff_b_13TeV_JPCorrection" : ["btag_bc_cb_up", "btag_bc_cb_down"],
    "CMS_eff_b_13TeV_CharmFragmentation" : ["btag_bc_cfragmentation_up", "btag_bc_cfragmentation_down"],
    "CMS_eff_b_13TeV_CharmTemplateCorrection" : ["btag_bc_cjets_up", "btag_bc_cjets_down"],
    "CMS_eff_b_13TeV_CharmToMuonBR" : ["btag_bc_dmux_up", "btag_bc_dmux_down"],
    "CMS_eff_b_13TeV_GluonSplitting" : ["btag_bc_gluonsplitting_up", "btag_bc_gluonsplitting_down"],
    "CMS_eff_b_13TeV_AwayJetTag" : ["btag_bc_jetaway_up", "btag_bc_jetaway_down"],
    "CMS_eff_b_13TeV_VzeroParticles" : ["btag_bc_ksl_up", "btag_bc_ksl_down"],
    "CMS_eff_b_13TeV_LightCharmRatio" : ["btag_bc_l2c_up", "btag_bc_l2c_down"],
    "CMS_eff_b_13TeV_LifetimeOthers" : ["btag_bc_ltothers_up", "btag_bc_ltothers_down"],
    "CMS_eff_b_13TeV_MuonDeltaR" : ["btag_bc_mudr_up", "btag_bc_mudr_down"],
    "CMS_eff_b_13TeV_MuonPt" : ["btag_bc_mupt_up", "btag_bc_mupt_down"],
    "CMS_eff_b_13TeV_MuonRelativePt" : ["btag_bc_ptrel_up", "btag_bc_ptrel_down"],
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

outdir = os.path.join(plot_outdir, jobid, f"Compare_BtagSourcesVariations")
if not os.path.isdir(outdir):
    os.makedirs(outdir)

era_colors = {
    "2016APV" : "#e41a1c", # red
    "2016"    : "#377eb8", # blue
    "2017"    : "#4daf4a", # green
    "2018"    : "k", # black
}

data_lumi_dict = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())

for jmult in sorted(hdicts["2016APV"].keys()):
    #set_trace()

    for sysname, (sys_up, sys_dw) in sys_groups.items(): # get systematic type name, nominal application name, and jet-specific name
        fig, ax = plt.subplots(figsize=(15.0, 10.0))
        fig.subplots_adjust(hspace=.07)

        for year in hdicts.keys():
            nominal = hdicts[year][jmult]["TT_nosys"].copy()
            x_lims = (0, nominal.dense_axes()[0].centers().size)

            histo_up = hdicts[year][jmult][f"TT_{sys_up}"].copy()
            histo_dw = hdicts[year][jmult][f"TT_{sys_dw}"].copy()

            if np.any(~np.isnan(histo_up.values()[()])):
                up_masked_vals, up_masked_bins = Plotter.get_ratio_arrays(num_vals=histo_up.values()[()], denom_vals=nominal.values()[()], input_bins=histo_up.dense_axes()[0].edges())
                ax.step(up_masked_bins, up_masked_vals, where="post", **{"color": era_colors[year], "linestyle": "-", "label": f"{cfeatures.year_labels[year]} Up"})
            if np.any(~np.isnan(histo_dw.values()[()])):
                dw_masked_vals, dw_masked_bins = Plotter.get_ratio_arrays(num_vals=histo_dw.values()[()], denom_vals=nominal.values()[()], input_bins=histo_dw.dense_axes()[0].edges())
                ax.step(dw_masked_bins, dw_masked_vals, where="post", **{"color": era_colors[year], "linestyle": "--", "label": f"{cfeatures.year_labels[year]} Down"})


        ax.legend(loc="upper right", ncol=2)
        ax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
        ax.autoscale()
        #set_trace()
        
        ax.set_xlim(x_lims)
        ax.set_xlabel("$m_{t\\bar{t}}$ [GeV]")
        ax.set_ylabel("Ratio to Nominal")
        #ax.set_ylim(max(ax.get_ylim()[0], 0.65), ax.get_ylim()[1]*1.005)
        
            # add lepton/jet multiplicity label
        ax.text(
            0.02, 0.92, f"TT {sysname}, {cfeatures.channel_labels[f'Lepton_{jmult}']}",
            fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
        )
            ## draw vertical lines for distinguishing different ctstar bins
        vlines = [x_lims[1]*ybin/5 for ybin in range(1, 5)]
        [ax.axvline(vline, color="k", linestyle="--") for vline in vlines]

        for idx, label in enumerate(ctstar_binlabels):
            ax.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                xytext=(0, 25), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)

        ax.set_xticks(mtt_bin_inds_to_plot)
        ax.set_xticklabels(mtt_bins_to_plot)

        hep.cms.label(ax=ax, data=False, label="Preliminary")
        
        pltdir = os.path.join(outdir, jmult)
        if not os.path.isdir(pltdir):
            os.makedirs(pltdir)

        figname = os.path.join(pltdir, f"TT_{jmult}_{sysname}_YearComp")
        fig.savefig(figname)
        print(f"{figname} written")
        plt.close(fig)

#set_trace()
## compare uncorrelated unc to statistic unc
sys_comps = {
    "Statistic" : ["btag_bc_statistic_up", "btag_bc_statistic_down", "#e41a1c"],
    "uncorr" : ["btag_bc_up_uncorrelated", "btag_bc_down_uncorrelated", "k"],
}
for jmult in sorted(hdicts["2016APV"].keys()):
    #set_trace()

    for year in hdicts.keys():
        year_to_use = cfeatures.year_labels[year]
        lumi_to_use = (data_lumi_dict[year]["Muons"]+data_lumi_dict[year]["Electrons"])/2000

        nominal = hdicts[year][jmult]["TT_nosys"].copy()
        x_lims = (0, nominal.dense_axes()[0].centers().size)

        fig, ax = plt.subplots(figsize=(15.0, 10.0))
        fig.subplots_adjust(hspace=.07)
        for sysname, (sys_up, sys_dw, sys_col) in sys_comps.items(): # get systematic type name, nominal application name, and jet-specific name

            histo_up = hdicts[year][jmult][f"TT_{sys_up}"].copy()
            histo_dw = hdicts[year][jmult][f"TT_{sys_dw}"].copy()

            if np.any(~np.isnan(histo_up.values()[()])):
                up_masked_vals, up_masked_bins = Plotter.get_ratio_arrays(num_vals=histo_up.values()[()], denom_vals=nominal.values()[()], input_bins=histo_up.dense_axes()[0].edges())
                ax.step(up_masked_bins, up_masked_vals, where="post", **{"color": sys_col, "linestyle": "-", "label": f"{sysname} Up"})
            if np.any(~np.isnan(histo_dw.values()[()])):
                dw_masked_vals, dw_masked_bins = Plotter.get_ratio_arrays(num_vals=histo_dw.values()[()], denom_vals=nominal.values()[()], input_bins=histo_dw.dense_axes()[0].edges())
                ax.step(dw_masked_bins, dw_masked_vals, where="post", **{"color": sys_col, "linestyle": "--", "label": f"{sysname} Down"})


        ax.legend(loc="upper right", ncol=2)
        ax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
        ax.autoscale()
        
        ax.set_xlim(x_lims)
        ax.set_xlabel("$m_{t\\bar{t}}$ [GeV]")
        ax.set_ylabel("Ratio to Nominal")
        ax.set_ylim(max(ax.get_ylim()[0], 0.65), ax.get_ylim()[1]*1.005)
        
            # add lepton/jet multiplicity label
        ax.text(
            0.02, 0.92, f"TT CMS_eff_b_13TeV, {cfeatures.channel_labels[f'Lepton_{jmult}']}",
            fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
        )
            ## draw vertical lines for distinguishing different ctstar bins
        vlines = [x_lims[1]*ybin/5 for ybin in range(1, 5)]
        [ax.axvline(vline, color="k", linestyle="--") for vline in vlines]

        for idx, label in enumerate(ctstar_binlabels):
            ax.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                xytext=(0, 25), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)

        ax.set_xticks(mtt_bin_inds_to_plot)
        ax.set_xticklabels(mtt_bins_to_plot)

        hep.cms.label(ax=ax, data=False, label="Preliminary", year=year_to_use, lumi=round(lumi_to_use, 1))

        #set_trace()        
        pltdir = os.path.join(outdir, jmult, "Stat_Uncorr_Comp")
        if not os.path.isdir(pltdir):
            os.makedirs(pltdir)

        figname = os.path.join(pltdir, f"TT_{jmult}_Stat_Uncorr_Comp_{year}")
        fig.savefig(figname)
        print(f"{figname} written")
        plt.close(fig)

toc = time.time()
print("Total time: %.1f" % (toc - tic))
