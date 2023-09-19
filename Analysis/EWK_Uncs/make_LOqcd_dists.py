#!/usr/bin/env python

# matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend("agg")
from matplotlib import rcParams
rcParams["font.size"] = 20
rcParams["savefig.format"] = "pdf"
rcParams["savefig.bbox"] = "tight"

import Utilities.Plotter as Plotter
from coffea import hist
from coffea.util import load, save
from pdb import set_trace
import os
import Utilities.prettyjson as prettyjson
import Utilities.plot_tools as plt_tools
import numpy as np
from coffea.lookup_tools.root_converters import convert_histo_root_file
from coffea.lookup_tools.dense_lookup import dense_lookup
from fnmatch import fnmatch
import Utilities.common_features as cfeatures

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--save_ratios", action="store_true", help="Save NNLO/tune ratios for all years (default is to not save them).")
args = parser.parse_args()

#set_trace()
proj_dir = os.environ["PROJECT_DIR"]
base_jobid = os.environ["base_jobid"]
jobid = os.environ["jobid"]
analyzer = "make_nnlo_dists"
plots_dir = os.environ["plots_dir"]
eos_dir = os.environ["eos_dir"]

f_ext = "TOT.coffea"
outdir = os.path.join(plots_dir, base_jobid, "LOqcd_dists")
if not os.path.isdir(outdir):
    os.makedirs(outdir)
    
style_dict = {
    "2016APV" : ("2016preVFP", "#984ea3"), # purple
    "2016" : ("2016postVFP", "#377eb8"), ## blue
    "2017" : ("2017", "#e41a1c"), ## red
    "2018" : ("2018", "#4daf4a"), ## green
}

tune_var = "mtt_vs_top_ctstar"
png_ext = "StatUncs"
save_dict = {} #if args.combine_years else {year: {tune_var: {}} for year in style_dict.keys()}

axis_labels_dict = {
    "mtt" : cfeatures.variable_names_to_labels["mtt"],
    "ctstar" : cfeatures.variable_names_to_labels["ctstar"],
}
    
data_lumi_dict = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())
lumi_correction = load(os.path.join(proj_dir, "Corrections", jobid, "MC_LumiWeights.coffea"))


def make_powheg_dists(mtt_binning, ctstar_binning):
        years_to_run = ["2016APV", "2016", "2017", "2018"]
        histos_dict = {}
        for year in years_to_run:
            input_dir = os.path.join(eos_dir, "results", f"{year}_{base_jobid}", analyzer)
            fnames = sorted(["%s/%s" % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)])
            hdict = plt_tools.add_coffea_files(fnames) if len(fnames) > 1 else load(fnames[0])
        
            tune_histo = hdict[tune_var]
            scale_dict = {tt : lumi_correction[year]["Muons"][tt] for tt in ["ttJetsSL", "ttJetsHad", "ttJetsDiLep"]}
            tune_histo.scale(scale_dict, axis="dataset")
    
                # orig
            tune_histo = tune_histo.integrate("dataset")
                # rebin
            tune_histo = tune_histo.rebin("mtt", hist.Bin("mtt", "mtt", mtt_binning))
            tune_histo = tune_histo.rebin("ctstar", hist.Bin("ctstar", "ctstar", ctstar_binning))
    
            histos_dict[year] = tune_histo
    
        tune_histo = histos_dict["2016APV"].add(histos_dict["2016"]).add(histos_dict["2017"]).add(histos_dict["2018"])
        tune_histo.scale(1./data_lumi_dict["TOT"]["Muons"]) # normalizes to NNLO cross section
        
        return tune_histo


lo_ztitle = "$\\sigma_{LO\ QCD}$ [pb]"
lo_diff_ztitle = "$\dfrac{d^{2} \\sigma_{LO}}{d m_{t\\bar{t}} d cos(\\theta^{*}_{t})}$"
powheg_ztitle = "$\\sigma_{Powheg}$ [pb]"
powheg_diff_ztitle = "$\dfrac{d^{2} \\sigma_{Powheg}}{d m_{t\\bar{t}} d cos(\\theta^{*}_{t})}$"

## Madgraph LO qcd
    ## get LO qcd values from root file
fname = "HATHORscaled.root"
#set_trace()
ewk_dists = convert_histo_root_file(os.path.join(proj_dir, "EWK_Uncs", fname))

# LO QCD dist is same for all directories
lo_xsec_vals, (tot_mtt_binning, tot_ctstar_binning) = ewk_dists[("EWno/tot", "dense_lookup")]
tot_mtt_binning, tot_ctstar_binning = np.around(tot_mtt_binning, decimals=0), np.around(tot_ctstar_binning, decimals=2)
rebinned_lo_xsec_vals, (rebinned_ctstar_binning, rebinned_mtt_binning) = ewk_dists[("EWno/tot_rebinned", "dense_lookup")]
rebinned_mtt_binning, rebinned_ctstar_binning = np.around(rebinned_mtt_binning, decimals=0), np.around(rebinned_ctstar_binning, decimals=2)

        # histo plotting params
max_mtt_val = 2000.
mtt_title, ctstar_title = cfeatures.variable_names_to_labels["mtt"], cfeatures.variable_names_to_labels["ctstar"]
tot_mtt_lims = (np.min(tot_mtt_binning), min(np.max(tot_mtt_binning), max_mtt_val))
tot_ctstar_lims = (np.min(tot_ctstar_binning), np.max(tot_ctstar_binning))
rebinned_mtt_lims = tot_mtt_lims
rebinned_ctstar_lims = (np.min(rebinned_ctstar_binning), np.max(rebinned_ctstar_binning))

        # only get mtt values < 2000 GeV
masked_lo_vals = lo_xsec_vals[np.where(tot_mtt_binning <= max_mtt_val)[0][0]:np.where(tot_mtt_binning <= max_mtt_val)[0][-1], :]
masked_tot_mtt_binning = tot_mtt_binning[np.where(tot_mtt_binning <= max_mtt_val)[0]]
masked_rebinned_lo_vals = rebinned_lo_xsec_vals[:, np.where(rebinned_mtt_binning <= max_mtt_val)[0][0]:np.where(rebinned_mtt_binning <= max_mtt_val)[0][-1]]
masked_rebinned_mtt_binning = rebinned_mtt_binning[np.where(rebinned_mtt_binning <= max_mtt_val)[0]]

        # get binning for differential plots
max_tot_mtt_val_bin = np.argwhere(masked_tot_mtt_binning == max_mtt_val)[0][0]
tot_mtt_binwidths = np.array([masked_tot_mtt_binning[i+1] - masked_tot_mtt_binning[i] for i in range(max_tot_mtt_val_bin)])
tot_ctstar_binwidths = np.array([tot_ctstar_binning[i+1] - tot_ctstar_binning[i] for i in range(len(tot_ctstar_binning) - 1)])
tiled_tot_mtt_binwidths = np.tile(tot_mtt_binwidths, (tot_ctstar_binwidths.size, 1)).T
tiled_tot_ctstar_binwidths = np.tile(tot_ctstar_binwidths, (tot_mtt_binwidths.size, 1))
            # rebinned
max_rebinned_mtt_val_bin = np.argwhere(masked_rebinned_mtt_binning == max_mtt_val)[0][0]
rebinned_mtt_binwidths = np.array([masked_rebinned_mtt_binning[i+1] - masked_rebinned_mtt_binning[i] for i in range(max_rebinned_mtt_val_bin)])
rebinned_ctstar_binwidths = np.array([rebinned_ctstar_binning[i+1] - rebinned_ctstar_binning[i] for i in range(len(rebinned_ctstar_binning) - 1)])
tiled_rebinned_mtt_binwidths = np.tile(rebinned_mtt_binwidths, (rebinned_ctstar_binwidths.size, 1)).T
tiled_rebinned_ctstar_binwidths = np.tile(rebinned_ctstar_binwidths, (rebinned_mtt_binwidths.size, 1))

        # get differential values
masked_lo_diff_vals = np.where(masked_lo_vals < 1e-10, .0, masked_lo_vals)/(tiled_tot_ctstar_binwidths*tiled_tot_mtt_binwidths)
masked_lo_diff_rebinned_vals = np.where(masked_rebinned_lo_vals.T < 1e-10, .0, masked_rebinned_lo_vals.T)/(tiled_rebinned_ctstar_binwidths*tiled_rebinned_mtt_binwidths)


    ### plot original LO QCD hist
fig, ax = plt.subplots()
fig.subplots_adjust(hspace=.07)
Plotter.plot_2d_norm(values=np.ma.masked_where(masked_lo_vals <= 0.0, masked_lo_vals), # mask nonzero probabilities for plotting
    xbins=masked_tot_mtt_binning, ybins=tot_ctstar_binning,
    xlimits=tot_mtt_lims, ylimits=tot_ctstar_lims, xlabel=mtt_title, ylabel=ctstar_title,
    ax=ax, **{"cmap_label" : lo_ztitle})
ax.text(
    0.98, 0.90, "$t\\bar{t}$\nparton level",
    fontsize=rcParams["font.size"], horizontalalignment="right", verticalalignment="bottom", transform=ax.transAxes, color="w",
)
hep.label.exp_label(ax=ax, exp="MADGRAPH", rlabel="")
figname = os.path.join(outdir, f"LO_QCD_xsec_mtt_vs_ctstar")
fig.savefig(figname)
print(f"{figname} written")
plt.close(fig)

    ### plot original LO QCD hist, differential
fig_diff, ax_diff = plt.subplots()
fig_diff.subplots_adjust(hspace=.07)
Plotter.plot_2d_norm(values=np.ma.masked_where(masked_lo_diff_vals <= 0.0, masked_lo_diff_vals), # mask nonzero probabilities for plotting
    xbins=masked_tot_mtt_binning, ybins=tot_ctstar_binning,
    xlimits=tot_mtt_lims, ylimits=tot_ctstar_lims, xlabel=mtt_title, ylabel=ctstar_title,
    ax=ax_diff, **{"cmap_label" : lo_diff_ztitle})
ax_diff.text(
    0.98, 0.90, "$t\\bar{t}$\nparton level",
    fontsize=rcParams["font.size"], horizontalalignment="right", verticalalignment="bottom", transform=ax_diff.transAxes, color="w",
)
hep.label.exp_label(ax=ax_diff, exp="MADGRAPH", rlabel="")
diff_figname = os.path.join(outdir, f"LO_QCD_xsec_mtt_vs_ctstar_Differential")
fig_diff.savefig(diff_figname)
print(f"{diff_figname} written")
plt.close(fig_diff)


    # plot rebinneded LO QCD hist
fig_rebinned, ax_rebinned = plt.subplots()
fig_rebinned.subplots_adjust(hspace=.07)
Plotter.plot_2d_norm(values=np.ma.masked_where(masked_rebinned_lo_vals.T <= 0.0, masked_rebinned_lo_vals.T), # mask nonzero probabilities for plotting
    xbins=masked_rebinned_mtt_binning, ybins=rebinned_ctstar_binning,
    xlimits=rebinned_mtt_lims, ylimits=rebinned_ctstar_lims, xlabel=mtt_title, ylabel=ctstar_title,
    ax=ax_rebinned, **{"cmap_label" : lo_ztitle})
ax_rebinned.text(
    0.98, 0.90, "$t\\bar{t}$\nparton level",
    fontsize=rcParams["font.size"], horizontalalignment="right", verticalalignment="bottom", transform=ax_rebinned.transAxes, color="w",
)
hep.label.exp_label(ax=ax_rebinned, exp="MADGRAPH", rlabel="")
figname_rebinned = os.path.join(outdir, f"LO_QCD_xsec_mtt_vs_ctstar_rebinned")
fig_rebinned.savefig(figname_rebinned)
print(f"{figname_rebinned} written")
plt.close(fig_rebinned)

##set_trace()
    ### plot rebinned differential LO QCD hist, normalized
fig_rebinned_diff, ax_rebinned_diff = plt.subplots()
fig_rebinned_diff.subplots_adjust(hspace=.07)
Plotter.plot_2d_norm(values=np.ma.masked_where(masked_lo_diff_rebinned_vals <= 0.0, masked_lo_diff_rebinned_vals), # mask nonzero probabilities for plotting
    xbins=masked_rebinned_mtt_binning, ybins=rebinned_ctstar_binning,
    xlimits=rebinned_mtt_lims, ylimits=rebinned_ctstar_lims, xlabel=mtt_title, ylabel=ctstar_title,
    ax=ax_rebinned_diff, **{"cmap_label" : lo_diff_ztitle})
ax_rebinned_diff.text(
    0.98, 0.90, "$t\\bar{t}$\nparton level",
    fontsize=rcParams["font.size"], horizontalalignment="right", verticalalignment="bottom", transform=ax_rebinned_diff.transAxes, color="w",
)
hep.label.exp_label(ax=ax_rebinned_diff, exp="MADGRAPH", rlabel="")
diff_rebinned_figname = os.path.join(outdir, f"LO_QCD_xsec_mtt_vs_ctstar_rebinned_Differential")
fig_rebinned_diff.savefig(diff_rebinned_figname)
print(f"{diff_rebinned_figname} written")
plt.close(fig_rebinned_diff)


    ## plot powheg dists using rebinned binning
powheg_rebinned_histo = make_powheg_dists(masked_rebinned_mtt_binning, rebinned_ctstar_binning)

    # plot powheg hist
fig, ax = plt.subplots()
fig.subplots_adjust(hspace=.07)
tune_vals = powheg_rebinned_histo.values()[()]        
tune_vals = np.where(abs(tune_vals) < 1e-10, 0, tune_vals) # set values that are less than 1e-10 to 0
Plotter.plot_2d_norm(values=np.ma.masked_where(tune_vals <= 0.0, tune_vals), # mask nonzero probabilities for plotting
    xbins=masked_rebinned_mtt_binning, ybins=rebinned_ctstar_binning,
    xlimits=rebinned_mtt_lims, ylimits=rebinned_ctstar_lims, xlabel=mtt_title, ylabel=ctstar_title,
    ax=ax, **{"cmap_label" : powheg_ztitle}
)
ax.text(
    0.98, 0.90, "$t\\bar{t}$\nparton level",
    fontsize=rcParams["font.size"], horizontalalignment="right", verticalalignment="bottom", transform=ax.transAxes, color="w",
)
hep.cms.label(ax=ax, data=False, year="Run 2")
figname = os.path.join(outdir, f"Powheg_AllYears_{base_jobid}_ttJets_xsec_{tune_var}_{png_ext}_rebinned")
fig.savefig(figname)
print(f"{figname} written")


        # plot powheg dist, differential
fig_diff, ax_diff = plt.subplots()
fig_diff.subplots_adjust(hspace=.07)
tune_diff_vals = tune_vals/(tiled_rebinned_mtt_binwidths*tiled_rebinned_ctstar_binwidths)
Plotter.plot_2d_norm(values=np.ma.masked_where(tune_diff_vals <= 0.0, tune_diff_vals), # mask nonzero probabilities for plotting
    xbins=masked_rebinned_mtt_binning, ybins=rebinned_ctstar_binning,
    xlimits=rebinned_mtt_lims, ylimits=rebinned_ctstar_lims, xlabel=mtt_title, ylabel=ctstar_title,
    ax=ax_diff, **{"cmap_label" : powheg_diff_ztitle}
)
ax_diff.text(
    0.98, 0.90, "$t\\bar{t}$\nparton level",
    fontsize=rcParams["font.size"], horizontalalignment="right", verticalalignment="bottom", transform=ax_diff.transAxes, color="w",
)
hep.cms.label(ax=ax_diff, data=False, year="Run 2")
figname_diff = os.path.join(outdir, f"Powheg_AllYears_{base_jobid}_ttJets_xsec_{tune_var}_{png_ext}_rebinned_Differential")
fig_diff.savefig(figname_diff)
print(f"{figname_diff} written")

#set_trace()
    # plot LO/powheg for combined years
fig_ratio, ax_ratio = plt.subplots()
fig_ratio.subplots_adjust(hspace=.07)
LO_to_Powheg_KFactor = masked_lo_diff_rebinned_vals/tune_diff_vals
Plotter.plot_2d_norm(values=np.ma.masked_where(LO_to_Powheg_KFactor <= 0.0, LO_to_Powheg_KFactor), # mask nonzero probabilities for plotting
    xbins=masked_rebinned_mtt_binning, ybins=rebinned_ctstar_binning,
    xlimits=rebinned_mtt_lims, ylimits=rebinned_ctstar_lims, xlabel=mtt_title, ylabel=ctstar_title,
    ax=ax_ratio, **{"cmap_label" : "Madgraph/Powheg"}
)
ax_ratio.text(
    0.98, 0.90, "$t\\bar{t}$\nparton level",
    fontsize=rcParams["font.size"], horizontalalignment="right", verticalalignment="bottom", transform=ax_ratio.transAxes, color="w",
)
hep.cms.label(ax=ax_ratio, data=False, year="Run 2")
figname_ratio = os.path.join(outdir, f"{base_jobid}_LO_to_POWHEG_Ratio_xsec_{tune_var}_rebinned_KFactor")
fig_ratio.savefig(figname_ratio)
print(f"{figname_ratio} written")

if args.save_ratios:
    #set_trace()
    LO_to_Powheg_KFactor_wts = np.where(LO_to_Powheg_KFactor == 0., 1, LO_to_Powheg_KFactor)
    save_dict["Rebinned_LO_to_Powheg_KFactor"] = dense_lookup(LO_to_Powheg_KFactor_wts, (masked_rebinned_mtt_binning, rebinned_ctstar_binning))
    ratios_fname = os.path.join(proj_dir, "EWK_Uncs", f"LO_to_Powheg_Ratios_{base_jobid}_CombineYears.coffea")# if args.combine_years else f"NNLO_to_Tune_Ratios_{base_jobid}_IndivYears.coffea")
    save(save_dict, ratios_fname)
    print(f"\n{ratios_fname} written")
