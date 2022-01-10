#!/usr/bin/env python

# matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend("agg")
from matplotlib import rcParams
rcParams["font.size"] = 20
rcParams["savefig.format"] = "png"
rcParams["savefig.bbox"] = "tight"

import Utilities.Plotter as Plotter
from coffea.hist import plot
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

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--combine_years", action="store_true", help="Combine all eras together.")
parser.add_argument("--save_ratios", action="store_true", help="Save NNLO/tune ratios for all years (default is to not save them).")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
base_jobid = os.environ["base_jobid"]
jobid = os.environ["jobid"]
analyzer = "make_nnlo_dists"

f_ext = "TOT.coffea"
outdir = os.path.join(proj_dir, "plots", base_jobid, analyzer)
if not os.path.isdir(outdir):
    os.makedirs(outdir)
    
style_dict = {
    "2016APV" : ("2016 pre-VFP", "#984ea3"), # purple
    "2016" : ("2016 post-VFP", "#377eb8"), ## blue
    "2017" : ("2017", "#e41a1c"), ## red
    "2018" : ("2018", "#4daf4a"), ## green
}

tune_var = "mtt_vs_top_ctstar"
save_dict = {} if args.combine_years else {year: {tune_var: {}} for year in style_dict.keys()}


    ## get values from NNLO root file
nnlo_fname = "MATRIX_ttmVStheta.root" # "xsec_central" dist has only statistical uncs
nnlo_file = convert_histo_root_file(os.path.join(proj_dir, "NNLO_files", nnlo_fname))
#set_trace()
nnlo_var = "xsec_central"
#nnlo_dict = Plotter.root_converters_dict_to_hist(nnlo_file, vars=[nnlo_var],
nnlo_dict = Plotter.root_converters_dict_to_hist(nnlo_file, vars=["xsec_central", "xsec_down", "xsec_up"],
    sparse_axes_list=[{"name": "dataset", "label" : "Event Process", "fill" : "nnlo"}],
    dense_axes_list=[{"name" : "ctstar", "idx" : 0}, {"name": "mtt", "idx" : 1}],
)
nnlo_histo = nnlo_dict[nnlo_var].integrate("dataset")
   # normalized
nnlo_normed_histo = nnlo_histo.copy()
nnlo_normed_histo.scale(1./nnlo_normed_histo.values(overflow="all")[()].sum())

nnlo_down_normed_histo = nnlo_dict["xsec_down"].integrate("dataset").copy()
nnlo_down_normed_histo.scale(1./nnlo_down_normed_histo.values(overflow="all")[()].sum())

nnlo_up_normed_histo = nnlo_dict["xsec_up"].integrate("dataset").copy()
nnlo_up_normed_histo.scale(1./nnlo_up_normed_histo.values(overflow="all")[()].sum())

mtt_binning = np.around(nnlo_dict[nnlo_var].axis("mtt").edges(), decimals=0)
ctstar_binning = np.around(nnlo_dict[nnlo_var].axis("ctstar").edges(), decimals=2)
#set_trace()
vlines = [(len(mtt_binning)-1)*ybin for ybin in range(1, len(ctstar_binning)-1)]
png_ext = "StatUncs"
nnlo_leg = "(Stat.)"

axis_labels_dict = {
    "mtt" : "$m_{t\\bar{t}}$",
    "ctstar" : "cos($\\theta^{*}_{t}$)"
}

        # histo plotting params
xtitle, ytitle = axis_labels_dict[nnlo_histo.dense_axes()[0].name], axis_labels_dict[nnlo_histo.dense_axes()[1].name]
ztitle = "$\dfrac{d^{2} \\sigma}{d m_{t\\bar{t}} d cos(\\theta^{*}_{t})}$"
opts = {"cmap_label" : ztitle}
norm_opts = {"cmap_label" : "$\dfrac{1}{\\sigma}$%s [$GeV^{-1}$]" % ztitle}
x_lims = (np.min(nnlo_histo.dense_axes()[0].edges()), np.max(nnlo_histo.dense_axes()[0].edges()))
y_lims = (np.min(nnlo_histo.dense_axes()[1].edges()), np.max(nnlo_histo.dense_axes()[1].edges()))

#set_trace()

if args.save_ratios:
    from Utilities import HistExport
    import uproot3
    outrname = os.path.join(proj_dir, "NNLO_files", f"NNLO_to_Tune_Ratios_{base_jobid}_CombineYears.root" if args.combine_years else f"NNLO_to_Tune_Ratios_{base_jobid}_IndivYears.root")
    upfout = uproot3.recreate(outrname, compression=uproot3.ZLIB(4)) if os.path.isfile(outrname) else uproot3.create(outrname)


data_lumi_dict = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())
lumi_correction = load(os.path.join(proj_dir, "Corrections", jobid, "MC_LumiWeights.coffea"))

if args.combine_years:
    years_to_run = ["2016APV", "2016", "2017", "2018"]
    histos_dict = {}
    #set_trace()
    for year in years_to_run:
        input_dir = os.path.join(proj_dir, "results", f"{year}_{base_jobid}", analyzer)
        fnames = sorted(["%s/%s" % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)])
        hdict = plt_tools.add_coffea_files(fnames) if len(fnames) > 1 else load(fnames[0])
    
        tune_histo = hdict[tune_var]
        scale_dict = {tt : lumi_correction[year]["Muons"][tt] for tt in ["ttJetsSL", "ttJetsHad", "ttJetsDiLep"]}
        tune_histo.scale(scale_dict, axis="dataset")
    
            # orig
        tune_histo = hdict[tune_var].integrate("dataset")
            # rebin
        tune_histo = tune_histo.rebin("mtt", hist.Bin("mtt", "mtt", mtt_binning))
        tune_histo = tune_histo.rebin("ctstar", hist.Bin("ctstar", "ctstar", ctstar_binning))

        histos_dict[year] = tune_histo

    tune_histo = histos_dict["2016APV"].add(histos_dict["2016"]).add(histos_dict["2017"]).add(histos_dict["2018"])
    tune_histo.scale(1./data_lumi_dict["TOT"]["Muons"])
    #set_trace()

        # save integral to make normalized hist
    tune_integral = tune_histo.values(overflow="all")[()].sum()
    
        # make normalized hist
    tune_normed_histo = tune_histo.copy()
    tune_normed_histo.scale(1./tune_integral)
    
        # plot yield 2d version of MC hist for each year
    fig, ax = plt.subplots()
    fig.subplots_adjust(hspace=.07)
    
    tune_vals = tune_histo.values()[()]        
    tune_vals = np.where(abs(tune_vals) < 1e-10, 0, tune_vals) # set values that are less than 1e-10 to 0
    Plotter.plot_2d_norm(tune_histo, xaxis_name=tune_histo.dense_axes()[0].name, yaxis_name=tune_histo.dense_axes()[1].name,
        values=np.ma.masked_where(tune_vals <= 0.0, tune_vals), # mask nonzero probabilities for plotting
        xlimits=y_lims, ylimits=x_lims, xlabel=axis_labels_dict[tune_histo.dense_axes()[0].name], ylabel=axis_labels_dict[tune_histo.dense_axes()[1].name],
        ax=ax, **opts)
    
    #set_trace()
    ax.text(
        0.98, 0.90, "$t\\bart$\nparton level",
        fontsize=rcParams["font.size"], horizontalalignment="right", verticalalignment="bottom", transform=ax.transAxes, color="w",
    )
    hep.cms.label(ax=ax, data=False, year="All Eras")
    
    figname = os.path.join(outdir, f"AllYears_{base_jobid}_ttJets_xsec_{tune_var}_{png_ext}")
    fig.savefig(figname)
    print(f"{figname} written")
    
    
        # plot normalized 2d version of MC hists for each year
    fig_norm, ax_norm = plt.subplots()
    fig_norm.subplots_adjust(hspace=.07)
    
    tune_norm_vals = tune_normed_histo.values()[()]
    tune_norm_vals = np.where(abs(tune_norm_vals) < 1e-10, 0, tune_norm_vals) # set values that are less than 1e-10 to 0
    Plotter.plot_2d_norm(tune_normed_histo, xaxis_name=tune_normed_histo.dense_axes()[0].name, yaxis_name=tune_normed_histo.dense_axes()[1].name,
        values=np.ma.masked_where(tune_norm_vals <= 0.0, tune_norm_vals), # mask nonzero probabilities for plotting
        xlimits=y_lims, ylimits=x_lims, xlabel=axis_labels_dict[tune_normed_histo.dense_axes()[0].name], ylabel=axis_labels_dict[tune_normed_histo.dense_axes()[1].name],
        ax=ax_norm, **norm_opts)
    
    ax_norm.text(
        0.98, 0.90, "$t\\bart$\nparton level",
        fontsize=rcParams["font.size"], horizontalalignment="right", verticalalignment="bottom", transform=ax_norm.transAxes, color="w",
    )
    hep.cms.label(ax=ax_norm, data=False, year="All Eras")
    
    figname_norm = os.path.join(outdir, f"AllYears_{base_jobid}_ttJets_xsec_{tune_var}_{png_ext}_Norm")
    fig_norm.savefig(figname_norm)
    print(f"{figname_norm} written")
    
    
        # plot NNLO/year (weights) each year
    fig_ratio, ax_ratio = plt.subplots()
    fig_ratio.subplots_adjust(hspace=.07)
    
    #set_trace()
    normed_ratio_vals = nnlo_normed_histo.values()[()].T/tune_norm_vals
    normed_ratio_vals = np.where(normed_ratio_vals <= 0, 1, normed_ratio_vals)
    normed_ratio_vals = np.where(normed_ratio_vals == np.inf, 1, normed_ratio_vals)
    Plotter.plot_2d_norm(tune_normed_histo, xaxis_name=tune_normed_histo.dense_axes()[0].name, yaxis_name=tune_normed_histo.dense_axes()[1].name,
        values=np.ma.masked_where(normed_ratio_vals <= 0.0, normed_ratio_vals), # mask nonzero probabilities for plotting
        xlimits=y_lims, ylimits=x_lims, xlabel=axis_labels_dict[tune_normed_histo.dense_axes()[0].name], ylabel=axis_labels_dict[tune_normed_histo.dense_axes()[1].name],
        ax=ax_ratio, **{"cmap_label" : "NNLO/Tune"})
    
    ax_ratio.text(
        0.98, 0.90, "$t\\bart$\nparton level",
        fontsize=rcParams["font.size"], horizontalalignment="right", verticalalignment="bottom", transform=ax_ratio.transAxes, color="w",
    )
    hep.cms.label(ax=ax_ratio, data=False, year="All Eras")
    
    figname_ratio = os.path.join(outdir, f"AllYears_{base_jobid}_NNLO_to_Tune_Ratio_xsec_{tune_var}")
    fig_ratio.savefig(figname_ratio)
    print(f"{figname_ratio} written")
    
    if args.save_ratios:
        #set_trace()
        save_dict[tune_var] = dense_lookup(normed_ratio_vals, (mtt_binning, ctstar_binning))
            # central
        tmp_central_hist = Plotter.np_array_TO_hist(sumw=normed_ratio_vals, sumw2=np.zeros(normed_ratio_vals.shape), hist_template=tune_normed_histo)
        upfout[f"MATRIX_xsec_central_over_AllYears_{base_jobid}_Normalized_Ratios"] = HistExport.export2d(tmp_central_hist)
        #    # down
        #normed_down_ratio_vals = nnlo_down_normed_histo.values()[()].T/tune_norm_vals
        #normed_down_ratio_vals = np.where(normed_down_ratio_vals <= 0, 1, normed_down_ratio_vals)
        #normed_down_ratio_vals = np.where(normed_down_ratio_vals == np.inf, 1, normed_down_ratio_vals)
        #tmp_down_hist = Plotter.np_array_TO_hist(sumw=normed_down_ratio_vals, sumw2=np.zeros(normed_down_ratio_vals.shape), hist_template=tune_normed_histo)
        #upfout[f"MATRIX_xsec_down_to_AllYeras_Normalized_Ratios_{base_jobid}"] = HistExport.export2d(tmp_down_hist)
        #    # up
        #normed_up_ratio_vals = nnlo_up_normed_histo.values()[()].T/tune_norm_vals
        #normed_up_ratio_vals = np.where(normed_up_ratio_vals <= 0, 1, normed_up_ratio_vals)
        #normed_up_ratio_vals = np.where(normed_up_ratio_vals == np.inf, 1, normed_up_ratio_vals)
        #tmp_up_hist = Plotter.np_array_TO_hist(sumw=normed_up_ratio_vals, sumw2=np.zeros(normed_up_ratio_vals.shape), hist_template=tune_normed_histo)
        #upfout[f"MATRIX_xsec_up_to_AllYears_Normalized_Ratios_{base_jobid}"] = HistExport.export2d(tmp_up_hist)
    
            # tune histo
        upfout[f"AllYears_{base_jobid}_xsec"] = HistExport.export2d(tune_histo)

else:
    #years_to_run = ["2018"]
    #years_to_run = ["2017"]
    years_to_run = ["2016APV", "2016", "2017", "2018"]    
    for year in years_to_run:
        input_dir = os.path.join(proj_dir, "results", f"{year}_{base_jobid}", analyzer)
        fnames = sorted(["%s/%s" % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)])
        hdict = plt_tools.add_coffea_files(fnames) if len(fnames) > 1 else load(fnames[0])
    
        tune_histo = hdict[tune_var]
        scale_dict = {tt : lumi_correction[year]["Muons"][tt]/data_lumi_dict[year]["Muons"] for tt in ["ttJetsSL", "ttJetsHad", "ttJetsDiLep"]}
        tune_histo.scale(scale_dict, axis="dataset")
    
            # orig
        tune_histo = hdict[tune_var].integrate("dataset")
    
        #set_trace()
            # rebin
        tune_histo = tune_histo.rebin("mtt", hist.Bin("mtt", "mtt", mtt_binning))
        tune_histo = tune_histo.rebin("ctstar", hist.Bin("ctstar", "ctstar", ctstar_binning))
    
            # save integral to make normalized hist
        tune_integral = tune_histo.values(overflow="all")[()].sum()
    
            # make normalized hist
        tune_normed_histo = tune_histo.copy()
        tune_normed_histo.scale(1./tune_integral)
    
            # plot yield 2d version of MC hist for each year
        fig, ax = plt.subplots()
        fig.subplots_adjust(hspace=.07)
    
        tune_vals = tune_histo.values()[()]        
        tune_vals = np.where(abs(tune_vals) < 1e-10, 0, tune_vals) # set values that are less than 1e-10 to 0
        Plotter.plot_2d_norm(tune_histo, xaxis_name=tune_histo.dense_axes()[0].name, yaxis_name=tune_histo.dense_axes()[1].name,
            values=np.ma.masked_where(tune_vals <= 0.0, tune_vals), # mask nonzero probabilities for plotting
            xlimits=y_lims, ylimits=x_lims, xlabel=axis_labels_dict[tune_histo.dense_axes()[0].name], ylabel=axis_labels_dict[tune_histo.dense_axes()[1].name],
            ax=ax, **opts)
        
        #set_trace()
        ax.text(
            0.98, 0.90, "$t\\bart$\nparton level",
            fontsize=rcParams["font.size"], horizontalalignment="right", verticalalignment="bottom", transform=ax.transAxes, color="w",
        )
        hep.cms.label(ax=ax, data=False, year=year)
        
        figname = os.path.join(outdir, f"{year}_{base_jobid}_ttJets_xsec_{tune_var}_{png_ext}")
        fig.savefig(figname)
        print(f"{figname} written")
    
    
            # plot normalized 2d version of MC hists for each year
        fig_norm, ax_norm = plt.subplots()
        fig_norm.subplots_adjust(hspace=.07)
        
        tune_norm_vals = tune_normed_histo.values()[()]
        tune_norm_vals = np.where(abs(tune_norm_vals) < 1e-10, 0, tune_norm_vals) # set values that are less than 1e-10 to 0
        Plotter.plot_2d_norm(tune_normed_histo, xaxis_name=tune_normed_histo.dense_axes()[0].name, yaxis_name=tune_normed_histo.dense_axes()[1].name,
            values=np.ma.masked_where(tune_norm_vals <= 0.0, tune_norm_vals), # mask nonzero probabilities for plotting
            xlimits=y_lims, ylimits=x_lims, xlabel=axis_labels_dict[tune_normed_histo.dense_axes()[0].name], ylabel=axis_labels_dict[tune_normed_histo.dense_axes()[1].name],
            ax=ax_norm, **norm_opts)
        
        ax_norm.text(
            0.98, 0.90, "$t\\bart$\nparton level",
            fontsize=rcParams["font.size"], horizontalalignment="right", verticalalignment="bottom", transform=ax_norm.transAxes, color="w",
        )
        hep.cms.label(ax=ax_norm, data=False, year=year)
        
        figname_norm = os.path.join(outdir, f"{year}_{base_jobid}_ttJets_xsec_{tune_var}_{png_ext}_Norm")
        fig_norm.savefig(figname_norm)
        print(f"{figname_norm} written")
    
    
            # plot NNLO/year (weights) each year
        fig_ratio, ax_ratio = plt.subplots()
        fig_ratio.subplots_adjust(hspace=.07)
    
        #set_trace()
        normed_ratio_vals = nnlo_normed_histo.values()[()].T/tune_norm_vals
        normed_ratio_vals = np.where(normed_ratio_vals <= 0, 1, normed_ratio_vals)
        normed_ratio_vals = np.where(normed_ratio_vals == np.inf, 1, normed_ratio_vals)
        Plotter.plot_2d_norm(tune_normed_histo, xaxis_name=tune_normed_histo.dense_axes()[0].name, yaxis_name=tune_normed_histo.dense_axes()[1].name,
            values=np.ma.masked_where(normed_ratio_vals <= 0.0, normed_ratio_vals), # mask nonzero probabilities for plotting
            xlimits=y_lims, ylimits=x_lims, xlabel=axis_labels_dict[tune_normed_histo.dense_axes()[0].name], ylabel=axis_labels_dict[tune_normed_histo.dense_axes()[1].name],
            ax=ax_ratio, **{"cmap_label" : "NNLO/Tune"})
        
        ax_ratio.text(
            0.98, 0.90, "$t\\bart$\nparton level",
            fontsize=rcParams["font.size"], horizontalalignment="right", verticalalignment="bottom", transform=ax_ratio.transAxes, color="w",
        )
        hep.cms.label(ax=ax_ratio, data=False, year=year)
        
        figname_ratio = os.path.join(outdir, f"{year}_{base_jobid}_NNLO_to_Tune_Ratio_xsec_{tune_var}")
        fig_ratio.savefig(figname_ratio)
        print(f"{figname_ratio} written")
    
        if args.save_ratios:
            #set_trace()
            save_dict[year][tune_var] = dense_lookup(normed_ratio_vals, (mtt_binning, ctstar_binning))
                # central
            tmp_central_hist = Plotter.np_array_TO_hist(sumw=normed_ratio_vals, sumw2=np.zeros(normed_ratio_vals.shape), hist_template=tune_normed_histo)
            upfout[f"MATRIX_xsec_central_to_{year}_Normalized_Ratios_{base_jobid}"] = HistExport.export2d(tmp_central_hist)
            #    # down
            #normed_down_ratio_vals = nnlo_down_normed_histo.values()[()].T/tune_norm_vals
            #normed_down_ratio_vals = np.where(normed_down_ratio_vals <= 0, 1, normed_down_ratio_vals)
            #normed_down_ratio_vals = np.where(normed_down_ratio_vals == np.inf, 1, normed_down_ratio_vals)
            #tmp_down_hist = Plotter.np_array_TO_hist(sumw=normed_down_ratio_vals, sumw2=np.zeros(normed_down_ratio_vals.shape), hist_template=tune_normed_histo)
            #upfout[f"MATRIX_xsec_down_to_{year}_Normalized_Ratios_{base_jobid}"] = HistExport.export2d(tmp_down_hist)
            #    # up
            #normed_up_ratio_vals = nnlo_up_normed_histo.values()[()].T/tune_norm_vals
            #normed_up_ratio_vals = np.where(normed_up_ratio_vals <= 0, 1, normed_up_ratio_vals)
            #normed_up_ratio_vals = np.where(normed_up_ratio_vals == np.inf, 1, normed_up_ratio_vals)
            #tmp_up_hist = Plotter.np_array_TO_hist(sumw=normed_up_ratio_vals, sumw2=np.zeros(normed_up_ratio_vals.shape), hist_template=tune_normed_histo)
            #upfout[f"MATRIX_xsec_up_to_{year}_Normalized_Ratios_{base_jobid}"] = HistExport.export2d(tmp_up_hist)
    
                # tune histo
            upfout[f"{year}_{base_jobid}_xsec"] = HistExport.export2d(tune_histo)
            #save_dict.update({year: dense_lookup(normed_ratio_vals, (mtt_binning, ctstar_binning))})

# save weights
if args.save_ratios:
    ratios_fname = os.path.join(proj_dir, "NNLO_files", f"NNLO_to_Tune_Ratios_{base_jobid}_CombineYears.coffea" if args.combine_years else f"NNLO_to_Tune_Ratios_{base_jobid}_IndivYears.coffea")
    save(save_dict, ratios_fname)
    print(f"\n{ratios_fname} written")

    upfout["MATRIX_xec_central"] = HistExport.export2d(nnlo_histo)
    upfout["MATRIX_xec_down"] = HistExport.export2d(nnlo_dict["xsec_down"].integrate("dataset"))
    upfout["MATRIX_xec_up"] = HistExport.export2d(nnlo_dict["xsec_up"].integrate("dataset"))
    upfout.close()
    print(f"{outrname} written") 
