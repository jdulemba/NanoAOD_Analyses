# matplotlib
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
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
import Utilities.plot_tools as plt_tools
import Utilities.prettyjson as prettyjson
from coffea import hist
import numpy as np
import fnmatch
import Utilities.Plotter as Plotter
from Utilities.styles import styles as hstyles
import Utilities.final_analysis_binning as final_binning
from scipy.signal import peak_widths
import Utilities.common_features as cfeatures

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2018"], help="What year is the ntuple from.")
#parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("lepton", choices=["Electron", "Muon"], help="Choose which lepton to make plots for")
parser.add_argument("plots", choices=["RECO", "GEN", "RESO", "REL_RESO", "All"], help="Choose which hists to make plots for")
parser.add_argument("binning", choices=["Final", "Default"], help="Choose which mtt binning to use for histograms.")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
base_jobid = os.environ["base_jobid"]
jobid = os.environ["jobid"]
analyzer = "binning_check"
plot_outdir = os.environ["plots_dir"]
eos_dir = os.environ["eos_dir"]

input_dir = os.path.join(eos_dir, "results", f"{args.year}_{jobid}", analyzer)
fnames = [os.path.join(input_dir, fname) for fname in fnmatch.filter(os.listdir(input_dir), "*TOT.coffea")]
hdict = load(fnames[0])

outdir = os.path.join(plot_outdir, f"{args.year}_{jobid}", analyzer)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

rewt_style_dict = {
    "TT" : {"label" : "SM $\mathrm{t\\bar t}_{\ell j}$", "facecolor" : "none", "hatch" : "//", "edgecolor" : "r"},
    "Res" : {"facecolor" : "#377eb8", "edgecolor" : "#377eb8"}, ## blue
    "Int_neg" : {"color" : "k", "linestyle" : "-"},
    "Int_pos" : {"color" : "k", "linestyle" : "--"},
}


reco_variables = {
    "RECO_mtt_vs_tlep_ctstar" : (cfeatures.variable_names_to_labels["mtt"], cfeatures.variable_names_to_labels["tlep_ctstar"], final_binning.mtt_binning if (args.binning == "Final") else 1, 1, (200., 2000.),  (-1., 1.)),
    "RECO_mtt_vs_tlep_ctstar_abs" : (cfeatures.variable_names_to_labels["mtt"], cfeatures.variable_names_to_labels["tlep_ctstar_abs"], final_binning.mtt_binning if (args.binning == "Final") else 1, final_binning.ctstar_abs_binning if (args.binning == "Final") else 1, (200., 2000.),  (0., 1.)),
    "RECO_mtt_vs_topRapidity" : (cfeatures.variable_names_to_labels["mtt"], cfeatures.variable_names_to_labels["rapidity_top"], 1, 1, (200., 2000.), (-3., 3.)),
    "RECO_mtt_vs_deltaYtt" : (cfeatures.variable_names_to_labels["mtt"], "$\\Delta y_{t\\bar{t}}$", 1, 1, (200., 2000.), (-5., 5.)),
    "RECO_mtt_vs_deltaYtt_abs" : (cfeatures.variable_names_to_labels["mtt"], "|$\\Delta y_{t\\bar{t}}$|", 1, 1, (200., 2000.), (0., 5.)),
}
gen_variables = {
    "GEN_mtt_vs_tlep_ctstar" : (cfeatures.variable_names_to_labels["mtt"], cfeatures.variable_names_to_labels["tlep_ctstar"], final_binning.mtt_binning if (args.binning == "Final") else 1, 1, (200., 2000.),  (-1., 1.)),
    "GEN_mtt_vs_tlep_ctstar_abs" : (cfeatures.variable_names_to_labels["mtt"], cfeatures.variable_names_to_labels["tlep_ctstar_abs"], final_binning.mtt_binning if (args.binning == "Final") else 1, final_binning.ctstar_abs_binning if (args.binning == "Final") else 1, (200., 2000.),  (0., 1.)),
    "GEN_mtt_vs_topRapidity" : (cfeatures.variable_names_to_labels["mtt"], cfeatures.variable_names_to_labels["rapidity_top"], 1, 1, (200., 2000.), (-3., 3.)),
    "GEN_mtt_vs_deltaYtt" : (cfeatures.variable_names_to_labels["mtt"], "$\\Delta y_{t\\bar{t}}$", 1, 1, (200., 2000.), (-5., 5.)),
    "GEN_mtt_vs_deltaYtt_abs" : (cfeatures.variable_names_to_labels["mtt"], "|$\\Delta y_{t\\bar{t}}$|", 1, 1, (200., 2000.), (0., 5.)),
}
reso_variables = {
    "RESO_mtt_vs_mtt" : (cfeatures.variable_names_to_labels["mtt"], f"Gen-Reco {cfeatures.variable_names_to_labels['mtt']}", final_binning.mtt_binning if (args.binning == "Final") else 1, 1, (200., 2000.),  (-500., 500.)),
    "RESO_ctstar_vs_mtt" : (cfeatures.variable_names_to_labels["mtt"], f"Gen-Reco {cfeatures.variable_names_to_labels['tlep_ctstar']}", final_binning.mtt_binning if (args.binning == "Final") else 1, 1, (200., 2000.),  (-2., 2.)),
    "RESO_ctstar_abs_vs_mtt" : (cfeatures.variable_names_to_labels["mtt"], f"Gen-Reco {cfeatures.variable_names_to_labels['tlep_ctstar_abs']}", final_binning.mtt_binning if (args.binning == "Final") else 1, 1, (200., 2000.),  (-1., 1.)),
    ##"RESO_mtt_vs_tlep_ctstar" : ("$m_{t\\bar{t}}$", "cos($\\theta^{*}_{t_{l}}$)", 1, 1, (-500., 500.),  (-2., 2.)),
    ##"RESO_mtt_vs_tlep_ctstar_abs" : ("$m_{t\\bar{t}}$", "|cos($\\theta^{*}_{t_{l}}$)|", 1, 1, (-500., 500.),  (-1., 1.)),
}
rel_reso_variables = {
    "REL_RESO_mtt_vs_mtt" : ("$m^{reco}_{t\\bar{t}}$ [GeV]", "($m^{reco}_{t\\bar{t}}$-$m^{parton}_{t\\bar{t}}$)/$m^{parton}_{t\\bar{t}}$", final_binning.mtt_binning if (args.binning == "Final") else 1, 1, (200., 2000.),  (-2., 2.)),
}


    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())[args.year]
lumi_correction = load(os.path.join(proj_dir, "Corrections", base_jobid, "MC_LumiWeights.coffea"))[args.year][f"{args.lepton}s"]
        # scale ttJets events, split by reconstruction type, by normal ttJets lumi correction
ttJets_permcats = ["*right", "*matchable", "*unmatchable", "*sl_tau", "*other"]
#set_trace()
names = [dataset for dataset in sorted(set([key[0] for key in hdict["RECO_mtt_vs_tlep_ctstar_abs"].values().keys()]))] # get dataset names in hists
ttJets_cats = [name for name in names if any([fnmatch.fnmatch(name, cat) for cat in ttJets_permcats])] # gets ttJets(_PS)_other, ...
if len(ttJets_cats) > 0:
    for tt_cat in ttJets_cats:
        ttJets_lumi_topo = "_".join(tt_cat.split("_")[:-2]) if "sl_tau" in tt_cat else "_".join(tt_cat.split("_")[:-1]) # gets ttJets[SL, Had, DiLep] or ttJets_PS
        ttJets_eff_lumi = lumi_correction[ttJets_lumi_topo]
        lumi_correction.update({tt_cat: ttJets_eff_lumi})

## make groups based on process
process = hist.Cat("process", "Process", sorting="placement")
process_cat = "dataset"
process_groups = plt_tools.make_dataset_groups(args.lepton, args.year, samples=names, bkgdict="dataset" if (args.plots == "REL_RESO") else "templates")
#set_trace()

    # scale and group hists by process
for hname in hdict.keys():
    if "cutflow" in hname: continue
    #set_trace()
    hdict[hname].scale(lumi_correction, axis="dataset") # scale hists
    hdict[hname] = hdict[hname].group(process_cat, process, process_groups) # group by process
    hdict[hname] = hdict[hname][:, :, args.lepton].integrate("leptype") # only pick out specified lepton



#set_trace()
    ## make reco plots
if (args.plots == "RECO") or (args.plots == "All"):
    for hname in reco_variables.keys():
        if hname not in hdict.keys():
            raise ValueError(f"{hname} not found in file")
        #set_trace()
        histo = hdict[hname].copy()
        #set_trace()
    
        if histo.dense_dim() == 2:
            xtitle, ytitle, xrebinning, yrebinning, x_lims, y_lims = reco_variables[hname]
    
            xaxis_name = histo.dense_axes()[0].name
            yaxis_name = histo.dense_axes()[1].name
                ## rebin x axis
            if isinstance(xrebinning, np.ndarray):
                new_xbins = hist.Bin(xaxis_name, xaxis_name, xrebinning)
            elif isinstance(xrebinning, float) or isinstance(xrebinning, int):
                new_xbins = xrebinning
            histo = histo.rebin(xaxis_name, new_xbins)
                ## rebin y axis
            if isinstance(yrebinning, np.ndarray):
                new_ybins = hist.Bin(yaxis_name, yaxis_name, yrebinning)
            elif isinstance(yrebinning, float) or isinstance(yrebinning, int):
                new_ybins = yrebinning
            histo = histo.rebin(yaxis_name, new_ybins)

            ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
            for jmult in sorted(set([key[1] for key in histo.values().keys()])):
                print(*[jmult, hname], sep=", ")
                pltdir = os.path.join(outdir, args.lepton, jmult, "RECO")
                if not os.path.isdir(pltdir):
                    os.makedirs(pltdir)
    
                hslice = histo[:, jmult].integrate("jmult")
    
                    # make 1D projection along dense axes
                for dax in range(2):
                    if dax == 0:
                        xlabel = ytitle
                        xlimits = y_lims
                        figname = os.path.join(pltdir, "_".join([jmult, args.lepton, "FinalBinning", "RECO", yaxis_name, "Norm"])) if (args.binning == "Final") \
                            else os.path.join(pltdir, "_".join([jmult, args.lepton, "RECO", yaxis_name, "Norm"]))
                        figname_ratio = os.path.join(pltdir, "_".join([jmult, args.lepton, "FinalBinning", "RECO", yaxis_name, "REStoTTbar"])) if (args.binning == "Final") \
                            else os.path.join(pltdir, "_".join([jmult, args.lepton, "RECO", yaxis_name, "REStoTTbar"]))
                    else:
                        xlabel = xtitle
                        xlimits = x_lims
                        figname = os.path.join(pltdir, "_".join([jmult, args.lepton, "FinalBinning", "RECO", xaxis_name, "Norm"])) if (args.binning == "Final") \
                            else os.path.join(pltdir, "_".join([jmult, args.lepton, "RECO", xaxis_name, "Norm"]))
                        figname_ratio = os.path.join(pltdir, "_".join([jmult, args.lepton, "FinalBinning", "RECO", xaxis_name, "REStoTTbar"])) if (args.binning == "Final") \
                            else os.path.join(pltdir, "_".join([jmult, args.lepton, "RECO", xaxis_name, "REStoTTbar"]))
    
                    hproj = hslice.integrate(hslice.dense_axes()[dax].name)
                    bin_edges = hproj.dense_axes()[0].edges()
    
                        # normalized plots
                    fig, ax = plt.subplots()
                    fig.subplots_adjust(hspace=.07)

                    #set_trace()    
                    for proc in [("AtoTTJetsSL_M500_W10p0_Res",), ("AtoTTJetsSL_M500_W10p0_Int_neg",), ("AtoTTJetsSL_M500_W10p0_Int_pos",), ("TT",)]:
                        norm_vals = np.r_[hproj.values()[proc], hproj.values()[proc][-1]]/np.sum(hproj.values(overflow="all")[proc])
                        if proc[0] == "TT":
                            tt_vals = norm_vals
                            ax.fill_between(bin_edges, norm_vals, step="post", **rewt_style_dict["TT"])
                        elif "Res" in proc[0]:
                            res_vals = norm_vals
                            rewt_style_dict["Res"]["label"] = plt_tools.get_label(proc[0], hstyles)
                            ax.fill_between(bin_edges, norm_vals, step="post", **rewt_style_dict["Res"])
                        elif "Int_neg" in proc[0]:
                            rewt_style_dict["Int_neg"]["label"] = plt_tools.get_label(proc[0], hstyles)
                            ax.step(bin_edges, norm_vals, where="post", **rewt_style_dict["Int_neg"])
                        elif "Int_pos" in proc[0]:
                            rewt_style_dict["Int_pos"]["label"] = plt_tools.get_label(proc[0], hstyles)
                            ax.step(bin_edges, norm_vals, where="post", **rewt_style_dict["Int_pos"])
                            
                    ax.legend(loc="upper right")
                    ax.autoscale()
                    ax.set_ylim(0, ax.get_ylim()[1]*1.15)
                    ax.set_xlim(xlimits)
                    ax.set_xlabel(xlabel)
                    ax.set_ylabel("Probability Density [a.u.]")
    
                        # add lepton/jet multiplicity label
                    ax.text(
                        0.02, 0.88, f"{cfeatures.channel_labels[f'{args.lepton}_{jmult}']}\ndetector level",
                        fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                    )
                    hep.cms.label(ax=ax, data=False, label="Preliminary", year=cfeatures.year_labels[args.year])
    
                    #set_trace()
                    fig.savefig(figname)
                    print(f"{figname} written")
                    plt.close(fig)
    
                    #set_trace()
                        # plot resonant/ttbar
                    fig_ratio, ax_ratio = plt.subplots()
                    fig_ratio.subplots_adjust(hspace=.07)
    
                    ax_ratio.step(bin_edges, res_vals/tt_vals, where="post", **{"color" : "k"})
                            
                    #ax_ratio.legend(loc="upper right")
                    ax_ratio.autoscale()
                    ax_ratio.set_ylim(0, ax_ratio.get_ylim()[1]*1.15)
                    ax_ratio.set_xlim(xlimits)
                    ax_ratio.set_xlabel(xlabel)
                    ax_ratio.set_ylabel("Resonant/SM $t\\bar{t}$ [a.u.]")
                    ax_ratio.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
    
                        # add lepton/jet multiplicity label
                    ax_ratio.text(
                        0.02, 0.88, f"{cfeatures.channel_labels[f'{args.lepton}_{jmult}']}\ndetector level",
                        fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax_ratio.transAxes
                    )
                    hep.cms.label(ax=ax_ratio, data=False, label="Preliminary", year=cfeatures.year_labels[args.year])
    
                    #set_trace()
                    fig_ratio.savefig(figname_ratio)
                    print(f"{figname_ratio} written")
                    plt.close(fig_ratio)


                    #set_trace()
                        # plott TT yields with final binning overlaid on top
                    if (args.binning == "Default") and (hslice.dense_axes()[dax].name != "mtt"):
                        #set_trace()
                            # plots
                        fig, ax = plt.subplots()
                        fig.subplots_adjust(hspace=.07)
   
                        bin_edges = hproj.dense_axes()[0].edges()
                        vals = hproj.values()[("TT",)]
                        ax.fill_between(bin_edges, np.r_[vals, vals[-1]], step="post", **{"label": "$\\mathrm{t\\bar t}_{\\ell j}$", "facecolor": "r", "edgecolor": "k"})

                        ax.legend(loc="upper right")
                        ax.autoscale()
                        #set_trace()
                        ymax = ax.get_ylim()[1]*1.01
                        [ax.vlines(vline, ymin=0., ymax=ymax, color="k", linestyle="--") for vline in final_binning.mtt_binning]

                        ax.set_xlabel(xlabel)
                        ax.set_ylabel("Events")
                        ax.set_xlim((final_binning.mtt_binning[0], final_binning.mtt_binning[-1]))
                        ax.set_ylim(0, ax.get_ylim()[1]*1.15)
                                
    
                            # add lepton/jet multiplicity label
                        ax.text(
                            0.02, 0.88, f"{cfeatures.channel_labels[f'{args.lepton}_{jmult}']}\ndetector level",
                            fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                        )
                        hep.cms.label(ax=ax, data=False, label="Preliminary", year=cfeatures.year_labels[args.year], lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))
    
                        #set_trace()
                        figname = os.path.join(pltdir, "_".join([jmult, args.lepton, "RECO", xaxis_name]))
                        fig.savefig(figname)
                        print(f"{figname} written")
                        plt.close(fig)
    

    ## make gen plots
if (args.plots == "GEN") or (args.plots == "All"):
    for hname in gen_variables.keys():
        if hname not in hdict.keys():
            raise ValueError(f"{hname} not found in file")
        #set_trace()
        histo = hdict[hname].copy()
        #set_trace()
    
        if histo.dense_dim() == 2:
            xtitle, ytitle, xrebinning, yrebinning, x_lims, y_lims = gen_variables[hname]
    
            xaxis_name = histo.dense_axes()[0].name
            yaxis_name = histo.dense_axes()[1].name
                ## rebin x axis
            if isinstance(xrebinning, np.ndarray):
                new_xbins = hist.Bin(xaxis_name, xaxis_name, xrebinning)
            elif isinstance(xrebinning, float) or isinstance(xrebinning, int):
                new_xbins = xrebinning
            histo = histo.rebin(xaxis_name, new_xbins)
                ## rebin y axis
            if isinstance(yrebinning, np.ndarray):
                new_ybins = hist.Bin(yaxis_name, yaxis_name, yrebinning)
            elif isinstance(yrebinning, float) or isinstance(yrebinning, int):
                new_ybins = yrebinning
            histo = histo.rebin(yaxis_name, new_ybins)
    
            ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
            for jmult in sorted(set([key[1] for key in histo.values().keys()])):
                print(*[jmult, hname], sep=", ")
                pltdir = os.path.join(outdir, args.lepton, jmult, "GEN")
                if not os.path.isdir(pltdir):
                    os.makedirs(pltdir)

                hslice = histo[:, jmult].integrate("jmult").integrate("process") # only TT events filled at gen-level
   
                    # make 1D projection along dense axes
                for dax in range(2):
                    if dax == 0:
                        xlabel = ytitle
                        xlimits = y_lims
                        figname = os.path.join(pltdir, "_".join([jmult, args.lepton, "FinalBinning", "GEN", yaxis_name])) if (args.binning == "Final") \
                            else os.path.join(pltdir, "_".join([jmult, args.lepton, "GEN", yaxis_name]))
                    else:
                        xlabel = xtitle
                        xlimits = x_lims
                        figname = os.path.join(pltdir, "_".join([jmult, args.lepton, "FinalBinning", "GEN", xaxis_name])) if (args.binning == "Final") \
                            else os.path.join(pltdir, "_".join([jmult, args.lepton, "GEN", xaxis_name]))
    
                    hproj = hslice.integrate(hslice.dense_axes()[dax].name)
                    bin_edges = hproj.dense_axes()[0].edges()
                    vals = hproj.values()[()]
    
                        # plots
                    fig, ax = plt.subplots()
                    fig.subplots_adjust(hspace=.07)
   
                    Plotter.plot_1D(vals, bin_edges, xlimits=xlimits, xlabel=xlabel, ax=ax, label=rewt_style_dict["TT"]["label"])
                            
                    ax.legend(loc="upper right")
                    ax.set_ylim(0, ax.get_ylim()[1]*1.15)
    
                        # add lepton/jet multiplicity label
                    ax.text(
                        0.02, 0.88, f"{cfeatures.channel_labels[f'{args.lepton}_{jmult}']}\nparton level",
                        fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                    )
                    hep.cms.label(ax=ax, data=False, label="Preliminary", year=cfeatures.year_labels[args.year], lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))
    
                    #set_trace()
                    fig.savefig(figname)
                    print(f"{figname} written")
                    plt.close(fig)


    ## make resolution plots
if (args.plots == "RESO") or (args.plots == "All"):
    for hname in reso_variables.keys():
        if hname not in hdict.keys():
            raise ValueError(f"{hname} not found in file")
        #set_trace()
        histo = hdict[hname].copy()
        #set_trace()
    
        if histo.dense_dim() == 2:
            xtitle, ytitle, xrebinning, yrebinning, x_lims, y_lims = reso_variables[hname]
    
            xaxis_name = histo.dense_axes()[0].name
            yaxis_name = histo.dense_axes()[1].name
                ## rebin x axis
            if isinstance(xrebinning, np.ndarray):
                new_xbins = hist.Bin(xaxis_name, xaxis_name, xrebinning)
            elif isinstance(xrebinning, float) or isinstance(xrebinning, int):
                new_xbins = xrebinning
            histo = histo.rebin(xaxis_name, new_xbins)
                ## rebin y axis
            if isinstance(yrebinning, np.ndarray):
                new_ybins = hist.Bin(yaxis_name, yaxis_name, yrebinning)
            elif isinstance(yrebinning, float) or isinstance(yrebinning, int):
                new_ybins = yrebinning
            histo = histo.rebin(yaxis_name, new_ybins)
    
            ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
            for jmult in sorted(set([key[1] for key in histo.values().keys()])):
                print(*[jmult, hname], sep=", ")
                pltdir = os.path.join(outdir, args.lepton, jmult, "RESO")
                if not os.path.isdir(pltdir):
                    os.makedirs(pltdir)

                hslice = histo[:, jmult].integrate("jmult").integrate("process") # only TT filled at gen-level
    
                    # make 1D projection along dense axes
                for dax in range(2):
                    if hslice.dense_axes()[dax].name != "mtt": continue # skip normal mtt dist
                    if dax == 0:
                        xlabel = ytitle
                        xlimits = y_lims
                        figname = os.path.join(pltdir, "_".join([jmult, args.lepton, "FinalBinning", "RESO", yaxis_name])) if (args.binning == "Final") \
                            else os.path.join(pltdir, "_".join([jmult, args.lepton, "RESO", yaxis_name]))
                    else:
                        xlabel = xtitle
                        xlimits = x_lims
                        figname = os.path.join(pltdir, "_".join([jmult, args.lepton, "FinalBinning", "RESO", xaxis_name])) if (args.binning == "Final") \
                            else os.path.join(pltdir, "_".join([jmult, args.lepton, "RESO", xaxis_name]))
   
                    #set_trace() 
                    hproj = hslice.integrate(hslice.dense_axes()[dax].name)
                    bin_edges = hproj.dense_axes()[0].edges()
                    bin_centers = hproj.dense_axes()[0].centers()
                    vals = hproj.values()[()]

                    mean = np.average(bin_centers, weights=vals)
                    stdev = np.sqrt(np.average((bin_centers-mean)**2, weights=vals))                    
    
                        # plots
                    fig, ax = plt.subplots()
                    fig.subplots_adjust(hspace=.07)
  
                    #set_trace() 
                    Plotter.plot_1D(vals, bin_edges, xlimits=xlimits, xlabel=xlabel, ax=ax, label=rewt_style_dict["TT"]["label"])
                            
                    ax.legend(loc="upper right")
                    ax.set_ylim(0, ax.get_ylim()[1]*1.15)
    
                        # add lepton/jet multiplicity label
                    ax.text(
                        0.02, 0.84, "%s\nmean=%.4f\nstd=%.4f" % (cfeatures.channel_labels[f"{args.lepton}_{jmult}"], mean, stdev),
                        fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                    )
                    hep.cms.label(ax=ax, data=False, label="Preliminary", year=cfeatures.year_labels[args.year], lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))
    
                    #set_trace()
                    fig.savefig(figname)
                    print(f"{figname} written")
                    plt.close(fig)


    ## make relative resolution plots
if (args.plots == "REL_RESO") or (args.plots == "All"):
    #from hepstats.modeling import bayesian_blocks
    for hname in rel_reso_variables.keys():
        if hname not in hdict.keys():
            raise ValueError(f"{hname} not found in file")
        #set_trace()
        histo = hdict[hname].copy()
        #set_trace()
    
        if histo.dense_dim() == 2:
            xtitle, ytitle, xrebinning, yrebinning, x_lims, y_lims = rel_reso_variables[hname]
    
            xaxis_name = histo.dense_axes()[0].name
            yaxis_name = histo.dense_axes()[1].name
                ## rebin x axis
            if isinstance(xrebinning, np.ndarray):
                new_xbins = hist.Bin(xaxis_name, xaxis_name, xrebinning)
            elif isinstance(xrebinning, float) or isinstance(xrebinning, int):
                new_xbins = xrebinning
            histo = histo.rebin(xaxis_name, new_xbins)
                ## rebin y axis
            if isinstance(yrebinning, np.ndarray):
                new_ybins = hist.Bin(yaxis_name, yaxis_name, yrebinning)
            elif isinstance(yrebinning, float) or isinstance(yrebinning, int):
                new_ybins = yrebinning
            histo = histo.rebin(yaxis_name, new_ybins)
    
            ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
            for jmult in sorted(set([key[1] for key in histo.values().keys()])):
                print(*[jmult, hname], sep=", ")
                pltdir = os.path.join(outdir, args.lepton, jmult, "REL_RESO")
                if not os.path.isdir(pltdir):
                    os.makedirs(pltdir)
   
                #set_trace() 
                hslice = histo[:, jmult].integrate("jmult")
   
                    # make 1D projection along dense axes
                for dax in range(1):
                    if hslice.dense_axes()[dax].name != "mtt": continue # skip normal mtt dist
                    if dax == 0:
                        xlabel = ytitle
                        xlimits = y_lims
                        base_figname = os.path.join(pltdir, "_".join([jmult, args.lepton, "FinalBinning", "REL_RESO", yaxis_name])) if (args.binning == "Final") \
                            else os.path.join(pltdir, "_".join([jmult, args.lepton, "REL_RESO", yaxis_name]))
                    else:
                        xlabel = xtitle
                        xlimits = x_lims
                        base_figname = os.path.join(pltdir, "_".join([jmult, args.lepton, "FinalBinning", "REL_RESO", xaxis_name])) if (args.binning == "Final") \
                            else os.path.join(pltdir, "_".join([jmult, args.lepton, "REL_RESO", xaxis_name]))
   
                    #set_trace() 
                    hproj = hslice.integrate(hslice.dense_axes()[dax].name)
                    bin_edges = hproj.dense_axes()[0].edges()
                    bin_centers = hproj.dense_axes()[0].centers()
                    bin_numbers = np.arange(0, bin_edges.size) 
                        # make plots for individual categories
                    for proc in hproj.values().keys():
                        proc_histo = hproj[proc].integrate("process")
                        vals = proc_histo.values()[()]

                        mean = np.average(bin_centers, weights=vals)
                        stdev = np.sqrt(np.average((bin_centers-mean)**2, weights=vals))                    
   
                            # find FWHM of distribution
                        peak = np.array([np.argmax(vals)])
                        results_half_width, results_half_height, results_half_leftpos, results_half_rightpos = peak_widths(vals, peak, rel_height=0.5)
                        #if results_half_height.size != 1: raise ValueError(f"{results_half_height.size} peaks have been found")
                                # interpolate bin widths into bin values
                        interp_bin_minval, interp_bin_maxval = np.interp(results_half_leftpos, bin_numbers, bin_edges), np.interp(results_half_rightpos, bin_numbers, bin_edges)
                        interp_fwhm = abs(interp_bin_maxval- interp_bin_minval)[0]
 
                            # plots
                        fig, ax = plt.subplots()
                        fig.subplots_adjust(hspace=.07)
  
                        #set_trace()
                        Plotter.plot_1D(vals, bin_edges, xlimits=xlimits, xlabel=xlabel, ax=ax, label=plt_tools.get_label(proc[0], hstyles))
                        ax.hlines(results_half_height, interp_bin_minval, interp_bin_maxval, color="r")
                                
                        ax.legend(loc="upper right")
                        ax.set_ylim(0, ax.get_ylim()[1]*1.15)
    
                            # add lepton/jet multiplicity label
                        ax.text(
                            0.02, 0.82, "%s\nmean=%.4f\nstd=%.4f\nFWHM=%.4f" % (cfeatures.channel_labels[f"{args.lepton}_{jmult}"], mean, stdev, interp_fwhm),
                            fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                        )
                        hep.cms.label(ax=ax, data=False, label="Preliminary", year=cfeatures.year_labels[args.year], lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))

                        figname = "_".join([base_figname, proc[0]]) 
                        #set_trace()
                        fig.savefig(figname)
                        print(f"{figname} written")
                        plt.close(fig)

                        # make plots for all categories
                    proc_histo = hproj.integrate("process")
                    vals = proc_histo.values()[()]

                    mean = np.average(bin_centers, weights=vals)
                    stdev = np.sqrt(np.average((bin_centers-mean)**2, weights=vals))                    
    
                        # find FWHM of distribution
                    peak = np.array([np.argmax(vals)])
                    results_half_width, results_half_height, results_half_leftpos, results_half_rightpos = peak_widths(vals, peak, rel_height=0.5)
                    #if results_half_height.size != 1: raise ValueError(f"{results_half_height.size} peaks have been found")
                            # interpolate bin widths into bin values
                    interp_bin_minval, interp_bin_maxval = np.interp(results_half_leftpos, bin_numbers, bin_edges), np.interp(results_half_rightpos, bin_numbers, bin_edges)
                    interp_fwhm = abs(interp_bin_maxval- interp_bin_minval)[0]
 
                        # plots
                    fig, ax = plt.subplots()
                    fig.subplots_adjust(hspace=.07)
  
                    #set_trace()
                    Plotter.plot_1D(vals, bin_edges, xlimits=xlimits, xlabel=xlabel, ax=ax, label=plt_tools.get_label("ttJetsSL", hstyles))
                    ax.hlines(results_half_height, interp_bin_minval, interp_bin_maxval, color="r")
                            
                    ax.legend(loc="upper right")
                    ax.set_ylim(0, ax.get_ylim()[1]*1.15)
    
                        # add lepton/jet multiplicity label
                    ax.text(
                        0.02, 0.82, "%s\nmean=%.4f\nstd=%.4f\nFWHM=%.4f" % (cfeatures.channel_labels[f"{args.lepton}_{jmult}"], mean, stdev, interp_fwhm),
                        fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                    )
                    hep.cms.label(ax=ax, data=False, label="Preliminary", year=cfeatures.year_labels[args.year], lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))

                    figname = "_".join([base_figname, "ttJetsSL"]) 
                    #set_trace()
                    fig.savefig(figname)
                    print(f"{figname} written")
                    plt.close(fig)


                # make plots of relative bin width and resolution for each bin in mtt
                    # calculate bin widths
                #set_trace()
                mtt_edges = hslice.axis("mtt").edges()
                mtt_centers = hslice.axis("mtt").centers()
                #mtt_edges = mtt_edges[np.argwhere(mtt_edges == 300.)[0][0]:] # start mtt bins at 300 GeV
                mtt_min, mtt_max = mtt_edges[0], mtt_edges[-1]
                rel_mtt_bin_widths = np.array([(mtt_edges[idx+1]-mtt_edges[idx])/mtt_centers[idx] for idx in range(mtt_centers.size)])
                #rel_mtt_bin_widths = np.array([(mtt_edges[idx+1]-mtt_edges[idx])/(mtt_max-mtt_min) for idx in range(mtt_edges.size-1)])

                rel_res_bin_centers = hslice.axis("rel_reso_mtt").centers()
                rel_res_bin_edges = hslice.axis("rel_reso_mtt").edges()
                rel_res_bin_numbers = np.arange(0, rel_res_bin_edges.size)

                ### initialize plots comparing correct and all l+jets ttbar
                fig_comp, ax_comp = plt.subplots()
                fig_comp.subplots_adjust(hspace=.07)

                #set_trace()
                    # make plots for individual categories
                for proc in hproj.values().keys():
                    proc_histo = hslice[proc].integrate("process")
                    rel_res_array = np.zeros(mtt_edges.size-1)
                    rel_fwhm_array = np.zeros(mtt_edges.size-1)
                        # loop over mtt bin values and calculate mean, stdev of relative resolution
                    #set_trace()
                    for idx in range(mtt_edges.size-1):
                        mtt_low, mtt_hi = mtt_edges[idx], mtt_edges[idx+1]
                        rel_res_vals = proc_histo[mtt_low:mtt_hi, :].integrate("mtt").values()[()]

                        if np.sum(rel_res_vals) > 0:
                            bin_mean = np.average(rel_res_bin_centers, weights=rel_res_vals)
                            bin_stdev = np.sqrt(np.average((rel_res_bin_centers-bin_mean)**2, weights=rel_res_vals))
                            rel_res_array[idx] = bin_stdev

                            #set_trace()
                                # find FWHM of distribution
                            peak = np.array([np.argmax(rel_res_vals)])
                            results_half_width, results_half_height, results_half_leftpos, results_half_rightpos = peak_widths(rel_res_vals, peak, rel_height=0.5)
                            #if results_half_height.size != 1: raise ValueError(f"{results_half_height.size} peaks have been found")
                                    # interpolate bin widths into bin values
                            interp_bin_minval, interp_bin_maxval = np.interp(results_half_leftpos, rel_res_bin_numbers, rel_res_bin_edges), np.interp(results_half_rightpos, rel_res_bin_numbers, rel_res_bin_edges)
                            interp_fwhm = abs(interp_bin_maxval- interp_bin_minval)[0]
                            rel_fwhm_array[idx] = interp_fwhm
 
                        else:
                            rel_res_array[idx] = np.nan


                        # plot correct ttbar in comp figure
                    #if proc[0] == "ttJets_right": set_trace()
                    if proc[0] == "ttJets_right": Plotter.plot_1D(rel_fwhm_array, mtt_edges, xlimits=x_lims, xlabel=xtitle, ylabel="a.u.", ax=ax_comp, label=f"{plt_tools.get_label(proc[0], hstyles)} FWHM", color="b")
                    #if proc[0] == "ttJets_right": Plotter.plot_1D(rel_res_array, mtt_edges, xlimits=x_lims, xlabel=xtitle, ylabel="a.u.", ax=ax_comp, label=f"{plt_tools.get_label(proc[0], hstyles)}"+" $m_{t\\bar{t}}$ resolution", color="m")

                    #set_trace()
                    fig, ax = plt.subplots()
                    fig.subplots_adjust(hspace=.07)
  
                    Plotter.plot_1D(rel_res_array, mtt_edges, xlimits=x_lims, xlabel=xtitle, ylabel="a.u.", ax=ax, label="$m_{t\\bar{t}}$ resolution", color="r")
                    Plotter.plot_1D(rel_fwhm_array, mtt_edges, xlimits=x_lims, xlabel=xtitle, ylabel="a.u.", ax=ax, label="$m_{t\\bar{t}}$ FWHM", color="b")
                    Plotter.plot_1D(rel_mtt_bin_widths, mtt_edges, xlimits=x_lims, xlabel=xtitle, ylabel="a.u.", ax=ax, label="Relative bin width", color="k")
                            
                    ax.legend(loc="upper right", title=plt_tools.get_label(proc[0], hstyles))
                    ax.set_ylim(0, ax.get_ylim()[1]*1.15)
    
                        # add lepton/jet multiplicity label
                    ax.text(
                        0.02, 0.92, cfeatures.channel_labels[f"{args.lepton}_{jmult}"],
                        fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                    )
                    hep.cms.label(ax=ax, data=False, label="Preliminary", year=cfeatures.year_labels[args.year], lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))

                    #set_trace()
                    figname = os.path.join(pltdir, "_".join([jmult, args.lepton, "FinalBinning", "RelativeResolutionBinWidth", "vs", xaxis_name, proc[0]])) if (args.binning == "Final") \
                        else os.path.join(pltdir, "_".join([jmult, args.lepton, "RelativeResolutionBinWidth", "vs", xaxis_name, proc[0]]))
                    #set_trace()
                    fig.savefig(figname)
                    print(f"{figname} written")
                    plt.close(fig)

                    # make plots for combined categories
                proc_histo = hslice.integrate("process")
                rel_res_array = np.zeros(mtt_edges.size-1)
                rel_fwhm_array = np.zeros(mtt_edges.size-1)
                    # loop over mtt bin values and calculate mean, stdev of relative resolution
                #set_trace()
                for idx in range(mtt_edges.size-1):
                    mtt_low, mtt_hi = mtt_edges[idx], mtt_edges[idx+1]
                    rel_res_vals = proc_histo[mtt_low:mtt_hi, :].integrate("mtt").values()[()]

                    if np.sum(rel_res_vals) > 0:
                        bin_mean = np.average(rel_res_bin_centers, weights=rel_res_vals)
                        bin_stdev = np.sqrt(np.average((rel_res_bin_centers-bin_mean)**2, weights=rel_res_vals))
                        rel_res_array[idx] = bin_stdev

                            # find FWHM of distribution
                        peak = np.array([np.argmax(rel_res_vals)])
                        results_half_width, results_half_height, results_half_leftpos, results_half_rightpos = peak_widths(rel_res_vals, peak, rel_height=0.5)
                        #if results_half_height.size != 1: raise ValueError(f"{results_half_height.size} peaks have been found")
                                # interpolate bin widths into bin values
                        interp_bin_minval, interp_bin_maxval = np.interp(results_half_leftpos, rel_res_bin_numbers, rel_res_bin_edges), np.interp(results_half_rightpos, rel_res_bin_numbers, rel_res_bin_edges)
                        interp_fwhm = abs(interp_bin_maxval- interp_bin_minval)[0]
                        rel_fwhm_array[idx] = interp_fwhm
 
                    else:
                        rel_res_array[idx] = np.nan

                # plot total l+jets result in comparison figure
                #Plotter.plot_1D(rel_res_array, mtt_edges, xlimits=x_lims, xlabel=xtitle, ylabel="a.u.", ax=ax_comp, label=f"{plt_tools.get_label('ttJetsSL', hstyles)}"+" $m_{t\\bar{t}}$ resolution", color="g")
                Plotter.plot_1D(rel_fwhm_array, mtt_edges, xlimits=x_lims, xlabel=xtitle, ylabel="a.u.", ax=ax_comp, label=f"{plt_tools.get_label('ttJetsSL', hstyles)} FWHM", color="r")
                Plotter.plot_1D(rel_mtt_bin_widths, mtt_edges, xlimits=x_lims, xlabel=xtitle, ylabel="a.u.", ax=ax_comp, label="Relative bin width", color="k")
                        
                ax_comp.legend(loc="upper right")#, title=plt_tools.get_label("ttJetsSL", hstyles))
                ax_comp.set_ylim(0, ax_comp.get_ylim()[1]*1.15)
    
                    # add lepton/jet multiplicity label
                ax_comp.text(
                    0.02, 0.92, cfeatures.channel_labels[f"{args.lepton}_{jmult}"],
                    fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax_comp.transAxes
                )
                hep.cms.label(ax=ax_comp, data=False, label="Preliminary", year=cfeatures.year_labels[args.year], lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))

                #set_trace()
                #figname = os.path.join(pltdir, "_".join([jmult, args.lepton, "FinalBinning", "RelativeResolutionBinWidth", "vs", xaxis_name, "CorrectVSttJetsSL_FWHMvsSigma_Comp"])) if (args.binning == "Final") \
                #    else os.path.join(pltdir, "_".join([jmult, args.lepton, "RelativeResolutionBinWidth", "vs", xaxis_name, "CorrectVSttJetsSL_Comp"]))
                figname = os.path.join(pltdir, "_".join([jmult, args.lepton, "FinalBinning", "RelativeResolutionBinWidth", "vs", xaxis_name, "CorrectVSttJetsSL_Comp"])) if (args.binning == "Final") \
                    else os.path.join(pltdir, "_".join([jmult, args.lepton, "RelativeResolutionBinWidth", "vs", xaxis_name, "CorrectVSttJetsSL_Comp"]))
                #set_trace()
                fig_comp.savefig(figname)
                print(f"{figname} written")
                plt.close(fig_comp)


                #set_trace()
                fig, ax = plt.subplots()
                fig.subplots_adjust(hspace=.07)
  
                Plotter.plot_1D(rel_res_array, mtt_edges, xlimits=x_lims, xlabel=xtitle, ylabel="a.u.", ax=ax, label="$m_{t\\bar{t}}$ resolution", color="r")
                Plotter.plot_1D(rel_fwhm_array, mtt_edges, xlimits=x_lims, xlabel=xtitle, ylabel="a.u.", ax=ax, label="$m_{t\\bar{t}}$ FWHM", color="b")
                Plotter.plot_1D(rel_mtt_bin_widths, mtt_edges, xlimits=x_lims, xlabel=xtitle, ylabel="a.u.", ax=ax, label="Relative bin width", color="k")
                        
                ax.legend(loc="upper right", title=plt_tools.get_label("ttJetsSL", hstyles))
                ax.set_ylim(0, ax.get_ylim()[1]*1.15)
    
                    # add lepton/jet multiplicity label
                ax.text(
                    0.02, 0.92, cfeatures.channel_labels[f"{args.lepton}_{jmult}"],
                    fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                )
                hep.cms.label(ax=ax, data=False, label="Preliminary", year=cfeatures.year_labels[args.year], lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))

                #set_trace()
                figname = os.path.join(pltdir, "_".join([jmult, args.lepton, "FinalBinning", "RelativeResolutionBinWidth", "vs", xaxis_name, "ttJetsSL"])) if (args.binning == "Final") \
                    else os.path.join(pltdir, "_".join([jmult, args.lepton, "RelativeResolutionBinWidth", "vs", xaxis_name, "ttJetsSL"]))
                #set_trace()
                fig.savefig(figname)
                print(f"{figname} written")
                plt.close(fig)
