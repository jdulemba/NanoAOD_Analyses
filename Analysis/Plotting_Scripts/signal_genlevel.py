import time
tic = time.time()

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
import fnmatch, re
import Utilities.Plotter as Plotter
from Utilities.styles import styles as hstyles
import Utilities.final_analysis_binning as final_binning
import Utilities.common_features as cfeatures

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2018"], help="What year is the ntuple from.")
#parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
base_jobid = os.environ["base_jobid"]
jobid = os.environ["jobid"]
analyzer = "signal_genlevel"
plot_outdir = os.environ["plots_dir"]
eos_dir = os.environ["eos_dir"]

input_dir = os.path.join(eos_dir, "results", f"{args.year}_{base_jobid}", analyzer)
fnames = [os.path.join(input_dir, fname) for fname in fnmatch.filter(os.listdir(input_dir), "*TOT.coffea")]
hdict = load(fnames[0])

outdir = os.path.join(plot_outdir, f"{args.year}_{base_jobid}", analyzer)
if not os.path.isdir(outdir):
    os.makedirs(outdir)


linearize_binning = (
    final_binning.mtt_binning,
    final_binning.ctstar_abs_binning
)

#set_trace()
tlep_ctstar_binlabels = [r"%s $\in [%s, %s)$" % (cfeatures.variable_names_to_labels["tlep_ctstar_abs"], linearize_binning[1][bin], linearize_binning[1][bin+1]) for bin in range(len(linearize_binning[1])-1)]
tlep_ctstar_binlabels[-1] = r"%s $\in [%s, %s]$" % (cfeatures.variable_names_to_labels["tlep_ctstar_abs"], linearize_binning[1][-2], linearize_binning[1][-1])
tlep_ctstar_bin_locs = np.linspace((len(linearize_binning[0])-1)/2, (len(linearize_binning[0])-1)*(len(linearize_binning[1])-1) - (len(linearize_binning[0])-1)/2, len(linearize_binning[1])-1)
top_ctstar_binlabels = [r"%s $\in [%s, %s)$" % (cfeatures.variable_names_to_labels["top_ctstar_abs"], linearize_binning[1][bin], linearize_binning[1][bin+1]) for bin in range(len(linearize_binning[1])-1)]
top_ctstar_binlabels[-1] = r"%s $\in [%s, %s]$" % (cfeatures.variable_names_to_labels["top_ctstar_abs"], linearize_binning[1][-2], linearize_binning[1][-1])
top_ctstar_bin_locs = np.linspace((len(linearize_binning[0])-1)/2, (len(linearize_binning[0])-1)*(len(linearize_binning[1])-1) - (len(linearize_binning[0])-1)/2, len(linearize_binning[1])-1)
mtt_vals_to_plot = np.array([400, 600, 1000])
mtt_tiled_labels = np.tile(linearize_binning[0][:-1], linearize_binning[1].size-1)
mtt_bin_inds_to_plot = np.where(np.in1d(mtt_tiled_labels, mtt_vals_to_plot))[0]
mtt_bins_to_plot = np.tile(mtt_vals_to_plot, linearize_binning[1].size-1)

variables = {
    "mtt" : (cfeatures.variable_names_to_labels["mtt"], 1, (340., 1100.)),
    #"mtop" : (cfeatures.variable_names_to_labels["mtop"], 2, (0., 300.)),
    #"pt_top" : (cfeatures.variable_names_to_labels["pt_top"], 2, (0., 500.)),
    #"pt_tt" : (cfeatures.variable_names_to_labels["pt_tt"], 2, (0., 500.)),
    #"eta_top" : (cfeatures.variable_names_to_labels["eta_top"], 2, (-4., 4.)),
    #"eta_tt" : (cfeatures.variable_names_to_labels["eta_tt"], 2, (-4., 4.)),
    "top_ctstar" : (cfeatures.variable_names_to_labels["top_ctstar"], 1, (-1., 1.)),
    "top_ctstar_abs" : (cfeatures.variable_names_to_labels["top_ctstar_abs"], 1, (0., 1.)),
    #"mtt_vs_top_ctstar_abs" : (cfeatures.variable_names_to_labels["mtt"], cfeatures.variable_names_to_labels["top_ctstar_abs"], linearize_binning[0], linearize_binning[1],
    #    (linearize_binning[0][0], linearize_binning[0][-1]), (linearize_binning[1][0], linearize_binning[1][-1]), True, None, (top_ctstar_binlabels, top_ctstar_bin_locs)),
    #"tlep_ctstar" : (cfeatures.variable_names_to_labels["tlep_ctstar"], 1, (-1., 1.)),
    #"tlep_ctstar_abs" : (cfeatures.variable_names_to_labels["tlep_ctstar_abs"], 1, (0., 1.)),
    #"mtt_vs_tlep_ctstar_abs" : (cfeatures.variable_names_to_labels["mtt"], cfeatures.variable_names_to_labels["tlep_ctstar_abs"], linearize_binning[0], linearize_binning[1],
    #    (linearize_binning[0][0], linearize_binning[0][-1]), (linearize_binning[1][0], linearize_binning[1][-1]), True, None, (tlep_ctstar_binlabels, tlep_ctstar_bin_locs)),
}


    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())[args.year]
lumi_correction = load(os.path.join(proj_dir, "Corrections", base_jobid, "MC_LumiWeights.coffea"))[args.year]["Muons"]
##set_trace()
## make groups based on process
names = [dataset for dataset in sorted(set([key[0] for key in hdict["mtt"].values().keys()]))] # get dataset names in hists
process = hist.Cat("process", "Process", sorting="placement")
process_cat = "dataset"
process_groups = plt_tools.make_dataset_groups("Muon", args.year, samples=names, bkgdict="templates", sigdict="MC_indiv")
#set_trace()

    # scale and group hists by process
for hname in hdict.keys():
    if "cutflow" in hname: continue
    #set_trace()
    hdict[hname].scale(lumi_correction, axis="dataset") # scale hists
    hdict[hname] = hdict[hname].group(process_cat, process, process_groups) # group by process
    hdict[hname] = hdict[hname][:, "nosys"].integrate("sys") # only pick out specified lepton


masses = [365, 400, 500, 600, 800, 1000]
widths = ["W2p5", "W10p0", "W25p0"]
parities = ["A", "H"]
couplings = ["Res", "Int", "Total"]

## make plots
for hname in variables.keys():
    if hname not in hdict.keys():
        raise ValueError(f"{hname} not found in file")
    #set_trace()
    histo = hdict[hname].copy()

    #is2d = histo.dense_dim() == 2
    #axes_to_sum = (histo.dense_axes()[0].name,) if not is2d else (histo.dense_axes()[0].name, histo.dense_axes()[1].name)

    if histo.dense_dim() == 1:
        xtitle, rebinning, x_lims = variables[hname]
        xaxis_name = histo.dense_axes()[0].name

        if rebinning != 1:
            axes_to_sum = (histo.dense_axes()[0].name,)
            histo = histo.rebin(*axes_to_sum, rebinning)

    if histo.dense_dim() == 2:
        print("2D hists not supported yet")
        set_trace()
        xtitle, ytitle, xrebinning, yrebinning, x_lims, y_lims = variables[hname]

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

    bin_edges = histo.axis(xaxis_name).edges()

    # make plots comparing different masses for each width
    pltdir = os.path.join(outdir, "Comp_Masses")
    if not os.path.isdir(pltdir):
        os.makedirs(pltdir)

    for width in widths:
        for par in parities:
            for coupling in couplings:
                fig, ax = plt.subplots()
                fig.subplots_adjust(hspace=.07)
                fig_norm, ax_norm = plt.subplots()
                fig_norm.subplots_adjust(hspace=.07)

                for mass in masses:
                    if coupling == "Int":
                        #set_trace()
                        pos_vals = histo[f"{par}toTTJetsSL_M{mass}_{width}_{coupling}_pos"].integrate("process").values()[()]
                        neg_vals = histo[f"{par}toTTJetsSL_M{mass}_{width}_{coupling}_neg"].integrate("process").values()[()]
                        vals = pos_vals + neg_vals
                        norm_vals = vals/np.sum(np.abs(pos_vals) + np.abs(neg_vals))
                    elif coupling == "Res":
                        vals = histo[f"{par}toTTJetsSL_M{mass}_{width}_{coupling}"].integrate("process").values()[()]
                        norm_vals = vals/np.sum(vals)
                    elif coupling == "Total":
                        #set_trace()
                        pos_vals = histo[f"{par}toTTJetsSL_M{mass}_{width}_Int_pos"].integrate("process").values()[()]
                        neg_vals = histo[f"{par}toTTJetsSL_M{mass}_{width}_Int_neg"].integrate("process").values()[()]
                        res_vals = histo[f"{par}toTTJetsSL_M{mass}_{width}_Res"].integrate("process").values()[()]
                        vals = pos_vals + neg_vals + res_vals
                        norm_vals = vals/np.sum(np.abs(pos_vals) + np.abs(neg_vals) + np.abs(res_vals))
                    else:
                        raise ValueError(f"{coupling} not supported")

                    if coupling == "Res":
                        scale = 1
                        label = str(mass)
                        #hep.plot.histplot(vals, bin_edges, ax=ax, label=str(mass), linestyle="-", color=hstyles[f"M{mass}"]["color"], histtype="step")
                    elif coupling == "Int":
                        if width == "W2p5":
                            #if mass >= 800:
                            if mass == 800:
                                scale = 10
                                label = f"{str(mass)}x10"
                            elif mass == 1000:
                                scale = 50
                                label = f"{str(mass)}x50"
                            else:
                                scale = 1
                                label = str(mass)
                        else:
                            #if mass == 800:
                            if mass >= 800:
                                scale = 10
                                label = f"{str(mass)}x10"
                            #elif mass == 1000:
                            #    scale = 50
                            #    label = f"{str(mass)}x50"
                            else:
                                scale = 1
                                label = str(mass)
                    else:
                        if width == "W2p5":
                            if mass == 800:
                                scale = 10
                                label = f"{str(mass)}x10"
                            elif mass == 1000:
                                scale = 50
                                label = f"{str(mass)}x50"
                            else:
                                scale = 1
                                label = str(mass)
                        else:
                            #if mass == 800:
                            if mass >= 800:
                                scale = 10
                                label = f"{str(mass)}x10"
                            #elif mass == 1000:
                            #    scale = 50
                            #    label = f"{str(mass)}x50"
                            else:
                                scale = 1
                                label = str(mass)
                    hep.plot.histplot(vals*scale, bin_edges, ax=ax, label=label, linestyle="-", color=hstyles[f"M{mass}"]["color"], histtype="step")
                    #hep.plot.histplot(vals, bin_edges, ax=ax, label=str(mass), linestyle="-", color=hstyles[f"M{mass}"]["color"], histtype="step")
                    hep.plot.histplot(norm_vals, bin_edges, ax=ax_norm, label=str(mass), linestyle="-", color=hstyles[f"M{mass}"]["color"], histtype="step")

                if coupling == "Res": ax.set_yscale("log")                            
                ax.legend(loc="upper right", title="$\mathrm{m_%s}$ [GeV]" % par, ncol=2)
                ax.set_ylim(1., ax.get_ylim()[1]*5) if coupling == "Res" else ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.3)
                ax.set_xlim(x_lims)
                ax.set_xlabel(xtitle)
                ax.set_ylabel("Events")
                ax.axhline(0, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})

                    # add lepton/jet multiplicity label
                ax.text(
                    0.02, 0.90, "$\mathrm{\Gamma}$/$\mathrm{m_%s}$=%s%%, %s\nparton level" % (par, width.strip("W").replace("p", "."), coupling),
                    fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                )
                hep.cms.label(ax=ax, data=False, label="Preliminary", year=cfeatures.year_labels[args.year], lumi=round(data_lumi_year["Muons"]/1000., 1))

                figname = os.path.join(pltdir, f"{par}toTTJetsSL_{width}_{coupling}_MassComp_{hname}")
                fig.savefig(figname)
                print(f"{figname} written")
                plt.close(fig)
                #set_trace()

                # plot normalized figures            
                ax_norm.legend(loc="upper right", title="$\mathrm{m_%s}$ [GeV]" % par, ncol=2)
                ax_norm.set_ylim(ax_norm.get_ylim()[0], ax_norm.get_ylim()[1]*1.10 if coupling == "Res" else ax_norm.get_ylim()[1]*1.3)
                #ax_norm.set_ylim(ax_norm.get_ylim()[0], ax_norm.get_ylim()[1]*1.3 if coupling == "Int" else ax_norm.get_ylim()[1]*1.10)
                ax_norm.set_xlim(x_lims)
                ax_norm.set_xlabel(xtitle)
                ax_norm.set_ylabel("a.u.")
                ax_norm.axhline(0, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})

                    # add lepton/jet multiplicity label
                ax_norm.text(
                    0.02, 0.90, "$\mathrm{\Gamma}$/$\mathrm{m_%s}$=%s%%, %s\nparton level" % (par, width.strip("W").replace("p", "."), coupling),
                    fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax_norm.transAxes
                )
                hep.cms.label(ax=ax_norm, data=False, label="Preliminary", year=cfeatures.year_labels[args.year])

                figname_norm = os.path.join(pltdir, f"{par}toTTJetsSL_{width}_{coupling}_MassComp_{hname}_Norm")
                fig_norm.savefig(figname_norm)
                print(f"{figname_norm} written")
                plt.close(fig_norm)
                #if coupling == "Total": set_trace()
            
toc = time.time()
print("Total time: %.1f" % (toc - tic))
