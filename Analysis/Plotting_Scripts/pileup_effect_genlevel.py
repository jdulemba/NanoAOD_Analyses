#!/usr/bin/env python
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

import Utilities.Plotter as Plotter
from coffea.hist import plot
from coffea.util import load#, save
from pdb import set_trace
import os
import Utilities.prettyjson as prettyjson
import Utilities.plot_tools as plt_tools
import numpy as np
from coffea.lookup_tools.root_converters import convert_histo_root_file
import Utilities.common_features as cfeatures

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="Specify which year to run over")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
base_jobid = os.environ["base_jobid"]
plot_outdir = os.environ["plots_dir"]
eos_dir = os.environ["eos_dir"]
analyzer = "PileupEffect_GenLevel"

mc_input_dir = os.path.join(eos_dir, "results", base_jobid, analyzer, args.year)
f_ext = "TOT.coffea"
mc_fnames = sorted([os.path.join(mc_input_dir, fname) for fname in os.listdir(mc_input_dir) if fname.endswith(f_ext)])
hdict = plt_tools.add_coffea_files(mc_fnames) if len(mc_fnames) > 1 else load(mc_fnames[0])

outdir = os.path.join(plot_outdir, base_jobid, analyzer, args.year)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

    # get lumi info
data_lumi_dict = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())
lumi_to_use = (data_lumi_dict[args.year]["Muons"]+data_lumi_dict[args.year]["Electrons"])/2000.

year_to_use = cfeatures.year_labels[args.year]

variables = {
    "nTrueInt" : ("Number of Pileup Interactions", 1, (0, 80), True),
    #"nTrueInt" : ("Number of Pileup Interactions", 1, (0, 100), True),
    "nPV" : ("Number of Primary Vertices", 1, (0, 50), False),
    "rho" : ("Fastjet $\\rho$", 10, (0, 50), False),
}

rewt_style_dict = {
    "nosys" : {"label" : "Simulation, Uncorrected", "color" : "#e41a1c", "linestyle" : "-"}, ## red
    "Pileup" : {"label" : "Simulation, Corrected", "color" : "k", "linestyle" : "-"}, ## black
    #"PileupUp" : {"label" : "Corrected, Up", "color" : "#e41a1c", "linestyle" : "-"}, ## red
    #"PileupDown" : {"label" : "Corrected, Down", "color" : "#377eb8", "linestyle" : "-"}, ## blue
}
data_style_dict = {
        "data" : {"label" : "data, 69.2mb", "facecolor": "#4daf4a"}, ## green
}

central_files = {
    "2016APV" : {
        "Down" : "PileupHistogram-goldenJSON-13tev-2016-preVFP-66000ub-99bins.root",
        "Cen"  : "PileupHistogram-goldenJSON-13tev-2016-preVFP-69200ub-99bins.root",
        "Up"   : "PileupHistogram-goldenJSON-13tev-2016-preVFP-72400ub-99bins.root",
    },
    "2016" : {
        "Down" : "PileupHistogram-goldenJSON-13tev-2016-postVFP-66000ub-99bins.root",
        "Cen"  : "PileupHistogram-goldenJSON-13tev-2016-postVFP-69200ub-99bins.root",
        "Up"   : "PileupHistogram-goldenJSON-13tev-2016-postVFP-72400ub-99bins.root",
    },
    "2017" : {
        "Down" : "PileupHistogram-goldenJSON-13tev-2017-66000ub-99bins.root",
        "Cen"  : "PileupHistogram-goldenJSON-13tev-2017-69200ub-99bins.root",
        "Up"   : "PileupHistogram-goldenJSON-13tev-2017-72400ub-99bins.root",
    },
    "2018" : {
        "Down" : "PileupHistogram-goldenJSON-13tev-2018-66000ub-99bins.root",
        "Cen"  : "PileupHistogram-goldenJSON-13tev-2018-69200ub-99bins.root",
        "Up"   : "PileupHistogram-goldenJSON-13tev-2018-72400ub-99bins.root",
    },
}

data_input_dir = os.path.join(proj_dir, "inputs", "data", base_jobid, "Pileup")


#set_trace()
for hname in variables.keys():
    xtitle, rebinning, x_lims, withData = variables[hname]
    histo = hdict[hname]["ttJetsSL"].integrate("dataset").copy()

    axes_to_sum = (histo.dense_axes()[0].name,)
    if isinstance(rebinning, np.ndarray):
        new_xbins = hist.Bin(axes_to_sum[0], axes_to_sum[0], rebinning)
    elif isinstance(rebinning, float) or isinstance(rebinning, int):
        new_xbins = rebinning
    histo = histo.rebin(*axes_to_sum, new_xbins)

        ## plot simulation and nominal data
        ## plot data
    if withData:
        fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
        fig.subplots_adjust(hspace=.07)
                    # central
        data_pu_central = convert_histo_root_file(os.path.join(data_input_dir, central_files[args.year]["Cen"]))
        data_pu_dict = Plotter.root_converters_dict_to_hist(data_pu_central, vars=["pileup"],
            dense_axes_list=[{"name" : "pu_nTrueInt", "idx" : 0}],
            sparse_axes_list=[{"name": "dataset", "label" : "Event Process", "fill" : "data"}],
        )
        data_pu_histos = data_pu_dict["pileup"].copy()
            # plot data
        hep.plot.histplot(data_pu_histos.values()[("data",)], data_pu_histos.dense_axes()[0].edges(), ax=ax, histtype="fill", density=True, **data_style_dict["data"])

        for rewt, style_dict in rewt_style_dict.items():
                # plot yields
            hep.plot.histplot(histo.values()[(rewt,)], histo.dense_axes()[0].edges(), ax=ax, histtype="step", density=True, **style_dict)

            #set_trace()
                # plot data/MC ratio 
            ratio_masked_vals, ratio_masked_bins = Plotter.get_ratio_arrays(num_vals=np.r_[data_pu_histos.values()[("data",)]/np.sum(data_pu_histos.values()[("data",)]), (data_pu_histos.values()[("data",)]/np.sum(data_pu_histos.values()[("data",)]))[-1]],
                denom_vals=histo.values()[(rewt,)]/np.sum(histo.values()[(rewt,)]), input_bins=histo.dense_axes()[0].edges())
            rax.step(ratio_masked_bins, ratio_masked_vals, where='post', **style_dict)

        ax.autoscale()
        ax.set_xlim(x_lims)
        ax.set_ylabel("Probability Density")
        ax.set_xlabel(None)

        rax.autoscale()
        rax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
        rax.set_ylabel("data/MC")
        rax.set_xlabel(xtitle)
        rax.set_xlim(x_lims)
        rax.set_ylim(0.5, 1.5)

        ## set legend and corresponding colors
        ax.legend(loc="upper right")

            # add cms label
        hep.cms.label(ax=ax, data=withData, label="Preliminary", year=year_to_use, lumi=round(lumi_to_use, 1))
        
        figname = os.path.join(outdir, f"{args.year}_{analyzer}_{hname}_dataMC_Comp_Norm")
        fig.savefig(figname)
        print(f"{figname} written")
        plt.close(fig)

    else:
        fig, ax = plt.subplots()
        fig.subplots_adjust(hspace=.07)

        for rewt, style_dict in rewt_style_dict.items():
                # plot yields
            hep.plot.histplot(histo.values()[(rewt,)], histo.dense_axes()[0].edges(), ax=ax, histtype="step", density=True, **style_dict)


        ax.autoscale()
        ax.set_xlim(x_lims)
        ax.set_ylabel("Probability Density")
        ax.set_xlabel(xtitle)

        ## set legend and corresponding colors
        ax.legend(loc="upper right")

            # add cms label
        hep.cms.label(ax=ax, data=withData, label="Preliminary", year=year_to_use, lumi=round(lumi_to_use, 1))
        
        figname = os.path.join(outdir, f"{args.year}_{analyzer}_{hname}_Norm")
        fig.savefig(figname)
        print(f"{figname} written")
        plt.close(fig)


toc = time.time()
print("Total time: %.1f" % (toc - tic))
