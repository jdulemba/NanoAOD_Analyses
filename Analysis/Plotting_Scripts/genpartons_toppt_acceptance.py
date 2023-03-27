#! /bin/env python

import time
tic = time.time()

# matplotlib
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import mplhep as hep
plt.style.use(hep.style.CMS)
plt.switch_backend("agg")
from matplotlib import rcParams
rcParams["font.size"] = 20
rcParams["savefig.format"] = "pdf"
rcParams["savefig.bbox"] = "tight"

from coffea.hist import plot
from coffea.util import load
from pdb import set_trace
import os
import Utilities.plot_tools as plt_tools
import Utilities.prettyjson as prettyjson
import re
from coffea import hist
import numpy as np
import fnmatch
import Utilities.Plotter as Plotter
import Utilities.systematics as systematics
from copy import deepcopy
from Utilities.styles import styles as styles
import Utilities.final_analysis_binning as final_binning
import Utilities.common_features as cfeatures

base_jobid = os.environ["base_jobid"]

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2018"], help="What year is the ntuple from.")
#parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
analyzer = "genpartons_toppt_acceptance"
eos_dir = os.environ["eos_dir"]
plots_dir = os.environ["plots_dir"]

#set_trace()
base_ext = "*TOT.coffea"
input_dir = os.path.join(eos_dir, "results", jobid, analyzer, args.year)
fnames = fnmatch.filter(os.listdir(input_dir), base_ext)
fnames = [os.path.join(input_dir, fname) for fname in fnames]
hdict = load(fnames[0])
#set_trace()

outdir = os.path.join(plots_dir, jobid, analyzer, args.year)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

year_to_use = cfeatures.year_labels[args.year]


    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())[args.year]
lumi_correction = load(os.path.join(proj_dir, "Corrections", base_jobid, "MC_LumiWeights.coffea"))[args.year]["Muons"]
        # scale ttJets events, split by reconstruction type, by normal ttJets lumi correction
ttJets_permcats = ["*right", "*matchable", "*unmatchable", "*sl_tau", "*other"]
names = [dataset for dataset in sorted(set([key[0] for key in hdict["Acceptance"].values().keys()]))] # get dataset names in hists
ttJets_cats = [name for name in names if any([fnmatch.fnmatch(name, cat) for cat in ttJets_permcats])] # gets ttJets(_PS)_other, ...
if len(ttJets_cats) > 0:
    for tt_cat in ttJets_cats:
        ttJets_lumi_topo = "_".join(tt_cat.split("_")[:-2]) if "sl_tau" in tt_cat else "_".join(tt_cat.split("_")[:-1]) # gets ttJets[SL, Had, DiLep] or ttJets_PS
        ttJets_eff_lumi = lumi_correction[ttJets_lumi_topo]
        lumi_correction.update({tt_cat: ttJets_eff_lumi})

## make groups based on process
process_groups = plt_tools.make_dataset_groups("Muon", args.year, samples=names, bkgdict="templates", sigdict="MC_indiv")
tt_LHEscale_wts_name_dict = {
    "FACTORDown" : "uF_down",
    "FACTORUp"   : "uF_up",
    "RENORMDown" : "uR_down",
    "RENORMUp"   : "uR_up",
}
if sorted(set(process_groups.keys())) != ["TT"]:
    raise ValueError("Expected only ttbar events in this file")

    # scale and group hists by process
for hname in hdict.keys():
    if hname not in hdict.keys():
        print(f"{hname} not found in file")
        continue
    if "cutflow" in hname: continue
    #set_trace()
    hdict[hname].scale(lumi_correction, axis="dataset") # scale hists
    for nom_tt in process_groups["TT"]:
        for cat in ttJets_permcats:
            # rescale LHEscale systematics correctly
            for sysname, dname in tt_LHEscale_wts_name_dict.items():
                if f"{nom_tt}_{cat[1:]}" not in ttJets_cats: continue
                #print(f"{nom_tt}_{cat[1:]}_{dname}")
                lhe_scale = lumi_correction[f"{nom_tt}_{dname}"]/lumi_correction[f"{nom_tt}_{cat[1:]}"]
                hdict[hname].scale({(f"{nom_tt}_{cat[1:]}", sysname) : lhe_scale}, axis=("dataset", "sys"))
    hdict[hname] = hdict[hname][:].integrate("dataset")

styles_dict = {
    "nosys" : {"color":"k", "linestyle":"-", "label":"Nominal"},
    "TopPtUp" : {"color":"r", "linestyle":"-", "label":"$p_{T}$($t$) Up"},
    "TopPtDown" : {"color":"r", "linestyle":"--", "label":"$p_{T}$($t$) Down"},
    "RENORMUp" : {"color":"b", "linestyle":"-", "label":"$\\mu_{R}$ Up"},
    "RENORMDown" : {"color":"b", "linestyle":"--", "label":"$\\mu_{R}$ Down"},
    "FACTORUp" : {"color":"g", "linestyle":"-", "label":"$\\mu_{F}$ Up"},
    "FACTORDown" : {"color":"g", "linestyle":"--", "label":"$\\mu_{F}$ Down"},
}

    ## make plots
#set_trace()
histo = hdict["Acceptance"]

## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
for jmult in sorted(set([key[1] for key in histo.values().keys()])):
    print(jmult)
    pltdir = os.path.join(outdir, jmult)
    if not os.path.isdir(pltdir):
        os.makedirs(pltdir)

    fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
    fig.subplots_adjust(hspace=.07)

    hslice = histo[:, jmult].integrate("jmult")
    x_lims = (min(hslice.dense_axes()[0].edges()), max(hslice.dense_axes()[0].edges()))

    rows = [("Lumi: %s fb^-1" % format(round(data_lumi_year["Muons"]/1000., 1), ".1f"), "", "", "")]
    rows += [("Sys", "All", "Kin. Cuts", "Frac")]

        # plot nominal
    hep.plot.histplot(hslice.values()[("nosys",)], hslice.dense_axes()[0].edges(), ax=ax, histtype="step", **styles_dict["nosys"])
    rows += [("Nominal", format(hslice.values()[("nosys",)][0], ".1f"), format(hslice.values()[("nosys",)][1], ".1f"), format(hslice.values()[("nosys",)][1]/hslice.values()[("nosys",)][0], "0.4f"))]

        # plot systematics
    for sys in sorted(set([key[0] for key in hslice.values().keys()])):
        if sys == "nosys": continue
        hep.plot.histplot(hslice.values()[(sys,)], hslice.dense_axes()[0].edges(), ax=ax, histtype="step", **styles_dict[sys])

        masked_vals, masked_bins = Plotter.get_ratio_arrays(num_vals=hslice.values()[(sys,)], denom_vals=hslice.values()[("nosys",)], input_bins=hslice.dense_axes()[0].edges())
        rax.step(masked_bins, masked_vals, where="post", **styles_dict[sys])
        #set_trace()
        rows += [(sys, format(hslice.values()[(sys,)][0], ".1f"), format(hslice.values()[(sys,)][1], ".1f"), format(hslice.values()[(sys,)][1]/hslice.values()[(sys,)][0], "0.4f"))]

    #set_trace()
        # save and print acceptance values
    frac_name = os.path.join(pltdir, f"{jmult}_yields_and_fracs.txt")
    plt_tools.print_table(rows, filename=frac_name, print_output=True, header_line=1)
    print(f"{frac_name} written")

    ax.set_yscale("log")
    ax.legend(loc="upper right", ncol=2)
    ax.autoscale()
    ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.5)
    #ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.15)
    ax.set_xlim(x_lims)
    ax.set_xlabel(None)
    ax.set_ylabel("Events")
    ax.set_xticks(np.array([0.5, 1.5]))

    rax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
    rax.autoscale()
    rax.set_xlim(x_lims)
    rax.set_xlabel(None)
    rax.set_ylabel("Ratio to Nominal")
    rax.set_xticks(np.array([0.5, 1.5]))
    rax.set_xticklabels(["All", "Kin. Cuts"])#, rotation='vertical', fontsize=18)

        # add lepton/jet multiplicity label
    ##set_trace()
    ax.text(
        0.02, 0.86, f"{styles['ttJetsSL']['name']}\nparton level",
        fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
    )
    hep.cms.label(ax=ax, data=False, label="Preliminary", year=year_to_use, lumi=round(data_lumi_year["Muons"]/1000., 1))

    #set_trace()
    figname = os.path.join(pltdir, "_".join([jmult, "Acceptance"]))
    fig.savefig(figname)
    print(f"{figname} written")
    plt.close(fig)

toc = time.time()
print("Total time: %.1f" % (toc - tic))
