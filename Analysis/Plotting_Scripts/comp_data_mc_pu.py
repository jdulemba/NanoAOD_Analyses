#!/usr/bin/env python

# matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend("agg")
from matplotlib import rcParams
rcParams["font.size"] = 20
rcParams["savefig.format"] = "pdf"
#rcParams["savefig.format"] = "png"
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

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
base_jobid = os.environ["base_jobid"]
analyzer = "meta_info"
plot_outdir = os.environ["plots_dir"]
eos_dir = os.environ["eos_dir"]

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="Specify which year to run over")
args = parser.parse_args()

mc_input_dir = os.path.join(eos_dir, "results", f"{args.year}_{base_jobid}", analyzer)
f_ext = "TOT.coffea"
mc_fnames = sorted([os.path.join(mc_input_dir, fname) for fname in os.listdir(mc_input_dir) if fname.endswith(f_ext)])
mc_hdict = plt_tools.add_coffea_files(mc_fnames) if len(mc_fnames) > 1 else load(mc_fnames[0])

    # get hists
mc_nTrueInt_histo = mc_hdict["PU_nTrueInt"]

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

data_input_dir = os.path.join(proj_dir, "inputs", "data", base_jobid, "Pileup", args.year)
            # central
data_pu_central = convert_histo_root_file(os.path.join(data_input_dir, central_files[args.year]["Cen"]))
data_pu_dict = Plotter.root_converters_dict_to_hist(data_pu_central, vars=["pileup"],
    dense_axes_list=[{"name" : "pu_nTrueInt", "idx" : 0}],
    sparse_axes_list=[{"name": "dataset", "label" : "Event Process", "fill" : "data"}],
)
#set_trace()
            # up
data_pu_up = convert_histo_root_file(os.path.join(data_input_dir, central_files[args.year]["Up"]))
data_pu_up_dict = Plotter.root_converters_dict_to_hist(data_pu_up, vars=["pileup"],
    dense_axes_list=[{"name" : "pu_nTrueInt", "idx" : 0}],
    sparse_axes_list=[{"name": "dataset", "label" : "Event Process", "fill" : "data_up"}]
)
            # down
data_pu_dw = convert_histo_root_file(os.path.join(data_input_dir, central_files[args.year]["Down"]))
data_pu_dw_dict = Plotter.root_converters_dict_to_hist(data_pu_dw, vars=["pileup"],
    dense_axes_list=[{"name" : "pu_nTrueInt", "idx" : 0}],
    sparse_axes_list=[{"name": "dataset", "label" : "Event Process", "fill" : "data_down"}]
)
    # combine nominal, up, and down into one histo
data_pu_histos = data_pu_dict["pileup"].copy()
data_pu_histos.add(data_pu_up_dict["pileup"].copy())
data_pu_histos.add(data_pu_dw_dict["pileup"].copy())

outdir = os.path.join(plot_outdir, f"{args.year}_{base_jobid}", "comp_data_mc_pu")
if not os.path.isdir(outdir):
    os.makedirs(outdir)

    # get lumi info
data_lumi_dict = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())
lumi_to_use = (data_lumi_dict[args.year]["Muons"]+data_lumi_dict[args.year]["Electrons"])/2000.

year_to_use = cfeatures.year_labels[args.year]

## plot histograms
    ## plot simulation and nominal data
fig, ax = plt.subplots()
fig.subplots_adjust(hspace=.07)

# plot MC
mc_pu = mc_nTrueInt_histo["ttJetsSL"]
plot.plot1d(mc_pu, ax=ax, clear=False, density=True, fill_opts={"facecolor": "#4daf4a"}) # green

# plot data
plot.plot1d(data_pu_histos[("data",)], ax=ax, clear=False, line_opts={"linestyle" : "-", "color" : "k"}, density=True)

ax.set_ylabel("Probability Density")
ax.set_xlabel("Number of Pileup Interactions")
ax.set_xlim((0, 100))
ax.autoscale()#axis="x", tight=True)
ax.set_ylim(0, ax.get_ylim()[1]*1.15)

    ## set legend and corresponding colors
handles, labels = ax.get_legend_handles_labels()
for idx, sample in enumerate(labels):
    if sample == "data":
        continue
    else:
        labels[idx] = "Simulation"
# call ax.legend() with the new values
order = [1, 0]
ax.legend([handles[i] for i in order], [labels[i] for i in order], loc="upper right")

    # add cms label
hep.cms.label(ax=ax, data=True, label="Preliminary", year=year_to_use, lumi=round(lumi_to_use, 1))

figname = os.path.join(outdir, f"{args.year}_data_mc_pileup_comp_Norm_New")
#figname = os.path.join(outdir, f"{args.year}_data_mc_pileup_comp_Norm")
fig.savefig(figname)
print(f"{figname} written")
plt.close(fig)


    ## plot simulation and all data
fig_all, ax_all = plt.subplots()
fig_all.subplots_adjust(hspace=.07)

# plot MC
plot.plot1d(mc_pu, ax=ax_all, clear=False, density=True, fill_opts={"facecolor": "#4daf4a"}) # green

# plot nominal data
plot.plot1d(data_pu_histos[("data",)], ax=ax_all, clear=False, line_opts={"linestyle" : "-", "color" : "k"}, density=True,)

# plot up data
plot.plot1d(data_pu_histos[("data_up",)], ax=ax_all, clear=False, line_opts={"linestyle" : "-", "color" : "r"}, density=True,)

# plot down data
plot.plot1d(data_pu_histos[("data_down",)], ax=ax_all, clear=False, line_opts={"linestyle" : "-", "color" : "b"}, density=True,)

ax_all.set_ylabel("Probability Density")
ax_all.set_xlabel("Number of Pileup Interactions")
ax_all.set_xlim((0, 100))
ax_all.autoscale()#axis="x", tight=True)
ax_all.set_ylim(0, ax_all.get_ylim()[1]*1.15)

    ## set legend and corresponding colors
handles, labels = ax_all.get_legend_handles_labels()
for idx, sample in enumerate(labels):
    if sample == "data":
        continue
    elif sample == "data_up":
        labels[idx] = "data, Up"
    elif sample == "data_down":
        labels[idx] = "data, Down"
    else:
        labels[idx] = "Simulation"
# call ax.legend() with the new values
order = [1, 2, 3, 0]
ax_all.legend([handles[i] for i in order], [labels[i] for i in order], loc="upper right")

    # add cms label
hep.cms.label(ax=ax_all, data=True, label="Preliminary", year=year_to_use, lumi=round(lumi_to_use, 1))

figname_all = os.path.join(outdir, f"{args.year}_data_variations_mc_pileup_comp_Norm_New")
#figname_all = os.path.join(outdir, f"{args.year}_data_variations_mc_pileup_comp_Norm")
fig_all.savefig(figname_all)
print(f"{figname_all} written")
plt.close(fig_all)
