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
import Utilities.prettyjson as prettyjson
import Utilities.common_features as cfeatures
import numpy as np
from coffea.lookup_tools.root_converters import convert_histo_root_file
from coffea.lookup_tools.dense_lookup import dense_lookup
from pdb import set_trace
import os

proj_dir = os.environ["PROJECT_DIR"]
base_jobid = os.environ["base_jobid"]
eos_dir = os.environ["eos_dir"]
plot_outdir = os.environ["plots_dir"]

outdir = os.path.join(plot_outdir, base_jobid, "Compare_PileupProfiles") 
if not os.path.isdir(outdir):
    os.makedirs(outdir)

pu_path = os.path.join(proj_dir, "inputs", "data", base_jobid, "Pileup")

    # get lumi info
data_lumi_dict = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())

years_to_run = ["2016APV", "2016", "2017", "2018"]

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

for year in years_to_run:

    lumi_to_use = (data_lumi_dict[year]["Muons"]+data_lumi_dict[year]["Electrons"])/2000.

    # files that I made
        # get nominal distribution
    ur_central = convert_histo_root_file(os.path.join(pu_path, f"{year}_data.meta.pu.root"))
    ur_cen_vals, ur_bins = ur_central[("pileup", "dense_lookup")]

        # get up variation
    ur_up = convert_histo_root_file(os.path.join(pu_path, f"{year}_data.meta.pu_up.root"))
    ur_up_vals, _ = ur_up[("pileup", "dense_lookup")]

        # get down variation
    ur_dw = convert_histo_root_file(os.path.join(pu_path, f"{year}_data.meta.pu_down.root"))
    ur_dw_vals, _ = ur_dw[("pileup", "dense_lookup")]
    
    # central files
        # get nominal distribution
    central_central = convert_histo_root_file(os.path.join(pu_path, year, central_files[year]["Cen"]))
    central_cen_vals, central_bins = central_central[("pileup", "dense_lookup")]

        # get up variation
    central_up = convert_histo_root_file(os.path.join(pu_path, year, central_files[year]["Up"]))
    central_up_vals, _ = central_up[("pileup", "dense_lookup")]

        # get down variation
    central_dw = convert_histo_root_file(os.path.join(pu_path, year, central_files[year]["Down"]))
    central_dw_vals, _ = central_dw[("pileup", "dense_lookup")]


    fig, ax = plt.subplots()
    fig.subplots_adjust(hspace=.07)

        ## plot nominal vals
    ax.step(*central_bins, np.r_[central_cen_vals, central_cen_vals[-1]], where="post", **{"color" : "k", "label" : "CMS central"})
    ax.step(*ur_bins, np.r_[ur_cen_vals, ur_cen_vals[-1]], where="post", **{"color" : "k", "linestyle" : "--", "label" : "UR central"})
    ax.step(*central_bins, np.r_[central_up_vals, central_up_vals[-1]], where="post", **{"color" : "r", "label" : "CMS up"})
    ax.step(*ur_bins, np.r_[ur_up_vals, ur_up_vals[-1]], where="post", **{"color" : "r", "linestyle" : "--", "label" : "UR up"})
    ax.step(*central_bins, np.r_[central_dw_vals, central_dw_vals[-1]], where="post", **{"color" : "b", "label" : "CMS down"})
    ax.step(*ur_bins, np.r_[ur_dw_vals, ur_dw_vals[-1]], where="post", **{"color" : "b", "linestyle" : "--", "label" : "UR down"})

    ax.set_ylabel("Events")
    ax.set_xlabel("Number of Pileup Interactions")
    ax.autoscale()
    ax.legend(loc="upper right")
    ax.set_xlim(0, 100.)
    ax.set_ylim(0, ax.get_ylim()[1]*1.15)

    hep.cms.label(ax=ax, data=True, label="Preliminary", year=cfeatures.year_labels[year], lumi=round(lumi_to_use, 1))

    figname = os.path.join(outdir, f"{year}_pileup_profile")
    fig.savefig(figname)
    print(f"{figname} written")
    plt.close(fig)

    fig, ax = plt.subplots()
    fig.subplots_adjust(hspace=.07)

        ## plot nominal vals
    ax.step(*central_bins, np.r_[central_cen_vals/np.sum(central_cen_vals), (central_cen_vals/np.sum(central_cen_vals))[-1]], where="post", **{"color" : "k", "label" : "CMS central"})
    ax.step(*ur_bins, np.r_[ur_cen_vals/np.sum(ur_cen_vals), (ur_cen_vals/np.sum(ur_cen_vals))[-1]], where="post", **{"color" : "k", "linestyle" : "--", "label" : "UR central"})
    ax.step(*central_bins, np.r_[central_up_vals/np.sum(central_up_vals), (central_up_vals/np.sum(central_up_vals))[-1]], where="post", **{"color" : "r", "label" : "CMS up"})
    ax.step(*ur_bins, np.r_[ur_up_vals/np.sum(ur_up_vals), (ur_up_vals/np.sum(ur_up_vals))[-1]], where="post", **{"color" : "r", "linestyle" : "--", "label" : "UR up"})
    ax.step(*central_bins, np.r_[central_dw_vals/np.sum(central_dw_vals), (central_dw_vals/np.sum(central_dw_vals))[-1]], where="post", **{"color" : "b", "label" : "CMS down"})
    ax.step(*ur_bins, np.r_[ur_dw_vals/np.sum(ur_dw_vals), (ur_dw_vals/np.sum(ur_dw_vals))[-1]], where="post", **{"color" : "b", "linestyle" : "--", "label" : "UR down"})

    ax.set_ylabel("Probability Density")
    ax.set_xlabel("Number of Pileup Interactions")
    ax.autoscale()
    ax.legend(loc="upper right")
    ax.set_xlim(0, 100.)
    ax.set_ylim(0, ax.get_ylim()[1]*1.15)

    hep.cms.label(ax=ax, data=True, label="Preliminary", year=cfeatures.year_labels[year])# lumi=round(lumi_to_use, 1))

    figname = os.path.join(outdir, f"{year}_pileup_profile_norm")
    fig.savefig(figname)
    print(f"{figname} written")
    plt.close(fig)


toc = time.time()
print("Total time: %.1f" % (toc - tic))
