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
import Utilities.prettyjson as prettyjson
import numpy as np
import Utilities.Plotter as Plotter
import Utilities.final_analysis_binning as final_binning
import Utilities.common_features as cfeatures
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
args = parser.parse_args()

year_to_use = cfeatures.year_labels[args.year]

base_jobid = os.environ["base_jobid"]
proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
analyzer = "htt_btag_sb_regions"
plot_outdir = os.environ["plots_dir"]
eos_dir = os.environ["eos_dir"]

#set_trace()
print("Choose version")
version = "V31"

outdir = os.path.join(plot_outdir, jobid, f"Templates_{analyzer}", version, args.year,  "Singlet_Over_TTbar")#, args.lepton)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

linearize_binning = (
    final_binning.mtt_binning,
    final_binning.ctstar_abs_binning
)

ctstar_binlabels = [r"|$\cos (\theta^{*}_{t_{l}})$| $\in [%s, %s)$" % (linearize_binning[1][bin], linearize_binning[1][bin+1]) for bin in range(len(linearize_binning[1])-1)]
ctstar_binlabels[-1] = r"|$\cos (\theta^{*}_{t_{l}})$| $\in [%s, %s]$" % (linearize_binning[1][-2], linearize_binning[1][-1])
ctstar_bin_locs = np.linspace((len(linearize_binning[0])-1)/2, (len(linearize_binning[0])-1)*(len(linearize_binning[1])-1) - (len(linearize_binning[0])-1)/2, len(linearize_binning[1])-1)
mtt_vals_to_plot = np.array([400, 600, 1000])
mtt_tiled_labels = np.tile(linearize_binning[0][:-1], linearize_binning[1].size-1)
mtt_bin_inds_to_plot = np.where(np.in1d(mtt_tiled_labels, mtt_vals_to_plot))[0]
mtt_bins_to_plot = np.tile(mtt_vals_to_plot, linearize_binning[1].size-1)

input_dir = os.path.join(eos_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}")
final_hdict = load(os.path.join(input_dir, "FINAL", f"final_templates_lj_bkg_{args.year}_{jobid}_{version}.coffea"))

data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())[args.year]
lumi_to_use = (data_lumi_year["Muons"] + data_lumi_year["Electrons"])/2000.

ratio_fig, ratio_ax = plt.subplots(figsize=(15.0, 10.0))
ratio_fig.subplots_adjust(hspace=.07)

colors = {"Muon_3Jets" : "k", "Muon_4PJets" : "r", "Electron_3Jets" : "b", "Electron_4PJets" : "g"}

    ## make plots for background templates
for jmult in final_hdict.keys():
    for lep, hdict in final_hdict[jmult].items():
        #set_trace()
        TQ_nosys, TW_nosys, TB_nosys = hdict["TQ_nosys"][0].copy(), hdict["TW_nosys"][0].copy(), hdict["TB_nosys"][0].copy()
        ST_nosys_vals = TQ_nosys.values()[()] + TW_nosys.values()[()] + TB_nosys.values()[()]
        TT_nosys = hdict["TT_nosys"][0].copy()
        x_lims = (0, TT_nosys.dense_axes()[0].centers().size)

        print(lep, jmult)
        #set_trace()
        masked_ratio_vals, masked_ratio_bins = Plotter.get_ratio_arrays(num_vals=ST_nosys_vals, denom_vals=TT_nosys.values()[()], input_bins=TT_nosys.dense_axes()[0].edges())
        ratio_ax.step(masked_ratio_bins, masked_ratio_vals, where="post", **{"color": colors[f"{lep}_{jmult}"], "linestyle":"-", "label": cfeatures.channel_labels[f"{lep}_{jmult}"]})

#set_trace()
ratio_ax.legend(loc="upper right", title="Channel", ncol=2)                
ratio_ax.autoscale()
ratio_ax.set_xlim(x_lims)
ratio_ax.set_xlabel("$m_{t\\bar{t}}$ [GeV]")
ratio_ax.set_ylabel("single t/$t\\bar{t}$")

    ## draw vertical lines for distinguishing different ctstar bins
vlines = [x_lims[1]*ybin/5 for ybin in range(1, 5)]
[ratio_ax.axvline(vline, color="k", linestyle="--") for vline in vlines]

for idx, label in enumerate(ctstar_binlabels):
    ratio_ax.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
        xytext=(0, 10), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)

ratio_ax.set_xticks(mtt_bin_inds_to_plot)
ratio_ax.set_xticklabels(mtt_bins_to_plot)

#set_trace()
hep.cms.label(ax=ratio_ax, data=False, label="Preliminary", year=year_to_use, lumi=round(lumi_to_use, 1))

ratio_figname = os.path.join(outdir, f"Singlet_To_TTbar_Template_Ratio_{args.year}")
ratio_fig.savefig(ratio_figname)
print(f"{ratio_figname} written")
plt.close(ratio_fig)

toc = time.time()
print("Total time: %.1f" % (toc - tic))
