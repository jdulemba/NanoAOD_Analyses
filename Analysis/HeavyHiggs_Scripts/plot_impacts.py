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
from matplotlib import ticker

import numpy as np
from pdb import set_trace
import uproot
import os

import argparse
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("parity", choices=["A", "H"], help="Choose which parity to make plots for")
parser.add_argument("mass", choices=["400", "800"], help="Choose which mass to make plots for")
parser.add_argument("width", choices=["5.0"], help="Choose which width to make plots for")
parser.add_argument("--nvals", default=20, help="Choose the number of NPs to plot per figure")
parser.add_argument("--unblind", action="store_true", help="Make plots with best fit coupling value")
args = parser.parse_args()

jobid = os.environ["jobid"]
proj_dir = os.environ["PROJECT_DIR"]
plot_outdir = os.environ["plots_dir"]

sig_label = "$%s(%s\ GeV, %s\%%)$" % (args.parity, args.mass, args.width)

def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


basedir = "root://cmseos.fnal.gov//store/user/jdulemba/HeavyHiggsFitter/Summer20UL_DeepJet_Unblinded_Toponium_v35"
rname = f"IMPACTSv35_lj_data_1_{args.parity}_m{args.mass}_relw{args.width.replace('.', 'p')}_sig1.root"
rfile = uproot.open(os.path.join(basedir, rname))

outdir = os.path.join(plot_outdir, jobid, "Limits_htt_btag_sb_regions", os.path.basename(basedir).split(f"{jobid}_")[-1], f"{args.parity}_m{args.mass}_relw{args.width.replace('.', 'p')}")
if not os.path.isdir(outdir):
    os.makedirs(outdir)


hup   = rfile["IMPACTup"]
hdown = rfile["IMPACTdown"]
huncdown = rfile["UNCSDOWN"]
huncup = rfile["UNCSUP"]
hunccentral = rfile["UNCSCENTRAL"]

unc_labels = hup.axes[0].labels()

sigbin = np.argmax(np.abs(hup.values()))

    ## find best fit g
mu_cen = abs(hunccentral.values()[sigbin])
mu_up = huncup.values()[sigbin]
mu_down = huncdown.values()[sigbin]

## make list of [ [nuisance par name, impact_up, impact_dw, post-fit mean, post-fit unc down, post-fit unc up] ]
    # get indices corresponding to NPs that aren't defined (dilepton NPs), "double-counted", plus the sigbin)
inds_to_ignore = list(np.where(hunccentral.values() == 0.)[0]) + [idx for idx, name  in enumerate(unc_labels) if "tmass_AH" in name] + [sigbin]
info_list = [[name, hup.values()[idx], hdown.values()[idx], hunccentral.values()[idx], huncdown.values()[idx], huncup.values()[idx]] for idx, name in enumerate(unc_labels) if idx not in inds_to_ignore]
sorted_info_list = sorted(info_list, key = lambda x:  max(abs(x[1]), abs(x[2])), reverse = True)
    # get impact names, up values, and down values
impact_labels, impact_up, impact_dw = [val[0] for val in sorted_info_list], np.array([val[1] for val in sorted_info_list]), np.array([val[2] for val in sorted_info_list])
impact_labels = [label.replace("UNCORR", "").replace("_CMS_", "_").replace("CMS_", "").replace("_13TeV_", "_").replace("13TeV_", "").replace("tmass_TT", "tmass_AHTT") for label in impact_labels]
impact_labels = ["%s    (%i)" % (label, idx+1) for idx, label in enumerate(impact_labels)]

    ## calculate pulls
pulls = []
for val in sorted_info_list:
    if "tmass_TT" in val[0]:
        val[3] *= 3.
        val[4] *= 3.
        val[5] *= 3.
    pulls.append( val[3] / np.sqrt( 1. - 0.5 * ( np.square( val[4]) + np.square( val[5])) ) )
pulls = np.array(pulls)

    ## calculate post-fit values
pf_mean, pf_unc_dw, pf_unc_up = np.array([val[3] for val in sorted_info_list]), np.array([val[4] for val in sorted_info_list]), np.array([val[5] for val in sorted_info_list])


nplots = list(chunks(range(len(impact_labels)), int(args.nvals)))
## make plots
for idx, nplot in enumerate(nplots):
    #if idx > 0: continue
    idx_min, idx_max = nplot[0], nplot[-1]+1

    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
    fig.subplots_adjust(hspace=0.0, wspace=0.0)
    
    ## plot up and down impacts
    ax1.barh(
        range(len(impact_up[idx_min:idx_max])), impact_up[idx_min:idx_max][::-1],
        alpha=0.5, linestyle="solid", facecolor="r", height=1., align="edge", label="+1$\\sigma$ Impact",
    )
    ax1.barh(
        range(len(impact_dw[idx_min:idx_max])), impact_dw[idx_min:idx_max][::-1],
        alpha=0.5, linestyle="solid", facecolor="b", height=1., align="edge", label="-1$\\sigma$ Impact",
    )
    ax1.autoscale()
    ax1.grid(linestyle="dotted", color="k", axis="x", linewidth=2)
    ax1.set_ylim(0., len(impact_up[idx_min:idx_max]))
    ax1.set_xlabel("$\\Delta \\hat{g}_{%s}$" % args.parity, horizontalalignment="center")
    ax1.set_yticks(np.arange(0.5, len(impact_up[idx_min:idx_max]), step=1.), impact_labels[idx_min:idx_max][::-1], minor=True, ha="right")
    ax1.yaxis.set_major_locator(ticker.MultipleLocator(1.))
    ax1.yaxis.set_major_formatter(ticker.NullFormatter()) # removes numbers
    ax1.tick_params(axis="y", which="minor", length=0)
    ax1.text(-0.05, 1., sig_label, fontsize=rcParams["font.size"]*1.3, horizontalalignment="right", verticalalignment="bottom", transform=ax1.transAxes) ## align title above yaxis values

    hep.cms.label(ax=ax1, data=True, label="Preliminary", rlabel="")
    
    ## plot pulls and post-fit uncertainty
    ax2.scatter(x=pulls[idx_min:idx_max][::-1], y=np.arange(0.3, len(pulls[idx_min:idx_max]), step=1.), s=100, color="b", marker="x", label="Pull")
    ax2.errorbar(x=pf_mean[idx_min:idx_max][::-1], y=np.arange(0.5, len(pf_mean[idx_min:idx_max]), step=1.), xerr=[np.abs(pf_unc_dw)[idx_min:idx_max][::-1], np.abs(pf_unc_up)[idx_min:idx_max][::-1]],
        fmt="ok", capsize=7, ecolor="k", label="Post-fit $\\hat{\\theta}$ and $\\hat{\\sigma}$"
    )
        ## add axis labels
    for label, color, position in zip(["$(\\hat{\\theta} - \\theta_{i})/\\sigma_{i}$", "$(\\hat{\\theta} - \\theta_{i})/\\sqrt{\\sigma_{i}^2 - \\hat{\\sigma}^2}$"], ["k", "b"], [0, 0.4]):
        ax2.text(position, -0.1, label, color=color, transform=ax2.transAxes, verticalalignment="bottom", horizontalalignment="left")

    [ax2.axvline(vline, color="k", linestyle="dotted", linewidth=2) for vline in np.arange(np.floor(ax2.get_xlim()[0]), np.ceil(ax2.get_xlim()[1])+1, step=1)]
    ax2.tick_params(axis="y", which="minor", length=0)

        ## add best fit g value
    if args.unblind:
        ax2.text(0.3, 1.01, "$\\hat{g}_{%s} = %.3f^{+%.3f}_{%.3f}$" % (args.parity, mu_cen, mu_up, mu_down), horizontalalignment="left", verticalalignment="bottom", transform=ax2.transAxes)

        ## combine labels from both axes into one legend
    ax1_handles, ax1_labels = ax1.get_legend_handles_labels()
    ax2_handles, ax2_labels = ax2.get_legend_handles_labels()
    lgd = fig.legend(ax1_handles+ax2_handles, ax1_labels+ax2_labels, loc="lower center", ncol=2, bbox_to_anchor=(-0.05, 0.))
    [hand.set_ec("k") for idx, hand in enumerate(lgd.legend_handles) if idx < 2] ## only set edgecolor for impacts


    outfname = os.path.join(outdir, f"{rname.replace('.root', '')}_{idx}")
    fig.savefig(outfname)
    print(f"{outfname} written")
    plt.close(fig)

toc = time.time()
print("Total time: %.1f" % (toc - tic))
