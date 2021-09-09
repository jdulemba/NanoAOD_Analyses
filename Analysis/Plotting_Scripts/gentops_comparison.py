# matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend("agg")
from matplotlib import rcParams
rcParams["font.size"] = 20
rcParams["savefig.format"] = "png"
rcParams["savefig.bbox"] = "tight"

from coffea.hist import plot
from coffea.util import load
from pdb import set_trace
import os
import re, fnmatch
import sys
import Utilities.plot_tools as plt_tools
import Utilities.prettyjson as prettyjson
import Utilities.Plotter as Plotter
from Utilities.styles import styles as styles

proj_dir = os.environ["PROJECT_DIR"]
base_jobid = os.environ["base_jobid"]
analyzer = "gentops_comparison"

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016", "2017", "2018"] if base_jobid == "NanoAODv6" else ["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("--indiv_int", action="store_true", help="Plot individual positive and negative weights for Interference samples.")
parser.add_argument("--wildcard", default="*", help="Make plots for selected datasets using wildcard given ('*' is default).")
args = parser.parse_args()

input_dir = os.path.join(proj_dir, "results", "%s_%s" % (args.year, base_jobid), analyzer)
f_ext = "TOT.coffea"
outdir = os.path.join(proj_dir, "plots", "%s_%s" % (args.year, base_jobid), analyzer)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

#set_trace()
fnames = sorted(["%s/%s" % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)])
hdict = plt_tools.add_coffea_files(fnames) if len(fnames) > 1 else load(fnames[0])

objects = {
    "TTbar" : ("t\\bar{t}", (200., 2000.)),
    "Top" : ("t", (150., 200.)),
    "Tbar" : ("\\bar{t}", (150., 200.)),
}

isSignal = lambda x : (x.startswith("AtoTT") or x.startswith("HtoTT"))

variables = {
    "pt" : ("$p_{T}$($obj$) [GeV]", 2, (0., 500.)),
    "eta": ("$\\eta$($obj$)", 2, (-2.6, 2.6)),
    "phi": ("$\\phi$($obj$)", 2, (-4, 4)),
    "energy": ("$E_{obj}$ [GeV]", 2, (0., 2000.)),
    "mass": ("$m_{obj}$ [GeV]", 1, (0., 300.)),
    "ctstar" : ("cos($\\theta^{*}_{obj}$)", 1, (-1., 1.)),
    "ctstar_abs" : ("|cos($\\theta^{*}_{obj}$)|", 1, (0., 1.)),
}


    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, "inputs", "%s_lumis_data.json" % base_jobid)).read())[args.year]
lumi_to_use = (data_lumi_year["Muons"]+data_lumi_year["Electrons"])/2000.
lumi_correction = load(os.path.join(proj_dir, "Corrections", base_jobid, "MC_LumiWeights.coffea"))[args.year]

        # scale events by lumi correction
for hname in hdict.keys():
    if hname == "cutflow": continue
    hdict[hname].scale(lumi_correction["Muons"], axis="dataset")


    ## make bp plots
for hname in variables.keys():
    if hname not in hdict.keys():
        raise ValueError(f"{hname} not found in file")
    print(f"{hname}")
    h_tot = hdict[hname]
    xtitle, rebinning, x_lims = variables[hname]
    xaxis_name = h_tot.dense_axes()[0].name
    if rebinning != 1:
        h_tot = h_tot.rebin(xaxis_name, rebinning)

    all_possible_datasets = sorted(set([key[0] for key in h_tot.values().keys()]))
    datasets_to_use = all_possible_datasets if args.indiv_int else \
        ["_".join(signal.split("_")[:-1]) for signal in all_possible_datasets if "neg" in signal] + [signal for signal in all_possible_datasets if "Res" in signal] + [signal for signal in all_possible_datasets if "ttJets" in signal]

        # choose which datasets to make plots for based on wildcard
    datasets_to_use = [sample for sample in datasets_to_use if fnmatch.fnmatch(sample, args.wildcard)]
    if not datasets_to_use:
        print(f"No matches found for wildcard {args.wildcard}")
        sys.exit()

    for dataset in sorted(datasets_to_use):
        print(f"\t{dataset}")
        pltdir = os.path.join(outdir, "_".join(dataset.split("_")[:-1]), dataset.split("_")[-1]) if (("Int" in dataset) and (args.indiv_int)) else os.path.join(outdir, dataset)
        if not os.path.isdir(pltdir):
            os.makedirs(pltdir)

        decay_label = plt_tools.get_label(dataset, styles)
        dataset_mask = re.compile(r"((?:%s*))" % dataset)
        histo = h_tot[dataset_mask]
        #set_trace()
        for genobj, (objlabel, mass_range) in objects.items():
            if genobj not in sorted(set([key[1] for key in histo.values().keys()])): continue
            new_xtitle = xtitle.replace("obj", objlabel)
            if hname == "mass":
                x_lims = mass_range
            tt_histo = histo[:, genobj].integrate("objtype")

            fig, ax = plt.subplots()
            fig.subplots_adjust(hspace=.07)

            if ("Int" in dataset) and (not args.indiv_int):
                    # plot negative weights first
                plot.plot1d(
                    tt_histo[re.compile(r".*(neg)")].integrate("dataset"),
                    overlay=tt_histo.integrate("dataset").axes()[0].name,
                    line_opts={"color" : ["#e41a1c", "#e41a1c"], "linestyle" : ["-", "--"]}, # red
                    ax=ax, clear=False, stack=False,
                )
                    # plot positive weights
                plot.plot1d(
                    tt_histo[re.compile(r".*(pos)")].integrate("dataset"),
                    overlay=tt_histo.integrate("dataset").axes()[0].name,
                    line_opts={"color" : ["#377eb8", "#377eb8"], "linestyle" : ["-", "--"]}, # blue
                    ax=ax, clear=False, stack=False,
                )
                    # plot combined
                plot.plot1d(
                    tt_histo.integrate("dataset"),
                    overlay=tt_histo.integrate("dataset").axes()[0].name,
                    line_opts={"color" : ["k", "k"], "linestyle" : ["-", "--"]}, # black
                    ax=ax, clear=False, stack=False,
                )

                    # set legend labels by hand !! BE CAREFUL !!
                handles, labels = ax.get_legend_handles_labels()
                for idx, sample in enumerate(labels):
                    if (idx == 0) or (idx == 1):
                        labels[idx] = f"w$<$0, {sample}"
                    if (idx == 2) or (idx == 3):
                        labels[idx] = f"w$>$0, {sample}"
                    if (idx == 4) or (idx == 5):
                        labels[idx] = f"Combined, {sample}"
                ax.legend(handles, labels, title=decay_label, loc="upper right", title_fontsize=24)
                ax.axhline(y=0., color="k", linestyle="-")

            else:
                plot.plot1d(
                    tt_histo.integrate("dataset"),
                    overlay=tt_histo.integrate("dataset").axes()[0].name,
                    line_opts={"color" : ["#e41a1c", "k"]}, # red and black
                    ax=ax, clear=False, stack=False,
                )

                ax.legend(title=decay_label, loc="upper right", title_fontsize=24)

            ax.autoscale()
            ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.3)
            ax.set_xlabel(new_xtitle)
            ax.set_xlim(x_lims)
            hep.cms.label(ax=ax, data=False, paper=False, year=args.year, lumi=round(lumi_to_use, 1))

            #set_trace()
            figname = os.path.join(pltdir, "_".join([dataset, genobj, hname, "statusFlagsComp"]))
            fig.savefig(figname)
            print(f"{figname} written")
            plt.close()
