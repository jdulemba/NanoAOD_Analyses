# matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend("agg")
from matplotlib import rcParams
rcParams["font.size"] = 20
rcParams["savefig.format"] = "png"
rcParams["savefig.bbox"] = "tight"
from coffea.util import load
from pdb import set_trace
import os
import Utilities.plot_tools as plt_tools
import Utilities.prettyjson as prettyjson
import numpy as np
import Utilities.Plotter as Plotter
import Utilities.systematics as systematics
from Utilities.styles import styles as hstyles
import uproot

from argparse import ArgumentParser
parser = ArgumentParser()
#parser.add_argument("year", choices=["2018"], help="What year is the ntuple from.")
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("lepton", choices=["Electron", "Muon"], help="Choose which lepton to make plots for")
parser.add_argument("type", choices=["MC", "Folded", "Folded_LO"], help="Choose between signal produced via MC or Folding.")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
base_jobid = os.environ["base_jobid"]
jobid = os.environ["jobid"]
analyzer = "htt_btag_sb_regions"

#template_version = "v11"
template_version = "v12"
template_dname = "root://cmseos.fnal.gov//store/user/lpcbtagging/UR_ntuples/heavyhiggsinputs"
if args.type == "MC":
    output_dir = os.path.join(proj_dir, "plots", f"{args.year}_{jobid}", f"Templates_{analyzer}", "MC_SIG", template_version, args.lepton)
    template_fname = os.path.join(template_dname, f"{template_version}_smoothed/templates_lj_sig_{args.year}.root")
elif args.type == "Folded":
    output_dir = os.path.join(proj_dir, "plots", f"{args.year}_{jobid}", f"Templates_{analyzer}", "FOLDED_SIG", template_version, args.lepton)
    template_fname = os.path.join(template_dname, f"{template_version}_folded/templates_lj_sig_{args.year}.root")
elif args.type == "Folded_LO":
    output_dir = os.path.join(proj_dir, "plots", f"{args.year}_{jobid}", f"Templates_{analyzer}", "FOLDED_LO_SIG", template_version, args.lepton)
    template_fname = os.path.join(template_dname, f"{template_version}_folded_LO/templates_lj_sig_{args.year}.root")

#set_trace()
if not os.path.isdir(output_dir):
    os.makedirs(output_dir)

rfile = uproot.open(template_fname)

jet_mults = {
    "3Jets" : "3 jets",
    "4PJets" : "4+ jets"
}

leptypes = {
    "Muon" : "$\\mu$",
    "Electron" : "$e$",
}
lepdirs_dict = {
    "Muon, 3Jets" : "mu3jets",
    "Muon, 4PJets" : "mu4pjets",
    "Electron, 3Jets" : "e3jets",
    "Electron, 4PJets" : "e4pjets",
}

data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())[args.year]

def proc_to_name(proc):
    new_proc = proc.replace("A", "AtoTTJetsSL_").replace("H", "HtoTTJetsSL_").replace("relw", "W").replace("int", "Int").replace("res", "Res") # get almost correct format except for mass
    mass = new_proc.split("_")[1]
    new_proc = new_proc.replace(mass, "M"+mass)
    return new_proc
    

#set_trace()
    ## make plots for background templates
#for jmult in ["4PJets"]:
for jmult in ["3Jets", "4PJets"]:
    dirname = lepdirs_dict[f"{args.lepton}, {jmult}"]
        # find possible signal points (ignoring error distributions)
    sig_dists = sorted(set([(key.split("/")[1]).strip(";1") for key in rfile.keys() if ((key.split("/")[0] == dirname))]))# and ("error" not in key[0].split("/")[1]))]))

    #set_trace()
    systypes = systematics.sys_groups[args.year].keys()
    #systypes = ["BTAG_L_CORR"]
    #systypes = ["PILEUP"]
    for sys in systypes:
            # find histograms of associated systematics and their processes
        up_sysname = systematics.template_sys_to_name[args.year][systematics.sys_groups[args.year][sys][0]] if systematics.sys_groups[args.year][sys][0] in systematics.template_sys_to_name[args.year].keys() else None
        dw_sysname = systematics.template_sys_to_name[args.year][systematics.sys_groups[args.year][sys][1]] if systematics.sys_groups[args.year][sys][1] in systematics.template_sys_to_name[args.year].keys() else None

        if (not up_sysname) and (not dw_sysname):
            print(f"No match found for {sys}. Skipping")
            continue

        baseSysname = up_sysname.split("Up")[0] if up_sysname else dw_sysname.split("Down")[0]
        procs_sys = sorted(set([key.split(f"_{baseSysname}")[0] for key in sig_dists if baseSysname in key]))

        if not procs_sys: continue

        for proc in procs_sys:
            if ("A500_relw10p0" not in proc): continue
            #if (proc != "A1000_relw25p0_res"): continue

            pltdir = os.path.join(output_dir, jmult, "Comp", sys)
            if not os.path.isdir(pltdir):
                os.makedirs(pltdir)

            print(jmult, sys, proc)
            #set_trace()

            bins = [ax.edges() for ax in rfile[f"{dirname}/{proc}_{up_sysname}_orig"].axes][0]
            orig_up = rfile[f"{dirname}/{proc}_{up_sysname}_orig"].values()
            orig_dw = rfile[f"{dirname}/{proc}_{dw_sysname}_orig"].values()
            smooth_up = rfile[f"{dirname}/{proc}_{up_sysname}SMOOTH"].values()
            smooth_dw = rfile[f"{dirname}/{proc}_{dw_sysname}SMOOTH"].values()
            final_up = rfile[f"{dirname}/{proc}_{up_sysname}FINAL"].values()
            final_dw = rfile[f"{dirname}/{proc}_{dw_sysname}FINAL"].values()

            up_histos = [
                    (orig_up-1, {"color": "r", "linestyle": "--", "label": "Original Up"}, True),
                    (smooth_up-1, {"color": "r", "linestyle": "-", "label": "Smooth Up"}, False),
                    (final_up-1, {"color": "r", "linestyle": "--", "label": "Final Up"}, False),
            ]
            dw_histos = [
                    (orig_dw-1, {"color": "b", "linestyle": "--", "label": "Original Down"}, True),
                    (smooth_dw-1, {"color": "b", "linestyle": "-", "label": "Smooth Down"}, False),
                    (final_dw-1, {"color": "b", "linestyle": "--", "label": "Final Down"}, False),
            ]

            x_lims = (0, bins[-1])
    
            fig, ax = plt.subplots()
            fig.subplots_adjust(hspace=.07)

                ## plot relative deviations
            for up_histo, up_style, use_fill_between in up_histos:
                    # there is at least one actual value
                if np.any(~np.isnan(up_histo)):
                    ax.fill_between(bins, np.r_[up_histo, up_histo[-1]], facecolor=up_style["color"], step="post", alpha=0.5) if use_fill_between else ax.step(bins, np.r_[up_histo, up_histo[-1]], where="post", **up_style)

            for dw_histo, dw_style, use_fill_between in dw_histos:
                    # there is at least one actual value
                if np.any(~np.isnan(dw_histo)):
                    ax.fill_between(bins, np.r_[dw_histo, dw_histo[-1]], facecolor=dw_style["color"], step="post", alpha=0.5) if use_fill_between else ax.step(bins, np.r_[dw_histo, dw_histo[-1]], where="post", **dw_style)

            #set_trace()
            ax.legend(loc="upper right", title=f"{sys}, {plt_tools.get_label(proc_to_name(proc), hstyles)}")
            ax.axhline(0, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
            ax.autoscale()
            ax.set_ylim(max(ax.get_ylim()[0], -0.2), min(ax.get_ylim()[1], 0.2))
            ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.15)
            ax.set_xlim(x_lims)
            ax.set_xlabel("$m_{t\\bar{t}}$ $\otimes$ |cos($\\theta^{*}_{t_{l}}$)|")
            ax.set_ylabel("Rel. Deviaton from Nominal")
            
                # add lepton/jet multiplicity label
            ax.text(
                0.02, 0.88, f"{leptypes[args.lepton]}, {jet_mults[jmult]}\n{(args.type).replace('_', ' ')}",
                fontsize=rcParams["font.size"]*0.9, horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
            )
                ## draw vertical lines for distinguishing different ctstar bins
            vlines = [x_lims[1]*ybin/5 for ybin in range(1, 5)]
            for vline in vlines:
                ax.axvline(vline, color="k", linestyle="--")
            hep.cms.label(ax=ax, data=False, year=args.year, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))
            
            #set_trace()
            figname = os.path.join(pltdir, "_".join([jmult, args.lepton, sys, proc, "SysTemplates_Comp"]))
            fig.savefig(figname)
            print(f"{figname} written")
            plt.close()
