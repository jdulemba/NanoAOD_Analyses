#! /bin/env python
import time
tic = time.time()

import math
from pdb import set_trace
import numpy as np

# matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend("agg")
from matplotlib import rcParams
rcParams["font.size"] = 26
#rcParams["font.size"] = 20
rcParams["savefig.format"] = "pdf"
rcParams["savefig.bbox"] = "tight"

import matplotlib
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.ticker as ticker
import os
import uproot

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("channel", choices=["Run2", "2016", "2017", "2018", "3j", "4pj"], help="Choose which channel to make plots for.")
parser.add_argument("--widths", type=str, default="All", help="Choose which width points to compute limits for, multiple options can be input as ':' separated strings.")
parser.add_argument("--parity", type=str, default="All", help="Choose which parity to compute limits for, multiple options can be input as ':' separated strings.")
args = parser.parse_args()


# Stolen from Andrey
def max_g(cp, phi_mass, rel_width):
	"""Compute maximal allowed value of the coupling scale factor.

	Computed value corresponds to a 100% branching ratio of H->tt.

	Arguments:
		cp:  CP state, 'A' or 'H'.
		phi_mass:  Mass of the Higgs boson, GeV.
		width:  Relative width of the Higgs boson.
	"""

	gF = 1.1663787e-5  # GeV^(-2)
	mt = 172.5  # GeV

	if phi_mass <= 2 * mt:
		return 0.

	w = 3 * gF * mt ** 2 * phi_mass / (4 * math.pi * math.sqrt(2))
	beta = math.sqrt(1 - (2 * mt / phi_mass) ** 2)

	if cp == 'A':
		width_g1 = w * beta
	elif cp == 'H':
		width_g1 = w * beta ** 3
	else:
		raise RuntimeError('Cannot recognize CP state "{}".'.format(cp))

	return math.sqrt(rel_width * phi_mass / width_g1)



def compute_limit(parity, relw):
    fname = os.path.join(basedir, f"limit_{rfile_channel_name}_{parity}_{val2name(relw)}.root")
    #fname = os.path.join(basedir, f"limit_ASYMCRABv35_lj_data_5_{parity}_{val2name(relw)}.root")
    try:
        rfile = uproot.open(fname)
    except:
        print(f"{fname} not found")
        #continue
        #raise ValueError(f"{fname} not found")

    masses = np.copy(rfile["central"].values()[0])
    x_min, x_max = masses.min(), masses.max()
    maxg_values = [(mass, max_g(parity, mass, relw/100.)) for mass in np.arange(min(masses), max(masses)+5., 5.)]

    obs_data = np.array(np.copy(rfile["data"].values()[1]), dtype="f4")
    obs_exp = np.array(np.copy(rfile["central"].values()[1]), dtype="f4")
    pos1_exp = np.array(obs_exp + np.copy(rfile["one_sigma"].errors("high")[1]), dtype="f4")
    pos2_exp = np.array(obs_exp + np.copy(rfile["two_sigma"].errors("high")[1]), dtype="f4")
    neg1_exp = np.array(obs_exp - np.copy(rfile["one_sigma"].errors("low")[1]), dtype="f4")
    neg2_exp = np.array(obs_exp - np.copy(rfile["two_sigma"].errors("low")[1]), dtype="f4")

    max_coupling = np.max(np.concatenate((obs_data, pos2_exp)))
    max_plot_coupling = np.round(max_coupling*1.1, decimals=1)

    # obs_color = (103./255., 203./255., 123./255., 0.4)
    obs_color = (135./255., 206./255., 250./255., 0.5)
    
    fig = plt.figure(figsize=(10, 10), dpi=80, facecolor="w", edgecolor="k")
    #fig = plt.figure(figsize=(15, 10), dpi=80, facecolor="w", edgecolor="k")
    ax = fig.add_subplot(111)

    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.05))
    handles = []

        # plot expected    
    center = ax.plot(masses, obs_exp, color="k", linestyle="-")
    handles.append( (mlines.Line2D([], [], color="k", linestyle="-"), "Expected") )

    #set_trace()
    twosig = ax.fill_between(masses, neg2_exp, pos2_exp, color=twosigma)
    #handles.append( (mpatches.Patch(color=twosigma), r"95% Expected") )

    onesig = ax.fill_between(masses, pos1_exp, neg1_exp, color=onesigma)
    handles.append( (mpatches.Patch(color=onesigma), "$\\pm 1 \\sigma$") )
    handles.append( (mpatches.Patch(color=twosigma), "$\\pm 2 \\sigma$") )
    #handles.append( (mpatches.Patch(color=onesigma), r"68% expected") )
    
     # plot data
    observed = ax.plot(masses, obs_data, color="r", linestyle="-")
    #set_trace()
    #ax.fill_between(masses, obs_data, [val+0.05 for val in obs_data], color="none", hatch="xx", edgecolor="r", linewidth=0.)
    handles.append((mpatches.Patch(facecolor="none", hatch="xx", edgecolor="r", linewidth=1.), "data"))
    

        # set labels
    ax.set_xlabel(addenda["mass"] % (parity, parity, relw)+", "+xlabels["mass"] % parity, fontsize=32)
    ax.set_ylabel("$g_{%s \\rightarrow t\\bar t}$" % parity, fontsize=32)

        # plot unphysical region    
    #set_trace()
    handles.append((mpatches.Patch(color=obs_color), "$\Gamma_{%s \\rightarrow t\\bar t} > \Gamma_{%s}$" % (parity, parity)))
    maxg_xvalues_todraw = [val[0]for val in maxg_values if val[1] < max_plot_coupling]
    maxg_values_todraw = [val[1] for val in maxg_values if val[1] < max_plot_coupling]
    unphys_region = ax.plot(maxg_xvalues_todraw, maxg_values_todraw, color=obs_color, linestyle="-")
    ax.fill_between(maxg_xvalues_todraw, maxg_values_todraw, max_plot_coupling, color=obs_color)
    
    #legend
    ax.axhline(y=max_plot_coupling, color="k", linestyle="-")
    ax.autoscale()
    ax.set_ylim(0., ax.get_ylim()[1]*1.35)
    #ax.set_ylim(0., ax.get_ylim()[1]*1.4)
    ax.set_xlim(x_min, x_max)
    lheight = max_plot_coupling/ax.get_ylim()[1]
    lmargin = 0.02

        # plot hatching for data
    ax.fill_between(masses, obs_data, [val + max_plot_coupling*0.05 for val in obs_data], color="none", hatch="xx", edgecolor="r", linewidth=0.)

        # plot legend
    legend = ax.legend(
        x(handles), y(handles), loc="lower left", ncol=2, fontsize=29, frameon=False,
        bbox_to_anchor=(lmargin, lheight),
        borderpad=0., labelspacing=0.2,
    )
    fontP = matplotlib.font_manager.FontProperties()
    fontP.set_size(29)
    legend.set_title(title=type_label, prop=fontP)


    #CMS blurb
    hep.cms.label(ax=ax, data=True, label="Preliminary", year=year_label, lumi=lumi_to_use)
    
    ax.tick_params(axis="both", labelsize=29, which="both")

    #set_trace()
    figname = os.path.join(outdir, f"limit_{parity}_{val2name(relw)}")
    fig.savefig(figname)
    print(f"{figname} written")
    plt.close(fig)



if __name__ == "__main__":	
    proj_dir = os.environ["PROJECT_DIR"]
    base_jobid = os.environ["base_jobid"]
    jobid = os.environ["jobid"]
    analyzer = "htt_btag_sb_regions"
    user = os.environ["USER"]
    eos_dir = os.environ["eos_dir"]
    plot_dir = os.environ["plots_dir"]

    #set_trace()
    basedir = f"root://cmseos.fnal.gov//store/user/jdulemba/HeavyHiggsFitter/Thesis_constrainedYt_EtaT_v1/limits_{args.channel}"
    outdir = os.path.join(plot_dir, jobid, "Limits_htt_btag_sb_regions", basedir.split("/")[-2], "Limits", args.channel)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    #possible_widths = [25.0]
    ##possible_widths = [1.0, 25.0]
    possible_widths = [1.0, 2.0, 2.5, 3.0, 4.0, 5.0, 8.0, 10.0, 13.0, 15.0, 18.0, 21.0, 25.0]
    possible_widths = [str(width) for width in possible_widths]

    widths_to_run = args.widths
    if widths_to_run == "All":
        available_widths = possible_widths
    else:
        available_widths = widths_to_run.split(":")
    available_widths = [float(width) for width in available_widths if width in possible_widths]
    
    possible_parities = ["A", "H"]
    parities_to_run = args.parity
    if parities_to_run == "All":
        available_parities = possible_parities
    else:
        available_parities = parities_to_run.split(":")
        available_parities = [parity for parity in available_parities if parity in possible_parities]

    xlabels = {
    	"mass" : "$m_{%s}$ [GeV]",
    	"width": "$\Gamma_{%s}$/$m_{%s}$\%%",
    }
    
    addenda = {
    	"width": "$m_{%s}$ = %d [GeV]",
    	"mass" : r"$\Gamma_{%s}$/$m_{%s}$ = %.1f$\mathsf{\%%}$",
    	}
    
    vartoadd = {
    	"width" : "mass",
    	"mass" : "width"
    }
    val2name = lambda x: ("%.1f" % x).replace(".", "p")

    x = lambda vv: [i for i, _ in vv]
    y = lambda vv: [i for _, i in vv]

    if (args.channel == "Run2") or (args.channel == "3j") or (args.channel == "4pj"):
        year_label = "Run 2"
        lumi_to_use = 138
        if args.channel == "Run2":
            rfile_channel_name = "ASYMCRABv35_lj_data_5"
            type_label = "95% $\mathsf{CL_s}$ limits, $\\ell$+jets"
        elif args.channel == "3j":
            rfile_channel_name = "ASYMCRABv35_lj3j_data_1"
            type_label = "95% $\mathsf{CL_s}$ limits, $\\ell$/3 jets"
        elif args.channel == "4pj":
            rfile_channel_name = "ASYMCRABv35_lj4j_data_1"
            type_label = "95% $\mathsf{CL_s}$ limits, $\\ell$/4+ jets"
    else:
        import Utilities.prettyjson as prettyjson
        data_lumi_file = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())
        if args.channel == "2016":
            lumi_to_use = round((data_lumi_file["2016APV"]["Muons"] + data_lumi_file["2016APV"]["Electrons"] + data_lumi_file["2016"]["Muons"] + data_lumi_file["2016"]["Electrons"])/2000., 1)
            rfile_channel_name = "ASYMCRABv35_lj16_data_1"
        elif args.channel == "2017":
            rfile_channel_name = "ASYMCRABv35_lj17_data_1"
            lumi_to_use = round((data_lumi_file[args.channel]["Muons"] + data_lumi_file[args.channel]["Electrons"])/2000., 1)
        elif args.channel == "2018":
            rfile_channel_name = "ASYMCRABv35_lj18_data_1"
            lumi_to_use = round((data_lumi_file[args.channel]["Muons"] + data_lumi_file[args.channel]["Electrons"])/2000., 1)

        type_label = "95% $\mathsf{CL_s}$ limits, $\\ell$+jets"
        year_label = args.channel

    onesigma = "#00cc00" # kGreen + 1
    twosigma = "#ffcc00" # kOrange
    #onesigma = '#00f847'
    #twosigma = '#fffc4d'
    line = "#ff1521"

    for parity in available_parities:
        for relw in available_widths:
            compute_limit(parity, relw)

toc = time.time()
print("Total time: %.1f" % (toc - tic))
