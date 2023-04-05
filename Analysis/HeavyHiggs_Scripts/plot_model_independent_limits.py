#! /bin/env python
import math
from pdb import set_trace
import numpy as np

# matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend("agg")
from matplotlib import rcParams
rcParams["font.size"] = 20
rcParams["savefig.format"] = "pdf"
rcParams["savefig.bbox"] = "tight"

import matplotlib
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.ticker as ticker
import os
from coffea.util import load

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("type", choices=["MC", "MEreweighting", "Morphed"], help="Choose between signal produced via MC, ME reweighting, or Folding.")
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


def make_plot(subset, parity, masses, width, maxg_values=None):
    #set_trace()
    pos2_exp = np.array([val[0] for val in subset.values()], dtype="f4")
    pos1_exp = np.array([val[1] for val in subset.values()], dtype="f4")
    obs_exp = np.array([val[2] for val in subset.values()], dtype="f4")
    neg1_exp = np.array([val[3] for val in subset.values()], dtype="f4")
    neg2_exp = np.array([val[4] for val in subset.values()], dtype="f4")
    x_min, x_max = masses.min(), masses.max()
    
    # obs_color = (103./255., 203./255., 123./255., 0.4)
    obs_color = (135./255., 206./255., 250./255., 0.5)
    
    fig = plt.figure(figsize=(15, 10), dpi=80, facecolor="w", edgecolor="k")
    ax = fig.add_subplot(111)

    ax.xaxis.set_major_formatter(ticker.FormatStrFormatter("%d"))
    ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.1f"))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.05))
    handles = []
    
    # #FIXME: add observed
    #observed = ax.plot(masses, obs_exp, color="k", linestyle="None")
    #handles.append( (mpatches.Patch(color=obs_color), "Observed") )
    
    center = ax.plot(masses, obs_exp, color="k", linestyle="-")
    handles.append( (mlines.Line2D([], [], color="k", linestyle="-"), "Expected") )

    twosig = ax.fill_between(masses, neg2_exp, pos2_exp, color=twosigma)
    handles.append( (mpatches.Patch(color=twosigma), r"95% Expected") )

    onesig = ax.fill_between(masses, pos1_exp, neg1_exp, color=onesigma)
    handles.append( (mpatches.Patch(color=onesigma), r"68% Expected") )
    
    
    ##version bug, the opacity is not handled in mpatches, therefore we make it lighter
    ## alpha = obs_color[3]
    ## patch_color = [i*alpha+1.*(1-alpha) for i in obs_color]
    ## patch_color[3] = 1.
    ##set_trace()
    #upper_contour =  np.array([val if val < 3. else y_leg_cutoff for val in subset['obslower']])
    #observed_contour = ax.fill_between(subset[xvar], subset['obs'], upper_contour, color=obs_color)
    #observed_contour = ax.fill_between(subset[xvar], subset['obsupper'], [y_leg_cutoff for n in xrange(len(subset['obsupper']))], color=obs_color)
    ## observed_contour = ax.fill_between(subset[xvar], subset['obslower'], subset['obsupper'], color=obs_color)
    #lower_to_draw =  np.array([val if val < 3. else np.nan for val in subset['obslower']])
    #observed_lower = ax.plot(
    #	subset[xvar], lower_to_draw, 
    #	color='k', linestyle='None'#, markersize=10, marker='.'
    #	)
    #observed_upper = ax.plot(
    #	subset[xvar], subset['obsupper'], 
    #	color='k', linestyle='None'#, markersize=10, marker='.'
    #	)
    
    #Fake observed, just to check it works
    ## ax.fill_between(
    ## 	xs, [0]*len(xs), [8, 6, 5, 4.3, 3.5, 3, 2, 1.5], 
    ## 	facecolor=obs_color, edgecolor='k', linewidth=1
    ## )

    ax.set_xlabel(addenda["mass"] % (parity, parity, relw)+", "+xlabels["mass"] % parity, fontsize=32)
    ax.set_ylabel("$g_{%s \\rightarrow t\\bar t}$" % parity, fontsize=32)
    ax.autoscale()
    ax.set_ylim(0., min(ax.get_ylim()[1]*1.4, 4.0))
    ax.set_xlim(x_min, x_max)
    #set_trace()
    
    if maxg_values is not None:
    	handles.append((mpatches.Patch(facecolor="none", hatch="||", edgecolor="gray", linewidth=1.), "$\Gamma_{%s \\rightarrow t\\bar t} > \Gamma_{%s}$" % (parity, parity)))
    	maxg_xvalues_todraw = [val[0]for val in maxg_values ]
    	maxg_values_todraw = [val[1] for val in maxg_values]
    	unphys_region = ax.plot(maxg_xvalues_todraw, maxg_values_todraw, color="gray", linestyle="-")
    	ax.fill_between(maxg_xvalues_todraw, maxg_values_todraw, [val+0.05 for val in maxg_values_todraw], color="none", hatch="||", edgecolor="gray", linewidth=0.)
    
    #legend
    legend = ax.legend(
        x(handles), y(handles),
        loc="upper left",
        ncol=2,
        fontsize=29,
        frameon=False
    )
    fontP = matplotlib.font_manager.FontProperties()
    fontP.set_size(29)
    legend.set_title(title="95% $\mathsf{CL_s}$ limits", prop=fontP)


        # add label for signal type
    ax.text(
        0.96, 0.92, type_label,
        fontsize=29, horizontalalignment="right", verticalalignment="bottom", transform=ax.transAxes
    )
    
    #CMS blurb
    hep.cms.label("Preliminary", data=True, lumi=lumi, ax=ax)
    
    ax.tick_params(axis="both", labelsize=29, which="both")

    #set_trace()
    figname = os.path.join(output_dir, f"{args.type}_limit_{parity}_{wname}")
    fig.savefig(figname)
    print(f"{figname} written")
    plt.close(fig)



if __name__ == "__main__":	
    #set_trace()
    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]
    analyzer = "htt_btag_sb_regions"
    user = os.environ["USER"]
    eos_dir = os.environ["eos_dir"]
    plot_dir = os.environ["plots_dir"]

    #eos_base_input_dir = os.path.join("/eos", "user", user[0], user, "NanoAOD_Analyses", "results", jobid, f"Limits_{analyzer}")
    #eos_base_output_dir = os.path.join("/eos", "user", user[0], user, "NanoAOD_Analyses", "plots", jobid, f"Limits_{analyzer}")
    #eos_base_input_dir = os.path.join(eos_dir, "results", jobid, f"Limits_{analyzer}")
    #eos_base_output_dir = os.path.join(plot_dir, jobid, f"Limits_{analyzer}")
    eos_base_input_dir = os.path.join(eos_dir, "results", jobid, f"Limits_{analyzer}")
    eos_base_output_dir = os.path.join(plot_dir, jobid, f"Limits_{analyzer}")
    if args.type == "MEreweighting":
        version = "V21"
    
           # All uncs
        base_lim_dirname = os.path.join(eos_base_input_dir, version, "ASYMv21_lj_nosig_1")
        possible_widths = [1.0, 2.0, 2.5, 3.0, 4.0, 5.0, 8.0, 10.0, 13.0, 15.0, 18.0, 21.0, 25.0]
        #set_trace()
        #   # All uncs
        #base_lim_dirname = os.path.join(eos_base_input_dir, "ASYMv13mewsmoothed_nosig_2")
        #possible_widths = [2.5, 3.0, 4.0, 5.0, 8.0, 10.0, 13.0, 15.0, 18.0, 21.0, 25.0]    
        #   # Some uncs
        #base_lim_dirname = os.path.join(eos_base_input_dir, "ASYMv14mewsmoothed_nosig_4uncpart")
        #possible_widths = [2.5, 10.0, 15.0, 21.0, 25.0]
        #   # No uncs
        #base_lim_dirname = os.path.join(eos_base_input_dir, "ASYMv13mewsmoothed_nosig_3nouncs")
        #possible_widths = [2.5, 10.0, 15.0, 25.0]

        output_dir = os.path.join(eos_base_output_dir, version, os.path.basename(base_lim_dirname))
        #output_dir = os.path.join(eos_base_output_dir, os.path.basename(base_lim_dirname))
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

        possible_widths = [str(width) for width in possible_widths]

        type_label = "$e/\\mu$+jets channel"
        #type_label = "ME reweighted"

    elif args.type == "Morphed":
                #set_trace()
        #   # All uncs
        #base_lim_dirname = os.path.join(eos_base_input_dir, "ASYMv13mcsigmorphedsmoothed_nosig_4")
        #possible_widths = [2.5, 5.0, 10.0, 15.0, 20.0, 25.0]
        #   # Some uncs
        #base_lim_dirname = os.path.join(eos_base_input_dir, "ASYMv13mcsigmorphedsmoothed_nosig_4uncpart")
        #possible_widths = [2.5, 5.0, 10.0, 15.0, 21.0, 25.0]
           # Some uncs
        base_lim_dirname = os.path.join(eos_base_input_dir, "ASYMv13mcsigsmooth_nosig_3nouncs")
        possible_widths = [2.5, 10.0, 15.0, 25.0]

        output_dir = os.path.join(eos_base_output_dir, os.path.basename(base_lim_dirname))
        if not os.path.isdir(output_dir):
            os.makedirs(output_dir)

        possible_widths = [str(width) for width in possible_widths]

        type_label = "Morphed MC"

    
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

    lumi = 138
    onesigma = "#00cc00" # kGreen + 1
    twosigma = "#ffcc00" # kOrange
    #onesigma = '#00f847'
    #twosigma = '#fffc4d'
    line = "#ff1521"

    #set_trace()
    for boson in available_parities:
        for relw in available_widths:
            wname = val2name(relw)
            fname = os.path.join(base_lim_dirname, f"{boson}_relw{wname}.coffea")
            if not os.path.isfile(fname): raise ValueError(f"{fname} not found")
            limits_dict = load(fname)

            masses = np.array(list(limits_dict.keys()), dtype=int)
            maxg_values = [(mass, max_g(boson, mass, relw/100.)) for mass in np.arange(min(masses), max(masses)+5., 5.)]

            make_plot(subset=limits_dict, parity=boson, masses=masses, width=relw, maxg_values=maxg_values)
