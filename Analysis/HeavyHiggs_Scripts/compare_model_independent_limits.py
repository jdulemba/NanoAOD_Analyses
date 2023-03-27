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
rcParams["font.size"] = 20
rcParams["savefig.format"] = "png"
rcParams["savefig.bbox"] = "tight"

import matplotlib
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.ticker as ticker
import os
from coffea.util import load

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("uncs_set", choices=["All", "Some", "No"], help="Choose between limits computed from different sets of uncertainties.")
parser.add_argument("--only_exp", action="store_true", help="Make comparision plots only comparing expected limits, not sigma bands.")
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


def make_plot(limits_dict, parity, masses, width, maxg_values=None):
    #set_trace()
    
    # obs_color = (103./255., 203./255., 123./255., 0.4)
    obs_color = (135./255., 206./255., 250./255., 0.5)
    
    fig = plt.figure(figsize=(15, 10), dpi=80, facecolor="w", edgecolor="k")
    ax = fig.add_subplot(111)
    ax.xaxis.set_major_formatter(
    	ticker.FormatStrFormatter("%d")
    	)
    ax.yaxis.set_major_formatter(
    	ticker.FormatStrFormatter("%.1f")
    	)
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.05))
    handles = []

    x_min, x_max = masses.min(), masses.max()

    #set_trace()
    idx = 0
    for lim_name, limits in limits_dict.items():
        pos2_exp = np.array([val[0] for val in limits.values()], dtype="f4")
        pos1_exp = np.array([val[1] for val in limits.values()], dtype="f4")
        obs_exp = np.array([val[2] for val in limits.values()], dtype="f4")
        neg1_exp = np.array([val[3] for val in limits.values()], dtype="f4")
        neg2_exp = np.array([val[4] for val in limits.values()], dtype="f4")

        if args.only_exp:
            # expected        
            center = ax.plot(np.array(list(limits.keys()), dtype=int), obs_exp, color="k" if idx==0 else "r", linestyle="-")
            handles.append( (mlines.Line2D([], [], color="k" if idx==0 else "r", linestyle="-"), f"Expected {lim_name}") )
        else:
            # expected        
            center = ax.plot(np.array(list(limits.keys()), dtype=int), obs_exp, color="k", linestyle="-" if idx == 0  else "--")
            if idx == 0: handles.append( (mlines.Line2D([], [], color="k", linestyle="-" if idx == 0  else "--"), "Expected") )

                # +/-1 sigma
            neg1_sig = ax.plot(np.array(list(limits.keys()), dtype=int), neg1_exp, color=onesigma, linestyle="-" if idx == 0  else "--")
            if idx == 0: handles.append( (mlines.Line2D([], [], color=onesigma, linestyle="-" if idx == 0  else "--"), "68% Expected") )
            pos1_sig = ax.plot(np.array(list(limits.keys()), dtype=int), pos1_exp, color=onesigma, linestyle="-" if idx == 0  else "--")

                # +/-2 sigma
            neg2_sig = ax.plot(np.array(list(limits.keys()), dtype=int), neg2_exp, color=twosigma, linestyle="-" if idx == 0  else "--")
            if idx == 0: handles.append( (mlines.Line2D([], [], color=twosigma, linestyle="-" if idx == 0  else "--"), "95% Expected") )
            pos2_sig = ax.plot(np.array(list(limits.keys()), dtype=int), pos2_exp, color=twosigma, linestyle="-" if idx == 0  else "--")

            handles.append( (mlines.Line2D([], [], color="k", linestyle="-" if idx == 0  else "--"), lim_name) )

        idx += 1
    
    
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
        numpoints=1,
        fontsize=29,
        frameon=False
    )
    fontP = matplotlib.font_manager.FontProperties()
    fontP.set_size(29)
    #set_trace()
    legend.set_title(title="95% $\mathsf{CL_s}$ limits"+f", {uncs_name}", prop=fontP)
    #legend.set_title(title="95% $\mathsf{CL_s}$ limits", prop=fontP)

    #CMS blurb
    hep.cms.label("Preliminary", data=True, lumi=lumi, ax=ax)
    
    ax.tick_params(axis="both", labelsize=29, which="both")

    #set_trace()
    figname = os.path.join(output_dir, f"compare_limits_{parity}_{wname}_Expected" if args.only_exp else f"compare_limits_{parity}_{wname}")
    fig.savefig(figname)
    print(f"{figname} written")
    plt.close(fig)



if __name__ == "__main__":	
    #set_trace()
    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]
    analyzer = "htt_btag_sb_regions"
    user = os.environ["USER"]
    
    eos_base_input_dir = os.path.join("/eos", "user", user[0], user, "NanoAOD_Analyses", "results", jobid, f"Limits_{analyzer}")
    eos_base_output_dir = os.path.join("/eos", "user", user[0], user, "NanoAOD_Analyses", "plots", jobid, f"Limits_{analyzer}", "Comparisons")

        # All uncs
    if args.uncs_set == "All":
            # get ME reweight results
        MErewt_base_lim_dirname = os.path.join(eos_base_input_dir, "ASYMv13mewsmoothed_nosig_2")
        MErewt_possible_widths = [2.5, 3.0, 4.0, 5.0, 8.0, 10.0, 13.0, 15.0, 18.0, 21.0, 25.0]
            # get morphed mc results
        Morph_base_lim_dirname = os.path.join(eos_base_input_dir, "ASYMv13mcsigmorphedsmoothed_nosig_4")
        Morph_possible_widths = [2.5, 5.0, 10.0, 15.0, 20.0, 25.0]

        uncs_name = "All Uncertainties"

        # Some uncs
    if args.uncs_set == "Some":
           # get ME reweight results
        MErewt_base_lim_dirname = os.path.join(eos_base_input_dir, "ASYMv14mewsmoothed_nosig_4uncpart")
        MErewt_possible_widths = [2.5, 10.0, 15.0, 21.0, 25.0]
           # get morphed mc results
        Morph_base_lim_dirname = os.path.join(eos_base_input_dir, "ASYMv13mcsigmorphedsmoothed_nosig_4uncpart")
        Morph_possible_widths = [2.5, 5.0, 10.0, 15.0, 21.0, 25.0]

        uncs_name = "Some Uncertainties"

        # No uncs
    if args.uncs_set == "No":
            # get ME reweight results
        MErewt_base_lim_dirname = os.path.join(eos_base_input_dir, "ASYMv13mewsmoothed_nosig_3nouncs")
        MErewt_possible_widths = [2.5, 10.0, 15.0, 25.0]
            # get morphed mc results
        Morph_base_lim_dirname = os.path.join(eos_base_input_dir, "ASYMv13mcsigsmooth_nosig_3nouncs")
        Morph_possible_widths = [2.5, 10.0, 15.0, 25.0]

        uncs_name = "No Uncertainties"

    output_dir = os.path.join(eos_base_output_dir, f"Comp_{os.path.basename(MErewt_base_lim_dirname)}_{os.path.basename(Morph_base_lim_dirname)}")
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)


    possible_widths = sorted(set(MErewt_possible_widths) - (set(MErewt_possible_widths) - set(Morph_possible_widths)))
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

    lumi = 137    
    onesigma = "#00cc00" # kGreen + 1
    twosigma = "#ffcc00" # kOrange
    line = "#ff1521"

    #set_trace()
    for boson in available_parities:
        for relw in available_widths:
            wname = val2name(relw)
                # get ME reweight files
            ME_fname = os.path.join(MErewt_base_lim_dirname, f"{boson}_relw{wname}.coffea")
            if not os.path.isfile(ME_fname): raise ValueError(f"{ME_fname} not found")
            ME_limits_dict = load(ME_fname)

                # get morphed MC files
            Morph_fname = os.path.join(Morph_base_lim_dirname, f"{boson}_relw{wname}.coffea")
            if not os.path.isfile(Morph_fname): raise ValueError(f"{Morph_fname} not found")
            Morph_limits_dict = load(Morph_fname)

                # get available mass points between the two methods
            ME_masses = np.array(list(ME_limits_dict.keys()), dtype=int)
            Morph_masses = np.array(list(Morph_limits_dict.keys()), dtype=int)

            masses = np.unique(np.concatenate((np.arange(min(ME_masses), max(ME_masses)+5., 5.), np.arange(min(Morph_masses), max(Morph_masses)+5., 5.))))
            maxg_values = [(mass, max_g(boson, mass, relw/100.)) for mass in np.arange(min(masses), max(masses)+5., 5.)]

            make_plot(limits_dict={"ME reweighted" : ME_limits_dict, "Morphed MC" : Morph_limits_dict}, parity=boson, masses=masses, width=relw, maxg_values=maxg_values)

    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
