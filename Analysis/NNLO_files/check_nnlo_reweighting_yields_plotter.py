#!/usr/bin/env python

# matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend("agg")
from matplotlib import rcParams
rcParams["font.size"] = 20
rcParams["savefig.format"] = "png"
rcParams["savefig.bbox"] = "tight"

import Utilities.Plotter as Plotter
from coffea.hist import plot
from coffea import hist
from coffea.util import load
from pdb import set_trace
import os
import Utilities.prettyjson as prettyjson
import Utilities.plot_tools as plt_tools
import numpy as np
from fnmatch import fnmatch

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
base_jobid = os.environ["base_jobid"]
analyzer = "check_nnlo_reweighting_yields"

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016", "2017", "2018"] if base_jobid == "NanoAODv6" else ["2016APV", "2016", "2017", "2018"], help="Specify which year to run over")
args = parser.parse_args()

f_ext = "TOT.coffea"
input_dir = os.path.join(proj_dir, "results", "%s_%s" % (args.year, base_jobid), analyzer)
fnames = sorted(["%s/%s" % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)])
hdict = plt_tools.add_coffea_files(fnames) if len(fnames) > 1 else load(fnames[0])
    
outdir = os.path.join(proj_dir, "plots", base_jobid, analyzer, args.year)
if not os.path.isdir(outdir):
    os.makedirs(outdir)
    
rewt_style_dict = {
    "Nominal" : ("Nominal", "k"),
    "thad_pt" : ("$W_{Orig}$($p_{T}$($t_{h}$))", "#e42a2c"), ## red
    "thad_pt_Interp" : ("$W_{Int}$($p_{T}$($t_{h}$))", "#4daf4a"), ## green
    "mtt_vs_thad_ctstar" : ("$W_{Orig}$($m_{t\\bar{t}}$, cos($\\theta^{*}_{t_{h}}$))", "#377eb8"), ## blue
    "mtt_vs_thad_ctstar_Interp" : ("$W_{Int}$($m_{t\\bar{t}}$, cos($\\theta^{*}_{t_{h}}$))", "#ff7f00"), ## orange
}

variables = {
        # thad
    "pt_thad" : ("$p_{T}$($t_{h}$) [GeV]", 2, (0., 500.)),
    "eta_thad" : ("$\eta$($t_{h}$)", 2, (-4., 4.)),
    "phi_thad" : ("$\phi$($t_{h}$)", 1, (-3.2, 3.2)),
    "mass_thad" : ("m($t_{h}$) [GeV]", 1, (150., 200.)),
    "ctstar_thad" : ("cos($\\theta^{*}_{t_{h}}$)", 2, (-1., 1.)),
    "ctstar_abs_thad" : ("|cos($\\theta^{*}_{t_{h}}$)|", 1, (0., 1.)),
        # tlep
    "pt_tlep" : ("$p_{T}$($t_{l}$) [GeV]", 2, (0., 500.)),
    "eta_tlep" : ("$\eta$($t_{l}$)", 2, (-4., 4.)),
    "phi_tlep" : ("$\phi$($t_{l}$)", 1, (-3.2, 3.2)),
    "mass_tlep" : ("m($t_{l}$) [GeV]", 1, (150., 200.)),
    "ctstar_tlep" : ("cos($\\theta^{*}_{t_{l}}$)", 2, (-1., 1.)),
    "ctstar_abs_tlep" : ("|cos($\\theta^{*}_{t_{l}}$)|", 1, (0., 1.)),
        # ttbar
    "pt_tt" : ("$p_{T}$($t\\bar{t}$) [GeV]", 2, (0., 500.)),
    "eta_tt" : ("$\\eta$($t\\bar{t}$)", 2, (-4., 4.)),
    "phi_tt" : ("$\phi$($t\\bar{t}$)", 1, (-3.2, 3.2)),
    "mtt" : ("m($t\\bar{t}$) [GeV]", 10, (200., 4000.)),
    "mtt_vs_thad_ctstar" : ("m($t\\bar{t}$) $\otimes$ cos($\\theta^{*}_{t_{h}}$)", 1, None),
}


# make dict to check that yields for all dists are the same
yields_dict = {}

    ## get data lumi and scale MC by lumi
ttSL = "ttJets_PS" if ((args.year == "2016") and (base_jobid == "NanoAODv6")) else "ttJetsSL"
data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, "inputs", "%s_lumis_data.json" % base_jobid)).read())[args.year]
lumi_to_use = (data_lumi_year["Muons"]+data_lumi_year["Electrons"])/2000.

lumi_correction_dict = load(os.path.join(proj_dir, "Corrections", jobid, "MC_LumiWeights.coffea"))[args.year]
avg_lumi_corr = (lumi_correction_dict["Muons"][ttSL]+lumi_correction_dict["Electrons"][ttSL])/2.

process_axis = hist.Cat("process", "Process")
reweighting_axis = hist.Cat("rewt", "Reweighting Type")

for hname in variables.keys():
    if hname not in hdict.keys():
        raise ValueError("%s not found in file" % hname)
    histo = hdict[hname].integrate("dataset")
    histo.scale(avg_lumi_corr)

    is2d = histo.dense_dim() == 2
    axes_to_sum = (histo.dense_axes()[0].name,) if not is2d else (histo.dense_axes()[0].name, histo.dense_axes()[1].name)
    yields_dict[hname] = histo.sum(*axes_to_sum, overflow="all").values()

    #set_trace()
    if is2d:
            # vars only for linearized mtt ctstar hist
        mtt_binning, ctstar_binning = histo.axis(axes_to_sum[0]).edges(), histo.axis(axes_to_sum[1]).edges()
        mtt_binlabels = ["%s $\leq$ m($t\\bar{t}$) $\leq$ %s" % (mtt_binning[bin], mtt_binning[bin+1]) for bin in range(len(mtt_binning)-1)]
        ctstar_binlabels = ["%s $\leq$ cos($\\theta^{*}_{t_{h}}$) $\leq$ %s" % (ctstar_binning[bin], ctstar_binning[bin+1]) for bin in range(len(ctstar_binning)-1)]*len(mtt_binlabels)
        mtt_bin_locs = np.linspace((len(ctstar_binning)-1)/2, (len(ctstar_binning)-1)*(len(mtt_binning)-1)-(len(ctstar_binning)-1)/2, len(mtt_binning)-1)
        vlines = [len(mtt_binlabels)*ybin for ybin in range(1, len(mtt_binlabels))] # last element is hardcoded from binning

            # easier to rename sparse axis than change linearize()
        tmp_histo = histo.copy()
        tmp_histo = tmp_histo.group(histo.sparse_axes()[0].name, process_axis, {key[0]:key[0] for key in histo.values().keys()})
        hline = Plotter.linearize_hist(tmp_histo, no_transpose=True)
            # revert sparse axis name to original
        hline = hline.group(hline.sparse_axes()[0].name, reweighting_axis, {key[0]:key[0] for key in histo.values().keys()})
        histo = hline

    xtitle, rebinning, x_lims = variables[hname]
    if rebinning != 1:
        histo = histo.rebin(*axes_to_sum, rebinning)

        # orig
    fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
    fig.subplots_adjust(hspace=.07)

    for rewt_key in sorted(histo.values().keys()):
            # plot yields    
        ax = plot.plot1d(histo[rewt_key].integrate("rewt"),
            ax=ax, clear=False,
            line_opts={"linestyle" : "-", "color" : rewt_style_dict[rewt_key[0]][1]},
        )
            ## plot ratios
        if rewt_key == ("Nominal",): continue
                # orig
        ratio_masked_vals, ratio_bins = Plotter.get_ratio_arrays(num_vals=histo.values()[rewt_key], denom_vals=histo.values()[("Nominal",)], input_bins=histo.dense_axes()[0].edges())
        rax.step(ratio_bins, ratio_masked_vals, where="post", **{"linestyle" : "-", "color" : rewt_style_dict[rewt_key[0]][1]})

    #set_trace()
    # plot yields
        # format axes
    ax.autoscale(axis="x", tight=True)
    ax.set_xlim(x_lims)
    ax.set_ylim(None)
    ax.set_xlabel(None)
    ax.set_ylabel("Events")
    
    rax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
    rax.set_ylabel("W(var)/Nominal")
    if not is2d: rax.set_xlabel(xtitle)
    
    ## set plotting styles
        ## set legend and corresponding colors
    handles, labels = ax.get_legend_handles_labels()
    for idx, handle in enumerate(handles):
        hstyle = [val for val in rewt_style_dict.values() if handle.get_c() in val[1]]
        labels[idx] = hstyle[0][0]
    # call ax.legend() with the new values
    ax.legend(handles,labels, loc="upper right")

    ax.text(
        0.02, 0.86,
        "$t\\bart \\rightarrow e/\mu + jets$\nparton level",
        fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
    )

    if is2d:
        # plot ctstar labels
        rax.set_xticks(np.arange(len(ctstar_binlabels)))
        rax.set_xticklabels(ctstar_binlabels)
        ax.tick_params(which="minor", bottom=False, top=False)
        rax.tick_params(which="minor", bottom=False, top=False)

        plt.setp(rax.get_xticklabels(), rotation=90, ha="left", fontsize=rcParams["font.size"]*0.5)

        # plot vertical lines for mtt vs ctstar
        for vline in vlines:
                # orig
            ax.axvline(vline, color="k", linestyle="--")
            rax.axvline(vline, color="k", linestyle="--")
    
        # plot mtt labels
        for idx, label in enumerate(mtt_binlabels):
            rax.annotate(label, xy=(mtt_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                xytext=(0, -120), textcoords="offset points", va="top", ha="center", fontsize=rcParams["font.size"]*0.5, rotation=45)

    hep.cms.label(ax=ax, data=False, paper=False, year=args.year, lumi=round(lumi_to_use, 1))    

    figname = os.path.join(outdir, "SemilepTTbar_RewtComp_%s" % hname)
    fig.savefig(figname)
    #set_trace()
    print(f"{figname} written")

# save yields dict    
dec_tolerance = 7
nom_same, pt_same, mtt_ct_same, pt_interp_same, mtt_ct_interp_same = (False,)*5
    # lower decimal tolerance until all three yields have only one value
while not (nom_same and pt_same and mtt_ct_same and pt_interp_same and mtt_ct_interp_same):
    nom_yields = np.around(np.array([yields_dict[key][("Nominal",)] for key in yields_dict.keys()]), decimals=dec_tolerance)
    nom_same = np.all(nom_yields == nom_yields[0]) # values are all equal up to specified decimal place
    if nom_same:
        nom_yield = nom_yields[0]
    
    pt_yields = np.around(np.array([yields_dict[key][("thad_pt",)] for key in yields_dict.keys()]), decimals=dec_tolerance)
    pt_same = np.all(pt_yields == pt_yields[0]) # values are all equal up to specified decimal place
    if pt_same:
        pt_yield = pt_yields[0]
    
    mtt_ct_yields = np.around(np.array([yields_dict[key][("mtt_vs_thad_ctstar",)] for key in yields_dict.keys()]), decimals=dec_tolerance)
    mtt_ct_same = np.all(mtt_ct_yields == mtt_ct_yields[0]) # values are all equal up to specified decimal place
    if mtt_ct_same:
        mtt_ct_yield = mtt_ct_yields[0]

    pt_interp_yields = np.around(np.array([yields_dict[key][("thad_pt_Interp",)] for key in yields_dict.keys()]), decimals=dec_tolerance)
    pt_interp_same = np.all(pt_interp_yields == pt_interp_yields[0]) # values are all equal up to specified decimal place
    if pt_interp_same:
        pt_interp_yield = pt_interp_yields[0]
    
    mtt_ct_interp_yields = np.around(np.array([yields_dict[key][("mtt_vs_thad_ctstar_Interp",)] for key in yields_dict.keys()]), decimals=dec_tolerance)
    mtt_ct_interp_same = np.all(mtt_ct_interp_yields == mtt_ct_interp_yields[0]) # values are all equal up to specified decimal place
    if mtt_ct_interp_same:
        mtt_ct_interp_yield = mtt_ct_interp_yields[0]

    if not (nom_same and pt_same and mtt_ct_same and pt_interp_same and mtt_ct_interp_same): dec_tolerance -= 1

rows = [("Lumi: %s fb^-1" % format(lumi_to_use, ".1f"), "", "Values are the same up to %i decimal places" % dec_tolerance)]
rows += [("Reweighting", "Yield", "Yield/Nominal")]
rows += [("Nominal", format(nom_yield, ".%if" % dec_tolerance), format(nom_yield/nom_yield, ".%if" % dec_tolerance))]
rows += [("thad_pt", format(pt_yield, ".%if" % dec_tolerance), format(pt_yield/nom_yield, ".%if" % dec_tolerance))]
rows += [("mtt_vs_thad_ctstar", format(mtt_ct_yield, ".%if" % dec_tolerance), format(mtt_ct_yield/nom_yield, ".%if" % dec_tolerance))]
rows += [("thad_pt_interp", format(pt_interp_yield, ".%if" % dec_tolerance), format(pt_interp_yield/nom_yield, ".%if" % dec_tolerance))]
rows += [("mtt_vs_thad_ctstar_interp", format(mtt_ct_interp_yield, ".%if" % dec_tolerance), format(mtt_ct_interp_yield/nom_yield, ".%if" % dec_tolerance))]

yields_fname = os.path.join(outdir, "reweighting_yields.txt")
plt_tools.print_table(rows, filename=yields_fname, print_output=True, header_line=1)
print(f"{yields_fname} written")
