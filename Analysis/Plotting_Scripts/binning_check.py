# matplotlib
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
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
import Utilities.plot_tools as plt_tools
import Utilities.prettyjson as prettyjson
import re
from coffea import hist
import numpy as np
import fnmatch
import Utilities.Plotter as Plotter
import uproot3
from Utilities.HistExport import export2d

base_jobid = os.environ["base_jobid"]

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("lepton", choices=["Electron", "Muon"], help="Choose which lepton to make plots for")
parser.add_argument("plots", choices=["RECO", "GEN", "RESO", "REL_RESO", "All"], help="Choose which hists to make plots for")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
analyzer = "binning_check"

input_dir = os.path.join(proj_dir, "results", f"{args.year}_{jobid}", analyzer)
f_ext = "TOT.coffea"
#outdir = os.path.join(proj_dir, "plots", f"{args.year}_{jobid}", analyzer, "NoNNLOqcdNLOew")
outdir = os.path.join(proj_dir, "plots", f"{args.year}_{jobid}", analyzer)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

fnames = sorted(["%s/%s" % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)])

#set_trace()
hdict = plt_tools.add_coffea_files(fnames) if len(fnames) > 1 else load(fnames[0])

jet_mults = {
    "3Jets" : "3 jets",
    "4PJets" : "4+ jets"
}

objtypes = {
    "Jets" : "jets",
    "Lep" :  {
        "Muon" : "$\\mu$",
        "Electron" : "$e$",
    }
}

rewt_style_dict = {
    "TT" : {"label" : "SM $t\\bar{t}$", "facecolor" : "none", "hatch" : "//", "edgecolor" : "r"},
    "Res" : {"label" : "$A_{400}^{10\\%}$, res.", "facecolor" : "#377eb8", "edgecolor" : "#377eb8"}, ## blue
    "Int_neg" : {"label" : "$A_{400}^{10\\%}$, int. w $<$ 0 ", "color" : "k", "linestyle" : "-"},
    "Int_pos" : {"label" : "$A_{400}^{10\\%}$, int. w $>$ 0 ", "color" : "k", "linestyle" : "--"},
}


reco_variables = {
    #"RECO_mtt_vs_tlep_ctstar" : ("$m_{t\\bar{t}}$ [GeV]", "cos($\\theta^{*}_{t_{l}}$)", 1, 1, (200., 2000.),  (-1., 1.)),
    #"RECO_mtt_vs_tlep_ctstar_abs" : ("$m_{t\\bar{t}}$ [GeV]", "|cos($\\theta^{*}_{t_{l}}$)|", 1, 1, (200., 2000.),  (0., 1.)),
    "RECO_mtt_vs_topRapidity" : ("$m_{t\\bar{t}}$ [GeV]", "$y_t$", 1, 1, (200., 2000.), (-5., 5.)),
    "RECO_mtt_vs_deltaYtt" : ("$m_{t\\bar{t}}$ [GeV]", "$\\delta y_{t\\bar{t}}$", 1, 1, (200., 2000.), (-10., 10.)),
}
gen_variables = {
    #"GEN_mtt_vs_tlep_ctstar" : ("$m_{t\\bar{t}}$ [GeV]", "|cos($\\theta^{*}_{t_{l}}$)|", 1, 1, (200., 2000.),  (-1., 1.)),
    #"GEN_mtt_vs_tlep_ctstar_abs" : ("$m_{t\\bar{t}}$ [GeV]", "|cos($\\theta^{*}_{t_{l}}$)|", 1, 1, (200., 2000.),  (0., 1.)),
    "GEN_mtt_vs_topRapidity" : ("$m_{t\\bar{t}}$ [GeV]", "$y_t$", 1, 1, (200., 2000.), (-5., 5.)),
    "GEN_mtt_vs_deltaYtt" : ("$m_{t\\bar{t}}$ [GeV]", "$\\delta y_{t\\bar{t}}$", 1, 1, (200., 2000.), (-10., 10.)),
}
reso_variables = {
    "RESO_mtt_vs_mtt" : ("$m_{t\\bar{t}}$ [GeV]", "Gen-Reco $m_{t\\bar{t}}$ [GeV]", 1, 1, (200., 2000.),  (-500., 500.)),
    "RESO_ctstar_vs_mtt" : ("$m_{t\\bar{t}}$ [GeV]", "Gen-Reco cos($\\theta^{*}_{t_{l}}$)", 1, 1, (200., 2000.),  (-2., 2.)),
    "RESO_ctstar_abs_vs_mtt" : ("$m_{t\\bar{t}}$ [GeV]", "Gen-Reco |cos($\\theta^{*}_{t_{l}}$)|", 1, 1, (200., 2000.),  (-1., 1.)),
    #"RESO_mtt_vs_tlep_ctstar" : ("$m_{t\\bar{t}}$ [GeV]", "cos($\\theta^{*}_{t_{l}}$)", 1, 1, (-500., 500.),  (-2., 2.)),
    #"RESO_mtt_vs_tlep_ctstar_abs" : ("$m_{t\\bar{t}}$ [GeV]", "|cos($\\theta^{*}_{t_{l}}$)|", 1, 1, (-500., 500.),  (-1., 1.)),
}
rel_reso_variables = {
    "REL_RESO_mtt_vs_mtt" : ("$m^{reco}_{t\\bar{t}}$ [GeV]", "($m^{reco}_{t\\bar{t}}$-$m^{parton}_{t\\bar{t}}$)/$m^{parton}_{t\\bar{t}}$", 1, 1, (200., 2000.),  (-2., 2.)),
}


    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())[args.year]
lumi_correction = load(os.path.join(proj_dir, "Corrections", base_jobid, "MC_LumiWeights.coffea"))[args.year][f"{args.lepton}s"]
        # scale ttJets events, split by reconstruction type, by normal ttJets lumi correction
ttJets_permcats = ["*right", "*matchable", "*unmatchable", "*sl_tau", "*other"]
#set_trace()
names = [dataset for dataset in sorted(set([key[0] for key in hdict["RECO_mtt_vs_tlep_ctstar_abs"].values().keys()]))] # get dataset names in hists
ttJets_cats = [name for name in names if any([fnmatch.fnmatch(name, cat) for cat in ttJets_permcats])] # gets ttJets(_PS)_other, ...
if len(ttJets_cats) > 0:
    for tt_cat in ttJets_cats:
        ttJets_lumi_topo = "_".join(tt_cat.split("_")[:-2]) if "sl_tau" in tt_cat else "_".join(tt_cat.split("_")[:-1]) # gets ttJets[SL, Had, DiLep] or ttJets_PS
        ttJets_eff_lumi = lumi_correction[ttJets_lumi_topo]
        lumi_correction.update({tt_cat: ttJets_eff_lumi})

## make groups based on process
process = hist.Cat("process", "Process", sorting="placement")
process_cat = "dataset"
#process_groups = plt_tools.make_dataset_groups(args.lepton, args.year, samples=names, gdict="templates")
process_groups = plt_tools.make_dataset_groups(args.lepton, args.year, samples=names)
#set_trace()

    # scale and group hists by process
for hname in hdict.keys():
    if "cutflow" in hname: continue
    #set_trace()
    hdict[hname].scale(lumi_correction, axis="dataset") # scale hists
    hdict[hname] = hdict[hname].group(process_cat, process, process_groups) # group by process
    hdict[hname] = hdict[hname][:, :, args.lepton].integrate("leptype") # only pick out specified lepton



#set_trace()
    ## make reco plots
if (args.plots == "RECO") or (args.plots == "All"):
    for hname in reco_variables.keys():
        if hname not in hdict.keys():
            raise ValueError(f"{hname} not found in file")
        #set_trace()
        histo = hdict[hname].copy()
        #set_trace()
    
        if histo.dense_dim() == 2:
            xtitle, ytitle, xrebinning, yrebinning, x_lims, y_lims = reco_variables[hname]
    
            xaxis_name = histo.dense_axes()[0].name
            yaxis_name = histo.dense_axes()[1].name
                ## rebin x axis
            if isinstance(xrebinning, np.ndarray):
                new_xbins = hist.Bin(xaxis_name, xaxis_name, xrebinning)
            elif isinstance(xrebinning, float) or isinstance(xrebinning, int):
                new_xbins = xrebinning
            histo = histo.rebin(xaxis_name, new_xbins)
                ## rebin y axis
            if isinstance(yrebinning, np.ndarray):
                new_ybins = hist.Bin(yaxis_name, yaxis_name, yrebinning)
            elif isinstance(yrebinning, float) or isinstance(yrebinning, int):
                new_ybins = yrebinning
            histo = histo.rebin(yaxis_name, new_ybins)
    
            ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
            for jmult in sorted(set([key[1] for key in histo.values().keys()])):
                print(*[jmult, hname], sep=", ")
                pltdir = os.path.join(outdir, args.lepton, jmult)
                if not os.path.isdir(pltdir):
                    os.makedirs(pltdir)
    
                hslice = histo[:, jmult].integrate("jmult")
    
                    # make 1D projection along dense axes
                for dax in range(2):
                    if dax == 0:
                        xlabel = ytitle
                        xlimits = y_lims
                        figname = os.path.join(pltdir, "_".join([jmult, args.lepton, "RECO", yaxis_name, "Norm"]))
                        figname_ratio = os.path.join(pltdir, "_".join([jmult, args.lepton, "RECO", yaxis_name, "REStoTTbar"]))
                    else:
                        xlabel = xtitle
                        xlimits = x_lims
                        figname = os.path.join(pltdir, "_".join([jmult, args.lepton, "RECO", xaxis_name, "Norm"]))
                        figname_ratio = os.path.join(pltdir, "_".join([jmult, args.lepton, "RECO", xaxis_name, "REStoTTbar"]))
    
                    hproj = hslice.integrate(hslice.dense_axes()[dax].name)
                    bin_edges = hproj.dense_axes()[0].edges()
    
                        # normalized plots
                    fig, ax = plt.subplots()
                    fig.subplots_adjust(hspace=.07)
    
                    for proc in [("AtoTTJetsSL_M400_W10p0_Res",), ("AtoTTJetsSL_M400_W10p0_Int_neg",), ("AtoTTJetsSL_M400_W10p0_Int_pos",), ("TT",)]:
                        norm_vals = np.r_[hproj.values()[proc], hproj.values()[proc][-1]]/np.sum(hproj.values(overflow="all")[proc])
                        if proc[0] == "TT":
                            tt_vals = norm_vals
                            ax.fill_between(bin_edges, norm_vals, step="post", **rewt_style_dict["TT"])
                        elif "Res" in proc[0]:
                            res_vals = norm_vals
                            ax.fill_between(bin_edges, norm_vals, step="post", **rewt_style_dict["Res"])
                        elif "Int_neg" in proc[0]:
                            ax.step(bin_edges, norm_vals, where="post", **rewt_style_dict["Int_neg"])
                        elif "Int_pos" in proc[0]:
                            ax.step(bin_edges, norm_vals, where="post", **rewt_style_dict["Int_pos"])
                            
                    ax.legend(loc="upper right")
                    ax.autoscale()
                    ax.set_ylim(0, ax.get_ylim()[1]*1.15)
                    ax.set_xlim(xlimits)
                    ax.set_xlabel(xlabel)
                    ax.set_ylabel("Probability Density [a.u.]")
    
                        # add lepton/jet multiplicity label
                    ax.text(
                        0.02, 0.88, "%s, %s\ndetector level" % (objtypes["Lep"][args.lepton], jet_mults[jmult]),
                        fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                    )
                    hep.cms.label(ax=ax, data=False, year=args.year)
    
                    #set_trace()
                    fig.savefig(figname)
                    print(f"{figname} written")
                    plt.close()
    
                    #set_trace()
                        # plot resonant/ttbar
                    fig_ratio, ax_ratio = plt.subplots()
                    fig_ratio.subplots_adjust(hspace=.07)
    
                    ax_ratio.step(bin_edges, res_vals/tt_vals, where="post", **{"color" : "k"})
                            
                    #ax_ratio.legend(loc="upper right")
                    ax_ratio.autoscale()
                    ax_ratio.set_ylim(0, ax_ratio.get_ylim()[1]*1.15)
                    ax_ratio.set_xlim(xlimits)
                    ax_ratio.set_xlabel(xlabel)
                    ax_ratio.set_ylabel("Resonant/SM $t\\bar{t}$ [a.u.]")
                    ax_ratio.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
    
                        # add lepton/jet multiplicity label
                    ax_ratio.text(
                        0.02, 0.88, "%s, %s\ndetector level" % (objtypes["Lep"][args.lepton], jet_mults[jmult]),
                        fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax_ratio.transAxes
                    )
                    hep.cms.label(ax=ax_ratio, data=False, year=args.year)
    
                    #set_trace()
                    fig_ratio.savefig(figname_ratio)
                    print(f"{figname_ratio} written")
                    plt.close()


    ## make gen plots
if (args.plots == "GEN") or (args.plots == "All"):
    for hname in gen_variables.keys():
        if hname not in hdict.keys():
            raise ValueError(f"{hname} not found in file")
        #set_trace()
        histo = hdict[hname].copy()
        #set_trace()
    
        if histo.dense_dim() == 2:
            xtitle, ytitle, xrebinning, yrebinning, x_lims, y_lims = gen_variables[hname]
    
            xaxis_name = histo.dense_axes()[0].name
            yaxis_name = histo.dense_axes()[1].name
                ## rebin x axis
            if isinstance(xrebinning, np.ndarray):
                new_xbins = hist.Bin(xaxis_name, xaxis_name, xrebinning)
            elif isinstance(xrebinning, float) or isinstance(xrebinning, int):
                new_xbins = xrebinning
            histo = histo.rebin(xaxis_name, new_xbins)
                ## rebin y axis
            if isinstance(yrebinning, np.ndarray):
                new_ybins = hist.Bin(yaxis_name, yaxis_name, yrebinning)
            elif isinstance(yrebinning, float) or isinstance(yrebinning, int):
                new_ybins = yrebinning
            histo = histo.rebin(yaxis_name, new_ybins)
    
            ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
            for jmult in sorted(set([key[1] for key in histo.values().keys()])):
                print(*[jmult, hname], sep=", ")
                pltdir = os.path.join(outdir, args.lepton, jmult)
                if not os.path.isdir(pltdir):
                    os.makedirs(pltdir)
    
                hslice = histo[:, jmult].integrate("jmult").integrate("process")
   
                    # make 1D projection along dense axes
                for dax in range(2):
                    if dax == 0:
                        xlabel = ytitle
                        xlimits = y_lims
                        figname = os.path.join(pltdir, "_".join([jmult, args.lepton, "GEN", yaxis_name]))
                    else:
                        xlabel = xtitle
                        xlimits = x_lims
                        figname = os.path.join(pltdir, "_".join([jmult, args.lepton, "GEN", xaxis_name]))
    
                    hproj = hslice.integrate(hslice.dense_axes()[dax].name)
                    bin_edges = hproj.dense_axes()[0].edges()
                    vals = hproj.values()[()]
    
                        # plots
                    fig, ax = plt.subplots()
                    fig.subplots_adjust(hspace=.07)
   
                    Plotter.plot_1D(vals, bin_edges, xlimits=xlimits, xlabel=xlabel, ax=ax, label=rewt_style_dict["TT"]["label"])
                            
                    ax.legend(loc="upper right")
                    ax.set_ylim(0, ax.get_ylim()[1]*1.15)
    
                        # add lepton/jet multiplicity label
                    ax.text(
                        0.02, 0.88, "%s, %s\nparton level" % (objtypes["Lep"][args.lepton], jet_mults[jmult]),
                        fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                    )
                    hep.cms.label(ax=ax, data=False, year=args.year, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))
    
                    #set_trace()
                    fig.savefig(figname)
                    print(f"{figname} written")
                    plt.close()


    ## make resolution plots
if (args.plots == "RESO") or (args.plots == "All"):
    for hname in reso_variables.keys():
        if hname not in hdict.keys():
            raise ValueError(f"{hname} not found in file")
        #set_trace()
        histo = hdict[hname].copy()
        #set_trace()
    
        if histo.dense_dim() == 2:
            xtitle, ytitle, xrebinning, yrebinning, x_lims, y_lims = reso_variables[hname]
    
            xaxis_name = histo.dense_axes()[0].name
            yaxis_name = histo.dense_axes()[1].name
                ## rebin x axis
            if isinstance(xrebinning, np.ndarray):
                new_xbins = hist.Bin(xaxis_name, xaxis_name, xrebinning)
            elif isinstance(xrebinning, float) or isinstance(xrebinning, int):
                new_xbins = xrebinning
            histo = histo.rebin(xaxis_name, new_xbins)
                ## rebin y axis
            if isinstance(yrebinning, np.ndarray):
                new_ybins = hist.Bin(yaxis_name, yaxis_name, yrebinning)
            elif isinstance(yrebinning, float) or isinstance(yrebinning, int):
                new_ybins = yrebinning
            histo = histo.rebin(yaxis_name, new_ybins)
    
            ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
            for jmult in sorted(set([key[1] for key in histo.values().keys()])):
                print(*[jmult, hname], sep=", ")
                pltdir = os.path.join(outdir, args.lepton, jmult)
                if not os.path.isdir(pltdir):
                    os.makedirs(pltdir)
    
                hslice = histo[:, jmult].integrate("jmult").integrate("process")
    
                    # make 1D projection along dense axes
                for dax in range(2):
                    if dax == 0:
                        xlabel = ytitle
                        #xlabel = f"Gen-Reco {ytitle}"
                        xlimits = y_lims
                        figname = os.path.join(pltdir, "_".join([jmult, args.lepton, "RESO", yaxis_name]))
                    else:
                        xlabel = xtitle
                        #xlabel = f"Gen-Reco {xtitle}"
                        xlimits = x_lims
                        figname = os.path.join(pltdir, "_".join([jmult, args.lepton, "RESO", xaxis_name]))
   
                    #set_trace() 
                    hproj = hslice.integrate(hslice.dense_axes()[dax].name)
                    bin_edges = hproj.dense_axes()[0].edges()
                    bin_centers = hproj.dense_axes()[0].centers()
                    vals = hproj.values()[()]

                    mean = np.average(bin_centers, weights=vals)
                    stdev = np.sqrt(np.average((bin_centers-mean)**2, weights=vals))                    
    
                        # plots
                    fig, ax = plt.subplots()
                    fig.subplots_adjust(hspace=.07)
  
                    #set_trace() 
                    Plotter.plot_1D(vals, bin_edges, xlimits=xlimits, xlabel=xlabel, ax=ax, label=rewt_style_dict["TT"]["label"])
                            
                    ax.legend(loc="upper right")
                    ax.set_ylim(0, ax.get_ylim()[1]*1.15)
    
                        # add lepton/jet multiplicity label
                    ax.text(
                        0.02, 0.86, "%s, %s\nmean=%.4f\nstd=%.4f" % (objtypes["Lep"][args.lepton], jet_mults[jmult], mean, stdev),
                        fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                    )
                    hep.cms.label(ax=ax, data=False, year=args.year, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))
    
                    #set_trace()
                    fig.savefig(figname)
                    print(f"{figname} written")
                    plt.close()


    ## make relative resolution plots
if (args.plots == "REL_RESO") or (args.plots == "All"):
    tmp_rname = "output_NNLO_NLOwts_relreso_splitprocs.root"
    fout = uproot3.recreate(tmp_rname, compression=uproot3.ZLIB(4)) if os.path.isfile(tmp_rname) else uproot3.create(tmp_rname)
    for hname in rel_reso_variables.keys():
        if hname not in hdict.keys():
            raise ValueError(f"{hname} not found in file")
        #set_trace()
        histo = hdict[hname].copy()
        #set_trace()
    
        if histo.dense_dim() == 2:
            xtitle, ytitle, xrebinning, yrebinning, x_lims, y_lims = rel_reso_variables[hname]
    
            xaxis_name = histo.dense_axes()[0].name
            yaxis_name = histo.dense_axes()[1].name
                ## rebin x axis
            if isinstance(xrebinning, np.ndarray):
                new_xbins = hist.Bin(xaxis_name, xaxis_name, xrebinning)
            elif isinstance(xrebinning, float) or isinstance(xrebinning, int):
                new_xbins = xrebinning
            histo = histo.rebin(xaxis_name, new_xbins)
                ## rebin y axis
            if isinstance(yrebinning, np.ndarray):
                new_ybins = hist.Bin(yaxis_name, yaxis_name, yrebinning)
            elif isinstance(yrebinning, float) or isinstance(yrebinning, int):
                new_ybins = yrebinning
            histo = histo.rebin(yaxis_name, new_ybins)
    
            ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
            for jmult in sorted(set([key[1] for key in histo.values().keys()])):
                print(*[jmult, hname], sep=", ")
                pltdir = os.path.join(outdir, args.lepton, jmult)
                if not os.path.isdir(pltdir):
                    os.makedirs(pltdir)
    
                #hslice = histo[:, jmult].integrate("jmult").integrate("process")
                hslice = histo[:, jmult].integrate("jmult")
   
                    # make 1D projection along dense axes
                for dax in range(1):
                    if dax == 0:
                        xlabel = ytitle
                        xlimits = y_lims
                        figname = os.path.join(pltdir, "_".join([jmult, args.lepton, "REL_RESO", yaxis_name]))
                    else:
                        xlabel = xtitle
                        xlimits = x_lims
                        figname = os.path.join(pltdir, "_".join([jmult, args.lepton, "REL_RESO", xaxis_name]))
   
                    #set_trace() 
                    hproj = hslice.integrate(hslice.dense_axes()[dax].name)
                    bin_edges = hproj.dense_axes()[0].edges()
                    bin_centers = hproj.dense_axes()[0].centers()
                    #vals = hproj.values()[()]
                    for proc in hproj.values().keys():
                        proc_histo = hproj[proc].integrate("process")
                        vals = proc_histo.values()[()]

                        mean = np.average(bin_centers, weights=vals)
                        stdev = np.sqrt(np.average((bin_centers-mean)**2, weights=vals))                    
    
                            # plots
                        fig, ax = plt.subplots()
                        fig.subplots_adjust(hspace=.07)
  
                        Plotter.plot_1D(vals, bin_edges, xlimits=xlimits, xlabel=xlabel, ax=ax, label=proc[0])
                        #Plotter.plot_1D(vals, bin_edges, xlimits=xlimits, xlabel=xlabel, ax=ax, label=rewt_style_dict["TT"]["label"])
                                
                        ax.legend(loc="upper right")
                        ax.set_ylim(0, ax.get_ylim()[1]*1.15)
    
                            # add lepton/jet multiplicity label
                        ax.text(
                            0.02, 0.86, "%s, %s\nmean=%.4f\nstdev=%.4f" % (objtypes["Lep"][args.lepton], jet_mults[jmult], mean, stdev),
                            fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                        )
                        hep.cms.label(ax=ax, data=False, year=args.year, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))

                        figname = "_".join([figname, proc[0]]) 
                        #set_trace()
                        fig.savefig(figname)
                        print(f"{figname} written")
                        plt.close()

                        fout[f"{hname}_{jmult}_{proc[0]}"] = export2d(hslice[proc].integrate("process"))
                #set_trace()
                #fout[f"{hname}_{jmult}"] = export2d(hslice)

                # make plots of relative bin width and resolution for each bin in mtt
                #test_mtt_bins = np.array([300., 360.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 625.0, 650.0, 675.0, 
                #    700.0, 730.0, 760.0, 800.0, 850.0, 900.0, 950., 1000.0, 1050., 1100.0, 1150., 1200., 1300., 1500., 2000.])
                #hslice = hslice.rebin(xaxis_name, hist.Bin(xaxis_name, xaxis_name, test_mtt_bins))

                    # calculate bin widths
                mtt_edges = hslice.axis("mtt").edges()
                #mtt_edges = mtt_edges[np.argwhere(mtt_edges == 300.)[0][0]:] # start mtt bins at 300 GeV
                mtt_min, mtt_max = mtt_edges[0], mtt_edges[-1]
                rel_mtt_bin_widths = np.array([(mtt_edges[idx+1]-mtt_edges[idx])/(mtt_max-mtt_min) for idx in range(mtt_edges.size-1)])

                rel_res_bin_centers = hslice.axis("rel_reso_mtt").centers()
                #rel_res_list = []
                for proc in hproj.values().keys():
                    proc_histo = hslice[proc].integrate("process")
                    rel_res_array = np.zeros(mtt_edges.size-1)
                        # loop over mtt bin values and calculate mean, stdev of relative resolution
                    #set_trace()
                    for idx in range(mtt_edges.size-1):
                        mtt_low, mtt_hi = mtt_edges[idx], mtt_edges[idx+1]
                        rel_res_vals = proc_histo[mtt_low:mtt_hi, :].integrate("mtt").values()[()]
                        #rel_res_vals = hslice[mtt_low:mtt_hi, :].integrate("mtt").values()[()]

                        
                        if np.sum(rel_res_vals) > 0:
                            bin_mean = np.average(rel_res_bin_centers, weights=rel_res_vals)
                            bin_stdev = np.sqrt(np.average((rel_res_bin_centers-bin_mean)**2, weights=rel_res_vals))
                            rel_res_array[idx] = bin_stdev
                        else:
                            rel_res_array[idx] = np.nan
    
                    #set_trace()
                    fig, ax = plt.subplots()
                    fig.subplots_adjust(hspace=.07)
  
                    Plotter.plot_1D(rel_res_array, mtt_edges, xlimits=x_lims, xlabel=xtitle, ylabel="a.u.", ax=ax, label="Relative resolution", color="r")
                    Plotter.plot_1D(rel_mtt_bin_widths, mtt_edges, xlimits=x_lims, xlabel=xtitle, ylabel="a.u.", ax=ax, label="Relative bin width", color="k")
                            
                    ax.legend(loc="upper right", title=proc[0])
                    #ax.legend(loc="upper right")
                    ax.set_ylim(0, ax.get_ylim()[1]*1.15)
    
                        # add lepton/jet multiplicity label
                    ax.text(
                        0.02, 0.92, "%s, %s" % (objtypes["Lep"][args.lepton], jet_mults[jmult]),
                        fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                    )
                    hep.cms.label(ax=ax, data=False, year=args.year, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))

                    #figname = os.path.join(pltdir, "_".join([jmult, args.lepton, "RelativeResolutionBinWidth", "vs", xaxis_name]))
                    figname = os.path.join(pltdir, "_".join([jmult, args.lepton, "RelativeResolutionBinWidth", "vs", xaxis_name, proc[0]]))
                    #set_trace()
                    fig.savefig(figname)
                    print(f"{figname} written")
                    plt.close()
    fout.close()
    print(f"{tmp_rname} written")
