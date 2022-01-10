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
from sklearn.decomposition import PCA
import itertools

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("lepton", choices=["Electron", "Muon"], help="Choose which lepton to make plots for")
parser.add_argument("--plot", action="store_true", help="Make plots")
parser.add_argument("--pca", action="store_true", help="Perform PCA")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
base_jobid = os.environ["base_jobid"]
jobid = os.environ["jobid"]
analyzer = "htt_pdfUncs"

input_dir = os.path.join(proj_dir, "results", f"{args.year}_{jobid}", analyzer)
f_ext = "TOT.coffea"
outdir = os.path.join(proj_dir, "plots", f"{args.year}_{jobid}", analyzer)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

fnames = sorted(["%s/%s" % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)])
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

linearize_binning = (
    #np.array([300.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 625.0, 650.0, 675.0, 700.0, 730.0, 760.0, 800.0, 850.0, 900.0, 1000.0, 1200.0]), # orig
    #np.array([300.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 625.0, 650.0, 675.0,
    #            700.0, 730.0, 760.0, 800.0, 850.0, 900.0, 950., 1000.0, 1050., 1100.0, 1150., 1200., 1300., 1500., 2000.]),
    ##np.array([300.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 620.0, 650., 700.0, 750.0, 800.0, 850.0, 900.0, 2000.0]),
    np.array([360.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 625.0, 650.0, 675.0,
        700.0, 730.0, 760.0, 800.0, 850.0, 900.0, 950., 1000.0, 1050., 1100.0, 1150., 1200., 1300., 1500., 2000.]),
    np.array([0.0, 0.4, 0.6, 0.75, 0.9, 1.0])
)

ctstar_binlabels = ["%s $\leq$ |cos($\\theta^{*}_{t_{l}}$)| $\leq$ %s" % (linearize_binning[1][bin], linearize_binning[1][bin+1]) for bin in range(len(linearize_binning[1])-1)]
ctstar_bin_locs = np.linspace((len(linearize_binning[0])-1)/2, (len(linearize_binning[0])-1)*(len(linearize_binning[1])-1) - (len(linearize_binning[0])-1)/2, len(linearize_binning[1])-1)
mtt_binlabels = ["%s $\leq$ m($t\\bar{t}$) $\leq$ %s" % (linearize_binning[0][bin], linearize_binning[0][bin+1]) for bin in range(len(linearize_binning[0])-1)]*len(ctstar_binlabels)


variables = {
    "mtt_vs_tlep_ctstar_abs" : ("$m_{t\\bar{t}}$", "|cos($\\theta^{*}_{t_{l}}$)|", linearize_binning[0], linearize_binning[1], (linearize_binning[0][0], linearize_binning[0][-1]),
        (linearize_binning[1][0], linearize_binning[1][-1]), True, None, (ctstar_binlabels, ctstar_bin_locs)),
}


    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())[args.year]
lumi_correction = load(os.path.join(proj_dir, "Corrections", base_jobid, "MC_LumiWeights.coffea"))[args.year][f"{args.lepton}s"]
        # scale ttJets events, split by reconstruction type, by normal ttJets lumi correction
ttJets_permcats = ["*right", "*matchable", "*unmatchable", "*sl_tau", "*other"]
names = [dataset for dataset in sorted(set([key[0] for key in hdict["mtt_vs_tlep_ctstar_abs"].values().keys()]))] # get dataset names in hists
ttJets_cats = [name for name in names if any([fnmatch.fnmatch(name, cat) for cat in ttJets_permcats])] # gets ttJets(_PS)_other, ...
if len(ttJets_cats) > 0:
    for tt_cat in ttJets_cats:
        ttJets_lumi_topo = "_".join(tt_cat.split("_")[:-2]) if "sl_tau" in tt_cat else "_".join(tt_cat.split("_")[:-1]) # gets ttJets[SL, Had, DiLep] or ttJets_PS
        ttJets_eff_lumi = lumi_correction[ttJets_lumi_topo]
        lumi_correction.update({tt_cat: ttJets_eff_lumi})

## make groups based on process
process = hist.Cat("process", "Process", sorting="placement")
process_cat = "dataset"
process_groups = {"TT" : names}

    ## make plots
if args.plot:
        # scale and group hists by process
    for hname in hdict.keys():
        if "cutflow" in hname: continue
        #set_trace()
        hdict[hname].scale(lumi_correction, axis="dataset") # scale hists
        hdict[hname] = hdict[hname].group(process_cat, process, process_groups) # group by process
        #hdict[hname] = hdict[hname][:, :, :, args.lepton].integrate("leptype") # only pick out specified lepton
        hdict[hname] = hdict[hname][:, :, :, args.lepton].integrate("leptype").integrate("process") # only pick out specified lepton

    for hname in variables.keys():
        if hname not in hdict.keys():
            raise ValueError(f"{hname} not found in file")
        #set_trace()
        histo = hdict[hname].copy() # process, sys, jmult
    
        xtitle, ytitle, xrebinning, yrebinning, x_lims, y_lims, withData, plot_xlabels, plot_ylabels = variables[hname]
    
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
            #print(*[jmult, hname], sep=", ")
            #pltdir = os.path.join(outdir, args.lepton, jmult)
            #if not os.path.isdir(pltdir):
            #    os.makedirs(pltdir)
    
            #hslice = histo[:, jmult].integrate("jmult")
            #set_trace()
            nosys_histo = histo["nosys", jmult].integrate("jmult").integrate("sys")
            for pdf_sys in sorted(set([key[0] for key in histo.values().keys()])):
                if pdf_sys == "nosys": continue
                print(f"{args.year}, {args.lepton}, {jmult}, {pdf_sys}")
                #print(*[jmult, hname, pdf_sys], sep=", ")
                pltdir = os.path.join(outdir, args.lepton, jmult, pdf_sys)
                if not os.path.isdir(pltdir):
                    os.makedirs(pltdir)
                hslice = histo[pdf_sys, jmult].integrate("jmult").integrate("sys")
    
                    # make 1D projection along dense axes
                for dax in range(2):
                    if dax == 0:
                        xlabel = ytitle
                        xlimits = y_lims
                        figname = os.path.join(pltdir, "_".join([jmult, args.lepton, yaxis_name, "proj"]))
                    else:
                        xlabel = f"{xtitle} [GeV]"
                        xlimits = x_lims
                        figname = os.path.join(pltdir, "_".join([jmult, args.lepton, xaxis_name, "proj"]))
    
                    fig, ax = plt.subplots()
                    fig.subplots_adjust(hspace=.07)
    
                    nosys_hproj = nosys_histo.integrate(nosys_histo.dense_axes()[dax].name)
                    hproj = hslice.integrate(hslice.dense_axes()[dax].name)
                    hproj_masked_vals, hproj_masked_bins = Plotter.get_ratio_arrays(num_vals=hproj.values()[()]-nosys_hproj.values()[()], denom_vals=nosys_hproj.values()[()], input_bins=nosys_hproj.dense_axes()[0].edges())
                    ax.step(hproj_masked_bins, hproj_masked_vals, where="post", **{"color": "r", "linestyle": "-", "label": pdf_sys})
                    ax.legend(loc="upper right", title="PDF")
                    ax.axhline(0, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
                    ax.autoscale()
                    ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.15)
                    ax.set_xlim(xlimits)
                    ax.set_xlabel(xlabel)
                    ax.set_ylabel("Rel. Deviaton from Nominal")
    
                        # add lepton/jet multiplicity label
                    ax.text(
                        0.02, 0.92, "%s, %s" % (objtypes["Lep"][args.lepton], jet_mults[jmult]),
                        fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                    )
                    hep.cms.label(ax=ax, data=False, year=args.year, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))
    
                    #set_trace()
                    fig.savefig(figname)
                    print(f"{figname} written")
                    plt.close()
    
                # plot linearized view 
                fig, ax = plt.subplots()
                fig.subplots_adjust(hspace=.07)
    
                nosys_hline = Plotter.linearize_hist(nosys_histo)
                hline = Plotter.linearize_hist(hslice)
                hline_masked_vals, hline_masked_bins = Plotter.get_ratio_arrays(num_vals=hline.values()[()]-nosys_hline.values()[()], denom_vals=nosys_hline.values()[()], input_bins=nosys_hline.dense_axes()[0].edges())
                ax.step(hline_masked_bins, hline_masked_vals, where="post", **{"color": "r", "linestyle": "-", "label": pdf_sys})
                ax.legend(loc="upper right", title="PDF")
                ax.axhline(0, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
                ax.autoscale()
                ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.15)
                ax.set_xlim((0, hline_masked_bins.size-1))
                ax.set_xlabel(f"{xtitle} $\otimes$ {ytitle}")
                ax.set_ylabel("Rel. Deviaton from Nominal")
    
                    # draw vertical lines separating ctstar bins
                bin_sep_lines = [hslice.values()[()].shape[0]*ybin for ybin in range(1, hslice.values()[()].shape[1])]
                for binline in bin_sep_lines:
                    ax.axvline(binline, color="k", linestyle="--")
                    #if rax is not None: rax.axvline(binline, color="k", linestyle="--")
    
                    # add lepton/jet multiplicity label
                ax.text(
                    0.02, 0.92, "%s, %s" % (objtypes["Lep"][args.lepton], jet_mults[jmult]),
                    fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                )
                hep.cms.label(ax=ax, data=False, year=args.year, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))
    
                #set_trace()
                figname = os.path.join(pltdir, "_".join([jmult, args.lepton, hname]))
                fig.savefig(figname)
                print(f"{figname} written")
                plt.close()
                set_trace()


if args.pca:
        # scale and group hists by process
    for hname in hdict.keys():
        if "cutflow" in hname: continue
        #set_trace()
        hdict[hname].scale(lumi_correction, axis="dataset") # scale hists
        hdict[hname] = hdict[hname].group(process_cat, process, process_groups) # group by process
        hdict[hname] = hdict[hname][:, :, :, args.lepton].integrate("leptype") # only pick out specified lepton

    hname = "mtt_vs_tlep_ctstar_abs"
    if hname not in hdict.keys():
        raise ValueError(f"{hname} not found in file")
    histo = hdict[hname].copy()
    
    xtitle, ytitle, xrebinning, yrebinning, x_lims, y_lims, withData, plot_xlabels, plot_ylabels = variables[hname]
    
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
    for jmult in sorted(set([key[2] for key in histo.values().keys()])):
        print(f"{args.year}, {args.lepton}, {jmult}")
        pltdir = os.path.join(outdir, args.lepton, jmult, "PCA")
        if not os.path.isdir(pltdir):
            os.makedirs(pltdir)
    
        hline = Plotter.linearize_hist(histo[:, :, jmult].integrate("jmult")).integrate("process")
        nosys_raw = hline["nosys"].integrate("sys").values()[()]

        pdf_raw_yields = {pdf: hline[pdf].integrate("sys").values()[()] for pdf in [key[0] for key in hline.values().keys() if "nosys" not in key[0]]}
            # remove 0th and last 2 pdf replicas 0th is base set compatible with 1, last two sets are variations in alpha_S
        del pdf_raw_yields["pdf_0"]
        del pdf_raw_yields["pdf_101"]
        del pdf_raw_yields["pdf_102"]

            # find RMS of total pdf relative deviation 
        set_trace()
        rel_distance = {pdf : (pdf_raw_yields[pdf]-nosys_raw)/nosys_raw for pdf in pdf_raw_yields.keys()} # d_n = r_n - c
        tot_rel_uncs = np.sqrt(np.sum(np.square(np.array([rel_distance[pdf] for pdf in rel_distance.keys()])), axis=0))
        Dsquared_dict = {pdf : np.sum(np.square(pdf_raw_yields[pdf]-nosys_raw)) for pdf in pdf_raw_yields.keys()} 
        epsilon_dict = {pdf : (pdf_raw_yields[pdf]-nosys_raw)/abs(np.sum((pdf_raw_yields[pdf]-nosys_raw))) for pdf in pdf_raw_yields.keys()}
        new_tot_rel_uncs = np.sqrtnp.sum(np.array([Dsquared_dict[pdf]*np.square(epsilon_dict[pdf]) for pdf in pdf_raw_yields.keys()]), axis=0))
        ##tot_rel_uncs = np.sqrt(np.sum(np.square(np.array([rel_distance[pdf] for pdf in rel_distance.keys()])), axis=0)/len(sorted(rel_distance.keys())))
        #rel_distance_X = np.array([rel_distance[pdf] for pdf in rel_distance.keys()])
        #set_trace()

        #    # use PCA to reduce number of pdfs
        #pdf_rel_norm_yields = {pdf: (pdf_raw_yields[pdf]-nosys_raw)/(nosys_raw) for pdf in pdf_raw_yields.keys()}
        ##pdf_norm_yields = {pdf: (pdf_raw_yields[pdf]-nosys_raw)/(abs(sum(pdf_raw_yields[pdf]-nosys_raw))) for pdf in pdf_raw_yields.keys()}
        ##norm_X = np.array([pdf_norm_yields[pdf] for pdf in pdf_norm_yields.keys()])
        #rel_norm_X = np.array([pdf_rel_norm_yields[pdf] for pdf in pdf_rel_norm_yields.keys()])

        ###pca = PCA(n_components=2)
        ##norm_pca = PCA()
        ##norm_pca.fit(norm_X.T)
        ##norm_eigenvectors = norm_pca.components_
        ##norm_eigenvalues = norm_pca.explained_variance_
        ##norm_exp_ratio = norm_pca.explained_variance_ratio_
        ##print(norm_eigenvalues)

        #rel_distance_pca = PCA()
        #rel_distance_pca.fit(rel_distance_X)
        ##rel_distance_pca.fit(rel_distance_X.T)
        #rel_distance_eigenvectors = rel_distance_pca.components_
        #rel_distance_eigenvalues = rel_distance_pca.explained_variance_
        #rel_distance_exp_ratio = rel_distance_pca.explained_variance_ratio_
        #print(rel_distance_eigenvalues)

        #rel_norm_pca = PCA()
        #rel_norm_pca.fit(rel_norm_X.T)
        #rel_norm_eigenvectors = rel_norm_pca.components_
        #rel_norm_eigenvalues = rel_norm_pca.explained_variance_
        #rel_norm_exp_ratio = rel_norm_pca.explained_variance_ratio_
        #print(rel_norm_eigenvalues)

            # try Yiting's method for reducing number of pdfs from AN2016_479 Appendix B
        pdf_pairs = list(itertools.combinations(list(pdf_raw_yields.keys()), 2))
        C_dict = {}
        for pdf_k, pdf_l in pdf_pairs:
            D_k = np.sqrt(np.sum(np.square(rel_distance[pdf_k])))
            D_l = np.sqrt(np.sum(np.square(rel_distance[pdf_l])))
            C_kl = np.sum(rel_distance[pdf_k]*rel_distance[pdf_l])/(D_k*D_l)
            C_dict.update({f"{pdf_k}_{pdf_l}" : C_kl})
        set_trace()

        # plot linearized view 
        fig, ax = plt.subplots()
        fig.subplots_adjust(hspace=.07)
   
        ax.step(hline.dense_axes()[0].edges(), np.r_[tot_rel_uncs, tot_rel_uncs[-1]], where="post", **{"color": "r", "linestyle": "-", "label": "Original 100 replicas"})
        ax.legend(loc="upper right", title="PDF")
        #ax.axhline(0, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
        ax.autoscale()
        ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.15)
        ax.set_xlim((0, hline.dense_axes()[0].edges().size-1))
        ax.set_xlabel(f"{xtitle} $\otimes$ {ytitle}")
        ax.set_ylabel("Rel. Deviaton from Nominal")
    
        #set_trace()
            # draw vertical lines separating ctstar bins
        bin_sep_lines = [(xrebinning.size-1)*ybin for ybin in range(1, yrebinning.size)]
        for binline in bin_sep_lines:
            ax.axvline(binline, color="k", linestyle="--")
    
            # add lepton/jet multiplicity label
        ax.text(
            0.02, 0.92, "%s, %s" % (objtypes["Lep"][args.lepton], jet_mults[jmult]),
            fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
        )
        hep.cms.label(ax=ax, data=False, year=args.year, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))
    
        #set_trace()
        figname = os.path.join(pltdir, "_".join([jmult, args.lepton, "PDF_PCA_comp"]))
        fig.savefig(figname)
        print(f"{figname} written")
        plt.close()
        set_trace()

