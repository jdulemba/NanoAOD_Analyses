#! /bin/env python

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
import Utilities.plot_tools as plt_tools
from coffea import hist
import numpy as np
import fnmatch
import Utilities.Plotter as Plotter
import coffea.processor as processor    
import Utilities.final_analysis_binning as final_binning
import Utilities.btag_sideband_regions as btag_sidebands
import Utilities.prettyjson as prettyjson
import uproot
import Utilities.common_features as cfeatures

import argparse
class ParseKwargs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, dict())
        for value in values:
            key, value = value.split("=")
            getattr(namespace, self.dest)[key] = value


from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--nomSMTTxsec", action="store_true", help="Apply nominal SM cross sections to top mass and LHE scale weights")
parser.add_argument("--plot", action="store_true", help="Make plots of sys uncs")
args = parser.parse_args()


def get_bkg_templates(year):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """
    #set_trace()
    hdict = load(input_fname)

        # get correct hist and rebin
    hname_to_use = "mtt_vs_tlep_ctstar_abs"
    if hname_to_use not in hdict.keys():
        raise ValueError(f"{hname_to_use} not found in file")
    xrebinning, yrebinning = linearize_binning
    histo = hdict[hname_to_use][Plotter.nonsignal_samples] # process, sys, jmult, leptype, btag, lepcat
    
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
    rebin_histo = histo.rebin(yaxis_name, new_ybins)

    #set_trace()    
        ## scale ttJets events, split by reconstruction type, by normal ttJets lumi correction
    ttJets_permcats = ["*right", "*matchable", "*unmatchable", "*sl_tau", "*other"]
    names = [dataset for dataset in sorted(set([key[0] for key in histo.values().keys()]))] # get dataset names in hists
    ttJets_cats = [name for name in names if any([fnmatch.fnmatch(name, cat) for cat in ttJets_permcats])] # gets ttJets(_PS)_other, ...

        ## make groups based on process
    process = hist.Cat("process", "Process", sorting="placement")
    process_cat = "dataset"

    histo_dict = processor.dict_accumulator({njets : {"Muon" : {}, "Electron" :{}} for njets in njets_to_run})

    for lep in ["Muon", "Electron"]:
        #set_trace()    
        ## make groups based on process
        process_groups = plt_tools.make_dataset_groups(lep, year, samples=names, bkgdict="templates")
        
        lumi_correction = lumi_corr_dict[year][f"{lep}s"]
                # scale ttJets events, split by reconstruction type, by normal ttJets lumi correction
        if len(ttJets_cats) > 0:
            for tt_cat in ttJets_cats:
                ttJets_lumi_topo = "_".join(tt_cat.split("_")[:-2]) if "sl_tau" in tt_cat else "_".join(tt_cat.split("_")[:-1]) # gets ttJets[SL, Had, DiLep] or ttJets_PS
                ttJets_eff_lumi = lumi_correction[ttJets_lumi_topo]
                lumi_correction.update({tt_cat: ttJets_eff_lumi})
    
        histo = rebin_histo.copy()
        #set_trace()
        histo.scale(lumi_correction, axis="dataset")
        for nom_tt in ["ttJetsSL", "ttJetsDiLep", "ttJetsHad"]:
            for cat in ttJets_permcats:
                # rescale LHEscale systematics correctly
                for hname, dname in tt_LHEscale_wts_name_dict.items():
                    if f"{nom_tt}_{cat[1:]}" not in ttJets_cats: continue
                    #print(f"{nom_tt}_{cat[1:]}_{dname}")
                    lhe_scale = lumi_correction[f"{nom_tt}_{dname}"]/lumi_correction[f"{nom_tt}_{cat[1:]}"]
                    histo.scale({(f"{nom_tt}_{cat[1:]}", hname) : lhe_scale}, axis=("dataset", "sys"))

        histo = histo.group(process_cat, process, process_groups)[:, :, :, lep, :, :].integrate("leptype")

        systs = sorted(set([key[1] for key in histo.values().keys()]))
        #set_trace()

            # loop over each jet multiplicity
        for jmult in njets_to_run:
                ## open and save 'nosys' histo from existing templates file
            histo_dict[jmult][lep]["TT_nosys"] = existing_template_file[jmult][lep]["TT_nosys"].copy()

            linearized_histo = Plotter.linearize_hist(histo[:, :, jmult, btag_reg_names_dict["Signal"]["reg"]].integrate("jmult").integrate("btag"))
                # loop over each systematic
            for sys in systs:
                sys_histo = (linearized_histo[:, sys].integrate("sys")).copy()

                    ## write nominal and systematic variations for each topology to file
                for proc in sorted(set([key[0] for key in sys_histo.values().keys()])):
                    if proc != "TT": raise ValueError(f"Only valid proc is TT! Not {proc}")
                    if not sys_histo[proc].values().keys():
                        print(f"Systematic {sys} for {lep} {jmult} {proc} not found, skipping")
                        continue

                    print(year, lep, jmult, sys, proc)
                    #set_trace()
                    template_histo = (sys_histo[proc].integrate("process")).copy()
                    sumw, sumw2 = template_histo.values(sumw2=True)[()]
                    rel_err = np.sqrt(sumw2)/np.abs(sumw)
                    rel_err_mask = rel_err > 10
                    if np.any(rel_err_mask):
                        set_trace()

                        ## save template histos to coffea dict
                    histo_dict[jmult][lep][f"{proc}_{sys}"] = template_histo.copy()

    #set_trace()
    return histo_dict


def make_output_dir(analyzer, year):
    outdir = os.path.join(eos_dir, "results", f"{year}_{jobid}", f"Templates_{analyzer}")
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    return outdir
  

def combine_templates(input_dict):
        ## combine templates across years and lepton flavour channels
    #set_trace()
    combined_acc = processor.dict_accumulator({njets : {} for njets in njets_to_run})
    for jmult in njets_to_run:
            # find common names across lepton channels
        mu_tnames_dict = {year : sorted(input_dict[year][jmult]["Muon"].keys()) for year in input_dict.keys() }
        el_tnames_dict = {year : sorted(input_dict[year][jmult]["Electron"].keys()) for year in input_dict.keys() }
        tnames = sorted( set.intersection(*map(set, sorted(mu_tnames_dict.values()))) & set.intersection(*map(set, sorted(el_tnames_dict.values()))) )
        for tname in tnames:
            print(jmult, tname)

            combined_year_template = None
            idx = 0
            for year in input_dict.keys():
                for lep in sorted(["Muon", "Electron"]):
                    combined_year_template = input_dict[year][jmult][lep][tname].copy() if idx == 0 else combined_year_template.add(input_dict[year][jmult][lep][tname].copy())
                    idx += 1

                ## save template histos to coffea dict
            combined_acc[jmult][tname] = combined_year_template.copy()

    if args.plot:
        make_plots(combined_acc)

        ## save templates to root file for each distribution for each year
    #set_trace()
    rname = "test.root"
    rfile = uproot.recreate(rname, compression=uproot.ZLIB(4)) if os.path.isfile(rname) else uproot.create(rname)

    for jmult in combined_acc.keys():
        for hname, template in combined_acc[jmult].items():
            for year in input_dict.keys():
                if year == "2016APV": year_to_use = "2016pre"
                elif year == "2016": year_to_use = "2016post"
                else: year_to_use = year

                for lep in input_dict[year][jmult].keys():
                    orig_lepdir = "muNJETS" if lep == "Muon" else "eNJETS"
                    lepdir = orig_lepdir.replace("NJETS", jmult.lower())

                    dirname = f"{lepdir}_{year_to_use}"
                    rfile.mkdir(dirname)

                    #set_trace()
                    if "nosys" in hname:
                        continue
                        rfile[dirname]["TT"] = (input_dict[year][jmult][lep][hname].copy()).to_hist()
                        continue

                    elif ("NNLOqcd" in hname) or ("LOqcd" in hname):
                    #elif "NNLOqcd" in hname:
                        #set_trace()
                        rfile[dirname][hname] = (input_dict[year][jmult][lep][hname].copy()).to_hist()
                        continue

                        # full EWK_scheme unc wts applied
                    elif "EWschemeYt" in hname:
                        continue
                        #set_trace()
                        output_histo = template.copy()
                        sys_vals = np.copy(output_histo.values()[()]/(combined_acc[jmult]["TT_nosys"].copy()).values()[()] * (input_dict[year][jmult][lep]["TT_nosys"].copy()).values()[()])
                        variances = np.zeros_like(sys_vals)
                        output_histo.values(sumw2=True)[()][0][:], output_histo.values(sumw2=True)[()][1][:] = sys_vals, variances

                        # individual dQCD and dEW wts applied
                    else:
                        continue
                        #set_trace()
                        output_histo = template.copy()
                        sys_vals = np.copy(output_histo.values()[()]/(combined_acc[jmult]["TT_nosys"].copy()).values()[()])
                        variances = np.zeros_like(sys_vals)
                        output_histo.values(sumw2=True)[()][0][:], output_histo.values(sumw2=True)[()][1][:] = sys_vals, variances


                    rfile[dirname][hname] = (output_histo.copy()).to_hist()

    rfile.close()
    print(f"{rname} written")


## make plots of systematic uncs
def make_plots(hdict):
    outdir = os.path.join(plot_outdir, jobid, "Yukawa_EWKscheme_Uncs")
    if not os.path.isdir(outdir):
        os.makedirs(outdir)    

    year_to_use = cfeatures.year_labels["Total"]

    ctstar_binlabels = [r"%s $\in [%s, %s)$" % (cfeatures.variable_names_to_labels["tlep_ctstar_abs"], linearize_binning[1][bin], linearize_binning[1][bin+1]) for bin in range(len(linearize_binning[1])-1)]
    ctstar_binlabels[-1] = r"%s $\in [%s, %s]$" % (cfeatures.variable_names_to_labels["tlep_ctstar_abs"], linearize_binning[1][-2], linearize_binning[1][-1])
    ctstar_bin_locs = np.linspace((len(linearize_binning[0])-1)/2, (len(linearize_binning[0])-1)*(len(linearize_binning[1])-1) - (len(linearize_binning[0])-1)/2, len(linearize_binning[1])-1)
    mtt_vals_to_plot = np.array([400, 600, 1000])
    mtt_tiled_labels = np.tile(linearize_binning[0][:-1], linearize_binning[1].size-1)
    mtt_bin_inds_to_plot = np.where(np.in1d(mtt_tiled_labels, mtt_vals_to_plot))[0]
    mtt_bins_to_plot = np.tile(mtt_vals_to_plot, linearize_binning[1].size-1)

    x_lims = (0, (linearize_binning[1].size - 1)* (linearize_binning[0].size - 1))

    EWscheme_style_dict = {
        "EWschemeYt0p88" : {"color" : "b", "label" : "$y_t$ = 0.88"},
        "EWschemeYt1p0"  : {"color" : "k", "label" : "$y_t$ = 1.0"},
        "EWschemeYt1p11" : {"color" : "r", "label" : "$y_t$ = 1.11"},
    }
    dEW_style_dict = {
        "dEWyt0p88" : {"color" : "b", "label" : "$y_t$ = 0.88"},
        "dEWyt1p0"  : {"color" : "k", "label" : "$y_t$ = 1.0"},
        "dEWyt1p11" : {"color" : "r", "label" : "$y_t$ = 1.11"},
    }
    dQCD_style_dict = {
        "dQCD" : {"color" : "r"},
    }
    QCDorders_style_dict = {
        "NNLOqcd" : {"color" : "r", "label" : "NNLO QCD"},
        "LOqcd"   : {"color" : "b", "label" : "LO QCD"},
    }

    sys_groups = {
        "EWK_Scheme" : [["EWschemeYt0p88", "EWschemeYt1p0", "EWschemeYt1p11"], EWscheme_style_dict],
        "deltaEW" : [["dEWyt0p88", "dEWyt1p0", "dEWyt1p11"], dEW_style_dict],
        "deltaQCD" : [["dQCD"], dQCD_style_dict],
        "QCDorders" : [["NNLOqcd", "LOqcd"], QCDorders_style_dict],
    }
    #set_trace()
        # loop over jmults
    for jmult in hdict.keys():
        pltdir = os.path.join(outdir, jmult)
        if not os.path.isdir(pltdir):
            os.makedirs(pltdir)

        nosys = hdict[jmult]["TT_nosys"].copy()

        for sysname, (sys_list, style_dict) in sys_groups.items():
            #set_trace()
            fig, ax = plt.subplots(figsize=(15.0, 10.0))
            fig.subplots_adjust(hspace=.07)

            for sys in sys_list:
                template = hdict[jmult][f"TT_{sys}"].copy() if (("NNLOqcd" in sys) or ("LOqcd" in sys)) else hdict[jmult][f"TT_{sys}Up"].copy()
                masked_vals, masked_bins = Plotter.get_ratio_arrays(num_vals=template.values()[()], denom_vals=nosys.values()[()], input_bins=template.dense_axes()[0].edges())
                ax.step(masked_bins, masked_vals, where="post", **style_dict[sys])

            ax.legend(loc="upper right", title="TT", ncol=1)

            if sysname == "EWK_Scheme":
                ax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
                ax.autoscale()
                ax.set_ylabel("Ratio to Nominal")
                ax.set_ylim(max(ax.get_ylim()[0], 0.65), min(ax.get_ylim()[1], 1.25))
            else:
                if sysname == "deltaEW" :
                    ax.axhline(0, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
                    ax.set_ylabel("Relative Deviation")
                if sysname == "QCDorders" : 
                    ax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
                    ax.set_ylabel("Ratio to Nominal")
                ax.autoscale()

            ax.set_xlim(x_lims)
            ax.set_xlabel("$m_{t\\bar{t}}$ [GeV]")

                # add lepton/jet multiplicity label
            ax.text(
                0.02, 0.92, sysname,
                fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
            )
                ## draw vertical lines for distinguishing different ctstar bins
            vlines = [x_lims[1]*ybin/5 for ybin in range(1, 5)]
            [ax.axvline(vline, color="k", linestyle="--") for vline in vlines]

            for idx, label in enumerate(ctstar_binlabels):
                ax.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                    xytext=(0, 25), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)

            ax.set_xticks(mtt_bin_inds_to_plot)
            ax.set_xticklabels(mtt_bins_to_plot)
            hep.cms.label(ax=ax, data=False, year=year_to_use)

            figname = os.path.join(pltdir, "_".join(["TT", jmult, sysname, "Yt_Comp" if "EW" in sysname  else "Comp"]))
            fig.savefig(figname)
            print(f"{figname} written")
            plt.close(fig)




if __name__ == "__main__":
    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]
    base_jobid = os.environ["base_jobid"]
    eos_dir = os.environ["eos_dir"]
    plot_outdir = os.environ["plots_dir"]

    njets_to_run = sorted(["3Jets", "4PJets"])
    years_to_run = ["2016APV", "2016", "2017", "2018"]

        ## initialize lumi scaling files
    lumi_name = "MC_LumiWeights"
    if args.nomSMTTxsec:
        lumi_name += "_nomSMTTxsec"
    lumi_corr_dict = load(os.path.join(proj_dir, "Corrections", base_jobid, f"{lumi_name}.coffea"))
  
    tt_LHEscale_wts_name_dict = {
        "FACTORDown" : "uF_down",
        "FACTORUp"   : "uF_up",
        "RENORMDown" : "uR_down",
        "RENORMUp"   : "uR_up",
    }

    jet_pars = prettyjson.loads(open(os.path.join(proj_dir, "cfg_files", f"cfg_pars_{jobid}.json")).read())["Jets"]
    btagger = jet_pars["btagger"]
    wps_to_use = list(set([jet_pars["permutations"]["tightb"],jet_pars["permutations"]["looseb"]]))
    if not( len(wps_to_use) == 1):
        raise IOError("Only 1 unique btag working point supported now")
    btag_wp = btagger+wps_to_use[0]

    linearize_binning = (
        final_binning.mtt_binning,
        final_binning.ctstar_abs_binning
    ) 

    #set_trace()
    analyzer = "htt_btag_sb_regions"

    tot_hdict = {year : {} for year in years_to_run}

    # loop over each year to lumi scale systematics
    for year in years_to_run:
            ## open existing template file to copy nosys
        outdir = make_output_dir(analyzer, year)
        existing_template_file = load(os.path.join(outdir, f"raw_templates_lj_bkg_nomSMTTxsec_{year}_{jobid}.coffea" if args.nomSMTTxsec else f"raw_templates_lj_bkg_{year}_{jobid}.coffea"))

        btag_reg_names_dict = btag_sidebands.btag_reg_names_dict[year][btag_wp]
        input_fname = os.path.join(eos_dir, "results",f"{year}_{jobid}", analyzer, f"BATCH_EWKscheme_LOqcd_NNLOqcd_Variations_23August2023_{year}_{jobid}_TOT.coffea")
        #input_fname = os.path.join(eos_dir, "results",f"{year}_{jobid}", analyzer, f"BATCH_EWKschemeVariations_07August2023_{year}_{jobid}_TOT.coffea")
        if not os.path.isfile(input_fname): raise ValueError("No input file found.")
        tot_hdict[year] = get_bkg_templates(year)

    #set_trace()
    combine_templates(tot_hdict)
    #set_trace()

    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
