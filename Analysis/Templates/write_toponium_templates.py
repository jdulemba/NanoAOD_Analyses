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

from coffea.util import load, save
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
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("template_type", choices=["Raw", "Combined_Era_Lepton", "Combined_Lepton"], help="What type of template do you want to make?")
parser.add_argument("--plot", action="store_true", help="Make plots of sys uncs")
parser.add_argument("--write", action="store_true", help="Write toponium templates to a root file.")
args = parser.parse_args()


def make_year_raw_templates(input_fname, year):
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

        ## make groups based on process
    process = hist.Cat("process", "Process", sorting="placement")
    process_cat = "dataset"
    process_groups = {"Toponium" : ["ToponiumSL", "ToponiumDiLep"]}

    btag_reg_names_dict = btag_sidebands.btag_reg_names_dict[year][btag_wp]

    histo_dict = processor.dict_accumulator({njets : {"Muon" : {}, "Electron" :{}} for njets in njets_to_run})

    for lep in ["Muon", "Electron"]:
        ## make groups and normalize based on processes
        lumi_correction = lumi_corr_dict[year][f"{lep}s"]
            # scale toponium events by mass window cut for lumi correction
        histo = rebin_histo.copy()
            # scale {(ToponiumSL, nosys) : ToponiumSL_nosys}
        histo.scale({key : lumi_correction["_".join(key)] for key in sorted(set([(key[0], key[1]) for key in histo.values().keys()]))}, axis=("dataset", "sys"))
        histo = histo.group(process_cat, process, process_groups)[:, :, :, lep, :, :].integrate("leptype")

        systs = sorted(set([key[1] for key in histo.values().keys()]))

            # loop over each jet multiplicity
        for jmult in njets_to_run:
            linearized_histo = Plotter.linearize_hist(histo[:, :, jmult, btag_reg_names_dict["Signal"]["reg"]].integrate("jmult").integrate("btag"))
                # loop over each systematic
            for sys in systs:
                sys_histo = (linearized_histo[:, sys].integrate("sys")).copy()

                    ## write nominal and systematic variations for each topology to file
                for proc in sorted(set([key[0] for key in sys_histo.values().keys()])):
                    if proc != "Toponium": raise ValueError(f"Only valid proc is Toponium! Not {proc}")
                    if not sys_histo[proc].values().keys():
                        print(f"Systematic {sys} for {lep} {jmult} {proc} not found, skipping")
                        continue

                    print(year, lep, jmult, sys, proc)
                    template_histo = (sys_histo[proc].integrate("process")).copy()
                    sumw, sumw2 = template_histo.values(sumw2=True)[()]
                    rel_err = np.sqrt(sumw2)/np.abs(sumw)
                    rel_err_mask = rel_err > 10
                    if np.any(rel_err_mask):
                        set_trace()

                        ## save template histos to coffea dict
                    histo_dict[jmult][lep][f"{proc}_{sys}"] = template_histo.copy()

    return histo_dict


def make_output_dir(analyzer):
    outdir = os.path.join(eos_dir, "results", jobid, analyzer, "Templates")
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    return outdir
  

def combine_era_lepton_templates(input_dict):
        ## combine templates across years and lepton flavour channels
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

    coffea_out = os.path.join(outdir, f"raw_combined_era_lepton_templates_lj_toponium_{jobid}.coffea")
    save(combined_acc, coffea_out)
    print(f"{coffea_out} written")

    if args.plot:
        make_combined_era_lepton_plots(combined_acc)


## make plots of systematic uncs
def make_combined_era_lepton_plots(hdict):
    plt_outdir = os.path.join(plot_outdir, jobid, "Toponium_HBSR", "Combined_Era_Lepton")
    if not os.path.isdir(plt_outdir):
        os.makedirs(plt_outdir)    

    year_to_use = cfeatures.year_labels["Total"]
    lumi_to_use = (data_lumi["TOT"]["Electrons"] + data_lumi["TOT"]["Muons"])/2000.


        # loop over jmults
    for jmult in hdict.keys():
        nosys = hdict[jmult]["Toponium_nosys"].copy()

        fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True, figsize=(15.0, 10.0))
        fig.subplots_adjust(hspace=.07)
        for sysname in sys_style_dict.keys():
            template = hdict[jmult][f"Toponium_{sysname}"].copy()
            masked_vals, masked_bins = Plotter.get_ratio_arrays(num_vals=template.values()[()], denom_vals=np.ones_like(template.values()[()]), input_bins=template.dense_axes()[0].edges())
            ax.step(masked_bins, masked_vals, where="post", **sys_style_dict[sysname])
            ratio_masked_vals, ratio_masked_bins = Plotter.get_ratio_arrays(num_vals=template.values()[()], denom_vals=nosys.values()[()], input_bins=template.dense_axes()[0].edges())
            rax.step(ratio_masked_bins, ratio_masked_vals, where="post", **sys_style_dict[sysname])

        ax.legend(loc="upper right", title="$\\eta_{t}$", ncol=1)
        ax.autoscale()
        ax.set_ylabel("Events")
        ax.set_ylim(0, ax.get_ylim()[1]*1.3)
        ax.set_xlabel(None)
        ax.set_xlim(x_lims)

        rax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
        rax.autoscale()
        rax.set_ylabel("Ratio to Nominal")
        #rax.set_ylim(max(rax.get_ylim()[0], 0.65), min(rax.get_ylim()[1], 1.25))
        rax.set_xlim(x_lims)
        rax.set_xlabel("$m_{t\\bar{t}}$ [GeV]")

            # add lepton/jet multiplicity label
        ax.text(
            0.02, 0.92, cfeatures.channel_labels[f"Lepton_{jmult}"],
            fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
        )
            ## draw vertical lines for distinguishing different ctstar bins
        vlines = [x_lims[1]*ybin/5 for ybin in range(1, 5)]
        [ax.axvline(vline, color="k", linestyle="--") for vline in vlines]
        [rax.axvline(vline, color="k", linestyle="--") for vline in vlines]

        for idx, label in enumerate(ctstar_binlabels):
            ax.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                xytext=(0, 25), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)

        ax.set_xticks(mtt_bin_inds_to_plot)
        ax.set_xticklabels(mtt_bins_to_plot)
        rax.set_xticks(mtt_bin_inds_to_plot)
        rax.set_xticklabels(mtt_bins_to_plot)
        hep.cms.label(ax=ax, data=False, label="Preliminary", year=year_to_use, lumi=int(round(lumi_to_use, 0)))

        figname = os.path.join(plt_outdir, "_".join(["Toponium", jmult, "Combined_Era_Lepton", "SysComp"]))
        fig.savefig(figname)
        print(f"{figname} written")
        plt.close(fig)


def make_indiv_plots(hdict, year):
    plt_outdir = os.path.join(plot_outdir, jobid, "Toponium_HBSR", "Indiv", year)
    if not os.path.isdir(plt_outdir):
        os.makedirs(plt_outdir)    

    year_to_use = cfeatures.year_labels[year]

        # loop over jmults
    for jmult in hdict.keys():
            # loop over leptons
        for lep in hdict[jmult].keys():
            lumi_to_use = (data_lumi[year][f"{lep}s"])/1000.
    
            nosys = hdict[jmult][lep]["Toponium_nosys"].copy()
    
            fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True, figsize=(15.0, 10.0))
            fig.subplots_adjust(hspace=.07)
            for sysname in sys_style_dict.keys():
                template = hdict[jmult][lep][f"Toponium_{sysname}"].copy()
                masked_vals, masked_bins = Plotter.get_ratio_arrays(num_vals=template.values()[()], denom_vals=np.ones_like(template.values()[()]), input_bins=template.dense_axes()[0].edges())
                ax.step(masked_bins, masked_vals, where="post", **sys_style_dict[sysname])
                ratio_masked_vals, ratio_masked_bins = Plotter.get_ratio_arrays(num_vals=template.values()[()], denom_vals=nosys.values()[()], input_bins=template.dense_axes()[0].edges())
                rax.step(ratio_masked_bins, ratio_masked_vals, where="post", **sys_style_dict[sysname])
    
            ax.legend(loc="upper right", title="$\\eta_{t}$", ncol=1)
            ax.autoscale()
            ax.set_ylabel("Events")
            ax.set_ylim(0, ax.get_ylim()[1]*1.3)
            ax.set_xlabel(None)
            ax.set_xlim(x_lims)
    
            rax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
            rax.autoscale()
            rax.set_ylabel("Ratio to Nominal")
            #rax.set_ylim(max(rax.get_ylim()[0], 0.65), min(rax.get_ylim()[1], 1.25))
            rax.set_xlim(x_lims)
            rax.set_xlabel("$m_{t\\bar{t}}$ [GeV]")
    
                # add lepton/jet multiplicity label
            ax.text(
                0.02, 0.92, f"{cfeatures.channel_labels[f'{lep}_{jmult}']}",
                fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
            )
                ## draw vertical lines for distinguishing different ctstar bins
            vlines = [x_lims[1]*ybin/5 for ybin in range(1, 5)]
            [ax.axvline(vline, color="k", linestyle="--") for vline in vlines]
            [rax.axvline(vline, color="k", linestyle="--") for vline in vlines]
    
            for idx, label in enumerate(ctstar_binlabels):
                ax.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                    xytext=(0, 25), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)
    
            ax.set_xticks(mtt_bin_inds_to_plot)
            ax.set_xticklabels(mtt_bins_to_plot)
            rax.set_xticks(mtt_bin_inds_to_plot)
            rax.set_xticklabels(mtt_bins_to_plot)
            hep.cms.label(ax=ax, data=False, label="Preliminary", year=year_to_use, lumi=round(lumi_to_use, 1))
    
            figname = os.path.join(plt_outdir, "_".join(["Toponium", lep, jmult, "Indiv", year, "SysComp"]))
            fig.savefig(figname)
            print(f"{figname} written")
            plt.close(fig)



def make_raw_tot_templates():
    tot_hdict = {year : {} for year in years_to_run}
    
    # loop over each year to lumi scale systematics
    for year in years_to_run:
        input_dir = os.path.join(eos_dir, "results", jobid, "Toponium_HBSR", year)
        if os.path.isdir(input_dir):
            fnames = fnmatch.filter(os.listdir(input_dir), "*TOT.coffea")
            fnames = [os.path.join(input_dir, fname) for fname in fnames]
            if len(fnames) > 1: raise ValueError("Multiple input files found")
        else: raise ValueError("No input file found.")

        toponium_fname = fnames[0]
        tot_hdict[year] = make_year_raw_templates(toponium_fname, year)
        if args.plot:
            make_indiv_plots(tot_hdict[year], year)
    
    coffea_out = os.path.join(outdir, f"raw_templates_lj_toponium_{jobid}.coffea")
    save(tot_hdict, coffea_out)
    print(f"{coffea_out} written")


def write_rootfile(dist_type):
        ## save templates to root file for each distribution for each year
    if dist_type == "Raw":
        input_tfile = os.path.join(outdir, f"raw_templates_lj_toponium_{jobid}.coffea")
        if not os.path.isfile(input_tfile): raise ValueError(f"{input_tfile} not found")
        tdict = load(input_tfile)
    else:
        raise ValueError("Only 'Raw' allowed as distribution type at this point")

    output_rname = input_tfile.replace(".coffea", ".root")
    rfile = uproot.recreate(output_rname, compression=uproot.ZLIB(4)) if os.path.isfile(output_rname) else uproot.create(output_rname)

    for year in tdict.keys():
        if year == "2016APV": year_to_use = "2016pre"
        elif year == "2016": year_to_use = "2016post"
        else: year_to_use = year

        for jmult in tdict[year].keys():
            for lep in tdict[year][jmult].keys():
                orig_lepdir = "muNJETS" if lep == "Muon" else "eNJETS"
                lepdir = orig_lepdir.replace("NJETS", jmult.lower())

                dirname = f"{lepdir}_{year_to_use}"
                rfile.mkdir(dirname)

                for hname, template in tdict[year][jmult][lep].items():
                    if "nosys" in hname:
                        rfile[dirname]["EtaT"] = (template.copy()).to_hist()
                    elif "BindingEnergy" in hname:
                        rfile[dirname][hname.replace("Toponium", "EtaT").replace("BindingEnergy", "Eb")] = (template.copy()).to_hist()
                    elif "TopMass" in hname:
                        rfile[dirname][hname.replace("Toponium", "EtaT").replace("TopMass", "tmass_EtaT")] = (template.copy()).to_hist()
                    else:
                        print(hname)
                        set_trace()

    rfile.close()
    print(f"{output_rname} written")




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
    lumi_corr_dict = load(os.path.join(proj_dir, "Corrections", base_jobid, f"{lumi_name}.coffea"))

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

    ctstar_binlabels = [r"%s $\in [%s, %s)$" % (cfeatures.variable_names_to_labels["tlep_ctstar_abs"], linearize_binning[1][bin], linearize_binning[1][bin+1]) for bin in range(len(linearize_binning[1])-1)]
    ctstar_binlabels[-1] = r"%s $\in [%s, %s]$" % (cfeatures.variable_names_to_labels["tlep_ctstar_abs"], linearize_binning[1][-2], linearize_binning[1][-1])
    ctstar_bin_locs = np.linspace((len(linearize_binning[0])-1)/2, (len(linearize_binning[0])-1)*(len(linearize_binning[1])-1) - (len(linearize_binning[0])-1)/2, len(linearize_binning[1])-1)
    mtt_vals_to_plot = np.array([400, 600, 1000])
    mtt_tiled_labels = np.tile(linearize_binning[0][:-1], linearize_binning[1].size-1)
    mtt_bin_inds_to_plot = np.where(np.in1d(mtt_tiled_labels, mtt_vals_to_plot))[0]
    mtt_bins_to_plot = np.tile(mtt_vals_to_plot, linearize_binning[1].size-1)

    x_lims = (0, (linearize_binning[1].size - 1)* (linearize_binning[0].size - 1))

    sys_style_dict = {
        "nosys"             : {"color" : "k", "label" : "Nominal"},
        "BindingEnergyUp"   : {"color" : "#e41a1c", "label" : "Binding Energy Up"}, # red
        "BindingEnergyDown" : {"color" : "#377eb8", "label" : "Binding Energy Down"}, #blue 
        "TopMassUp"         : {"color" : "#ff7f00", "label" : "Top Mass Up"}, # orange
        "TopMassDown"       : {"color" : "#4daf4a", "label" : "Top Mass Down"}, # green
    }

    if args.plot: data_lumi = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())

    outdir = make_output_dir("Toponium_HBSR")

    if args.template_type == "Raw":
        make_raw_tot_templates()
        write_rootfile("Raw")

    if args.template_type == "Combined_Era_Lepton":
        #set_trace()
            # make sure raw template file exists
        raw_tfile = os.path.join(outdir, f"raw_templates_lj_toponium_{jobid}.coffea")
        if not os.path.isfile(raw_tfile):
            make_raw_tot_templates()

        combine_era_lepton_templates(load(raw_tfile))
        #combine_era_lepton_templates(tot_hdict)

    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
