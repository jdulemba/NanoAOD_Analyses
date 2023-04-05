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
import numpy as np
import Utilities.Plotter as Plotter
import Utilities.final_analysis_binning as final_binning
import Utilities.common_features as cfeatures
import Utilities.prettyjson as prettyjson
import Utilities.systematics as systematics
from Utilities.styles import styles as hstyles

import argparse
class ParseKwargs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, dict())
        for value in values:
            key, value = value.split("=")
            getattr(namespace, self.dest)[key] = value

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("templates_to_run", type=str, help="Choose which type of templates to run, multiple options can be input as ':' separated strings.")
parser.add_argument("combination", choices=["lepton", "era_lepton", "indiv"], help="What type of combination has been performed.")
parser.add_argument("--year", choices=["2016APV", "2016", "2017", "2018"], help="Specify which year to use if combination is only across leptons.")
parser.add_argument("--MEopts", nargs="*", action=ParseKwargs, help="Options to specify which mass and width points to produce for ME reweighted signal.")
args = parser.parse_args()

def plot_bkg_dists():
    if args.combination == "era_lepton":
        outdir = os.path.join(plot_outdir, jobid, "Mtop3GeV_Uncs", "ERA_LEPTON")
    
        year_to_use = cfeatures.year_labels["Total"]
        lumi_to_use = round((data_lumi_dict["TOT"]["Muons"] + data_lumi_dict["TOT"]["Electrons"])/2000.,0)
        systs = [
            "FSRUp", "FSRDown",
        ]
        mtop3gev_hdict = load(os.path.join(mtop3gev_dir, "raw_combined_year_and_lepton_templates.coffea"))
        nominal_hdict = load(os.path.join(eos_dir, "results", jobid, "Templates_htt_btag_sb_regions", f"raw_combined_year_and_lepton_templates_lj_bkg_{jobid}.coffea"))
    
    
    if args.combination == "lepton":
        #set_trace()
        year_to_use = cfeatures.year_labels[args.year]
        lumi_to_use = round((data_lumi_dict[args.year]["Muons"] + data_lumi_dict[args.year]["Electrons"])/2000.,1)
        systs = [
            "JES_FlavorQCDOnlyLightJets_DW", "JES_FlavorQCDOnlyLightJets_UP",
        ]
    
        outdir = os.path.join(plot_outdir, jobid, "Mtop3GeV_Uncs", "LEPTON", args.year)
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
    
        mtop3gev_hdict = load(os.path.join(mtop3gev_dir, f"raw_combined_lep_templates_{args.year}.coffea"))
        nominal_hdict = load(os.path.join(eos_dir, "results", f"{args.year}_{jobid}", "Templates_htt_btag_sb_regions", f"raw_combined_lep_templates_lj_bkg_{args.year}_{jobid}.coffea"))
    
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    
    
    #set_trace()
    for sys in systs:
        syslabel = systematics.combine_template_sys_to_name[args.year][sys] if args.combination == "lepton" else systematics.combine_template_sys_to_name["2017"][sys]
        for jmult in ["3Jets", "4PJets"]:
            fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True, figsize=(15.0, 10.0))
            fig.subplots_adjust(hspace=.07)
    
            mt_up_nosys, mt_up_sys = mtop3gev_hdict[jmult][f"TT_mtop1755_nosys"].copy(), mtop3gev_hdict[jmult][f"TT_mtop1755_{sys}"].copy()
            mt_dw_nosys, mt_dw_sys = mtop3gev_hdict[jmult][f"TT_mtop1695_nosys"].copy(), mtop3gev_hdict[jmult][f"TT_mtop1695_{sys}"].copy()
            mt_cen_nosys, mt_cen_sys = nominal_hdict[jmult][f"TT_nosys"].copy(), nominal_hdict[jmult][f"TT_{sys}"].copy()
    
            # calculate sys/nosys ratios for each mass
                # mt=172.5
            cen_masked_vals, cen_masked_bins = Plotter.get_ratio_arrays(num_vals=mt_cen_sys.values()[()], denom_vals=mt_cen_nosys.values()[()], input_bins=mt_cen_nosys.dense_axes()[0].edges())
            ax.step(cen_masked_bins, cen_masked_vals, where="post", **mtop_styles["TT"])
                # mt=175.5
            up_masked_vals, up_masked_bins = Plotter.get_ratio_arrays(num_vals=mt_up_sys.values()[()], denom_vals=mt_up_nosys.values()[()], input_bins=mt_up_nosys.dense_axes()[0].edges())
            ax.step(up_masked_bins, up_masked_vals, where="post", **mtop_styles["TT_mtop1755"])
            rax.step(up_masked_bins, up_masked_vals/cen_masked_vals, where="post", **mtop_styles["TT_mtop1755"])
                # mt=169.5
            dw_masked_vals, dw_masked_bins = Plotter.get_ratio_arrays(num_vals=mt_dw_sys.values()[()], denom_vals=mt_dw_nosys.values()[()], input_bins=mt_dw_nosys.dense_axes()[0].edges())
            ax.step(dw_masked_bins, dw_masked_vals, where="post", **mtop_styles["TT_mtop1695"])
            rax.step(dw_masked_bins, dw_masked_vals/cen_masked_vals, where="post", **mtop_styles["TT_mtop1695"])
    
            #ax.set_yscale("log")
            ax.legend(loc="upper right", title=mtop_label)
            ax.autoscale()
            #ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.01)
            ax.set_xlim(x_lims)
            ax.set_xlabel(None)
            ax.set_ylabel("Ratio to Nominal")
            #ax.set_ylabel("Events")
            
            rax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
            rax.autoscale()
            rax.set_xlim(x_lims)
            rax.set_xlabel("$m_{t\\bar{t}}$ [GeV]")
            rax.set_ylabel("Ratio to $m_{t}=172.5$")
    
                # add lepton/jet multiplicity label
            ax.text(
                0.02, 0.86, f"{syslabel}\n{cfeatures.channel_labels[f'Lepton_{jmult}']}",
                fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
            )
                ## draw vertical lines for distinguishing different ctstar bins
            vlines = [x_lims[1]*ybin/5 for ybin in range(1, 5)]
            [ax.axvline(vline, color="k", linestyle="--") for vline in vlines]
            [rax.axvline(vline, color="k", linestyle="--") for vline in vlines]
    
            [ax.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                    xytext=(0, 25), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0) for idx, label in enumerate(ctstar_binlabels)]
    
            ax.set_xticks(mtt_bin_inds_to_plot)
            ax.set_xticklabels(mtt_bins_to_plot)
    
            hep.cms.label(ax=ax, label="Preliminary", data=False, year=year_to_use, lumi=lumi_to_use)
            
            pltdir = os.path.join(outdir, jmult)
            if not os.path.isdir(pltdir):
                os.makedirs(pltdir)
    
            figname = os.path.join(pltdir, f"{jmult}_{sys}_MtopComp")
            fig.savefig(figname)
            print(f"{figname} written")
            plt.close(fig)


def plot_sig_dists():
    outdir = os.path.join(plot_outdir, jobid, "Mtop3GeV_Uncs", "MEsig", args.year, "Indiv")
    year_to_use = cfeatures.year_labels[args.year]
    lumi_to_use = round((data_lumi_dict[args.year]["Muons"] + data_lumi_dict[args.year]["Electrons"])/2000.,1)

    orig_mtop_histo = load(os.path.join(input_dir, "BATCH_MEsig_MTopUncs_relw0p5TO25p0_16March2023_2018_Summer20UL_DeepJet.coffea"))["mtt_vs_tlep_ctstar_abs"]
    nominal_histo = load(os.path.join(eos_dir, "results", f"{args.year}_{jobid}", "Templates_htt_btag_sb_regions", f"raw_templates_lj_MEsig_{args.year}_{jobid}.coffea"))
    
    #set_trace()
        # get signal names that are allowed based on inputs
    names = sorted(set([id.name for id in orig_mtop_histo.axis("dataset").identifiers()]))
    names = [name for name in names for mass in allowed_masses for width in allowed_widths if f"{mass}_{width}" in name]
    process_groups = plt_tools.make_dataset_groups("Muon", args.year, samples=names, bkgdict="templates", sigdict="MEreweight_combined")

    mtop_histo = orig_mtop_histo.group("dataset", hist.Cat("process", "Process", sorting="placement"), process_groups)
    signals = sorted(set([key[0] for key in mtop_histo.values().keys()]))
    
    for lep in ["Muon", "Electron"]:
        for jmult in ["3Jets", "4PJets"]:
            for sig in signals:
                fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True, figsize=(15.0, 10.0))
                fig.subplots_adjust(hspace=.07)
   
                nom = nominal_histo[jmult][lep][f"{sig}_nosys"].copy()
                mt_up = Plotter.linearize_hist(mtop_histo[sig, "mtop173p5", jmult, lep].integrate("jmult").integrate("process").integrate("sys").integrate("leptype"))
                mt_dw = Plotter.linearize_hist(mtop_histo[sig, "mtop171p5", jmult, lep].integrate("jmult").integrate("process").integrate("sys").integrate("leptype"))
    
                # calculate sys/nosys ratios for each mass
                    # mt=172.5
                hep.plot.histplot(nom.values()[()], nom.dense_axes()[0].edges(), ax=ax, histtype="step", **mtop_styles["nominal"])
                    # mt=173.5
                hep.plot.histplot(mt_up.values()[()], mt_up.dense_axes()[0].edges(), ax=ax, histtype="step", **mtop_styles["mtop173p5"])
                up_masked_vals, up_masked_bins = Plotter.get_ratio_arrays(num_vals=mt_up.values()[()], denom_vals=nom.values()[()], input_bins=mt_up.dense_axes()[0].edges())
                rax.step(up_masked_bins, up_masked_vals, where="post", **mtop_styles["mtop173p5"])
                    # mt=171.5
                hep.plot.histplot(mt_dw.values()[()], mt_dw.dense_axes()[0].edges(), ax=ax, histtype="step", **mtop_styles["mtop171p5"])
                dw_masked_vals, dw_masked_bins = Plotter.get_ratio_arrays(num_vals=mt_dw.values()[()], denom_vals=nom.values()[()], input_bins=mt_dw.dense_axes()[0].edges())
                rax.step(dw_masked_bins, dw_masked_vals, where="post", **mtop_styles["mtop171p5"])
    
                #ax.set_yscale("log")
                ax.legend(loc="upper right", title=hstyles[sig]["name"])
                ax.autoscale()
                #ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.01)
                ax.set_xlim(x_lims)
                ax.set_xlabel(None)
                ax.set_ylabel("Events")
                
                rax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
                rax.autoscale()
                rax.set_xlim(x_lims)
                rax.set_xlabel("$m_{t\\bar{t}}$ [GeV]")
                rax.set_ylabel("Ratio to Nominal")
                rax.set_ylim(0.5, 1.5)
    
                    # add lepton/jet multiplicity label
                ax.text(
                    0.02, 0.86, f"tmass_AH\n{cfeatures.channel_labels[f'{lep}_{jmult}']}",
                    fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                )
                    ## draw vertical lines for distinguishing different ctstar bins
                vlines = [x_lims[1]*ybin/5 for ybin in range(1, 5)]
                [ax.axvline(vline, color="k", linestyle="--") for vline in vlines]
                [rax.axvline(vline, color="k", linestyle="--") for vline in vlines]
    
                [ax.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                        xytext=(0, 25), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0) for idx, label in enumerate(ctstar_binlabels)]
    
                ax.set_xticks(mtt_bin_inds_to_plot)
                ax.set_xticklabels(mtt_bins_to_plot)
    
                hep.cms.label(ax=ax, label="Preliminary", data=False, year=year_to_use, lumi=lumi_to_use)
                
                pltdir = os.path.join(outdir, lep, jmult, sig.split("_Int")[0] if "Int" in sig else sig.split("_Res")[0])
                if not os.path.isdir(pltdir):
                    os.makedirs(pltdir)

                figname = os.path.join(pltdir, "_".join([sig, jmult, lep, "MtopComp"]))
                fig.savefig(figname)
                print(f"{figname} written")
                plt.close(fig)
                #set_trace()




if __name__ == "__main__":
    allowed_template_options = ["bkg", "MEsig"]
    templates_to_run = [template for template in (args.templates_to_run).split(":") if template in allowed_template_options]
    templates_to_not_run = [template for template in (args.templates_to_run).split(":") if template not in allowed_template_options]
    if templates_to_not_run:
        print(f"{templates_to_not_run} are not valid options for making templates, will be skipped")

    proj_dir = os.environ["PROJECT_DIR"]
    base_jobid = os.environ["base_jobid"]
    jobid = os.environ["jobid"]
    analyzer = "htt_btag_sb_regions"
    plot_outdir = os.environ["plots_dir"]
    eos_dir = os.environ["eos_dir"]
    
    linearize_binning = (
        final_binning.mtt_binning,
        final_binning.ctstar_abs_binning
    )
    
    nbins = (len(linearize_binning[0]) - 1) * (len(linearize_binning[1]) - 1)
    x_lims = (0, nbins)
    
    ctstar_binlabels = [r"%s $\in [%s, %s)$" % (cfeatures.variable_names_to_labels["tlep_ctstar_abs"], linearize_binning[1][bin], linearize_binning[1][bin+1]) for bin in range(len(linearize_binning[1])-1)]
    ctstar_binlabels[-1] = r"%s $\in [%s, %s]$" % (cfeatures.variable_names_to_labels["tlep_ctstar_abs"], linearize_binning[1][-2], linearize_binning[1][-1])
    ctstar_bin_locs = np.linspace((len(linearize_binning[0])-1)/2, (len(linearize_binning[0])-1)*(len(linearize_binning[1])-1) - (len(linearize_binning[0])-1)/2, len(linearize_binning[1])-1)
    mtt_vals_to_plot = np.array([400, 600, 1000])
    mtt_tiled_labels = np.tile(linearize_binning[0][:-1], linearize_binning[1].size-1)
    mtt_bin_inds_to_plot = np.where(np.in1d(mtt_tiled_labels, mtt_vals_to_plot))[0]
    mtt_bins_to_plot = np.tile(mtt_vals_to_plot, linearize_binning[1].size-1)
    
    mtop_label = "$m_{t}$ [GeV]"
    data_lumi_dict = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())
        

    if "bkg" in templates_to_run:
        if (args.combination == "lepton") and (not args.year):
            raise ValueError("Year must be specified when the combination is only across leptons")
    
        mtop_styles = {
            "TT" : {"label" : "172.5", "color" : "k", "linestyle" : "-"},
            "TT_mtop1755" : {"label" : "175.5", "color" : "#e41a1c", "linestyle" : "-"}, # red
            "TT_mtop1695" : {"label" : "169.5", "color" : "#377eb8", "linestyle" : "-"}, # blue
        }

        mtop3gev_dir = os.path.join(eos_dir, "results", jobid, "Mtop3GeV_Uncs", "Templates")

        plot_bkg_dists()    
    
    if "MEsig" in templates_to_run:
        from coffea import hist
        import Utilities.plot_tools as plt_tools

        possible_masses = [365, 380, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000]
        possible_widths = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 8.0, 10.0, 13.0, 15.0, 18.0, 21.0, 25.0]
        possible_masses = [str(mass) for mass in possible_masses]
        possible_widths = [str(width) for width in possible_widths]

        masses_to_run = "All" if not args.MEopts else args.MEopts.get("allowed_masses", "All")
        widths_to_run = "All" if not args.MEopts else args.MEopts.get("allowed_widths", "All")

        if masses_to_run == "All":
            allowed_masses = possible_masses
        else:
            allowed_masses = masses_to_run.split(":")
            allowed_masses = [mass for mass in allowed_masses if mass in possible_masses]

        if widths_to_run == "All":
            allowed_widths = possible_widths
        else:
            allowed_widths = widths_to_run.split(":")
            allowed_widths = [width for width in allowed_widths if width in possible_widths]

        #set_trace()
        allowed_masses = ["M"+mass for mass in allowed_masses]
        allowed_widths = ["W"+width.replace(".", "p") for width in allowed_widths]
        print(f"Making target signal points for {allowed_masses} and {allowed_widths}")

        mtop_styles = {
            "nominal" : {"label" : "172.5 [GeV]", "color" : "k", "linestyle" : "-"},
            "mtop173p5" : {"label" : "173.5 [GeV]", "color" : "#e41a1c", "linestyle" : "-"}, # red
            "mtop171p5" : {"label" : "171.5 [GeV]", "color" : "#377eb8", "linestyle" : "-"}, # blue
        }
        input_dir = os.path.join(eos_dir, "results", f"{args.year}_{jobid}", "signal_ME_evtReweighting", "RecoLevel")
        plot_sig_dists()    

    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
