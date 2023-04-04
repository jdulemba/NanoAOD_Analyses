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
import Utilities.prettyjson as prettyjson
from coffea import hist
import numpy as np
import Utilities.Plotter as Plotter
import Utilities.systematics as systematics
from coffea.hist import plot
import Utilities.final_analysis_binning as final_binning
import Utilities.common_features as cfeatures
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("lepton", choices=["Electron", "Muon"], help="Choose which lepton to make plots for")
parser.add_argument("process", choices=["bkg", "PDF"], help="Specify which process to use.")
args = parser.parse_args()

year_to_use = cfeatures.year_labels[args.year]

base_jobid = os.environ["base_jobid"]
proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
analyzer = "htt_btag_sb_regions"
plot_outdir = os.environ["plots_dir"]
eos_dir = os.environ["eos_dir"]

#set_trace()
print("Choose version")
version = "V29"
#version = "V27"

outdir = os.path.join(plot_outdir, jobid, f"Templates_{analyzer}", version, args.year,  "FINAL", (args.process).upper(), args.lepton)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

linearize_binning = (
    final_binning.mtt_binning,
    final_binning.ctstar_abs_binning
)

ctstar_binlabels = [r"|$\cos (\theta^{*}_{t_{l}})$| $\in [%s, %s)$" % (linearize_binning[1][bin], linearize_binning[1][bin+1]) for bin in range(len(linearize_binning[1])-1)]
ctstar_binlabels[-1] = r"|$\cos (\theta^{*}_{t_{l}})$| $\in [%s, %s]$" % (linearize_binning[1][-2], linearize_binning[1][-1])
ctstar_bin_locs = np.linspace((len(linearize_binning[0])-1)/2, (len(linearize_binning[0])-1)*(len(linearize_binning[1])-1) - (len(linearize_binning[0])-1)/2, len(linearize_binning[1])-1)
mtt_vals_to_plot = np.array([400, 600, 1000])
mtt_tiled_labels = np.tile(linearize_binning[0][:-1], linearize_binning[1].size-1)
mtt_bin_inds_to_plot = np.where(np.in1d(mtt_tiled_labels, mtt_vals_to_plot))[0]
mtt_bins_to_plot = np.tile(mtt_vals_to_plot, linearize_binning[1].size-1)


baseSys = lambda sys : "_".join(sys.split("_")[:-1])

input_dir = os.path.join(eos_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}")
if args.process == "bkg":
    orig_hdict = load(os.path.join(input_dir, f"raw_templates_lj_bkg_{args.year}_{jobid}.coffea"))
    final_hdict = load(os.path.join(input_dir, "FINAL", f"final_templates_lj_bkg_{args.year}_{jobid}_{version}.coffea"))

if args.process == "PDF":
    input_dir = os.path.join(eos_dir, "results", f"{args.year}_{jobid}", "Templates_htt_pdfUncs")
    orig_hdict = load(os.path.join(input_dir, f"raw_pdf_templates_lj_bkg_{args.year}_{jobid}.coffea"))
    final_hdict = load(os.path.join(input_dir, "FINAL", f"final_pdf_templates_lj_bkg_{args.year}_{jobid}_{version}.coffea"))

data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())[args.year]

nom_styles = {"color":"k", "linestyle":"-", "label":"Nominal"}
orig_styles = {"color":"b", "linestyle":"-", "label":"Original"}


    ## make plots for PDF templates
if args.process == "PDF":
    #set_trace()
    for jmult in final_hdict.keys():
        for lep, hdict in final_hdict[jmult].items():
            if lep != args.lepton: continue
            final_nominal, _ = hdict["TT_nosys"].copy()
            orig_nominal = orig_hdict[jmult][lep]["TT_nosys"].copy()
            x_lims = (0, final_nominal.dense_axes()[0].centers().size)
            vlines = [x_lims[1]*ybin/5 for ybin in range(1, 5)]

            pltdir = os.path.join(outdir, jmult)
            if not os.path.isdir(pltdir):
                os.makedirs(pltdir)

            ## init dict of for highest variations
            pdf_vars_dict = {}

            ### plot alphaS variations ratios
            fig, ax = plt.subplots(figsize=(15.0, 10.0))
            fig.subplots_adjust(hspace=.07)
    
                # get original variations
            orig_alphaSup = orig_hdict[jmult][lep]["TT_alphaSUp"].copy()
            orig_alphaSdw = orig_hdict[jmult][lep]["TT_alphaSDown"].copy()
            orig_alphaSUp_vals, orig_alphaSUp_bins = Plotter.get_ratio_arrays(num_vals=orig_alphaSup.values()[()], denom_vals=orig_nominal.values()[()], input_bins=orig_nominal.dense_axes()[0].edges())
            orig_alphaSDown_vals, orig_alphaSDown_bins = Plotter.get_ratio_arrays(num_vals=orig_alphaSdw.values()[()], denom_vals=orig_nominal.values()[()], input_bins=orig_nominal.dense_axes()[0].edges())
                    # plot
            ax.fill_between(orig_alphaSUp_bins, orig_alphaSUp_vals, y2=1., facecolor="r", step="post", alpha=0.5)
            ax.fill_between(orig_alphaSDown_bins, orig_alphaSDown_vals, y2=1., facecolor="b", step="post", alpha=0.5)

                # get final variations
            final_alphaSup, _ = hdict["TT_alphaSUp"].copy()
            final_alphaSdw, _ = hdict["TT_alphaSDown"].copy()
            final_alphaSUp_vals, final_alphaSUp_bins = Plotter.get_ratio_arrays(num_vals=final_alphaSup.values()[()], denom_vals=final_nominal.values()[()], input_bins=final_nominal.dense_axes()[0].edges())
            final_alphaSDown_vals, final_alphaSDown_bins = Plotter.get_ratio_arrays(num_vals=final_alphaSdw.values()[()], denom_vals=final_nominal.values()[()], input_bins=final_nominal.dense_axes()[0].edges())
                    # plot                
            ax.step(final_alphaSUp_bins, final_alphaSUp_vals, where="post", **{"color": "r", "linestyle": "-", "label": f"Final Up"})
            ax.step(final_alphaSDown_bins, final_alphaSDown_vals, where="post", **{"color": "b", "linestyle": "-", "label": f"Final Down"})

            ax.legend(loc="upper right", title="TT", ncol=2)
            ax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
            ax.autoscale()
            ax.set_xlim(x_lims)
            ax.set_xlabel("$m_{t\\bar{t}}$ [GeV]")
            ax.set_ylabel("Ratio to Nominal")

                # draw vertical lines separating ctstar bins
            [ax.axvline(vline, color="k", linestyle="--") for vline in vlines]

            [ax.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                xytext=(0, 400), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0) for idx, label in enumerate(ctstar_binlabels)]

            #set_trace()
            ax.set_xticks(mtt_bin_inds_to_plot)
            ax.set_xticklabels(mtt_bins_to_plot)

                # add lepton/jet multiplicity label
            ax.text(
                0.02, 0.90, "$\\alpha_{S}$\n%s" % cfeatures.channel_labels[f"{args.lepton}_{jmult}"],
                fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
            )
            hep.cms.label(ax=ax, label="Preliminary", data=False, year=year_to_use, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))

            figname = os.path.join(pltdir, "_".join([jmult, args.lepton, "alphaS", "FinalSysTemplates_RatioOnly"]))
            fig.savefig(figname)
            print(f"{figname} written")
            plt.close(fig)

            #set_trace()
            for tname in hdict.keys():
                if ("nosys" in tname) or ("alphaS" in tname): continue
                pdf_sys = tname.split("TT_")[-1]
                print(f"{args.year}, {args.lepton}, {jmult}, {pdf_sys}")
   
                fig, ax = plt.subplots(figsize=(15.0, 10.0))
                fig.subplots_adjust(hspace=.07)

                    # get original variations
                orig_pdf = orig_hdict[jmult][lep][f"TT_{pdf_sys}"].copy()
                orig_pdf_vals, orig_pdf_bins = Plotter.get_ratio_arrays(num_vals=orig_pdf.values()[()], denom_vals=orig_nominal.values()[()], input_bins=orig_nominal.dense_axes()[0].edges())
                        # plot
                ax.fill_between(orig_pdf_bins, orig_pdf_vals, y2=1., facecolor="r", step="post", alpha=0.5)

                    # get final variations
                final_pdf, _ = hdict[f"TT_{pdf_sys}"].copy()
                final_pdf_vals, final_pdf_bins = Plotter.get_ratio_arrays(num_vals=final_pdf.values()[()], denom_vals=final_nominal.values()[()], input_bins=final_nominal.dense_axes()[0].edges())
                        # plot                
                ax.step(final_pdf_bins, final_pdf_vals, where="post", **{"color": "r", "linestyle": "-", "label": f"Final"})

                    # save sum of distance from nominal and ratio-to-nominal array
                pdf_vars_dict[np.sum(np.sqrt(np.square(final_pdf.values()[()]-final_nominal.values()[()])))] = [pdf_sys, final_pdf_vals]

                ax.legend(loc="upper right", title="TT")
                ax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
                ax.autoscale()
                ax.set_xlim(x_lims)
                ax.set_xlabel("$m_{t\\bar{t}}$ [GeV]")
                ax.set_ylabel("Ratio to Nominal")

                    # draw vertical lines separating ctstar bins
                [ax.axvline(vline, color="k", linestyle="--") for vline in vlines]

                [ax.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                    xytext=(0, 400), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0) for idx, label in enumerate(ctstar_binlabels)]

                #set_trace()
                ax.set_xticks(mtt_bin_inds_to_plot)
                ax.set_xticklabels(mtt_bins_to_plot)

                    # add lepton/jet multiplicity label
                ax.text(
                    0.02, 0.90, f"{pdf_sys}\n{cfeatures.channel_labels[f'{args.lepton}_{jmult}']}",
                    fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                )
                hep.cms.label(ax=ax, label="Preliminary", data=False, year=year_to_use, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))

                figname = os.path.join(pltdir, "_".join([jmult, args.lepton, pdf_sys, "FinalSysTemplates_RatioOnly"]))
                fig.savefig(figname)
                print(f"{figname} written")
                plt.close(fig)
                #set_trace()

            #set_trace()    
            ### make plots comparing top 5 highest impact pdf uncs + alphaS variations
            top_vars = np.sort(np.array(list(pdf_vars_dict.keys())))[::-1][:5]
            fig, ax = plt.subplots(figsize=(15.0, 10.0))
            fig.subplots_adjust(hspace=.07)
    
            ax.step(final_alphaSUp_bins, final_alphaSUp_vals, where="post", **{"linestyle": "-", "label": "$\\alpha_{S}$ Up"})
            ax.step(final_alphaSDown_bins, final_alphaSDown_vals, where="post", **{"linestyle": "-", "label": "$\\alpha_{S}$ Down"})
    
            for pdf_var in top_vars:
                pdf_sys, pdf_vals = pdf_vars_dict[pdf_var]
                ax.step(final_alphaSDown_bins, pdf_vals, where="post", **{"linestyle": "-", "label": pdf_sys})
    
            ax.legend(loc="upper right", ncol=3)
            ax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
            ax.autoscale()
            #set_trace()
            ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.01)
            ax.set_xlim(x_lims)
            ax.set_xlabel("$m_{t\\bar{t}}$ [GeV]")
            ax.set_ylabel("Ratio to Nominal")
    
                # draw vertical lines separating ctstar bins
            [ax.axvline(vline, color="k", linestyle="--") for vline in vlines]
    
            [ax.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                xytext=(0, 50), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0) for idx, label in enumerate(ctstar_binlabels)]
    
            #set_trace()
            ax.set_xticks(mtt_bin_inds_to_plot)
            ax.set_xticklabels(mtt_bins_to_plot)
    
                # add lepton/jet multiplicity label
            ax.text(
                0.02, 0.92, cfeatures.channel_labels[f'{args.lepton}_{jmult}'],
                fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
            )
            hep.cms.label(ax=ax, label="Preliminary", data=False, year=year_to_use, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))
    
            #set_trace()
            figname = os.path.join(pltdir, "_".join([jmult, args.lepton, "PDF_Comp", "FinalSysTemplates_RatioOnly"]))
            fig.savefig(figname)
            print(f"{figname} written")
            plt.close(fig)


    ## make plots for background templates
else:
    for jmult in final_hdict.keys():
        for lep, hdict in final_hdict[jmult].items():
            if lep != args.lepton: continue
            #set_trace()
            systypes = sorted(systematics.sys_groups[args.year].keys())
            for sys in systypes:
                    # find histograms of associated systematics and their processes
                up_sysname = systematics.sys_groups[args.year][sys][0]
                dw_sysname = systematics.sys_groups[args.year][sys][1]
                if (sys == "CR1") or (sys == "CR2") or (sys == "erdON"):
                    dw_sysname = up_sysname + "Down"
                    up_sysname += "Up"
                procs_sys = sorted(set([key.split(f"_{up_sysname}")[0] for key in hdict.keys() if up_sysname in key])) if not dw_sysname \
                    else sorted(set([key.split(f"_{up_sysname}")[0] for key in hdict.keys() if up_sysname in key] + [key.split(f"_{dw_sysname}")[0] for key in hdict.keys() if dw_sysname in key]))
    
                if not procs_sys: continue
    
                combine_sysname = systematics.combine_template_sys_to_name[args.year][up_sysname].replace("Up", "")
                if "LEP" in combine_sysname: combine_sysname = combine_sysname.replace("LEP", args.lepton[0].lower())
                if (args.year == "2016APV") and ("2016APV" in combine_sysname): combine_sysname = combine_sysname.replace("2016APV", "2016pre")
                if (args.year == "2016") and ("2016" in combine_sysname): combine_sysname = combine_sysname.replace("2016", "2016post")
                if "CHAN_" in combine_sysname: combine_sysname = combine_sysname.replace("CHAN_", "")
    
                for proc in procs_sys:
                    #print(lep, jmult, sys, proc)
                    orig_exists = (f"{proc}_{up_sysname}" in orig_hdict[jmult][lep].keys()) & (f"{proc}_{dw_sysname}" in orig_hdict[jmult][lep].keys())
                    final_exists = (f"{proc}_{up_sysname}" in final_hdict[jmult][lep].keys()) & (f"{proc}_{dw_sysname}" in final_hdict[jmult][lep].keys())
                    if (not orig_exists) and (not final_exists): continue
    
                    print(lep, jmult, combine_sysname, proc)
                    #set_trace()
                        # make sure output directory exists
                    pltdir = os.path.join(outdir, jmult, combine_sysname) if args.process == "bkg" else os.path.join(outdir, jmult, proc, combine_sysname)
                    if not os.path.isdir(pltdir):
                        os.makedirs(pltdir)
    
                    up_histos, dw_histos = [], []
    
                    #set_trace()
                    nominal, _ = hdict[f"{proc}_nosys"]
                    x_lims = (0, nominal.dense_axes()[0].centers().size)
    
                    if (sys == "CR1") or (sys == "CR2") or (sys == "erdON"):
                        orig_up = orig_hdict[jmult][lep][f"{proc}_{sys}"].copy()
                        up_histos.append((orig_up, {"color": "r", "linestyle": "--", "label": f"Original Up"}, True))

                    if orig_exists:
                        orig_up, orig_dw = orig_hdict[jmult][lep][f"{proc}_{up_sysname}"].copy(), orig_hdict[jmult][lep][f"{proc}_{dw_sysname}"].copy()
                        up_histos.append((orig_up, {"color": "r", "linestyle": "--", "label": f"Original Up"}, True))
                        dw_histos.append((orig_dw, {"color": "b", "linestyle": "--", "label": f"Original Down"}, True))
    
                    if final_exists:
                        final_up, treatment_up = hdict[f"{proc}_{up_sysname}"]
                        final_dw, treatment_dw = hdict[f"{proc}_{dw_sysname}"]
                        up_histos.append((final_up, {"color": "r", "linestyle": "-", "label": f"Final Up"}, False))
                        if not ((sys == "CR1") or (sys == "CR2") or (sys == "erdON")):
                            dw_histos.append((final_dw, {"color": "b", "linestyle": "-", "label": f"Final Down"}, False))
    
                    if (not nominal.values()): continue
    
                    fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True, figsize=(15.0, 10.0))
                    fig.subplots_adjust(hspace=.07)
    
                    ratio_fig, ratio_ax = plt.subplots(figsize=(15.0, 10.0))
                    ratio_fig.subplots_adjust(hspace=.07)
    
                    #set_trace()
                    hep.plot.histplot(nominal.values()[()], nominal.dense_axes()[0].edges(), ax=ax, histtype="step", **{"color":"k", "linestyle":"-", "label":"Nominal"})
    
                        ## plot relative deviations
                    for up_histo, up_style, use_fill_between in up_histos:
                        #set_trace()
                            # there is at least one actual value
                        if np.any(~np.isnan(up_histo.values()[()])):
                            up_masked_vals, up_masked_bins = Plotter.get_ratio_arrays(num_vals=up_histo.values()[()], denom_vals=nominal.values()[()], input_bins=nominal.dense_axes()[0].edges())
                            if use_fill_between:
                                rax.fill_between(up_masked_bins, up_masked_vals, y2=1., facecolor=up_style["color"], step="post", alpha=0.5)
                                ratio_ax.fill_between(up_masked_bins, up_masked_vals, y2=1., facecolor=up_style["color"], step="post", alpha=0.5)
                            else:
                                hep.plot.histplot(up_histo.values()[()], nominal.dense_axes()[0].edges(), ax=ax, histtype="step", **up_style)
                                rax.step(up_masked_bins, up_masked_vals, where="post", **up_style)
                                ratio_ax.step(up_masked_bins, up_masked_vals, where="post", **up_style)
    
                    for dw_histo, dw_style, use_fill_between in dw_histos:
                            # there is at least one actual value
                        if np.any(~np.isnan(dw_histo.values()[()])):
                            dw_masked_vals, dw_masked_bins = Plotter.get_ratio_arrays(num_vals=dw_histo.values()[()], denom_vals=nominal.values()[()], input_bins=nominal.dense_axes()[0].edges())
                            if use_fill_between:
                                rax.fill_between(dw_masked_bins, dw_masked_vals, y2=1., facecolor=dw_style["color"], step="post", alpha=0.5)
                                ratio_ax.fill_between(dw_masked_bins, dw_masked_vals, y2=1., facecolor=dw_style["color"], step="post", alpha=0.5)
                            else:
                                hep.plot.histplot(dw_histo.values()[()], nominal.dense_axes()[0].edges(), ax=ax, histtype="step", **dw_style)
                                rax.step(dw_masked_bins, dw_masked_vals, where="post", **dw_style)
                                ratio_ax.step(dw_masked_bins, dw_masked_vals, where="post", **dw_style)
    
                    ax.set_yscale("log")
                    ax.legend(loc="upper right", title=proc, ncol=2)
                    ax.autoscale()
                    ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.5)
                    ax.set_xlim(x_lims)
                    ax.set_xlabel(None)
                    ax.set_ylabel("Events")
    
                    rax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
                    rax.autoscale()
                    rax.set_xlim(x_lims)
                    rax.set_xlabel("$m_{t\\bar{t}}$ [GeV]")
                    rax.set_ylabel("Ratio to Nominal")
    
                    ratio_ax.legend(loc="upper right", title=proc, ncol=2)                
                    ratio_ax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
                    ratio_ax.autoscale()
                    ratio_ax.set_xlim(x_lims)
                    ratio_ax.set_xlabel("$m_{t\\bar{t}}$ [GeV]")
                    ratio_ax.set_ylabel("Ratio to Nominal")
                    
                        # add lepton/jet multiplicity label
                    ax.text(
                        0.02, 0.86, f"{combine_sysname}\n{cfeatures.channel_labels[f'{args.lepton}_{jmult}']}",
                        fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                    )
                    ratio_ax.text(
                        0.02, 0.90, f"{combine_sysname}\n{cfeatures.channel_labels[f'{args.lepton}_{jmult}']}",
                        fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ratio_ax.transAxes
                    )
                        ## draw vertical lines for distinguishing different ctstar bins
                    vlines = [x_lims[1]*ybin/5 for ybin in range(1, 5)]
                    [ax.axvline(vline, color="k", linestyle="--") for vline in vlines]
                    [rax.axvline(vline, color="k", linestyle="--") for vline in vlines]
                    [ratio_ax.axvline(vline, color="k", linestyle="--") for vline in vlines]
    
                    for idx, label in enumerate(ctstar_binlabels):
                        ax.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                            xytext=(0, 10), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)
                        ratio_ax.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                            xytext=(0, 10), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)
    
                    #set_trace()
                    rax.set_xticks(mtt_bin_inds_to_plot)
                    rax.set_xticklabels(mtt_bins_to_plot)
    
                    ratio_ax.set_xticks(mtt_bin_inds_to_plot)
                    ratio_ax.set_xticklabels(mtt_bins_to_plot)
    
                    hep.cms.label(ax=ax, data=False, label="Preliminary", year=year_to_use, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))
                    hep.cms.label(ax=ratio_ax, data=False, label="Preliminary", year=year_to_use, lumi=round(data_lumi_year[f"{args.lepton}s"]/1000., 1))
                    
                    figname = os.path.join(pltdir, "_".join([proc, jmult, args.lepton, combine_sysname, "FinalSysTemplates"]))
                    fig.savefig(figname)
                    print(f"{figname} written")
                    plt.close(fig)
                    
                    ratio_figname = os.path.join(pltdir, "_".join([proc, jmult, args.lepton, combine_sysname, "FinalSysTemplates_RatioOnly"]))
                    ratio_fig.savefig(ratio_figname)
                    print(f"{ratio_figname} written")
                    plt.close(ratio_fig)
                    #set_trace()

toc = time.time()
print("Total time: %.1f" % (toc - tic))
