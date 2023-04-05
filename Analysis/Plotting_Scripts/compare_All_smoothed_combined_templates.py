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
from coffea import hist
import numpy as np
import Utilities.Plotter as Plotter
import Utilities.systematics as systematics
from coffea.hist import plot
from Utilities.styles import styles as hstyles
import Utilities.final_analysis_binning as final_binning
import Utilities.common_features as cfeatures

base_jobid = os.environ["base_jobid"]

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("combination", choices=["lepton", "era_lepton", "indiv"], help="What type of combination has been performed.")
#parser.add_argument("process", choices=["bkg"], help="Specify which process to use.")
parser.add_argument("process", choices=["bkg", "MEreweight_sig", "sig"], help="Specify which process to use.")
parser.add_argument("--year", choices=["2016APV", "2016", "2017", "2018"], help="Specify which year to use if combination is only across leptons.")
parser.add_argument("--nomSMTTxsec", action="store_true", help="Apply nominal SM cross sections to top mass and LHE scale weights")
args = parser.parse_args()

if (args.combination == "lepton") and (not args.year):
    raise ValueError("Year must be specified when the combination is only across leptons")

if (args.combination == "indiv") and (not args.year):
    raise ValueError("Year must be specified when the individual channel templates are used")

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
analyzer = "htt_btag_sb_regions"
plot_outdir = os.environ["plots_dir"]
eos_dir = os.environ["eos_dir"]

outdir = os.path.join(plot_outdir, jobid if args.combination == "era_lepton" else f"{args.year}_{jobid}", f"Compare_All_Smoothed_Combined_Templates_{analyzer}", (args.process).upper(), (args.combination).upper())

year_to_use = cfeatures.year_labels[args.year] if args.year else cfeatures.year_labels["Total"]

#set_trace()
if args.process == "bkg":
    if args.combination == "lepton":
        import Utilities.prettyjson as prettyjson
        data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())[args.year]
        lumi_to_use = (data_lumi_year["Muons"]+data_lumi_year["Electrons"])/2000

        orig_hdict = load(os.path.join(eos_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}", f"raw_combined_lep_templates_lj_bkg_{args.year}_{jobid}.coffea"))
        smoothed_hdict = load(os.path.join(eos_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}", f"smoothed_combined_lep_templates_lj_bkg_{args.year}_{jobid}.coffea"))

        systypes = [
            "JES_Absolute", f"JES_Absolute_{args.year}",
            "JES_BBEC1", 
            f"JES_BBEC1_{args.year}", 
            "JES_RelativeBal", f"JES_RelativeSample_{args.year}",
            "JES_FlavorQCD",
            "JES_FlavorQCDOnlyLightJets", "JES_FlavorPureBottomOnlyBottomJets", "JES_FlavorPureBottom",
            "JES_FlavorPureCharmOnlyCharmJets", "JES_FlavorPureCharm", "JES_FlavorPureQuarkOnlyQuarkJets", "JES_FlavorPureQuark",
            "JES_FlavorPureGluonOnlyGluonJets", "JES_FlavorPureGluon", 
            "JER", "MET", "SHAPE",
            "PILEUP",
            "BTAG_BC_JES", "BTAG_BC_PILEUP", "BTAG_BC_STATISTIC", "BTAG_BC_TYPE3",
            "BTAG_BC_CORR", "BTAG_BC_UNCORR",
            "BTAG_L_CORR",
            "BTAG_L_UNCORR",
            "BTAG_BC_BFRAGMENTATION", "BTAG_BC_BTEMPCORR", "BTAG_BC_CB", "BTAG_BC_CFRAGMENTATION",
            "BTAG_BC_CJETS", "BTAG_BC_DMUX", "BTAG_BC_GLUONSPLITTING", "BTAG_BC_JETAWAY",
            "BTAG_BC_KSL", "BTAG_BC_L2C", "BTAG_BC_LTOTHERS", "BTAG_BC_MUDR", "BTAG_BC_MUPT", "BTAG_BC_PTREL"
        ]
        if args.year != "2018": systypes.append("PREFIRE")

        year_use = args.year

    if args.combination == "era_lepton":
        orig_hdict = load(os.path.join(eos_dir, "results", jobid, f"Templates_{analyzer}", f"raw_combined_year_and_lepton_templates_lj_bkg_nomSMTTxsec_{jobid}.coffea" if args.nomSMTTxsec else f"raw_combined_year_and_lepton_templates_lj_bkg_{jobid}.coffea"))
        smoothed_hdict = load(os.path.join(eos_dir, "results", jobid, f"Templates_{analyzer}", f"smoothed_combined_year_and_lepton_templates_lj_bkg_nomSMTTxsec_{jobid}.coffea" if args.nomSMTTxsec else f"smoothed_combined_year_and_lepton_templates_lj_bkg_{jobid}.coffea"))
        smoothed_1d_hdict = load(os.path.join(eos_dir, "results", jobid, f"Templates_{analyzer}", f"smoothed_1D_combined_year_and_lepton_templates_lj_bkg_{jobid}.coffea"))

        #set_trace()
        if args.nomSMTTxsec:
            outdir = os.path.join(outdir, "NominalTTbarXsec")
            systypes = ["FACTOR", "RENORM", "MTOP3GEV", "MTOP1GEV"]
            #systypes = ["EWK_scheme", "EWK_yukawa", "ISR", "FSR", "FACTOR", "RENORM", "HDAMP", "UE", "MTOP3GEV", "MTOP1GEV", "CR1", "CR2", "erdON"]
        else:
            systypes = ["EWK_scheme", "EWK_yukawa", "ISR", "FSR", "FACTOR", "RENORM", "HDAMP", "UE", "MTOP3GEV", "MTOP1GEV", "CR1", "CR2", "erdON"]
        year_use = "2017"

    if args.combination == "indiv":
        import Utilities.prettyjson as prettyjson
        data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())[args.year]
        lumi_to_use = (data_lumi_year["Muons"]+data_lumi_year["Electrons"])/2000

        orig_hdict = load(os.path.join(eos_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}", f"raw_templates_lj_bkg_{args.year}_{jobid}.coffea"))
        smoothed_hdict = load(os.path.join(eos_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}", f"smoothed_templates_lj_bkg_{args.year}_{jobid}.coffea"))

        systypes = systematics.sys_groups[args.year].keys()
        year_use = args.year

#if args.process == "sig":
#    orig_hdict = load(os.path.join(input_dir, f"raw_templates_lj_sig_{args.year}_{jobid}.coffea"))
#    comb_lep_hdict = load(os.path.join(input_dir, f"raw_combined_lep_templates_lj_sig_{args.year}_{jobid}.coffea"))
#    #comb_era_hdict = load(os.path.join(proj_dir, "results", jobid, f"Templates_{analyzer}", f"combined_year_templates_lj_sig_{jobid}.coffea"))
#    comb_era_lep_hdict = load(os.path.join(proj_dir, "results", jobid, f"Templates_{analyzer}", f"raw_combined_year_and_lepton_templates_lj_sig_{jobid}.coffea"))
#
#    systypes = ["AH_ISR", "AH_FSR", "AH_RENORM", "AH_FACTOR"]
#

#if args.process == "MEreweight_sig":
#    templates_names = {
#        "Orig" : load(os.path.join(input_dir, f"raw_templates_lj_MEsig_{args.year}_{jobid}.coffea")),
#        #"Smooth" : load(os.path.join(input_dir, f"smoothed_templates_lj_MEsig_{args.year}_{jobid}.coffea")),
#        #"Flat" : load(os.path.join(input_dir, f"flattened_templates_lj_MEsig_{args.year}_{jobid}.coffea")),
#        #"Symm" : load(os.path.join(input_dir, f"symmetrized_templates_lj_MEsig_{args.year}_{jobid}.coffea")),
#    }

if not os.path.isdir(outdir):
    os.makedirs(outdir)

linearize_binning = (
    final_binning.mtt_binning,
    final_binning.ctstar_abs_binning
)
ctstar_binlabels = [r"%s $\in [%s, %s)$" % (cfeatures.variable_names_to_labels["tlep_ctstar_abs"], linearize_binning[1][bin], linearize_binning[1][bin+1]) for bin in range(len(linearize_binning[1])-1)]
ctstar_binlabels[-1] = r"%s $\in [%s, %s]$" % (cfeatures.variable_names_to_labels["tlep_ctstar_abs"], linearize_binning[1][-2], linearize_binning[1][-1])
ctstar_bin_locs = np.linspace((len(linearize_binning[0])-1)/2, (len(linearize_binning[0])-1)*(len(linearize_binning[1])-1) - (len(linearize_binning[0])-1)/2, len(linearize_binning[1])-1)
#set_trace()
mtt_vals_to_plot = np.array([400, 600, 1000])
mtt_tiled_labels = np.tile(linearize_binning[0][:-1], linearize_binning[1].size-1)
mtt_bin_inds_to_plot = np.where(np.in1d(mtt_tiled_labels, mtt_vals_to_plot))[0]
mtt_bins_to_plot = np.tile(mtt_vals_to_plot, linearize_binning[1].size-1)

#set_trace()

if args.combination == "indiv":
    for jmult in sorted(orig_hdict.keys()):
        for lep in ["Muon", "Electron"]:
            orig_dict = orig_hdict[jmult][lep]
            smooth_dict = smoothed_hdict[jmult][lep]  
  
                # get all keys from both files to make sure they"re the same    
            orig_keys = sorted(orig_dict.keys())
    
            #set_trace()
            for sys in systypes:
                    # find histograms of associated systematics and their processes
                up_sysname = systematics.sys_groups[year_use][sys][0]
                dw_sysname = systematics.sys_groups[year_use][sys][1] if len(systematics.sys_groups[year_use][sys]) > 1 else None
                procs_sys = sorted(set([key.split(f"_{up_sysname}")[0] for key in orig_keys if up_sysname in key])) if not dw_sysname \
                    else sorted(set([key.split(f"_{up_sysname}")[0] for key in orig_keys if up_sysname in key] + [key.split(f"_{dw_sysname}")[0] for key in orig_keys if dw_sysname in key]))
    
                if not procs_sys: continue
    
                combine_sysname = systematics.combine_template_sys_to_name[year_use][up_sysname].replace("Up", "")
                if "LEP" in combine_sysname: combine_sysname = combine_sysname.replace("LEP", lep[0].lower())
                if (args.year == "2016APV") and ("2016APV" in combine_sysname): combine_sysname = combine_sysname.replace("2016APV", "2016pre")
                if (args.year == "2016") and ("2016" in combine_sysname): combine_sysname = combine_sysname.replace("2016", "2016post")
                if "CHAN_" in combine_sysname: combine_sysname = combine_sysname.replace("CHAN_", "")
    
                for proc in procs_sys:
                    # check if same hists exist in other files
                    smooth_exists = (f"{proc}_{up_sysname}" in smooth_dict.keys()) & (f"{proc}_{dw_sysname}" in smooth_dict.keys())
                    if not smooth_exists: continue
    
                    if not ( (f"{proc}_{up_sysname}" in orig_dict.keys() and f"{proc}_{dw_sysname}" in orig_dict.keys()) and f"{proc}_nosys" in orig_dict.keys() ): continue
    
                    if ((args.process == "sig") or (args.process == "MEreweight_sig")) and ("M500_W10p0" not in proc): continue
                    
                    print(jmult, lep, combine_sysname, proc)
    
                    #set_trace()
                        # make sure output directory exists
                    pltdir = os.path.join(outdir, lep, jmult, combine_sysname) if args.process == "bkg" else os.path.join(outdir, lep, jmult, proc, combine_sysname)
                    if not os.path.isdir(pltdir):
                        os.makedirs(pltdir)
    
                    up_histos, dw_histos = [], []
    
                    #set_trace()
                    orig_nominal, orig_up, orig_dw = orig_dict[f"{proc}_nosys"].copy(), orig_dict[f"{proc}_{up_sysname}"].copy(), orig_dict[f"{proc}_{dw_sysname}"].copy()
                    up_histos.append((orig_nominal, orig_up, {"color": "r", "linestyle": "-", "label": f"Original Up, {cfeatures.channel_labels[f'{lep}_{jmult}']}"}, True))
                    dw_histos.append((orig_nominal, orig_dw, {"color": "b", "linestyle": "-", "label": f"Original Down, {cfeatures.channel_labels[f'{lep}_{jmult}']}"}, True))
    
                    if smooth_exists:
                        smooth_up, smooth_dw = smooth_dict[f"{proc}_{up_sysname}"].copy(), smooth_dict[f"{proc}_{dw_sysname}"].copy()
                        up_histos.append((orig_nominal, smooth_up, {"color": "g", "linestyle": "-", "label": f"Smooth Up, {cfeatures.channel_labels[f'{lep}_{jmult}']}"}, False))
                        dw_histos.append((orig_nominal, smooth_dw, {"color": "k", "linestyle": "-", "label": f"Smooth Down, {cfeatures.channel_labels[f'{lep}_{jmult}']}"}, False))
    
                    x_lims = (0, orig_nominal.dense_axes()[0].centers().size)
            
                    fig, ax = plt.subplots(figsize=(15.0, 10.0))
                    fig.subplots_adjust(hspace=.07)
    
                        ## plot relative deviations
                    for nominal, up_histo, up_style, use_fill_between in up_histos:
                            # there is at least one actual value
                        if np.any(~np.isnan(up_histo.values()[()])):
                            up_masked_vals, up_masked_bins = Plotter.get_ratio_arrays(num_vals=up_histo.values()[()], denom_vals=nominal.values()[()] if nominal else np.ones(up_histo.values()[()].size), input_bins=up_histo.dense_axes()[0].edges())
                            ax.fill_between(up_masked_bins, up_masked_vals, y2=1., facecolor=up_style["color"], label=up_style["label"], step="post", alpha=0.5) if use_fill_between\
                                else ax.step(up_masked_bins, up_masked_vals, where="post", **up_style)
    
                    for nominal, dw_histo, dw_style, use_fill_between in dw_histos:
                            # there is at least one actual value
                        if np.any(~np.isnan(dw_histo.values()[()])):
                            dw_masked_vals, dw_masked_bins = Plotter.get_ratio_arrays(num_vals=dw_histo.values()[()], denom_vals=nominal.values()[()] if nominal else np.ones(dw_histo.values()[()].size), input_bins=dw_histo.dense_axes()[0].edges())
                            ax.fill_between(dw_masked_bins, dw_masked_vals, y2=1., facecolor=dw_style["color"], label=dw_style["label"], step="post", alpha=0.5) if use_fill_between \
                                else ax.step(dw_masked_bins, dw_masked_vals, where="post", **dw_style)
    
                    #set_trace()
                    ax.legend(loc="upper right", title=proc, ncol=2)
                    ax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
                    ax.autoscale()
                    #set_trace()
                    
                    ax.set_xlim(x_lims)
                    ax.set_xlabel("$m_{t\\bar{t}}$ [GeV]")
                    ax.set_ylabel("Ratio to Nominal")
                    ax.set_ylim(max(ax.get_ylim()[0], 0.65), min(ax.get_ylim()[1], 1.25))
                    
                        # add lepton/jet multiplicity label
                    ax.text(
                        0.02, 0.92, combine_sysname,
                        fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                    )
                        ## draw vertical lines for distinguishing different ctstar bins
                    vlines = [x_lims[1]*ybin/5 for ybin in range(1, 5)]
                    [ax.axvline(vline, color="k", linestyle="--") for vline in vlines]
    
                    for idx, label in enumerate(ctstar_binlabels):
                        ax.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                            xytext=(0, 25), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)
    
                    #set_trace()
                    ax.set_xticks(mtt_bin_inds_to_plot)
                    ax.set_xticklabels(mtt_bins_to_plot)
    
                    hep.cms.label(ax=ax, data=False, year=year_to_use, lumi=round(lumi_to_use, 1))
                    
                    #set_trace()
                    figname = os.path.join(pltdir, "_".join([proc, jmult, lep, combine_sysname, "Indiv_Smoothed_Templates_Comp"]))
                    fig.savefig(figname)
                    print(f"{figname} written")
                    plt.close(fig)


else:
    for jmult in sorted(orig_hdict.keys()):
        orig_dict = orig_hdict[jmult]
    
            # get all keys from both files to make sure they"re the same    
        orig_keys = sorted(orig_dict.keys())
    
        #set_trace()
        for sys in systypes:
                # find histograms of associated systematics and their processes
            up_sysname = systematics.sys_groups[year_use][sys][0]
            dw_sysname = systematics.sys_groups[year_use][sys][1]
            procs_sys = sorted(set([key.split(f"_{up_sysname}")[0] for key in orig_keys if up_sysname in key])) if not dw_sysname \
                else sorted(set([key.split(f"_{up_sysname}")[0] for key in orig_keys if up_sysname in key] + [key.split(f"_{dw_sysname}")[0] for key in orig_keys if dw_sysname in key]))
    
            if not procs_sys: continue

            #set_trace()    
            combine_sysname = systematics.combine_template_sys_to_name[year_use][up_sysname+"Up" if ((sys == "CR1") or (sys == "CR2") or (sys == "erdON")) else up_sysname].replace("Up", "")
            if "LEP" in combine_sysname: combine_sysname = combine_sysname.replace("LEP", args.lepton[0].lower())
            if (args.year == "2016APV") and ("2016APV" in combine_sysname): combine_sysname = combine_sysname.replace("2016APV", "2016pre")
            if (args.year == "2016") and ("2016" in combine_sysname): combine_sysname = combine_sysname.replace("2016", "2016post")
            if "CHAN_" in combine_sysname: combine_sysname = combine_sysname.replace("CHAN_", "")
    
            for proc in procs_sys:
                # check if same hists exist in other files
                if (up_sysname is not None) & (dw_sysname is not None):
                    smooth_exists = (f"{proc}_{up_sysname}" in smoothed_hdict[jmult].keys()) & (f"{proc}_{dw_sysname}" in smoothed_hdict[jmult].keys())
                    orig_exists = (f"{proc}_{up_sysname}" in orig_dict.keys() and f"{proc}_{dw_sysname}" in orig_dict.keys()) and f"{proc}_nosys" in orig_dict.keys()
                elif (up_sysname is not None) & (dw_sysname is None):
                    smooth_exists = (f"{proc}_{up_sysname}" in smoothed_hdict[jmult].keys())
                    orig_exists = f"{proc}_{up_sysname}" in orig_dict.keys() and f"{proc}_nosys" in orig_dict.keys()
                if not smooth_exists: continue
                if not orig_exists: continue
    
                if ((args.process == "sig") or (args.process == "MEreweight_sig")) and ("M500_W10p0" not in proc): continue
                
                print(jmult, combine_sysname, proc)
    
                #set_trace()
                    # make sure output directory exists
                pltdir = os.path.join(outdir, jmult, combine_sysname) if args.process == "bkg" else os.path.join(outdir, jmult, proc, combine_sysname)
                if not os.path.isdir(pltdir):
                    os.makedirs(pltdir)
    
                up_histos, dw_histos = [], []
    
                orig_nominal, orig_up = orig_dict[f"{proc}_nosys"].copy(), orig_dict[f"{proc}_{up_sysname}"].copy()
                up_histos.append((orig_nominal, orig_up, {"color": "r", "linestyle": "-", "label": f"Original Up, {cfeatures.channel_labels[f'Lepton_{jmult}']}"}, True))
                if dw_sysname is not None:
                    orig_dw = orig_dict[f"{proc}_{dw_sysname}"].copy()
                    dw_histos.append((orig_nominal, orig_dw, {"color": "b", "linestyle": "-", "label": f"Original Down, {cfeatures.channel_labels[f'Lepton_{jmult}']}"}, True))
    
                if smooth_exists:
                    #set_trace()
                    smooth_up = smoothed_hdict[jmult][f"{proc}_{up_sysname}"].copy()
                    up_histos.append((None, smooth_up, {"color": "g", "linestyle": "-", "label": f"Smooth Up, {cfeatures.channel_labels[f'Lepton_{jmult}']}"}, False))
                        # 1d hists
                    smooth_1d_up = smoothed_1d_hdict[jmult][f"{proc}_{up_sysname}"].copy()
                    up_histos.append((None, smooth_1d_up, {"color": "r", "linestyle": "-", "label": f"Smooth 1D Up, {cfeatures.channel_labels[f'Lepton_{jmult}']}"}, False))
                    if dw_sysname is not None:
                        smooth_dw = smoothed_hdict[jmult][f"{proc}_{dw_sysname}"].copy()
                        dw_histos.append((None, smooth_dw, {"color": "k", "linestyle": "-", "label": f"Smooth Down, {cfeatures.channel_labels[f'Lepton_{jmult}']}"}, False))
                            # 1d hists
                        smooth_1d_dw = smoothed_1d_hdict[jmult][f"{proc}_{dw_sysname}"].copy()
                        dw_histos.append((None, smooth_1d_dw, {"color": "b", "linestyle": "-", "label": f"Smooth 1D Down, {cfeatures.channel_labels[f'Lepton_{jmult}']}"}, False))
    
                x_lims = (0, orig_nominal.dense_axes()[0].centers().size)
        
                fig, ax = plt.subplots(figsize=(15.0, 10.0))
                fig.subplots_adjust(hspace=.07)
    
                    ## plot relative deviations
                for nominal, up_histo, up_style, use_fill_between in up_histos:
                        # there is at least one actual value
                    if np.any(~np.isnan(up_histo.values()[()])):
                        up_masked_vals, up_masked_bins = Plotter.get_ratio_arrays(num_vals=up_histo.values()[()], denom_vals=nominal.values()[()] if nominal else np.ones(up_histo.values()[()].size), input_bins=up_histo.dense_axes()[0].edges())
                        ax.fill_between(up_masked_bins, up_masked_vals, y2=1., facecolor=up_style["color"], label=up_style["label"], step="post", alpha=0.5) if use_fill_between\
                            else ax.step(up_masked_bins, up_masked_vals, where="post", **up_style)
    
                for nominal, dw_histo, dw_style, use_fill_between in dw_histos:
                        # there is at least one actual value
                    if np.any(~np.isnan(dw_histo.values()[()])):
                        dw_masked_vals, dw_masked_bins = Plotter.get_ratio_arrays(num_vals=dw_histo.values()[()], denom_vals=nominal.values()[()] if nominal else np.ones(dw_histo.values()[()].size), input_bins=dw_histo.dense_axes()[0].edges())
                        ax.fill_between(dw_masked_bins, dw_masked_vals, y2=1., facecolor=dw_style["color"], label=dw_style["label"], step="post", alpha=0.5) if use_fill_between \
                            else ax.step(dw_masked_bins, dw_masked_vals, where="post", **dw_style)
    
                #set_trace()
                ax.legend(loc="upper right", title=proc, ncol=2)
                ax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
                ax.autoscale()
                #set_trace()
                
                ax.set_xlim(x_lims)
                ax.set_xlabel("$m_{t\\bar{t}}$ [GeV]")
                ax.set_ylabel("Ratio to Nominal")
                ax.set_ylim(max(ax.get_ylim()[0], 0.65), min(ax.get_ylim()[1], 1.25))
                
                    # add lepton/jet multiplicity label
                ax.text(
                    0.02, 0.92, combine_sysname,
                    fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                )
                    ## draw vertical lines for distinguishing different ctstar bins
                vlines = [x_lims[1]*ybin/5 for ybin in range(1, 5)]
                [ax.axvline(vline, color="k", linestyle="--") for vline in vlines]
    
                for idx, label in enumerate(ctstar_binlabels):
                    ax.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                        xytext=(0, 25), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)
    
                #set_trace()
                ax.set_xticks(mtt_bin_inds_to_plot)
                ax.set_xticklabels(mtt_bins_to_plot)
    
                hep.cms.label(ax=ax, data=False, year=year_to_use) if args.combination == "era_lepton" else hep.cms.label(ax=ax, data=False, year=year_to_use, lumi=round(lumi_to_use, 1))
                
                #set_trace()
                figname = os.path.join(pltdir, "_".join([proc, jmult, combine_sysname, "Combined_Smoothed_Templates_Comp"]))
                fig.savefig(figname)
                print(f"{figname} written")
                plt.close(fig)

toc = time.time()
print("Total time: %.1f" % (toc - tic))
