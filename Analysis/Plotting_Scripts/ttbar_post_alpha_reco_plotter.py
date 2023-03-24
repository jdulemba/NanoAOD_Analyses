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
import Utilities.styles as styles
import Utilities.plot_tools as plt_tools
import Utilities.prettyjson as prettyjson
import re
from coffea import hist
import numpy as np
import fnmatch
import Utilities.Plotter as Plotter
from coffea.hist import plot
import Utilities.common_features as cfeatures

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
base_jobid = os.environ["base_jobid"]
analyzer = "ttbar_post_alpha_reco"
plot_outdir = os.environ["plots_dir"]
eos_dir = os.environ["eos_dir"]

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("--plot", default="all", choices=["nosys", "uncs", "all"], help="Make plots for no systematics, variations of JES/JER systematics, or both.")
args = parser.parse_args()

input_dir = os.path.join(eos_dir, "results", f"{args.year}_{jobid}", analyzer)
f_ext = "TOT.coffea"
fnames = sorted(["%s/%s" % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)])
hdict = plt_tools.add_coffea_files(fnames) if len(fnames) > 1 else load(fnames[0])

outdir = os.path.join(plot_outdir, f"{args.year}_{jobid}", analyzer)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

alpha_corrections = {
    "E_All_1D" : ("$\\alpha_{E}$ All Lin.", "#377eb8", 2), ## blue
    "E_All_2D" : ("$\\alpha_{E}$ All Quad.", "#e42a2c", 2), ## red
    "E_Mtt" : ("$\\alpha_{E}$ Mtt", "#4daf4a", 2), ## green
    "P_All_1D" : ("$\\alpha_{P}$ All Lin.", "#984ea3", 2), ## purple
    "P_All_2D" : ("$\\alpha_{P}$ All Quad.", "#ff7f00", 2), ## orange
    "P_Mtt" : ("$\\alpha_{P}$ Mtt", "#a65628", 2), ## brown
    "Uncorrected" : ("Uncorrected", "k", 2),
}

corr_to_use = "E_All_2D"
comp_3j_mask = re.compile(r"((?:%s))" % "|".join(["E_All_2D*", "Uncorrected*"]))
alpha_corrections_mask = {
    "E_All_2D" : ("$\\alpha_{E}$", "#e42a2c", 2), ## red
    "Uncorrected" : ("Uncorrected", "k", 2),
}


systematics = {
    "nosys" : ("Nominal", "k", 2),
    "JES_Total_UP" : ("JES Total Up", "#e42a2c", 2), ## red
    "JES_Total_DW" : ("JES Total Down", "#377eb8", 2), ## blue
    "JER_UP" : ("JER Up", "#4daf4a", 2), ## green
    "JER_DW" : ("JER Down", "#984ea3", 2), ## purple
}

variables = {
    ###"Reco_mtt": ("$m_{t\\bar{t}}$ [GeV]", 2, (200., 2000.)),
    "Reco_mtt": (cfeatures.variable_names_to_labels["mtt"], 2, (200., 1000.)),
    "Reco_mthad": (cfeatures.variable_names_to_labels["mthad"], 2, (0., 250.)),
    "Reco_thad_ctstar": (cfeatures.variable_names_to_labels["thad_ctstar"], 2, (-1., 1.)),
    "Reco_thad_ctstar_abs": (cfeatures.variable_names_to_labels["thad_ctstar_abs"], 2, (0., 1.)),
    "Reso_mtt": ("$m_{t\\bar{t}}$ Resolution [GeV]", 5, (-300., 300.)),
    "Reso_mthad": ("$m_{t_{h}}$ Resolution [GeV]", 2, (-50., 150.)),
    "Reso_thad_ctstar": ("$\cos(\\theta^{*}_{t_{h}})$ Resolution", 5, (-1., 1.)),
    "Reso_thad_ctstar_abs": ("|$\cos(\\theta^{*}_{t_{h}})$| Resolution", 5, (-1., 1.)),
    ##"Reso_mtt": ("$m_{t\\bar{t}}$ Resolution [GeV]", 1, (-300., 300.)),
    ##"Reso_mthad": ("$m_{t_{h}}$ Resolution [GeV]", 2, (-200., 200.)),
    ##"Reso_thad_ctstar": ("cos($\\theta^{*}_{t_{h}}$) Resolution", 2, (-1., 1.)),
    ##"Reso_thad_ctstar_abs": ("|cos($\\theta^{*}_{t_{h}}$)| Resolution", 2, (-1., 1.)),
}


    ## get plotting colors/settings
hstyles = styles.styles

    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())[args.year]
lumi_correction = load(os.path.join(proj_dir, "Corrections", base_jobid, "MC_LumiWeights.coffea"))[args.year]
lumi_to_use = (data_lumi_year["Muons"]+data_lumi_year["Electrons"])/2000.

        # scale ttJets events, split by reconstruction type, by normal ttJets lumi correction
ttJets_permcats = ["*right", "*matchable", "*unmatchable", "*sl_tau", "*other"]
names = [dataset for dataset in sorted(set([key[0] for key in hdict[sorted(variables.keys())[0]].values().keys()]))] # get dataset names in hists

ttJets_cats = [name for name in names if any([fnmatch.fnmatch(name, cat) for cat in ttJets_permcats])] # gets ttJets(_PS)_other, ...
if len(ttJets_cats) > 0:
    for tt_cat in ttJets_cats:
        ttJets_lumi_topo = "_".join(tt_cat.split("_")[:-2]) if "sl_tau" in tt_cat else "_".join(tt_cat.split("_")[:-1]) # gets ttJets[SL, Had, DiLep] or ttJets_PS
        mu_lumi = lumi_correction["Muons"][ttJets_lumi_topo]
        el_lumi = lumi_correction["Electrons"][ttJets_lumi_topo]
        lumi_correction["Muons"].update({tt_cat: mu_lumi})
        lumi_correction["Electrons"].update({tt_cat: el_lumi})


## make groups based on process
process = hist.Cat("process", "Process", sorting="placement")
process_cat = "dataset"
process_groups = plt_tools.make_dataset_groups("Muon", args.year, samples=names)
    

    ## make plots
for hname in variables.keys():
    if hname not in hdict.keys():
        raise ValueError(f"{hname} not found in file")

    histo = hdict[hname]
    h_mu = histo[:, :, :, "Muon", :].integrate("leptype")
    h_mu.scale(lumi_correction["Muons"], axis="dataset")
    h_el = histo[:, :, :, "Electron", :].integrate("leptype")
    h_el.scale(lumi_correction["Electrons"], axis="dataset")
    h_tot = h_mu+h_el
    h_tot = h_tot.group(process_cat, process, process_groups)
    #set_trace()    

    if (args.plot == "nosys") or (args.plot == "all"):
        xtitle, rebinning, x_lims = variables[hname]
        if rebinning != 1:
            xaxis_name = h_tot.dense_axes()[0].name
            h_tot = h_tot.rebin(xaxis_name, rebinning)
    
        nosys_histo = h_tot[:, "nosys", :, :].integrate("sys")
            # make plot for each jet multiplicity
        #for jmult in ["3Jets"]:
        for jmult in sorted(set([key[1] for key in nosys_histo.values().keys()])):
            for cat in sorted(set([key[0] for key in nosys_histo.values().keys()])):
                pltdir = os.path.join(outdir, jmult, cat)
                if not os.path.isdir(pltdir):
                    os.makedirs(pltdir)

                hslice = nosys_histo[cat, jmult].integrate("process").integrate("jmult")
    
                fig, ax = plt.subplots()
                fig.subplots_adjust(hspace=.07)
    
                plot.plot1d(hslice, overlay=hslice.axes()[0].name,
                    ax=ax, clear=False, line_opts={"linestyle" : "-"},
                )
                ax.autoscale()
                ax.set_ylim(0, ax.get_ylim()[1]*1.15)
                ax.set_xlabel(xtitle)
                ax.set_xlim(x_lims)
            
                #set_trace()
                    ## set legend and corresponding colors
                handles, labels = ax.get_legend_handles_labels()
                for idx, label in enumerate(labels):
                    labels[idx] = alpha_corrections[label][0]
                    handles[idx].set_color(alpha_corrections[label][1])
                    handles[idx].set_linewidth(alpha_corrections[label][2])
                # call ax.legend() with the new values
                ax.legend(handles,labels, loc="upper right", ncol=2)
    
                    # add perm category 
                #set_trace()
                ax.text(
                    0.02, 0.90, f"{cfeatures.channel_labels[f'Lepton_{jmult}']}\n{hstyles[cat]['name']}",
                    horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                )
                hep.cms.label(ax=ax, data=False, label="Preliminary", year=cfeatures.year_labels[args.year], lumi=round(lumi_to_use, 1))
    
                figname = os.path.join(pltdir, "_".join([args.year, jobid, jmult, cat, hname]))
                fig.savefig(figname)
                print(f"{figname} written")
                plt.close()
    
            # compare plots for each jet multiplicity 
        for cat in sorted(set([key[0] for key in nosys_histo.values().keys()])):
            pltdir = os.path.join(outdir, "Comp_JMults", cat)
            if not os.path.isdir(pltdir):
                os.makedirs(pltdir)
   
                # yields 
            fig, ax = plt.subplots()
            fig.subplots_adjust(hspace=.07)
    
                # norm 
            fig_norm, ax_norm = plt.subplots()
            fig_norm.subplots_adjust(hspace=.07)
    
            for jmult in sorted(set([key[1] for key in nosys_histo.values().keys()])):
                hslice = nosys_histo[cat, jmult, :].integrate("process").integrate("jmult") if jmult == "3Jets" else nosys_histo[cat, jmult, :].integrate("process").integrate("corrtype")
                if jmult == "3Jets":
                    hslice = hslice[comp_3j_mask]

                    # yields 
                plot.plot1d(hslice, overlay=hslice.axes()[0].name,
                    ax=ax, clear=False, line_opts={"linestyle" : "-"},
                )

                    # norm
                #set_trace()
                for corr in sorted(hslice.values().keys()):
                    Plotter.plot_1D(values=hslice.values()[corr]/np.sum(hslice.values()[corr]), bins=hslice.dense_axes()[0].edges(),
                        ax=ax_norm, xlimits=x_lims, xlabel=xtitle, ylabel="Probability Density", label=corr[0], histtype="step")

                ## set legend and corresponding colors
            handles, labels = ax.get_legend_handles_labels()
            for idx, label in enumerate(labels):
                if label == "4Jets":
                    labels[idx] = cfeatures.channel_labels[f"Lepton_{label}"]
                    handles[idx].set_color("b")
                    handles[idx].set_linewidth(2)
                elif label == "5PJets":
                    labels[idx] = cfeatures.channel_labels[f"Lepton_{label}"]
                    handles[idx].set_color("g")
                    handles[idx].set_linewidth(2)
                else:
                    labels[idx] = f"{cfeatures.channel_labels['Lepton_3Jets']}, {alpha_corrections_mask[label][0]}"
                    handles[idx].set_color(alpha_corrections_mask[label][1])
                    handles[idx].set_linewidth(alpha_corrections_mask[label][2])

            # set axes and call ax.legend() with the new values
            ax.autoscale()
            ax.set_ylim(0, ax.get_ylim()[1]*1.15)
            ax.set_xlabel(xtitle)
            ax.set_xlim(x_lims)
                
            ax.legend(handles,labels, loc="upper right", ncol=2)
    
                # add perm category 
            ax.text(
                0.02, 0.93, hstyles[cat]["name"],
                horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
            )
            hep.cms.label(ax=ax, data=False, label="Preliminary", year=cfeatures.year_labels[args.year], lumi=round(lumi_to_use, 1))

            figname = os.path.join(pltdir, "_".join([args.year, jobid, "Comp_JMults", cat, hname]))
            fig.savefig(figname)
            print(f"{figname} written")
            plt.close(fig)
    
                # for norm hists
            handles, labels = ax_norm.get_legend_handles_labels()
            for idx, label in enumerate(labels):
                if label == "4Jets":
                    labels[idx] = cfeatures.channel_labels[f"Lepton_{label}"]
                    handles[idx].set_color("b")
                    handles[idx].set_linewidth(2)
                elif label == "5PJets":
                    labels[idx] = cfeatures.channel_labels[f"Lepton_{label}"]
                    handles[idx].set_color("g")
                    handles[idx].set_linewidth(2)
                else:
                    labels[idx] = f"{cfeatures.channel_labels['Lepton_3Jets']}, {alpha_corrections_mask[label][0]}"
                    handles[idx].set_color(alpha_corrections_mask[label][1])
                    handles[idx].set_linewidth(alpha_corrections_mask[label][2])

            # set axes and call ax.legend() with the new values
            ax_norm.autoscale()
            ax_norm.set_ylim(0, ax_norm.get_ylim()[1]*1.15)
            ax_norm.set_xlabel(xtitle)
            ax_norm.set_ylabel("Probability Density")
            ax_norm.set_xlim(x_lims)
                
            ax_norm.legend(handles,labels, loc="upper right", ncol=2)
    
                # add perm category 
            ax_norm.text(
                0.02, 0.93, hstyles[cat]["name"],
                horizontalalignment="left", verticalalignment="bottom", transform=ax_norm.transAxes
            )
            hep.cms.label(ax=ax_norm, data=False, label="Preliminary", year=cfeatures.year_labels[args.year], lumi=round(lumi_to_use, 1))
            #set_trace()
    
            figname_norm = os.path.join(pltdir, "_".join([args.year, jobid, "Comp_JMults", cat, hname, "Norm"]))
            fig_norm.savefig(figname_norm)
            print(f"{figname_norm} written")
            plt.close(fig_norm)
    
    
    if (args.plot == "uncs") or (args.plot == "all"):
        xtitle, rebinning, x_lims = variables[hname]
        if rebinning != 1:
            xaxis_name = h_tot.dense_axes()[0].name
            h_tot = h_tot.rebin(xaxis_name, rebinning)
   
        uncs_histo = h_tot[:, :, "3Jets", :].integrate("jmult")
            # make plot for each category
        for cat in sorted(set([key[0] for key in uncs_histo.values().keys()])):
            pltdir = os.path.join(outdir, "3Jets", "Sys_Vars", cat)
            if not os.path.isdir(pltdir):
                os.makedirs(pltdir)
            for corr in sorted(set([key[2] for key in uncs_histo.values().keys()])):

                hslice = uncs_histo[cat, :, corr].integrate("process").integrate("corrtype")

                fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                fig.subplots_adjust(hspace=.07)
    
                    # plot yields
                plot.plot1d(hslice, overlay=hslice.axes()[0].name,
                    ax=ax, clear=False, line_opts={"linestyle" : "-"},
                )
                ax.autoscale()
                ax.set_ylim(0, ax.get_ylim()[1]*1.15)
                ax.set_xlabel(None)
                ax.set_xlim(x_lims)

                    # plot ratios
                for sys in sorted([key[0] for key in hslice.values().keys()]):
                    if sys == "nosys": continue
                    nom_histo = hslice["nosys"].integrate("sys")
                    sys_histo = hslice[sys].integrate("sys")

                    ratio_masked_vals, ratio_bins = Plotter.get_ratio_arrays(num_vals=sys_histo.values()[()], denom_vals=nom_histo.values()[()], input_bins=nom_histo.dense_axes()[0].edges())
                    rax.step(ratio_bins, ratio_masked_vals, where="post", **{"linestyle" : "-", "color" : systematics[sys][1], "linewidth" :  systematics[sys][2]})

                rax.set_xlabel(xtitle)
                rax.set_ylabel("Sys/Nominal")
                rax.set_xlim(x_lims)
                rax.set_ylim(0.9, 1.1)
                #rax.set_ylim(0.8, 1.2)
                #rax.set_ylim(0.5, 1.5)
                rax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
            
                    ## set legend and corresponding colors
                handles, labels = ax.get_legend_handles_labels()
                for idx, label in enumerate(labels):
                    labels[idx] = systematics[label][0]
                    handles[idx].set_color(systematics[label][1])
                    handles[idx].set_linewidth(systematics[label][2])
                # call ax.legend() with the new values
                ax.legend(handles,labels, loc="upper right", title=alpha_corrections[corr][0], ncol=2)
    
                    # add lepton/jet mult, and tt perm category 
                ax.text(
                    0.02, 0.86, f"{cfeatures.channel_labels['Lepton_3Jets']}\n{hstyles[cat]['name']}",
                    fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                )
                hep.cms.label(ax=ax, data=False, label="Preliminary", year=cfeatures.year_labels[args.year], lumi=round(lumi_to_use, 1))
    
                #set_trace()
                figname = os.path.join(pltdir, "_".join([args.year, jobid, "3Jets", cat, corr, hname, "Sys_Comp"]))
                fig.savefig(figname)
                print(f"{figname} written")
                plt.close(fig)
