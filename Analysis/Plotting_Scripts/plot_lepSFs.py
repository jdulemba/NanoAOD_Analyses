#!/usr/bin/env python

# matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend("agg")
from matplotlib import rcParams
rcParams["font.size"] = 20
rcParams["savefig.format"] = "pdf"
rcParams["savefig.bbox"] = "tight"
color_cycler = ["k", "r", "b", "#008000", "#984ea3", "#ff7f00", "y"] # dark green, purple, orange

import coffea
from coffea.util import load
from pdb import set_trace
import os
import Utilities.prettyjson as prettyjson
import numpy as np

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("lepton", choices=["Electrons", "Muons", "Both"], help="Choose which lepton to make plots for")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
base_jobid = os.environ["base_jobid"]
analyzer = "Plot_LepSFs"
plot_outdir = os.environ["plots_dir"]

leptons_to_run = ["Muons", "Electrons"] if args.lepton == "Both" else [args.lepton]

titles = {
    "Muons" : {
        "ID" : "Muon ID: TightID/TrackerMuons",
        "ISO" : "Muon ISO: TightRelIso/TightIDandIPCut",
        "TRIG" : "Muon Trigger: IsoMuX/TightIDandTightPFIso",
        "RECO" : "Muon Reconstruction: TrackerMuons/genTracks",
    },
    "Electrons" : {
        "ID" : "Electron ID: Tight cut-based ID Fall17 V2",
        "TRIG": "Electron Trigger",
        "RECO" : "Electron Reconstruction",
        #"RECO" : "Electron Reconstruction: GsfTrack and supercluster $H/E<0.5$",
    }
}


era_styles = {
    "2016APV" : {"color" : "#e41a1c", "label" : "2016preVFP"}, # red
    "2016"    : {"color" : "#377eb8", "label" : "2016postVFP"}, # blue
    "2017"    : {"color" : "#4daf4a", "label" : "2017"}, # green
    "2018"    : {"color" : "#ff7f00", "label" : "2018"}, # orange
}

#set_trace()
cfg_pars = prettyjson.loads(open(os.path.join(proj_dir, "cfg_files", f"cfg_pars_{jobid}.json")).read())
lepSFs = load(os.path.join(proj_dir, "Corrections", base_jobid, cfg_pars["corrections"]["lepton"]))

outdir = os.path.join(plot_outdir, jobid, analyzer)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

for lep in leptons_to_run:
    for sf_type, sf_title in titles[lep].items():
        print(lep, sf_type)

            # get central lookup tables for all eras
        cen_lookups = {era : lepSFs[era][lep][sf_type]["Central"] for era in ["2016APV", "2016", "2017", "2018"]}
        errTOT_lookups = {era : lepSFs[era][lep][sf_type]["Error_tot"] for era in ["2016APV", "2016", "2017", "2018"]}
        #errSTAT_lookups = {era : lepSFs[era][lep][sf_type]["Error_stat"] for era in ["2016APV", "2016", "2017", "2018"]} if "Error_stat" in lepSFs["2018"][lep][sf_type].keys() else None
        #errSYST_lookups = {era : lepSFs[era][lep][sf_type]["Error_syst"] for era in ["2016APV", "2016", "2017", "2018"]} if "Error_syst" in lepSFs["2018"][lep][sf_type].keys() else None

        if lep == "Muons":
            #set_trace()
            if lepSFs["2018"][lep][sf_type]["Central"]._dimension == 1:
                #set_trace()
                eta_binedges = lepSFs["2018"][lep][sf_type]["Central"]._axes
                eta_bincenters = np.array([(eta_binedges[idx]+eta_binedges[idx+1])/2 for idx in range(len(eta_binedges)-1)])
                eta_xerrs = np.array([eta_bincenters[idx]-eta_binedges[idx] for idx in range(len(eta_bincenters))])

                cen_etavals = {era : cen_lookups[era]._values for era in cen_lookups.keys()}
                errTOT_etavals = {era : errTOT_lookups[era]._values for era in errTOT_lookups.keys()}

                fig, ax = plt.subplots()
                fig.subplots_adjust(hspace=.07)

                [ax.errorbar(eta_bincenters, values, xerr=eta_xerrs, yerr=errTOT_etavals[year], fmt="none", capsize=4, **era_styles[year]) for year, values in cen_etavals.items()]

                ax.set_ylabel("Scale Factor")
                ax.set_xlabel("$\\left| \\eta \\right|$" if lepSFs["2018"][lep][sf_type]["isAbsEta"] else "$\\eta$")
                ax.legend(loc="upper right")
                ax.autoscale()
                ax.set_xlim(0., eta_binedges[-1])
                ax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
                ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.01)

                    # add title
                plt.title(sf_title, loc="center", fontsize=22)

                #set_trace()                
                figname = os.path.join(outdir, f"{lep}SFs_{sf_type}_AllEta")
                fig.savefig(figname)
                print(f"{figname} written")
                plt.close(fig)

            elif lepSFs["2018"][lep][sf_type]["Central"]._dimension == 2:
                for era in lepSFs.keys():
                        # get central lookup tables for all eras
                    cen_lookup = lepSFs[era][lep][sf_type]["Central"]
                    errTOT_lookup = lepSFs[era][lep][sf_type]["Error_tot"]

                    eta_binedges = cen_lookup._axes[0]
                    eta_bincenters = np.array([(eta_binedges[idx]+eta_binedges[idx+1])/2 for idx in range(len(eta_binedges)-1)])
                    eta_xerrs = np.array([eta_bincenters[idx]-eta_binedges[idx] for idx in range(len(eta_bincenters))])
    
                    #set_trace()
                    fig, ax = plt.subplots()
                    fig.subplots_adjust(hspace=.07)
   
                    col_idx = 0
                    for pt_idx in range(cen_lookup._axes[1].size - 1):
                        ptmin, ptmax = cen_lookup._axes[1][pt_idx], cen_lookup._axes[1][pt_idx+1]
                        if ptmax <= 30.: continue
                        pt_cen_vals = cen_lookup._values[:, pt_idx]
                        pt_err_vals = errTOT_lookup._values[:, pt_idx]
                        ax.errorbar(eta_bincenters, pt_cen_vals, xerr=eta_xerrs, yerr=pt_err_vals, fmt="none", color=color_cycler[col_idx], capsize=4, label=f"${ptmin} \\leq p_T \\leq {ptmax}$ GeV")
                        col_idx += 1

                    ax.set_ylabel("Scale Factor")
                    ax.set_xlabel("$\\left| \\eta \\right|$" if lepSFs[era][lep][sf_type]["isAbsEta"] else "$\\eta$")
                    ax.legend(loc="upper right", title=era_styles[era]["label"])
                    ax.autoscale()
                    ax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
                    ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.01)
                    ax.set_xlim(eta_binedges[0], eta_binedges[-1])
    
                        # add title
                    plt.title(sf_title, loc="center", fontsize=22)
    
                    figname = os.path.join(outdir, f"{lep}SFs_{sf_type}_{era_styles[era]['label']}")
                    fig.savefig(figname)
                    print(f"{figname} written")
                    plt.close(fig)
                    #set_trace()
            else:
                raise ValueError("Only 1D and 2D scale factors are supported.")

        if lep == "Electrons":
            for era in lepSFs.keys():
                    # get central lookup tables for all eras
                cen_lookup = lepSFs[era][lep][sf_type]["Central"]
                errTOT_lookup = lepSFs[era][lep][sf_type]["Error_tot"]

                if isinstance(cen_lookup, dict):
                    #set_trace()
                    eta_inds = np.arange(len(lepSFs[era][lep][sf_type]["eta_ranges"]))
                    for eta_idx in range(int(len(lepSFs[era][lep][sf_type]["eta_ranges"])/2)):
                        neg_eta_idx, pos_eta_idx = eta_inds[eta_idx], eta_inds[-(eta_idx+1)]
                        if not np.array_equal(cen_lookup[f"eta_bin{neg_eta_idx}"]._axes, cen_lookup[f"eta_bin{pos_eta_idx}"]._axes): raise ValueError(f"Pt bins are not the same for {neg_eta_idx} and {pos_eta_idx}")
                        pt_binedges = cen_lookup[f"eta_bin{neg_eta_idx}"]._axes
                        pt_bincenters = np.array([(pt_binedges[idx]+pt_binedges[idx+1])/2 for idx in range(len(pt_binedges)-1)])
                        pt_xerrs = np.array([pt_bincenters[idx]-pt_binedges[idx] for idx in range(len(pt_bincenters))])

                        neg_etamin, neg_etamax = lepSFs[era][lep][sf_type]["eta_ranges"][neg_eta_idx]
                        neg_eta_cen_vals = cen_lookup[f"eta_bin{neg_eta_idx}"]._values
                        neg_eta_err_vals = errTOT_lookup[f"eta_bin{neg_eta_idx}"]._values

                        pos_etamin, pos_etamax = lepSFs[era][lep][sf_type]["eta_ranges"][pos_eta_idx]
                        pos_eta_cen_vals = cen_lookup[f"eta_bin{pos_eta_idx}"]._values
                        pos_eta_err_vals = errTOT_lookup[f"eta_bin{pos_eta_idx}"]._values
                        
                        fig, ax = plt.subplots()
                        fig.subplots_adjust(hspace=.07)

                        #if era == "2018": set_trace()
                        ax.errorbar(pt_bincenters[np.where(pt_binedges >= 30.)[0][0]:], neg_eta_cen_vals[np.where(pt_binedges >= 30.)[0][0]:], xerr=pt_xerrs[np.where(pt_binedges >= 30.)[0][0]:], yerr=neg_eta_err_vals[np.where(pt_binedges >= 30.)[0][0]:], fmt="none", capsize=4, color="r", label=f"$\\eta_{{SC}}$: [{neg_etamin}, {neg_etamax}]")
                        ax.errorbar(pt_bincenters[np.where(pt_binedges >= 30.)[0][0]:], pos_eta_cen_vals[np.where(pt_binedges >= 30.)[0][0]:], xerr=pt_xerrs[np.where(pt_binedges >= 30.)[0][0]:], yerr=pos_eta_err_vals[np.where(pt_binedges >= 30.)[0][0]:], fmt="none", capsize=4, color="b", label=f"$\\eta_{{SC}}$: [{pos_etamin}, {pos_etamax}]")
                  
                        ax.set_ylabel("Scale Factor")
                        ax.set_xlabel("$p_{T}$ [GeV]")
                        ax.legend(loc="upper right", title=era_styles[era]["label"])
                        ax.autoscale()
                        ax.set_xlim(pt_binedges[np.where(pt_binedges >= 30.)[0][0]], pt_binedges[-1])
                        if era != "2018": ax.set_xticks(np.r_[ax.get_xlim()[0], ax.get_xticks()[1:]]) # make sure first bin value is shown
    
                            # add title
                        plt.title(sf_title, loc="center", fontsize=22)
    
                        figname = os.path.join(outdir, f"{era_styles[era]['label']}_{lep}SFs_{sf_type}_etabin{eta_idx}")
                        fig.savefig(figname)
                        print(f"{figname} written")
                        plt.close(fig)

                elif isinstance(cen_lookup, coffea.lookup_tools.dense_lookup.dense_lookup):
                    eta_binedges = cen_lookup._axes[0]
                    eta_bincenters = np.array([(eta_binedges[idx]+eta_binedges[idx+1])/2 for idx in range(len(eta_binedges)-1)])
                    ecal_gap_mask = np.where((np.abs(eta_bincenters) < 1.4442) | (np.abs(eta_bincenters) > 1.5660))
                    eta_xerrs = np.array([eta_bincenters[idx]-eta_binedges[idx] for idx in range(len(eta_bincenters))])
    
                    if lepSFs[era][lep][sf_type]["isAbsEta"]:
                        set_trace()
                    else:
                        #set_trace()
                        fig, ax = plt.subplots()
                        fig.subplots_adjust(hspace=.07)
   
                        #set_trace()
                        col_idx = 0
                        for pt_idx in range(cen_lookup._axes[1].size - 1):
                            ptmin, ptmax = cen_lookup._axes[1][pt_idx], cen_lookup._axes[1][pt_idx+1]
                            if ptmax < 30.: continue
                            pt_cen_vals = cen_lookup._values[:, pt_idx]
                            pt_err_vals = errTOT_lookup._values[:, pt_idx]
                            #set_trace()
    
                            ax.errorbar(eta_bincenters[ecal_gap_mask], pt_cen_vals[ecal_gap_mask], xerr=eta_xerrs[ecal_gap_mask], yerr=pt_err_vals[ecal_gap_mask], fmt="none", color=color_cycler[col_idx], capsize=4, label=f"${ptmin} \\leq p_T \\leq {ptmax}$ GeV")
                            col_idx += 1

                        ax.set_ylabel("Scale Factor")
                        ax.set_xlabel("$\\eta_{SC}$")
                        ax.legend(loc="upper right", title=era_styles[era]["label"])
                        ax.autoscale()
                        ax.axhline(1, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
                        ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.05)
                        ax.set_xlim(eta_binedges[0], eta_binedges[-1])
    
                            # add title
                        plt.title(sf_title, loc="center", fontsize=22)
    
                        #set_trace()
                        figname = os.path.join(outdir, f"{lep}SFs_{sf_type}_{era_styles[era]['label']}")
                        fig.savefig(figname)
                        print(f"{figname} written")
                        plt.close(fig)
                else:
                    raise ValueError("Only dictionaries and dense_lookup tables are supported as SF formats")
