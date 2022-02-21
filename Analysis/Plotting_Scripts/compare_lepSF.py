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
from coffea.util import load#, save
from pdb import set_trace
import os
import Utilities.prettyjson as prettyjson
import Utilities.plot_tools as plt_tools
import numpy as np
from coffea.lookup_tools.root_converters import convert_histo_root_file

import numpy as np
from coffea.lookup_tools.dense_lookup import dense_lookup

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
base_jobid = os.environ["base_jobid"]
analyzer = "Compare_LepSFs"

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="Specify which year to run over")
#parser.add_argument("sf", choices=["ID", "ISO", "ALL"], help="Specify which scale factor you want to plot")
parser.add_argument("sf", choices=["RECO", "ID", "ISO", "TRIGGER", "ALL"], help="Specify which scale factor you want to plot")
args = parser.parse_args()


## pt and eta lists are hardcoded for now
POG_leptons = {
    "Muons" : {
        "2016APV" : {
            "RECO" : ["Efficiency_muon_generalTracks_Run2016preVFP_UL_trackerMuon.root", "NUM_TrackerMuons_DEN_genTracks"],
            "ID" : ["Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ID.root", "NUM_TightID_DEN_TrackerMuons_abseta_pt"],
            "ISO" : ["Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ISO.root", "NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt"],
            "TRIGGER" : ["Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_SingleMuonTriggers.root", "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight_eta_pt"],
        },
        "2016" : {
            "ID" : ["Efficiencies_muon_generalTracks_Z_Run2016_UL_ID.root", "NUM_TightID_DEN_TrackerMuons_abseta_pt"],
            "ISO" : ["Efficiencies_muon_generalTracks_Z_Run2016_UL_ISO.root", "NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt"],
            "TRIGGER" : ["Efficiencies_muon_generalTracks_Z_Run2016_UL_SingleMuonTriggers.root", "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight_eta_pt"],
            "RECO" : ["Efficiency_muon_generalTracks_Run2016postVFP_UL_trackerMuon.root", "NUM_TrackerMuons_DEN_genTracks"],
        },
        "2017" : {
            "ID" : ["Efficiencies_muon_generalTracks_Z_Run2017_UL_ID.root", "NUM_TightID_DEN_TrackerMuons_abseta_pt"],
            "ISO" : ["Efficiencies_muon_generalTracks_Z_Run2017_UL_ISO.root", "NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt"],
            "TRIGGER" : ["Efficiencies_muon_generalTracks_Z_Run2017_UL_SingleMuonTriggers.root", "NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight_eta_pt"],
            "RECO" : ["Efficiency_muon_generalTracks_Run2017_UL_trackerMuon.root", "NUM_TrackerMuons_DEN_genTracks"],
        },
        "2018" : {
            "ID" : ["Efficiencies_muon_generalTracks_Z_Run2018_UL_ID.root", "NUM_TightID_DEN_TrackerMuons_abseta_pt"],
            "ISO" : ["Efficiencies_muon_generalTracks_Z_Run2018_UL_ISO.root", "NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt"],
            "TRIGGER" : ["Efficiencies_muon_generalTracks_Z_Run2018_UL_SingleMuonTriggers.root", "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight_eta_pt"],
            "RECO" : ["Efficiency_muon_generalTracks_Run2018_UL_trackerMuon.root", "NUM_TrackerMuons_DEN_genTracks"],
        },
    },
    "Electrons" : {
        "2016APV" : {
            "ID" : ["egammaEffi.txt_Ele_Tight_preVFP_EGM2D.root", "EGamma_SF2D"],
            ##"ISO" : ["", ""],
            ##"TRIGGER" : ["", ""],
            #"RECO" : ["egammaEffi_ptAbove20.txt_EGM2D_UL2016preVFP.root", "EGamma_SF2D"],
        },
        "2016" : {
            "ID" : ["egammaEffi.txt_Ele_Tight_postVFP_EGM2D.root", "EGamma_SF2D"],
            ##"ISO" : ["", ""],
            ##"TRIGGER" : ["", ""],
            #"RECO" : ["egammaEffi_ptAbove20.txt_EGM2D_UL2016postVFP.root", "EGamma_SF2D"],
        },
        "2017" : {
            "ID" : ["egammaEffi.txt_EGM2D_Tight_UL17.root", "EGamma_SF2D"],
            ##"ISO" : ["", ""],
            ##"TRIGGER" : ["", ""],
            #"RECO" : ["egammaEffi_ptAbove20.txt_EGM2D_UL2017.root", "EGamma_SF2D"],
        },
        "2018" : {
            "ID" : ["egammaEffi.txt_Ele_Tight_EGM2D.root", "EGamma_SF2D"],
            ##"ISO" : ["", ""],
            ##"TRIGGER" : ["", ""],
            #"RECO" : ["egammaEffi_ptAbove20.txt_EGM2D_UL2018.root", "EGamma_SF2D"],
        },
    }
}
Otto_leptons = {
    "Muons" : {
        "eta" : "eta_MuTOT",
        "2016APV" : "SF_Mu_V161400.root",
        "2016" : "SF_Mu_V162400.root",
        "2017" : "SF_Mu_V170400.root",
        "2018" : "SF_Mu_V180400.root",
        "RECO" : "SF_MuTRK",
        "ID" : "SF_MuSA_MuID",
        #"ID" : "SF_MuID",
        #"ISO" : "SF_MuTOT", # parts of muon reco+ID+ISO sf
        "ISO" : "SF_MuISO", # ID+ISO sf
        #"ISO" : "SF_MuID_MuISO", # ID+ISO sf
        "TRIGGER" : "SF_MuISOTRG"
    },
    "Electrons" : {
        "eta" : "eta_ElTOT",
        "2016APV" : "SF_El_V161400.root",
        "2016" : "SF_El_V162400.root",
        "2017" : "SF_El_V170400.root",
        "2018" : "SF_El_V180400.root",
        "ID" : "SF_ElTOT",
        #"TRIGGER" : "SF_ElISOTRG",
        "RECO" : "SF_ElReco",
    }
}


SFs_to_run = ["RECO", "ID", "ISO", "TRIGGER"] if args.sf == "ALL" else [args.sf]
#SFs_to_run = ["ID", "ISO"] if args.sf == "ALL" else [args.sf]

titles = {
    "Muons" : {
        "ID" : "Muon Tight ID",
        "ISO" : "Muon TightRelIso + Tight ID",
        "TRIGGER" : "Muon Trigger",
        "RECO" : "Muon Reconstruction",
    },
    "Electrons" : {
        "ID" : "Electron Tight cut-based ID",
        "RECO" : "Electron Reconstruction",
    }
}

pog_abseta_styles = {"color":"k", "linestyle":"-", "label":"POG, $\\left| \\eta \\right|$: [LOW, HI]"}
pog_eta_styles = {"color":"k", "linestyle":"-", "label":"POG, $\\eta$: [LOW, HI]"}

absEta_to_eta_bins = {
    "0" : ["3", "4"],
    "1" : ["2", "5"],
    "2" : ["1", "6"],
    "3" : ["0", "7"],
}

otto_dir_to_use = "leptonSF_noTAGvariation"
input_dir = os.path.join(proj_dir, "inputs", "data", base_jobid, "lepSFs")
outdir = os.path.join(proj_dir, "plots", base_jobid, "Compare_LepSFs", otto_dir_to_use)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

#set_trace()
for leptype in POG_leptons.keys():
    for sf_type, (pog_sf_fname, pog_sf_dist) in POG_leptons[leptype][args.year].items():
        if sf_type not in SFs_to_run: continue
            # get values and errors for POG SFs
        pog_sf_file = convert_histo_root_file(os.path.join(input_dir, "POG_SFs", args.year, leptype, pog_sf_fname))
        pog_cen_vals, (pog_eta_binning, pog_pt_binning) = pog_sf_file[(pog_sf_dist, "dense_lookup")]
        pog_err_vals = pog_sf_file[(f"{pog_sf_dist}_error", "dense_lookup")][0]
        isAbsEta = np.all(pog_eta_binning >= 0)

            # get values and errors for Otto's SFs
        otto_sf_file = convert_histo_root_file(os.path.join(input_dir, otto_dir_to_use, Otto_leptons[leptype][args.year]))
        otto_eta_binning = otto_sf_file[(Otto_leptons[leptype]["eta"], "dense_lookup")][1][0]

        #set_trace()        
        # make plots in terms of eta
        if isAbsEta:
            for eta_bin_idx in range(len(pog_eta_binning)-1):
                # set style for pog vals based on eta
                pog_eta_style = pog_abseta_styles.copy()
                pog_eta_style["label"] = pog_eta_style["label"].replace("LOW", str(pog_eta_binning[eta_bin_idx])).replace("HI", str(pog_eta_binning[eta_bin_idx+1]))

                neg_eta_sf_cen_vals, neg_eta_sf_pt_bins = otto_sf_file[(f"{Otto_leptons[leptype][sf_type]}_{absEta_to_eta_bins[str(eta_bin_idx)][0]}", "dense_lookup")]
                neg_eta_sf_err_vals = otto_sf_file[(f"{Otto_leptons[leptype][sf_type]}_{absEta_to_eta_bins[str(eta_bin_idx)][0]}_error", "dense_lookup")][0]
                pos_eta_sf_cen_vals, pos_eta_sf_pt_bins = otto_sf_file[(f"{Otto_leptons[leptype][sf_type]}_{absEta_to_eta_bins[str(eta_bin_idx)][1]}", "dense_lookup")]
                pos_eta_sf_err_vals = otto_sf_file[(f"{Otto_leptons[leptype][sf_type]}_{absEta_to_eta_bins[str(eta_bin_idx)][1]}_error", "dense_lookup")][0]

                #set_trace()
                    # convert sf values to arrays of equal length since pT binning is different
                unique_pt_binvals = np.unique(np.concatenate((pog_pt_binning, neg_eta_sf_pt_bins, pos_eta_sf_pt_bins), axis=None))
                min_pt_binval = np.min(unique_pt_binvals)
                #min_pt_binval = np.min(np.concatenate((pog_pt_binning, neg_eta_sf_pt_bins, pos_eta_sf_pt_bins), axis=None))
                max_pt_binval = unique_pt_binvals[-2]+10
                #max_pt_binval = np.max(np.concatenate((pog_pt_binning, neg_eta_sf_pt_bins, pos_eta_sf_pt_bins), axis=None))
                pt_binvals = np.arange(min_pt_binval, max_pt_binval+1)

                    # set POG values
                pog_cen_vals_to_plot, pog_err_vals_to_plot = np.ones(pt_binvals.size-1)*pog_cen_vals[eta_bin_idx][0], np.ones(pt_binvals.size-1)*pog_err_vals[eta_bin_idx][0]
                for pog_pt_binIdx in range(pog_pt_binning.size-1):
                    inds_to_set = np.where(pt_binvals >= pog_pt_binning[pog_pt_binIdx])[0][:-1]
                    pog_cen_vals_to_plot[inds_to_set] = pog_cen_vals[eta_bin_idx][pog_pt_binIdx]
                    pog_err_vals_to_plot[inds_to_set] = pog_err_vals[eta_bin_idx][pog_pt_binIdx]

                #set_trace()
                    # set Otto values
                neg_eta_cen_vals_to_plot, neg_eta_err_vals_to_plot = np.ones(pt_binvals.size-1)*neg_eta_sf_cen_vals[0], np.ones(pt_binvals.size-1)*neg_eta_sf_err_vals[0]
                pos_eta_cen_vals_to_plot, pos_eta_err_vals_to_plot = np.ones(pt_binvals.size-1)*pos_eta_sf_cen_vals[0], np.ones(pt_binvals.size-1)*pos_eta_sf_err_vals[0]
                for otto_pt_binIdx in range(neg_eta_sf_pt_bins[0].size-1):
                    inds_to_set = np.where(pt_binvals >= neg_eta_sf_pt_bins[0][otto_pt_binIdx])[0][:-1]
                    neg_eta_cen_vals_to_plot[inds_to_set] = neg_eta_sf_cen_vals[otto_pt_binIdx]
                    neg_eta_err_vals_to_plot[inds_to_set] = neg_eta_sf_err_vals[otto_pt_binIdx]
                    pos_eta_cen_vals_to_plot[inds_to_set] = pos_eta_sf_cen_vals[otto_pt_binIdx]
                    pos_eta_err_vals_to_plot[inds_to_set] = pos_eta_sf_err_vals[otto_pt_binIdx]

                #set_trace()
                fig, ax = plt.subplots()
                fig.subplots_adjust(hspace=.07)

                    ## plot POG vals
                ax.step(pt_binvals, np.r_[pog_cen_vals_to_plot,pog_cen_vals_to_plot[-1]], where="post", **pog_eta_style)
                ax.fill_between(pt_binvals, np.r_[pog_cen_vals_to_plot+pog_err_vals_to_plot,(pog_cen_vals_to_plot+pog_err_vals_to_plot)[-1]],
                    np.r_[pog_cen_vals_to_plot-pog_err_vals_to_plot,(pog_cen_vals_to_plot-pog_err_vals_to_plot)[-1]], facecolor="k", step="post", alpha=0.5)
                    ## plot ottos neg eta vals
                ax.step(pt_binvals, np.r_[neg_eta_cen_vals_to_plot, neg_eta_cen_vals_to_plot[-1]], where="post",
                    **{"color":"r", "label":f"UR, $\\eta$: [{otto_eta_binning[int(absEta_to_eta_bins[str(eta_bin_idx)][0])]}, {otto_eta_binning[int(absEta_to_eta_bins[str(eta_bin_idx)][0])+1]}]"})
                ax.fill_between(pt_binvals, np.r_[neg_eta_cen_vals_to_plot+neg_eta_err_vals_to_plot, (neg_eta_cen_vals_to_plot+neg_eta_err_vals_to_plot)[-1]],
                    np.r_[neg_eta_cen_vals_to_plot-neg_eta_err_vals_to_plot, (neg_eta_cen_vals_to_plot-neg_eta_err_vals_to_plot)[-1]], facecolor="r", step="post", alpha=0.5)
                    ## plot ottos pos eta vals
                ax.step(pt_binvals, np.r_[pos_eta_cen_vals_to_plot, pos_eta_cen_vals_to_plot[-1]], where="post",
                    **{"color":"b", "label":f"UR, $\\eta$: [{otto_eta_binning[int(absEta_to_eta_bins[str(eta_bin_idx)][1])]}, {otto_eta_binning[int(absEta_to_eta_bins[str(eta_bin_idx)][1])+1]}]"})
                ax.fill_between(pt_binvals, np.r_[pos_eta_cen_vals_to_plot+pos_eta_err_vals_to_plot, (pos_eta_cen_vals_to_plot+pos_eta_err_vals_to_plot)[-1]],
                    np.r_[pos_eta_cen_vals_to_plot-pos_eta_err_vals_to_plot, (pos_eta_cen_vals_to_plot-pos_eta_err_vals_to_plot)[-1]], facecolor="b", step="post", alpha=0.5)
          
                ax.set_ylabel("SF")
                ax.set_xlabel("$p_{T}$ [GeV]")
                ax.autoscale()
                ax.legend(loc="upper right")
                ax.set_xlim(min_pt_binval, max_pt_binval)

                    # add title        
                plt.title(f"{titles[leptype][sf_type]}, {args.year}", loc="center", fontsize=18)
                
                figname = os.path.join(outdir, f"LepSFs_{args.year}_{leptype}_{sf_type}_etabin{eta_bin_idx}")
                fig.savefig(figname)
                print(f"{figname} written")
                plt.close(fig)

        else:
            for eta_bin_idx in range(len(pog_eta_binning)-1):
                #set_trace()
                if pog_eta_binning.size != otto_eta_binning.size: raise ValueError(f"POG eta binning is {pog_eta_binning}, Otto's is {otto_eta_binning}")

                pog_eta_style = pog_eta_styles.copy()
                pog_eta_style["label"] = pog_eta_style["label"].replace("LOW", str(pog_eta_binning[eta_bin_idx])).replace("HI", str(pog_eta_binning[eta_bin_idx+1]))

                otto_eta_sf_cen_vals, otto_eta_sf_pt_bins = otto_sf_file[(f"{Otto_leptons[leptype][sf_type]}_{eta_bin_idx}", "dense_lookup")]
                otto_eta_sf_err_vals = otto_sf_file[(f"{Otto_leptons[leptype][sf_type]}_{eta_bin_idx}_error", "dense_lookup")][0]

                #set_trace()
                    # convert sf values to arrays of equal length since pT binning is different
                unique_pt_binvals = np.unique(np.concatenate((pog_pt_binning, otto_eta_sf_pt_bins), axis=None))
                min_pt_binval = np.min(unique_pt_binvals)
                max_pt_binval = unique_pt_binvals[-2]+10
                pt_binvals = np.arange(min_pt_binval, max_pt_binval+1)

                    # set POG values to plot
                pog_cen_vals_to_plot, pog_err_vals_to_plot = np.ones(pt_binvals.size-1)*pog_cen_vals[eta_bin_idx][0], np.ones(pt_binvals.size-1)*pog_err_vals[eta_bin_idx][0]
                for pog_pt_binIdx in range(pog_pt_binning.size-1):
                    inds_to_set = np.where(pt_binvals >= pog_pt_binning[pog_pt_binIdx])[0][:-1]
                    pog_cen_vals_to_plot[inds_to_set] = pog_cen_vals[eta_bin_idx][pog_pt_binIdx]
                    pog_err_vals_to_plot[inds_to_set] = pog_err_vals[eta_bin_idx][pog_pt_binIdx]

                #set_trace()
                    # set Otto values to plot
                otto_eta_cen_vals_to_plot, otto_eta_err_vals_to_plot = np.ones(pt_binvals.size-1)*otto_eta_sf_cen_vals[0], np.ones(pt_binvals.size-1)*otto_eta_sf_err_vals[0]
                for otto_pt_binIdx in range(otto_eta_sf_pt_bins[0].size-1):
                    inds_to_set = np.where(pt_binvals >= otto_eta_sf_pt_bins[0][otto_pt_binIdx])[0][:-1]
                    otto_eta_cen_vals_to_plot[inds_to_set] = otto_eta_sf_cen_vals[otto_pt_binIdx]
                    otto_eta_err_vals_to_plot[inds_to_set] = otto_eta_sf_err_vals[otto_pt_binIdx]

                #set_trace()
                fig, ax = plt.subplots()
                fig.subplots_adjust(hspace=.07)

                    ## plot POG vals
                ax.step(pt_binvals, np.r_[pog_cen_vals_to_plot,pog_cen_vals_to_plot[-1]], where="post", **pog_eta_style)
                ax.fill_between(pt_binvals, np.r_[pog_cen_vals_to_plot+pog_err_vals_to_plot,(pog_cen_vals_to_plot+pog_err_vals_to_plot)[-1]],
                    np.r_[pog_cen_vals_to_plot-pog_err_vals_to_plot,(pog_cen_vals_to_plot-pog_err_vals_to_plot)[-1]], facecolor="k", step="post", alpha=0.5)
                    ## plot ottos eta vals
                ax.step(pt_binvals, np.r_[otto_eta_cen_vals_to_plot, otto_eta_cen_vals_to_plot[-1]], where="post",
                    **{"color":"r", "label":f"UR, $\\eta$: [{otto_eta_binning[eta_bin_idx]}, {otto_eta_binning[eta_bin_idx+1]}]"})
                ax.fill_between(pt_binvals, np.r_[otto_eta_cen_vals_to_plot+otto_eta_err_vals_to_plot, (otto_eta_cen_vals_to_plot+otto_eta_err_vals_to_plot)[-1]],
                    np.r_[otto_eta_cen_vals_to_plot-otto_eta_err_vals_to_plot, (otto_eta_cen_vals_to_plot-otto_eta_err_vals_to_plot)[-1]], facecolor="r", step="post", alpha=0.5)
          
                ax.set_ylabel("SF")
                ax.set_xlabel("$p_{T}$ [GeV]")
                ax.autoscale()
                ax.legend(loc="upper right")
                ax.set_xlim(min_pt_binval, max_pt_binval)

                    # add title        
                plt.title(f"{titles[leptype][sf_type]}, {args.year}", loc="center", fontsize=18)
                
                figname = os.path.join(outdir, f"LepSFs_{args.year}_{leptype}_{sf_type}_etabin{eta_bin_idx}")
                fig.savefig(figname)
                print(f"{figname} written")
                plt.close(fig)

