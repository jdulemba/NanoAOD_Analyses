import numpy as np
from coffea.util import save, load
from coffea.lookup_tools.root_converters import convert_histo_root_file
from coffea.lookup_tools.dense_lookup import dense_lookup
from pdb import set_trace
import os

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
base_jobid = os.environ["base_jobid"]

indir = os.path.join(proj_dir, "inputs", "data", base_jobid, "lepSFs", "POG_SFs")
outdir = os.path.join(proj_dir, "Corrections", jobid)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

leptons = {
    "2016APV" : {
        "Muons" : {
            "ID" : ["Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ID.root", "NUM_TightID_DEN_TrackerMuons_abseta_pt"],
            "ISO" : ["Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ISO.root", "NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt"],
            "TRIG" : ["Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_SingleMuonTriggers.root", "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight_eta_pt"],
            "RECO" : ["Efficiency_muon_generalTracks_Run2016preVFP_UL_trackerMuon.root", "NUM_TrackerMuons_DEN_genTracks"],
        },
        "Electrons" : {
            "ID" : ["egammaEffi.txt_Ele_Tight_preVFP_EGM2D.root", "EGamma_SF2D"],
            "TRIG" : ["", ""],
            "RECO" : ["egammaEffi_ptAbove20.txt_EGM2D_UL2016preVFP.root", "EGamma_SF2D"],
        },
    },
    "2016" : {
        "Muons" : {
            "ID" : ["Efficiencies_muon_generalTracks_Z_Run2016_UL_ID.root", "NUM_TightID_DEN_TrackerMuons_abseta_pt"],
            "ISO" : ["Efficiencies_muon_generalTracks_Z_Run2016_UL_ISO.root", "NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt"],
            "TRIG" : ["Efficiencies_muon_generalTracks_Z_Run2016_UL_SingleMuonTriggers.root", "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight_eta_pt"],
            "RECO" : ["Efficiency_muon_generalTracks_Run2016postVFP_UL_trackerMuon.root", "NUM_TrackerMuons_DEN_genTracks"],
        },
        "Electrons" : {
            "ID" : ["egammaEffi.txt_Ele_Tight_postVFP_EGM2D.root", "EGamma_SF2D"],
            "TRIG" : ["", ""],
            "RECO" : ["egammaEffi_ptAbove20.txt_EGM2D_UL2016postVFP.root", "EGamma_SF2D"],
        },
    },
    "2017" : {
        "Muons" : {
            "ID" : ["Efficiencies_muon_generalTracks_Z_Run2017_UL_ID.root", "NUM_TightID_DEN_TrackerMuons_abseta_pt"],
            "ISO" : ["Efficiencies_muon_generalTracks_Z_Run2017_UL_ISO.root", "NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt"],
            "TRIG" : ["Efficiencies_muon_generalTracks_Z_Run2017_UL_SingleMuonTriggers.root", "NUM_IsoMu27_DEN_CutBasedIdTight_and_PFIsoTight_eta_pt"],
            "RECO" : ["Efficiency_muon_generalTracks_Run2017_UL_trackerMuon.root", "NUM_TrackerMuons_DEN_genTracks"],
        },
        "Electrons" : {
            "ID" : ["egammaEffi.txt_EGM2D_Tight_UL17.root", "EGamma_SF2D"],
            "TRIG" : ["", ""],
            "RECO" : ["egammaEffi_ptAbove20.txt_EGM2D_UL2017.root", "EGamma_SF2D"],
        },
    },
    "2018" : {
        "Muons" : {
            "ID" : ["Efficiencies_muon_generalTracks_Z_Run2018_UL_ID.root", "NUM_TightID_DEN_TrackerMuons_abseta_pt"],
            "ISO" : ["Efficiencies_muon_generalTracks_Z_Run2018_UL_ISO.root", "NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt"],
            "TRIG" : ["Efficiencies_muon_generalTracks_Z_Run2018_UL_SingleMuonTriggers.root", "NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight_eta_pt"],
            "RECO" : ["Efficiency_muon_generalTracks_Run2018_UL_trackerMuon.root", "NUM_TrackerMuons_DEN_genTracks"],
        },
        "Electrons" : {
            "ID" : ["egammaEffi.txt_Ele_Tight_EGM2D.root", "EGamma_SF2D"],
            "TRIG" : ["", ""],
            "RECO" : ["egammaEffi_ptAbove20.txt_EGM2D_UL2018.root", "EGamma_SF2D"],
        },
    },
}

sf_output = {
    year : {
        lep : {
            sf_type : {
                "Central" : {},
                "Error_tot" : {},
                "Error_stat" : {},
                "Error_syst" : {},
            } for sf_type in leptons[year][lep].keys()
        } for lep in leptons[year].keys()
    } for year in leptons.keys()
}

    # load UR lepton SFs to add electron triggers
UR_SFs = load(os.path.join(proj_dir, "Corrections", jobid, "leptonSF_noTAGvariation.coffea"))


def average_sfs(sf_file, sf_dist):
    central_vals = dense_lookup(*sf_file[(sf_dist, "dense_lookup")])
    error_tot_vals = dense_lookup(*sf_file[(f"{sf_dist}_error", "dense_lookup")])

        # find bins where pt >= 5 GeV
    valid_ptbins = np.where(central_vals._axes[1] >= 5)[0][:-1]
    #print(f"Central: {central_vals._values[:, valid_ptbins]}")
    #print(f"Error: {error_tot_vals._values[:, valid_ptbins]}")
        # find average values and errors weighted by 1./variance for each eta bin
    sumeta_numerator = np.sum(np.where(error_tot_vals._values[:, valid_ptbins] == 0., 0., central_vals._values[:, valid_ptbins]/(error_tot_vals._values[:, valid_ptbins]**2)), axis=1)
    sumeta_denom = np.sum(np.where(error_tot_vals._values[:, valid_ptbins] == 0., 0., 1./(error_tot_vals._values[:, valid_ptbins]**2)), axis=1)
    average_eta_vals, average_eta_errors = sumeta_numerator/sumeta_denom, 1./np.sqrt(sumeta_denom)
    central_lookup = dense_lookup(average_eta_vals, central_vals._axes[0])
    toterror_lookup = dense_lookup(average_eta_errors, error_tot_vals._axes[0])
    if (f"{sf_dist}_stat_error", "dense_lookup") in sf_file.keys():
        set_trace()
    if (f"{sf_dist}_syst_error", "dense_lookup") in sf_file.keys():
        set_trace()

    return central_lookup, toterror_lookup


for year in leptons.keys():
    for lep in leptons[year].keys():
        for sf_type, (sf_fname, sf_dist) in leptons[year][lep].items():
            print(f"Making SFs for {year}, {lep}, {sf_type}")

                # save UR electron trigger SFs
            if (lep == "Electrons") and (sf_type == "TRIG"):
                    # make SFs in 2018 start at 34 GeV instead of 35
                if year == "2018":
                    orig_cen_lookup = UR_SFs[year][lep]["Trig"]["Central"]
                    orig_err_lookup = UR_SFs[year][lep]["Trig"]["Error"]
                        # set dense lookups with updated bins
                    sf_output[year][lep][sf_type]["Central"] = {eta_bin : dense_lookup(lookup._values[np.where(lookup._axes >= 35.)[0][0]:], np.r_[np.array([34.0]), lookup._axes[np.where(lookup._axes >= 35.)[0][0]+1:]]) for eta_bin, lookup in orig_cen_lookup.items()}
                    sf_output[year][lep][sf_type]["Error_tot"] = {eta_bin : dense_lookup(lookup._values[np.where(lookup._axes >= 35.)[0][0]:], np.r_[np.array([34.0]), lookup._axes[np.where(lookup._axes >= 35.)[0][0]+1:]]) for eta_bin, lookup in orig_err_lookup.items()}
                else:                
                    sf_output[year][lep][sf_type]["Central"] = UR_SFs[year][lep]["Trig"]["Central"]
                    sf_output[year][lep][sf_type]["Error_tot"] = UR_SFs[year][lep]["Trig"]["Error"]
                sf_output[year][lep][sf_type]["eta_ranges"] = UR_SFs[year][lep]["eta_ranges"]
                sf_output[year][lep][sf_type]["isAbsEta"] = UR_SFs[year][lep]["eta_ranges"][0][0] >= 0
                del sf_output[year][lep][sf_type]["Error_stat"]
                del sf_output[year][lep][sf_type]["Error_syst"]
            else:
                sf_file = convert_histo_root_file(os.path.join(indir, year, lep, sf_fname))

                sf_output[year][lep][sf_type]["isAbsEta"] = np.all(sf_file[(sf_dist, "dense_lookup")][1][0] >= 0)
                if (lep == "Muons") and (sf_type == "RECO"):
                    central_lookup, tot_error_lookup = average_sfs(sf_file, sf_dist)
                    #set_trace()
                    sf_output[year][lep][sf_type]["Central"] = central_lookup
                    sf_output[year][lep][sf_type]["Error_tot"] = tot_error_lookup
                else:
                    sf_output[year][lep][sf_type]["Central"] = dense_lookup(*sf_file[(sf_dist, "dense_lookup")])
                    sf_output[year][lep][sf_type]["Error_tot"] = dense_lookup(*sf_file[(f"{sf_dist}_error", "dense_lookup")])

                if (sf_type != "RECO") and (lep == "Muons"):
                    sf_output[year][lep][sf_type]["Error_stat"] = dense_lookup(*sf_file[(f"{sf_dist}_stat_error", "dense_lookup")])
                    sf_output[year][lep][sf_type]["Error_syst"] = dense_lookup(*sf_file[(f"{sf_dist}_syst_error", "dense_lookup")])
                else:
                    del sf_output[year][lep][sf_type]["Error_stat"]
                    del sf_output[year][lep][sf_type]["Error_syst"]

#set_trace()
lepSF_name = os.path.join(outdir, f"lepton_{jobid}.coffea")
#lepSF_name = os.path.join(outdir, "lepton_POG_SFs.coffea")
save(sf_output, lepSF_name)    
print(f"{lepSF_name} written")
