import numpy as np
from coffea.util import save
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

## pt and eta lists are hardcoded for now
    #"Electrons" : {
    #    "pt" : {
    #        "reco_id" : ["SF_ElReco_ElTOT_0", "SF_ElReco_ElTOT_1", "SF_ElReco_ElTOT_2", "SF_ElReco_ElTOT_3", "SF_ElReco_ElTOT_4", "SF_ElReco_ElTOT_5", "SF_ElReco_ElTOT_6", "SF_ElReco_ElTOT_7", "SF_ElReco_ElTOT_8", "SF_ElReco_ElTOT_9"],
    #        "trig" : ["SF_ElISOTRG_0", "SF_ElISOTRG_1", "SF_ElISOTRG_2", "SF_ElISOTRG_3", "SF_ElISOTRG_4", "SF_ElISOTRG_5", "SF_ElISOTRG_6", "SF_ElISOTRG_7", "SF_ElISOTRG_8", "SF_ElISOTRG_9"],
    #    },
    #    "eta" : "eta_ElTOT",
    #    "fnames" : {
    #        #"2016APV" : "SF_El_V161400b.root",
    #        #"2016" : "SF_El_V162400b.root", 
    #        #"2017" : "SF_El_V170400b.root",
    #        #"2018" : "SF_El_V180400b.root",
    #        "2016APV" : "SF_El_V161400.root",
    #        "2016" : "SF_El_V162400.root",
    #        "2017" : "SF_El_V170400.root",
    #        "2018" : "SF_El_V180400.root",
    #    }
    #},
leptons = {
    "Muons" : {
        "2016APV" : {
            "ID" : ["Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ID.root", "NUM_TightID_DEN_TrackerMuons_abseta_pt"],
            "ISO" : ["Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ISO.root", "NUM_TightRelIso_DEN_TightIDandIPCut_abseta_pt"],
            "TRIGGER" : ["Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_SingleMuonTriggers.root", "NUM_IsoMu24_or_IsoTkMu24_DEN_CutBasedIdTight_and_PFIsoTight_eta_pt"],
            "RECO" : ["Efficiency_muon_generalTracks_Run2016preVFP_UL_trackerMuon.root", "NUM_TrackerMuons_DEN_genTracks"],
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
    }
}

leps_to_run = sorted(leptons.keys())
years_to_run = sorted(leptons[leps_to_run[0]].keys())
sfs_to_run = sorted(leptons[leps_to_run[0]][years_to_run[0]].keys())
sf_output = {
    year : {
        lep : {
            sf_type : {
                "Central" : {},
                "Error_tot" : {},
                "Error_stat" : {},
                "Error_syst" : {},
            } for sf_type in sfs_to_run
        } for lep in leps_to_run
    } for year in years_to_run
}

#set_trace()
for lep in leptons.keys():
    for year in leptons[lep].keys():
        for sf_type, (sf_fname, sf_dist) in leptons[lep][year].items():
            print(f"Making SFs for {lep}, {year}, {sf_type}")
            sf_file = convert_histo_root_file(os.path.join(indir, year, sf_fname))

            sf_output[year][lep][sf_type]["Central"] = dense_lookup(*sf_file[(sf_dist, "dense_lookup")])
            sf_output[year][lep][sf_type]["Error_tot"] = dense_lookup(*sf_file[(f"{sf_dist}_error", "dense_lookup")])
            if sf_type != "RECO":
                sf_output[year][lep][sf_type]["Error_stat"] = dense_lookup(*sf_file[(f"{sf_dist}_stat_error", "dense_lookup")])
                sf_output[year][lep][sf_type]["Error_syst"] = dense_lookup(*sf_file[(f"{sf_dist}_syst_error", "dense_lookup")])

set_trace()
lepSF_name = os.path.join(outdir, "lepton_POG_SFs.coffea")
save(sf_output, lepSF_name)    
print(f"{lepSF_name} written")
