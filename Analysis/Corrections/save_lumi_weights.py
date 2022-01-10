from coffea.util import load, save
from pdb import set_trace
import os
from fnmatch import fnmatch
import Utilities.prettyjson as prettyjson

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--kfactors", action="store_true", help="Apply signal k-factors to signal")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
base_jobid = os.environ["base_jobid"]

outdir = os.path.join(proj_dir, "Corrections", base_jobid)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

data_lumi = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read()) # file with integrated luminosity for all three years
signal_xsecs = prettyjson.loads(open(os.path.join(proj_dir, "inputs", "signal_xsecs.json")).read()) # file with signal cross sections
signal_xabs = prettyjson.loads(open(os.path.join(proj_dir, "inputs", "signal_xabs.json")).read()) # file with signal cross sections
signal_kfactors = prettyjson.loads(open(os.path.join(proj_dir, "inputs", "signal_kfactors.json")).read()) # file with signal NNLO/LO k-factors

LHEScale_wts_dict = {
    "uF_up" : 5,
    "uF_down" : 3,
    "uR_up" : 7,
    "uR_down" : 1,
    "uF_up_uR_up" : 8,
    "uF_down_uR_down" : 0,
    #"uF_up_uR_down" : 6,
    #"uF_down_uR_up" : 2,
}
# ttbar cross section values and scale uncertainties for 13 TeV mt = 172.5 https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO
tt_NNLO_xsec = {
    "central" : 831.76,
    "up" : 19.77,
    "down" : -29.20,
}
tt_brs = {"ttJetsDiLep" : 0.105, "ttJetsHad" : 0.457, "ttJetsSL" : 0.438}
#PS_wts_dict = {
#    "ISRUp" : 0,
#    "ISRDown" : 2,
#    "FSRUp" : 1,
#    "FSRDown" : 3,
#}

nominal_ttJets = ["ttJetsDiLep", "ttJetsHad", "ttJetsSL"]

years_to_run = ["2016APV", "2016", "2017", "2018"]
lumi_weights = {year:{"Electrons" : {}, "Muons" : {}} for year in years_to_run}

# for each year, read sumGenWeights from all meta.json files
for year in years_to_run:
    print(year)
    xsec_file = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"samples_{year}_{base_jobid}.json")).read()) # file with cross sections
    datasets = list(filter(lambda x: fnmatch(x["name"], "*"), xsec_file))
    for dataset in datasets:
        sample = dataset["name"]
        if sample.startswith("data_Single"): continue
        if dataset["DBSName"] == "NOT PRESENT":
            #set_trace()
            if os.path.isfile(os.path.join(proj_dir, "inputs", f"{year}_{base_jobid}", f"{sample}.meta.json")):
                #set_trace()
                print(f"\t{sample}")
                meta_json = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{year}_{base_jobid}", f"{sample}.meta.json")).read())
                sumGenWeights = meta_json["sumGenWeights"]
                xsec = signal_xsecs[sample] if (sample.startswith("AtoTT") or sample.startswith("HtoTT")) else dataset["xsection"]
                for lep in ["Electrons", "Muons"]:
                    lumi_weights[year][lep][sample] = data_lumi[year][lep]/(sumGenWeights/xsec)
            else:
                print(f"Dataset {sample} not present, will be skipped")
            continue

        if "Int" in sample:
            #set_trace()
            print(f"\t{sample}")
                # get pos wts
            if not os.path.isfile(os.path.join(proj_dir, "inputs", f"{year}_{base_jobid}", f"{sample}_pos.meta.json")):
                print(f"No meta.json file found for dataset {sample}_pos, skipping")
                continue
            meta_json_pos = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{year}_{base_jobid}", f"{sample}_pos.meta.json")).read())
            sumGenWeights_pos = meta_json_pos["sumGenWeights"]
                # get neg wts
            if not os.path.isfile(os.path.join(proj_dir, "inputs", f"{year}_{base_jobid}", f"{sample}_neg.meta.json")):
                print(f"No meta.json file found for dataset {sample}_neg, skipping")
                continue
            meta_json_neg = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{year}_{base_jobid}", f"{sample}_neg.meta.json")).read())
            sumGenWeights_neg = meta_json_neg["sumGenWeights"]

            xabs = signal_xabs[sample]
            int_scale = xabs/(abs(sumGenWeights_pos) + abs(sumGenWeights_neg))
            if args.kfactors:
                int_scale = int_scale*signal_kfactors[sample]
            for lep in ["Electrons", "Muons"]:
                for wt_type in ["pos", "neg"]:
                    lumi_weights[year][lep][f"{sample}_{wt_type}"] = data_lumi[year][lep]*int_scale

                # add LHE scale weights
            #set_trace()
            for lhe_wt_name, idx in LHEScale_wts_dict.items():
                lhe_wts_pos = meta_json_pos["sumLHEscaleWeights"][idx]
                lhe_wts_neg = meta_json_neg["sumLHEscaleWeights"][idx]
                lhe_wts_xabs = signal_xabs[f"{sample}_{lhe_wt_name}"]
                lhe_int_scale = lhe_wts_xabs/(abs(lhe_wts_pos) + abs(lhe_wts_neg))
                if args.kfactors:
                    lhe_int_scale = lhe_int_scale*signal_kfactors[f"{sample}_{lhe_wt_name}"]
                for lep in ["Electrons", "Muons"]:
                    for wt_type in ["pos", "neg"]:
                        lumi_weights[year][lep][f"{sample}_{wt_type}_{lhe_wt_name}"] = data_lumi[year][lep]*lhe_int_scale
                
        else:
            if not os.path.isfile(os.path.join(proj_dir, "inputs", f"{year}_{base_jobid}", f"{sample}.meta.json")):
                print(f"No meta.json file found for dataset {sample}, skipping")
                continue
            print(f"\t{sample}")
            meta_json = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{year}_{base_jobid}", f"{sample}.meta.json")).read())
            sumGenWeights = meta_json["sumGenWeights"]
            if (sample.startswith("AtoTT") or sample.startswith("HtoTT")):
                #set_trace()
                xsec = signal_xabs[sample]
                if args.kfactors:
                    xsec = xsec*signal_kfactors[sample]

                    # add LHE scale weights for signal
                for lhe_wt_name, idx in LHEScale_wts_dict.items():
                    lhe_wts = meta_json["sumLHEscaleWeights"][idx]
                    lhe_wts_xabs = signal_xabs[f"{sample}_{lhe_wt_name}"]
                    lhe_res_scale = lhe_wts_xabs/abs(lhe_wts)
                    if args.kfactors:
                        lhe_res_scale = lhe_res_scale*signal_kfactors[f"{sample}_{lhe_wt_name}"]
                    for lep in ["Electrons", "Muons"]:
                        lumi_weights[year][lep][f"{sample}_{lhe_wt_name}"] = data_lumi[year][lep]*lhe_res_scale
            else:
                xsec = dataset["xsection"]
                if sample in nominal_ttJets:
                    #set_trace()
                        # add LHE scale weights for nominal ttJets
                    for lhe_wt_name, idx in LHEScale_wts_dict.items():
                        lhe_wts = meta_json["sumLHEscaleWeights"][idx]
                        xsec_to_use = (tt_NNLO_xsec["central"]+tt_NNLO_xsec["up"])*tt_brs[sample] if "up" in lhe_wt_name else (tt_NNLO_xsec["central"]+tt_NNLO_xsec["down"])*tt_brs[sample]
                        lhe_wts_scale = xsec_to_use/lhe_wts
                        for lep in ["Electrons", "Muons"]:
                            lumi_weights[year][lep][f"{sample}_{lhe_wt_name}"] = data_lumi[year][lep]*lhe_wts_scale

                    #    # add PS weights for nominal ttJets
                    #for ps_wt_name, idx in PS_wts_dict.items():
                    #    ps_wts = meta_json["sumPSWeights"][idx]
                    #    ps_wts_scale = xsec/ps_wts
                    #    for lep in ["Electrons", "Muons"]:
                    #        lumi_weights[year][lep][f"{sample}_{ps_wt_name}"] = data_lumi[year][lep]*ps_wts_scale

            for lep in ["Electrons", "Muons"]:
                lumi_weights[year][lep][sample] = data_lumi[year][lep]*(xsec/sumGenWeights)

    print(f"{year} calculated")

#set_trace()
    # save files
mcweights_name = os.path.join(outdir, "MC_LumiWeights_kfactors.coffea" if args.kfactors else "MC_LumiWeights.coffea")
save(lumi_weights, mcweights_name)
print(f"\n{mcweights_name} written")
