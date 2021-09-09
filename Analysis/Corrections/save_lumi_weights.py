from coffea.util import load, save
from pdb import set_trace
import os
from fnmatch import fnmatch
import Utilities.prettyjson as prettyjson

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
base_jobid = os.environ["base_jobid"]

outdir = os.path.join(proj_dir, "Corrections", base_jobid)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

data_lumi = prettyjson.loads(open(os.path.join(proj_dir, "inputs", "%s_lumis_data.json" % base_jobid)).read()) # file with integrated luminosity for all three years
signal_xsecs = prettyjson.loads(open(os.path.join(proj_dir, "inputs", "signal_xsecs.json")).read()) # file with signal cross sections

years_to_run = ["2016", "2017", "2018"] if base_jobid == "NanoAODv6" else ["2016APV", "2016", "2017", "2018"]
lumi_weights = {year:{"Electrons" : {}, "Muons" : {}} for year in years_to_run}

# for each year, read sumGenWeights from all meta.json files
for year in years_to_run:
    print(year)
    xsec_file = prettyjson.loads(open(os.path.join(proj_dir, "inputs", "samples_%s_%s.json" % (year, base_jobid))).read()) # file with cross sections
    datasets = list(filter(lambda x: fnmatch(x["name"], "*"), xsec_file))
    for dataset in datasets:
        sample = dataset["name"]
        if sample.startswith("data_Single"): continue
        if dataset["DBSName"] == "NOT PRESENT":
            print(f"Dataset {sample} not present, will be skipped")
            continue

        if "Int" in sample:
            #set_trace()
            for wt_type in ["pos", "neg"]:
                if not os.path.isfile(os.path.join(proj_dir, "inputs", "%s_%s" % (year, base_jobid), "%s_%s.meta.json" % (sample, wt_type))):
                    print(f"No meta.json file found for dataset {sample}_{wt_type}, skipping")
                    continue
                print(f"\t{sample}_{wt_type}")
                meta_json = prettyjson.loads(open(os.path.join(proj_dir, "inputs", "%s_%s" % (year, base_jobid), "%s_%s.meta.json" % (sample, wt_type))).read())
                sumGenWeights = meta_json["sumGenWeights"]
                xsec = signal_xsecs[sample]
                for lep in ["Electrons", "Muons"]:
                    lumi_weights[year][lep]["%s_%s" % (sample, wt_type)] = data_lumi[year][lep]/abs(sumGenWeights/xsec)

        else:
            if not os.path.isfile(os.path.join(proj_dir, "inputs", "%s_%s" % (year, base_jobid), "%s.meta.json" % sample)):
                print(f"No meta.json file found for dataset {sample}, skipping")
                continue
            print(f"\t{sample}")
            meta_json = prettyjson.loads(open(os.path.join(proj_dir, "inputs", "%s_%s" % (year, base_jobid), "%s.meta.json" % sample)).read())
            sumGenWeights = meta_json["sumGenWeights"]
            xsec = signal_xsecs[sample] if (sample.startswith("AtoTT") or sample.startswith("HtoTT")) else dataset["xsection"]
            for lep in ["Electrons", "Muons"]:
                lumi_weights[year][lep][sample] = data_lumi[year][lep]/(sumGenWeights/xsec)

    print(f"{year} calculated")

#set_trace()
    # save files
mcweights_name = os.path.join(outdir, "MC_LumiWeights.coffea")
save(lumi_weights, mcweights_name)
print(f"\n{mcweights_name} written")
