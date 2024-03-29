import time
tic = time.time()

import numpy as np
from coffea.lookup_tools.root_converters import convert_histo_root_file
from coffea.lookup_tools.dense_lookup import dense_lookup
from coffea.util import load, save
from pdb import set_trace
import os

proj_dir = os.environ["PROJECT_DIR"]
base_jobid = os.environ["base_jobid"]
eos_dir = os.environ["eos_dir"]
analyzer = "meta_info"

outdir = os.path.join(proj_dir, "Corrections", base_jobid) 
if not os.path.isdir(outdir):
    os.makedirs(outdir)

pu_files = {
    "2016APV" : {
        "Down" : "PileupHistogram-goldenJSON-13tev-2016-preVFP-66000ub-99bins.root",
        "Cen"  : "PileupHistogram-goldenJSON-13tev-2016-preVFP-69200ub-99bins.root",
        "Up"   : "PileupHistogram-goldenJSON-13tev-2016-preVFP-72400ub-99bins.root",
    },
    "2016" : {
        "Down" : "PileupHistogram-goldenJSON-13tev-2016-postVFP-66000ub-99bins.root",
        "Cen"  : "PileupHistogram-goldenJSON-13tev-2016-postVFP-69200ub-99bins.root",
        "Up"   : "PileupHistogram-goldenJSON-13tev-2016-postVFP-72400ub-99bins.root",
    },
    "2017" : {
        "Down" : "PileupHistogram-goldenJSON-13tev-2017-66000ub-99bins.root",
        "Cen"  : "PileupHistogram-goldenJSON-13tev-2017-69200ub-99bins.root",
        "Up"   : "PileupHistogram-goldenJSON-13tev-2017-72400ub-99bins.root",
    },
    "2018" : {
        "Down" : "PileupHistogram-goldenJSON-13tev-2018-66000ub-99bins.root",
        "Cen"  : "PileupHistogram-goldenJSON-13tev-2018-69200ub-99bins.root",
        "Up"   : "PileupHistogram-goldenJSON-13tev-2018-72400ub-99bins.root",
    },
}

pu_path = os.path.join(proj_dir, "inputs", "data", base_jobid, "Pileup")

years_to_run = ["2016APV", "2016", "2017", "2018"]
mc_pu_weights = {year: {} for year in years_to_run}
data_pu_dists = {year: {} for year in years_to_run}

for year in years_to_run:
    input_dir = os.path.join(eos_dir, "results", f"{year}_{base_jobid}", analyzer)
    fnames = [os.path.join(input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith("TOT.coffea")]

        # get nominal distribution
    data_pu_central = convert_histo_root_file(os.path.join(pu_path, year, pu_files[year]["Cen"]))
    central_hist = dense_lookup(*data_pu_central[("pileup", "dense_lookup")])
    central_hist._values = central_hist._values/sum(central_hist._values) # normalize values
    data_pu_dists[year]["central"] = central_hist
        # get up variation
    data_pu_up = convert_histo_root_file(os.path.join(pu_path, year, pu_files[year]["Up"]))
    up_hist = dense_lookup(*data_pu_up[("pileup", "dense_lookup")])
    up_hist._values = up_hist._values/sum(up_hist._values) # normalize values
    data_pu_dists[year]["up"] = up_hist
        # get down variation
    data_pu_dw = convert_histo_root_file(os.path.join(pu_path, year, pu_files[year]["Down"]))
    dw_hist = dense_lookup(*data_pu_dw[("pileup", "dense_lookup")])
    dw_hist._values = dw_hist._values/sum(dw_hist._values) # normalize values
    data_pu_dists[year]["down"] = dw_hist
    
    for fname in fnames:
        if not os.path.isfile(fname):
            raise IOError(f"{fname} not found")
        hists = load(fname)
    
        pu_histo = hists["PU_nTrueInt"]
        for dataset in sorted(pu_histo.values().keys()):
            if "data_Single" in dataset[0]: continue
            histo = pu_histo[dataset].integrate("dataset")
            mc_vals = histo.values()[()]
            mc_vals = mc_vals/sum(mc_vals) # normalize values
            edges = histo.axes()[-1].edges()
            mc_pu_weights[year][dataset[0]] = {}
            for sys_var in ["central", "up", "down"]:
                #set_trace()
                mc_weights = data_pu_dists[year][sys_var]._values/mc_vals[:data_pu_dists[year][sys_var]._values.size]
                mc_weights[mc_weights == np.inf] = np.nan
                mc_weights = np.nan_to_num(mc_weights)
                mc_pu_weights[year][dataset[0]][sys_var] = dense_lookup(mc_weights, edges)

    # save files
mcweights_name = os.path.join(outdir, "MC_PU_Weights.coffea")
save(mc_pu_weights, mcweights_name)
print(f"\n{mcweights_name} written")

data_pu_name = os.path.join(outdir, "data_PU_dists.coffea")
save(data_pu_dists, data_pu_name)
print(f"\n{data_pu_name} written")

toc = time.time()
print("Total time: %.1f" % (toc - tic))
