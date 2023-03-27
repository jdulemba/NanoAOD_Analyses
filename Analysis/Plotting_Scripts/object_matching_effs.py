#!/usr/bin/env python

import time
tic = time.time()

from coffea.util import load
from pdb import set_trace
import os
import Utilities.prettyjson as prettyjson
import Utilities.plot_tools as plt_tools
import numpy as np

proj_dir = os.environ["PROJECT_DIR"]
base_jobid = os.environ["base_jobid"]
jobid = os.environ["jobid"]
eos_dir = os.environ["eos_dir"]
plot_outdir = os.environ["plots_dir"]
analyzer = "object_matching_effs"

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="Specify which year to run over")
args = parser.parse_args()

input_dir = os.path.join(eos_dir, "results", f"{args.year}_{jobid}", analyzer)
f_ext = "TOT.coffea"
fnames = sorted(["%s/%s" % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)])
hdict = plt_tools.add_coffea_files(fnames) if len(fnames) > 1 else load(fnames[0])
#set_trace()

outdir = os.path.join(plot_outdir, jobid, analyzer, args.year)
if not os.path.isdir(outdir):
    os.makedirs(outdir)


#histograms = ["nBestPerms", "Object_Identifying", "ValidBestPerms"]
histograms = ["ValidBestPerms"]
#histograms = ["Object_Identifying"]
jmults = sorted(["3Jets", "4PJets"])
leptypes = sorted(["Muon", "Electron"])


for hname in histograms:
    if hname == "nBestPerms":
        BestPerms_yields_and_effs_json = {}
        bp_rows = [(f"{args.year} Channel", "Likelihood Solutions", "Total", "Efficiency")]
        for jmult in jmults:
            for lepton in leptypes:
                passing_vals = hdict[f"{hname}_Solutions"][:, jmult, lepton].integrate("dataset").integrate("jmult").integrate("leptype").values()[()][0]
                total_vals = hdict[f"{hname}_Total"][:, jmult, lepton].integrate("dataset").integrate("jmult").integrate("leptype").values()[()][0]

                bp_rows += [(f"{jmult} {lepton}", format(passing_vals, ".1f"), format(total_vals, ".1f"), format(passing_vals/total_vals, ".4f"))]
                BestPerms_yields_and_effs_json[f"{lepton}_{jmult}"] = {"Passing" : passing_vals, "Total" : total_vals, "Eff" : passing_vals/total_vals}

        ## save best perms values
        bp_fname = os.path.join(outdir, f"BestPerms_Yields_and_Effs_{args.year}.txt")
        plt_tools.print_table(bp_rows, filename=bp_fname, print_output=True)
        print(f"{bp_fname} written")
        with open(bp_fname.replace(".txt", ".json"), "w") as out:
            out.write(prettyjson.dumps(BestPerms_yields_and_effs_json))
        print(f"{bp_fname.replace('.txt', '.json')} written")

    elif hname == "Object_Identifying":
        obj_rows = [(f"{args.year} Object", "Correctly ID", "Total", "Efficiency")]
        object_yields_and_effs_json = {}
        for jmult in jmults:
            for lepton in leptypes:
                object_yields_and_effs_json[f"{jmult}_{lepton}"] = {}
                obj_rows += [(jmult, lepton, "", "")]
                for decay_obj in sorted(set([key[0] for key in hdict[hname][:, jmult, lepton].integrate("dataset").integrate("jmult").integrate("leptype").values().keys()])):
                    hvals = hdict[hname][:, jmult, lepton, decay_obj].integrate("dataset").integrate("jmult").integrate("leptype").integrate("objtype").values()[()]
                    obj_rows += [(decay_obj, format(hvals[1], ".1f"), format(np.sum(hvals), ".1f"), format(hvals[1]/np.sum(hvals), ".4f"))]
                    object_yields_and_effs_json[f"{jmult}_{lepton}"].update({decay_obj : {"Passing" : hvals[1], "Total" : np.sum(hvals), "Eff" : hvals[1]/np.sum(hvals)}})
                obj_rows += [("", "", "", "")]

        #set_trace()
        ## save best perm object identification efficiencies
        obj_fname = os.path.join(outdir, f"BPobject_Identification_Yields_and_Effs_{args.year}.txt")
        plt_tools.print_table(obj_rows, filename=obj_fname, print_output=True)
        print(f"{obj_fname} written")
        with open(obj_fname.replace(".txt", ".json"), "w") as out:
            out.write(prettyjson.dumps(object_yields_and_effs_json))
        print(f"{obj_fname.replace('.txt', '.json')} written")

    elif hname == "ValidBestPerms":
        reco_cats_dict = {
            "Correcct" : 1,
            "Matchable" : 2,
            "Unmatchable" : 3,
            "SL_Tau" : 4,
        }
        #set_trace()
        ttcat_rows = [(f"{args.year} lj ttbar reconstruction categories", "", "", "")]
        ttcat_rows += [("Correct", "Matchable", "Unmatchable", "Tau")]
        ttcat_yields_and_effs_json = {}
        for jmult in jmults:
            for lepton in leptypes:
                ttcat_rows += [(f"{jmult} {lepton}", "", "", "")]
                vals = hdict[hname][:, jmult, lepton].integrate("dataset").integrate("jmult").integrate("leptype").values()[()]

                tot_val = np.sum(vals)
                effs = vals/tot_val

                ttcat_rows += [tuple(format(val, ".1f")+" ("+format(effs[idx], ".4f")+")" for idx, val in enumerate(vals) if idx > 0)]
                ttcat_yields_and_effs_json[f"{lepton}_{jmult}"] = {key: {"Yield"  : vals[idx] , "Eff" : effs[idx]} for key, idx in reco_cats_dict.items()}

        #set_trace()
        ## save reconstruction categories
        ttcat_fname = os.path.join(outdir, f"ttSL_Reconstruction_Categories_Yields_and_Effs_{args.year}.txt")
        plt_tools.print_table(ttcat_rows, filename=ttcat_fname, header_line=1, print_output=True)
        print(f"{ttcat_fname} written")
        with open(ttcat_fname.replace(".txt", ".json"), "w") as out:
            out.write(prettyjson.dumps(ttcat_yields_and_effs_json))
        print(f"{ttcat_fname.replace('.txt', '.json')} written")

    else:
        print(f"{hname} not supported")

toc = time.time()
print("Total time: %.1f" % (toc - tic))
