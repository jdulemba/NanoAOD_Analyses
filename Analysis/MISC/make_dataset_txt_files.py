#! /bin/env python

"""
This script searches DBS for the datasets from the input json file and dumps the associated rootfiles into txt files

Created by Joseph Dulemba
30 October 2020
"""
from pdb import set_trace
import os, sys
import subprocess
from fnmatch import fnmatch
from pdb import set_trace
import Utilities.prettyjson as prettyjson
import Utilities.das as das
from site_mapping import site_name_to_xrootd

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("json", help="json file containing the samples definition")
parser.add_argument("--sample", help="Sample to run on, POSIX regex allowed")
parser.add_argument("--test", action="store_true", help="Output txt file named 'tmp.txt'")
parser.add_argument("--options", help="command-line arguments"
                    " to be passed to the configuration", default="")

args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]

if not os.path.isfile(args.json):
   raise ValueError(f"file {args.json} does not exist")

outdir = os.path.join(proj_dir, "inputs", "_".join(os.path.basename(args.json).split(".")[0].split("_")[1:])) # get name of json file except for "samples_"
if not os.path.isdir(outdir): os.makedirs(outdir)

all_samples = prettyjson.loads(open(args.json).read())
samples_to_run = list(filter(
   lambda x: fnmatch(x["name"], args.sample if args.sample else "*"),
     all_samples
))
if not len(samples_to_run):
    raise RuntimeError("Could not find any sample matching the pattern")

#set_trace()
analyzer_inputs = []
for sample in samples_to_run:
    #set_trace()
    
    if "DBSName" in sample:
        if sample["DBSName"] == "NOT PRESENT": continue
        if "Ext" in sample["name"]: print("Must combine %s with non-extenstion dataset!" % sample["name"])

        txtname = os.path.join(outdir, f"{sample['name']}_tmp.txt") if args.test else os.path.join(outdir, f"{sample['name']}.txt")

        flist = das.query(f"file dataset={sample['DBSName']} instance={sample['tier']}") if "tier" in sample else das.query(f"file dataset={sample['DBSName']}")#, True)
        #set_trace()
        for idx, fname in enumerate(flist):
            site_list = das.query(f"site file={fname}")
            already_changed = False
            for site in site_list:
                if already_changed: continue
                if site in site_name_to_xrootd.keys():
                    flist[idx] = fname.replace("/store", f"root://{site_name_to_xrootd[site]}//store")
                    #flist[idx] = fname.replace("/store", "root://%s//store" % site_name_to_xrootd[site])
                    already_changed = True

            if not already_changed:
                set_trace()
                print(f"No site found for {fname}")
                flist[idx] = fname.replace("/store", f"root://{site_name_to_xrootd['NEED']}//store")

        fnames = "\n".join(sorted(flist))

        txt_out = open(txtname, "w")
        txt_out.write(fnames)
        txt_out.close()
        print(f"{txtname} written")
        analyzer_inputs.append("\n%s" % sample["name"])
    else:
        raise ValueError(f"No DBSName found for sample: {sample['name']}")

if args.test: set_trace()
#if args.sample: set_trace()
analyzer_inputs_name = os.path.join(outdir, "analyzer_inputs.txt")
analyzer_inputs_out = open(analyzer_inputs_name, "a") if os.path.isfile(analyzer_inputs_name) else open(analyzer_inputs_name, "w")
analyzer_inputs_out.write("\n".join(analyzer_inputs))
analyzer_inputs_out.close()
print(f"{analyzer_inputs_name} written")
