#!/usr/bin/env python

import os
from pdb import set_trace
import Run_Jobs.tools as tools

import argparse
class ParseKwargs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, dict())
        for value in values:
            key, value = value.split("=")
            getattr(namespace, self.dest)[key] = value

parser = argparse.ArgumentParser("submit analyzer to the batch queues")
parser.add_argument("analyzer", help="Analyzer to use.")
parser.add_argument("frange", type=str, help="Specify start:stop indices for files, inclusive. 0:1 means valid files in indices 0 and 1 will be used")
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="Specify which year to run over")
parser.add_argument("sample", type=str, help="Specify which sample to run over")
parser.add_argument("--opts", nargs="*", action=ParseKwargs, help="Options to pass to analyzers.")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
analyzer=args.analyzer

# define dictionary of options to pass
opts_dict = args.opts
if not opts_dict: opts_dict = {}
base_jobid = opts_dict.get("base_jobid", os.environ["base_jobid"])
jobid = opts_dict.get("jobid", os.environ["jobid"])

    ## get samples to use
#set_trace()
indir = os.path.join(proj_dir, "inputs", f"{args.year}_{base_jobid}")
samples_to_use = tools.get_sample_list(indir=indir, sample=args.sample)

fileset = {}
for sample in samples_to_use:
    if not os.path.isfile(sample):
        raise IOError(f"Sample file {sample} not found")

    sample_name = sample.split("/")[-1].split(".")[0]
    sfiles = open(sample, "r")
    valid_inds_fnames = [(idx, fname.strip("\n")) for idx, fname in enumerate(sfiles) if not fname.startswith("#")]
    valid_inds = [x[0] for x in valid_inds_fnames]
    valid_files = [x[1] for x in valid_inds_fnames]

    if ":" in args.frange:
        file_start, file_stop = int((args.frange).split(":")[0]), int((args.frange).split(":")[1])
    else:
        file_start = 0
        file_stop = len(valid_files)-1 if (args.frange).lower() == "all" else int(args.frange)-1

    if file_start >= 0 and file_stop < len(valid_files):
        valid_files = valid_files[file_start:file_stop+1]
        inds_to_use = valid_inds[file_start:file_stop+1]
    else:
        raise IOError(f"The number of root files available for the {sample} sample is %i. args.frange must be less than this." % len(valid_files) )

    fileset[sample_name] = valid_files

print(fileset)

    ## save output to coffea pkl file
#set_trace()
if (args.frange).lower() == "all":
    if "outfname" in opts_dict.keys():
        cfname = opts_dict["outfname"]
    elif args.sample:
        cfname = f"{args.sample}.coffea"
    else:
        cfname = f"test_{analyzer}.coffea"
else:
    if ":" in args.frange:
        if "outfname" in opts_dict.keys():
            cfname = opts_dict["outfname"]
        elif args.sample:
            cfname = f"{args.sample}_{file_start}to{file_stop}.coffea"
        else:
            cfname = f"test_{file_start}to{file_stop}.coffea"
    else:
        if "outfname" in opts_dict.keys():
            cfname = opts_dict["outfname"]
        elif args.sample:
            cfname = f"{args.sample}_{args.year}_{analyzer}.test.coffea"
        else:
            cfname = f"{args.year}_{analyzer}.test.coffea"


run_cmd = """python {PROJDIR}/bin/{ANALYZER}.py "{FSET}" {YEAR} {OUTFNAME} "{OPTS}" """.format(
        PROJDIR=proj_dir,
        ANALYZER=analyzer,
        FSET=fileset,
        YEAR=args.year,
        OUTFNAME=cfname,
        OPTS=opts_dict
)
print(f"Running command: {run_cmd}")
os.system(run_cmd)
