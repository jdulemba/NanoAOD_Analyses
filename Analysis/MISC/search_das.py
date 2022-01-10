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

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("dataset_name", help="Dataset string name to search for in DAS")
parser.add_argument("--status", default="VALID", help="Output txt file named 'tmp_dataset.txt'")
parser.add_argument("--txtname", default="tmp_dataset", help="Output txt file named 'tmp_dataset.txt'")

args = parser.parse_args()

jobid = os.environ["jobid"]
proj_dir = os.environ["PROJECT_DIR"]

dlist = das.query(f"dataset status={args.status} dataset={args.dataset_name}")
dnames = "\n".join(dlist)

txtname = os.path.join(proj_dir, f"{args.txtname}.txt")
if os.path.isfile(txtname):
   raise ValueError(f"File {txtname} already exists!")

#set_trace()
txt_out = open(txtname, "w")
txt_out.write(dnames)
txt_out.close()
print(f"{txtname} written")
