#! /bin/env python

import time
tic = time.time()
from pdb import set_trace
import os
import argparse
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--nomSMTTxsec", action="store_true", help="Apply nominal SM cross sections to top mass and LHE scale weights")
args = parser.parse_args()

eos_dir = os.environ["eos_dir"]
proj_dir = os.environ["PROJECT_DIR"]
if "srv" in proj_dir:
    raise ValueError("This can't be run within a singularity environment.")

jobid = os.environ["jobid"]

outdir = os.path.join(eos_dir, "results", jobid, "Templates_FINAL")
files_not_found = []

version = "V36"
#version = "V33"
#version = "V32"
#version = "V31"
#version = "V30"
#version = "V29"

bkg_rname = os.path.join(outdir, f"final_templates_lj_bkg_nomSMTTxsec_{jobid}_{version}.root" if args.nomSMTTxsec else f"final_templates_lj_bkg_{jobid}_{version}.root")
if not os.path.isfile(bkg_rname): files_not_found.appen(bkg_rname)

pdf_rname = os.path.join(outdir, f"final_pdf_templates_lj_bkg_{jobid}_V29.root")
#pdf_rname = os.path.join(outdir, f"final_pdf_templates_lj_bkg_{jobid}_{version}.root")
if not os.path.isfile(pdf_rname): files_not_found.appen(pdf_rname)

if files_not_found:
    raise ValueError(f"{files_not_found} not found to be combined")

target_rname = os.path.join(outdir, f"final_BKG_PDF_templates_lj_nomSMTTxsec_{jobid}_{version}.root" if args.nomSMTTxsec else f"final_BKG_PDF_templates_lj_{jobid}_{version}.root")

hadd_cmd = f"hadd -f {target_rname} {bkg_rname} {pdf_rname}"
print(f"Executing: {hadd_cmd}")
set_trace()
os.system(hadd_cmd)

toc = time.time()
print("Total time: %.1f" % (toc - tic))
