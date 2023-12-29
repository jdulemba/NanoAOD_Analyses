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

bkg_rname = os.path.join(outdir, f"final_templates_lj_bkg_nomSMTTxsec_{jobid}_{version}.root" if args.nomSMTTxsec else f"final_templates_lj_bkg_{jobid}_{version}.root")
if not os.path.isfile(bkg_rname): files_not_found.append(bkg_rname)

pdf_rname = os.path.join(outdir, f"final_pdf_templates_lj_bkg_{jobid}_V29.root")
if not os.path.isfile(pdf_rname): files_not_found.append(pdf_rname)

toponium_rname = os.path.join(eos_dir, "results", jobid, "Toponium_HBSR", "Templates", "V1", f"raw_templates_lj_toponium_{jobid}_V1.root")
if not os.path.isfile(toponium_rname): files_not_found.append(toponium_rname)

qcdorder_rname = os.path.join(outdir, f"NNLOqcd_LOqcd_templates_lj_{jobid}.root")
if not os.path.isfile(qcdorder_rname): files_not_found.append(qcdorder_rname)

if files_not_found:
    raise ValueError(f"{files_not_found} not found to be combined")

target_rname = os.path.join(outdir, f"final_BKG_PDF_ETAT_templates_lj_nomSMTTxsec_{jobid}_{version}.root" if args.nomSMTTxsec else f"final_BKG_PDF_ETAT_templates_lj_{jobid}_{version}.root")

hadd_cmd = f"hadd -f {target_rname} {bkg_rname} {pdf_rname} {toponium_rname} {qcdorder_rname}"
print(f"Executing: {hadd_cmd}")
set_trace()
os.system(hadd_cmd)

toc = time.time()
print("Total time: %.1f" % (toc - tic))
