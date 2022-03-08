#! /bin/env python

import time
tic = time.time()

from coffea.util import load
from pdb import set_trace
import os
import numpy as np
import Utilities.systematics as systematics
from coffea import hist
import uproot
    
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("--combine", action="store_true", help="Specify if output file is supposed to be for Combine or not.")
args = parser.parse_args()


def final_bkg_templates(bkg_dict):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """

    rname = os.path.join(input_dir, f"final_combine_pdf_templates_lj_bkg_{args.year}_{jobid}.root") if args.combine else os.path.join(input_dir, f"final_pdf_templates_lj_bkg_{args.year}_{jobid}.root")
    upfout = uproot.recreate(rname, compression=uproot.ZLIB(4)) if os.path.isfile(rname) else uproot.create(rname)

    #set_trace()
    for jmult, hdict in bkg_dict.items():
        for lep, histo in hdict.items():
            orig_lepdir = "muNJETS" if lep == "Muon" else "eNJETS"
            lepdir = orig_lepdir.replace("NJETS", jmult.lower())

            if args.combine:
                if args.year == "2016APV": year_to_use = "2016pre"
                elif args.year == "2016": year_to_use = "2016post"
                else: year_to_use = args.year
                dirname = f"{lepdir}_{year_to_use}"
            else:
                dirname = lepdir
            upfout.mkdir(dirname)

            #set_trace()
            systs = sorted(set(["_".join(key.split("_")[1:]) for key in histo.keys() if not ("data_obs" in key or len(key.split("_")) == 1 or "shape" in key)]))
            for sys in systs:
                if sys == "nosys": continue
                #set_trace()
                    # find histograms of associated systematics and their processes
                procs = sorted(set([key.split(f"_{sys}")[0] for key in histo.keys() if sys in key]))
                for proc in procs:
                    print(lep, jmult, sys, proc)
                    if sys == "nosys":
                        template, treatment = histo[f"{proc}_{sys}"]
                        upfout[dirname][proc] = template.to_hist()

                    else:
                        template, treatment = histo[f"{proc}_{sys}"]
                        outhname = "_".join([proc, "CMS", "PDF", sys.split("_")[-1]]) if args.combine else f"{proc}_{sys}"
                        upfout[dirname][outhname] = template.to_hist()

    upfout.close()
    print(f"{rname} written")


if __name__ == "__main__":
    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]

    analyzer = "htt_pdfUncs"
    base_bkg_template_name = f"final_pdf_templates_lj_bkg_{args.year}_{jobid}.coffea"

    input_dir = os.path.join(proj_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}", "FINAL")
    bkg_fname = os.path.join(input_dir, base_bkg_template_name)
    if not os.path.isfile(bkg_fname): raise ValueError("No background file found.")
    bkg_hdict = load(bkg_fname)

    print("Creating final background templates")
    final_bkg_templates(bkg_hdict)
    
    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
