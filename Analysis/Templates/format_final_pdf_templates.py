#! /bin/env python

import time
tic = time.time()

from pdb import set_trace
import os
from rootpy.io import root_open
    
base_jobid = os.environ["base_jobid"]
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
args = parser.parse_args()


def write_correct_template_format(in_fname, isSignal=None):
    """
    Opens temporary root file where template distributions are and then saves them with the correct structure/naming
    """
    if isSignal is None:
        raise ValueError("isSignal needs to be set to True to write signal templates, False for background")

    rfile = root_open(in_fname) if in_fname.endswith(".root") else root_open("%s.root" % in_fname)

    out_rfile = os.path.join(outdir, os.path.basename(in_fname).replace("temp_", "final_"))

    #set_trace()
    mu_3j_keys = sorted(set([key.name for key in rfile.keys() if "mu3jets" in key.name]))
    el_3j_keys = sorted(set([key.name for key in rfile.keys() if "e3jets" in key.name]))
    mu_4pj_keys = sorted(set([key.name for key in rfile.keys() if "mu4pjets" in key.name]))
    el_4pj_keys = sorted(set([key.name for key in rfile.keys() if "e4pjets" in key.name]))

    with root_open(out_rfile, "w") as rout:
            # mu, 3jets
        mu_3j_dir = rout.mkdir("mu3jets")
        for key in mu_3j_keys:
            hname = key.split("_mu3jets")[0]
            histo = rfile.Get(key)
            histo.name = hname
            if (hname == "data_obs") and (args.maskData):
                histo.Reset()
            mu_3j_dir.WriteTObject(histo, hname)

            # el, 3jets
        el_3j_dir = rout.mkdir("e3jets")
        for key in el_3j_keys:
            hname = key.split("_e3jets")[0]
            histo = rfile.Get(key)
            histo.name = hname
            if (hname == "data_obs") and (args.maskData):
                histo.Reset()
            el_3j_dir.WriteTObject(histo, hname)

            # mu, 4pjets
        mu_4pj_dir = rout.mkdir("mu4pjets")
        for key in mu_4pj_keys:
            hname = key.split("_mu4pjets")[0]
            histo = rfile.Get(key)
            histo.name = hname
            if (hname == "data_obs") and (args.maskData):
                histo.Reset()
            mu_4pj_dir.WriteTObject(histo, hname)

            # el, 4pjets
        el_4pj_dir = rout.mkdir("e4pjets")
        for key in el_4pj_keys:
            hname = key.split("_e4pjets")[0]
            histo = rfile.Get(key)
            histo.name = hname
            if (hname == "data_obs") and (args.maskData):
                histo.Reset()
            el_4pj_dir.WriteTObject(histo, hname)

    print(f"{out_rfile} written")



if __name__ == "__main__":
    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]

    analyzer = "htt_pdfUncs"
    #set_trace()
    outdir = os.path.join(proj_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}", "FINAL")
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    
    print("Creating final background templates")
    tmp_bkg_rname = os.path.join(outdir, f"temp_pdf_templates_lj_bkg_{args.year}_{jobid}.root")
    write_correct_template_format(tmp_bkg_rname, isSignal=False)

    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
