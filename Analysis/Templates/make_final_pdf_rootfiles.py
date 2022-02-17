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
            upfout.mkdir(lepdir)

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
                        upfout[lepdir][proc] = template.to_hist()

                    else:
                        template, treatment = histo[f"{proc}_{sys}"]
                        outhname = "_".join([proc, "CMS", "PDF", sys.split("_")[-1]]) if args.combine else f"{proc}_{sys}"
                        upfout[lepdir][outhname] = template.to_hist()
                
                        # add extra 'chi2' histogram for variations that were smoothed/flattened
                        if treatment != "raw":
                            set_trace()
                            #print("\t", sys, treatment)
                            #if treatment == "flat": set_trace()
                            treated_up_val = up_template.values(overflow="over")[()][-1]
                            treated_dw_val = dw_template.values(overflow="over")[()][-1]
            
                            chi2_outhname = "_".join(list(filter(None, [proc, systematics.new_template_sys_to_name[args.year][dw_sysname].split("Down")[0], "chi2", lepdir])))
                            if "LEP" in up_outhname: up_outhname = up_outhname.replace("LEP", lep[0].lower())

                            sumw = np.array([0., 0., 0., 0., treated_up_val, treated_dw_val])
                            tmp_chi2_histo = orig_chi2_histo.copy()
                                ## fill bins
                            for xbin in range(sumw.size):
                                tmp_chi2_histo.values()[()][xbin] = sumw[xbin]
                            
                            #set_trace()
                            upfout[lepdir][chi2_outhname] = tmp_chi2_histo.to_hist()

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

    orig_chi2_histo = hist.Hist("Events", hist.Bin("x_y", "x_y", np.arange(7)))
    orig_chi2_histo.fill(x_y=np.zeros(0))

    print("Creating final background templates")
    final_bkg_templates(bkg_hdict)
    
    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
