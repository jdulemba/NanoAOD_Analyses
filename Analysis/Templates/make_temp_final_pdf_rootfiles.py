#! /bin/env python

import time
tic = time.time()

from coffea.util import load
from pdb import set_trace
import os
import numpy as np
import fnmatch
import Utilities.systematics as systematics
import uproot3
from coffea import hist
    
base_jobid = os.environ["base_jobid"]
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
args = parser.parse_args()


def final_bkg_templates(bkg_dict):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """

    tmp_rname = os.path.join(outdir, f"temp_pdf_templates_lj_bkg_{args.year}_{jobid}.root")
    upfout = uproot3.recreate(tmp_rname, compression=uproot3.ZLIB(4)) if os.path.isfile(tmp_rname) else uproot3.create(tmp_rname)

    #set_trace()
    for jmult, hdict in bkg_dict.items():
        for lep, histo in hdict.items():
            orig_lepdir = "muNJETS" if lep == "Muon" else "eNJETS"
            lepdir = orig_lepdir.replace("NJETS", jmult.lower())

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
                        upfout[f"{proc}_{lepdir}"] = hist.export1d(template.copy())

                    else:
                        template, treatment = histo[f"{proc}_{sys}"]
                        upfout[f"{proc}_{sys}_{lepdir}"] = hist.export1d(template.copy())
                
                        # add extra 'chi2' histogram for variations that were smoothed/flattened
                        if treatment != "raw":
                            set_trace()
                            #print("\t", sys, treatment)
                            #if treatment == "flat": set_trace()
                            treated_up_val = up_template.values(overflow="over")[()][-1]
                            treated_dw_val = dw_template.values(overflow="over")[()][-1]
            
                            chi2_outhname = "_".join(list(filter(None, [proc, systematics.template_sys_to_name[args.year][dw_sysname].split("Down")[0], "chi2", lepdir])))
                            if "LEP" in chi2_outhname: chi2_outhname = chi2_outhname.replace("LEP", "muon") if lep == "Muon" else chi2_outhname.replace("LEP", "electron")

                            sumw = np.array([0., 0., 0., 0., treated_up_val, treated_dw_val])
                            tmp_chi2_histo = orig_chi2_histo.copy()
                                ## fill bins
                            for xbin in range(sumw.size):
                                tmp_chi2_histo.values()[()][xbin] = sumw[xbin]
                            
                            #set_trace()
                            upfout[chi2_outhname] = hist.export1d(tmp_chi2_histo.copy())

    upfout.close()
    print(f"{tmp_rname} written")



if __name__ == "__main__":
    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]

    analyzer = "htt_pdfUncs"
    base_bkg_template_name = f"final_pdf_templates_lj_NJETS_bkg_{args.year}_{jobid}"

    #set_trace()
    input_dir = os.path.join(proj_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}", "FINAL")
    outdir = os.path.join(proj_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}", "FINAL")
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    
    if os.path.isdir(input_dir):
            # find files for 3 jets
        bkg_3j_fnmatch = "%s.coffea" % base_bkg_template_name.replace("NJETS", "3Jets")
        bkg_3j_fnames = fnmatch.filter(os.listdir(input_dir), bkg_3j_fnmatch)
        bkg_3j_fnames = [os.path.join(input_dir, fname) for fname in bkg_3j_fnames]
        if len(bkg_3j_fnames) > 1: raise ValueError("More than one input file found for 3 jets templates.")
            # find files for 4+ jets
        bkg_4pj_fnmatch = "%s.coffea" % base_bkg_template_name.replace("NJETS", "4PJets")
        bkg_4pj_fnames = fnmatch.filter(os.listdir(input_dir), bkg_4pj_fnmatch)
        bkg_4pj_fnames = [os.path.join(input_dir, fname) for fname in bkg_4pj_fnames]
        if len(bkg_4pj_fnames) > 1: raise ValueError("More than one input file found for 4+ jets templates.")
    
        bkg_dict = {"3Jets" : load(bkg_3j_fnames[0]), "4PJets" : load(bkg_4pj_fnames[0])}
    else: print("No background file found.")
    
    orig_chi2_histo = hist.Hist("Events", hist.Bin("x_y", "x_y", np.arange(7)))
    orig_chi2_histo.fill(x_y=np.zeros(0))
    
    baseSys = lambda sys : "_".join(sys.split("_")[:-1])
    
    print("Creating final background templates")
    final_bkg_templates(bkg_dict)

    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
