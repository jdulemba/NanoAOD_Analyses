#! /bin/env python

import time
tic = time.time()

from coffea.util import load
from pdb import set_trace
import os
import numpy as np
import Utilities.systematics as systematics
import uproot
    
import argparse
class ParseKwargs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, dict())
        for value in values:
            key, value = value.split("=")
            getattr(namespace, self.dest)[key] = value

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("years_to_run", type=str, help="Choose which year to run, multiple options can be input as ':' separated strings.")
parser.add_argument("templates_to_run", type=str, help="Choose which type of templates to run, multiple options can be input as ':' separated strings.")
parser.add_argument("--unblind", action="store_true", help="Include unblinded data in templates, default is False.")
parser.add_argument("--nomSMTTxsec", action="store_true", help="Apply nominal SM cross sections to top mass and LHE scale weights")
args = parser.parse_args()

output_version = "V33"
input_version = "V33"
#output_version = "V32"
#input_version = "V31"
#version = "V31"
#version = "V30"
#version = "V29"
#version = "V27"
#version = "V26"

def final_bkg_templates(years_to_run):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """

    #set_trace()
    if len(years_to_run) == 4:
        outdir = global_outdir
        rname = os.path.join(outdir, f"final_templates_lj_bkg_nomSMTTxsec_{jobid}_{output_version}.root" if args.nomSMTTxsec else f"final_templates_lj_bkg_{jobid}_{output_version}.root")
        upfout = uproot.recreate(rname, compression=uproot.ZLIB(4)) if os.path.isfile(rname) else uproot.create(rname)

    for year in years_to_run:
        if len(years_to_run) < 4:
            outdir = os.path.join(eos_dir, "results", f"{year}_{jobid}", f"Templates_{analyzer}", "FINAL")
            rname = os.path.join(outdir, f"final_templates_lj_bkg_nomSMTTxsec_{year}_{jobid}_{output_version}.root" if args.nomSMTTxsec else f"final_templates_lj_bkg_{year}_{jobid}_{output_version}.root")
            upfout = uproot.recreate(rname, compression=uproot.ZLIB(4)) if os.path.isfile(rname) else uproot.create(rname)

        base_bkg_template_name = f"final_templates_lj_bkg_nomSMTTxsec_{year}_{jobid}_{input_version}.coffea" if args.nomSMTTxsec else f"final_templates_lj_bkg_{year}_{jobid}_{input_version}.coffea"
        input_dir = os.path.join(eos_dir, "results", f"{year}_{jobid}", f"Templates_{analyzer}", "FINAL")
        bkg_fname = os.path.join(input_dir, base_bkg_template_name)
        if not os.path.isfile(bkg_fname): raise ValueError("No background file found.")
        hdict = load(bkg_fname)

        for jmult in hdict.keys():
            for lep, histo in hdict[jmult].items():
                orig_lepdir = "muNJETS" if lep == "Muon" else "eNJETS"
                lepdir = orig_lepdir.replace("NJETS", jmult.lower())

                if year == "2016APV": year_to_use = "2016pre"
                elif year == "2016": year_to_use = "2016post"
                else: year_to_use = year
                dirname = f"{lepdir}_{year_to_use}"

                upfout.mkdir(dirname)

                systypes = ["nosys"] + sorted(systematics.final_systypes[year])
                for sys in systypes:
                    if (jobid == "Summer20UL_POG_lepSFs") or (jobid == "Summer20UL_DeepJet"):
                        if (lep == "Muon") and (("IDtot" in sys) or ("ISOtot" in sys) or ("TRIGtot" in sys)):
                            continue

                    if jobid == "Summer20UL_regroupedJECs":
                        if (sys == "nosys") or (sys == "JES_FlavorQCD") or (sys == "JES_RelativeBal"): continue

                        # find histograms of associated systematics and their processes
                    if sys == "nosys":
                        procs = sorted(set([key.split(f"_{sys}")[0] for key in histo.keys() if sys in key]))
                    else:
                        if (sys == "CR1") or (sys == "CR2") or (sys == "erdON"):
                            up_sysname = f"{sys}Up"
                            dw_sysname = f"{sys}Down"
                        else:
                            up_sysname = systematics.sys_groups[year][sys][0]
                            dw_sysname = systematics.sys_groups[year][sys][1]
                        procs = sorted(set([key.split(f"_{up_sysname}")[0] for key in sorted(histo.keys()) if f"_{up_sysname}" in key])) if not dw_sysname \
                            else sorted(set([key.split(f"_{up_sysname}")[0] for key in sorted(histo.keys()) if f"_{up_sysname}" in key] + [key.split(f"_{dw_sysname}")[0] for key in sorted(histo.keys()) if f"_{dw_sysname}" in key]))

                        ## make sure only 1 process is actually included in procs
                    procs = [proc for proc in procs if (len(proc.split("_")) == 1) or (proc == "data_obs")]

                    if not procs: continue

                    for proc in procs:
                        print(year, lep, jmult, sys, proc)
                        if sys == "nosys":
                            template, treatment = histo[f"{proc}_{sys}"]
                            if np.any(np.isnan(template.values()[()])):
                                print(f"\tSome bins contain nan values!")
                                set_trace()
                            if (proc == "data_obs") and (not args.unblind): template.clear()
                            upfout[dirname][proc] = template.to_hist()

                        else:
                            up_template, treatment = histo[f"{proc}_{up_sysname}"]
                            dw_template, treatment = histo[f"{proc}_{dw_sysname}"]

                            up_outhname = "_".join(list(filter(None, [proc, systematics.combine_template_sys_to_name[year][up_sysname]])))
                            if "LEP" in up_outhname: up_outhname = up_outhname.replace("LEP", lep[0].lower())
                            if "_tot" in up_outhname: up_outhname = up_outhname.replace("_tot", "")
                            dw_outhname = "_".join(list(filter(None, [proc, systematics.combine_template_sys_to_name[year][dw_sysname]])))
                            if "LEP" in dw_outhname: dw_outhname = dw_outhname.replace("LEP", lep[0].lower())
                            if "_tot" in dw_outhname: dw_outhname = dw_outhname.replace("_tot", "")
                            if (sys == "SHAPE") or (sys == "EWQCD_TTsub"):
                                up_outhname = up_outhname.replace("CHAN", lepdir)
                                dw_outhname = dw_outhname.replace("CHAN", lepdir)

                                # replace '2016APV' with '2016pre' and '2016' with '2016post'
                            if (year == "2016APV") and ("2016APV" in up_outhname):
                                #set_trace()
                                up_outhname = up_outhname.replace("2016APV", "2016pre")
                                dw_outhname = dw_outhname.replace("2016APV", "2016pre")
                            if (year == "2016") and ("2016" in up_outhname):
                                up_outhname = up_outhname.replace("2016", "2016post")
                                dw_outhname = dw_outhname.replace("2016", "2016post")

                                # replace 'ST' in single top ME/PS uncs with individual proc
                            if sys.startswith("ST_"):
                                up_outhname = up_outhname.replace("ST", proc)
                                dw_outhname = dw_outhname.replace("ST", proc)

                            if np.any(np.isnan(up_template.values()[()])):
                                print(f"\tSome bins in up_template contain nan values!")
                                set_trace()
                            if np.any(np.isnan(dw_template.values()[()])):
                                print(f"\tSome bins in dw_template contain nan values!")
                                set_trace()

                            upfout[dirname][up_outhname] = up_template.to_hist()
                            upfout[dirname][dw_outhname] = dw_template.to_hist()
                    
    upfout.close()
    print(f"{rname} written")


def final_pdf_templates(years_to_run):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """

    #set_trace()
    if len(years_to_run) == 4:
        outdir = global_outdir
        rname = os.path.join(outdir, f"final_pdf_templates_lj_bkg_{jobid}_{output_version}.root")
        upfout = uproot.recreate(rname, compression=uproot.ZLIB(4)) if os.path.isfile(rname) else uproot.create(rname)

    for year in years_to_run:
        if len(years_to_run) < 4:
            outdir = os.path.join(eos_dir, "results", f"{year}_{jobid}", f"Templates_{analyzer}", "FINAL")
            rname = os.path.join(outdir, f"final_pdf_templates_lj_bkg_{year}_{jobid}.root")
            upfout = uproot.recreate(rname, compression=uproot.ZLIB(4)) if os.path.isfile(rname) else uproot.create(rname)

        base_bkg_template_name = f"final_pdf_templates_lj_bkg_{year}_{jobid}_{input_version}.coffea"
        input_dir = os.path.join(eos_dir, "results", f"{year}_{jobid}", f"Templates_{analyzer}", "FINAL")
        bkg_fname = os.path.join(input_dir, base_bkg_template_name)
        if not os.path.isfile(bkg_fname): raise ValueError("No background file found.")
        hdict = load(bkg_fname)

        for jmult in hdict.keys():
            for lep, histo in hdict[jmult].items():
                orig_lepdir = "muNJETS" if lep == "Muon" else "eNJETS"
                lepdir = orig_lepdir.replace("NJETS", jmult.lower())

                if year == "2016APV": year_to_use = "2016pre"
                elif year == "2016": year_to_use = "2016post"
                else: year_to_use = year
                dirname = f"{lepdir}_{year_to_use}"
                upfout.mkdir(dirname)

                systs = sorted(set(["_".join(key.split("_")[1:]) for key in histo.keys() if not ("data_obs" in key or len(key.split("_")) == 1 or "shape" in key)]))
                for sys in systs:
                    if sys == "nosys": continue
                        # find histograms of associated systematics and their processes
                    procs = sorted(set([key.split(f"_{sys}")[0] for key in histo.keys() if sys in key]))
                    for proc in procs:
                        print(year, lep, jmult, sys, proc)
                        if sys == "nosys":
                            template, treatment = histo[f"{proc}_{sys}"]
                            upfout[dirname][proc] = template.to_hist()

                        else:
                            template, treatment = histo[f"{proc}_{sys}"]
                            if "alphaS" in sys:
                                outhname = "_".join([proc, "CMS", "PDF", sys.split("_")[-1]])
                            else:
                                pdf_idx = int(sys.split("pdf_")[-1])
                                outhname =  "_".join([proc, "PDF", f"{pdf_idx - 1}Up"])
                            upfout[dirname][outhname] = template.to_hist()

    upfout.close()
    print(f"{rname} written")


if __name__ == "__main__":
    allowed_template_options = ["bkg", "PDF"]
    templates_to_run = (args.templates_to_run).split(":")
    templates_to_run = [template for template in templates_to_run if template in allowed_template_options]

    allowed_year_options = ["2016APV", "2016", "2017", "2018"]
    years_to_run = (args.years_to_run).split(":")
    years_to_run = [year for year in years_to_run if year in allowed_year_options]

    jobid = os.environ["jobid"]
    eos_dir = os.environ["eos_dir"]

    global_outdir = os.path.join(eos_dir, "results", jobid, f"Templates_FINAL")
    if not os.path.isdir(global_outdir):
        os.makedirs(global_outdir)

    if "bkg" in templates_to_run:
        analyzer = "htt_btag_sb_regions"

        print("Creating final background templates")
        final_bkg_templates(years_to_run)

    if "PDF" in templates_to_run:
        analyzer = "htt_pdfUncs"

        print("Creating final background templates")
        final_pdf_templates(years_to_run)


    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
