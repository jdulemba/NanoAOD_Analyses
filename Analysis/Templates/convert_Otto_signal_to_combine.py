#! /bin/env python

import time
tic = time.time()

from pdb import set_trace
import os, fnmatch
import Utilities.systematics as systematics
import uproot

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("sig_type", choices=["MC", "Folded", "Folded_LO"], help="Choose between signal produced via MC or Folding.")
parser.add_argument("--dir", type=str, default="*", help="Specify which directories to use.")
parser.add_argument("--test", action="store_true", help="Only choose specific signal to use for testing.")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
base_jobid = os.environ["base_jobid"]
jobid = os.environ["jobid"]
analyzer = "htt_btag_sb_regions"

output_dir = os.path.join(proj_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}", "FINAL")
if not os.path.isdir(output_dir):
    os.makedirs(output_dir)

template_version = "v12"
template_dname = "root://cmseos.fnal.gov//store/user/lpcbtagging/UR_ntuples/heavyhiggsinputs"
if args.sig_type == "MC":
    template_fname = os.path.join(template_dname, f"{template_version}_smoothed/templates_lj_sig_{args.year}.root")
elif args.sig_type == "Folded":
    template_fname = os.path.join(template_dname, f"{template_version}_folded/templates_lj_sig_{args.year}.root")
elif args.sig_type == "Folded_LO":
    template_fname = os.path.join(template_dname, f"{template_version}_folded_LO/templates_lj_sig_{args.year}.root")

rfile = uproot.open(template_fname)

#set_trace()
if args.test:
    out_rname = os.path.join(proj_dir, "tmp.root")
else:
    out_rname = os.path.join(output_dir, f"final_templates_lj_{template_version}_{args.sig_type.lower()}_sig_{args.year}_{jobid}_{args.dir}.root") if args.dir != "*" \
        else os.path.join(output_dir, f"final_templates_lj_{template_version}_{args.sig_type.lower()}_sig_{args.year}_{jobid}_ALL.root")
upfout = uproot.recreate(out_rname, compression=uproot.ZLIB(4)) if os.path.isfile(out_rname) else uproot.create(out_rname)

dir_choices = ["mu3jets", "mu4pjets", "e3jets", "e4pjets"]
dirs_to_use = fnmatch.filter(dir_choices, args.dir)
if not dirs_to_use:
    raise ValueError(f"Not match found for '{args.dir}' in {dir_choices}")

if args.year == "2016APV": year_to_use = "2016pre"
elif args.year == "2016": year_to_use = "2016post"
else: year_to_use = args.year

#set_trace()
valid_dists = [key for key in rfile.keys() if ( ("FINAL" not in key) and ("SMOOTH" not in key) and ("orig" not in key) and (len(key.split("/")) > 1) )] # get names of histograms
for key in valid_dists:
    #set_trace()
    dname, signal_hname = key.split("/")
    if dname not in dirs_to_use: continue
    if args.test:
        if "A365" not in signal_hname: continue
    if f"{dname}_{year_to_use};1" not in upfout.keys(): upfout.mkdir(f"{dname}_{year_to_use}")

    signal = signal_hname.strip(";1")
    if "res" in signal:
        sysname = signal.split("res_")[-1]
    elif "int_neg" in signal:
        sysname = signal.split("int_neg_")[-1]
    elif "int_pos" in signal:
        sysname = signal.split("int_pos_")[-1]
    else:
        raise ValueError(f"Treatment for {key} not found.")

    signame = signal.split(f"_{sysname}")[0]
    if signame == sysname: sysname = "nosys"

    print(dname, signame, sysname)
    combine_signame = signame.replace("relw", "w").replace("A", "A_").replace("H", "H_")
    mass = combine_signame.split("_")[1]
    combine_signame = combine_signame.replace(mass, "m"+mass)

    if sysname == "nosys":
        hname = combine_signame
    else:
        if ("muonsf" in sysname) or ("electronsf" in sysname):
            #set_trace()
            combine_sysname = systematics.combine_template_sys_to_name[args.year][[key for key, val in systematics.template_sys_to_name[args.year].items() \
                if sysname.replace("muon", "LEP").replace("electron", "LEP") == val][0]]
            combine_sysname = combine_sysname.replace("LEP", "m" if "muon" in sysname else "e")
        else:
            combine_sysname = systematics.combine_template_sys_to_name[args.year][[key for key, val in systematics.template_sys_to_name[args.year].items() if sysname == val][0]]

        combine_sysname = combine_sysname.replace(args.year, year_to_use)
        hname = f"{combine_signame}_{combine_sysname}"

    upfout[f"{dname}_{year_to_use}"][hname] = rfile[key].to_hist()

upfout.close()
print(f"{out_rname} written")

#set_trace()
toc = time.time()
print("Total time: %.1f" % (toc - tic))
