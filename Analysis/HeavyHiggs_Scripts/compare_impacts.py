#! /bin/env python

import time
tic = time.time()

print("Do not run this within singularity")

import json
import os
import numpy as np
from pdb import set_trace
from tqdm import tqdm
import Utilities.plot_tools as plt_tools

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("parity", choices=["A", "H"], help="Choose which signal parity")
parser.add_argument("mass", choices=["400", "800"], help="Choose which signal mass")
parser.add_argument("--param", default="tmass_3GeV_TT", help="Choose which parameter to print out")
#parser.add_argument("width", choices=["5"], help="Choose which signal width")
args = parser.parse_args()

eos_dir = os.environ["eos_dir"]
jobid = os.environ["jobid"]

#input_dir = os.path.join(eos_dir, "results", jobid, "ahtt_run2ul", "misc", "impact_smtt_shapeu_230121") 
input_dir = os.path.join(eos_dir, "results", jobid, "ahtt_run2ul", "misc", "impact_smtt_shapeu_230124") 

width = "5p0" # args.width
input_jsons = [os.path.join(input_dir, fname) for fname in os.listdir(input_dir) if ((f"{args.parity}_m{args.mass}_w{width}" in fname) and fname.endswith(".json") )]
#set_trace()

    # configurations/naming based on parameter choices
str_to_float = lambda x : x.replace("p", ".")
multipliers = {"Gaussian" : 3000, "Unconstrained" : 3000}
#multipliers = {"Gaussian" : 500, "Unconstrained" : 3000}
channels = {"lx" : "lj+ll", "lj" : "lj", "ll" : "ll"}
prior_type = {"g" : "Gaussian", "f" : "flat"}
def mt_shift(mt_str):
    shift_val = mt_str.split("mt")[-1]
    if shift_val == "p0":
        ret_str = "0"
    else:
        ret_str = shift_val.replace("m", "-").replace("p", "+")
    return ret_str


rows = [(f"{args.parity} {args.mass} {width}", args.param, "", "", "", "", "", ""),
    ("channel", "prior type", "mt fit range (GeV)", "Asimov g", "fitted g", "mt Asimov shift (MeV)", "pull (MeV)", "impact")]

for idx in tqdm(range(len(input_jsons))):
    fname = os.path.basename(input_jsons[idx]).strip(".json")
    if len(fname.split("_"))  != 11: continue #set_trace()

    parity, mass, width, channel_prior_mtrange, asimov_g, mt_shift_str, _, _, _, _, _ = fname.split("_")
    channel, prior, mtrange = channel_prior_mtrange[:2], channel_prior_mtrange[2], channel_prior_mtrange[-1]

    fit_json = json.load(open(input_jsons[idx]))

        # get overall fitted value of POI
    gvals = [g_dict for idx, g_dict in enumerate(fit_json["POIs"]) if g_dict["name"] == "r"][0]["fit"]

        # get impact and constraints on wanted parameter
    #set_trace()
    param_dict_list = [p_dict for idx, p_dict in enumerate(fit_json["params"]) if p_dict["name"] == args.param]
    if len(param_dict_list) != 1:
        print(os.path.basename(input_jsons[idx]))
        #set_trace()
        rows += [(channels[channel], prior_type[prior], mtrange, str_to_float(asimov_g[1:]), "-", mt_shift(mt_shift_str), "-", "-")]
        continue

    param_dict = param_dict_list[0]

    c_vals = param_dict["fit"] # get constraint vals
    c_multiplier = multipliers[param_dict["type"]] # get multiplier for constraint to convert value into MeV
    impact = param_dict["impact_r"]
    
    #set_trace()
    rows += [(channels[channel], prior_type[prior], mtrange, str_to_float(asimov_g[1:]),
            format(gvals[1], ".3f")+"  +"+format(gvals[2]-gvals[1], "0.3f")+"/"+format(gvals[0]-gvals[1], "0.3f"),
            mt_shift(mt_shift_str), format(c_multiplier*c_vals[1], "0.1f")+"  +"+format(c_multiplier*(c_vals[2]-c_vals[1]), "0.1f")+"/"+format(c_multiplier*(c_vals[0]-c_vals[1]), "0.1f"),
             format(impact, ".4f"))
    ]
    

#set_trace()
frac_name = os.path.join(input_dir, f"{args.parity}_m{args.mass}_w{width}_{args.param}_impacts_comp.txt")
plt_tools.print_table(rows, filename=frac_name, print_output=True, header_line=1)
print(f"{frac_name} written")

toc = time.time()
print("Total time: %.1f" % (toc - tic))
