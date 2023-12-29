import time
tic = time.time()

"""
This script checks to see if any bins from the input root file have one-sided variations for the up and down uncertainties.
"""

import uproot
import numpy as np
import os
from pdb import set_trace

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("inputfile")
args = parser.parse_args()

if not os.path.isfile(args.inputfile):
    raise ValueError(f"{args.inputfile} not found")

rfile = uproot.open(args.inputfile)

dirnames = sorted(set([key.split("/")[0].split(";")[0] for key in rfile.keys()]))
mc_procs = sorted(set([key.split(";")[0] for key in rfile[dirnames[0]].keys() if (("data_obs" not in key) and ("Up" not in key) and ("Down" not in key))]))
for dirname in dirnames:
    print(f"Checking {dirname}")
    for proc in mc_procs:
        print(f"  {proc}")
        systypes = sorted(set([key.split(";")[0].replace(f"{proc}_", "").replace("Up", "").replace("Down", "")  for key in rfile[dirname].keys() if key.startswith(f"{proc}_")]))
        nom_vals = np.copy(rfile[dirname][proc].values())
        for sys in systypes:
            try:
                up_vals = np.copy(rfile[dirname][f"{proc}_{sys}Up"].values())
            except:
                up_vals = None
            try:
                dw_vals = np.copy(rfile[dirname][f"{proc}_{sys}Down"].values())
            except:
                dw_vals = None

                # make sure up and down variations actually exist
            if (up_vals is None) or (dw_vals is None): continue

                # check one-sidedness by looking at sign of value after multiplying relative values
            up_rel, dw_rel = np.copy(1. - up_vals/nom_vals), np.copy(1. - dw_vals/nom_vals)
            if np.any(np.nan_to_num(up_rel * dw_rel) > 0):
                #set_trace()
                print("    %s %s" % (sys, np.where(np.nan_to_num(up_rel * dw_rel) > 0)[0]))
    

toc = time.time()
print("Total time: %.1f" % (toc - tic))
