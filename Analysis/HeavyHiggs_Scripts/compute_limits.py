#! /bin/env python

import time
tic = time.time()

# matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend("agg")
from matplotlib import rcParams
rcParams["font.size"] = 20
rcParams["savefig.format"] = "png"
rcParams["savefig.bbox"] = "tight"
from coffea.util import load
from pdb import set_trace
import os
import numpy as np
import Utilities.Plotter as Plotter
import Utilities.systematics as systematics
from Utilities.styles import styles as hstyles
import uproot
import fnmatch
from scipy.interpolate import interp1d
from coffea.util import save

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("type", choices=["MC", "MEreweighting", "Morphed"], help="Choose between signal produced via MC, ME reweighting, or Folding.")
#parser.add_argument("type", choices=["MC", "Folded_LO", "MEreweighting_LO", "Morphed"], help="Choose between signal produced via MC, ME reweighting, or Folding.")
parser.add_argument("--widths", type=str, default="All", help="Choose which width points to compute limits for, multiple options can be input as ':' separated strings.")
parser.add_argument("--parity", type=str, default="All", help="Choose which parity to compute limits for, multiple options can be input as ':' separated strings.")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
base_jobid = os.environ["base_jobid"]
jobid = os.environ["jobid"]
analyzer = "htt_btag_sb_regions"

eos_input_dir = "root://cmseos.fnal.gov//store/user/jdulemba/Limit_Results"


if args.type == "MC":
    set_trace()
    output_dir = os.path.join("/eos", "user", os.environ["USER"][0], os.environ["USER"], "NanoAOD_Analyses", "plots", f"{args.year}_{jobid}", f"Templates_{analyzer}", "MC_SIG", template_version, args.lepton)
    template_fname = os.path.join(template_dname, f"{template_version}_mcsig_smoothed/templates_lj_sig.root") if template_version == "v13" else os.path.join(template_dname, f"{template_version}_smoothed/templates_lj_sig_{args.year}.root")

    rfile = uproot.open(template_fname)

elif args.type == "Folded":
    set_trace()
    output_dir = os.path.join("/eos", "user", os.environ["USER"][0], os.environ["USER"], "NanoAOD_Analyses", "plots", f"{args.year}_{jobid}", f"Templates_{analyzer}", "FOLDED_LO_SIG", template_version, args.lepton)
    template_fname = os.path.join(template_dname, f"{template_version}_folded_LO/templates_lj_sig_{args.year}.root")
    rfile = uproot.open(template_fname)

elif args.type == "MEreweighting":
    #set_trace()

    version = "V21"
       # All uncs
    base_lim_dirname = os.path.join(eos_input_dir, version, "ASYMv21_lj_nosig_1")
    possible_widths = [1.0, 2.0, 2.5, 3.0, 4.0, 5.0, 8.0, 10.0, 13.0, 15.0, 18.0, 21.0, 25.0]
    #   # All uncs
    #base_lim_dirname = os.path.join(eos_input_dir, "ASYMv13mewsmoothed_nosig_2")
    #possible_widths = [2.5, 3.0, 4.0, 5.0, 8.0, 10.0, 13.0, 15.0, 18.0, 21.0, 25.0]
    #   # Some uncs
    #base_lim_dirname = os.path.join(eos_input_dir, "ASYMv14mewsmoothed_nosig_4uncpart")
    #possible_widths = [2.5, 10.0, 15.0, 21.0, 25.0]
    #   # No uncs
    #base_lim_dirname = os.path.join(eos_input_dir, "ASYMv13mewsmoothed_nosig_3nouncs")
    #possible_widths = [2.5, 10.0, 15.0, 25.0]

    output_dir = os.path.join("/eos", "user", os.environ["USER"][0], os.environ["USER"], "NanoAOD_Analyses", "results", jobid, f"Limits_{analyzer}", version, os.path.basename(base_lim_dirname))
    #output_dir = os.path.join("/eos", "user", os.environ["USER"][0], os.environ["USER"], "NanoAOD_Analyses", "results", jobid, f"Limits_{analyzer}", os.path.basename(base_lim_dirname))
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    available_masses = [365, 380, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000]
    available_masses = [str(mass) for mass in available_masses]
    possible_widths = [str(width) for width in possible_widths]

elif args.type == "Morphed":
    #set_trace()
    #   # All uncs
    #base_lim_dirname = os.path.join(eos_input_dir, "ASYMv13mcsigmorphedsmoothed_nosig_4")
    #possible_widths = [2.5, 5.0, 10.0, 15.0, 20.0, 25.0]
    #   # Some uncs
    #base_lim_dirname = os.path.join(eos_input_dir, "ASYMv13mcsigmorphedsmoothed_nosig_4uncpart")
    #possible_widths = [2.5, 5.0, 10.0, 15.0, 21.0, 25.0]
       # No uncs
    base_lim_dirname = os.path.join(eos_input_dir, "ASYMv13mcsigsmooth_nosig_3nouncs")
    possible_widths = [2.5, 10.0, 15.0, 25.0]

    output_dir = os.path.join("/eos", "user", os.environ["USER"][0], os.environ["USER"], "NanoAOD_Analyses", "results", jobid, f"Limits_{analyzer}", os.path.basename(base_lim_dirname))
    if not os.path.isdir(output_dir):
        os.makedirs(output_dir)

    available_masses = sorted(np.arange(370., 1010., 10., dtype=int))
    available_masses = [str(mass) for mass in available_masses]
    possible_widths = [str(width) for width in possible_widths]


#set_trace()
widths_to_run = args.widths
if widths_to_run == "All":
    available_widths = possible_widths
else:
    available_widths = widths_to_run.split(":")
    available_widths = [width for width in available_widths if width in possible_widths]

possible_parities = ["A", "H"]
parities_to_run = args.parity
if parities_to_run == "All":
    available_parities = possible_parities
else:
    available_parities = parities_to_run.split(":")
    available_parities = [parity for parity in available_parities if parity in possible_parities]

#set_trace()
for boson in available_parities:
    for relw in available_widths:
        lim_dirname = f"{base_lim_dirname}_{boson}_relw{relw.replace('.', 'p')}"

        limit_results_dict = {}
        for mass_idx in range(len(available_masses)):
            fname = os.path.join(lim_dirname, f"fitresults_{mass_idx}_{len(available_masses)}.root")
            try:
                rfile = uproot.open(fname)
            except:
                print(f"{fname} not found")
                continue
                #raise ValueError(f"{fname} not found")
            #rfile = uproot.open(fname)

            print(f"Finding limits for {boson} {relw} {available_masses[mass_idx]}")
            g_crossing_per_mass = []
            for cls_ind in sorted([str(ind) for ind in np.arange(5)]):
                cls_vals = rfile[f"CLs_{cls_ind}"].values()
                g_vals = rfile[f"CLs_{cls_ind}"].to_hist().axes[0].edges

                    # find index where CLs values cross 0.95
                crossings = np.where(np.diff(np.sign(cls_vals - 0.95)))[0] # finds index right before crossing
                if crossings.size == 0:
                    print(f"No CLs values cross 0.95 for {fname}, need to investigate...")
                    continue

                elif crossings.size == 1:
                        # get points before and after crossing to be used in interpolation
                    cls_cross = cls_vals[crossings[0]:crossings[0]+2]
                    g_cross = g_vals[crossings[0]:crossings[0]+2]
                        # calculate interpolation and estimated new CLs values between the known points
                    f_interp = interp1d(g_cross, cls_cross)
                    g_new = np.linspace(g_cross[0], g_cross[1], 101)
                    cls_new = f_interp(g_new)
                        # find where closest estimated point to g = 0.05 is
                    estimated_g_crossing = g_new[np.where(np.abs(cls_new-0.95) == min(np.abs(cls_new-0.95)))[0]][0]

                    g_crossing_per_mass.append(estimated_g_crossing)

                else:
                    print(f"There are {crossings.size} crossings for {fname}, need to investigate...")
                    #set_trace()
                    # use first crossing for -1/-2 sigma (inds 3 and 4)
                    if int(cls_ind) > 2:
                            # get points before and after crossing to be used in interpolation
                        cls_cross = cls_vals[crossings[0]:crossings[0]+2]
                        g_cross = g_vals[crossings[0]:crossings[0]+2]
                            # calculate interpolation and estimated new CLs values between the known points
                        f_interp = interp1d(g_cross, cls_cross)
                        g_new = np.linspace(g_cross[0], g_cross[1], 101)
                        cls_new = f_interp(g_new)
                            # find where closest estimated point to g = 0.05 is
                        estimated_g_crossing = g_new[np.where(np.abs(cls_new-0.95) == min(np.abs(cls_new-0.95)))[0]][0]

                        g_crossing_per_mass.append(estimated_g_crossing)

                    # use last crossing for expected, +1/+2 sigma (inds 0, 1, 2)
                    else:
                            # get points before and after crossing to be used in interpolation
                        cls_cross = cls_vals[crossings[-1]:crossings[-1]+2]
                        g_cross = g_vals[crossings[-1]:crossings[-1]+2]
                            # calculate interpolation and estimated new CLs values between the known points
                        f_interp = interp1d(g_cross, cls_cross)
                        g_new = np.linspace(g_cross[0], g_cross[1], 101)
                        cls_new = f_interp(g_new)
                            # find where closest estimated point to g = 0.05 is
                        estimated_g_crossing = g_new[np.where(np.abs(cls_new-0.95) == min(np.abs(cls_new-0.95)))[0]][0]

                        g_crossing_per_mass.append(estimated_g_crossing)

                    #continue
            print(f"\t{np.array(g_crossing_per_mass)}")
            limit_results_dict[str(available_masses[mass_idx])] = np.array(g_crossing_per_mass)

        #set_trace()
        out_fname = os.path.join(output_dir, f"{boson}_relw{relw.replace('.', 'p')}.coffea")
        save(limit_results_dict, out_fname)
        print(f"{out_fname} written")


toc = time.time()
print("Total time: %.1f" % (toc - tic))
