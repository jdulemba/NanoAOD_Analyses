#! /bin/env python

import time
tic = time.time()

from pdb import set_trace
import os
import numpy as np
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
parser.add_argument("input_dir", type=str, help="Specify the directory name from which to choose the signal files.")
parser.add_argument("--MEopts", nargs="*", action=ParseKwargs, help="Options to specify which mass and width points to produce for ME reweighted signal.")
args = parser.parse_args()

if __name__ == "__main__":
    possible_masses = [365, 380, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000]
    possible_widths = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 8.0, 10.0, 13.0, 15.0, 18.0, 21.0, 25.0]
    possible_masses = [str(mass) for mass in possible_masses]
    possible_widths = [str(width) for width in possible_widths]

    masses_to_run = "All" if args.MEopts is None else args.MEopts.get("allowed_masses", "All")
    widths_to_run = "All" if args.MEopts is None else args.MEopts.get("allowed_widths", "All")

    if masses_to_run == "All":
        allowed_masses = possible_masses
    else:
        allowed_masses = masses_to_run.split(":")
        allowed_masses = [mass for mass in allowed_masses if mass in possible_masses]

    if widths_to_run == "All":
        allowed_widths = possible_widths
    else:
        allowed_widths = widths_to_run.split(":")
        allowed_widths = [width for width in allowed_widths if width in possible_widths]

    allowed_masses = ["m"+mass for mass in allowed_masses]
    allowed_widths = ["w"+width.replace(".", "p") for width in allowed_widths]
    print(f"Making target signal points for {allowed_masses} and {allowed_widths}")

    input_dir = f"root://cmseos.fnal.gov//store/user/jdulemba/Htt_Templates/{args.input_dir}/FINAL"
    rnames = [
        os.path.join(input_dir, "templates_lj_sig_2016pre.root"),
        os.path.join(input_dir, "templates_lj_sig_2016post.root"),
        os.path.join(input_dir, "templates_lj_sig_2017.root"),
        os.path.join(input_dir, "templates_lj_sig_2018.root")
    ]

    for mass in allowed_masses:
        print(f"Mass point: {mass}")

            ## make name of output root file
        rname = f"templates_lj_sig_{mass}.root" if (widths_to_run == "All") else f"templates_lj_sig_{mass}_{''.join(allowed_widths).lower()}.root"
        upfout = uproot.recreate(rname, compression=uproot.ZLIB(4)) if os.path.isfile(rname) else uproot.create(rname)

        for rfname in rnames:
            try:
                rfile = uproot.open(rfname)
            except:
                raise ValueError(f"Could not open {rfname}")

            dirnames = sorted(rfile.keys(recursive=False, cycle=False))
            for dirname in dirnames:
                print(f"\tCopying {dirname}")
                upfout.mkdir(dirname)
                #set_trace()
                for width in allowed_widths:
                    print(f"\t\t{width}")
                    upfout[dirname].copy_from(rfile[dirname], filter_name=f"*{mass}*{width}*", require_matches=True)

            rfile.close()

        upfout.close()
        print(f"{rname} written")

            ## copy output file to input dir and then delete local copy
        exec_str = f"xrdcp {rname} {input_dir} && rm {rname}"
        try:
            print(f"Executing {exec_str}")
            os.system(exec_str)
        except:
            raise ValueError(f"Could not copy {rname} to {input_dir}")

    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
