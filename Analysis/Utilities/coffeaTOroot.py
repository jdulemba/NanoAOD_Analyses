#! /bin/env python

import time
tic = time.time()

from coffea.util import load
from pdb import set_trace
import os
from coffea import hist
import uproot
    
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("input_file", type=str, help="Input file to convert into a root file")
parser.add_argument("--output_fname", type=str, help="Name of output file with file extension, default is same name (and directory) as input but with '.root' extension")
args = parser.parse_args()


def walk_cdict(input):
    """
    Appends all keys from nested dictionary to 'path'
    except most nested key and then makes rootfile dir
    out of this path, with most nested key set as 
    histogram name.
    """
    global path
    for key, value in iter(input.items()):
        if isinstance(value, dict):
            path.append(key)
            walk_cdict(value)
            path.pop()
        elif isinstance(value, hist.Hist):
            dname = "_".join(path)
            upfout.mkdir(dname)
            upfout[dname][key] = value.to_hist()
            print(dname, key)
        else:
            set_trace()


if __name__ == "__main__":

    cfile = load(args.input_file)

    path = []
    rname = args.output_fname if args.output_fname else (args.input_file).replace(".coffea", ".root")
    upfout = uproot.recreate(rname, compression=uproot.ZLIB(4)) if os.path.isfile(rname) else uproot.create(rname)
    walk_cdict(cfile)

    upfout.close()
    print(f"{rname} written")
                    
    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
