#! /bin/env python

import time
tic = time.time()

from coffea.util import load
from pdb import set_trace
import os
import Utilities.plot_tools as plt_tools
import Utilities.prettyjson as prettyjson
import numpy as np
import fnmatch

base_jobid = os.environ["base_jobid"]

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2018"], help="What year is the ntuple from.")
#parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
analyzer = "jet_puID_effs"
eos_dir = os.environ["eos_dir"]
plots_dir = os.environ["plots_dir"]

base_ext = "*TOT.coffea"
input_dir = os.path.join(eos_dir, "results", f"{args.year}_{jobid}", analyzer)
fnames = fnmatch.filter(os.listdir(input_dir), base_ext)
fnames = [os.path.join(input_dir, fname) for fname in fnames]
hdict = load(fnames[0])

outdir = os.path.join(plots_dir, jobid, analyzer, args.year)
if not os.path.isdir(outdir):
    os.makedirs(outdir)


jet_bin_mapping = ["hiPt", "lowPt", "gen match", "loosePU", "mediumPU", "tightPU"]

rows = [("Channel", "Total tightID jets", "njets pT >= 50", "njets pT < 50", "gen match", "loosePU", "mediumPU", "tightPU")]

    # individual channels
for lep in ["Muon", "Electron"]:
    for jmult in ["3Jets", "4PJets"]:
        vals = hdict["NJets"].integrate("dataset").values()[(jmult, lep)]
        yields_dict = {key : vals[idx] for idx, key in enumerate(jet_bin_mapping)}
        yields_dict.update({"tot_njets" : np.sum(vals[:2])})
        rows += [( f"{lep}/{jmult}", f"{yields_dict['tot_njets']}",
            "%.1f (%.3f)" % (yields_dict["hiPt"], yields_dict["hiPt"]/yields_dict["tot_njets"]),
            "%.1f (%.3f)" % (yields_dict["lowPt"], yields_dict["lowPt"]/yields_dict["tot_njets"]),
            "%.1f (%.3f)" % (yields_dict["gen match"], (yields_dict["gen match"]+yields_dict["hiPt"])/yields_dict["tot_njets"]),
            "%.1f (%.3f)" % (yields_dict["loosePU"], (yields_dict["loosePU"]+yields_dict["hiPt"])/yields_dict["tot_njets"]),
            "%.1f (%.3f)" % (yields_dict["mediumPU"], (yields_dict["mediumPU"]+yields_dict["hiPt"])/yields_dict["tot_njets"]),
            "%.1f (%.3f)" % (yields_dict["tightPU"], (yields_dict["tightPU"]+yields_dict["hiPt"])/yields_dict["tot_njets"])
        )]
        rows += [("", "", "", "", "", "", "", "")]

#set_trace()
    # only jmult channels (integrate lepton)
for jmult in ["3Jets", "4PJets"]:
    vals = hdict["NJets"].integrate("dataset").integrate("leptype").values()[(jmult,)]
    yields_dict = {key : vals[idx] for idx, key in enumerate(jet_bin_mapping)}
    yields_dict.update({"tot_njets" : np.sum(vals[:2])})
    rows += [( jmult, f"{yields_dict['tot_njets']}",
        "%.1f (%.3f)" % (yields_dict["hiPt"], yields_dict["hiPt"]/yields_dict["tot_njets"]),
        "%.1f (%.3f)" % (yields_dict["lowPt"], yields_dict["lowPt"]/yields_dict["tot_njets"]),
        "%.1f (%.3f)" % (yields_dict["gen match"], (yields_dict["gen match"]+yields_dict["hiPt"])/yields_dict["tot_njets"]),
        "%.1f (%.3f)" % (yields_dict["loosePU"], (yields_dict["loosePU"]+yields_dict["hiPt"])/yields_dict["tot_njets"]),
        "%.1f (%.3f)" % (yields_dict["mediumPU"], (yields_dict["mediumPU"]+yields_dict["hiPt"])/yields_dict["tot_njets"]),
        "%.1f (%.3f)" % (yields_dict["tightPU"], (yields_dict["tightPU"]+yields_dict["hiPt"])/yields_dict["tot_njets"])
    )]
    rows += [("", "", "", "", "", "", "", "")]


frac_name = os.path.join(outdir, "TTbar_NJets_PU_yields_and_fracs.txt")
plt_tools.print_table(rows, filename=frac_name, print_output=True)#, header_line=1)
print(f"{frac_name} written")

toc = time.time()
print("Total time: %.1f" % (toc - tic))
