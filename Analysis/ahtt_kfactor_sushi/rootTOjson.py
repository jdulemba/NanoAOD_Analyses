#!/usr/bin/env python

from pdb import set_trace
import os
import numpy as np
import Utilities.prettyjson as prettyjson
from coffea.lookup_tools.root_converters import convert_histo_root_file
from coffea.lookup_tools.dense_lookup import dense_lookup

proj_dir = os.environ["PROJECT_DIR"]

rname = os.path.join(proj_dir, "ahtt_kfactor_sushi", "Final_ULkfactor_03Sep2021.root")
if not os.path.isfile(rname): raise ValueError(f"{rname} not found")

bosons = ["A", "H"]
masses = [365, 380, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000] # available mass points (GeV)
widths = [0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 8, 10, 13, 15, 18, 21, 25] # available widths (%)
scales = ["nominal", "uF_up", "uF_down", "uR_up", "uR_down", "uF_up_uR_up", "uF_down_uR_down"] # available scale
channels = {"ll" : "DiLep", "lj" : "SL"}
procs = {"res" : "Res", "int" : "Int"}
types = ["xsec", "xabs"] # <type> can be either xsec for actual cross section or xabs for the absolute cross section i.e. sum of the magnitude of positive and negative parts of the cross section

pos_evt_fraction_name = "int_mg5_pdf_325500_scale_dyn_0p5mtt_positive_event_fraction" # <A/H>_int_mg5_pdf_325500_scale_dyn_0p5mtt_positive_event_fraction contain the fraction of positive events for each signal points
kfactors_name = "sushi_nnlo_mg5_lo_kfactor_pdf_325500" # <A/H>_<res/int>_sushi_nnlo_mg5_lo_kfactor_pdf_325500_<scale> contain the NNLO to LO k-factors
mg5_LO_xsecs_name = "mg5_pdf_325500_scale_dyn_0p5mtt" # <A/H>_<res/int>_mg5_pdf_325500_scale_dyn_0p5mtt_<scale>_<type>_<channel> contain the MG5 LO cross section for the relevant process and channel

fdict = convert_histo_root_file(rname)

widthTOname = lambda width : str(float(width)).replace('.', 'p')

# create dict of the fraction of positive events for each signal points
#outdict = {}
#for boson in bosons:
#    sig_dname = "_".join([boson, pos_evt_fraction_name]) # name of signal dist
#    vals = dense_lookup(*fdict[(sig_dname, "dense_lookup")])
#    errs = dense_lookup(*fdict[(f"{sig_dname}_error", "dense_lookup")])

    # create dict of LO xsection values and errors
#set_trace()
LO_outdict = {}
for boson in bosons:
    for proc, proc_name in procs.items():
        for scale in scales:
            for xsec_type in types:
                for channel, chan_name in channels.items():
                    sig_dname = "_".join([boson, proc, mg5_LO_xsecs_name, scale, xsec_type, channel])
                    vals = dense_lookup(*fdict[(sig_dname, "dense_lookup")])
                    errs = dense_lookup(*fdict[(f"{sig_dname}_error", "dense_lookup")])
                    for mtt in masses:
                        for width in widths:
                            outname = f"{boson}toTTJets{chan_name}_M{mtt}_W{widthTOname(width)}_{proc_name}"
                            if scale != "nominal": outname = f"{outname}_{scale}"
                            LO_outdict[outname] = vals(mtt, width)

lo_xsec_fname = os.path.join(proj_dir, "inputs", "signal_xsecs.json")
with open(lo_xsec_fname, "w") as out:
    out.write(prettyjson.dumps(LO_outdict))
print(f"{lo_xsec_fname} written")
