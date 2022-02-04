#!/usr/bin/env python

from pdb import set_trace
import os
import numpy as np
import Utilities.prettyjson as prettyjson
import Utilities.root_converters as root_conv

proj_dir = os.environ["PROJECT_DIR"]

kfactor_fname = "ulkfactor_final_220129"
rname = os.path.join(proj_dir, "ahtt_kfactor_sushi", f"{kfactor_fname}.root")
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

fdict = root_conv.convert_TGraph_root_file(rname)

widthTOname = lambda width : str(float(width)).replace('.', 'p')

# create dict of the fraction of positive events for each signal points
#outdict = {}
#for boson in bosons:
#    sig_dname = "_".join([boson, pos_evt_fraction_name]) # name of signal dist
#    vals = dense_lookup(*fdict[(sig_dname, "dense_lookup")])
#    errs = dense_lookup(*fdict[(f"{sig_dname}_error", "dense_lookup")])

#set_trace()

    # create json for NNLO/LO cross section k-factors
kfactors_outdict = {}
for boson in bosons:
    for proc, proc_name in procs.items():
        for scale in scales:
            sig_dname = "_".join([boson, proc, kfactors_name, scale])
            vals_dict = fdict[sig_dname]
            for mtt in masses:
                for width in widths:
                    for channel, chan_name in channels.items():
                        outname = f"{boson}toTTJets{chan_name}_M{mtt}_W{widthTOname(width)}_{proc_name}"
                        if scale != "nominal": outname = f"{outname}_{scale}"
                        kfactors_outdict[outname] = vals_dict[(mtt, width)]

#set_trace()
kfactors_fname = os.path.join(proj_dir, "inputs", f"signal_kfactors_{kfactor_fname}.json")
#kfactors_fname = os.path.join(proj_dir, "inputs", "signal_kfactors.json")
with open(kfactors_fname, "w") as out:
    out.write(prettyjson.dumps(kfactors_outdict))
print(f"{kfactors_fname} written")


    # create dict of LO xsection values and errors
LO_outdict_xsec = {}
for boson in bosons:
    for proc, proc_name in procs.items():
        for scale in scales:
            for xsec_type in ["xsec"]:
                for channel, chan_name in channels.items():
                    sig_dname = "_".join([boson, proc, mg5_LO_xsecs_name, scale, xsec_type, channel])
                    vals_dict = fdict[sig_dname]
                    for mtt in masses:
                        for width in widths:
                            outname = f"{boson}toTTJets{chan_name}_M{mtt}_W{widthTOname(width)}_{proc_name}"
                            if scale != "nominal": outname = f"{outname}_{scale}"
                            LO_outdict_xsec[outname] = vals_dict[(mtt, width)]

lo_xsec_fname = os.path.join(proj_dir, "inputs", f"signal_xsecs_{kfactor_fname}.json")
#lo_xsec_fname = os.path.join(proj_dir, "inputs", "signal_xsecs.json")
with open(lo_xsec_fname, "w") as out:
    out.write(prettyjson.dumps(LO_outdict_xsec))
print(f"{lo_xsec_fname} written")


    # create json for absolute cross section
LO_outdict_xabs = {}
for boson in bosons:
    for proc, proc_name in procs.items():
        for scale in scales:
            for xsec_type in ["xabs"]:
                for channel, chan_name in channels.items():
                    sig_dname = "_".join([boson, proc, mg5_LO_xsecs_name, scale, xsec_type, channel])
                    vals_dict = fdict[sig_dname]
                    for mtt in masses:
                        for width in widths:
                            outname = f"{boson}toTTJets{chan_name}_M{mtt}_W{widthTOname(width)}_{proc_name}"
                            if scale != "nominal": outname = f"{outname}_{scale}"
                            LO_outdict_xabs[outname] = vals_dict[(mtt, width)]

lo_xabs_fname = os.path.join(proj_dir, "inputs", f"signal_xabs_{kfactor_fname}.json")
#lo_xabs_fname = os.path.join(proj_dir, "inputs", "signal_xabs.json")
with open(lo_xabs_fname, "w") as out:
    out.write(prettyjson.dumps(LO_outdict_xabs))
print(f"{lo_xabs_fname} written")
