import time
tic = time.time()

import uproot
import numpy as np
from collections import defaultdict
import hist
from pdb import set_trace

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("outputfile")
args = parser.parse_args()

#proc = "TT"
nnlo_name = "NNLOqcd"
lo_name = "LOqcd"

channels = ["em", "ee", "mm"]
years = ["2016pre", "2016post", "2017", "2018"]

templates_to_copy = defaultdict(dict) ## make dict to store variations which are just copied

#set_trace()
#orig_rfilename = "/eos/user/j/jdulemba/NanoAOD_Analyses/results/Summer20UL_DeepJet/LL_Templates/AH_templates_021023/before_yt_ll.root"
ewk_yt_rfilename = "/eos/user/j/jdulemba/NanoAOD_Analyses/results/Summer20UL_DeepJet/LL_Templates/AH_templates_021023/after_yt_ll.root"
output_rfilename = "/eos/user/j/jdulemba/NanoAOD_Analyses/results/Summer20UL_DeepJet/LL_Templates/templates_ll_bkg.root"

with uproot.open(ewk_yt_rfilename) as rfile:
    for channel in channels:
        for year in years:
            cat = channel + "_" + year
            print("Loading templates for " + cat)
            if cat in rfile:
                rcat = rfile[cat]
                    ## save 'TT' from post yt file as TT_NNLOqcd
                templates_to_copy[cat][nnlo_name] = np.stack([np.copy(rcat["TT"].values()), np.copy(rcat["TT"].variances())], axis=-1)
                ewk_tt_templates = sorted(set([key.split(";")[0] for key in rcat.keys() if key.startswith("EWK_TT")]))
                for ewk_tt_temp in ewk_tt_templates:
                    templates_to_copy[cat][ewk_tt_temp] = np.copy(rcat[ewk_tt_temp].values())
                    # save toponium
                templates_to_copy[cat]["EtaT"] = np.stack([np.copy(rcat["EtaT"].values()), np.copy(rcat["EtaT"].variances())], axis=-1)

#set_trace()

    # add TT_NNLO dist to output rfile
with uproot.update(output_rfilename) as rfile:
    for cat in templates_to_copy.keys():
        print(f"Adding templates for {cat}")
        rfile[cat + "/TT_NNLOqcd"] = hist.Hist(hist.axis.Variable(np.arange(len(templates_to_copy[cat][nnlo_name])+1)), "Weight", data=templates_to_copy[cat][nnlo_name])
        rfile[cat + "/EtaT"] = hist.Hist(hist.axis.Variable(np.arange(len(templates_to_copy[cat]["EtaT"])+1)), "Weight", data=templates_to_copy[cat]["EtaT"])


output_templates = {}
#set_trace()
for cat in templates_to_copy.keys():
    nnloqcd_temp = templates_to_copy[cat][nnlo_name]
    bin_edges = np.arange(len(nnloqcd_temp)+1)

        # find individual coefficients
    b0_pos, b0_neg = templates_to_copy[cat]["EWK_TT_const_pos"], templates_to_copy[cat]["EWK_TT_const_neg"]
    b1_pos, b1_neg = templates_to_copy[cat]["EWK_TT_lin_pos"], templates_to_copy[cat]["EWK_TT_lin_neg"]
    b2_pos, b2_neg = templates_to_copy[cat]["EWK_TT_quad_pos"], templates_to_copy[cat]["EWK_TT_quad_neg"]

        # combine coefficients
    b0, b1, b2 = b0_pos - b0_neg, b1_pos - b1_neg, b2_pos - b2_neg 
    variances = np.zeros_like(b0)

        ## equation for TT_LOqcd dist when yt = 0
    tt_loqcd = np.copy((templates_to_copy[cat]["EWK_TT_const_pos_EWK_schemeUp"] - templates_to_copy[cat]["EWK_TT_const_neg_EWK_schemeUp"])/b0 * nnloqcd_temp[:,0])


    #output_templates[cat + "/TT_NNLOqcd"] = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=nnloqcd_temp)
    output_templates[cat + "/TT_LOqcd"] = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([tt_loqcd, variances], axis=-1))
    output_templates[cat + "/EWK_TT_const"] = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([b0, variances], axis=-1))
    output_templates[cat + "/EWK_TT_lin"]   = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([b1, variances], axis=-1))
    output_templates[cat + "/EWK_TT_quad"]  = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([b2, variances], axis=-1))

import time
named_tuple = time.localtime() # get struct_time
time_str = time.strftime("%d%m%Y", named_tuple)
output_rname = "_".join([(args.outputfile).split(".root")[0], time_str])+".root" ## append date onto root filename
with uproot.recreate(output_rname) as rfile:
    for key, hist in output_templates.items():
        channel, sys = key.split("/")
        print(f"Calculating templates for {channel} {sys}")
        rfile[key] = hist.copy()

print(f"\n\nFinished writing {output_rname} to disk")

toc = time.time()
print("Total time: %.1f" % (toc - tic))
