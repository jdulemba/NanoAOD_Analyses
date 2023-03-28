#! /bin/env python

import time
tic = time.time()

# matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend("agg")
from matplotlib import rcParams
rcParams["font.size"] = 20
rcParams["savefig.format"] = "pdf"
rcParams["savefig.bbox"] = "tight"
#import matplotlib
#import matplotlib.patches as mpatches
#import matplotlib.lines as mlines
#import matplotlib.ticker as ticker

from pdb import set_trace
import numpy as np
import os
import uproot
from coffea.lookup_tools.dense_lookup import dense_lookup
import Utilities.root_converters as root_conv
from itertools import product

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--parity", choices=["A", "H", "All"], default="All", help="Choose which parity to compute limits for.")
args = parser.parse_args()


def plot_ratios(lookup, boson, channel, sig_comp, sys, clear=True):
    #set_trace()
    pltdir = os.path.join(plot_dir, base_jobid, "KFactors_Ratios", f"{boson}toTT")
    if not os.path.isdir(pltdir): os.makedirs(pltdir)

    fig, ax = plt.subplots()
    fig.subplots_adjust(hspace=.07)

    opts = {}
    #opts = {"cmap" : "OrRd"}
    xedges, yedges = lookup._axes
    vals = lookup._values
    pc = ax.pcolormesh(xedges, yedges, vals.T, **opts)
    ax.add_collection(pc)
    if clear:
        fig.colorbar(pc, ax=ax, label=f"sushi/mg5: {sys_opts_dict[sys]['label']}, {sig_components[sig_comp]}, {channel}", pad=0.)
    ax.autoscale()
    #ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    #ax.xaxis.set_minor_formatter(mpl.ticker.ScalarFormatter())
    #ax.ticklabel_format(style="plain", axis="x",useOffset=False)

        ## set axes labels and titles
    ax.set_xlabel("$m_{%s}$ [GeV]" % boson)
    ax.set_ylabel(r"$\Gamma_{%s}$/$m_{%s}$ [%%]" % (boson, boson))

    #set_trace()
    figname = os.path.join(pltdir, f"{boson}toTTJets{channels_dict[channel]}_{sig_comp}_{sys_opts_dict[sys]['name']}_LO_ratios")
    fig.savefig(figname)
    print(f"{figname} written")
    plt.close(fig)


def get_kfactor(sigpnt):
    set_trace()
    khist = {
        syst : [
            dense_lookup(*nom_tgraph[(f"{parity}_res_sushi_nnlo_mg5_lo_kfactor_pdf_325500_{syst}", "dense_lookup")),
            dense_lookup(*nom_tgraph[(f"{parity}_int_sushi_nnlo_mg5_lo_kfactor_pdf_325500_{syst}", "dense_lookup")),
        ] for syst in scale_choices
    }
        (kfile.Get(sigpnt[0] + "_res_sushi_nnlo_mg5_lo_kfactor_pdf_325500_" + syst),
         kfile.Get(sigpnt[0] + "_int_sushi_nnlo_mg5_lo_kfactor_pdf_325500_" + syst))
        for syst in scale_choices
    ]
    #kfile = uproot.open(kfactor_file_name)
    ##kfile = TFile.Open(kfactor_file_name, "read")
    #khist = [
    #    (kfile.Get(sigpnt[0] + "_res_sushi_nnlo_mg5_lo_kfactor_pdf_325500_" + syst),
    #     kfile.Get(sigpnt[0] + "_int_sushi_nnlo_mg5_lo_kfactor_pdf_325500_" + syst))
    #    for syst in scale_choices
    #]
    #kvals = tuple([(syst[0].Interpolate(sigpnt[1], sigpnt[2]), syst[1].Interpolate(sigpnt[1], sigpnt[2])) for syst in khist])
    #kfile.Close()

    return kvals

def get_lo_ratio(parity, channel):
    # modeled after this https://github.com/afiqaize/CombineHarvester/blob/ahtt_run2ul_dev/ahtt/scripts/make_datacard.py#L42-L59
    xhist = { syst :
        [
            dense_lookup(*nom_tgraph[(f"{parity}_res_mg5_pdf_325500_scale_dyn_0p5mtt_{syst}_xsec_{channel}", "dense_lookup")]),
            dense_lookup(*nom_tgraph[(f"{parity}_int_mg5_pdf_325500_scale_dyn_0p5mtt_{syst}_xabs_{channel}", "dense_lookup")]),
            dense_lookup(*nom_tgraph[(f"{parity}_int_mg5_pdf_325500_scale_dyn_0p5mtt_{syst}_positive_event_fraction_{channel}", "dense_lookup")])
         ] for syst in scale_choices
    }
        # add top mass variations
    xhist["mt_up"] = [
            dense_lookup(*mt_up_tgraph[(f"{parity}_res_mg5_pdf_325500_scale_dyn_0p5mtt_nominal_xsec_{channel}", "dense_lookup")]),
            dense_lookup(*mt_up_tgraph[(f"{parity}_int_mg5_pdf_325500_scale_dyn_0p5mtt_nominal_xabs_{channel}", "dense_lookup")]),
            dense_lookup(*mt_up_tgraph[(f"{parity}_int_mg5_pdf_325500_scale_dyn_0p5mtt_nominal_positive_event_fraction_{channel}", "dense_lookup")])
    ]
    xhist["mt_dw"] = [
            dense_lookup(*mt_dw_tgraph[(f"{parity}_res_mg5_pdf_325500_scale_dyn_0p5mtt_nominal_xsec_{channel}", "dense_lookup")]),
            dense_lookup(*mt_dw_tgraph[(f"{parity}_int_mg5_pdf_325500_scale_dyn_0p5mtt_nominal_xabs_{channel}", "dense_lookup")]),
            dense_lookup(*mt_dw_tgraph[(f"{parity}_int_mg5_pdf_325500_scale_dyn_0p5mtt_nominal_positive_event_fraction_{channel}", "dense_lookup")])
    ]
    lookups = {
        syst : [
            dist[0],
            dense_lookup(*(dist[1]._values * dist[2]._values, dist[0]._axes)),
            dense_lookup(*(dist[1]._values * (1. - dist[2]._values), dist[0]._axes))
        ] for syst, dist in xhist.items()
    }

        # syst/nominal values for [res, pos, neg] signal components
    lookups_ratios = {
        syst : {
            "Res" : dense_lookup(*(lookups[syst][0]._values / lookups["nominal"][0]._values, lookups[syst][0]._axes)),
            "Int_pos" : dense_lookup(*(lookups[syst][1]._values / lookups["nominal"][1]._values, lookups[syst][0]._axes)),
            "Int_neg" : dense_lookup(*(lookups[syst][2]._values / lookups["nominal"][2]._values, lookups[syst][0]._axes))
        } for syst in lookups.keys()
    }
    return lookups_ratios



if __name__ == "__main__":	
    #set_trace()
    base_jobid = os.environ["base_jobid"]
    plot_dir = os.environ["plots_dir"]

    rfiles = {
        "Nominal" : "root://eosuser.cern.ch//eos/cms/store/user/afiqaize/ahtt_kfactor_sushi/ulkfactor_final_220129.root",
        "Mt_1715" : "root://eosuser.cern.ch//eos/cms/store/user/afiqaize/ahtt_kfactor_sushi/ulkfactor_final_mt171p5_230317.root",
        "Mt_1735" : "root://eosuser.cern.ch//eos/cms/store/user/afiqaize/ahtt_kfactor_sushi/ulkfactor_final_mt173p5_230317.root",
    }
    nom_tgraph = root_conv.convert_TGraph_root_file(rfiles["Nominal"])
    mt_up_tgraph = root_conv.convert_TGraph_root_file(rfiles["Mt_1735"])
    mt_dw_tgraph = root_conv.convert_TGraph_root_file(rfiles["Mt_1715"])


    masses = [365, 380, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000] # available mass points (GeV)
    widths = [0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 8, 10, 13, 15, 18, 21, 25] # available widths (%)
    channels_dict = {
        "lj" : "SL",
        "ll" : "DiLep"
    }
    sig_components = {"Res" : "resonant", "Int_neg" : "negative int.", "Int_pos" : "positive int."}

    scale_choices = ["nominal", "uF_up", "uF_down", "uR_up", "uR_down"]
    sys_opts_dict = {
        "nominal" : {"label" : "", "name" : "Nominal"},
        "uF_up" : {"label" : "$\mu_{F}$ up", "name" : "uF_up"},
        "uF_down" : {"label" : "$\mu_{F}$ down", "name" : "uF_down"},
        "uR_up" : {"label" : "$\mu_{R}$ up", "name" : "uR_up"},
        "uR_down" : {"label" : "$\mu_{R}$ down", "name" : "uR_down"},
        "mt_up" : {"label" : "$m_{t}$ up", "name" : "mt_up"},
        "mt_dw" : {"label" : "$m_{t}$ down", "name" : "mt_down"},
    }

    parities_to_run = ["A", "H"] if args.parity == "All" else [args.parity]
    #set_trace()
    for parity in parities_to_run:
        for channel in channels_dict.keys():
            LO_ratios = get_lo_ratio(parity, channel)
            for sys in LO_ratios.keys():
                if sys == "nominal": continue
                for sig_comp in LO_ratios[sys].keys():
                    plot_ratios(lookup=LO_ratios[sys][sig_comp], boson=parity, channel=channel, sig_comp=sig_comp, sys=sys)

    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
