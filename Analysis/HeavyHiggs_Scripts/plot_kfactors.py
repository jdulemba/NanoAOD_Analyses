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
rcParams["savefig.format"] = "pdf"
rcParams["savefig.bbox"] = "tight"

from pdb import set_trace
import numpy as np
import os
from coffea.lookup_tools.dense_lookup import dense_lookup
import Utilities.root_converters as root_conv

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--plot", choices=["LO", "KFactor", "All"], default="All", help="Choose which parity to compute limits for.")
parser.add_argument("--parity", choices=["A", "H", "All"], default="All", help="Choose which parity to compute limits for.")
args = parser.parse_args()


def make_plots(lookup, boson, channel, sig_comp, sys, isKFactor=False, clear=True):
    #set_trace()
    pltdir = os.path.join(plot_dir, base_jobid, "KFactors_Ratios", "KFactors" if isKFactor else "LOratios", f"{boson}toTT")
    if not os.path.isdir(pltdir): os.makedirs(pltdir)

    fig, ax = plt.subplots()
    fig.subplots_adjust(hspace=.07)

    zaxis_label = "k-factor" if isKFactor else "sushi/mg5"
    zaxis_name = "NNLO_kfactor" if isKFactor else "LO_ratios"

    xedges, yedges = lookup._axes
    vals = lookup._values
    pc = ax.pcolormesh(xedges, yedges, vals.T, **{})
    ax.add_collection(pc)
    if clear:
        if channel:
            zlabel = f"{zaxis_label}: {sys_opts_dict[sys]['label']}, {sig_components[sig_comp]}, {channel}"
        else:
            zlabel = f"{boson} {sig_components[sig_comp]} {zaxis_label}" if sys == "nominal" else f"{boson} {sig_components[sig_comp]} {zaxis_label}: {sys_opts_dict[sys]['label']}"
        fig.colorbar(pc, ax=ax, label=zlabel, pad=0.)
    ax.autoscale()

        ## set axes labels and titles
    ax.set_xlabel("$m_{%s}$ [GeV]" % boson)
    ax.set_ylabel(r"$\Gamma_{%s}$/$m_{%s}$ [%%]" % (boson, boson))

    figname = os.path.join(pltdir, f"{boson}toTTJets{channels_dict[channel]}_{sig_comp}_{sys_opts_dict[sys]['name']}_{zaxis_name}" if channel else f"{boson}toTT_{sig_comp}_{sys_opts_dict[sys]['name']}_{zaxis_name}")
    fig.savefig(figname)
    print(f"{figname} written")
    plt.close(fig)
    #set_trace()


def get_kfactor(parity):
    # modeled after this https://github.com/afiqaize/CombineHarvester/blob/ahtt_run2ul_dev/ahtt/scripts/make_datacard.py#L30-L40
    kfactors = {
        syst : {
            "Res" : dense_lookup(*nom_tgraph[(f"{parity}_res_sushi_nnlo_mg5_lo_kfactor_pdf_325500_{syst}", "dense_lookup")]),
            "Int" : dense_lookup(*nom_tgraph[(f"{parity}_int_sushi_nnlo_mg5_lo_kfactor_pdf_325500_{syst}", "dense_lookup")]),
        } for syst in scale_choices
    }
    kfactors["mt_up"] = {
            "Res" : dense_lookup(*mt_up_tgraph[(f"{parity}_res_sushi_nnlo_mg5_lo_kfactor_pdf_325500_nominal", "dense_lookup")]),
            "Int" : dense_lookup(*mt_up_tgraph[(f"{parity}_int_sushi_nnlo_mg5_lo_kfactor_pdf_325500_nominal", "dense_lookup")]),
    }
    kfactors["mt_dw"] = {
            "Res" : dense_lookup(*mt_dw_tgraph[(f"{parity}_res_sushi_nnlo_mg5_lo_kfactor_pdf_325500_nominal", "dense_lookup")]),
            "Int" : dense_lookup(*mt_dw_tgraph[(f"{parity}_int_sushi_nnlo_mg5_lo_kfactor_pdf_325500_nominal", "dense_lookup")]),
    }
    return kfactors


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
    base_jobid = os.environ["base_jobid"]
    plot_dir = os.environ["plots_dir"]

    nom_tgraph = root_conv.convert_TGraph_root_file("root://eosuser.cern.ch//eos/cms/store/user/afiqaize/ahtt_kfactor_sushi/ulkfactor_final_220129.root")
    mt_up_tgraph = root_conv.convert_TGraph_root_file("root://eosuser.cern.ch//eos/cms/store/user/afiqaize/ahtt_kfactor_sushi/ulkfactor_final_mt173p5_230317.root")
    mt_dw_tgraph = root_conv.convert_TGraph_root_file("root://eosuser.cern.ch//eos/cms/store/user/afiqaize/ahtt_kfactor_sushi/ulkfactor_final_mt171p5_230317.root")


    masses = [365, 380, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000] # available mass points (GeV)
    widths = [0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 8, 10, 13, 15, 18, 21, 25] # available widths (%)
    channels_dict = {
        "lj" : "SL",
        "ll" : "DiLep"
    }
    sig_components = {"Res" : "resonance", "Int_neg" : "negative int.", "Int_pos" : "positive int.", "Int" : "interference"}

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

    plots_to_make = ["LO", "KFactor"] if args.plot == "All" else [args.plot]
    parities_to_run = ["A", "H"] if args.parity == "All" else [args.parity]

    for parity in parities_to_run:
        if "KFactor" in plots_to_make:
            kfactors = get_kfactor(parity)
            for sys in kfactors.keys():
                for sig_comp in kfactors[sys].keys():
                    make_plots(lookup=kfactors[sys][sig_comp], boson=parity, channel=None, sig_comp=sig_comp, sys=sys, isKFactor=True)

        if "LO" in plots_to_make:
            for channel in channels_dict.keys():
                LO_ratios = get_lo_ratio(parity, channel)
                for sys in LO_ratios.keys():
                    if sys == "nominal": continue
                    for sig_comp in LO_ratios[sys].keys():
                        if np.any(LO_ratios[sys][sig_comp]._values > 2.0):
                            vals, edges = LO_ratios[sys][sig_comp]._values, LO_ratios[sys][sig_comp]._axes
                            failing_mass_inds, failing_width_inds = np.where(vals > 2.)
                            failing_points = [(masses[failing_mass_inds[idx]], widths[failing_width_inds[idx]]) for idx in range(len(failing_mass_inds))]
                            print(f"{parity} {sig_comp} {channel} {sys}: {failing_points}\n")
                            
                        make_plots(lookup=LO_ratios[sys][sig_comp], boson=parity, channel=channel, sig_comp=sig_comp, sys=sys, isKFactor=False)

    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
