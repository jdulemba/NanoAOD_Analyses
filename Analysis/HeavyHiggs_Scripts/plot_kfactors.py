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

__allowed_types__ = ["LOvals", "LOratios", "KFactorVals", "KFactorRatios", "NNLOvals", "NNLOratios"]

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("plot", choices=__allowed_types__, help="Choose which parity to compute limits for.")
parser.add_argument("--parity", choices=["A", "H", "All"], default="All", help="Choose which parity to compute limits for.")
args = parser.parse_args()


def make_plots(lookup, boson, channel, sig_comp, sys, clear=True):
    #set_trace()
    pltdir = os.path.join(plot_dir, base_jobid, "KFactors_LOratios", file_version, args.plot, f"{boson}toTT")
    if not os.path.isdir(pltdir): os.makedirs(pltdir)

    fig, ax = plt.subplots()
    fig.subplots_adjust(hspace=.07)

    if args.plot == "LOvals":
        zlabel = f"{boson} {channel} {sig_components[sig_comp]}" + " $\\sigma$ [pb]"
        if sys != "nominal": zlabel += f": {sys_opts_dict[sys]['label']}"
        hep.label.lumitext(ax=ax, text="Madgraph LO QCD")
        figname = os.path.join(pltdir, f"{boson}toTTJets{channels_dict[channel]}_{sig_comp}_{sys_opts_dict[sys]['name']}_{args.plot}")

    if args.plot == "LOratios":
        zlabel = f"{boson} {channel} {sig_components[sig_comp]}" + " $\\sigma_{%s}$/$\\sigma_{nominal}$" % sys_opts_dict[sys]["label"].replace("$", "").replace(" ", "\ ")
        hep.label.lumitext(ax=ax, text="Madgraph LO QCD")
        figname = os.path.join(pltdir, f"{boson}toTTJets{channels_dict[channel]}_{sig_comp}_{sys_opts_dict[sys]['name']}_{args.plot}")

    if args.plot == "KFactorVals":
        zlabel = f"{boson} {sig_components[sig_comp]} k-factor" if sys == "nominal" else f"{boson} {sig_components[sig_comp]} {sys_opts_dict[sys]['label']} k-factor"
        hep.label.lumitext(ax=ax, text="Sushi NNLO / Madgraph LO")
        figname = os.path.join(pltdir, f"{boson}toTT_{sig_comp}_{sys_opts_dict[sys]['name']}_{args.plot}")

    if args.plot == "KFactorRatios":
        zlabel = f"{boson} {sig_components[sig_comp]} {sys_opts_dict[sys]['label']}/nominal k-factor ratios"
        hep.label.lumitext(ax=ax, text="Sushi NNLO / Madgraph LO")
        figname = os.path.join(pltdir, f"{boson}toTT_{sig_comp}_{sys_opts_dict[sys]['name']}_{args.plot}")

    if args.plot == "NNLOvals":
        zlabel = f"{boson} {channel} {sig_components[sig_comp]}" + " $\\sigma$ [pb]"
        if sys != "nominal": zlabel += f": {sys_opts_dict[sys]['label']}"
        hep.label.lumitext(ax=ax, text="Sushi NNLO")
        figname = os.path.join(pltdir, f"{boson}toTTJets{channels_dict[channel]}_{sig_comp}_{sys_opts_dict[sys]['name']}_{args.plot}")

    if args.plot == "NNLOratios":
        zlabel = f"{boson} {channel} {sig_components[sig_comp]}" + " $\\sigma_{%s}$/$\\sigma_{nominal}$" % sys_opts_dict[sys]["label"].replace("$", "").replace(" ", "\ ")
        hep.label.lumitext(ax=ax, text="Sushi NNLO")
        figname = os.path.join(pltdir, f"{boson}toTTJets{channels_dict[channel]}_{sig_comp}_{sys_opts_dict[sys]['name']}_{args.plot}")


    xedges, yedges = lookup._axes
    vals = lookup._values
    pc = ax.pcolormesh(xedges, yedges, vals.T, **{})
    ax.add_collection(pc)
    if clear:
        fig.colorbar(pc, ax=ax, label=zlabel, pad=0.)
    ax.autoscale()

        ## set axes labels and titles
    ax.set_xlabel("$m_{%s}$ [GeV]" % boson)
    ax.set_ylabel(r"$\Gamma_{%s}$/$m_{%s}$ [%%]" % (boson, boson))

    fig.savefig(figname)
    print(f"{figname} written")
    plt.close(fig)
    #set_trace()


def get_kfactor(parity):
    # modeled after this https://github.com/afiqaize/CombineHarvester/blob/ahtt_run2ul_dev/ahtt/scripts/make_datacard.py#L30-L40
    kfactor_vals = {
        syst : {
            "Res" : dense_lookup(*nom_tgraph[(f"{parity}_res_sushi_nnlo_mg5_lo_kfactor_pdf_325500_{syst}", "dense_lookup")]),
            "Int" : dense_lookup(*nom_tgraph[(f"{parity}_int_sushi_nnlo_mg5_lo_kfactor_pdf_325500_{syst}", "dense_lookup")]),
        } for syst in scale_choices
    }
    kfactor_vals["mt_up"] = {
            "Res" : dense_lookup(*mt_up_tgraph[(f"{parity}_res_sushi_nnlo_mg5_lo_kfactor_pdf_325500_nominal", "dense_lookup")]),
            "Int" : dense_lookup(*mt_up_tgraph[(f"{parity}_int_sushi_nnlo_mg5_lo_kfactor_pdf_325500_nominal", "dense_lookup")]),
    }
    kfactor_vals["mt_dw"] = {
            "Res" : dense_lookup(*mt_dw_tgraph[(f"{parity}_res_sushi_nnlo_mg5_lo_kfactor_pdf_325500_nominal", "dense_lookup")]),
            "Int" : dense_lookup(*mt_dw_tgraph[(f"{parity}_int_sushi_nnlo_mg5_lo_kfactor_pdf_325500_nominal", "dense_lookup")]),
    }
        # syst/nominal values for signal components
    kfactor_ratios = {
        syst : {
            sig_comp : dense_lookup(*(kfactor_vals[syst][sig_comp]._values / kfactor_vals["nominal"][sig_comp]._values, kfactor_vals[syst][sig_comp]._axes)) for sig_comp in kfactor_vals[syst].keys()
        } for syst in kfactor_vals.keys()
    }
    return kfactor_vals, kfactor_ratios


def get_lo(parity, channel):
    # modeled after this https://github.com/afiqaize/CombineHarvester/blob/ahtt_run2ul_dev/ahtt/scripts/make_datacard.py#L42-L59
    lo_dict = { syst :
        [
            dense_lookup(*nom_tgraph[(f"{parity}_res_mg5_pdf_325500_scale_dyn_0p5mtt_{syst}_xsec_{channel}", "dense_lookup")]),
            dense_lookup(*nom_tgraph[(f"{parity}_int_mg5_pdf_325500_scale_dyn_0p5mtt_{syst}_xabs_{channel}", "dense_lookup")]),
            dense_lookup(*nom_tgraph[(f"{parity}_int_mg5_pdf_325500_scale_dyn_0p5mtt_{syst}_positive_event_fraction_{channel}", "dense_lookup")])
         ] for syst in scale_choices
    }
        # add top mass variations
    lo_dict["mt_up"] = [
            dense_lookup(*mt_up_tgraph[(f"{parity}_res_mg5_pdf_325500_scale_dyn_0p5mtt_nominal_xsec_{channel}", "dense_lookup")]),
            dense_lookup(*mt_up_tgraph[(f"{parity}_int_mg5_pdf_325500_scale_dyn_0p5mtt_nominal_xabs_{channel}", "dense_lookup")]),
            dense_lookup(*mt_up_tgraph[(f"{parity}_int_mg5_pdf_325500_scale_dyn_0p5mtt_nominal_positive_event_fraction_{channel}", "dense_lookup")])
    ]
    lo_dict["mt_dw"] = [
            dense_lookup(*mt_dw_tgraph[(f"{parity}_res_mg5_pdf_325500_scale_dyn_0p5mtt_nominal_xsec_{channel}", "dense_lookup")]),
            dense_lookup(*mt_dw_tgraph[(f"{parity}_int_mg5_pdf_325500_scale_dyn_0p5mtt_nominal_xabs_{channel}", "dense_lookup")]),
            dense_lookup(*mt_dw_tgraph[(f"{parity}_int_mg5_pdf_325500_scale_dyn_0p5mtt_nominal_positive_event_fraction_{channel}", "dense_lookup")])
    ]

        # sushi xsec values for [res, pos, neg] signal components
    lo_vals = {
        syst : {
            "Res" : dist[0],
            "Int_pos" : dense_lookup(*(dist[1]._values * dist[2]._values, dist[0]._axes)),
            "Int_neg" : dense_lookup(*(dist[1]._values * (1. - dist[2]._values), dist[0]._axes))
        } for syst, dist in lo_dict.items()
    }
        # syst/nominal values for [res, pos, neg] signal components
    lo_ratios = {
        syst : {
            sig_comp : dense_lookup(*(lo_vals[syst][sig_comp]._values / lo_vals["nominal"][sig_comp]._values, lo_vals[syst][sig_comp]._axes)) for sig_comp in lo_vals[syst].keys()
        } for syst in lo_vals.keys()
    }
    return lo_vals, lo_ratios


def get_nnlo(kfactor_vals, lo_vals):
    #set_trace()
        # (NNLO kfactor) * (LO mg5 xsec) for each systematic variation
    final_norms = {
        syst : {
            "Res" : dense_lookup(*(kfactor_vals[syst]["Res"]._values * lo_vals[syst]["Res"]._values, kfactor_vals[syst]["Res"]._axes)),
            "Int_pos" : dense_lookup(*(kfactor_vals[syst]["Int"]._values * lo_vals[syst]["Int_pos"]._values, kfactor_vals[syst]["Int"]._axes)),
            "Int_neg" : dense_lookup(*(kfactor_vals[syst]["Int"]._values * lo_vals[syst]["Int_neg"]._values, kfactor_vals[syst]["Int"]._axes)),
        } for syst in lo_vals.keys()
    }

        # (kfactor * (LO ratio))_sys / (kfactor * (LO ratio))_nominal for each systematic variation
    final_norm_ratios = {
        syst : {
            sig_comp : dense_lookup(*(final_norms[syst][sig_comp]._values / final_norms["nominal"][sig_comp]._values, final_norms[syst][sig_comp]._axes)) for sig_comp in final_norms[syst].keys()
        } for syst in final_norms.keys()
    }
    return final_norms, final_norm_ratios

if __name__ == "__main__":	
    base_jobid = os.environ["base_jobid"]
    plot_dir = os.environ["plots_dir"]

    #file_version = "230317"
    file_version = "230329"
    nom_tgraph = root_conv.convert_TGraph_root_file("root://eosuser.cern.ch//eos/cms/store/user/afiqaize/ahtt_kfactor_sushi/ulkfactor_final_220129.root")
    mt_up_tgraph = root_conv.convert_TGraph_root_file(f"root://eosuser.cern.ch//eos/cms/store/user/afiqaize/ahtt_kfactor_sushi/ulkfactor_final_mt173p5_{file_version}.root")
    mt_dw_tgraph = root_conv.convert_TGraph_root_file(f"root://eosuser.cern.ch//eos/cms/store/user/afiqaize/ahtt_kfactor_sushi/ulkfactor_final_mt171p5_{file_version}.root")


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

    parities_to_run = ["A", "H"] if args.parity == "All" else [args.parity]

    for parity in parities_to_run:
            # (kfactor vals) * (LO vals) = (NNLO sushi vals)
        if (args.plot == "NNLOvals") or (args.plot == "NNLOratios"):
            kfactor_vals, kfactor_ratios = get_kfactor(parity)
            for channel in channels_dict.keys():
                lo_vals, lo_ratios = get_lo(parity, channel)
                nnlo_vals, nnlo_ratios = get_nnlo(kfactor_vals, lo_vals)
                for sys in nnlo_vals.keys():
                    for sig_comp in nnlo_vals[sys].keys():
                        if args.plot == "NNLOvals":
                            make_plots(lookup=nnlo_vals[sys][sig_comp], boson=parity, channel=channel, sig_comp=sig_comp, sys=sys)

                        if (sys != "nominal") and (args.plot == "NNLOratios"):
                            make_plots(lookup=nnlo_ratios[sys][sig_comp], boson=parity, channel=channel, sig_comp=sig_comp, sys=sys)


        if (args.plot == "KFactorVals") or (args.plot == "KFactorRatios"):
            kfactor_vals, kfactor_ratios = get_kfactor(parity)
            for sys in kfactor_vals.keys():
                for sig_comp in kfactor_vals[sys].keys():
                    if args.plot == "KFactorVals":
                        make_plots(lookup=kfactor_vals[sys][sig_comp], boson=parity, channel=None, sig_comp=sig_comp, sys=sys)

                    if (sys != "nominal") and (args.plot == "KFactorRatios"):
                        make_plots(lookup=kfactor_ratios[sys][sig_comp], boson=parity, channel=None, sig_comp=sig_comp, sys=sys)


        if (args.plot == "LOvals") or (args.plot == "LOratios"):
            for channel in channels_dict.keys():
                LO_vals, LO_ratios = get_lo(parity, channel)
                for sys in LO_vals.keys():
                    for sig_comp in LO_vals[sys].keys():
                        if args.plot == "LOvals":
                            make_plots(lookup=LO_vals[sys][sig_comp], boson=parity, channel=channel, sig_comp=sig_comp, sys=sys)

                        if (sys != "nominal") and (args.plot == "LOratios"):
                            if np.any(LO_ratios[sys][sig_comp]._values > 2.0):
                                vals, edges = LO_ratios[sys][sig_comp]._values, LO_ratios[sys][sig_comp]._axes
                                failing_mass_inds, failing_width_inds = np.where(vals > 2.)
                                failing_points = [(masses[failing_mass_inds[idx]], widths[failing_width_inds[idx]]) for idx in range(len(failing_mass_inds))]
                                print(f"{parity} {sig_comp} {channel} {sys}: {failing_points}\n")
                                
                            make_plots(lookup=LO_ratios[sys][sig_comp], boson=parity, channel=channel, sig_comp=sig_comp, sys=sys)

    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
