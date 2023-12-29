import time
tic = time.time()

import uproot
import argparse
import numpy as np
from collections import defaultdict
import hist
from pdb import set_trace
import Utilities.common_features as cfeatures

parser = argparse.ArgumentParser()
parser.add_argument("inputtemplates")
parser.add_argument("outputfile")
parser.add_argument("format", choices=["LJ", "LL"], help="Choose which formatting of output root file to make. One used to give to dilepton or one for Otto.")
parser.add_argument("--plot", action="store_true", help="Make plots of sys uncs")
args = parser.parse_args()

proc = "TT"
ewk_name = "EWK_yukawa"
nnlo_name = "NNLOqcd"
lo_name = "LOqcd"

#channels = ["em", "ee", "mm", "e3jets", "mu3jets", "e4pjets", "mu4pjets"]
channels = ["e3jets", "mu3jets", "e4pjets", "mu4pjets"]
years = ["2016pre", "2016post", "2017", "2018"]

# Weight-based systematics to copy to the EWK templates
# It is assumed that EWK_yukawa and the other nuisance are uncorellated. Since they are both
# weight-based, the template for both systematics varied by 1sigma is then just
# (EW_yukawa up/down) * (other up/down) / (nominal)
# Then, the same quadratic equation is solved for the systematic as for the nominal to get
# the variation templates for the EWK signal.

if args.format == "LJ": systs_to_copy = []
else:
    systs_to_copy = [
        "EWK_scheme",
        "CMS_eff_e_reco",
        "CMS_eff_e_id",
        "CMS_eff_m_id_stat",
        "CMS_eff_m_id_syst",
        "CMS_eff_m_iso_stat",
        "CMS_eff_m_iso_syst",
        "CMS_eff_trigger_ee",
        "CMS_eff_trigger_em",
        "CMS_eff_trigger_mm",
        "CMS_L1_prefire",
        "CMS_eff_trigger_m_syst",
        "CMS_eff_trigger_m_stat",
        "CMS_eff_trigger_e",
    ]

dyt_up = 0.11
dyt_dw = -0.12

yt_nom = 1.0
yt_up = yt_nom + dyt_up
yt_dw = yt_nom + dyt_dw

TT_templates = defaultdict(dict)

templates_to_copy = defaultdict(dict) ## make dict to store variations which are just copied
TT_templates_to_ignore = ["TT_LOqcd", "TT_NNLO_dQCDUp", "TT_NNLO_dQCDDown"] if args.format == "LL" else ["TT_NNLO_dQCDUp", "TT_NNLO_dQCDDown"]

with uproot.open(args.inputtemplates) as rfile:
    for channel in channels:
        for year in years:
            cat = channel + "_" + year
            print("Loading templates for " + cat)
            if cat in rfile:
                rcat = rfile[cat]
                TT_templates[cat]["nominal"] = np.copy(rcat[proc].values())

                if args.format == "LJ":
                    TT_templates[cat][ewk_name + "Up"]   = np.copy(rcat[proc + "_" + ewk_name + "Up"].values())
                    TT_templates[cat][ewk_name + "Down"] = np.copy(rcat[proc + "_" + ewk_name + "Down"].values())
                    TT_templates[cat][nnlo_name] = np.stack([np.copy(rcat[proc + "_" + nnlo_name].values()), np.copy(rcat[proc + "_" + nnlo_name].variances())], axis=-1)
                    TT_templates[cat][lo_name] = np.stack([np.copy(rcat[proc + "_" + lo_name].values()), np.copy(rcat[proc + "_" + lo_name].variances())], axis=-1)
                    for sys in systs_to_copy:
                        for sysdir in ["Up", "Down"]:
                            syskey = proc + "_" + sys + sysdir
                            if syskey in rcat:
                                TT_templates[cat][sys + sysdir] = rcat[syskey].values()

                if args.format == "LL":
                    processes = sorted(set([key.split(";")[0].split("_")[0] for key in rcat.keys()]))
                    processes = [process.replace("data", "data_obs") for process in processes]
                        # loop over all processes
                    for process in processes:
                            # find all available distributions for each process
                        proc_dists = sorted(set([key.split(";")[0] for key in rcat.keys() if key.startswith(process)]))
                        for distname in proc_dists:
                            if process == "TT":
                                if distname == process: continue
                                if ("Up" in distname) or ("Down" in distname):
                                    TT_templates[cat][distname.replace(f"{proc}_", "")] = np.copy(rcat[distname].values())
                                else:
                                    TT_templates[cat][distname.replace(f"{proc}_", "")] = np.stack([np.copy(rcat[distname].values()), np.copy(rcat[distname].variances())], axis=-1)
                            else:
                                templates_to_copy[cat][distname.replace("Eb", "bindingEnergy_EtaT") if ((process == "EtaT") and ("Eb" in distname)) else distname] = np.stack([np.copy(rcat[distname].values()), np.copy(rcat[distname].variances())], axis=-1)

"""
The yields for TTbar distributions at a given yt value are given by 
 TT(yt) = TT_NNLOqcd + TT_ewk(yt) == TT_NNLOqcd + b_0 + b_1 * (yt) + b_2 * (yt)^2
where TT_NNLOqcd is the contribution from including NNLOqcd terms via the NNLO QCD correction and 
 TT_ewk(yt) is the contribution from the NLO EWK terms via the EWK corrections

 TT = TT_NNLOqcd + b_0 + b_1 * (yt_nom) + b_2 * (yt_nom)^2
 TT_YukawaUp = TT_NNLOqcd + b_0 + b_1 * (yt_up) + b_2 * (yt_up)^2
 TT_YukawaDown = TT_NNLOqcd + b_0 + b_1 * (yt_dw) + b_2 * (yt_dw)^2


We can solve for b_0, b_1, and b_2 through the matrices Ax = B

   | 1.00   yt_nom   yt_nom^2 |   | b_0 | = | TT            - TT_NNLOqcd |
   | 1.00   yt_up    yt_up^2  |   | b_1 | = | TT_YukawaUp   - TT_NNLOqcd |
   | 1.00   yt_dw    yt_dw^2  |   | b_2 | = | TT_YukawaDown - TT_NNLOqcd |

which is the following with numpy (per bin of each distribution)

a = np.array([[1., yt_nom, yt_nom^2], [1., yt_up, yt_up^2], [1., yt_dw, yt_dw^2]])
b = np.array([TT - TT_NNLOqcd, TT_YukawaUp - TT_NNLOqcd, TT_YukawaDown - TT_NNLOqcd])
try:
    x = np.linalg.solve(a, b)
except LinAlgError:
    x = np.linalg.lstsq(a, b)[0]
"""

    ## matrix of coefficients from systems of equations, same for all channels
coeffs_matrix = np.array([[1., yt_nom, yt_nom**2], [1., yt_up, yt_up**2], [1., yt_dw, yt_dw**2]])

output_templates = {}

for cat in TT_templates.keys():
    template_nnloqcd = np.copy(TT_templates[cat][nnlo_name])
    template_loqcd   = np.copy(TT_templates[cat][lo_name])
    deltaQCD = np.copy((template_nnloqcd[:, 0] - template_loqcd[:, 0])/template_nnloqcd[:, 0])

    # Add bin edges so that uproot saves it as a TH1D
    bin_edges = np.arange(len(template_loqcd[:, 0])+1, dtype=int)
    variances = np.zeros_like(template_loqcd[:, 0])

    template_yukawa_nominal, template_yukawa_up, template_yukawa_down = np.copy(TT_templates[cat]["nominal"]), np.copy(TT_templates[cat]["EWK_yukawaUp"]), np.copy(TT_templates[cat]["EWK_yukawaDown"])

        ## find the terms that the ewk corrections are responsible for, dependent on yt (RHS of matrix equation)
    ewk_dists = [[template_yukawa_nominal[bin] - template_nnloqcd[bin, 0], template_yukawa_up[bin] - template_nnloqcd[bin, 0], template_yukawa_down[bin] - template_nnloqcd[bin, 0]] 
        for bin in range(template_yukawa_nominal.size)]


        # solve for EWK term coefficients
    try:
        b0, b1, b2 = np.array([np.linalg.solve(coeffs_matrix, np.array(bin)) for bin in ewk_dists]).T # transpose to make resulting vector size (nbins x 3)
    except LinAlgError:
        raise ValueError("Unable to solve for coefficients!")

    # Check that this reproduces the input templates
    def predict_yukawa(yt):
        return template_nnloqcd[:, 0] + b0 + b1 * yt + b2 * yt**2

    assert np.all(np.isclose(predict_yukawa(yt_nom), template_yukawa_nominal)) # nominal
    assert np.all(np.isclose(predict_yukawa(yt_up), template_yukawa_up)) # up
    assert np.all(np.isclose(predict_yukawa(yt_dw), template_yukawa_down)) # down


    if args.format == "LJ":
        print(f"Calculating templates for {cat}")
        hist_loqcd   = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=TT_templates[cat][lo_name])
        output_templates[cat + "/TT_" + lo_name] = hist_loqcd
        output_templates[cat + "/EWK_TT_const"] = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([b0, variances], axis=-1))
        output_templates[cat + "/EWK_TT_lin"]   = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([b1, variances], axis=-1))
        output_templates[cat + "/EWK_TT_quad"]  = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([b2, variances], axis=-1))

    if args.format == "LL":
        for sys in TT_templates[cat].keys():
            print("Calculating templates for " + cat + " " + sys)

            sys_template = np.copy(TT_templates[cat][sys])
            if sys == "nominal":
                syskey = ""
                b0_pos = np.where( b0 >= 0.,  b0, 0. )
                b0_neg = np.where( b0 <= 0., -b0, 0. )
                b1_pos = np.where( b1 >= 0.,  b1, 0. )
                b1_neg = np.where( b1 <= 0., -b1, 0. )
                b2_pos = np.where( b2 >= 0.,  b2, 0. )
                b2_neg = np.where( b2 <= 0., -b2, 0. )
                output_templates[cat + "/EWK_TT_const_pos" + syskey] = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([np.copy(b0_pos), np.copy(variances)], axis=-1))
                output_templates[cat + "/EWK_TT_const_neg" + syskey] = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([np.copy(b0_neg), np.copy(variances)], axis=-1))
                output_templates[cat + "/EWK_TT_lin_pos" + syskey]   = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([np.copy(b1_pos), np.copy(variances)], axis=-1))
                output_templates[cat + "/EWK_TT_lin_neg" + syskey]   = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([np.copy(b1_neg), np.copy(variances)], axis=-1))
                output_templates[cat + "/EWK_TT_quad_pos" + syskey]  = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([np.copy(b2_pos), np.copy(variances)], axis=-1))
                output_templates[cat + "/EWK_TT_quad_neg" + syskey]  = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([np.copy(b2_neg), np.copy(variances)], axis=-1))

            else:
                syskey = "_" + sys
                if sys.replace("Up", "").replace("Down", "") in systs_to_copy:
                    if "EWK_scheme" in sys:
                        if "Up" in sys:
                            coeff0, coeff1, coeff2 = np.copy(b0 * deltaQCD), np.copy(b1 * deltaQCD), np.copy(b2 * deltaQCD)
                        else:
                            coeff0, coeff1, coeff2 = np.copy(b0 / deltaQCD), np.copy(b1 / deltaQCD), np.copy(b2 / deltaQCD)
                    else:
                        nnlo_sys, coeff0, coeff1, coeff2 = (np.copy(template_nnloqcd[:, 0]), np.copy(b0), np.copy(b1), np.copy(b2)) * np.copy(sys_template) / np.copy(template_yukawa_nominal)
                        output_templates[cat + "/TT" + syskey] = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([np.copy(nnlo_sys), np.copy(variances)], axis=-1))

                    c0_pos = np.where( coeff0 >= 0.,  coeff0, 0. )
                    c0_neg = np.where( coeff0 <= 0., -coeff0, 0. )
                    c1_pos = np.where( coeff1 >= 0.,  coeff1, 0. )
                    c1_neg = np.where( coeff1 <= 0., -coeff1, 0. )
                    c2_pos = np.where( coeff2 >= 0.,  coeff2, 0. )
                    c2_neg = np.where( coeff2 <= 0., -coeff2, 0. )
                    output_templates[cat + "/EWK_TT_const_pos" + syskey] = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([np.copy(c0_pos), np.copy(variances)], axis=-1))
                    output_templates[cat + "/EWK_TT_const_neg" + syskey] = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([np.copy(c0_neg), np.copy(variances)], axis=-1))
                    output_templates[cat + "/EWK_TT_lin_pos" + syskey]   = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([np.copy(c1_pos), np.copy(variances)], axis=-1))
                    output_templates[cat + "/EWK_TT_lin_neg" + syskey]   = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([np.copy(c1_neg), np.copy(variances)], axis=-1))
                    output_templates[cat + "/EWK_TT_quad_pos" + syskey]  = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([np.copy(c2_pos), np.copy(variances)], axis=-1))
                    output_templates[cat + "/EWK_TT_quad_neg" + syskey]  = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([np.copy(c2_neg), np.copy(variances)], axis=-1))
                else:
                    if ("Up" in sys) or ("Down" in sys):
                        sys_vals_to_save = np.copy(sys_template - (b0 + b1 + b2))
                        output_templates[cat + "/TT" + syskey] = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([np.copy(sys_vals_to_save), np.copy(variances)], axis=-1))
                    else:
                        output_templates[cat + "/TT" + syskey] = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=sys_template)

    if args.format == "LL":
        for sys, np_stack in templates_to_copy[cat].items():
            ## hard coded ignoring of some templates
            if sys in TT_templates_to_ignore: continue
            print("Calculating templates for " + cat + " " + sys)
            if sys == "TT_NNLOqcd":
                output_templates[f"{cat}/TT"] = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np_stack)
            else:
                output_templates[f"{cat}/{sys}"] = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np_stack)

if args.plot:
    # matplotlib
    import matplotlib.pyplot as plt
    import mplhep as hep
    plt.style.use(hep.cms.style.ROOT)
    plt.switch_backend("agg")
    from matplotlib import rcParams
    rcParams["font.size"] = 20
    rcParams["savefig.format"] = "pdf"
    rcParams["savefig.bbox"] = "tight"

    import Utilities.Plotter as Plotter
    import Utilities.final_analysis_binning as final_binning
    import Utilities.common_features as cfeatures
    import os
    import warnings
    warnings.filterwarnings("ignore", "invalid value encountered in divide")

    plot_outdir = os.environ["plots_dir"]
    jobid = os.environ["jobid"]
    outdir = os.path.join(plot_outdir, jobid, "Yukawa_Coefficients")
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    linearize_binning = (
        final_binning.mtt_binning,
        final_binning.ctstar_abs_binning
    )
    ctstar_binlabels = [r"%s $\in [%s, %s)$" % (cfeatures.variable_names_to_labels["tlep_ctstar_abs"], linearize_binning[1][bin], linearize_binning[1][bin+1]) for bin in range(len(linearize_binning[1])-1)]
    ctstar_binlabels[-1] = r"%s $\in [%s, %s]$" % (cfeatures.variable_names_to_labels["tlep_ctstar_abs"], linearize_binning[1][-2], linearize_binning[1][-1])
    ctstar_bin_locs = np.linspace((len(linearize_binning[0])-1)/2, (len(linearize_binning[0])-1)*(len(linearize_binning[1])-1) - (len(linearize_binning[0])-1)/2, len(linearize_binning[1])-1)
    mtt_vals_to_plot = np.array([400, 600, 1000])
    mtt_tiled_labels = np.tile(linearize_binning[0][:-1], linearize_binning[1].size-1)
    mtt_bin_inds_to_plot = np.where(np.in1d(mtt_tiled_labels, mtt_vals_to_plot))[0]
    mtt_bins_to_plot = np.tile(mtt_vals_to_plot, linearize_binning[1].size-1)

    x_lims = (0, (linearize_binning[1].size - 1)* (linearize_binning[0].size - 1))

    styles_dict = {
        "EWK_TT_const" : {"label" : "$b_{0}$"},
        "EWK_TT_lin"   : {"label" : "$b_{1}$"},
        "EWK_TT_quad"  : {"label" : "$b_{2}$"},
    }

    for cat in templates.keys():
        channel, year =  cat.split("_")
        lep = "Electron" if channel.startswith("e") else "Muon"
        jmult = "3Jets" if "3" in channel else "4PJets"
        pltdir = os.path.join(outdir, cat.replace("_", "/"))
        if not os.path.isdir(pltdir):
            os.makedirs(pltdir)

        for part in ["EWK_TT_const", "EWK_TT_lin", "EWK_TT_quad"]:

            fig, ax = plt.subplots(figsize=(15.0, 10.0))
            fig.subplots_adjust(hspace=.07)

            Plotter.plot_1D(output_templates[f"{cat}/{part}"].values(), output_templates[f"{cat}/{part}"].axes.edges[0],
                xlimits=x_lims, color="k", ax=ax, label=styles_dict[part]["label"])
            ax.legend(loc="upper right", title="Coefficient", ncol=1)
            ax.axhline(0, **{"linestyle": "--", "color": (0, 0, 0, 0.5), "linewidth": 1})
            ax.autoscale()
            ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.3)
            ax.set_xlabel("$m_{t\\bar{t}}$ [GeV]")
            ax.set_xlim(0, output_templates[f"{cat}/{part}"].axes.centers[0].size)

                # add lepton/jet multiplicity label
            ax.text(
                0.02, 0.92, cfeatures.channel_labels[f"{lep}_{jmult}"],
                fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
            )
                ## draw vertical lines for distinguishing different ctstar bins
            vlines = [x_lims[1]*ybin/5 for ybin in range(1, 5)]
            [ax.axvline(vline, color="k", linestyle="--") for vline in vlines]

            for idx, label in enumerate(ctstar_binlabels):
                ax.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                    xytext=(0, 250), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)

            ax.set_xticks(mtt_bin_inds_to_plot)
            ax.set_xticklabels(mtt_bins_to_plot)
            hep.cms.label(ax=ax, data=False, year=year)

            figname = os.path.join(pltdir, "_".join([cat, part, "Coeff"]))
            fig.savefig(figname)
            print(f"{figname} written")
            plt.close(fig)


import time
named_tuple = time.localtime() # get struct_time
time_str = time.strftime("%d%m%Y", named_tuple)
output_rname = "_".join([(args.outputfile).split(".root")[0], f"For{args.format}", time_str])+".root" ## append date onto root filename
output_rfile = uproot.recreate(output_rname, compression=uproot.ZLIB(4))
for key, hist in output_templates.items():
    channel, sys = key.split("/")
    if sys in TT_templates_to_ignore: continue
    print("Calculating templates for " + channel + " " + sys)
    if sys == "TT_NNLOqcd":
        output_rfile[f"{channel}/TT"] = hist.copy()
    else:
        output_rfile[key] = hist.copy()
output_rfile.close()

print(f"\n\nFinished writing {output_rname} to disk")

toc = time.time()
print("Total time: %.1f" % (toc - tic))
