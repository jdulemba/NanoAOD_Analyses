import uproot
import argparse
import numpy as np
from collections import defaultdict
import hist
from pdb import set_trace

parser = argparse.ArgumentParser()
parser.add_argument("inputtemplates")
parser.add_argument("outputfile")
parser.add_argument("--plot", action="store_true", help="Make plots of sys uncs")
args = parser.parse_args()

proc = "TT"
ewk_name = "EWK_yukawa"
nnlo_name = "NNLOqcd"

channels = ["em", "ee", "mm", "e3jets", "mu3jets", "e4pjets", "mu4pjets"]
years = ["2016pre", "2016post", "2017", "2018"]

# Weight-based systematics to copy to the EWK templates
# It is assumed that EWK_yukawa and the other nuisance are uncorellated. Since they are both
# weight-based, the template for both systematics varied by 1sigma is then just
# (EW_yukawa up/down) * (other up/down) / (nominal)
# Then, the same quadratic equation is solved for the systematic as for the nominal to get
# the variation templates for the EWK signal.

systs_to_copy = [
    #"EWK_scheme",
    #"CMS_eff_e_reco",
    #"CMS_eff_e_id",
    #"CMS_eff_m_id_stat",
    #"CMS_eff_m_id_syst",
    #"CMS_eff_m_iso_stat",
    #"CMS_eff_m_iso_syst",
    #"CMS_eff_trigger_ee",
    #"CMS_eff_trigger_em",
    #"CMS_eff_trigger_mm",
    #"CMS_L1_prefire",
    #"CMS_eff_trigger_m_syst",
    #"CMS_eff_trigger_m_stat",
    #"CMS_eff_trigger_e",
]

dyt_up = 0.11
dyt_dw = -0.12

yt_nom = 1.0
yt_up = yt_nom + dyt_up
yt_dw = yt_nom + dyt_dw

templates = defaultdict(dict)
templates_yukawa_up = {}
templates_yukawa_down = {}
templates_nnloqcd = {}

#set_trace()
with uproot.open(args.inputtemplates) as rfile:
    for channel in channels:
        for year in years:
            cat = channel + "_" + year
            print("Loading templates for " + cat)
            if cat in rfile:
                rcat = rfile[cat]
                templates[cat]["nominal"] = rcat[proc].values()
                templates_yukawa_up[cat] = rcat[proc + "_" + ewk_name + "Up"].values()
                templates_yukawa_down[cat] = rcat[proc + "_" + ewk_name + "Down"].values()
                templates_nnloqcd[cat] = rcat[proc + "_" + nnlo_name].values()
                for sys in systs_to_copy:
                    for sysdir in ["Up", "Down"]:
                        syskey = proc + "_" + sys + sysdir
                        if syskey in rcat:
                            templates[cat][sys + sysdir] = rcat[syskey].values()

"""
The yields for TTbar distributions at a given yt value are given by 
 TT(yt) = TT_NNLOqcd + TT_ewk(yt) == TT_NNLOqcd + C_ewk + A * (yt) + B * (yt)^2
where TT_NNLOqcd is the contribution from including NNLOqcd terms via the NNLO QCD correction and 
 TT_ewk(yt) is the contribution from the NLO EWK terms via the EWK corrections

 TT = TT_NNLOqcd + Cew + A * (yt_nom) + B * (yt_nom)^2
 TT_YukawaUp = TT_NNLOqcd + Cew + A * (yt_up) + B * (yt_up)^2
 TT_YukawaDown = TT_NNLOqcd + Cew + A * (yt_dw) + B * (yt_dw)^2


We can solve for Cew, A, and B through the matrices Ax = B

   | 1.00   yt_nom   yt_nom^2 |   | Cew | = | TT            - TT_NNLOqcd |
   | 1.00   yt_up    yt_up^2  |   | A    | = | TT_YukawaUp   - TT_NNLOqcd |
   | 1.00   yt_dw    yt_dw^2  |   | B    | = | TT_YukawaDown - TT_NNLOqcd |

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

#set_trace()
for cat in templates.keys():
    for sys in templates[cat].keys():
        print("Calculating templates for " + cat + " " + sys)

        template_yukawa_nominal = np.copy(templates[cat][sys])
        template_yukawa_up      = np.copy(templates_yukawa_up[cat] * templates[cat][sys] / templates[cat]["nominal"])
        template_yukawa_down    = np.copy(templates_yukawa_down[cat] * templates[cat][sys] / templates[cat]["nominal"])
        template_nnloqcd        = np.copy(templates_nnloqcd[cat])

            ## find the terms that the ewk corrections are responsible for, dependent on yt (RHS of matrix equation)
        ewk_dists = [[template_yukawa_nominal[bin] - template_nnloqcd[bin], template_yukawa_up[bin] - template_nnloqcd[bin], template_yukawa_down[bin] - template_nnloqcd[bin]] 
            for bin in range(template_yukawa_nominal.size)]


            # solve for EWK term coefficients
        try:
            Cew, Aew, Bew = np.array([np.linalg.solve(coeffs_matrix, np.array(bin)) for bin in ewk_dists]).T # transpose to make resulting vector size (nbins x 3)
        except LinAlgError:
            raise ValueError("Unable to solve for coefficients!")
            #x = np.linalg.lstsq(coeffs_matrix, ewk_dists[bin])[0]

        # Check that this reproduces the input templates
        def predict_yukawa(yt):
            return template_nnloqcd + Cew + Aew * yt + Bew * yt**2

        assert np.all(np.isclose(predict_yukawa(yt_nom), template_yukawa_nominal)) # nominal
        assert np.all(np.isclose(predict_yukawa(yt_up), template_yukawa_up)) # up
        assert np.all(np.isclose(predict_yukawa(yt_dw), template_yukawa_down)) # down

        #set_trace()

        Cew_pos = np.where(Cew>=0 , Cew, 0.)
        Cew_neg = np.where(Cew<=0 , -Cew, 0.)
        Aew_pos = np.where(Aew>=0 , Aew, 0.)
        Aew_neg = np.where(Aew<=0 , -Aew, 0.)
        Bew_pos = np.where(Bew>=0 , Bew, 0.)
        Bew_neg = np.where(Bew<=0 , -Bew, 0.)

        # Add bin edges so that uproot saves it as a TH1D
        bin_edges = np.arange(len(template_yukawa_nominal)+1)

        syskey = "" if sys == "nominal" else "_" + sys

        variances = np.zeros_like(Cew)

        hist_Cew     = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([Cew, variances], axis=-1))
        hist_Cew_pos = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([Cew_pos, variances], axis=-1))
        hist_Cew_neg = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([Cew_neg, variances], axis=-1))
        hist_Aew     = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([Aew, variances], axis=-1))
        hist_Aew_pos = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([Aew_pos, variances], axis=-1))
        hist_Aew_neg = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([Aew_neg, variances], axis=-1))
        hist_Bew     = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([Bew, variances], axis=-1))
        hist_Bew_pos = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([Bew_pos, variances], axis=-1))
        hist_Bew_neg = hist.Hist(hist.axis.Variable(bin_edges), "Weight", data=np.stack([Bew_neg, variances], axis=-1))

        output_templates[cat + "/EWK_TT_const" + syskey]     = hist_Cew
        output_templates[cat + "/EWK_TT_const_pos" + syskey] = hist_Cew_pos
        output_templates[cat + "/EWK_TT_const_neg" + syskey] = hist_Cew_neg
        output_templates[cat + "/EWK_TT_lin" + syskey]       = hist_Aew
        output_templates[cat + "/EWK_TT_lin_pos" + syskey]   = hist_Aew_pos
        output_templates[cat + "/EWK_TT_lin_neg" + syskey]   = hist_Aew_neg
        output_templates[cat + "/EWK_TT_quad" + syskey]      = hist_Bew
        output_templates[cat + "/EWK_TT_quad_pos" + syskey]  = hist_Bew_pos
        output_templates[cat + "/EWK_TT_quad_neg" + syskey]  = hist_Bew_neg


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
        "EWK_TT_const" : {"label" : "$C_{EW}$"},
        "EWK_TT_lin"   : {"label" : "$A_{EW}$"},
        "EWK_TT_quad"  : {"label" : "$B_{EW}$"},
    }

    for cat in templates.keys():
        #set_trace()
        channel, year =  cat.split("_")
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
            ax.set_xlabel("$m_{t\\bar{t}}$ [GeV]")

                # add lepton/jet multiplicity label
            ax.text(
                0.02, 0.92, channel,
                fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
            )
                ## draw vertical lines for distinguishing different ctstar bins
            vlines = [x_lims[1]*ybin/5 for ybin in range(1, 5)]
            [ax.axvline(vline, color="k", linestyle="--") for vline in vlines]

            for idx, label in enumerate(ctstar_binlabels):
                ax.annotate(label, xy=(ctstar_bin_locs[idx], 0), xycoords=("data", "axes fraction"),
                    xytext=(0, 25), textcoords="offset points", va="bottom", ha="center", fontsize=rcParams["font.size"]*0.70, rotation=0)

            ax.set_xticks(mtt_bin_inds_to_plot)
            ax.set_xticklabels(mtt_bins_to_plot)
            hep.cms.label(ax=ax, data=False, year=year)

            figname = os.path.join(pltdir, "_".join([cat, part, "Coeff"]))
            fig.savefig(figname)
            print(f"{figname} written")
            plt.close(fig)

with uproot.recreate(args.outputfile) as rfile:
    for key, hist in output_templates.items():
        rfile[key] = hist

print(f"\n\nFinished writing {args.outputfile} to disk")
