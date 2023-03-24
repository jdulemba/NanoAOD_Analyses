## for good colors search here http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6

styles = {
#   "[WZ][WZ]" : {
#      "legendstyle" : "f",
#      "drawstyle" : "hist",
#      "facecolor" : "#9fcc18",
#      "linecolor" : "black",
#            "linewidth" : 1,
#      "name" : "VV",
#      #"name" : "diboson",
#      "fillstyle": "solid",
#      },
    "ZJets" : {
        "facecolor" : "#984ea3", ## purple
        "name" : "Z+jets",
    },
    "singlet" : {
        "facecolor" : "#ff7f00", ## orange
        "name" : "single top",
    },
   #"data" : {
   #   "legendstyle" : "p",
   #   "drawstyle" : "E0 X0",
   #   "markerstyle" : 20,
   #   "name" : "Observed",
   # },
   "WJets" : {
        "facecolor" : "#ff7f00", ## orange
        "name" : "W+jets",
    },
    "ttJets" : {
        "facecolor" : "#e41a1c", ## red
        "name" : "$t\\bar{t}$",
    },
    "tt" : {
        "facecolor" : "#e41a1c", ## red
        "name" : "$t\\bar{t}$",
    },
#   "tt[WZ]" : {
#      "legendstyle" : "f",
#      "drawstyle" : "hist",
#      "facecolor" : "#cc8d18",
#      #"facecolor" : ROOT.kOrange + 1,
#      "linecolor" : "black",
#            "linewidth" : 1,
#      "name" : "ttV",
#      "fillstyle": "solid",
#      },
    "QCD" : {
        "facecolor" : "#377eb8", ## blue
        "name" : "QCD (Sim.)",
    },
    "EWK" : {
        "facecolor" : "#4daf4a", ## green
        "name" : "EW (Sim.)",
    },
    "BKG" : {
        "facecolor" : "#4daf4a", ## green
        "name" : "BKG (Est.)",
    },
    "EWQCD" : {
        "facecolor" : "#4daf4a", ## green
        "name" : "EWQCD (Data-Driven)",
    },
        ## default name for coffea statistical errors
    "Sum unc." : {
        "facecolor" : "none",
        "linewidth" : 0,
        "name" : "Stat. Unc.",
        "hatch" : "///",
    },
    "Total MEscale Unc." : {
        "facecolor" : "none",
        "linewidth" : 0,
        "name" : "$\mu_{R}$+$\mu_{F}$ Unc.",
        "hatch" : "///",
    },
        ## default name for coffea statistical errors
    "Total Unc." : {
        "facecolor" : "none",
        "linewidth" : 0,
        "name" : "Stat.+Syst. Unc.",
        "hatch" : "///",
    },
    "data_err_opts" : {
        "marker": ".",
        "markersize": 10.,
        "color":"k",
        "elinewidth": 1,
    },

        ## permutation categories
    "ttJets_right" : {
        "facecolor" : "#e41a1c", ## red
        "name" : "$t\\bar{t}$ correct",
    },
    "ttJets_matchable" : {
        "facecolor" : "#984ea3", ## purple
        "name" : "$t\\bar{t}$ matchable",
    },
    "ttJets_unmatchable" : {
        "facecolor" : "#ffff33", ## yellow
        "name" : "$t\\bar{t}$ unmatchable",
    },
    "ttJets_sl_tau" : {
        "facecolor" : "b",
        "name" : "$t\\bar{t}$, SL $\\tau$",
    },
    "ttJets_other" : {
        "facecolor" : "#a65628", ## brown
        "name" : "$t\\bar{t}$ other",
    },
    "Correct_BLep" : {
        "name" : "Correct $b_{\ell}$",
        "linestyle" : "-",
        "color" : "r",
        "elinewidth" : 1,
    },
    "Wrong_BLep" : {
        "name" : "Wrong $b_{\ell}$",
        "linestyle" : "-",
        "color" : "b",
        "elinewidth" : 1,
    },
    "Correct_THad" : {
        "name" : "Correct $t_{h}$",
        "linestyle" : "-",
        "color" : "r",
        "elinewidth" : 1,
    },
    "Wrong_THad" : {
        "name" : "Wrong $t_{h}$",
        "linestyle" : "-",
        "color" : "b",
        "elinewidth" : 1,
    },
    "Correct" : {
        "name" : "Correct",
        "linestyle" : "-",
        "color" : "r",
        "elinewidth" : 1,
    },
    "Wrong" : {
        "name" : "Wrong",
        "linestyle" : "-",
        "color" : "b",
        "elinewidth" : 1,
    },
    "3Jets" : {
        #"name" : "3 jets",
        "marker": ".",
        "markersize": 10.,
        #"linestyle" : "-",
        "color" : "#e41a1c", # red
        "elinewidth" : 1,
    },
    "4PJets" : {
        #"name" : "4+ jets",
        "marker": ".",
        "markersize": 10.,
        #"linestyle" : "-",
        "color" : "#377eb8", # blue
        "elinewidth" : 1,
    },
    "3PJets" : {
        #"name" : "3+ jets",
        "marker": ".",
        "markersize": 10.,
        #"linestyle" : "-",
        "color" : "#4daf4a", # green
        "elinewidth" : 1,
    },
    "iso_reg_0" : {
        "marker": ".",
        "markersize": 10.,
        "color" : "#e41a1c", # red
        "elinewidth" : 1,
    },
    "iso_reg_1" : {
        "marker": ".",
        "markersize": 10.,
        "color" : "#377eb8", # blue
        "elinewidth" : 1,
    },
    "iso_reg_2" : {
        "marker": ".",
        "markersize": 10.,
        "color" : "#4daf4a", # green
        "elinewidth" : 1,
    },

    "AtoTT" : {
        "color" : "#377eb8", # blue
        "linewidth" : 2,
        "name" : "$\mathsf{A \\rightarrow t\\bar {t}}$",
    },
    "HtoTT" : {
        "color" : "#e41a1c", # red
        "linewidth" : 2,
        "name" : "$\mathsf{H \\rightarrow t\\bar {t}}$",
    },
    "M365" : { 
        "color" : "#e41a1c", # red
        "linewidth" : 2,
        "name" : "m=365 GeV",
    },
    "M400" : { 
        "color" : "#377eb8", # blue
        "linewidth" : 2,
        "name" : "m=400 GeV",
    },
    "M500" : { 
        "color" : "#4daf4a", # green
        "linewidth" : 2,
        "name" : "m=500 GeV",
    },
    "M600" : { 
        "color" : "#984ea3", # purple
        "linewidth" : 2,
        "name" : "m=600 GeV",
    },
    "M800" : { 
        "color" : "#ff7f00", # orange
        "linewidth" : 2,
        "name" : "m=800 GeV",
    },
    "M1000" : { 
        "color" : "k", # black
        "linewidth" : 2,
        "name" : "m=1000 GeV",
    },
    "W2p5" : { 
        "color" : "#e41a1c", # red
        "linewidth" : 2,
        "name" : "$\mathsf{\Gamma}$/m=2.5%",
    },
    "W5" : { 
        "color" : "#377eb8", # blue
        "linewidth" : 2,
        "name" : "$\mathsf{\Gamma}$/m=5%",
    },
    "W10" : { 
        "color" : "k", # black
        "linewidth" : 2,
        "name" : "$\mathsf{\Gamma}$/m=10%",
    },
    "W25" : { 
        "color" : "#4daf4a", # green
        "linewidth" : 2,
        "name" : "$\mathsf{\Gamma}$/m=25%",
    },

        ## systematic variations
    "Nominal" : {
        "color" : "k",
        "name" : "Nominal",
    },
    "Up" : {
        "color" : "k",
        "linestyle" : "-",
        #"color" : "#e41a1c", # red
        "label" : "Up",
    },
    "Down" : {
        "color" : "k", # blue
        #"color" : "#377eb8", # blue
        "linestyle" : "--",
        "label" : "Down",
    },


        ## ttbar decays
    "ttJetsSL" : {
        "facecolor" : "k",
        "name" : "$\mathrm{t\\bar{t}}_{\ell j}$",
    },
    "ttJetsDiLep" : {
        "facecolor" : "k",
        "name" : "$\mathrm{t\\bar{t}}_{\ell \ell}$",
    },
    "ttJetsHad" : {
        "facecolor" : "k",
        "name" : "$\mathrm{t\\bar{t}}_{j j}$",
    },

    "SL e" : {
        "facecolor" : "#FFFF33", # yellow
        "name" : "$e$",
    },
    "SL mu" : {
        "facecolor" : "#FF8000", # orange
        "name" : "$\mu$",
    },
    "SL tau->l" : {
        "facecolor" : "#FF3333", # red
        "name" : "$\\tau \\rightarrow l$",
    },
    "SL tau->h" : {
        "facecolor" : "#990000", # dark red
        "name" : "$\\tau \\rightarrow h$",
    },
    "Had Total" : {
        "facecolor" : "#b15928", # brown
        "name" : "All Had",
    },
    "DL e e" : {
        "facecolor" : "#4C9900", # dark green
        "name" : "$e e$",
    },
    "DL e mu" : {
        "facecolor" : "#00FF00", # green
        "name" : "$e \mu$",
    },
    "DL e tau->l" : {
        "facecolor" : "#003366", # dark blue
        "name" : "$e \\tau \\rightarrow l$",
    },
    "DL e tau->h" : {
        "facecolor" : "#0000FF", # blue
        "name" : "$e \\tau \\rightarrow h$",
    },
    "DL mu mu" : {
        "facecolor" : "#33FFFF", # light blue
        "name" : "$\mu \mu$",
    },
    "DL mu tau->l" : {
        "facecolor" : "#4C0099", # purple
        "name" : "$\mu \\tau \\rightarrow l$",
    },
    "DL mu tau->h" : {
        "facecolor" : "#9933FF", # light purple
        "name" : "$\mu \\tau \\rightarrow h$",
    },
    "DL tau tau->ll" : {
        "facecolor" : "#990099", # violet
        "name" : "$\\tau \\tau \\rightarrow ll$",
    },
    "DL tau tau->lh" : {
        "facecolor" : "#FF66FF", # pink
        "name" : "$\\tau \\tau \\rightarrow lh$",
    },
    "DL tau tau->hh" : {
        "facecolor" : "#606060", # grey
        "name" : "$\\tau \\tau \\rightarrow hh$",
    },
}


from itertools import product
import numpy as np
masses = np.arange(365., 1005., 5.)
widths = np.arange(0.5, 25.5, 0.5)
for bundle in product(masses, [str(width).replace(".", "p") for width in sorted(widths)]):

    # All ttbar decays events
        # individual positive, negative weights for interference
    styles["AtoTT_M%d_W%s_Int_neg" % bundle] = {
        "facecolor" : "#e41a1c",
        "name" : "$A_{%d\ GeV}^{%s\%%}$, Int, w$<$0" % (bundle[0], bundle[1].replace("p", ".")),
    }
    styles["AtoTT_M%d_W%s_Int_pos" % bundle] = {
        "facecolor" : "#e41a1c",
        "name" : "$A_{%d\ GeV}^{%s\%%}$, Int, w$>$0" % (bundle[0], bundle[1].replace("p", ".")),
    }
        # combined Interference
    styles["AtoTT_M%d_W%s_Int" % bundle] = {
        "facecolor" : "#e41a1c",
        "name" : "$A_{%d\ GeV}^{%s\%%}$, Int" % (bundle[0], bundle[1].replace("p", ".")),
    }
    styles["AtoTT_M%d_W%s_Res" % bundle] = {
        "facecolor" : "#e41a1c",
        "name" : "$A_{%d\ GeV}^{%s\%%}$, Res" % (bundle[0], bundle[1].replace("p", ".")),
    }
        # individual positive, negative weights for interference
    styles["HtoTT_M%d_W%s_Int_neg" % bundle] = {
        "facecolor" : "#e41a1c",
        "name" : "$H_{%d\ GeV}^{%s\%%}$, Int, w$<$0" % (bundle[0], bundle[1].replace("p", ".")),
    }
    styles["HtoTT_M%d_W%s_Int_pos" % bundle] = {
        "facecolor" : "#e41a1c",
        "name" : "$H_{%d\ GeV}^{%s\%%}$, Int, w$>$0" % (bundle[0], bundle[1].replace("p", ".")),
    }
        # combined Interference
    styles["HtoTT_M%d_W%s_Int" % bundle] = {
        "facecolor" : "#e41a1c",
        "name" : "$H_{%d\ GeV}^{%s\%%}$, Int" % (bundle[0], bundle[1].replace("p", ".")),
    }
    styles["HtoTT_M%d_W%s_Res" % bundle] = {
        "facecolor" : "#e41a1c",
        "name" : "$H_{%d\ GeV}^{%s\%%}$, Res" % (bundle[0], bundle[1].replace("p", ".")),
    }
        # total signal res+int
    styles["AtoTT_M%d_W%s_Total" % bundle] = {
        "facecolor" : "#e41a1c",
        "name" : "$A_{%d\ GeV}^{%s\%%}$, Total" % (bundle[0], bundle[1].replace("p", ".")),
    }
    styles["HtoTT_M%d_W%s_Total" % bundle] = {
        "facecolor" : "#e41a1c",
        "name" : "$H_{%d\ GeV}^{%s\%%}$, Total" % (bundle[0], bundle[1].replace("p", ".")),
    }

    # DiLep events
        # individual positive, negative weights for interference
    styles["AtoTTJetsDiLep_M%d_W%s_Int_neg" % bundle] = {
        "facecolor" : "#e41a1c",
        "name" : "$A_{%d\ GeV}^{%s\%%}$ $\\rightarrow \mathrm{t\\bar {t}}_{\ell \ell}$, Int, w$<$0" % (bundle[0], bundle[1].replace("p", ".")),
    }
    styles["AtoTTJetsDiLep_M%d_W%s_Int_pos" % bundle] = {
        "facecolor" : "#e41a1c",
        "name" : "$A_{%d\ GeV}^{%s\%%}$ $\\rightarrow \mathrm{t\\bar {t}}_{\ell \ell}$, Int, w$>$0" % (bundle[0], bundle[1].replace("p", ".")),
    }
        # combined Interference
    styles["AtoTTJetsDiLep_M%d_W%s_Int" % bundle] = {
        "facecolor" : "#e41a1c",
        "name" : "$A_{%d\ GeV}^{%s\%%}$ $\\rightarrow \mathrm{t\\bar {t}}_{\ell \ell}$, Int" % (bundle[0], bundle[1].replace("p", ".")),
    }
    styles["AtoTTJetsDiLep_M%d_W%s_Res" % bundle] = {
        "facecolor" : "#e41a1c",
        "name" : "$A_{%d\ GeV}^{%s\%%}$ $\\rightarrow \mathrm{t\\bar {t}}_{\ell \ell}$, Res" % (bundle[0], bundle[1].replace("p", ".")),
    }
        # individual positive, negative weights for interference
    styles["HtoTTJetsDiLep_M%d_W%s_Int_neg" % bundle] = {
        "facecolor" : "#e41a1c",
        "name" : "$H_{%d\ GeV}^{%s\%%}$ $\\rightarrow \mathrm{t\\bar {t}}_{\ell \ell}$, Int, w$<$0" % (bundle[0], bundle[1].replace("p", ".")),
    }
    styles["HtoTTJetsDiLep_M%d_W%s_Int_pos" % bundle] = {
        "facecolor" : "#e41a1c",
        "name" : "$H_{%d\ GeV}^{%s\%%}$ $\\rightarrow \mathrm{t\\bar {t}}_{\ell \ell}$, Int, w$>$0" % (bundle[0], bundle[1].replace("p", ".")),
    }
        # combined Interference
    styles["HtoTTJetsDiLep_M%d_W%s_Int" % bundle] = {
        "facecolor" : "#e41a1c",
        "name" : "$H_{%d\ GeV}^{%s\%%}$ $\\rightarrow \mathrm{t\\bar {t}}_{\ell \ell}$, Int" % (bundle[0], bundle[1].replace("p", ".")),
    }
    styles["HtoTTJetsDiLep_M%d_W%s_Res" % bundle] = {
        "facecolor" : "#e41a1c",
        "name" : "$H_{%d\ GeV}^{%s\%%}$ $\\rightarrow \mathrm{t\\bar {t}}_{\ell \ell}$, Res" % (bundle[0], bundle[1].replace("p", ".")),
    }
        # total signal res+int
    styles["AtoTTJetsDiLep_M%d_W%s_Total" % bundle] = {
        "facecolor" : "#e41a1c",
        "name" : "$A_{%d\ GeV}^{%s\%%}$ $\\rightarrow \mathrm{t\\bar {t}}_{\ell \ell}$, Total" % (bundle[0], bundle[1].replace("p", ".")),
    }
    styles["HtoTTJetsDiLep_M%d_W%s_Total" % bundle] = {
        "facecolor" : "#e41a1c",
        "name" : "$H_{%d\ GeV}^{%s\%%}$ $\\rightarrow \mathrm{t\\bar {t}}_{\ell \ell}$, Total" % (bundle[0], bundle[1].replace("p", ".")),
    }

    # Semilep events
        # individual positive, negative weights for interference
    styles["AtoTTJetsSL_M%d_W%s_Int_neg" % bundle] = {
        "facecolor" : "#e41a1c",
        "name" : "$A_{%d\ GeV}^{%s\%%}$ $\\rightarrow \mathrm{t\\bar {t}}_{\ell j}$, Int, w$<$0" % (bundle[0], bundle[1].replace("p", ".")),
    }
    styles["AtoTTJetsSL_M%d_W%s_Int_pos" % bundle] = {
        "facecolor" : "#e41a1c",
        "name" : "$A_{%d\ GeV}^{%s\%%}$ $\\rightarrow \mathrm{t\\bar {t}}_{\ell j}$, Int, w$>$0" % (bundle[0], bundle[1].replace("p", ".")),
    }
        # combined Interference
    styles["AtoTTJetsSL_M%d_W%s_Int" % bundle] = {
        "facecolor" : "#e41a1c",
        "name" : "$A_{%d\ GeV}^{%s\%%}$ $\\rightarrow \mathrm{t\\bar {t}}_{\ell j}$, Int" % (bundle[0], bundle[1].replace("p", ".")),
    }
    styles["AtoTTJetsSL_M%d_W%s_Res" % bundle] = {
        "facecolor" : "#e41a1c",
        "name" : "$A_{%d\ GeV}^{%s\%%}$ $\\rightarrow \mathrm{t\\bar {t}}_{\ell j}$, Res" % (bundle[0], bundle[1].replace("p", ".")),
    }

        # individual positive, negative weights for interference
    styles["HtoTTJetsSL_M%d_W%s_Int_neg" % bundle] = {
        "facecolor" : "#e41a1c",
        "name" : "$H_{%d\ GeV}^{%s\%%}$ $\\rightarrow \mathrm{t\\bar {t}}_{\ell j}$, Int, w$<$0" % (bundle[0], bundle[1].replace("p", ".")),
    }
    styles["HtoTTJetsSL_M%d_W%s_Int_pos" % bundle] = {
        "facecolor" : "#e41a1c",
        "name" : "$H_{%d\ GeV}^{%s\%%}$ $\\rightarrow \mathrm{t\\bar {t}}_{\ell j}$, Int, w$>$0" % (bundle[0], bundle[1].replace("p", ".")),
    }
        # combined Interference
    styles["HtoTTJetsSL_M%d_W%s_Int" % bundle] = {
        "facecolor" : "#e41a1c",
        "name" : "$H_{%d\ GeV}^{%s\%%}$ $\\rightarrow \mathrm{t\\bar {t}}_{\ell j}$, Int" % (bundle[0], bundle[1].replace("p", ".")),
    }
    styles["HtoTTJetsSL_M%d_W%s_Res" % bundle] = {
        "facecolor" : "#e41a1c",
        "name" : "$H_{%d\ GeV}^{%s\%%}$ $\\rightarrow \mathrm{t\\bar {t}}_{\ell j}$, Res" % (bundle[0], bundle[1].replace("p", ".")),
    }
        # total signal res+int
    styles["AtoTTJetsSL_M%d_W%s_Total" % bundle] = {
        "facecolor" : "#e41a1c",
        "name" : "$A_{%d\ GeV}^{%s\%%}$ $\\rightarrow \mathrm{t\\bar {t}}_{\ell j}$, Total" % (bundle[0], bundle[1].replace("p", ".")),
    }
    styles["HtoTTJetsSL_M%d_W%s_Total" % bundle] = {
        "facecolor" : "#e41a1c",
        "name" : "$H_{%d\ GeV}^{%s\%%}$ $\\rightarrow \mathrm{t\\bar {t}}_{\ell j}$, Total" % (bundle[0], bundle[1].replace("p", ".")),
    }
