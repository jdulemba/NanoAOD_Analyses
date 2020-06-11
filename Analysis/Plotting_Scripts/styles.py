## for good colors search here http://colorbrewer2.org/#type=qualitative&scheme=Set1&n=6

styles = {
#   '[WZ][WZ]' : {
#      'legendstyle' : 'f',
#      'drawstyle' : 'hist',
#      'facecolor' : '#9fcc18',
#      'linecolor' : 'black',
#            'linewidth' : 1,
#      'name' : "VV",
#      #'name' : "diboson",
#      'fillstyle': 'solid',
#      },
    'ZJets*' : {
        'facecolor' : '#984ea3', ## purple
        'name' : "Z+jets",
    },
    'single*' : {
        'facecolor' : '#ff7f00', ## orange
        'name' : "single top",
    },
   #'data*' : {
   #   'legendstyle' : 'p',
   #   'drawstyle' : 'E0 X0',
   #   'markerstyle' : 20,
   #   'name' : "Observed",
   # },
   'WJets*' : {
        'facecolor' : '#ff7f00', ## orange
        'name' : "W+jets",
    },
    'tt*' : {
        'facecolor' : '#e41a1c', ## red
        'name' : "$t\\bart$",
    },
#   'tt[WZ]*' : {
#      'legendstyle' : 'f',
#      'drawstyle' : 'hist',
#      'facecolor' : '#cc8d18',
#      #'facecolor' : ROOT.kOrange + 1,
#      'linecolor' : 'black',
#            'linewidth' : 1,
#      'name' : "ttV",
#      'fillstyle': 'solid',
#      },
    'QCD*' : {
        'facecolor' : '#377eb8', ## blue
        'name' : "QCD",
    },
    'EWK' : {
        'facecolor' : '#4daf4a', ## green
        'name' : "EW",
    },
        ## default name for coffea statistical errors
    'Sum unc.' : {
        'facecolor' : 'none',
        'linewidth' : 0,
        'name' : 'Stat. Unc.',
        'hatch' : '///',
    },
    'data_err_opts' : {
        'marker': '.',
        'markersize': 10.,
        'color':'k',
        'elinewidth': 1,
    },

        ## permutation categories
    'ttJets_right' : {
        'facecolor' : '#e41a1c', ## red
        'name' : "$t\\bart$ correct",
    },
    'ttJets_matchable' : {
        'facecolor' : '#984ea3', ## purple
        'name' : "$t\\bart$ matchable",
    },
    'ttJets_unmatchable' : {
        'facecolor' : '#ffff33', ## yellow
        'name' : "$t\\bart$ unmatchable",
    },
    'ttJets_other' : {
        'facecolor' : '#a65628', ## brown
        'name' : "$t\\bart$ other",
    },
    'Correct_BLep' : {
        'name' : 'Correct $b_{l}$',
        'linestyle' : '-',
        'color' : 'r',
        'elinewidth' : 1,
    },
    'Wrong_BLep' : {
        'name' : 'Wrong $b_{l}$',
        'linestyle' : '-',
        'color' : 'b',
        'elinewidth' : 1,
    },
    'Correct_THad' : {
        'name' : 'Correct $t_{h}$',
        'linestyle' : '-',
        'color' : 'r',
        'elinewidth' : 1,
    },
    'Wrong_THad' : {
        'name' : 'Wrong $t_{h}$',
        'linestyle' : '-',
        'color' : 'b',
        'elinewidth' : 1,
    },
    'Correct' : {
        'name' : 'Correct',
        'linestyle' : '-',
        'color' : 'r',
        'elinewidth' : 1,
    },
    'Wrong' : {
        'name' : 'Wrong',
        'linestyle' : '-',
        'color' : 'b',
        'elinewidth' : 1,
    },
    '3Jets' : {
        #'name' : '3 jets',
        'marker': '.',
        'markersize': 10.,
        #'linestyle' : '-',
        'color' : '#e41a1c', # red
        'elinewidth' : 1,
    },
    '4PJets' : {
        #'name' : '4+ jets',
        'marker': '.',
        'markersize': 10.,
        #'linestyle' : '-',
        'color' : '#377eb8', # blue
        'elinewidth' : 1,
    },
    '3PJets' : {
        #'name' : '3+ jets',
        'marker': '.',
        'markersize': 10.,
        #'linestyle' : '-',
        'color' : '#4daf4a', # green
        'elinewidth' : 1,
    },
    'iso_reg_0' : {
        'marker': '.',
        'markersize': 10.,
        'color' : '#e41a1c', # red
        'elinewidth' : 1,
    },
    'iso_reg_1' : {
        'marker': '.',
        'markersize': 10.,
        'color' : '#377eb8', # blue
        'elinewidth' : 1,
    },
    'iso_reg_2' : {
        'marker': '.',
        'markersize': 10.,
        'color' : '#4daf4a', # green
        'elinewidth' : 1,
    }
}
