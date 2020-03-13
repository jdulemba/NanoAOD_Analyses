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
        'facecolor' : '#4daf4a', ## green
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
        'facecolor' : '#ff7f00', ## orange
        'name' : "EW",
    },
        ## default name for coffea statistical errors
    'Sum unc.' : {
        'facecolor' : 'none',
        'linewidth' : 0,
        'name' : 'Stat. Unc.',
        'hatch' : '///',
    }

}
