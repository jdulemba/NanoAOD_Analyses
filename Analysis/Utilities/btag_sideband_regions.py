btag_sb_region_boundaries = {
    "2016APV" : {
        "DeepCSVMedium" : {
            "p00p20" : (0.00, 0.20),
            "p20p40" : (0.20, 0.40),
            "p40p60" : (0.40, 0.60),
        },
        "DeepJetMedium" : {
            "p00p0866" : (0.00, 0.0866),
            "p0866p1732" : (0.0866, 0.1732),
            "p1732p2598" : (0.1732, 0.2598),
        },
    },
    "2016" : {
        "DeepCSVMedium" : {
            "p00p1949" : (0.00, 0.1949),
            "p1949p3898" : (0.1949, 0.3898),
            "p3898p5847" : (0.3898, 0.5847),
        },
        "DeepJetMedium" : {
            "p00p0830" : (0.00, 0.0830),
            "p0830p1660" : (0.0830, 0.1660),
            "p1660p2489" : (0.1660, 0.2489),
        },
    },
    "2017" : {
        "DeepCSVMedium" : {
            "p00p1502" : (0.00, 0.1502),
            "p1502p3004" : (0.1502, 0.3004),
            "p3004p4506" : (0.3004, 0.4506),
        },
        "DeepJetMedium" : {
            "p00p1013" : (0.00, 0.1013),
            "p1013p2026" : (0.1013, 0.2026),
            "p2026p3040" : (0.2026, 0.3040),
        },
    },
    "2018" : {
        "DeepCSVMedium" : {
            "p00p1389" : (0.00, 0.1389),
            "p1389p2779" : (0.1389, 0.2779),
            "p2779p4168" : (0.2779, 0.4168),
        },
        "DeepJetMedium" : {
            "p00p0928" : (0.00, 0.0928),
            "p0928p1856" : (0.0928, 0.1856),
            "p1856p2783" : (0.1856, 0.2783),
        },
    },
}

btag_reg_names_dict = {
    "2016APV" : {
        "DeepCSVMedium" : {
            "Signal" : {"reg" : "btagPass"},
            "Down"   : {"reg" : "p00p20", "label" : "Down (0.0-0.2)", "color" : "b"},
            "Central": {"reg" : "p20p40", "label" : "Cen (0.2-0.4)", "color" : "k"},
            "Up"     : {"reg" : "p40p60", "label" : "Up (0.4-0.6)", "color" : "r"},
        },
        "DeepJetMedium" : {
            "Signal" : {"reg" : "btagPass"},
            "Down"   : {"reg" : "p00p0866", "label" : "Down (0.0-0.0866)", "color" : "b"},
            "Central": {"reg" : "p0866p1732", "label" : "Cen (0.0866-0.1732)", "color" : "k"},
            "Up"     : {"reg" : "p1732p2598", "label" : "Up (0.1732-0.2598)", "color" : "r"},
        },
    },
    "2016" : {
        "DeepCSVMedium" : {
            "Signal" : {"reg" : "btagPass"},
            "Down"   : {"reg" : "p00p1949", "label" : "Down (0.0-0.1949)", "color" : "b"},
            "Central": {"reg" : "p1949p3898", "label" : "Cen (0.1949-0.3898)", "color" : "k"},
            "Up"     : {"reg" : "p3898p5847", "label" : "Up (0.3898-0.5847)", "color" : "r"},
        },
        "DeepJetMedium" : {
            "Signal" : {"reg" : "btagPass"},
            "Down"   : {"reg" : "p00p0830", "label" : "Down (0.0-0.0830)", "color" : "b"},
            "Central": {"reg" : "p0830p1660", "label" : "Cen (0.0830-0.1660)", "color" : "k"},
            "Up"     : {"reg" : "p1660p2489", "label" : "Up (0.1660-0.2489)", "color" : "r"},
        },
    },
    "2017" : {
        "DeepCSVMedium" : {
            "Signal" : {"reg" : "btagPass"},
            "Down"   : {"reg" : "p00p1502", "label" : "Down (0.0-0.1502)", "color" : "b"},
            "Central": {"reg" : "p1502p3004", "label" : "Cen (0.1502-0.3004)", "color" : "k"},
            "Up"     : {"reg" : "p3004p4506", "label" : "Up (0.3004-0.4506)", "color" : "r"},
        },
        "DeepJetMedium" : {
            "Signal" : {"reg" : "btagPass"},
            "Down"   : {"reg" : "p00p1013", "label" : "Down (0.0-0.1013)", "color" : "b"},
            "Central": {"reg" : "p1013p2026", "label" : "Cen (0.1013-0.2026)", "color" : "k"},
            "Up"     : {"reg" : "p2026p3040", "label" : "Up (0.2026-0.3040)", "color" : "r"},
        },
    },
    "2018" : {
        "DeepCSVMedium" : {
            "Signal" : {"reg" : "btagPass"},
            "Down"   : {"reg" : "p00p1389", "label" : "Down (0.0-0.1389)", "color" : "b"},
            "Central": {"reg" : "p1389p2779", "label" : "Cen (0.1389-0.2779)", "color" : "k"},
            "Up"     : {"reg" : "p2779p4168", "label" : "Up (0.2779-0.4168)", "color" : "r"},
        },
        "DeepJetMedium" : {
            "Signal" : {"reg" : "btagPass"},
            "Down"   : {"reg" : "p00p0928", "label" : "Down (0.0-0.0928)", "color" : "b"},
            "Central": {"reg" : "p0928p1856", "label" : "Cen (0.0928-0.1856)", "color" : "k"},
            "Up"     : {"reg" : "p1856p2783", "label" : "Up (0.1856-0.2783)", "color" : "r"},
        },
    },
}

    # init dictionary of btag_groups and btag_cats
btag_cats = {year : {} for year in btag_sb_region_boundaries.keys()}
btag_groups = {year : {} for year in btag_sb_region_boundaries.keys()}

    # add btag sideband regions : region ranges from btag_sb_region_boundaries
    #   structure must be known beforehand
for year, year_dict in btag_sb_region_boundaries.items():
    for btagger, btag_dict in year_dict.items():
        btag_cats[year].update({key : f"{min} < max(btag discr) $\leq$ {max}" for key, (min, max) in btag_dict.items()})
    btag_cats[year]["btagPass"] ="$n_{btags} \geq$ 2"
    btag_groups[year].update({key : [key] for key in btag_cats[year].keys()})
