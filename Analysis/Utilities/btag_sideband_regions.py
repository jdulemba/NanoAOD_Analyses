btag_sb_region_boundaries = {
    "2016APV" : {
        "p00p20" : (0.00, 0.20),
        "p20p40" : (0.20, 0.40),
        "p40p60" : (0.40, 0.60),
    },
    "2016" : {
        "p00p1949" : (0.00, 0.1949),
        "p1949p3898" : (0.1949, 0.3898),
        "p3898p5847" : (0.3898, 0.5847),
    },
    "2017" : {
        "p00p1502" : (0.00, 0.1502),
        "p1502p3004" : (0.1502, 0.3004),
        "p3004p4506" : (0.3004, 0.4506),
    },
    "2018" : {
        "p00p1389" : (0.00, 0.1389),
        "p1389p2779" : (0.1389, 0.2779),
        "p2779p4168" : (0.2779, 0.4168),
    },
}

btag_reg_names_dict = {
    "2016APV" : {
        "Signal" : {"reg" : "btagPass"},
        "Down"   : {"reg" : "p00p20", "label" : "Down (0.0-0.2)", "color" : "b"},
        "Central": {"reg" : "p20p40", "label" : "Cen (0.2-0.4)", "color" : "k"},
        "Up"     : {"reg" : "p40p60", "label" : "Up (0.4-0.6)", "color" : "r"},
    },
    "2016" : {
        "Signal" : {"reg" : "btagPass"},
        "Down"   : {"reg" : "p00p1949", "label" : "Down (0.0-0.1949)", "color" : "b"},
        "Central": {"reg" : "p1949p3898", "label" : "Cen (0.1949-0.3898)", "color" : "k"},
        "Up"     : {"reg" : "p3898p5847", "label" : "Up (0.3898-0.5847)", "color" : "r"},
    },
    "2017" : {
        "Signal" : {"reg" : "btagPass"},
        "Down"   : {"reg" : "p00p1502", "label" : "Down (0.0-0.1502)", "color" : "b"},
        "Central": {"reg" : "p1502p3004", "label" : "Cen (0.1502-0.3004)", "color" : "k"},
        "Up"     : {"reg" : "p3004p4506", "label" : "Up (0.3004-0.4506)", "color" : "r"},
    },
    "2018" : {
        "Signal" : {"reg" : "btagPass"},
        "Down"   : {"reg" : "p00p1389", "label" : "Down (0.0-0.1389)", "color" : "b"},
        "Central": {"reg" : "p1389p2779", "label" : "Cen (0.1389-0.2779)", "color" : "k"},
        "Up"     : {"reg" : "p2779p4168", "label" : "Up (0.2779-0.4168)", "color" : "r"},
    },
}

btag_cats = {
    "btagPass" : "$n_{btags} \geq$ 2",
    "p00p20" : "0.0 < max(btag discr) $\leq$ 0.2",
    "p20p40" : "0.2 < max(btag discr) $\leq$ 0.4",
    "p40p60" : "0.4 < max(btag discr) $\leq$ 0.6",
    "p00p1949" : "0.0 < max(btag discr) $\leq$ 0.1949",
    "p1949p3898" : "0.1949 < max(btag discr) $\leq$ 0.3898",
    "p3898p5847" : "0.3898 < max(btag discr) $\leq$ 0.5847",
    "p00p1502" : "0.0 < max(btag discr) $\leq$ 0.1502",
    "p1502p3004" : "0.1502 < max(btag discr) $\leq$ 0.3004",
    "p3004p4506" : "0.3004 < max(btag discr) $\leq$ 0.4506",
    "p00p1389" : "0.0 < max(btag discr) $\leq$ 0.1389",
    "p1389p2779" : "0.1389 < max(btag discr) $\leq$ 0.2779",
    "p2779p4168" : "0.2779 < max(btag discr) $\leq$ 0.4168",
}

btag_groups = {
    "2016APV" : {
        "btagPass" : ["btagPass"],
        "p00p20" : ["p00p20"],
        "p20p40" : ["p20p40"],
        "p40p60" : ["p40p60"],
    },
    "2016" : {
        "btagPass" : ["btagPass"],
        "p00p1949" : ["p00p1949"],
        "p1949p3898" : ["p1949p3898"],
        "p3898p5847" : ["p3898p5847"],
    },
    "2017" : {
        "btagPass" : ["btagPass"],
        "p00p1502" : ["p00p1502"],
        "p1502p3004" : ["p1502p3004"],
        "p3004p4506" : ["p3004p4506"],
    },
    "2018" : {
        "btagPass" : ["btagPass"],
        "p00p1389" : ["p00p1389"],
        "p1389p2779" : ["p1389p2779"],
        "p2779p4168" : ["p2779p4168"],
    },
}
