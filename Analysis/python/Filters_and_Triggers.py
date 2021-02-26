import awkward as ak
from pdb import set_trace

single_el_trigger_paths = {
    '2016' : {
        'Iso' : ["Ele27_WPTight_Gsf"],
        'noIso' : ["Ele27_WPLoose_Gsf"],
    },
    '2017' : {
        'Iso' : ["Ele27_WPTight_Gsf", "Ele32_WPTight_Gsf"],
        'noIso' : ["Ele20_WPLoose_Gsf"],
    },
    '2018' : {
        'Iso' : ["Ele32_WPTight_Gsf"],
        'noIso' : ["Ele20_WPLoose_Gsf"],
    },
}

single_mu_trigger_paths = {
    '2016' : {
        'Iso' : ["IsoMu24", "IsoTkMu24"],
        'noIso' : ["Mu50"],
    },
    '2017' : {
        'Iso' : ["IsoMu27"],
        'noIso' : ["Mu50"],
    },
    '2018' : {
        'Iso' : ["IsoMu24"],
        'noIso' : ["Mu50"],
    },
}

def get_triggers(HLT, leptype, year, noIso=False, accumulator=None):
    ## event triggers to be used found here: https://twiki.cern.ch/twiki/bin/view/CMS/TopTriggerYear2016 or 2017, 2018...

    allowed_leptypes = ['Muon', 'Electron']
    if leptype not in allowed_leptypes:
        raise ValueError(f"Select 'leptype' from {allowed_leptypes}")

    iso_cat = 'noIso' if noIso else 'Iso'
    pass_triggers = ak.any((HLT[i] for i in single_mu_trigger_paths[year][iso_cat] if i in HLT.fields), axis=0) if leptype == 'Muon'\
            else ak.any((HLT[i] for i in single_el_trigger_paths[year][iso_cat] if i in HLT.fields), axis=0)

    if accumulator:
        accumulator['cutflow']['nEvts pass %s pass_triggers' % leptype] += ak.sum(pass_triggers)
        return pass_triggers, accumulator
    else:
        return pass_triggers


## Supported filters found here: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
met_filters = {}
met_filters["2016"] = [
    "goodVertices",
    "globalSuperTightHalo2016Filter",
    "HBHENoiseFilter",
    "HBHENoiseIsoFilter",
    "EcalDeadCellTriggerPrimitiveFilter",
    "BadPFMuonFilter"
]
met_filters["2017"] = [
    "goodVertices",
    "globalSuperTightHalo2016Filter",
    "HBHENoiseFilter",
    "HBHENoiseIsoFilter",
    "EcalDeadCellTriggerPrimitiveFilter",
    "BadPFMuonFilter",
    #"ecalBadCalibFilter", ## should not be used now
]
met_filters["2018"] = [
    "goodVertices",
    "globalSuperTightHalo2016Filter",
    "HBHENoiseFilter",
    "HBHENoiseIsoFilter",
    "EcalDeadCellTriggerPrimitiveFilter",
    "BadPFMuonFilter",
    #"ecalBadCalibFilter", ## should not be used now
]

def get_filters(Flags, year, accumulator=None):

    pass_filters = ak.all((Flags[i] for i in met_filters[year]), axis=0)

    if accumulator:
        accumulator['cutflow']['pass filters'] += ak.sum(pass_filters)
        return pass_filters, accumulator
    else:
        return pass_filters

