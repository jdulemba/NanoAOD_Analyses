from coffea.util import save
import numpy as np

from pdb import set_trace

single_el_trigger_paths = {
    '2016' : {
        'Iso' : ["HLT_Ele27_WPTight_Gsf"],
        'noIso' : ["HLT_Ele27_WPLoose_Gsf"],
    },
    '2017' : {
        'Iso' : ["HLT_Ele27_WPTight_Gsf", "HLT_Ele32_WPTight_Gsf"],
        'noIso' : ["HLT_Ele20_WPLoose_Gsf"],
    },
    '2018' : {
        'Iso' : ["HLT_Ele32_WPTight_Gsf"],
        'noIso' : ["HLT_Ele20_WPLoose_Gsf"],
    },
}

single_mu_trigger_paths = {
    '2016' : {
        'Iso' : ["HLT_IsoMu24", "HLT_IsoTkMu24"],
        'noIso' : ["HLT_Mu50"],
    },
    '2017' : {
        'Iso' : ["HLT_IsoMu27"],
        'noIso' : ["HLT_Mu50"],
    },
    '2018' : {
        'Iso' : ["HLT_IsoMu24"],
        'noIso' : ["HLT_Mu50"],
    },
}

def get_triggers(df, leptype, year, noIso=False, accumulator=None):
    ## event triggers to be used found here: https://twiki.cern.ch/twiki/bin/view/CMS/TopTriggerYear2016 or 2017, 2018...

    iso_cat = 'noIso' if noIso else 'Iso'
    if leptype == 'Muon':
        triggers = [df[i] for i in single_mu_trigger_paths[year][iso_cat] if i in df.columns]
        pass_triggers = np.stack(triggers, axis = 1).any(axis = 1)

    elif leptype == 'Electron':
        triggers = [df[i] for i in single_el_trigger_paths[year][iso_cat] if i in df.columns]
        pass_triggers = np.stack(triggers, axis = 1).any(axis = 1)

    else:
        raise ValueError("Only events analyzing muons OR electrons supported right now")

    if accumulator:
        accumulator['cutflow']['nEvts pass %s pass_triggers' % leptype] += pass_triggers.sum()
        return pass_triggers, accumulator
    else:
        return pass_triggers


## Supported filters found here: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
met_filters = {}
met_filters["2016"] = [
    "Flag_goodVertices",
    "Flag_globalSuperTightHalo2016Filter",
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter"
]
met_filters["2017"] = [
    "Flag_goodVertices",
    "Flag_globalSuperTightHalo2016Filter",
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter",
    #"Flag_ecalBadCalibFilter", ## should not be used now
]
met_filters["2018"] = [
    "Flag_goodVertices",
    "Flag_globalSuperTightHalo2016Filter",
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter",
    #"Flag_ecalBadCalibFilter", ## should not be used now
]

def get_filters(df, year, accumulator=None):

    filters = [df[i] for i in met_filters[year]]
    pass_filters = np.stack(filters, axis = 1).all(axis = 1)

    if accumulator:
        accumulator['cutflow']['pass filters'] += pass_filters.sum()
        return pass_filters, accumulator
    else:
        return pass_filters

