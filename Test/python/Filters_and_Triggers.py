from coffea.util import save
import numpy as np

from pdb import set_trace

year = "2016"

single_el_trigger_paths = {}
single_mu_trigger_paths = {}

single_el_trigger_paths["2016"] = ["HLT_Ele27_WPTight_Gsf"]
single_mu_trigger_paths["2016"] = ["HLT_IsoMu24", "HLT_IsoTkMu24"]

def get_triggers(df, leptype, accumulator=None):
    ## event triggers to be used found here: https://twiki.cern.ch/twiki/bin/view/CMS/TopTriggerYear2016 or 2017, 2018...

    if leptype == 'Muon':
        triggers = [df[i] for i in single_mu_trigger_paths[year]]
        pass_triggers = np.stack(triggers, axis = 1).any(axis = 1)        

    elif leptype == 'Electron':
        triggers = [df[i] for i in single_el_trigger_paths[year]]
        pass_triggers = np.stack(triggers, axis = 1).any(axis = 1)        

    else:
        raise ValueError("Only events analyzing muons OR electrons supported right now")

    if accumulator:
        accumulator['cutflow']['nEvts pass %s pass_triggers' % leptype] += pass_triggers.sum()
        return pass_triggers, accumulator
    else:
        return pass_triggers


met_filters = {}

met_filters["2016"] = [
    "Flag_goodVertices",
    "Flag_globalSuperTightHalo2016Filter",
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter"
]


def get_filters(df, accumulator=None):
    ## Supported filters found here: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2

    filters = [df[i] for i in met_filters[year]]
    pass_filters = np.stack(filters, axis = 1).all(axis = 1)

    if accumulator:
        accumulator['cutflow']['pass filters'] += pass_filters.sum()
        return pass_filters, accumulator
    else:
        return pass_filters
