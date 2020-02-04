from pdb import set_trace

def mu_triggers(dataframe):

    IsoMu24 = dataframe['HLT_IsoMu24']
    IsoTkMu24 = dataframe['HLT_IsoTkMu24']

    trigger = (IsoMu24 > 0) | (IsoTkMu24 > 0)
    
    return trigger

def el_triggers(dataframe):
    trigger = dataframe['HLT_Ele27_WPTight_Gsf']

    return trigger
