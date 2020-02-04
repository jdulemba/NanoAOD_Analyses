import coffea.nanoaod.nanoevents
import coffea.processor.dataframe
from pdb import set_trace

def mu_triggers(df):
    if isinstance(df, coffea.nanoaod.nanoevents.NanoEvents):
        trigger = (df['HLT']['IsoMu24']) | (df['HLT']['IsoTkMu24'])
    elif isinstance(df, coffea.processor.dataframe.LazyDataFrame):
        trigger = (df['HLT_IsoMu24'] > 0) | (df['HLT_IsoTkMu24'] > 0)
    else:
        raise ValueError("Only NanoEvents and LazyDataFrame formats supported right now")

    return trigger

def el_triggers(df):
    if isinstance(df, coffea.nanoaod.nanoevents.NanoEvents):
        trigger = (df['HLT']['Ele27_WPTight_Gsf'])
    elif isinstance(df, coffea.processor.dataframe.LazyDataFrame):
        trigger = df['HLT_Ele27_WPTight_Gsf']
    else:
        raise ValueError("Only NanoEvents and LazyDataFrame formats supported right now")

    return trigger
