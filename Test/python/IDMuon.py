#from coffea.analysis_objects import JaggedCandidateArray
import numpy as np
from pdb import set_trace


#def process_muons(dataframe):
#    muons = JaggedCandidateArray.candidatesfromcounts(
#        dataframe['nMuon'],
#        pt=dataframe['Muon_pt'],
#        eta=dataframe['Muon_eta'],
#        phi=dataframe['Muon_phi'],
#        mass=dataframe['Muon_mass'],
#        charge=dataframe['Muon_charge'],
#        #softId=dataframe['Muon_softId'],
#        tightId=dataframe['Muon_tightId'],
#        #IsoMutrig=dataframe['HLT_IsoMu24'],
#        #IsoTkMutrig=dataframe['HLT_IsoTkMu24'],
#    )
#    return muons
#


    ## Muon ID types
def fail(df):
    return df

def loose_12(df):
    eta_cut = (df['Muon'].eta <= 2.4)
    ID = (df['Muon'].looseId)
    Iso = (df['Muon'].tkRelIso < 0.10)

    return (eta_cut & ID & Iso)

def tight_12(df):
    eta_cut = (df['Muon'].eta <= 2.4)
    ID = (df['Muon'].tightId)
    Iso = (df['Muon'].tkRelIso < 0.05)

    return (eta_cut & ID & Iso)

def loose_12Db(df):
    eta_cut = (df['Muon'].eta <= 2.4)
    ID = (df['Muon'].looseId)
    Iso = (df['Muon'].pfRelIso04_all < 0.12)

    return (eta_cut & ID & Iso)

def tight_12Db(df):
    eta_cut = (df['Muon'].eta <= 2.4)
    ID = (df['Muon'].tightId)
    Iso = (df['Muon'].pfRelIso04_all < 0.12)

    return (eta_cut & ID & Iso)

def loose_15(df):
    ID = (df['Muon'].looseId)
    Iso = Iso = (df['Muon'].tkRelIso < 0.10)

    return (ID & Iso)

def tight_15(df):
    ID = (df['Muon'].tightId)
    Iso = Iso = (df['Muon'].tkRelIso < 0.05)

    return (ID & Iso)

def loose_15Db(df):
    ID = (df['Muon'].looseId)
    Iso = Iso = (df['Muon'].pfRelIso04_all < 0.25)

    return (ID & Iso)

def tight_15Db(df):
    ID = (df['Muon'].tightId)
    Iso = Iso = (df['Muon'].pfRelIso04_all < 0.15)

    return (ID & Iso)

def tight_noIso(df):
    return (df['Muon'].tightId)

def antiloose_15Db(df):
    ID = (df['Muon'].tightId)
    Iso = (df['Muon'].pfRelIso04_all >= 0.15) & (df['Muon'].pfRelIso04_all < 0.43)

    return (ID & Iso)

def make_muon_ids(df):

    mu_types = {
        'VETOMU' : {
            'id' : 'LOOSE_15Db',
            'ptmin' : 10.,
            'etamax' : 2.4
        },
        'LOOSEMU' : {
            'id' : 'ANTILOOSE_15Db',
            'ptmin' : 26.,
            'etamax' : 2.4
        },
        'TIGHTMU' : {
            'id' : 'TIGHT_15Db',
            'ptmin' : 26.,
            'etamax' : 2.4
        }
    }

    id_names = {
        'FAIL' : fail(df),
        'LOOSE_12' : loose_12(df),
        'TIGHT_12' : tight_12(df),
        'LOOSE_12Db' : loose_12Db(df),
        'TIGHT_12Db' : tight_12Db(df),
        'LOOSE_15' : loose_15(df),
        'TIGHT_15' : tight_15(df),
        'LOOSE_15Db' : loose_15Db(df),
        'TIGHT_15Db' : tight_15Db(df),
        'TIGHT_NOISO' : tight_noIso(df),
        'ANTILOOSE_15Db' : antiloose_15Db(df)
    }

    if mu_types['VETOMU']['id'] not in id_names.keys():
        raise IOError("veto Muon ID name not valid")
    if mu_types['LOOSEMU']['id'] not in id_names.keys():
        raise IOError("loose Muon ID name not valid")
    if mu_types['TIGHTMU']['id'] not in id_names.keys():
        raise IOError("tight Muon ID name not valid")

    #set_trace()
    for muID in mu_types.keys():
        pt_cut = (df['Muon'].pt >= mu_types[muID]['ptmin'])
        eta_cut = (np.abs(df['Muon'].eta) <= mu_types[muID]['etamax'])
        pass_id = id_names[mu_types[muID]['id']]
        df['Muon'][muID] = (pass_id) & (pt_cut) & (eta_cut)

    return df
