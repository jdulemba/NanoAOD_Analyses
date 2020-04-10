import numpy as np
from pdb import set_trace
import Utilities.prettyjson as prettyjson
import os

def process_muons(df, year):
    from coffea.analysis_objects import JaggedCandidateArray
    muons = JaggedCandidateArray.candidatesfromcounts(
        df['nMuon'],
        pt=df['Muon_pt'],
        eta=df['Muon_eta'],
        phi=df['Muon_phi'],
        mass=df['Muon_mass'],
        charge=df['Muon_charge'],
        looseId=df['Muon_looseId'],
        tightId=df['Muon_tightId'],
        pfRelIso=df['Muon_pfRelIso04_all'],
        tkRelIso=df['Muon_tkRelIso'],
    )

        ## add muon ID+Iso categories
    muons = make_muon_ids(muons, year)

    return muons



    ## Muon ID types
def fail(muons):
    return muons

def loose_12(muons):
    ID  = (muons.looseId)
    Iso = (muons.tkRelIso < 0.10)

    return (ID & Iso)

def tight_12(muons):
    ID  = (muons.tightId)
    Iso = (muons.tkRelIso < 0.05)

    return (ID & Iso)

def loose_12Db(muons):
    ID  = (muons.looseId)
    Iso = (muons.pfRelIso < 0.12)

    return (ID & Iso)

def tight_12Db(muons):
    ID  = (muons.tightId)
    Iso = (muons.pfRelIso < 0.12)

    return (ID & Iso)

def loose_15(muons):
    ID  = (muons.looseId)
    Iso = (muons.tkRelIso < 0.10)

    return (ID & Iso)

def tight_15(muons):
    ID  = (muons.tightId)
    Iso = (muons.tkRelIso < 0.05)

    return (ID & Iso)

def loose_15Db(muons):
    ID  = (muons.looseId)
    Iso = (muons.pfRelIso < 0.25)

    return (ID & Iso)

def tight_15Db(muons):
    ID  = (muons.tightId)
    Iso = (muons.pfRelIso < 0.15)

    return (ID & Iso)

def tight_noIso(muons):
    return (muons.tightId)

def antiloose_15Db(muons):
    ID  = (muons.tightId)
    Iso = (muons.pfRelIso >= 0.15) & (muons.pfRelIso < 0.43)

    return (ID & Iso)

def make_muon_ids(muons, year):

    mu_pars = prettyjson.loads(open('%s/cfg_files/cfg_pars_%s.json' % (os.environ['PROJECT_DIR'], os.environ['jobid'])).read())['Muons'][year]

    id_names = {
        'FAIL' : fail,
        'LOOSE_12' : loose_12,
        'TIGHT_12' : tight_12,
        'LOOSE_12Db' : loose_12Db,
        'TIGHT_12Db' : tight_12Db,
        'LOOSE_15' : loose_15,
        'TIGHT_15' : tight_15,
        'LOOSE_15Db' : loose_15Db,
        'TIGHT_15Db' : tight_15Db,
        'TIGHT_NOISO' : tight_noIso,
        'ANTILOOSE_15Db' : antiloose_15Db
    }

    if mu_pars['VETOMU']['id'] not in id_names.keys():
        raise IOError("veto Muon ID name not valid")
    if mu_pars['LOOSEMU']['id'] not in id_names.keys():
        raise IOError("loose Muon ID name not valid")
    if mu_pars['TIGHTMU']['id'] not in id_names.keys():
        raise IOError("tight Muon ID name not valid")

    for muID in mu_pars.keys():
        pt_cut = (muons.pt >= mu_pars[muID]['ptmin'])
        eta_cut = (np.abs(muons.eta) <= mu_pars[muID]['etamax'])
        pass_id = id_names[mu_pars[muID]['id']](muons)
        muons[muID] = (pass_id) & (pt_cut) & (eta_cut)

    return muons
