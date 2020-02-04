from pdb import set_trace
from python.IDMuon import process_muons as proc_mus
from python.IDElectron import process_electrons as proc_els
from python.IDJet import process_jets as proc_jets
import python.triggers as triggers

def select_muons(df):
    muons = proc_mus(df) # get muons

        # pass triggers
    trig_muons = triggers.mu_triggers(df)
    muons = muons[trig_muons]

        # single muon    
    onemuon = (muons.counts == 1)
    muons = muons[onemuon]

        # tight muon    
    tight_mu = (muons.tightId > 0)
    muons = muons[tight_mu]

    return muons


def select_jets(df):
    jets = proc_jets(df) # get jets

    return jets
