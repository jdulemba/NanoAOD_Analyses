from coffea.analysis_objects import JaggedCandidateArray
from pdb import set_trace

def process_muons(dataframe):
    muons = JaggedCandidateArray.candidatesfromcounts(
        dataframe['nMuon'],
        pt=dataframe['Muon_pt'],
        eta=dataframe['Muon_eta'],
        phi=dataframe['Muon_phi'],
        mass=dataframe['Muon_mass'],
        charge=dataframe['Muon_charge'],
        #softId=dataframe['Muon_softId'],
        tightId=dataframe['Muon_tightId'],
        #IsoMutrig=dataframe['HLT_IsoMu24'],
        #IsoTkMutrig=dataframe['HLT_IsoTkMu24'],
    )
    return muons

