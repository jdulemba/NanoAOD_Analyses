import numpy as np
from pdb import set_trace

def process_genParts(df):
    from coffea.analysis_objects import JaggedCandidateArray
    genParts = JaggedCandidateArray.candidatesfromcounts(
        df.nGenPart,
        pt=df.GenPart_pt,
        eta=df.GenPart_eta,
        phi=df.GenPart_phi,
        mass=df.GenPart_mass,
        momIdx=df.GenPart_genPartIdxMother,
        pdgId=df.GenPart_pdgId
    )

    return genParts

def process_lheParts(df):
    from coffea.analysis_objects import JaggedCandidateArray
    lheParts = JaggedCandidateArray.candidatesfromcounts(
        df.nLHEPart,
        pt=df.LHEPart_pt,
        eta=df.LHEPart_eta,
        phi=df.LHEPart_phi,
        mass=df.LHEPart_mass,
        pdgId=df.LHEPart_pdgId
    )

    return lheParts

def process_genJets(df):
    from coffea.analysis_objects import JaggedCandidateArray

    genJets = JaggedCandidateArray.candidatesfromcounts(
        df.nGenJet,
        pt=df.GenJet_pt,
        eta=df.GenJet_eta,
        phi=df.GenJet_phi,
        mass=df.GenJet_mass,
        pFlav=df.GenJet_partonFlavour,
        hFlav=df.GenJet_hadronFlavour,
    )

    return genJets