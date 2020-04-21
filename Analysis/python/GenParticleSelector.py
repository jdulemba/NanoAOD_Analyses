import numpy as np
from pdb import set_trace
import coffea.processor as processor
import awkward
import python.GenObjects as genobj
from coffea.analysis_objects import JaggedCandidateArray

def process_genParts(df):
    genParts = JaggedCandidateArray.candidatesfromcounts(
        df.nGenPart,
        pt=df.GenPart_pt,
        eta=df.GenPart_eta,
        phi=df.GenPart_phi,
        mass=df.GenPart_mass,
        momIdx=df.GenPart_genPartIdxMother,
        pdgId=df.GenPart_pdgId,
        status=df.GenPart_status,
        statusFlags=df.GenPart_statusFlags,
    )

    genParts['charge'] = genParts.pdgId.ones_like()*100.

    return genParts

def process_lheParts(df):
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

def select_lhe(df, w_decay_momid):
    if 'lheParts' not in df.columns:
        df['lheParts'] = process_lheParts(df)
    set_trace()

def select_normal(df, w_decay_momid):
    if 'genParts' not in df.columns:
        df['genParts'] = process_genParts(df)

    categories = processor.PackedSelection()
    categories.add('quarks', ( (df['genParts'].status > 21) & (df['genParts'].status < 30) & (df['genParts'].momIdx.counts > 0) ).flatten())
    categories.add('leptons', ( (df['genParts'].momIdx.counts > 0) & (np.abs(df['genParts'][df['genParts'].momIdx].pdgId) == w_decay_momid) ).flatten())
    categories.add('final_leptons', ( (df['genParts'].status == 1) & (df['genParts'].momIdx.counts > 0) & ( (np.abs(df['genParts'][df['genParts'].momIdx].pdgId) == 24) | (df['genParts'].pdgId == df['genParts'][df['genParts'].momIdx].pdgId) ) ).flatten())

    fermions = processor.PackedSelection()
    fermions.add('top',  categories.require(quarks=True) & (df['genParts'].pdgId == 6).flatten())
    fermions.add('tbar', categories.require(quarks=True) & (df['genParts'].pdgId == -6).flatten())
    fermions.add('b',    categories.require(quarks=True) & ( (df['genParts'].pdgId == 5) & (df['genParts'][df['genParts'].momIdx].pdgId != 24) ).flatten())
    fermions.add('bbar', categories.require(quarks=True) & ( (df['genParts'].pdgId == -5) & (df['genParts'][df['genParts'].momIdx].pdgId != -24) ).flatten())
    #fermions.add('wpartons', categories.require(quarks=True) & ( (np.abs(df['genParts'].pdgId) < 6) & (np.abs(df['genParts'][df['genParts'].momIdx].pdgId) == w_decay_momid) ).flatten())
    fermions.add('wpartons_up', categories.require(quarks=True) & ( (np.mod(df['genParts'].pdgId, 2) == 0) & (np.abs(df['genParts'].pdgId) < 6) & (np.abs(df['genParts'][df['genParts'].momIdx].pdgId) == w_decay_momid) ).flatten())
    fermions.add('wpartons_dw', categories.require(quarks=True) & ( (np.mod(df['genParts'].pdgId, 2) == 1) & (np.abs(df['genParts'].pdgId) < 6) & (np.abs(df['genParts'][df['genParts'].momIdx].pdgId) == w_decay_momid) ).flatten())

    fermions.add('charged_leps', categories.require(leptons=True) & ( (np.abs(df['genParts'].pdgId) == 11) | (np.abs(df['genParts'].pdgId) == 13) ).flatten())
    fermions.add('neutral_leps', categories.require(leptons=True) & ( (np.abs(df['genParts'].pdgId) == 12) | (np.abs(df['genParts'].pdgId) == 14) ).flatten())
    #fermions.add('tau_decay',    categories.require(leptons=True) & ( np.abs(df['genParts'].pdgId) == 16 ).flatten())
    fermions.add('final_charged_leps', categories.require(final_leptons=True) & ( (np.abs(df['genParts'].pdgId) == 11) | (np.abs(df['genParts'].pdgId) == 13) ).flatten())

    jagged_top_sel = awkward.JaggedArray.fromcounts(df['genParts'].counts, fermions.require(top=True))
    jagged_tbar_sel = awkward.JaggedArray.fromcounts(df['genParts'].counts, fermions.require(tbar=True))
    jagged_b_sel = awkward.JaggedArray.fromcounts(df['genParts'].counts, fermions.require(b=True))
    jagged_bbar_sel = awkward.JaggedArray.fromcounts(df['genParts'].counts, fermions.require(bbar=True))
    #jagged_wpartons_sel = awkward.JaggedArray.fromcounts(df['genParts'].counts, fermions.require(wpartons=True))
    jagged_wpartons_up_sel = awkward.JaggedArray.fromcounts(df['genParts'].counts, fermions.require(wpartons_up=True))
    jagged_wpartons_dw_sel = awkward.JaggedArray.fromcounts(df['genParts'].counts, fermions.require(wpartons_dw=True))

    jagged_chargedLeps_sel = awkward.JaggedArray.fromcounts(df['genParts'].counts, fermions.require(charged_leps=True))
    jagged_neutralLeps_sel = awkward.JaggedArray.fromcounts(df['genParts'].counts, fermions.require(neutral_leps=True))
    jagged_FinalChargedLeps_sel = awkward.JaggedArray.fromcounts(df['genParts'].counts, fermions.require(final_charged_leps=True))

    df['genParts']['charge'][jagged_top_sel] = df['genParts'][jagged_top_sel]['pdgId'].ones_like()*(2./3.)
    df['genParts']['charge'][jagged_tbar_sel] = df['genParts'][jagged_tbar_sel]['pdgId'].ones_like()*(-2./3.)
    df['genParts']['charge'][jagged_b_sel] = df['genParts'][jagged_b_sel]['pdgId'].ones_like()*(-1./3.)
    df['genParts']['charge'][jagged_bbar_sel] = df['genParts'][jagged_bbar_sel]['pdgId'].ones_like()*(1./3.)
    df['genParts']['charge'][jagged_wpartons_up_sel] = np.fmod(df['genParts'][jagged_wpartons_up_sel].pdgId+1, 2)*(2./3.)
    df['genParts']['charge'][jagged_wpartons_dw_sel] = np.fmod(df['genParts'][jagged_wpartons_dw_sel].pdgId, 2)*(-1./3.)
    df['genParts']['charge'][jagged_chargedLeps_sel] = np.fmod(df['genParts'][jagged_chargedLeps_sel].pdgId, 2)
    df['genParts']['charge'][jagged_neutralLeps_sel] = df['genParts'][jagged_neutralLeps_sel].pdgId.zeros_like()
    df['genParts']['charge'][jagged_FinalChargedLeps_sel] = np.fmod(df['genParts'][jagged_FinalChargedLeps_sel].pdgId, 2)

    tops = df['genParts'][jagged_top_sel]
    tbars = df['genParts'][jagged_tbar_sel]
    bs = df['genParts'][jagged_b_sel]
    bbars = df['genParts'][jagged_bbar_sel]
    wpartons_up = df['genParts'][jagged_wpartons_up_sel]
    wpartons_dw = df['genParts'][jagged_wpartons_dw_sel]
    #wpartons = df['genParts'][jagged_wpartons_sel]
    charged_leps = df['genParts'][jagged_chargedLeps_sel]
    neutral_leps = df['genParts'][jagged_neutralLeps_sel]
    final_charged_leps = df['genParts'][jagged_FinalChargedLeps_sel]

    #set_trace()
    GenObjects = awkward.Table(
        GenTop = tops,
        GenTbar = tbars,
        GenB = bs,
        GenBbar = bbars,
        #GenWParts = wpartons,
        GenWPartsUp = wpartons_up,
        GenWPartsDw = wpartons_dw,
        ChargedLeps = charged_leps,
        NeutralLeps = neutral_leps,
        FinalChargedLeps = final_charged_leps
    )

    return GenObjects

def select(df, mode='NORMAL'):
    modes_to_choose = {
        'NORMAL' : select_normal,
        'LHE' : select_lhe,
    }

    if mode not in modes_to_choose.keys():
        raise IOError("Gen Object mode %s not available" % mode)

    w_decay_momid = 6 if 'MADGRAPH' in mode else 24

    GenObjs = modes_to_choose[mode](df, w_decay_momid)

    set_trace()
    #ttbar = genobj.from_collections(wpartons=GenObjs['GenWParts'], charged_leps=GenObjs['ChargedLeps'], neutral_leps=GenObjs['NeutralLeps'], b=GenObjs['GenB'], bbar=GenObjs['GenBbar'], top=GenObjs['GenTop'], tbar=GenObjs['GenTbar'])
    ttbar = genobj.from_collections(wpartons_up=GenObjs['GenWPartsUp'], wpartons_dw=GenObjs['GenWPartsDw'], charged_leps=GenObjs['ChargedLeps'], neutral_leps=GenObjs['NeutralLeps'], b=GenObjs['GenB'], bbar=GenObjs['GenBbar'], top=GenObjs['GenTop'], tbar=GenObjs['GenTbar'])
    set_trace()
    #ttbar_final = genobj.from_collections(wpartons=GenObjs['GenWParts'], charged_leps=GenObjs['FinalChargedLeps'], neutral_leps=GenObjs['NeutralLeps'], b=GenObjs['GenB'], bbar=GenObjs['GenBbar'], top=GenObjs['GenTop'], tbar=GenObjs['GenTbar'])
    ttbar_final = genobj.from_collections(wpartons_up=GenObjs['GenWPartsUp'], wpartons_dw=GenObjs['GenWPartsDw'], charged_leps=GenObjs['FinalChargedLeps'], neutral_leps=GenObjs['NeutralLeps'], b=GenObjs['GenB'], bbar=GenObjs['GenBbar'], top=GenObjs['GenTop'], tbar=GenObjs['GenTbar'])



