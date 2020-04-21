import numpy as np
from pdb import set_trace
import coffea.processor as processor
import awkward
from Utilities.PDGID import PDGID
from coffea.analysis_objects import JaggedCandidateArray

def GenW(first, second):#, DecayType='INVALID'):
    # invalid decays == 0, leptonic  == 1, hadronic == 2

    set_trace()
    #genw_sel.add('LEPTONIC', np.abs(first.pdgId) >= 
    #genw_sel.add('LEPTONIC', ( (np.abs(first.pdgId) >= PDGID['e']) & (np.abs(first.pdgId) <= PDGID['nu_tau']) ))
    leptonic_decay = (np.abs(first.pdgId) >= PDGID['e']) & (np.abs(first.pdgId) <= PDGID['nu_tau'])
    hadronic_decay = (np.abs(first.pdgId) < PDGID['e']) | (np.abs(first.pdgId) > PDGID['nu_tau'])
    decay_type = awkward.JaggedArray.fromcounts(first.counts, np.zeros(first.size, dtype=int))
    decay_type[leptonic_decay] = decay_type[leptonic_decay].ones_like()
    decay_type[hadronic_decay] = decay_type[hadronic_decay].ones_like()*2

    charge = awkward.JaggedArray.fromcounts(first.counts, np.zeros(first.size, dtype=int))
    charge[leptonic_decay] = (np.fmod(first.pdgId, 2)+np.fmod(second.pdgId, 2))[leptonic_decay]

    GenW = JaggedCandidateArray.candidatesfromcounts(
        counts=first.counts,
        pt=(first.p4+second.p4).pt,
        eta=(first.p4+second.p4).eta,
        phi=(first.p4+second.p4).phi,
        mass=(first.p4+second.p4).mass,
        Charge=charge,
        DecayType=decay_type,
        First=first,
        Second=second,
    )
        

#def from_collections(wpartons=None, charged_leps=None, neutral_leps=None, b=None, bbar=None, top=None, tbar=None):
def from_collections(wpartons_up=None, wpartons_dw=None, charged_leps=None, neutral_leps=None, b=None, bbar=None, top=None, tbar=None):

    #set_trace()
        # events with W decays
    min_lep_counts = charged_leps.counts if charged_leps.counts.sum() < neutral_leps.counts.sum() else neutral_leps.counts # get array with counts for leptonic W events (charged and neutral leps don't have same # of events for some reason...

    lep_sel = ((charged_leps.counts > 0) & (charged_leps.counts == neutral_leps.counts))
    had_sel = ((wpartons_up.counts > 0) & (wpartons_dw.counts > 0) & (wpartons_up.counts == wpartons_dw.counts))
    dilep = (lep_sel & ~had_sel)
    dilep_array = np.repeat(dilep, (min_lep_counts+wpartons_up.counts))
    dihad = (had_sel & ~lep_sel)
    dihad_array = np.repeat(dihad, (min_lep_counts+wpartons_up.counts))
    semilep = (had_sel & lep_sel)
    sl_lep = np.repeat(semilep, (min_lep_counts+wpartons_up.counts))
    sl_had = np.repeat(semilep, (min_lep_counts+wpartons_up.counts))
    sl_lep_true = sl_lep[sl_lep == True]
    sl_lep_true[0::2] = False
    sl_lep[sl_lep == True] = sl_lep_true
    sl_had_true = sl_had[sl_had == True]
    sl_had_true[1::2] = False
    sl_had[sl_had == True] = sl_had_true
    leptonic_evts = awkward.JaggedArray.fromcounts( (min_lep_counts+wpartons_up.counts), (sl_lep | dilep_array))
    hadronic_evts = awkward.JaggedArray.fromcounts( (min_lep_counts+wpartons_up.counts), (sl_had | dihad_array))
    w_evts = awkward.JaggedArray.fromcounts( (min_lep_counts+wpartons_up.counts), (leptonic_evts.flatten() | hadronic_evts.flatten()) )

        # initialize variables to become attributes of GenW objecs
    pt = np.ones(w_evts.counts.sum())
    eta = np.ones(w_evts.counts.sum())
    phi = np.ones(w_evts.counts.sum())
    mass = np.ones(w_evts.counts.sum())
    Charge=np.zeros(w_evts.counts.sum())
    DecayType=np.zeros(w_evts.counts.sum())
    First=np.repeat( np.array([None]), w_evts.counts.sum())
    Second=np.repeat( np.array([None]), w_evts.counts.sum())
    Up=np.repeat( np.array([None]), w_evts.counts.sum())
    Down=np.repeat( np.array([None]), w_evts.counts.sum())
    #set_trace()

        # set values for leptonic Ws
    valid_charged_leps = charged_leps[leptonic_evts.sum() > 0]
    valid_neutral_leps = neutral_leps[leptonic_evts.sum() > 0]
    lep_id = valid_charged_leps.pdgId
    nu_id = -1*(lep_id + np.abs(lep_id)/lep_id)
    if (np.abs(valid_charged_leps.charge + valid_neutral_leps.charge) != 1).any().any():
        raise ValueError("Charged and neutral leptons don't have same corresponding indices!")
    lepW_p4 = (valid_charged_leps.p4 + valid_neutral_leps.p4)
    pt[leptonic_evts.flatten()] = lepW_p4.pt.flatten()
    eta[leptonic_evts.flatten()] = lepW_p4.eta.flatten()
    phi[leptonic_evts.flatten()] = lepW_p4.phi.flatten()
    mass[leptonic_evts.flatten()] = lepW_p4.mass.flatten()
    Charge[leptonic_evts.flatten()] = (valid_charged_leps.charge + valid_neutral_leps.charge).flatten()
    DecayType[leptonic_evts.flatten()] = np.ones(leptonic_evts.flatten().sum())
    First[leptonic_evts.flatten()] = valid_charged_leps.flatten()
    Second[leptonic_evts.flatten()] = valid_neutral_leps.flatten()
    Up[leptonic_evts.flatten()] = valid_neutral_leps.flatten()
    Down[leptonic_evts.flatten()] = valid_charged_leps.flatten()
    #set_trace()

        # set values for hadronic Ws
    if (np.abs(wpartons_up[hadronic_evts.sum() > 0].charge + wpartons_dw[hadronic_evts.sum() > 0].charge) != 1).any().any():
        raise ValueError("Up and Down-type wpartons don't have same corresponding indices!")
    hadW_p4 = (wpartons_up[hadronic_evts.sum() > 0].p4 + wpartons_dw[hadronic_evts.sum() > 0].p4)
    pt[hadronic_evts.flatten()] = hadW_p4.pt.flatten()
    eta[hadronic_evts.flatten()] = hadW_p4.eta.flatten()
    phi[hadronic_evts.flatten()] = hadW_p4.phi.flatten()
    mass[hadronic_evts.flatten()] = hadW_p4.mass.flatten()
    Charge[hadronic_evts.flatten()] = (wpartons_up[hadronic_evts.sum() > 0].charge + wpartons_dw[hadronic_evts.sum() > 0].charge).flatten()
    DecayType[hadronic_evts.flatten()] = np.ones(hadronic_evts.flatten().sum())*2
    up_isLeading = (wpartons_up[hadronic_evts.sum() > 0].pt > wpartons_dw[hadronic_evts.sum() > 0].pt).flatten()
    had_first = First[hadronic_evts.flatten()]
    had_first[up_isLeading] = wpartons_up[hadronic_evts.sum() > 0].flatten()[up_isLeading]
    had_first[~up_isLeading] = wpartons_dw[hadronic_evts.sum() > 0].flatten()[~up_isLeading]
    First[hadronic_evts.flatten()] = had_first
    had_second = Second[hadronic_evts.flatten()]
    had_second[up_isLeading] = wpartons_dw[hadronic_evts.sum() > 0].flatten()[up_isLeading]
    had_second[~up_isLeading] = wpartons_up[hadronic_evts.sum() > 0].flatten()[~up_isLeading]
    Second[hadronic_evts.flatten()] = had_second
    Up[hadronic_evts.flatten()] = wpartons_up[hadronic_evts.sum() > 0].flatten()
    Down[hadronic_evts.flatten()] = wpartons_dw[hadronic_evts.sum() > 0].flatten()

    set_trace()

    GenW = JaggedCandidateArray.candidatesfromcounts(
        counts=w_evts.counts,
        pt=pt,
        eta=eta,
        phi=phi,
        mass=mass,
        Charge=Charge,
        DecayType=DecayType,
        #First=First,
        #Second=Second,
        #Up=Up,
        #Down=Down,
    )



