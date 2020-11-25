import numpy as np
from pdb import set_trace
import awkward
from coffea.analysis_objects import JaggedCandidateArray

def from_collections(wpartons_up=None, wpartons_dw=None, charged_leps=None, neutral_leps=None, bs=None, bbars=None, tops=None, tbars=None):

    #set_trace()
    dilep_evts = ( (charged_leps.counts == 2) & (charged_leps.counts == neutral_leps.counts) & (charged_leps.charge.sum() == 0) )
    dihad_evts = ( (wpartons_up.counts == 2) & (wpartons_up.counts == wpartons_dw.counts) )
    semilep_evts = ( (charged_leps.counts == 1) & (neutral_leps.counts == 1) & (wpartons_up.counts == 1) & (wpartons_dw.counts == 1) )
    valid_evts = (dilep_evts | dihad_evts | semilep_evts)
    if not (valid_evts.sum() == tbars.size):
        raise ValueError("Not all ttbar events are valid!")
    dl_array = np.repeat(dilep_evts, valid_evts.astype(int)*2)
    dh_array = np.repeat(dihad_evts, valid_evts.astype(int)*2)
    sl_array = np.repeat(semilep_evts, valid_evts.astype(int)*2)

    if (bs[valid_evts].counts != 1).any():
        raise ValueError("Number of b partons in valid events != 1")
    if (bbars[valid_evts].counts != 1).any():
        raise ValueError("Number of bbar partons in valid events != 1")
    if (tops[valid_evts].counts != 1).any():
        raise ValueError("Number of top partons in valid events != 1")
    if (tbars[valid_evts].counts != 1).any():
        raise ValueError("Number of tbar partons in valid events != 1")

        # initialize variables to become attributes of GenW objecs
    Wpt       = np.ones(valid_evts.sum()*2)
    Weta      = np.ones(valid_evts.sum()*2)
    Wphi      = np.ones(valid_evts.sum()*2)
    Wmass     = np.ones(valid_evts.sum()*2)
    Wcharge   =np.zeros(valid_evts.sum()*2)
    Wdecaytype=np.zeros(valid_evts.sum()*2)

        # set values for fully leptonic evts
    dlW_p4 = (charged_leps[dilep_evts].p4 + neutral_leps[dilep_evts].p4)
    Wpt[dl_array]       = dlW_p4.pt.flatten()
    Weta[dl_array]      = dlW_p4.eta.flatten()
    Wphi[dl_array]      = dlW_p4.phi.flatten()
    Wmass[dl_array]     = dlW_p4.mass.flatten()
    Wcharge[dl_array]   = (charged_leps[dilep_evts].charge + neutral_leps[dilep_evts].charge).flatten()
    Wdecaytype[dl_array]= dlW_p4.pt.ones_like().flatten() # 0 is INVALID (shouldn't happen), 1 is LEPTONIC, 2 is HADRONIC
 
        # set values for fully hadronic evts
    dhW_p4 = (wpartons_up[dihad_evts].p4 + wpartons_dw[dihad_evts].p4)
    Wpt[dh_array]       = dhW_p4.pt.flatten()
    Weta[dh_array]      = dhW_p4.eta.flatten()
    Wphi[dh_array]      = dhW_p4.phi.flatten()
    Wmass[dh_array]     = dhW_p4.mass.flatten()
    Wcharge[dh_array]   = (wpartons_up[dihad_evts].charge + wpartons_dw[dihad_evts].charge).flatten()
    Wdecaytype[dh_array]= (dhW_p4.pt.ones_like()*2).flatten() # 0 is INVALID (shouldn't happen), 1 is LEPTONIC, 2 is HADRONIC
 
        # set values for semileptonic evts
    sl_inds = np.where(sl_array)[0]
    sl_had_inds = sl_inds[np.mod(sl_inds, 2) == 0] # hadronic decays are defined as first W
    sl_lep_inds = sl_inds[np.mod(sl_inds, 2) == 1] # leptonic decays are defined as second W
            ## leptonic decays
    slWlep_p4 = (charged_leps[semilep_evts].p4 + neutral_leps[semilep_evts].p4)
    Wpt[sl_lep_inds]       = slWlep_p4.pt.flatten()
    Weta[sl_lep_inds]      = slWlep_p4.eta.flatten()
    Wphi[sl_lep_inds]      = slWlep_p4.phi.flatten()
    Wmass[sl_lep_inds]     = slWlep_p4.mass.flatten()
    Wcharge[sl_lep_inds]   = (charged_leps[semilep_evts].charge + neutral_leps[semilep_evts].charge).flatten()
    Wdecaytype[sl_lep_inds]= (slWlep_p4.pt.ones_like()).flatten() # 0 is INVALID (shouldn't happen), 1 is LEPTONIC, 2 is HADRONIC
            ## hadronic decays
    slWhad_p4 = (wpartons_up[semilep_evts].p4 + wpartons_dw[semilep_evts].p4)
    Wpt[sl_had_inds]       = slWhad_p4.pt.flatten()
    Weta[sl_had_inds]      = slWhad_p4.eta.flatten()
    Wphi[sl_had_inds]      = slWhad_p4.phi.flatten()
    Wmass[sl_had_inds]     = slWhad_p4.mass.flatten()
    Wcharge[sl_had_inds]   = (wpartons_up[semilep_evts].charge + wpartons_dw[semilep_evts].charge).flatten()
    Wdecaytype[sl_had_inds]= (slWhad_p4.pt.ones_like()*2).flatten() # 0 is INVALID (shouldn't happen), 1 is LEPTONIC, 2 is HADRONIC

        # create W
    GenW = JaggedCandidateArray.candidatesfromcounts(
        counts=valid_evts.astype(int)*2,
        pt       = Wpt,
        eta      = Weta,
        phi      = Wphi,
        mass     = Wmass,
        charge   = Wcharge,
        decaytype= Wdecaytype,
    )

        # create leading/subleading objects from wpartons
            # init vars
    Lpt       = np.ones(wpartons_up.counts.sum())
    Leta      = np.ones(wpartons_up.counts.sum())
    Lphi      = np.ones(wpartons_up.counts.sum())
    Lmass     = np.ones(wpartons_up.counts.sum())
    Lcharge   = np.zeros(wpartons_up.counts.sum())
    LpdgId    = np.zeros(wpartons_up.counts.sum())
    Ldecaytype= np.ones(wpartons_up.counts.sum())*2
    Spt       = np.ones(wpartons_up.counts.sum())
    Seta      = np.ones(wpartons_up.counts.sum())
    Sphi      = np.ones(wpartons_up.counts.sum())
    Smass     = np.ones(wpartons_up.counts.sum())
    Scharge   = np.zeros(wpartons_up.counts.sum())
    SpdgId    = np.zeros(wpartons_up.counts.sum())

    up_isLeading = (wpartons_up.pt > wpartons_dw.pt).flatten()
            # Leading
    Lpt    = np.where(up_isLeading, wpartons_up.pt.flatten(), wpartons_dw.pt.flatten())
    Leta   = np.where(up_isLeading, wpartons_up.eta.flatten(), wpartons_dw.eta.flatten())
    Lphi   = np.where(up_isLeading, wpartons_up.phi.flatten(), wpartons_dw.phi.flatten())
    Lmass  = np.where(up_isLeading, wpartons_up.mass.flatten(), wpartons_dw.mass.flatten())
    Lcharge= np.where(up_isLeading, wpartons_up.charge.flatten(), wpartons_dw.charge.flatten())
    LpdgId = np.where(up_isLeading, wpartons_up.pdgId.flatten(), wpartons_dw.pdgId.flatten())
    First = JaggedCandidateArray.candidatesfromcounts(
        counts = wpartons_up.counts,
        pt       = Lpt,
        eta      = Leta,
        phi      = Lphi,
        mass     = Lmass,
        charge   = Lcharge,
        pdgId    = LpdgId,
        decaytype= Ldecaytype,
    )
            # Subleading
    Spt    = np.where(~up_isLeading, wpartons_up.pt.flatten(), wpartons_dw.pt.flatten())
    Seta   = np.where(~up_isLeading, wpartons_up.eta.flatten(), wpartons_dw.eta.flatten())
    Sphi   = np.where(~up_isLeading, wpartons_up.phi.flatten(), wpartons_dw.phi.flatten())
    Smass  = np.where(~up_isLeading, wpartons_up.mass.flatten(), wpartons_dw.mass.flatten())
    Scharge= np.where(~up_isLeading, wpartons_up.charge.flatten(), wpartons_dw.charge.flatten())
    SpdgId = np.where(~up_isLeading, wpartons_up.pdgId.flatten(), wpartons_dw.pdgId.flatten())
    Second = JaggedCandidateArray.candidatesfromcounts(
        counts = wpartons_up.counts,
        pt       = Spt,
        eta      = Seta,
        phi      = Sphi,
        mass     = Smass,
        charge   = Scharge,
        pdgId    = SpdgId,
        decaytype= Ldecaytype,
    )

        # make objects from top/tbar decays
    Wplus = GenW[GenW.charge == 1]
    GenB = JaggedCandidateArray.candidatesfromcounts(
        counts = valid_evts.astype(int),
        pt       = bs[valid_evts].pt.flatten(),
        eta      = bs[valid_evts].eta.flatten(),
        phi      = bs[valid_evts].phi.flatten(),
        mass     = bs[valid_evts].mass.flatten(),
        charge   = bs[valid_evts].charge.flatten(),
        decaytype= GenW[GenW.charge == 1].decaytype.flatten(),
    )
    GenTop = JaggedCandidateArray.candidatesfromcounts(
        counts = valid_evts.astype(int),
        pt       = tops[valid_evts].pt.flatten(),
        eta      = tops[valid_evts].eta.flatten(),
        phi      = tops[valid_evts].phi.flatten(),
        mass     = tops[valid_evts].mass.flatten(),
        charge   = tops[valid_evts].charge.flatten(),
        decaytype= GenW[GenW.charge == 1].decaytype.flatten(),
    )
    Wminus = GenW[GenW.charge == -1]
    GenBbar = JaggedCandidateArray.candidatesfromcounts(
        counts = valid_evts.astype(int),
        pt       = bbars[valid_evts].pt.flatten(),
        eta      = bbars[valid_evts].eta.flatten(),
        phi      = bbars[valid_evts].phi.flatten(),
        mass     = bbars[valid_evts].mass.flatten(),
        charge   = bbars[valid_evts].charge.flatten(),
        decaytype= GenW[GenW.charge == -1].decaytype.flatten(),
    )
    GenTbar = JaggedCandidateArray.candidatesfromcounts(
        counts = valid_evts.astype(int),
        pt       = tbars[valid_evts].pt.flatten(),
        eta      = tbars[valid_evts].eta.flatten(),
        phi      = tbars[valid_evts].phi.flatten(),
        mass     = tbars[valid_evts].mass.flatten(),
        charge   = tbars[valid_evts].charge.flatten(),
        decaytype= GenW[GenW.charge == -1].decaytype.flatten(),
    )

    ttbars = (GenTop.p4 + GenTbar.p4)
    GenTTbar = JaggedCandidateArray.candidatesfromcounts(
        counts = valid_evts.astype(int),
        pt       = ttbars.pt.flatten(),
        eta      = ttbars.eta.flatten(),
        phi      = ttbars.phi.flatten(),
        mass     = ttbars.mass.flatten(),
        decaytype= GenW[valid_evts].decaytype.sum(), # 0 is for INVALID, 2 for DILEP, 3 for SEMILEP, 4 for HADRONIC
    )    


    #set_trace()
        ## make tables of TTbar event objects
    DILEP_evts = awkward.Table(
        TTbar = awkward.JaggedArray.fromcounts(dilep_evts.astype(int), GenTTbar[dilep_evts].flatten()),
        Top   = awkward.JaggedArray.fromcounts(dilep_evts.astype(int), GenTop[dilep_evts].flatten()),
        Tbar  = awkward.JaggedArray.fromcounts(dilep_evts.astype(int), GenTbar[dilep_evts].flatten()),
        B     = awkward.JaggedArray.fromcounts(dilep_evts.astype(int), GenB[dilep_evts].flatten()),
        Bbar  = awkward.JaggedArray.fromcounts(dilep_evts.astype(int), GenBbar[dilep_evts].flatten()),
        Wplus = awkward.JaggedArray.fromcounts(dilep_evts.astype(int), GenW[GenW.charge == 1][dilep_evts].flatten()),
        Wminus = awkward.JaggedArray.fromcounts(dilep_evts.astype(int), GenW[GenW.charge == -1][dilep_evts].flatten()),
        First_plus = awkward.JaggedArray.fromcounts(dilep_evts.astype(int), charged_leps[dilep_evts][charged_leps[dilep_evts].charge > 0].flatten()), # charged lepton always made leading
        Second_plus = awkward.JaggedArray.fromcounts(dilep_evts.astype(int), neutral_leps[dilep_evts][charged_leps[dilep_evts].charge > 0].flatten()), # neutral lepton always made subleading
        First_minus = awkward.JaggedArray.fromcounts(dilep_evts.astype(int), charged_leps[dilep_evts][charged_leps[dilep_evts].charge < 0].flatten()), # charged lepton always made leading
        Second_minus = awkward.JaggedArray.fromcounts(dilep_evts.astype(int), neutral_leps[dilep_evts][charged_leps[dilep_evts].charge < 0].flatten()), # neutral lepton always made subleading
        Up_plus = awkward.JaggedArray.fromcounts(dilep_evts.astype(int), neutral_leps[dilep_evts][charged_leps[dilep_evts].charge > 0].flatten()),
        Down_plus = awkward.JaggedArray.fromcounts(dilep_evts.astype(int), charged_leps[dilep_evts][charged_leps[dilep_evts].charge > 0].flatten()),
        Up_minus = awkward.JaggedArray.fromcounts(dilep_evts.astype(int), neutral_leps[dilep_evts][charged_leps[dilep_evts].charge < 0].flatten()),
        Down_minus = awkward.JaggedArray.fromcounts(dilep_evts.astype(int), charged_leps[dilep_evts][charged_leps[dilep_evts].charge < 0].flatten()),
    )

    HAD_evts = awkward.Table(
        TTbar = awkward.JaggedArray.fromcounts(dihad_evts.astype(int), GenTTbar[dihad_evts].flatten()),
        Top   = awkward.JaggedArray.fromcounts(dihad_evts.astype(int), GenTop[dihad_evts].flatten()),
        Tbar  = awkward.JaggedArray.fromcounts(dihad_evts.astype(int), GenTbar[dihad_evts].flatten()),
        B     = awkward.JaggedArray.fromcounts(dihad_evts.astype(int), GenB[dihad_evts].flatten()),
        Bbar  = awkward.JaggedArray.fromcounts(dihad_evts.astype(int), GenBbar[dihad_evts].flatten()),
        Wplus = awkward.JaggedArray.fromcounts(dihad_evts.astype(int), GenW[GenW.charge == 1][dihad_evts].flatten()),
        Wminus= awkward.JaggedArray.fromcounts(dihad_evts.astype(int), GenW[GenW.charge == -1][dihad_evts].flatten()),
        First_plus   = awkward.JaggedArray.fromcounts(dihad_evts.astype(int), First[First.charge > 0][dihad_evts].flatten()),
        Second_plus  = awkward.JaggedArray.fromcounts(dihad_evts.astype(int), Second[Second.charge > 0][dihad_evts].flatten()),
        First_minus  = awkward.JaggedArray.fromcounts(dihad_evts.astype(int), First[First.charge < 0][dihad_evts].flatten()),
        Second_minus = awkward.JaggedArray.fromcounts(dihad_evts.astype(int), Second[Second.charge < 0][dihad_evts].flatten()),
        Up_plus      = awkward.JaggedArray.fromcounts(dihad_evts.astype(int), wpartons_up[wpartons_up.charge > 0][dihad_evts].flatten()),
        Down_plus    = awkward.JaggedArray.fromcounts(dihad_evts.astype(int), wpartons_dw[wpartons_dw.charge > 0][dihad_evts].flatten()),
        Up_minus     = awkward.JaggedArray.fromcounts(dihad_evts.astype(int), wpartons_up[wpartons_up.charge < 0][dihad_evts].flatten()),
        Down_minus   = awkward.JaggedArray.fromcounts(dihad_evts.astype(int), wpartons_dw[wpartons_dw.charge < 0][dihad_evts].flatten()),
    )

        # make Had/Lep decaying objects for SEMILEP events
    top_isLep = (GenTop.decaytype == 1)[semilep_evts].flatten()
    Lep_Top = JaggedCandidateArray.candidatesfromcounts(
        counts = semilep_evts.astype(int),
        pt       = np.where(top_isLep, GenTop[semilep_evts].pt.flatten(), GenTbar[semilep_evts].pt.flatten()),
        eta      = np.where(top_isLep, GenTop[semilep_evts].eta.flatten(), GenTbar[semilep_evts].eta.flatten()),
        phi      = np.where(top_isLep, GenTop[semilep_evts].phi.flatten(), GenTbar[semilep_evts].phi.flatten()),
        mass     = np.where(top_isLep, GenTop[semilep_evts].mass.flatten(), GenTbar[semilep_evts].mass.flatten()),
        charge   = np.where(top_isLep, GenTop[semilep_evts].charge.flatten(), GenTbar[semilep_evts].charge.flatten()),
        decaytype= np.where(top_isLep, GenTop[semilep_evts].decaytype.flatten(), GenTbar[semilep_evts].decaytype.flatten()),
    )
    Had_Top = JaggedCandidateArray.candidatesfromcounts(
        counts = semilep_evts.astype(int),
        pt       = np.where(~top_isLep, GenTop[semilep_evts].pt.flatten(), GenTbar[semilep_evts].pt.flatten()),
        eta      = np.where(~top_isLep, GenTop[semilep_evts].eta.flatten(), GenTbar[semilep_evts].eta.flatten()),
        phi      = np.where(~top_isLep, GenTop[semilep_evts].phi.flatten(), GenTbar[semilep_evts].phi.flatten()),
        mass     = np.where(~top_isLep, GenTop[semilep_evts].mass.flatten(), GenTbar[semilep_evts].mass.flatten()),
        charge   = np.where(~top_isLep, GenTop[semilep_evts].charge.flatten(), GenTbar[semilep_evts].charge.flatten()),
        decaytype= np.where(~top_isLep, GenTop[semilep_evts].decaytype.flatten(), GenTbar[semilep_evts].decaytype.flatten()),
    )
    Lep_B = JaggedCandidateArray.candidatesfromcounts(
        counts = semilep_evts.astype(int),
        pt       = np.where(top_isLep, GenB[semilep_evts].pt.flatten(), GenBbar[semilep_evts].pt.flatten()),
        eta      = np.where(top_isLep, GenB[semilep_evts].eta.flatten(), GenBbar[semilep_evts].eta.flatten()),
        phi      = np.where(top_isLep, GenB[semilep_evts].phi.flatten(), GenBbar[semilep_evts].phi.flatten()),
        mass     = np.where(top_isLep, GenB[semilep_evts].mass.flatten(), GenBbar[semilep_evts].mass.flatten()),
        charge   = np.where(top_isLep, GenB[semilep_evts].charge.flatten(), GenBbar[semilep_evts].charge.flatten()),
        decaytype= np.where(top_isLep, GenB[semilep_evts].decaytype.flatten(), GenBbar[semilep_evts].decaytype.flatten()),
    )
    Had_B = JaggedCandidateArray.candidatesfromcounts(
        counts = semilep_evts.astype(int),
        pt       = np.where(~top_isLep, GenB[semilep_evts].pt.flatten(), GenBbar[semilep_evts].pt.flatten()),
        eta      = np.where(~top_isLep, GenB[semilep_evts].eta.flatten(), GenBbar[semilep_evts].eta.flatten()),
        phi      = np.where(~top_isLep, GenB[semilep_evts].phi.flatten(), GenBbar[semilep_evts].phi.flatten()),
        mass     = np.where(~top_isLep, GenB[semilep_evts].mass.flatten(), GenBbar[semilep_evts].mass.flatten()),
        charge   = np.where(~top_isLep, GenB[semilep_evts].charge.flatten(), GenBbar[semilep_evts].charge.flatten()),
        decaytype= np.where(~top_isLep, GenB[semilep_evts].decaytype.flatten(), GenBbar[semilep_evts].decaytype.flatten()),
    )
    
    SEMILEP_evts = awkward.Table(
        TTbar     = awkward.JaggedArray.fromcounts(semilep_evts.astype(int), GenTTbar[semilep_evts].flatten()),
        THad   = awkward.JaggedArray.fromcounts(semilep_evts.astype(int), Had_Top.flatten()),
        TLep   = awkward.JaggedArray.fromcounts(semilep_evts.astype(int), Lep_Top.flatten()),
        BHad     = awkward.JaggedArray.fromcounts(semilep_evts.astype(int), Had_B.flatten()),
        BLep     = awkward.JaggedArray.fromcounts(semilep_evts.astype(int), Lep_B.flatten()),
        WHad     = awkward.JaggedArray.fromcounts(semilep_evts.astype(int), GenW[GenW.decaytype == 2][semilep_evts].flatten()),
        WLep     = awkward.JaggedArray.fromcounts(semilep_evts.astype(int), GenW[GenW.decaytype == 1][semilep_evts].flatten()),
        Lepton    = awkward.JaggedArray.fromcounts(semilep_evts.astype(int), charged_leps[semilep_evts].flatten()),
        Nu        = awkward.JaggedArray.fromcounts(semilep_evts.astype(int), neutral_leps[semilep_evts].flatten()),
        WJa = awkward.JaggedArray.fromcounts(semilep_evts.astype(int), First[semilep_evts].flatten()),
        WJb= awkward.JaggedArray.fromcounts(semilep_evts.astype(int), Second[semilep_evts].flatten()),
        Up_Had    = awkward.JaggedArray.fromcounts(semilep_evts.astype(int), wpartons_up[semilep_evts].flatten()),
        Down_Had  = awkward.JaggedArray.fromcounts(semilep_evts.astype(int), wpartons_dw[semilep_evts].flatten()),
    )


    GenTTbar = awkward.Table(
        SL = SEMILEP_evts,
        DL = DILEP_evts,
        Had= HAD_evts,
    )

    return GenTTbar
