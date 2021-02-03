import numpy as np
from pdb import set_trace
import coffea.processor as processor
import awkward
import python.GenObjects as genobj
from coffea.analysis_objects import JaggedCandidateArray
import numba

remove_fast = lambda x : x.split('fast_')[-1]

@numba.njit
def find_Wdecay_prod_inds(wminus_idx_array, wplus_idx_array, gp_offsets, gp_momInds, wdecay_prod_momIdx_array):
    ## finds particle Idx of first W boson for input W decay products
    wdecay_firstW_idx = np.zeros(wdecay_prod_momIdx_array.size)

    for wdecay_idx in range(wdecay_prod_momIdx_array.size):
        wdecay_momIdx = wdecay_prod_momIdx_array[wdecay_idx]
        while (wdecay_momIdx != wplus_idx_array[wdecay_idx]) and (wdecay_momIdx != wminus_idx_array[wdecay_idx]):
            wdecay_momIdx = gp_momInds[gp_offsets[wdecay_idx]+wdecay_momIdx]
        wdecay_firstW_idx[wdecay_idx] = wdecay_momIdx
    
    return wdecay_firstW_idx



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

        # add mother pdgId to genparts
    motherIdx = awkward.JaggedArray.fromcounts(df.nGenPart, df.GenPart_genPartIdxMother)
    pdgId = awkward.JaggedArray.fromcounts(df.nGenPart, df.GenPart_pdgId)
    hasmother = ((0 <= genParts.momIdx) & (genParts.momIdx < genParts.counts))
    motherpdgId = genParts.momIdx.zeros_like()
    motherpdgId[hasmother] = genParts.pdgId[genParts.momIdx[hasmother]]
    genParts['mompdgId'] = motherpdgId
    genParts['hasmother'] = hasmother
        # add parton indices to genparts
    inds_list = [np.arange(count) for count in df.nGenPart]
    part_inds = np.concatenate(inds_list, axis=None)
    genParts['Idx'] = awkward.JaggedArray.fromcounts(df.nGenPart, part_inds)
        # add decaytype (0 is INVALID, 1 is LEPTONIC, 2 is HADRONIC) as placeholder
    genParts['decaytype'] = genParts.momIdx.zeros_like()

    return genParts

def process_lheParts(df):
    lheParts = JaggedCandidateArray.candidatesfromcounts(
        df.nLHEPart,
        pt=df.LHEPart_pt,
        eta=df.LHEPart_eta,
        phi=df.LHEPart_phi,
        mass=df.LHEPart_mass,
        pdgId=df.LHEPart_pdgId,
        #status=df.LHEPart_status,
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

    genparts = df.genParts

        # get tops, defined as last copy
    is_last_copy = genparts.statusFlags >> 13 & 1 == 1
    gen_tops = genparts[(is_last_copy) & (genparts.pdgId == 6)]
    gen_tbars = genparts[(is_last_copy) & (genparts.pdgId == -6)]


        # get direct top decay products (will be first copy)
    is_first_copy = genparts.statusFlags >> 12 & 1 == 1
    is_hard_process = genparts.statusFlags >> 7 & 1 == 1
    hard_gps = genparts[is_first_copy & is_hard_process]
    abspdg = abs(hard_gps.pdgId)
    sgn = np.sign(hard_gps.pdgId)

    gen_bs = hard_gps[(hard_gps.pdgId == 5) & (hard_gps.mompdgId == 6)]
    gen_bbars = hard_gps[(hard_gps.pdgId == -5) & (hard_gps.mompdgId == -6)]
    gen_wplus = hard_gps[(hard_gps.pdgId == 24) & (hard_gps.mompdgId == 6)]
    gen_wminus = hard_gps[(hard_gps.pdgId == -24) & (hard_gps.mompdgId == -6)]

    gen_wpartons_up = hard_gps[(np.mod(hard_gps.pdgId, 2) == 0) & (abspdg < 6) & (hard_gps.mompdgId == sgn * 24)]
    gen_wpartons_dw = hard_gps[(np.mod(hard_gps.pdgId, 2) == 1) & (abspdg < 6) & (hard_gps.mompdgId == sgn * -24)]

    gen_charged_leps = hard_gps[( (abspdg == 11) | (abspdg == 13) | (abspdg == 15) ) & (hard_gps.mompdgId == sgn * -24)]
    gen_neutral_leps = hard_gps[( (abspdg == 12) | (abspdg == 14) | (abspdg == 16) ) & (hard_gps.mompdgId == sgn * 24)]
    gen_taus = hard_gps[(abspdg == 15) & (hard_gps.mompdgId == sgn * -24)]
    gen_nu_taus = hard_gps[(abspdg == 16) & (hard_gps.mompdgId == sgn * 24)]

    # find Idx of first W boson for up/dw partons and leptons for reconstructing gen W bosons
        # wparton up
    Wminus_WPup_idx_array = np.repeat(gen_wminus.Idx.content, gen_wpartons_up.counts)
    Wplus_WPup_idx_array = np.repeat(gen_wplus.Idx.content, gen_wpartons_up.counts)
    gp_offsets_WPup_array = np.repeat(genparts.offsets[:-1], gen_wpartons_up.counts)
    wpartons_up_momIdx_array = gen_wpartons_up.momIdx.content

        # wparton dw
    Wminus_WPdw_idx_array = np.repeat(gen_wminus.Idx.content, gen_wpartons_dw.counts)
    Wplus_WPdw_idx_array = np.repeat(gen_wplus.Idx.content, gen_wpartons_dw.counts)
    gp_offsets_WPdw_array = np.repeat(genparts.offsets[:-1], gen_wpartons_dw.counts)
    wpartons_dw_momIdx_array = gen_wpartons_dw.momIdx.content

        # charged lep
    Wminus_ChLep_idx_array = np.repeat(gen_wminus.Idx.content, gen_charged_leps.counts)
    Wplus_ChLep_idx_array = np.repeat(gen_wplus.Idx.content, gen_charged_leps.counts)
    gp_offsets_ChLep_array = np.repeat(genparts.offsets[:-1], gen_charged_leps.counts)
    ch_lep_momIdx_array = gen_charged_leps.momIdx.content

        # neutral lep
    Wminus_NuLep_idx_array = np.repeat(gen_wminus.Idx.content, gen_neutral_leps.counts)
    Wplus_NuLep_idx_array = np.repeat(gen_wplus.Idx.content, gen_neutral_leps.counts)
    gp_offsets_NuLep_array = np.repeat(genparts.offsets[:-1], gen_neutral_leps.counts)
    nu_lep_momIdx_array = gen_neutral_leps.momIdx.content

    gp_momInds_array = genparts.momIdx.content
    
    #set_trace()
    wpart_up_firstW_idx = find_Wdecay_prod_inds(Wminus_WPup_idx_array, Wplus_WPup_idx_array, gp_offsets_WPup_array, gp_momInds_array, wpartons_up_momIdx_array)
    wpart_dw_firstW_idx = find_Wdecay_prod_inds(Wminus_WPdw_idx_array, Wplus_WPdw_idx_array, gp_offsets_WPdw_array, gp_momInds_array, wpartons_dw_momIdx_array)
    charged_lep_firstW_idx = find_Wdecay_prod_inds(Wminus_ChLep_idx_array, Wplus_ChLep_idx_array, gp_offsets_ChLep_array, gp_momInds_array, ch_lep_momIdx_array)
    neutral_lep_firstW_idx = find_Wdecay_prod_inds(Wminus_NuLep_idx_array, Wplus_NuLep_idx_array, gp_offsets_NuLep_array, gp_momInds_array, nu_lep_momIdx_array)

    if not np.array_equal(wpart_up_firstW_idx, wpart_dw_firstW_idx):
        raise ValueError("Wpartons not matched to same W!")
    if not np.array_equal(charged_lep_firstW_idx, neutral_lep_firstW_idx):
        raise ValueError("Leptons not matched to same W!")



            # get direct tau decay products from hard processes (subset of gen_taus events above)
    isDirectHardProcessTauDecayProduct = genparts.statusFlags >> 10 & 1 == 1
    tau_decay_prods = genparts[isDirectHardProcessTauDecayProduct]

        # only need decays to leptons (e/mu, tau nu) for event classification
    tau_TO_tau_nu = tau_decay_prods[np.abs(tau_decay_prods.pdgId) == 16]
    tau_TO_charged_lep = tau_decay_prods[(np.abs(tau_decay_prods.pdgId) == 11) | (np.abs(tau_decay_prods.pdgId) == 13)]
    tau_TO_neutral_lep = tau_decay_prods[(np.abs(tau_decay_prods.pdgId) == 12) | (np.abs(tau_decay_prods.pdgId) == 14)]

        # set decaytype for gen taus
    charged_lep_decaytype_array = gen_charged_leps.decaytype.content

            # semilep evts
    lep_jets_evts = (gen_charged_leps.counts == 1) & (gen_taus.counts == 0) & (gen_wpartons_up.counts == 1) & (gen_wpartons_dw.counts == 1)
    tau_jets_evts = (gen_charged_leps.counts == 1) & (gen_taus.counts == 1) & (gen_wpartons_up.counts == 1) & (gen_wpartons_dw.counts == 1)
    semilep_evts = (lep_jets_evts | tau_jets_evts)
    sl_evts_mask = np.repeat(semilep_evts, gen_charged_leps.counts)
    semilep_decaytype_array = np.zeros(semilep_evts.sum(), dtype=int)
                # tau -> l
    semilep_tau_leptonic_decay = (tau_jets_evts) &  (tau_TO_charged_lep.counts == 1) & (tau_TO_neutral_lep.counts == 1) & (tau_TO_tau_nu.counts == 1)
    semilep_tau_hadronic_decay = (tau_jets_evts) & (~semilep_tau_leptonic_decay)
    semilep_decaytype_array[semilep_tau_leptonic_decay[semilep_evts]] = np.ones(semilep_tau_leptonic_decay.sum(), dtype=int)
    semilep_decaytype_array[semilep_tau_hadronic_decay[semilep_evts]] = np.ones(semilep_tau_hadronic_decay.sum(), dtype=int)*2

        # set charged_lep_decatype_array for semileptonic events
    charged_lep_decaytype_array[sl_evts_mask] = semilep_decaytype_array


            # dilep evts
    lep_lep_evts = (gen_charged_leps.counts == 2) & (gen_taus.counts == 0) & (gen_wpartons_up.counts == 0) & (gen_wpartons_dw.counts == 0)
    lep_tau_evts = (gen_charged_leps.counts == 2) & (gen_taus.counts == 1) & (gen_wpartons_up.counts == 0) & (gen_wpartons_dw.counts == 0)
    tau_tau_evts = (gen_charged_leps.counts == 2) & (gen_taus.counts == 2) & (gen_wpartons_up.counts == 0) & (gen_wpartons_dw.counts == 0)
    dilep_evts = (lep_lep_evts | lep_tau_evts | tau_tau_evts)
    dl_evts_mask = np.repeat(dilep_evts, gen_charged_leps.counts)
    dilep_decaytype_array = np.zeros((dilep_evts.sum(), 2), dtype=int)
                # tau + tau
                    # tau + tau -> ll
    dilep_TauTau_ll_decay = (tau_tau_evts) & (tau_TO_charged_lep.counts == 2) & (tau_TO_neutral_lep.counts == 2) & (tau_TO_tau_nu.counts == 2)
    dilep_decaytype_array[dilep_TauTau_ll_decay[dilep_evts]] = np.ones((dilep_TauTau_ll_decay.sum(), 2), dtype=int)
                    # tau + tau -> hh
    dilep_TauTau_hh_decay = (tau_tau_evts) & ((tau_TO_charged_lep.counts + tau_TO_neutral_lep.counts) == 0) & (tau_TO_tau_nu.counts == 2)
    dilep_decaytype_array[dilep_TauTau_hh_decay[dilep_evts]] = np.ones((dilep_TauTau_hh_decay.sum(), 2), dtype=int)*2
                    # tau + tau -> lh
    dilep_TauTau_lh_decay = (tau_tau_evts) & ~(dilep_TauTau_ll_decay | dilep_TauTau_hh_decay)
                        # set index corresponding to leptonically decaying tau to 1, default array is set to 2
    dl_TauTau_to_lh_decaytype_array = np.ones(dilep_TauTau_lh_decay.sum()*2, dtype=int)*2
    lep_tau_mask = (np.repeat(tau_TO_charged_lep[dilep_TauTau_lh_decay].mompdgId.flatten(), 2) == gen_charged_leps[dilep_TauTau_lh_decay].pdgId.flatten())
    dl_TauTau_to_lh_decaytype_array[lep_tau_mask] = np.ones(lep_tau_mask.sum(), dtype=int)
    dilep_decaytype_array[dilep_TauTau_lh_decay[dilep_evts]] = dl_TauTau_to_lh_decaytype_array.reshape(dilep_TauTau_lh_decay.sum(), 2)

                # lep + tau
                    # tau -> l
    dilep_LepTau_l_decay = (lep_tau_evts) & (tau_TO_charged_lep.counts == 1) & (tau_TO_neutral_lep.counts == 1) & (tau_TO_tau_nu.counts == 1)
                        # set index corresponding to tau to 1
    dl_LepTau_to_Lep_decaytype_array = np.zeros(dilep_LepTau_l_decay.sum()*2, dtype=int)
    dl_LepTau_to_Lep_decaytype_array[(np.abs(gen_charged_leps[dilep_LepTau_l_decay].pdgId) == 15).flatten()] = np.ones(dilep_LepTau_l_decay.sum(), dtype=int)
    dilep_decaytype_array[dilep_LepTau_l_decay[dilep_evts]] = dl_LepTau_to_Lep_decaytype_array.reshape(dilep_LepTau_l_decay.sum(), 2)
                    # tau -> h
    dilep_LepTau_h_decay = (lep_tau_evts) & ~(dilep_LepTau_l_decay)
                        # set index corresponding to tau to 2
    dl_LepTau_to_Had_decaytype_array = np.zeros(dilep_LepTau_h_decay.sum()*2, dtype=int)
    dl_LepTau_to_Had_decaytype_array[(np.abs(gen_charged_leps[dilep_LepTau_h_decay].pdgId) == 15).flatten()] = np.ones(dilep_LepTau_h_decay.sum(), dtype=int)*2
    dilep_decaytype_array[dilep_LepTau_h_decay[dilep_evts]] = dl_LepTau_to_Had_decaytype_array.reshape(dilep_LepTau_h_decay.sum(), 2)

        # set charged_lep_decatype_array for semileptonic events
    charged_lep_decaytype_array[dl_evts_mask] = dilep_decaytype_array.flatten()


    # reconstruct gen objects with new charges
        # tops
    gen_tops_dict = {remove_fast(key) : gen_tops[key].flatten() for key in gen_tops.columns if key != 'p4'}
    gen_tops_dict['charge'] = (gen_tops.charge.ones_like()*(2./3.)).flatten()
    gen_tops = JaggedCandidateArray.candidatesfromcounts(
        counts = gen_tops.counts,
        **gen_tops_dict
    )
        # tbars
    gen_tbars_dict = {remove_fast(key) : gen_tbars[key].flatten() for key in gen_tbars.columns if key != 'p4'}
    gen_tbars_dict['charge'] = (gen_tbars.charge.ones_like()*(-2./3.)).flatten()
    gen_tbars = JaggedCandidateArray.candidatesfromcounts(
        counts = gen_tbars.counts,
        **gen_tbars_dict
    )
        # wplus
    gen_wplus_dict = {remove_fast(key) : gen_wplus[key].flatten() for key in gen_wplus.columns if key != 'p4'}
    gen_wplus_dict['charge'] = (gen_wplus.charge.ones_like()).flatten()
    gen_wplus = JaggedCandidateArray.candidatesfromcounts(
        counts = gen_wplus.counts,
        **gen_wplus_dict
    )
        # wminus
    gen_wminus_dict = {remove_fast(key) : gen_wminus[key].flatten() for key in gen_wminus.columns if key != 'p4'}
    gen_wminus_dict['charge'] = (gen_wminus.charge.ones_like()*(-1.)).flatten()
    gen_wminus = JaggedCandidateArray.candidatesfromcounts(
        counts = gen_wminus.counts,
        **gen_wminus_dict
    )
        # bs
    gen_bs_dict = {remove_fast(key) : gen_bs[key].flatten() for key in gen_bs.columns if key != 'p4'}
    gen_bs_dict['charge'] = (gen_bs.charge.ones_like()*(-1./3.)).flatten()
    gen_bs = JaggedCandidateArray.candidatesfromcounts(
        counts = gen_bs.counts,
        **gen_bs_dict
    )
        # bbars
    gen_bbars_dict = {remove_fast(key) : gen_bbars[key].flatten() for key in gen_bbars.columns if key != 'p4'}
    gen_bbars_dict['charge'] = (gen_bbars.charge.ones_like()*(1./3.)).flatten()
    gen_bbars = JaggedCandidateArray.candidatesfromcounts(
        counts = gen_bbars.counts,
        **gen_bbars_dict
    )

        # wpartons up
    gen_wpartons_up_dict = {remove_fast(key) : gen_wpartons_up[key].flatten() for key in gen_wpartons_up.columns if key != 'p4'}
    gen_wpartons_up_dict['charge'] = (sgn[(np.mod(hard_gps.pdgId, 2) == 0) & (abspdg < 6) & (hard_gps.mompdgId == sgn * 24)]*(2./3.)).flatten()
    gen_wpartons_up_dict['firstW_momIdx'] = wpart_up_firstW_idx
    gen_wpartons_up = JaggedCandidateArray.candidatesfromcounts(
        counts = gen_wpartons_up.counts,
        **gen_wpartons_up_dict
    )

        # wpartons dw
    gen_wpartons_dw_dict = {remove_fast(key) : gen_wpartons_dw[key].flatten() for key in gen_wpartons_dw.columns if key != 'p4'}
    gen_wpartons_dw_dict['charge'] = (sgn[(np.mod(hard_gps.pdgId, 2) == 1) & (abspdg < 6) & (hard_gps.mompdgId == sgn * -24)]*(-1./3.)).flatten()
    gen_wpartons_dw_dict['firstW_momIdx'] = wpart_dw_firstW_idx
    gen_wpartons_dw = JaggedCandidateArray.candidatesfromcounts(
        counts = gen_wpartons_dw.counts,
        **gen_wpartons_dw_dict
    )

        # neutral leps
    gen_neutral_leps_dict = {remove_fast(key) : gen_neutral_leps[key].flatten() for key in gen_neutral_leps.columns if key != 'p4'}
    gen_neutral_leps_dict['charge'] = (gen_neutral_leps.charge.zeros_like()).flatten()
    gen_neutral_leps_dict['firstW_momIdx'] = neutral_lep_firstW_idx
    gen_neutral_leps = JaggedCandidateArray.candidatesfromcounts(
        counts = gen_neutral_leps.counts,
        **gen_neutral_leps_dict
    )

        # reconstruct gen_charged_leps with new decaytype values and charge
    gen_charged_leps_dict = {remove_fast(key) : gen_charged_leps[key].flatten() for key in gen_charged_leps.columns if key != 'p4'}
    gen_charged_leps_dict['decaytype'] = charged_lep_decaytype_array
    gen_charged_leps_dict['charge'] = (sgn[( (abspdg == 11) | (abspdg == 13) | (abspdg == 15) ) & (hard_gps.mompdgId == sgn * -24)]* -1).flatten()
    gen_charged_leps_dict['firstW_momIdx'] = charged_lep_firstW_idx
    gen_charged_leps = JaggedCandidateArray.candidatesfromcounts(
        counts = gen_charged_leps.counts,
        **gen_charged_leps_dict
    )

        # jagged array of (wplus, wminus)
        # initialize variables to become attributes of GenW objecs
    genWpt       = np.zeros(gen_wplus.counts.sum()+gen_wminus.counts.sum())
    genWeta      = np.zeros(gen_wplus.counts.sum()+gen_wminus.counts.sum())
    genWphi      = np.zeros(gen_wplus.counts.sum()+gen_wminus.counts.sum())
    genWmass     = np.zeros(gen_wplus.counts.sum()+gen_wminus.counts.sum())
    genWcharge   =np.zeros(gen_wplus.counts.sum()+gen_wminus.counts.sum())
    genWdecaytype=np.zeros(gen_wplus.counts.sum()+gen_wminus.counts.sum())

    wplus_inds = np.arange(genWpt.size)[(np.mod(np.arange(genWpt.size), 2) == 0)] # wplus are defined as first W
    wminus_inds = np.arange(genWpt.size)[(np.mod(np.arange(genWpt.size), 2) == 1)] # wminus are defined as second W
        # wplus
    genWpt[wplus_inds] = gen_wplus.pt.flatten()
    genWeta[wplus_inds] = gen_wplus.eta.flatten()
    genWphi[wplus_inds] = gen_wplus.phi.flatten()
    genWmass[wplus_inds] = gen_wplus.mass.flatten()
    genWcharge[wplus_inds] = gen_wplus.charge.flatten()
    #set_trace()
            ## set decaytype
    wplus_decaytype = np.zeros(gen_wplus.pt.content.size)
                # SL events
    SL_hadWplus = np.where(gen_wplus[semilep_evts].Idx.flatten() == gen_wpartons_up[semilep_evts].firstW_momIdx.flatten())[0]
    SL_lepWplus = np.where(gen_wplus[semilep_evts].Idx.flatten() == gen_charged_leps[semilep_evts].firstW_momIdx.flatten())[0]
    SL_wplus_decaytype = np.zeros(semilep_evts.sum())
    SL_wplus_decaytype[SL_hadWplus] = np.ones(SL_hadWplus.size)*2 # 0 is INVALID (shouldn't happen), 1 is LEPTONIC, 2 is HADRONIC
    SL_wplus_decaytype[SL_lepWplus] = np.ones(SL_lepWplus.size) # 0 is INVALID (shouldn't happen), 1 is LEPTONIC, 2 is HADRONIC
    wplus_decaytype[semilep_evts] = SL_wplus_decaytype
                # DiLep events
    wplus_decaytype[dilep_evts] = np.ones(dilep_evts.sum()) # 0 is INVALID (shouldn't happen), 1 is LEPTONIC, 2 is HADRONIC
                # Had events
    wplus_decaytype[~(dilep_evts | semilep_evts)] = np.ones((~(dilep_evts | semilep_evts)).sum())*2 # 0 is INVALID (shouldn't happen), 1 is LEPTONIC, 2 is HADRONIC
    genWdecaytype[wplus_inds] = wplus_decaytype
        # wminus
    genWpt[wminus_inds] = gen_wminus.pt.flatten()
    genWeta[wminus_inds] = gen_wminus.eta.flatten()
    genWphi[wminus_inds] = gen_wminus.phi.flatten()
    genWmass[wminus_inds] = gen_wminus.mass.flatten()
    genWcharge[wminus_inds] = gen_wminus.charge.flatten()
            ## set decaytype
    wminus_decaytype = np.zeros(gen_wminus.pt.content.size)
                # SL events
    SL_hadWminus = np.where(gen_wminus[semilep_evts].Idx.flatten() == gen_wpartons_up[semilep_evts].firstW_momIdx.flatten())[0]
    SL_lepWminus = np.where(gen_wminus[semilep_evts].Idx.flatten() == gen_charged_leps[semilep_evts].firstW_momIdx.flatten())[0]
    SL_wminus_decaytype = np.zeros(semilep_evts.sum())
    SL_wminus_decaytype[SL_hadWminus] = np.ones(SL_hadWminus.size)*2 # 0 is INVALID (shouldn't happen), 1 is LEPTONIC, 2 is HADRONIC
    SL_wminus_decaytype[SL_lepWminus] = np.ones(SL_lepWminus.size) # 0 is INVALID (shouldn't happen), 1 is LEPTONIC, 2 is HADRONIC
    wminus_decaytype[semilep_evts] = SL_wminus_decaytype
                # DiLep events
    wminus_decaytype[dilep_evts] = np.ones(dilep_evts.sum()) # 0 is INVALID (shouldn't happen), 1 is LEPTONIC, 2 is HADRONIC
                # Had events
    wminus_decaytype[~(dilep_evts | semilep_evts)] = np.ones((~(dilep_evts | semilep_evts)).sum())*2 # 0 is INVALID (shouldn't happen), 1 is LEPTONIC, 2 is HADRONIC
    genWdecaytype[wminus_inds] = wminus_decaytype

    Gen_W = JaggedCandidateArray.candidatesfromcounts(
        counts=gen_wplus.counts+gen_wminus.counts,
        pt       = genWpt,
        eta      = genWeta,
        phi      = genWphi,
        mass     = genWmass,
        charge   = genWcharge,
        decaytype= genWdecaytype,
    )


    #set_trace()
    GenObjects = awkward.Table(
        GenTop      = gen_tops,
        GenTbar     = gen_tbars,
        GenW        = Gen_W,
        GenB        = gen_bs,
        GenBbar     = gen_bbars,
        GenWPartsUp = gen_wpartons_up,
        GenWPartsDw = gen_wpartons_dw,
        ChargedLeps = gen_charged_leps,
        NeutralLeps = gen_neutral_leps,
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

    ttbar = genobj.from_collections(wpartons_up=GenObjs['GenWPartsUp'], wpartons_dw=GenObjs['GenWPartsDw'], charged_leps=GenObjs['ChargedLeps'], neutral_leps=GenObjs['NeutralLeps'], bs=GenObjs['GenB'], bbars=GenObjs['GenBbar'], wbosons=GenObjs['GenW'], tops=GenObjs['GenTop'], tbars=GenObjs['GenTbar'])

    #set_trace()
    return ttbar

