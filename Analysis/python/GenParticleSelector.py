import numpy as np
from pdb import set_trace
import awkward as ak


def select_normal(genparts, w_decay_momid):

    columns = ["pt", "eta", "phi", "mass", "pdgId", "charge", "decaytype"]

        # get tops, defined as last copy
    gen_tops = genparts[(genparts.hasFlags(['isLastCopy'])) & (genparts.pdgId == 6)]
    gen_tops['charge'] = ak.ones_like(gen_tops.pt)*(2./3.)
    gen_tbars = genparts[(genparts.hasFlags(['isLastCopy'])) & (genparts.pdgId == -6)]
    gen_tbars['charge'] = ak.ones_like(gen_tbars.pt)*(-2./3.)

        # get direct top decay products (children are "isHardProcess" and "isFirstCopy")
    gen_bs = ak.flatten(gen_tops.children[(gen_tops.children.pdgId == 5)], axis=2)
    gen_bs['charge'] = ak.ones_like(gen_bs.pt)*(-1./3.)
    gen_bbars = ak.flatten(gen_tbars.children[(gen_tbars.children.pdgId == -5)], axis=2)
    gen_bbars['charge'] = ak.ones_like(gen_bbars.pt)*(1./3.)
    gen_wplus = ak.flatten(gen_tops.children[(gen_tops.children.pdgId == w_decay_momid)], axis=2)
    gen_wplus['charge'] = ak.ones_like(gen_wplus.pt)
    gen_wminus = ak.flatten(gen_tbars.children[(gen_tbars.children.pdgId == -1*w_decay_momid)], axis=2)
    gen_wminus['charge'] = ak.ones_like(gen_wminus.pt)*(-1.)

        # get w decay products
            # get last copy of W bosons
    last_gen_wplus = ak.flatten(gen_tops.distinctChildren[(gen_tops.distinctChildren.pdgId == w_decay_momid) & (gen_tops.distinctChildren.hasFlags(['isLastCopy']))], axis=2)
    last_gen_wminus = ak.flatten(gen_tbars.distinctChildren[(gen_tbars.distinctChildren.pdgId == -1*w_decay_momid) & (gen_tbars.distinctChildren.hasFlags(['isLastCopy']))], axis=2)

    gen_wplus_decaytype = np.zeros(len(gen_wplus))
    gen_wminus_decaytype = np.zeros(len(gen_wminus))
            # up/down partons from last W+
    gen_wpartons_up_fromWplus = ak.flatten( last_gen_wplus.children[(np.mod(last_gen_wplus.children.pdgId, 2) == 0) & (np.abs(last_gen_wplus.children.pdgId) < 6)], axis=2)
    gen_wpartons_up_fromWplus['charge'] = np.sign(gen_wpartons_up_fromWplus.pdgId)*(2./3.)
    gen_wplus_decaytype[ak.num(gen_wpartons_up_fromWplus) > 0] = np.ones(ak.sum(ak.num(gen_wpartons_up_fromWplus) > 0))*2
    gen_wpartons_dw_fromWplus = ak.flatten( last_gen_wplus.children[(np.mod(last_gen_wplus.children.pdgId, 2) == 1) & (np.abs(last_gen_wplus.children.pdgId) < 6)], axis=2)
    gen_wpartons_dw_fromWplus['charge'] = np.sign(gen_wpartons_up_fromWplus.pdgId)*(-1./3.)
            # up/down partons from last W-
    gen_wpartons_up_fromWminus = ak.flatten( last_gen_wminus.children[(np.mod(last_gen_wminus.children.pdgId, 2) == 0) & (np.abs(last_gen_wminus.children.pdgId) < 6)], axis=2)
    gen_wpartons_up_fromWminus['charge'] = np.sign(gen_wpartons_up_fromWminus.pdgId)*(2./3.)
    gen_wminus_decaytype[ak.num(gen_wpartons_up_fromWminus) > 0] = np.ones(ak.sum(ak.num(gen_wpartons_up_fromWminus) > 0))*2
    gen_wpartons_dw_fromWminus = ak.flatten( last_gen_wminus.children[(np.mod(last_gen_wminus.children.pdgId, 2) == 1) & (np.abs(last_gen_wminus.children.pdgId) < 6)], axis=2)
    gen_wpartons_dw_fromWminus['charge'] = np.sign(gen_wpartons_up_fromWminus.pdgId)*(-1./3.)

            # charged leps from last W+
    gen_charged_leps_fromWplus = ak.flatten( last_gen_wplus.children[(np.abs(last_gen_wplus.children.pdgId) == 11) | (np.abs(last_gen_wplus.children.pdgId) == 13) | (np.abs(last_gen_wplus.children.pdgId) == 15)], axis=2)
    gen_charged_leps_fromWplus['charge'] = ak.ones_like(gen_charged_leps_fromWplus.pdgId)
    gen_wplus_decaytype[ak.num(gen_charged_leps_fromWplus) > 0] = np.ones(ak.sum(ak.num(gen_charged_leps_fromWplus) > 0))
        # add decaytype (0 is INVALID, 1 is LEPTONIC, 2 is HADRONIC)
    gen_wplus['decaytype'] = ak.unflatten(gen_wplus_decaytype, ak.num(gen_wplus))
    gen_tops['decaytype'] = gen_wplus['decaytype'] # set decaytype for tops
    gen_bs['decaytype'] = gen_wplus['decaytype'] # set decaytype for bs
            # neutral leps from last W+
    gen_neutral_leps_fromWplus = ak.flatten( last_gen_wplus.children[(np.abs(last_gen_wplus.children.pdgId) == 12) | (np.abs(last_gen_wplus.children.pdgId) == 14) | (np.abs(last_gen_wplus.children.pdgId) == 16)], axis=2)
    gen_neutral_leps_fromWplus['charge'] = ak.zeros_like(gen_neutral_leps_fromWplus.pt)

            # charged leps from last W-
    gen_charged_leps_fromWminus = ak.flatten( last_gen_wminus.children[(np.abs(last_gen_wminus.children.pdgId) == 11) | (np.abs(last_gen_wminus.children.pdgId) == 13) | (np.abs(last_gen_wminus.children.pdgId) == 15)], axis=2)
    gen_charged_leps_fromWminus['charge'] = ak.ones_like(gen_charged_leps_fromWminus.pdgId)*(-1.)
    gen_wminus_decaytype[ak.num(gen_charged_leps_fromWminus) > 0] = np.ones(ak.sum(ak.num(gen_charged_leps_fromWminus) > 0))
        # add decaytype (0 is INVALID, 1 is LEPTONIC, 2 is HADRONIC)
    gen_wminus['decaytype'] = ak.unflatten(gen_wminus_decaytype, ak.num(gen_wminus))
    gen_tbars['decaytype'] = gen_wminus['decaytype'] # set decaytype for tbars
    gen_bbars['decaytype'] = gen_wminus['decaytype'] # set decaytype for bbars
            # neutral leps from last W-
    gen_neutral_leps_fromWminus = ak.flatten( last_gen_wminus.children[(np.abs(last_gen_wminus.children.pdgId) == 12) | (np.abs(last_gen_wminus.children.pdgId) == 14) | (np.abs(last_gen_wminus.children.pdgId) == 16)], axis=2)
    gen_neutral_leps_fromWminus['charge'] = ak.zeros_like(gen_neutral_leps_fromWminus.pt)

    gen_wpartons_up = ak.Array({}, with_name="PtEtaPhiMLorentzVector")
    gen_wpartons_dw = ak.Array({}, with_name="PtEtaPhiMLorentzVector")
    gen_charged_leps = ak.Array({}, with_name="PtEtaPhiMLorentzVector")
    gen_neutral_leps = ak.Array({}, with_name="PtEtaPhiMLorentzVector")
    for column in columns:
        if column == 'decaytype': continue
        gen_wpartons_up[column] = ak.flatten(ak.concatenate( [gen_wpartons_up_fromWplus[column], gen_wpartons_up_fromWminus[column]] , axis=1)) # (up-type partons from W+, W-)
        gen_wpartons_dw[column] = ak.flatten(ak.concatenate( [gen_wpartons_dw_fromWplus[column], gen_wpartons_dw_fromWminus[column]] , axis=1)) # (dw-type partons from W+, W-)
        gen_charged_leps[column] = ak.flatten(ak.concatenate( [gen_charged_leps_fromWplus[column], gen_charged_leps_fromWminus[column]] , axis=1)) # (charged leps from W+, W-)
        gen_neutral_leps[column] = ak.flatten(ak.concatenate( [gen_neutral_leps_fromWplus[column], gen_neutral_leps_fromWminus[column]] , axis=1)) # (neutral leps from W+, W-)

    gen_wpartons_up = ak.unflatten(gen_wpartons_up, ak.num(gen_wpartons_up_fromWplus)+ak.num(gen_wpartons_up_fromWminus) )
    gen_wpartons_dw = ak.unflatten(gen_wpartons_dw, ak.num(gen_wpartons_dw_fromWplus)+ak.num(gen_wpartons_dw_fromWminus) )
    gen_charged_leps = ak.unflatten(gen_charged_leps, ak.num(gen_charged_leps_fromWplus)+ak.num(gen_charged_leps_fromWminus))
    gen_neutral_leps = ak.unflatten(gen_neutral_leps, ak.num(gen_neutral_leps_fromWplus)+ak.num(gen_neutral_leps_fromWminus))
    gen_taus = gen_charged_leps[np.abs(gen_charged_leps.pdgId) == 15]
    gen_nu_taus = gen_neutral_leps[np.abs(gen_neutral_leps.pdgId) == 16]

            # fully hadronic evts
    had_evts = (ak.num(gen_charged_leps) == 0) & (ak.num(gen_neutral_leps) == 0) & (ak.num(gen_wpartons_up) == 2) & (ak.num(gen_wpartons_dw) == 2)

    #set_trace()
            # get direct tau decay products from hard processes (subset of gen_taus events above)
    tau_decay_prods = genparts[genparts.hasFlags(['isDirectHardProcessTauDecayProduct'])]

        # only need decays to leptons (e/mu, tau nu) for event classification
    tau_TO_tau_nu = tau_decay_prods[np.abs(tau_decay_prods.pdgId) == 16]
    tau_TO_charged_lep = tau_decay_prods[(np.abs(tau_decay_prods.pdgId) == 11) | (np.abs(tau_decay_prods.pdgId) == 13)]
    tau_TO_neutral_lep = tau_decay_prods[(np.abs(tau_decay_prods.pdgId) == 12) | (np.abs(tau_decay_prods.pdgId) == 14)]

        # set decaytype for gen taus
    charged_lep_decaytype_array = ak.to_numpy(ak.flatten(ak.zeros_like(gen_charged_leps['pt'])))

            # semilep evts
    semilep_evts = (ak.num(gen_charged_leps) == 1) & (ak.num(gen_neutral_leps) == 1) & (ak.num(gen_wpartons_up) == 1) & (ak.num(gen_wpartons_dw) == 1)
    tau_jets_evts = semilep_evts & (ak.num(gen_taus) == 1)
    sl_evts_mask = np.repeat(ak.to_numpy(semilep_evts), ak.to_numpy(ak.num(gen_charged_leps)))
    semilep_decaytype_array = np.zeros(ak.to_numpy(semilep_evts).sum(), dtype=int)
                # tau -> l
    semilep_tau_leptonic_decay = (tau_jets_evts) &  (ak.num(tau_TO_charged_lep) == 1) & (ak.num(tau_TO_neutral_lep) == 1) & (ak.num(tau_TO_tau_nu) == 1)
    semilep_tau_hadronic_decay = (tau_jets_evts) & (~semilep_tau_leptonic_decay)
    semilep_decaytype_array[semilep_tau_leptonic_decay[semilep_evts]] = np.ones(ak.to_numpy(semilep_tau_leptonic_decay).sum(), dtype=int)
    semilep_decaytype_array[semilep_tau_hadronic_decay[semilep_evts]] = np.ones(ak.to_numpy(semilep_tau_hadronic_decay).sum(), dtype=int)*2

        # set charged_lep_decatype_array for semileptonic events
    charged_lep_decaytype_array[sl_evts_mask] = semilep_decaytype_array


            # dilep evts
    dilep_evts = (ak.num(gen_charged_leps) == 2) & (ak.num(gen_neutral_leps) == 2) & (ak.num(gen_wpartons_up) == 0) & (ak.num(gen_wpartons_dw) == 0)
    lep_lep_evts = dilep_evts & (ak.num(gen_taus) == 0)
    lep_tau_evts = dilep_evts & (ak.num(gen_taus) == 1)
    tau_tau_evts = dilep_evts & (ak.num(gen_taus) == 2)
    dl_evts_mask = np.repeat(ak.to_numpy(dilep_evts), ak.to_numpy(ak.num(gen_charged_leps)))
    dilep_decaytype_array = np.zeros((ak.to_numpy(dilep_evts).sum(), 2), dtype=int)
                # tau + tau
                    # tau + tau -> ll
    dilep_TauTau_ll_decay = (tau_tau_evts) & (ak.num(tau_TO_charged_lep) == 2) & (ak.num(tau_TO_neutral_lep) == 2) & (ak.num(tau_TO_tau_nu) == 2)
    dilep_decaytype_array[dilep_TauTau_ll_decay[dilep_evts]] = np.ones((ak.to_numpy(dilep_TauTau_ll_decay).sum(), 2), dtype=int)
                    # tau + tau -> hh
    dilep_TauTau_hh_decay = (tau_tau_evts) & ((ak.num(tau_TO_charged_lep) + ak.num(tau_TO_neutral_lep)) == 0) & (ak.num(tau_TO_tau_nu) == 2)
    dilep_decaytype_array[dilep_TauTau_hh_decay[dilep_evts]] = np.ones((ak.to_numpy(dilep_TauTau_hh_decay).sum(), 2), dtype=int)*2
                    # tau + tau -> lh
    dilep_TauTau_lh_decay = (tau_tau_evts) & ~(dilep_TauTau_ll_decay | dilep_TauTau_hh_decay)
                        # set index corresponding to leptonically decaying tau to 1, default array is set to 2
    dl_TauTau_to_lh_decaytype_array = np.ones(ak.to_numpy(dilep_TauTau_lh_decay).sum()*2, dtype=int)*2
    lep_tau_mask = (np.repeat(ak.to_numpy(ak.flatten(tau_TO_charged_lep[dilep_TauTau_lh_decay].parent.pdgId)), 2) == ak.flatten(gen_charged_leps[dilep_TauTau_lh_decay].pdgId))
    dl_TauTau_to_lh_decaytype_array[lep_tau_mask] = np.ones(lep_tau_mask.sum(), dtype=int)
    dilep_decaytype_array[dilep_TauTau_lh_decay[dilep_evts]] = dl_TauTau_to_lh_decaytype_array.reshape(ak.to_numpy(dilep_TauTau_lh_decay).sum(), 2)

                # lep + tau
                    # tau -> l
    dilep_LepTau_l_decay = (lep_tau_evts) & (ak.num(tau_TO_charged_lep) == 1) & (ak.num(tau_TO_neutral_lep) == 1) & (ak.num(tau_TO_tau_nu) == 1)
                        # set index corresponding to tau to 1
    dl_LepTau_to_Lep_decaytype_array = np.zeros(ak.to_numpy(dilep_LepTau_l_decay).sum()*2, dtype=int)
    dl_LepTau_to_Lep_decaytype_array[ak.flatten(np.abs(gen_charged_leps[dilep_LepTau_l_decay].pdgId) == 15)] = np.ones(ak.sum(dilep_LepTau_l_decay), dtype=int)
    dilep_decaytype_array[dilep_LepTau_l_decay[dilep_evts]] = dl_LepTau_to_Lep_decaytype_array.reshape(ak.sum(dilep_LepTau_l_decay), 2)
                    # tau -> h
    dilep_LepTau_h_decay = (lep_tau_evts) & ~(dilep_LepTau_l_decay)
                        # set index corresponding to tau to 2
    dl_LepTau_to_Had_decaytype_array = np.zeros(ak.sum(dilep_LepTau_h_decay)*2, dtype=int)
    dl_LepTau_to_Had_decaytype_array[ak.flatten(np.abs(gen_charged_leps[dilep_LepTau_h_decay].pdgId) == 15)] = np.ones(ak.sum(dilep_LepTau_h_decay), dtype=int)*2
    dilep_decaytype_array[dilep_LepTau_h_decay[dilep_evts]] = dl_LepTau_to_Had_decaytype_array.reshape(ak.sum(dilep_LepTau_h_decay), 2)

        # set charged_lep_decatype_array for dileptonic events
    charged_lep_decaytype_array[dl_evts_mask] = dilep_decaytype_array.flatten()

        #  set charged lepton decaytype (defined for taus only, e/mu are 0) (1 is LEPTONIC, 2 is HADRONIC)
    gen_charged_leps['decaytype'] = ak.unflatten(charged_lep_decaytype_array, ak.num(gen_charged_leps))

    # make awkward arrays of (top decay prods, tbar decay prods)
    Gen_Top_Pairs = ak.Array({}, with_name="PtEtaPhiMLorentzVector")
    Gen_B_Pairs = ak.Array({}, with_name="PtEtaPhiMLorentzVector")
    Gen_W_Pairs = ak.Array({}, with_name="PtEtaPhiMLorentzVector")
    Gen_Wparton_Pairs = ak.Array({}, with_name="PtEtaPhiMLorentzVector")
    for column in columns:
        Gen_Top_Pairs[column] = ak.flatten(ak.concatenate( [gen_tops[column], gen_tbars[column]] , axis=1)) # (top, tbar)
        Gen_B_Pairs[column] = ak.flatten(ak.concatenate( [gen_bs[column], gen_bbars[column]] , axis=1)) # (b, bbar)
        Gen_W_Pairs[column] = ak.flatten(ak.concatenate( [gen_wplus[column], gen_wminus[column]] , axis=1)) # (W+, W-)
        if column != "decaytype":
        #if column is not "decaytype":
            Gen_Wparton_Pairs[column] = ak.flatten(ak.concatenate( [ak.pad_none(gen_wpartons_up[column], 1, axis=1), ak.pad_none(gen_wpartons_dw[column], 1, axis=1)] , axis=1)) # (up-type wpartons, down-type wpartons)

    Gen_Top_Pairs = ak.unflatten(Gen_Top_Pairs, ak.num(gen_tops)+ak.num(gen_tbars))
    Gen_B_Pairs = ak.unflatten(Gen_B_Pairs, ak.num(gen_bs)+ak.num(gen_bbars))
    Gen_W_Pairs = ak.unflatten(Gen_W_Pairs, ak.num(gen_wplus)+ak.num(gen_wminus))
    Gen_Wparton_Pairs = ak.unflatten(Gen_Wparton_Pairs, ak.num(ak.pad_none(gen_wpartons_up, 1, axis=1))+ak.num(ak.pad_none(gen_wpartons_dw, 1, axis=1)))
    if ak.sum(ak.num(gen_wpartons_up)+ak.num(gen_wpartons_dw)) > 0:
        Gen_Wparton_Pairs = Gen_Wparton_Pairs[ak.argsort(Gen_Wparton_Pairs["pt"], ascending=False)] # sort by pt

    Gen_TTbar = ak.Array({
        "pt" : (gen_tops+gen_tbars).pt,
        "eta": (gen_tops+gen_tbars).eta,
        "phi": (gen_tops+gen_tbars).phi,
        "mass":(gen_tops+gen_tbars).mass,
        "decaytype": gen_tops["decaytype"]+gen_tbars["decaytype"], # 0 is for INVALID, 2 for DILEP, 3 for SEMILEP, 4 for HADRONIC
    }, with_name="PtEtaPhiMLorentzVector")

    ## make "table" of gen objects for certain decays
        ## DILEP
    DILEP_evts = ak.zip({
        "TTbar"       : Gen_TTbar[dilep_evts] if ak.any(dilep_evts) else ak.unflatten(Gen_TTbar[dilep_evts], ak.values_astype(dilep_evts, int)),
        "Top"         : Gen_Top_Pairs[np.sign(Gen_Top_Pairs.charge) == 1][dilep_evts] if ak.any(dilep_evts) else ak.unflatten(Gen_Top_Pairs[np.sign(Gen_Top_Pairs.charge) == 1][dilep_evts], ak.values_astype(dilep_evts, int)),
        "Tbar"        : Gen_Top_Pairs[np.sign(Gen_Top_Pairs.charge) == -1][dilep_evts] if ak.any(dilep_evts) else ak.unflatten(Gen_Top_Pairs[np.sign(Gen_Top_Pairs.charge) == -1][dilep_evts], ak.values_astype(dilep_evts, int)),
        "B"           : Gen_B_Pairs[np.sign(Gen_B_Pairs.charge) == -1][dilep_evts] if ak.any(dilep_evts) else ak.unflatten(Gen_B_Pairs[np.sign(Gen_B_Pairs.charge) == -1][dilep_evts], ak.values_astype(dilep_evts, int)),
        "Bbar"        : Gen_B_Pairs[np.sign(Gen_B_Pairs.charge) == 1][dilep_evts] if ak.any(dilep_evts) else ak.unflatten(Gen_B_Pairs[np.sign(Gen_B_Pairs.charge) == 1][dilep_evts], ak.values_astype(dilep_evts, int)),
        "Wplus"       : Gen_W_Pairs[Gen_W_Pairs.charge == 1][dilep_evts] if ak.any(dilep_evts) else ak.unflatten(Gen_W_Pairs[Gen_W_Pairs.charge == 1][dilep_evts], ak.values_astype(dilep_evts, int)),
        "Wminus"      : Gen_W_Pairs[Gen_W_Pairs.charge == -1][dilep_evts] if ak.any(dilep_evts) else ak.unflatten(Gen_W_Pairs[Gen_W_Pairs.charge == -1][dilep_evts], ak.values_astype(dilep_evts, int)),
        "First_plus"  : gen_charged_leps[gen_charged_leps.charge > 0][dilep_evts] if ak.any(dilep_evts) else ak.unflatten(gen_charged_leps[gen_charged_leps.charge > 0][dilep_evts], ak.values_astype(dilep_evts, int)), # charged lepton always made leading
        "Second_plus" : gen_neutral_leps[gen_charged_leps.charge > 0][dilep_evts] if ak.any(dilep_evts) else ak.unflatten(gen_neutral_leps[gen_charged_leps.charge > 0][dilep_evts], ak.values_astype(dilep_evts, int)), # neutral lepton always made subleading
        "First_minus" : gen_charged_leps[gen_charged_leps.charge < 0][dilep_evts] if ak.any(dilep_evts) else ak.unflatten(gen_charged_leps[gen_charged_leps.charge < 0][dilep_evts], ak.values_astype(dilep_evts, int)), # charged lepton always made leading
        "Second_minus": gen_neutral_leps[gen_charged_leps.charge < 0][dilep_evts] if ak.any(dilep_evts) else ak.unflatten(gen_neutral_leps[gen_charged_leps.charge < 0][dilep_evts], ak.values_astype(dilep_evts, int)), # neutral lepton always made subleading
        "Up_plus"     : gen_neutral_leps[gen_charged_leps.charge > 0][dilep_evts] if ak.any(dilep_evts) else ak.unflatten(gen_neutral_leps[gen_charged_leps.charge > 0][dilep_evts], ak.values_astype(dilep_evts, int)), # same as second plus
        "Down_plus"   : gen_charged_leps[gen_charged_leps.charge > 0][dilep_evts] if ak.any(dilep_evts) else ak.unflatten(gen_charged_leps[gen_charged_leps.charge > 0][dilep_evts], ak.values_astype(dilep_evts, int)), # same as first plus
        "Up_minus"    : gen_neutral_leps[gen_charged_leps.charge < 0][dilep_evts] if ak.any(dilep_evts) else ak.unflatten(gen_neutral_leps[gen_charged_leps.charge < 0][dilep_evts], ak.values_astype(dilep_evts, int)), # same as second minus
        "Down_minus"  : gen_charged_leps[gen_charged_leps.charge < 0][dilep_evts] if ak.any(dilep_evts) else ak.unflatten(gen_charged_leps[gen_charged_leps.charge < 0][dilep_evts], ak.values_astype(dilep_evts, int)), # same as first minus
    })

        ## HAD
    HAD_evts = ak.zip({
        "TTbar"       : Gen_TTbar[had_evts] if ak.any(had_evts) else ak.unflatten(Gen_TTbar[had_evts], ak.values_astype(had_evts, int)),
        "Top"         : Gen_Top_Pairs[np.sign(Gen_Top_Pairs.charge) == 1][had_evts] if ak.any(had_evts) else ak.unflatten(Gen_Top_Pairs[np.sign(Gen_Top_Pairs.charge) == 1][had_evts], ak.values_astype(had_evts, int)),
        "Tbar"        : Gen_Top_Pairs[np.sign(Gen_Top_Pairs.charge) == -1][had_evts] if ak.any(had_evts) else ak.unflatten(Gen_Top_Pairs[np.sign(Gen_Top_Pairs.charge) == -1][had_evts], ak.values_astype(had_evts, int)),
        "B"           : Gen_B_Pairs[np.sign(Gen_B_Pairs.charge) == -1][had_evts] if ak.any(had_evts) else ak.unflatten(Gen_B_Pairs[np.sign(Gen_B_Pairs.charge) == -1][had_evts], ak.values_astype(had_evts, int)),
        "Bbar"        : Gen_B_Pairs[np.sign(Gen_B_Pairs.charge) == 1][had_evts] if ak.any(had_evts) else ak.unflatten(Gen_B_Pairs[np.sign(Gen_B_Pairs.charge) == 1][had_evts], ak.values_astype(had_evts, int)),
        "Wplus"       : Gen_W_Pairs[Gen_W_Pairs.charge == 1][had_evts] if ak.any(had_evts) else ak.unflatten(Gen_W_Pairs[Gen_W_Pairs.charge == 1][had_evts], ak.values_astype(had_evts, int)),
        "Wminus"      : Gen_W_Pairs[Gen_W_Pairs.charge == -1][had_evts] if ak.any(had_evts) else ak.unflatten(Gen_W_Pairs[Gen_W_Pairs.charge == -1][had_evts], ak.values_astype(had_evts, int)),
        "First_plus"  : Gen_Wparton_Pairs[had_evts][Gen_Wparton_Pairs[had_evts].charge > 0][:, 0] if ak.any(had_evts) else ak.unflatten(Gen_Wparton_Pairs[had_evts][Gen_Wparton_Pairs[had_evts].charge > 0][:, 0], ak.values_astype(had_evts, int)), # leading positively-charged parton
        "Second_plus" : Gen_Wparton_Pairs[had_evts][Gen_Wparton_Pairs[had_evts].charge > 0][:, 1] if ak.any(had_evts) else ak.unflatten(Gen_Wparton_Pairs[had_evts][Gen_Wparton_Pairs[had_evts].charge > 0][:, 1], ak.values_astype(had_evts, int)), # subleading positively-charged parton
        "First_minus" : Gen_Wparton_Pairs[had_evts][Gen_Wparton_Pairs[had_evts].charge < 0][:, 0] if ak.any(had_evts) else ak.unflatten(Gen_Wparton_Pairs[had_evts][Gen_Wparton_Pairs[had_evts].charge < 0][:, 0], ak.values_astype(had_evts, int)), # leading negatively-charged parton
        "Second_minus": Gen_Wparton_Pairs[had_evts][Gen_Wparton_Pairs[had_evts].charge < 0][:, 1] if ak.any(had_evts) else ak.unflatten(Gen_Wparton_Pairs[had_evts][Gen_Wparton_Pairs[had_evts].charge < 0][:, 1], ak.values_astype(had_evts, int)), # subleading negatively-charged parton
        "Up_plus"     : gen_wpartons_up[gen_wpartons_up.charge > 0][had_evts] if ak.any(had_evts) else ak.unflatten(gen_wpartons_up[gen_wpartons_up.charge > 0][had_evts], ak.values_astype(had_evts, int)), # positively-charged up-type parton
        "Down_plus"   : gen_wpartons_dw[gen_wpartons_dw.charge > 0][had_evts] if ak.any(had_evts) else ak.unflatten(gen_wpartons_dw[gen_wpartons_dw.charge > 0][had_evts], ak.values_astype(had_evts, int)), # positively-charged down-type parton
        "Up_minus"    : gen_wpartons_up[gen_wpartons_up.charge < 0][had_evts] if ak.any(had_evts) else ak.unflatten(gen_wpartons_up[gen_wpartons_up.charge < 0][had_evts], ak.values_astype(had_evts, int)), # negatively-charged up-type parton
        "Down_minus"  : gen_wpartons_dw[gen_wpartons_dw.charge < 0][had_evts] if ak.any(had_evts) else ak.unflatten(gen_wpartons_dw[gen_wpartons_dw.charge < 0][had_evts], ak.values_astype(had_evts, int)), # negatively-charged down-type parton
    })

        ## SEMILEP
    SEMILEP_evts = ak.zip({
        "TTbar"   : Gen_TTbar[semilep_evts] if ak.any(semilep_evts) else ak.unflatten(Gen_TTbar[semilep_evts], ak.values_astype(semilep_evts, int)),
        "Top"     : Gen_Top_Pairs[np.sign(Gen_Top_Pairs.charge) == 1][semilep_evts] if ak.any(semilep_evts) else ak.unflatten(Gen_Top_Pairs[np.sign(Gen_Top_Pairs.charge) == 1][semilep_evts], ak.values_astype(semilep_evts, int)),
        "Tbar"    : Gen_Top_Pairs[np.sign(Gen_Top_Pairs.charge) == -1][semilep_evts] if ak.any(semilep_evts) else ak.unflatten(Gen_Top_Pairs[np.sign(Gen_Top_Pairs.charge) == -1][semilep_evts], ak.values_astype(semilep_evts, int)),
        "THad"    : Gen_Top_Pairs[Gen_Top_Pairs.decaytype == 2][semilep_evts] if ak.any(semilep_evts) else ak.unflatten(Gen_Top_Pairs[Gen_Top_Pairs.decaytype == 2][semilep_evts], ak.values_astype(semilep_evts, int)), 
        "TLep"    : Gen_Top_Pairs[Gen_Top_Pairs.decaytype == 1][semilep_evts] if ak.any(semilep_evts) else ak.unflatten(Gen_Top_Pairs[Gen_Top_Pairs.decaytype == 1][semilep_evts], ak.values_astype(semilep_evts, int)),
        "BHad"    : Gen_B_Pairs[Gen_B_Pairs.decaytype == 2][semilep_evts]if ak.any(semilep_evts) else ak.unflatten(Gen_B_Pairs[Gen_B_Pairs.decaytype == 2][semilep_evts], ak.values_astype(semilep_evts, int)),
        "BLep"    : Gen_B_Pairs[Gen_B_Pairs.decaytype == 1][semilep_evts] if ak.any(semilep_evts) else ak.unflatten(Gen_B_Pairs[Gen_B_Pairs.decaytype == 1][semilep_evts], ak.values_astype(semilep_evts, int)),
        "WHad"    : Gen_W_Pairs[Gen_W_Pairs.decaytype == 2][semilep_evts] if ak.any(semilep_evts) else ak.unflatten(Gen_W_Pairs[Gen_W_Pairs.decaytype == 2][semilep_evts], ak.values_astype(semilep_evts, int)),
        "WLep"    : Gen_W_Pairs[Gen_W_Pairs.decaytype == 1][semilep_evts] if ak.any(semilep_evts) else ak.unflatten(Gen_W_Pairs[Gen_W_Pairs.decaytype == 1][semilep_evts], ak.values_astype(semilep_evts, int)),
        "Lepton"  : gen_charged_leps[semilep_evts] if ak.any(semilep_evts) else ak.unflatten(gen_charged_leps[semilep_evts], ak.values_astype(semilep_evts, int)),
        "Nu"      : gen_neutral_leps[semilep_evts] if ak.any(semilep_evts) else ak.unflatten(gen_neutral_leps[semilep_evts], ak.values_astype(semilep_evts, int)),
        "WJa"     : Gen_Wparton_Pairs[:, 0][semilep_evts] if ak.any(semilep_evts) else ak.unflatten(Gen_Wparton_Pairs[:, 0][semilep_evts], ak.values_astype(semilep_evts, int)),
        "WJb"     : Gen_Wparton_Pairs[:, 1][semilep_evts] if ak.any(semilep_evts) else ak.unflatten(Gen_Wparton_Pairs[:, 1][semilep_evts], ak.values_astype(semilep_evts, int)),
        "Up_Had"  : gen_wpartons_up[semilep_evts] if ak.any(semilep_evts) else ak.unflatten(gen_wpartons_up[semilep_evts], ak.values_astype(semilep_evts, int)),
        "Down_Had": gen_wpartons_dw[semilep_evts] if ak.any(semilep_evts) else ak.unflatten(gen_wpartons_dw[semilep_evts], ak.values_astype(semilep_evts, int)),
    })

        # make dictionary to return
    GenObjects = dict({
        "SL" : SEMILEP_evts,
        "DL" : DILEP_evts,
        "Had" : HAD_evts,
    })

    return GenObjects

def select(events, mode='NORMAL'):
    modes_to_choose = {
        'NORMAL' : select_normal,
        #'LHE' : select_lhe,
    }

    if mode not in modes_to_choose.keys():
        raise IOError("Gen Object mode %s not available" % mode)

    w_decay_momid = 6 if 'MADGRAPH' in mode else 24

    GenObjs = modes_to_choose[mode](events['GenPart'], w_decay_momid)

    events["SL"] = GenObjs["SL"]
    events["DL"] = GenObjs["DL"]
    events["Had"] = GenObjs["Had"]
    
