from pdb import set_trace
import numpy as np
from coffea.analysis_tools import Weights
import Utilities.systematics as systematics
import awkward as ak
import Utilities.make_variables as make_vars

def get_event_weights(events, year: str, corrections, isTTbar=False, isSignal=False):
    weights = Weights(len(events), storeIndividual=True) # store individual variations

        ## only apply to MC
    if not events.metadata["dataset"].startswith("data_Single"):
            ## Prefire Corrections
        if (year != "2018") and (corrections["Prefire"] == True) and ("L1PreFiringWeight" in events.fields):
            weights.add("Prefire",
                ak.copy(events["L1PreFiringWeight"]["Nom"]),
                ak.copy(events["L1PreFiringWeight"]["Up"]),
                ak.copy(events["L1PreFiringWeight"]["Dn"])
            )

            ## Generator Weights
        weights.add("genweight", ak.copy(events["genWeight"]))
    
            ## Pileup Reweighting
        if "Pileup" in corrections.keys():
                # treat interference samples differently
            if (events.metadata["dataset"].startswith("AtoTT") or events.metadata["dataset"].startswith("HtoTT")) and ("Int" in events.metadata["dataset"]):
                central_pu_wt = ak.where(events["genWeight"] > 0, corrections["Pileup"]["%s_pos" % events.metadata["dataset"]]["central"](events["Pileup"]["nTrueInt"]),\
                    corrections["Pileup"]["%s_neg" % events.metadata["dataset"]]["central"](events["Pileup"]["nTrueInt"]))
                up_pu_wt = ak.where(events["genWeight"] > 0, corrections["Pileup"]["%s_pos" % events.metadata["dataset"]]["up"](events["Pileup"]["nTrueInt"]),\
                    corrections["Pileup"]["%s_neg" % events.metadata["dataset"]]["up"](events["Pileup"]["nTrueInt"]))
                down_pu_wt = ak.where(events["genWeight"] > 0, corrections["Pileup"]["%s_pos" % events.metadata["dataset"]]["down"](events["Pileup"]["nTrueInt"]),\
                    corrections["Pileup"]["%s_neg" % events.metadata["dataset"]]["down"](events["Pileup"]["nTrueInt"]))

                weights.add("Pileup",
                    ak.copy(central_pu_wt),
                    ak.copy(up_pu_wt),
                    ak.copy(down_pu_wt)
                )
            else:
                weights.add("Pileup",
                    ak.copy(corrections["Pileup"][events.metadata["dataset"]]["central"](events["Pileup"]["nTrueInt"])),
                    ak.copy(corrections["Pileup"][events.metadata["dataset"]]["up"](events["Pileup"]["nTrueInt"])),
                    ak.copy(corrections["Pileup"][events.metadata["dataset"]]["down"](events["Pileup"]["nTrueInt"]))
                )

        ## PS and LHE weights for ttbar events
        if isTTbar:
            ## PS Weight variations
            # PSWeight definitions can be found here: https://github.com/cms-sw/cmssw/blob/CMSSW_10_6_X/PhysicsTools/NanoAOD/plugins/GenWeightsTableProducer.cc#L543-L546
            psweights = events["PSWeight"]
            weights.add("ISR",
                np.ones(len(events)),
                ak.copy(psweights[:, 0]), # (ISR=2, FSR=1)
                ak.copy(psweights[:, 2]), # (ISR=0.5, FSR=1)
            )
            weights.add("FSR",
                np.ones(len(events)),
                ak.copy(psweights[:, 1]), # (ISR=1, FSR=2)
                ak.copy(psweights[:, 3]), # (ISR=1, FSR=0.5)
            )

            ## LHEScale Weight Variations
            # LHEScaleWeight definitions can be found here: https://cms-nanoaod-integration.web.cern.ch/integration/master/mc94X_doc.html#LHE
            lheweights = events["LHEScaleWeight"]
            weights.add("FACTOR",
                np.ones(len(events)),
                ak.copy(lheweights[:, 5]), # (muR=1, muF=2)
                ak.copy(lheweights[:, 3]), # (muR=1, muF=0.5)
            )
            weights.add("RENORM",
                np.ones(len(events)),
                ak.copy(lheweights[:, 7]), # (muR=2, muF=1)
                ak.copy(lheweights[:, 1]), # (muR=0.5, muF=1)
            )
            weights.add("RENORM_FACTOR_SAME",
                np.ones(len(events)),
                ak.copy(lheweights[:, 8]), # (muR=2, muF=2), RENORM_UP_FACTOR_UP
                ak.copy(lheweights[:, 0]), # (muR=0.5, muF=0.5), RENORM_DW_FACTOR_DW
            )
            weights.add("RENORM_FACTOR_DIFF",
                np.ones(len(events)),
                ak.copy(lheweights[:, 6]), # (muR=2, muF=0.5), RENORM_UP_FACTOR_DW
                ak.copy(lheweights[:, 2]), # (muR=0.5, muF=2), RENORM_DW_FACTOR_UP
            )

        ## PS and LHE weights for signal events
        if isSignal:
            ## PS Weight variations
            ## PSWeight definitions can be found here: https://github.com/cms-sw/cmssw/blob/CMSSW_10_6_X/PhysicsTools/NanoAOD/plugins/GenWeightsTableProducer.cc#L543-L546
            psweights = events["PSWeight"]
            weights.add("AH_ISR",
                np.ones(len(events)),
                ak.copy(psweights[:, 0]), # (ISR=2, FSR=1)
                ak.copy(psweights[:, 2]), # (ISR=0.5, FSR=1)
            )
            weights.add("AH_FSR",
                np.ones(len(events)),
                ak.copy(psweights[:, 1]), # (ISR=1, FSR=2)
                ak.copy(psweights[:, 3]), # (ISR=1, FSR=0.5)
            )

            ## LHEScale Weight Variations
            # LHEScaleWeight definitions can be found here: https://cms-nanoaod-integration.web.cern.ch/integration/master/mc94X_doc.html#LHE
            lheweights = events["LHEScaleWeight"]
            weights.add("AH_FACTOR",
                np.ones(len(events)),
                ak.copy(lheweights[:, 4]) if "W9p0" in events.metadata["dataset"] else ak.copy(lheweights[:, 5]), # (muR=1, muF=2)
                ak.copy(lheweights[:, 3]), # (muR=1, muF=0.5)
            )
            weights.add("AH_RENORM",
                np.ones(len(events)),
                ak.copy(lheweights[:, 6]) if "W9p0" in events.metadata["dataset"] else ak.copy(lheweights[:, 7]), # (muR=2, muF=1)
                ak.copy(lheweights[:, 1]), # (muR=0.5, muF=1)
            )
            weights.add("AH_RENORM_FACTOR_SAME",
                np.ones(len(events)),
                ak.copy(lheweights[:, 7]) if "W9p0" in events.metadata["dataset"] else ak.copy(lheweights[:, 8]), # (muR=2, muF=2), RENORM_UP_FACTOR_UP
                ak.copy(lheweights[:, 0]), # (muR=0.5, muF=0.5), RENORM_DW_FACTOR_DW
            )
            weights.add("AH_RENORM_FACTOR_DIFF",
                np.ones(len(events)),
                ak.copy(lheweights[:, 5]) if "W9p0" in events.metadata["dataset"] else ak.copy(lheweights[:, 6]), # (muR=2, muF=0.5), RENORM_UP_FACTOR_DW
                ak.copy(lheweights[:, 2]), # (muR=0.5, muF=2), RENORM_DW_FACTOR_UP
            )

    return weights    


def compute_btagSF_weights(constructor, jets, wp, mask_3j, mask_4pj, sysnames, evt_weights):
    btag_weights = {key : np.ones(len(jets)) for key in ["central"]+sysnames}

    btagger_3j_wts = constructor["3Jets"].get_scale_factor(jets=jets[mask_3j], passing_cut=wp)
    btagger_4pj_wts = constructor["4PJets"].get_scale_factor(jets=jets[mask_4pj], passing_cut=wp)

        # fll dict of btag weights
    for wt_name in btag_weights.keys():
        btag_weights[wt_name][mask_3j]  = ak.prod(btagger_3j_wts[wt_name], axis=1)
        btag_weights[wt_name][mask_4pj] = ak.prod(btagger_4pj_wts[wt_name], axis=1)

    #set_trace()
    btag_names = sorted(set([name.replace("_up", "") for name in sysnames if name.endswith("_up")] + [name.replace("_down", "") for name in sysnames if name.endswith("_down")]))

    evt_weights.add_multivariation(
        name="btag",
        weight=np.copy(btag_weights["central"]),
        modifierNames=btag_names,
        weightsUp=[np.copy(btag_weights[f"{name}_up"]) for name in btag_names],
        weightsDown=[np.copy(btag_weights[f"{name}_down"]) for name in btag_names],
    )

    return evt_weights
    


def get_pdf_weights(df):
    if len(sorted(set(df["nLHEPdfWeight"]))) > 1:
        raise ValueError("Differing number of PDF weights for events")
    pdfweights = ak.JaggedArray.fromcounts(df["nLHEPdfWeight"], df["LHEPdfWeight"])
    df["PDFWeights"] = pdfweights


def SF_pt(pt):
    SF = 0.103 * np.exp(-0.0118*pt) - 0.000134*pt + 0.973
    return SF

def get_TopPt_weights(events):
    # values and application of top pT reweighting taken from https://twiki.cern.ch/twiki/bin/view/CMS/TopPtReweighting#How_to_practically_apply_default
    #set_trace()
    genparts = events["GenPart"]
    gen_tops = genparts[(genparts.hasFlags(["isLastCopy"])) & (genparts.pdgId == 6)]
    gen_tbars = genparts[(genparts.hasFlags(["isLastCopy"])) & (genparts.pdgId == -6)]
    SFtop = SF_pt(gen_tops.pt)
    SFtbar = SF_pt(gen_tbars.pt)
    SFtot = np.sqrt(SFtop * SFtbar)
    return SFtot


def get_nnlo_weights(correction, events):
    #set_trace()
    var = correction["Var"]
    dist = correction["Correction"]
    if var != "mtt_vs_top_ctstar":
        raise ValueError(f"{var} not supported for NNLO kinematic reweighting")

    if "SL" in events.fields:
        wts = np.ones(len(events))

            # SL events
        sl_evts = ak.num(events["SL"]) > 0
        sl_top_ctstar, sl_tbar_ctstar = make_vars.ctstar(events["SL"][sl_evts]["Top"], events["SL"][sl_evts]["Tbar"])
        sl_top_ctstar = ak.flatten(sl_top_ctstar, axis=None)
        wts[sl_evts] = dist(ak.flatten(events["SL"][sl_evts]["TTbar"].mass, axis=None), sl_top_ctstar)

            # DL events
        dl_evts = ak.num(events["DL"]) > 0
        dl_top_ctstar, dl_tbar_ctstar = make_vars.ctstar(events["DL"][dl_evts]["Top"], events["DL"][dl_evts]["Tbar"])
        dl_top_ctstar = ak.flatten(dl_top_ctstar, axis=None)
        wts[dl_evts] = dist(ak.flatten(events["DL"][dl_evts]["TTbar"].mass, axis=None), dl_top_ctstar)

            # Had events
        had_evts = ak.num(events["Had"]) > 0
        had_top_ctstar, had_tbar_ctstar = make_vars.ctstar(events["Had"][had_evts]["Top"], events["Had"][had_evts]["Tbar"])
        had_top_ctstar = ak.flatten(had_top_ctstar, axis=None)
        wts[had_evts] = dist(ak.flatten(events["Had"][had_evts]["TTbar"].mass, axis=None), had_top_ctstar)
    else:
        genparts = events["GenPart"]
        gen_tops = genparts[(genparts.hasFlags(["isLastCopy"])) & (genparts.pdgId == 6)]
        gen_tbars = genparts[(genparts.hasFlags(["isLastCopy"])) & (genparts.pdgId == -6)]
        top_ctstar, tbar_ctstar = make_vars.ctstar(gen_tops, gen_tbars)
        mtt = (gen_tops+gen_tbars).mass
        wts = dist(ak.flatten(mtt, axis=None), ak.flatten(top_ctstar, axis=None))

    ##if "thad_pt" in var:
    #if "mtt_vs_top_ctstar" in var:
    #    wts = dist(ak.flatten(mtt, axis=None), ak.flatten(top_ctstar, axis=None))
    #else:
    #    raise ValueError(f"{var} not supported for NNLO kinematic reweighting")

    return wts

    

def get_ewk_weights(correction, events):
    event_status = -1*np.ones(len(events)) # initialize array to categorize events as gg (0), qqu (1), or qqd (2)
    
    lheparts = events["LHEPart"]
    incoming_lhe = lheparts[lheparts.status == -1]
    outgoing_lhe = lheparts[lheparts.status == 1]
    
    # find easiest event classifications
    gg_events = ak.all(incoming_lhe.pdgId == 21, axis=1) # events whose initial state partons are both gluons
    qqu_events = ak.all(abs(incoming_lhe.pdgId) < 7, axis=1) & ak.all(np.mod(incoming_lhe.pdgId, 2) == 0, axis=1) # events whose intial state partons are both up-type quarks
    qqd_events = ak.all(abs(incoming_lhe.pdgId) < 7, axis=1) & ak.all(np.mod(incoming_lhe.pdgId, 2) == 1, axis=1) # events whose intial state partons are both down-type quarks
    
    event_status[gg_events] = 0
    event_status[qqu_events] = 1
    event_status[qqd_events] = 2
    
    # find classification for qg initial state events
    qg_events = (ak.sum(incoming_lhe.pdgId == 21, axis=1) == 1) & (ak.sum(abs(incoming_lhe.pdgId) < 7, axis=1) == 1) # events that have one gluon and one quark as inital state partons
    qg_event_status = -1*np.ones(ak.sum(qg_events))
    
    qg_quarks = incoming_lhe[qg_events][incoming_lhe[qg_events].pdgId != 21] # quarks in qg initial states
    qg_gluons = incoming_lhe[qg_events][incoming_lhe[qg_events].pdgId == 21] # gluons in qg initial states
    
        # set lorentz vectors to correct values to find CM frame
    qg_quarks_lv = ak.Array({
        "x" : qg_quarks.px, "y" : qg_quarks.py, "z" : qg_quarks.incomingpz, "t" : abs(qg_quarks.incomingpz)
    }, with_name="LorentzVector")
    qg_gluons_lv = ak.Array({
        "x" : qg_gluons.px, "y" : qg_gluons.py, "z" : qg_gluons.incomingpz, "t" : abs(qg_gluons.incomingpz)
    }, with_name="LorentzVector")
    
        # find a boost incoming lhe partons CM system
    qg_system = (qg_quarks_lv+qg_gluons_lv)
    qg_system_boost = qg_system.boostvec
    
        # find additional outgoing parton from qg events and boost into CM system
    extra_out_quark = outgoing_lhe[:, 0][qg_events]
    rf_extra_out_quark = extra_out_quark.boost(qg_system_boost*-1)
    
    extra_out_quark_up = np.mod(abs(outgoing_lhe[qg_events].pdgId[:, 0]), 2) == 0
    extra_out_quark_dw = np.mod(abs(outgoing_lhe[qg_events].pdgId[:, 0]), 2) == 1
    rf_extra_out_quark_up = rf_extra_out_quark[extra_out_quark_up]
    rf_extra_out_quark_dw = rf_extra_out_quark[extra_out_quark_dw]
        # if sign(pz boosted outgoing up quark) == sign(pz incoming gluon) event is qqu else gg
    qg_event_status[extra_out_quark_up] = ak.to_numpy(ak.flatten(ak.where(np.sign(rf_extra_out_quark_up.pz) == np.sign(qg_gluons.incomingpz[extra_out_quark_up]), 1, 0)))
        # if sign(pz boosted outgoing down quark) == sign(pz incoming gluon) event is qqd else gg
    qg_event_status[extra_out_quark_dw] = ak.to_numpy(ak.flatten(ak.where(np.sign(rf_extra_out_quark_dw.pz) == np.sign(qg_gluons.incomingpz[extra_out_quark_dw]), 2, 0)))
    
        # set event status for qg events
    event_status[qg_events] = qg_event_status
    
    if "SL" in events.fields:
            # calculate cpTP
                # SL events
        sl_evts = ak.num(events["SL"]) > 0
        sl_gentops, sl_gentbars = events["SL"]["Top"], events["SL"]["Tbar"]
        sl_top_cpTP, sl_tbar_cpTP = make_vars.cpTP(sl_gentops, sl_gentbars)
        sl_top_cpTP, sl_tbar_cpTP = ak.flatten(sl_top_cpTP, axis=None), ak.flatten(sl_tbar_cpTP, axis=None)
                # DL events
        dl_evts = ak.num(events["DL"]) > 0
        dl_gentops, dl_gentbars = events["DL"]["Top"], events["DL"]["Tbar"]
        dl_top_cpTP, dl_tbar_cpTP = make_vars.cpTP(dl_gentops, dl_gentbars)
        dl_top_cpTP, dl_tbar_cpTP = ak.flatten(dl_top_cpTP, axis=None), ak.flatten(dl_tbar_cpTP, axis=None)
                # Had events
        had_evts = ak.num(events["Had"]) > 0
        had_gentops, had_gentbars = events["Had"]["Top"], events["Had"]["Tbar"]
        had_top_cpTP, had_tbar_cpTP = make_vars.cpTP(had_gentops, had_gentbars)
        had_top_cpTP, had_tbar_cpTP = ak.flatten(had_top_cpTP, axis=None), ak.flatten(had_tbar_cpTP, axis=None)

            # set values for mtt and top cpTP
        gen_mtt, gen_top_cpTP = -10*np.ones(len(events)), -10*np.ones(len(events))
        gen_mtt[sl_evts] = ak.flatten(events["SL"]["TTbar"].mass, axis=None)
        gen_mtt[dl_evts] = ak.flatten(events["DL"]["TTbar"].mass, axis=None)
        gen_mtt[had_evts] = ak.flatten(events["Had"]["TTbar"].mass, axis=None)
        gen_top_cpTP[sl_evts], gen_top_cpTP[dl_evts], gen_top_cpTP[had_evts] = sl_top_cpTP, dl_top_cpTP, had_top_cpTP

    else:
        genparts = events["GenPart"]
        gen_tops = genparts[(genparts.hasFlags(["isLastCopy"])) & (genparts.pdgId == 6)]
        gen_tbars = genparts[(genparts.hasFlags(["isLastCopy"])) & (genparts.pdgId == -6)]
        gen_top_cpTP, gen_tbar_cpTP = make_vars.cpTP(gen_tops, gen_tbars)
        gen_top_cpTP = ak.flatten(gen_top_cpTP, axis=None)
        gen_mtt = ak.flatten((gen_tops+gen_tbars).mass, axis=None)

    ewk_wts_dict = {val : np.ones(len(events)) for val in ["yt_0", "yt_0.5", "yt_1", "yt_1.5", "yt_2"]}
    for rewt_type, wts in ewk_wts_dict.items():
        wts[np.where(event_status == 0)] = correction[f"ratio_gg_{rewt_type}"](gen_mtt[np.where(event_status == 0)], gen_top_cpTP[np.where(event_status == 0)])
        wts[np.where(event_status == 1)] = correction[f"ratio_qqu_{rewt_type}"](gen_mtt[np.where(event_status == 1)], gen_top_cpTP[np.where(event_status == 1)])
        wts[np.where(event_status == 2)] = correction[f"ratio_qqd_{rewt_type}"](gen_mtt[np.where(event_status == 2)], gen_top_cpTP[np.where(event_status == 2)])

    return ewk_wts_dict

    
def get_Otto_ewk_weights(correction, events):
    if "SL" in events.fields:
            # calculate ctstar
                # SL events
        sl_evts = ak.num(events["SL"]) > 0
        sl_gentops, sl_gentbars = events["SL"]["Top"], events["SL"]["Tbar"]
        sl_top_ctstar, sl_tbar_ctstar = make_vars.ctstar(sl_gentops, sl_gentbars)
        sl_top_ctstar, sl_tbar_ctstar = ak.flatten(sl_top_ctstar, axis=None), ak.flatten(sl_tbar_ctstar, axis=None)
                # DL events
        dl_evts = ak.num(events["DL"]) > 0
        dl_gentops, dl_gentbars = events["DL"]["Top"], events["DL"]["Tbar"]
        dl_top_ctstar, dl_tbar_ctstar = make_vars.ctstar(dl_gentops, dl_gentbars)
        dl_top_ctstar, dl_tbar_ctstar = ak.flatten(dl_top_ctstar, axis=None), ak.flatten(dl_tbar_ctstar, axis=None)
                # Had events
        had_evts = ak.num(events["Had"]) > 0
        had_gentops, had_gentbars = events["Had"]["Top"], events["Had"]["Tbar"]
        had_top_ctstar, had_tbar_ctstar = make_vars.ctstar(had_gentops, had_gentbars)
        had_top_ctstar, had_tbar_ctstar = ak.flatten(had_top_ctstar, axis=None), ak.flatten(had_tbar_ctstar, axis=None)

            # set values for mtt and top ctstar
        gen_mtt, gen_top_ctstar = -10*np.ones(len(events)), -10*np.ones(len(events))
        gen_mtt[sl_evts] = ak.flatten(events["SL"]["TTbar"].mass, axis=None)
        gen_mtt[dl_evts] = ak.flatten(events["DL"]["TTbar"].mass, axis=None)
        gen_mtt[had_evts] = ak.flatten(events["Had"]["TTbar"].mass, axis=None)
        gen_top_ctstar[sl_evts], gen_top_ctstar[dl_evts], gen_top_ctstar[had_evts] = sl_top_ctstar, dl_top_ctstar, had_top_ctstar

    else:
        genparts = events["GenPart"]
        gen_tops = genparts[(genparts.hasFlags(["isLastCopy"])) & (genparts.pdgId == 6)]
        gen_tbars = genparts[(genparts.hasFlags(["isLastCopy"])) & (genparts.pdgId == -6)]
        gen_top_ctstar, gen_tbar_ctstar = make_vars.ctstar(gen_tops, gen_tbars)
        gen_top_ctstar = ak.flatten(gen_top_ctstar, axis=None)
        gen_mtt = ak.flatten((gen_tops+gen_tbars).mass, axis=None)

    #set_trace()
        # init dict of weights for different yukawa ratio values
    ewk_wts_dict = {}

        # find weight for delta QCD 
    deltaQCD_wt = correction["DeltaQCD"](gen_mtt, gen_top_ctstar)
    ewk_wts_dict["DeltaQCD"] = ak.copy(deltaQCD_wt)

    rewt_vals = ["0.0", "0.5", "0.88", "0.9", "1.0", "1.1", "1.11", "1.5", "2.0", "3.0", "4.0"]
    #rewt_vals = ["0.0", "0.5", "0.9", "1.0", "1.1", "1.5", "2.0", "3.0", "4.0"]
    for rewt_type in rewt_vals:
            # NLO EW kfactor weights
        ewk_wts_dict[f"KFactor_{rewt_type}"] = ak.copy(correction[f"EW_KFactor_{rewt_type}"](gen_mtt, gen_top_ctstar))
        ewk_wts_dict[f"Rebinned_KFactor_{rewt_type}"] = ak.copy(correction[f"Rebinned_EW_KFactor_{rewt_type}"](gen_mtt, gen_top_ctstar)) # coarser binning
            # delta EW weights
        ewk_wts_dict[f"DeltaEW_{rewt_type}"] = ak.copy(correction[f"DeltaEW_{rewt_type}"](gen_mtt, gen_top_ctstar))
        ewk_wts_dict[f"Rebinned_DeltaEW_{rewt_type}"] = ak.copy(correction[f"Rebinned_DeltaEW_{rewt_type}"](gen_mtt, gen_top_ctstar)) # coarser binning


    return ewk_wts_dict
    


def get_comb_lepSF(year: str, lepton: str, corrections, pt: np.ndarray, eta: np.ndarray, shift="Central"):
    if not (shift == "Central" or shift == "Error"):
        raise ValueError("Shift value %s not defined" % shift)
    #set_trace()
    sf_dict = corrections[lepton]
    eta_ranges = sf_dict["eta_ranges"]
    lepSFs = np.ones(pt.size)

    for idx, eta_range in enumerate(eta_ranges):
        mask = (eta >= eta_range[0]) & (eta < eta_range[1]) # find inds that are within given eta range
        if not mask.any(): continue # no values fall within eta range
        #set_trace()
        recoSF = sf_dict["Reco_ID"]["Central"]["eta_bin%i" % idx]
        trigSF = sf_dict["Trig"]["Central"]["eta_bin%i" % idx]
        lepSFs[mask] = recoSF(pt[mask])*trigSF(pt[mask])

    return lepSFs


def get_lepton_sf(sf_dict, pt, eta, tight_lep_mask, leptype):
    mu_schema = {
        "central" : ["ID_Central", "ISO_Central", "TRIG_Central", "RECO_Central"],
            # variations of ID
        "IDtotUp" : ["ID_Error_totUp", "ISO_Central", "TRIG_Central", "RECO_Central"],
        "IDtotDown" : ["ID_Error_totDown", "ISO_Central", "TRIG_Central", "RECO_Central"],
        "IDstatUp" : ["ID_Error_statUp", "ISO_Central", "TRIG_Central", "RECO_Central"],
        "IDstatDown" : ["ID_Error_statDown", "ISO_Central", "TRIG_Central", "RECO_Central"],
        "IDsystUp" : ["ID_Error_systUp", "ISO_Central", "TRIG_Central", "RECO_Central"],
        "IDsystDown" : ["ID_Error_systDown", "ISO_Central", "TRIG_Central", "RECO_Central"],
            # variations of ISO
        "ISOtotUp" : ["ID_Central", "ISO_Error_totUp", "TRIG_Central", "RECO_Central"],
        "ISOtotDown" : ["ID_Central", "ISO_Error_totDown", "TRIG_Central", "RECO_Central"],
        "ISOstatUp" : ["ID_Central", "ISO_Error_statUp", "TRIG_Central", "RECO_Central"],
        "ISOstatDown" : ["ID_Central", "ISO_Error_statDown", "TRIG_Central", "RECO_Central"],
        "ISOsystUp" : ["ID_Central", "ISO_Error_systUp", "TRIG_Central", "RECO_Central"],
        "ISOsystDown" : ["ID_Central", "ISO_Error_systDown", "TRIG_Central", "RECO_Central"],
            # variations of TRIG
        "TRIGtotUp" : ["ID_Central", "ISO_Central", "TRIG_Error_totUp", "RECO_Central"],
        "TRIGtotDown" : ["ID_Central", "ISO_Central", "TRIG_Error_totDown", "RECO_Central"],
        "TRIGstatUp" : ["ID_Central", "ISO_Central", "TRIG_Error_statUp", "RECO_Central"],
        "TRIGstatDown" : ["ID_Central", "ISO_Central", "TRIG_Error_statDown", "RECO_Central"],
        "TRIGsystUp" : ["ID_Central", "ISO_Central", "TRIG_Error_systUp", "RECO_Central"],
        "TRIGsystDown" : ["ID_Central", "ISO_Central", "TRIG_Error_systDown", "RECO_Central"],
            # variations of RECO
        "RECOtotUp" : ["ID_Central", "ISO_Central", "TRIG_Central", "RECO_Error_totUp"],
        "RECOtotDown" : ["ID_Central", "ISO_Central", "TRIG_Central", "RECO_Error_totDown"],
    }
    el_schema = {
        "central" : ["ID_Central", "TRIG_Central", "RECO_Central"],
            # variations of ID
        "IDtotUp" : ["ID_Error_totUp", "TRIG_Central", "RECO_Central"],
        "IDtotDown" : ["ID_Error_totDown", "TRIG_Central", "RECO_Central"],
            # variations of TRIG
        "TRIGtotUp" : ["ID_Central", "TRIG_Error_totUp", "RECO_Central"],
        "TRIGtotDown" : ["ID_Central", "TRIG_Error_totDown", "RECO_Central"],
            # variations of RECO
        "RECOtotUp" : ["ID_Central", "TRIG_Central", "RECO_Error_totUp"],
        "RECOtotDown" : ["ID_Central", "TRIG_Central", "RECO_Error_totDown"],
    }
    schema_to_use = mu_schema if leptype == "Muons" else el_schema

    #set_trace()
    indiv_SF_sources = {}
    for sf_type in sf_dict.keys():
        isAbsEta = sf_dict[sf_type]["isAbsEta"]
        if "eta_ranges" in sf_dict[sf_type].keys():
            #if leptype == "Electrons" : set_trace()
            eta_ranges = sf_dict[sf_type]["eta_ranges"]
            for sf_var in sf_dict[sf_type].keys():
                if ((sf_var == "eta_ranges") or (sf_var == "isAbsEta")): continue
                #print(leptype, sf_type, sf_var)
                if sf_var == "Central":
                    cen_wts = np.ones(len(pt))
                elif "Error" in sf_var:
                    errup_wts, errdw_wts = np.ones(len(pt)), np.ones(len(pt))
                else:
                    set_trace()

                for idx, eta_range in enumerate(eta_ranges):
                    mask = (np.abs(eta) >= eta_range[0]) & (np.abs(eta) < eta_range[1]) if isAbsEta else (eta >= eta_range[0]) & (eta < eta_range[1]) # find inds that are within given eta range
                    if not ak.any(mask): continue # no values fall within eta range
                    if sf_var == "Central":
                        cen_wts[mask] = sf_dict[sf_type][sf_var][f"eta_bin{idx}"](pt[mask])
                    else:
                        errup_wts[mask] = sf_dict[sf_type]["Central"][f"eta_bin{idx}"](pt[mask]) + sf_dict[sf_type][sf_var][f"eta_bin{idx}"](pt[mask])
                        errdw_wts[mask] = sf_dict[sf_type]["Central"][f"eta_bin{idx}"](pt[mask]) - sf_dict[sf_type][sf_var][f"eta_bin{idx}"](pt[mask])
                if sf_var == "Central":
                    indiv_SF_sources[f"{sf_type}_{sf_var}"] = cen_wts
                else:
                    indiv_SF_sources[f"{sf_type}_{sf_var}Up"] = errup_wts
                    indiv_SF_sources[f"{sf_type}_{sf_var}Down"] = errdw_wts
        else:
            for sf_var in sf_dict[sf_type].keys():
                if sf_var == "isAbsEta": continue
                #print(leptype, sf_type, sf_var)
                if sf_var == "Central":
                    if sf_dict[sf_type][sf_var]._dimension == 1:
                        indiv_SF_sources[f"{sf_type}_{sf_var}"] = sf_dict[sf_type][sf_var](np.abs(eta)) if isAbsEta else sf_dict[sf_type][sf_var](eta)
                    elif sf_dict[sf_type][sf_var]._dimension == 2:
                        indiv_SF_sources[f"{sf_type}_{sf_var}"] = sf_dict[sf_type][sf_var](np.abs(eta), pt) if isAbsEta else sf_dict[sf_type][sf_var](eta, pt)
                    else:
                        raise ValueError("Only 1D or 2D scale factors are supported!")
                else:
                    if sf_dict[sf_type][sf_var]._dimension == 1:
                        indiv_SF_sources[f"{sf_type}_{sf_var}Up"]   = sf_dict[sf_type]["Central"](np.abs(eta)) + sf_dict[sf_type][sf_var](np.abs(eta)) if isAbsEta else sf_dict[sf_type]["Central"](eta) + sf_dict[sf_type][sf_var](eta)
                        indiv_SF_sources[f"{sf_type}_{sf_var}Down"] = sf_dict[sf_type]["Central"](np.abs(eta)) - sf_dict[sf_type][sf_var](np.abs(eta)) if isAbsEta else sf_dict[sf_type]["Central"](eta) - sf_dict[sf_type][sf_var](eta)
                    elif sf_dict[sf_type][sf_var]._dimension == 2:
                        indiv_SF_sources[f"{sf_type}_{sf_var}Up"]   = sf_dict[sf_type]["Central"](np.abs(eta), pt) + sf_dict[sf_type][sf_var](np.abs(eta), pt) if isAbsEta else sf_dict[sf_type]["Central"](eta, pt) + sf_dict[sf_type][sf_var](eta, pt)
                        indiv_SF_sources[f"{sf_type}_{sf_var}Down"] = sf_dict[sf_type]["Central"](np.abs(eta), pt) - sf_dict[sf_type][sf_var](np.abs(eta), pt) if isAbsEta else sf_dict[sf_type]["Central"](eta, pt) - sf_dict[sf_type][sf_var](eta, pt)
                    else:
                        raise ValueError("Only 1D or 2D scale factors are supported!")

    output_SFs = {}
    for key, sf_list in schema_to_use.items():
            # initialize event weights to 1
        evt_wts = np.ones(len(tight_lep_mask))
            # make list of event weights from inidividual sources in sf_list
        arrays_list = [ak.to_numpy(indiv_SF_sources[sf_source]) for sf_source in sf_list]
            # set events that pass tight lepton mask equal to product of all event weights from inidividual sources in sf_list
        evt_wts[tight_lep_mask] = np.prod(np.vstack(arrays_list), axis=0)
        output_SFs[key] = np.copy(evt_wts)

    return output_SFs
