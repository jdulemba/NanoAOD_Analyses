from pdb import set_trace
import numpy as np
from coffea.analysis_tools import Weights
import Utilities.systematics as systematics
import awkward as ak
import Utilities.make_variables as make_vars

def get_event_weights(events, year: str, corrections, isTTbar=False):
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

    return weights    


def get_pdf_weights(df):
    if len(sorted(set(df["nLHEPdfWeight"]))) > 1:
        raise ValueError("Differing number of PDF weights for events")
    pdfweights = ak.JaggedArray.fromcounts(df["nLHEPdfWeight"], df["LHEPdfWeight"])
    df["PDFWeights"] = pdfweights
    

def get_nnlo_weights(correction, events):
    sl_evts = ak.num(events["SL"]) > 0
    dl_evts = ak.num(events["DL"]) > 0
    had_evts = ak.num(events["Had"]) > 0

        # choose seed for random integers in similar fashio to https://github.com/CoffeaTeam/coffea/blob/master/coffea/jetmet_tools/CorrectedJetsFactory.py#L246-L254
    seeds = np.array(ak.flatten(ak.concatenate([events["SL"][sl_evts]["TTbar"].pt, events["DL"][dl_evts]["TTbar"].pt, events["Had"][had_evts]["TTbar"].pt]), axis=None))[
        [0, -1]
    ].view("i4")
    randomstate = np.random.Generator(np.random.PCG64(seeds))
    which_top_to_use = randomstate.integers(low=0, high=2, size=len(events)) # top is 0, tbar is 1 

    var = correction["Var"]
    dist = correction["Correction"]

    wts = np.ones(len(events))
    if "thad_pt" in var:
            # set wts for semilep evts
        wts[sl_evts] = dist(ak.flatten(events["SL"][sl_evts]["THad"].pt, axis=None))
            # set wts for dilep evts
        dl_pt = np.where(which_top_to_use[dl_evts], ak.flatten(events["DL"][dl_evts]["Top"].pt, axis=None), ak.flatten(events["DL"][dl_evts]["Tbar"].pt, axis=None))
        wts[dl_evts] = dist(dl_pt)
            # set wts for had evts
        had_pt = np.where(which_top_to_use[had_evts], ak.flatten(events["Had"][had_evts]["Top"].pt, axis=None), ak.flatten(events["Had"][had_evts]["Tbar"].pt, axis=None))
        wts[had_evts] = dist(had_pt)

    elif "mtt_vs_thad_ctstar" in var:
            # set wts for semilep evts
        thad_ctstar, tlep_ctstar = make_vars.ctstar(events["SL"][sl_evts]["THad"], events["SL"][sl_evts]["TLep"])
        thad_ctstar, tlep_ctstar = ak.flatten(thad_ctstar, axis=None), ak.flatten(tlep_ctstar, axis=None)
        wts[sl_evts] = dist(ak.flatten(events["SL"][sl_evts]["TTbar"].mass, axis=None), thad_ctstar)
            # set wts for dilep evts
        dl_top_ctstar, dl_tbar_ctstar = make_vars.ctstar(events["DL"][dl_evts]["Top"], events["DL"][dl_evts]["Tbar"])
        dl_top_ctstar, dl_tbar_ctstar = ak.flatten(dl_top_ctstar, axis=None), ak.flatten(dl_tbar_ctstar, axis=None)
        dl_ctstar = np.where(which_top_to_use[dl_evts], dl_top_ctstar, dl_tbar_ctstar)
        wts[dl_evts] = dist(ak.flatten(events["DL"][dl_evts]["TTbar"].mass, axis=None), dl_ctstar)
            # set wts for had evts
        had_top_ctstar, had_tbar_ctstar = make_vars.ctstar(events["Had"][had_evts]["Top"], events["Had"][had_evts]["Tbar"])
        had_top_ctstar, had_tbar_ctstar = ak.flatten(had_top_ctstar, axis=None), ak.flatten(had_tbar_ctstar, axis=None)
        had_ctstar = np.where(which_top_to_use[had_evts], had_top_ctstar, had_tbar_ctstar)
        wts[had_evts] = dist(ak.flatten(events["Had"][had_evts]["TTbar"].mass, axis=None), had_ctstar)

    else:
        raise ValueError("%s not supported for NNLO kinematic reweighting" % var)

    return wts

    


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


def get_lepton_sf(year: str, lepton: str, corrections, pt, eta):
    sf_dict = corrections[lepton]
    eta_ranges = sf_dict["eta_ranges"]
    reco_cen = np.ones(len(pt))
    reco_err = np.ones(len(pt))
    trig_cen = np.ones(len(pt))
    trig_err = np.ones(len(pt))

    for idx, eta_range in enumerate(eta_ranges):
        mask = (eta >= eta_range[0]) & (eta < eta_range[1]) # find inds that are within given eta range
        if not ak.any(mask): continue # no values fall within eta range
        recoSF_cen = sf_dict["Reco_ID"]["Central"]["eta_bin%i" % idx]
        recoSF_err = sf_dict["Reco_ID"]["Error"]["eta_bin%i" % idx]
        trigSF_cen = sf_dict["Trig"]["Central"]["eta_bin%i" % idx]
        trigSF_err = sf_dict["Trig"]["Error"]["eta_bin%i" % idx]
        reco_cen[mask] = recoSF_cen(pt[mask])
        reco_err[mask] = recoSF_err(pt[mask])
        trig_cen[mask] = trigSF_cen(pt[mask])
        trig_err[mask] = trigSF_err(pt[mask])

    output_SFs = {
        "RECO_CEN" : reco_cen,
        "RECO_ERR" : reco_err,
        "TRIG_CEN" : trig_cen,
        "TRIG_ERR" : trig_err,
    }

    return output_SFs
