from pdb import set_trace
import numpy as np
from coffea.analysis_tools import Weights
import Utilities.systematics as systematics
import awkward as ak
import Utilities.make_variables as make_vars

def get_event_weights(events, year: str, corrections, BTagSFs=[], isTTbar=False):
    weights = Weights(len(events), storeIndividual=True) # store individual variations

        ## only apply to MC
    if not events.metadata['dataset'].startswith('data_Single'):
            ## Prefire Corrections
        if (year != '2018') and (corrections['Prefire'] == True):
            weights.add('Prefire',
                events['L1PreFiringWeight']['Nom'],
                events['L1PreFiringWeight']['Up'],
                events['L1PreFiringWeight']['Dn']
            )

            ## Generator Weights
        weights.add('genweight', events['genWeight'])
    
            ## Pileup Reweighting
        if 'Pileup' in corrections.keys():
            weights.add('Pileup',
                corrections['Pileup'][events.metadata['dataset']]['central'](events['Pileup']['nTrueInt']),
                corrections['Pileup'][events.metadata['dataset']]['up'](events['Pileup']['nTrueInt']),
                corrections['Pileup'][events.metadata['dataset']]['down'](events['Pileup']['nTrueInt'])
            )

        ## PS and LHE weights for ttbar events
        if isTTbar:
            ## PS Weight variations
            # PSWeight definitions can be found here: https://cms-nanoaod-integration.web.cern.ch/integration/master/mc94X_doc.html#PSWeight
            psweights = events['PSWeight']
            weights.add('ISR',
                np.ones(len(events)),
                psweights[:, 2], # (ISR=2, FSR=1)
                psweights[:, 0], # (ISR=0.5, FSR=1)
            )
            weights.add('FSR',
                np.ones(len(events)),
                psweights[:, 3], # (ISR=1, FSR=2)
                psweights[:, 1], # (ISR=1, FSR=0.5)
            )

            ## LHEScale Weight Variations
            # LHEScaleWeight definitions can be found here: https://cms-nanoaod-integration.web.cern.ch/integration/master/mc94X_doc.html#LHE
            lheweights = events['LHEScaleWeight']
            weights.add('FACTOR',
                np.ones(len(events)),
                lheweights[:, 5], # (muR=1, muF=2)
                lheweights[:, 3], # (muR=1, muF=0.5)
            )
            weights.add('RENORM',
                np.ones(len(events)),
                lheweights[:, 7], # (muR=2, muF=1)
                lheweights[:, 1], # (muR=0.5, muF=1)
            )
            weights.add('RENORM_FACTOR_SAME',
                np.ones(len(events)),
                lheweights[:, 8], # (muR=2, muF=2), RENORM_UP_FACTOR_UP
                lheweights[:, 0], # (muR=0.5, muF=0.5), RENORM_DW_FACTOR_DW
            )
            weights.add('RENORM_FACTOR_DIFF',
                np.ones(len(events)),
                lheweights[:, 6], # (muR=2, muF=0.5), RENORM_UP_FACTOR_DW
                lheweights[:, 2], # (muR=0.5, muF=2), RENORM_DW_FACTOR_UP
            )

    return weights    


def get_pdf_weights(df):
    if len(sorted(set(df['nLHEPdfWeight']))) > 1:
        raise ValueError('Differing number of PDF weights for events')
    pdfweights = ak.JaggedArray.fromcounts(df['nLHEPdfWeight'], df['LHEPdfWeight'])
    df['PDFWeights'] = pdfweights
    

def get_kin_weights(correction, genttbar):
    sl_evts = genttbar['SL']['TTbar'].counts > 0
    dl_evts = genttbar['DL']['TTbar'].counts > 0
    had_evts = genttbar['Had']['TTbar'].counts > 0

    np.random.seed(10) # sets seed so values from random distributions are reproducible (same as JER corrections)
    which_top_to_use = np.random.randint(2, size=genttbar.size) # top is 0, tbar is 1 
    var = correction['Var']
    dist = correction['Correction']

    wts = np.ones(genttbar.size)
    if 'thad_pt' in var:
            # set wts for semilep evts
        wts[sl_evts] = dist(genttbar['SL']['THad'].p4.pt[sl_evts].flatten())
            # set wts for dilep evts
        dl_pt = np.where(which_top_to_use[dl_evts], genttbar['DL']['Top'].p4.pt[dl_evts].flatten(), genttbar['DL']['Tbar'].p4.pt[dl_evts].flatten())
        wts[dl_evts] = dist(dl_pt)
            # set wts for had evts
        had_pt = np.where(which_top_to_use[had_evts], genttbar['Had']['Top'].p4.pt[had_evts].flatten(), genttbar['Had']['Tbar'].p4.pt[had_evts].flatten())
        wts[had_evts] = dist(had_pt)

    elif 'mtt_vs_thad_ctstar' in var:
            # set wts for semilep evts
        thad_p4, tlep_p4 = genttbar['SL']['THad'].p4[sl_evts].flatten(), genttbar['SL']['TLep'].p4[sl_evts].flatten()
        thad_ctstar, tlep_ctstar = make_vars.ctstar_flat(thad_p4, tlep_p4)
        wts[sl_evts] = dist(genttbar['SL']['TTbar'].p4.mass[sl_evts].flatten(), thad_ctstar)
            # set wts for dilep evts
        dl_top_p4, dl_tbar_p4 = genttbar['DL']['Top'].p4[dl_evts].flatten(), genttbar['DL']['Tbar'].p4[dl_evts].flatten()
        dl_top_ctstar, dl_tbar_ctstar = make_vars.ctstar_flat(dl_top_p4, dl_tbar_p4)
        dl_ctstar = np.where(which_top_to_use[dl_evts], dl_top_ctstar, dl_tbar_ctstar)
        wts[dl_evts] = dist(genttbar['DL']['TTbar'].p4.mass[dl_evts].flatten(), dl_ctstar)
            # set wts for had evts
        had_top_p4, had_tbar_p4 = genttbar['Had']['Top'].p4[had_evts].flatten(), genttbar['Had']['Tbar'].p4[had_evts].flatten()
        had_top_ctstar, had_tbar_ctstar = make_vars.ctstar_flat(had_top_p4, had_tbar_p4)
        had_ctstar = np.where(which_top_to_use[had_evts], had_top_ctstar, had_tbar_ctstar)
        wts[had_evts] = dist(genttbar['Had']['TTbar'].p4.mass[had_evts].flatten(), had_ctstar)

    else:
        raise ValueError("%s not supported for NNLO kinematic reweighting" % var)

    return wts

    


def get_comb_lepSF(year: str, lepton: str, corrections, pt: np.ndarray, eta: np.ndarray, shift='Central'):
    if not (shift == 'Central' or shift == 'Error'):
        raise ValueError('Shift value %s not defined' % shift)
    #set_trace()
    sf_dict = corrections[lepton]
    eta_ranges = sf_dict['eta_ranges']
    lepSFs = np.ones(pt.size)

    for idx, eta_range in enumerate(eta_ranges):
        mask = (eta >= eta_range[0]) & (eta < eta_range[1]) # find inds that are within given eta range
        if not mask.any(): continue # no values fall within eta range
        #set_trace()
        recoSF = sf_dict['Reco_ID']['Central']['eta_bin%i' % idx]
        trigSF = sf_dict['Trig']['Central']['eta_bin%i' % idx]
        lepSFs[mask] = recoSF(pt[mask])*trigSF(pt[mask])

    return lepSFs


def get_lepton_sf(year: str, lepton: str, corrections, pt, eta):
    sf_dict = corrections[lepton]
    eta_ranges = sf_dict['eta_ranges']
    reco_cen = np.ones(len(pt))
    reco_err = np.ones(len(pt))
    trig_cen = np.ones(len(pt))
    trig_err = np.ones(len(pt))

    for idx, eta_range in enumerate(eta_ranges):
        mask = (eta >= eta_range[0]) & (eta < eta_range[1]) # find inds that are within given eta range
        if not ak.any(mask): continue # no values fall within eta range
        recoSF_cen = sf_dict['Reco_ID']['Central']['eta_bin%i' % idx]
        recoSF_err = sf_dict['Reco_ID']['Error']['eta_bin%i' % idx]
        trigSF_cen = sf_dict['Trig']['Central']['eta_bin%i' % idx]
        trigSF_err = sf_dict['Trig']['Error']['eta_bin%i' % idx]
        reco_cen[mask] = recoSF_cen(pt[mask])
        reco_err[mask] = recoSF_err(pt[mask])
        trig_cen[mask] = trigSF_cen(pt[mask])
        trig_err[mask] = trigSF_err(pt[mask])

    output_SFs = {
        'RECO_CEN' : reco_cen,
        'RECO_ERR' : reco_err,
        'TRIG_CEN' : trig_cen,
        'TRIG_ERR' : trig_err,
    }

    return output_SFs
