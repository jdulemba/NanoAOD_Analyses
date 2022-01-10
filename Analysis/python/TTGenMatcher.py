from pdb import set_trace
import awkward as ak
import numpy as np
import compiled.pynusolver as pynusolver
from python.Permutations import make_perm_table
from numba import njit, objmode

@njit()
def find_nu(bleps, leptons, met, nu_array):
    for idx in range(bleps.shape[0]):
            # blep doesn't exist for this event
        if (bleps[idx, 3] == -999): continue
            ## run NS to get nschi2 to be used solver
        with objmode():
            pynusolver.run_nu_solver(leptons[idx], bleps[idx], met[idx], nu_array[idx]) # input order for nusolver is (lepton, jet, met, nu)

    return nu_array


def best_match(gen_hyp=None, jets=None, leptons=None, met=None):
    if gen_hyp is None:
        raise ValueError("Gen Objects gen_hyp needed for matching")
    if jets is None:
        raise ValueError("Reco jets needed for matching")
    if leptons is None:
        raise ValueError("Reco leptons needed for matching")
    if met is None:
        raise ValueError("Reco met needed for matching")

    if not ak.all(ak.num(gen_hyp) == 1):
        raise ValueError("Not all events for matching are semileptonic")

    jets_ak = ak.with_name(jets[["pt", "eta", "phi", "mass"]],"PtEtaPhiMLorentzVector")
    leps_ak = ak.with_name(leptons[["pt", "eta", "phi", "mass"]],"PtEtaPhiMLorentzVector")

        # init dict of objects
    matched_objects = {}

        # match jet closest to gen objects 
    for genobj in ['BHad', 'BLep', 'WJa', 'WJb']:
        genobj_ak = ak.with_name(gen_hyp[genobj][["pt", "eta", "phi", "mass"]],"PtEtaPhiMLorentzVector")
        jets_akc, genobj_akc = ak.unzip(ak.cartesian([jets_ak, genobj_ak], nested=False))
        deltaRs = jets_akc.delta_r(genobj_akc)  # find deltaRs between jets and gen object
        indexOfMin = ak.unflatten(ak.argmin(deltaRs, axis=1), ak.num(genobj_ak))
        passing_inds = deltaRs[indexOfMin] < 0.4

        matched_jets_inds = indexOfMin[passing_inds]
        matched_jets = jets[matched_jets_inds]

        ## add matched perm objects
        matched_objects[genobj] = ak.Array({
            'pt' : matched_jets.pt,
            'eta' : matched_jets.eta,
            'phi' : matched_jets.phi,
            'mass' : matched_jets.mass,
            'jetIdx' : matched_jets_inds, # index of jet that the gen object is matched to in the event
        }, with_name="PtEtaPhiMLorentzVector")
        
        # match lepton closest to gen lepton
    genlep_ak = ak.with_name(gen_hyp['Lepton'][["pt", "eta", "phi", "mass"]],"PtEtaPhiMLorentzVector")
    lep_akc, genlep_akc = ak.unzip(ak.cartesian([leps_ak, genlep_ak], nested=False))
    lepDRs = lep_akc.delta_r(genlep_akc)
    lepIdxOfMin = ak.unflatten(ak.argmin(lepDRs, axis=1), ak.num(genlep_ak))
    passing_inds = lepDRs[lepIdxOfMin] < 0.4
    matched_leps_inds = lepIdxOfMin[passing_inds]
    matched_leps = leptons[matched_leps_inds]
    
    ## add matched perm objects
    matched_objects['Lepton'] = ak.Array({key: matched_leps[key] for key in matched_leps.fields}, with_name="PtEtaPhiMLorentzVector")

        # solve for neutrino
    nu_array = np.zeros((len(ak.num(jets)), 4), dtype='float64')
            # convert all inputs into 2d numpy arrays of dtype=float64 (won't work if they're not float64)
    blep_inputs = np.stack((ak.to_numpy(ak.flatten(ak.fill_none(ak.pad_none(matched_objects['BLep'].px, 1), -999, axis=1))).astype('float64'), ak.to_numpy(ak.flatten(ak.fill_none(ak.pad_none(matched_objects['BLep'].py, 1), -999, axis=1))).astype('float64'),\
        ak.to_numpy(ak.flatten(ak.fill_none(ak.pad_none(matched_objects['BLep'].pz, 1), -999, axis=1))).astype('float64'), ak.to_numpy(ak.flatten(ak.fill_none(ak.pad_none(matched_objects['BLep'].energy, 1), -999, axis=1))).astype('float64')), axis=-1)
    lep_inputs = np.stack((ak.to_numpy(ak.flatten(ak.fill_none(ak.pad_none(matched_objects['Lepton'].px, 1), -999, axis=1))).astype('float64'), ak.to_numpy(ak.flatten(ak.fill_none(ak.pad_none(matched_objects['Lepton'].py, 1), -999, axis=1))).astype('float64'),\
        ak.to_numpy(ak.flatten(ak.fill_none(ak.pad_none(matched_objects['Lepton'].pz, 1), -999, axis=1))).astype('float64'), ak.to_numpy(ak.flatten(ak.fill_none(ak.pad_none(matched_objects['Lepton'].energy, 1), -999, axis=1))).astype('float64')), axis=-1)
    met_inputs = np.stack((ak.to_numpy(ak.fill_none(met.px, -999, axis=1)).astype('float64'), ak.to_numpy(ak.fill_none(met.py, -999, axis=1)).astype('float64')), axis=-1)
    nu_array = find_nu(bleps=blep_inputs, leptons=lep_inputs, met=met_inputs, nu_array=nu_array)

    valid_nu = ~((nu_array[:, 3] > 1e20) | (nu_array[:, 3] == 0)) # events that have a solution and matched blep

        # convert px, py, pz to pt, eta, phi
    nu_px, nu_py, nu_pz = nu_array[:, 0][valid_nu], nu_array[:, 1][valid_nu], nu_array[:, 2][valid_nu]
    nu_mom, nu_pt = np.sqrt(np.square(nu_px)+np.square(nu_py)+np.square(nu_pz)), np.sqrt(np.square(nu_px)+np.square(nu_py))
    nu_phi = np.arctan2(nu_py, nu_px)
    nu_eta = np.arcsinh(nu_pz/nu_pt)
    matched_objects['Nu'] = ak.Array({
        'pt' : ak.unflatten(nu_pt, valid_nu.astype(int)),
        'eta' : ak.unflatten(nu_eta, valid_nu.astype(int)),
        'phi' : ak.unflatten(nu_phi, valid_nu.astype(int)),
        'mass' : ak.zeros_like(ak.unflatten(nu_array[:, 0][valid_nu], valid_nu.astype(int))),
        'chi2' : ak.unflatten(nu_array[:, 3][valid_nu], valid_nu.astype(int)),
    }, with_name="PtEtaPhiMLorentzVector")

    matched_perm = make_perm_table(bhad=matched_objects['BHad'], blep=matched_objects['BLep'], wja=matched_objects['WJa'], wjb=matched_objects['WJb'], lepton=matched_objects['Lepton'], met=met, nu=matched_objects['Nu'])

    return matched_perm
