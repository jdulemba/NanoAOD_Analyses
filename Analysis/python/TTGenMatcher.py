from pdb import set_trace
import awkward
from coffea.analysis_objects import JaggedCandidateArray
import numpy as np
import compiled.pynusolver as pynusolver
from python.Permutations import make_perm_table

remove_fast = lambda x : x.split('fast_')[-1]

def best_match(gen_hyp=None, jets=None, leptons=None, met=None):
    if gen_hyp is None:
        raise ValueError("Gen Objects gen_hyp needed for matching")
    if jets is None:
        raise ValueError("Reco jets needed for matching")
    if leptons is None:
        raise ValueError("Reco leptons needed for matching")
    if met is None:
        raise ValueError("Reco met needed for matching")

    if (gen_hyp['SL']['TTbar'].counts != 1).any():
        raise ValueError("Not all events for matching are semileptonic")

        # init dicts of counts and variables
    matched_dict_counts = {}
    matched_dict_vars ={}


        # match jet closest to gen objects 
    for genobj in ['BHad', 'BLep', 'WJa', 'WJb']:
        deltaRs = jets.p4.delta_r(gen_hyp['SL'][genobj].p4.flatten())
        indexOfMin = deltaRs.argmin()
        passing_inds = deltaRs[indexOfMin] < 0.4
        matched_jets_inds = indexOfMin[passing_inds]
        matched_jets = jets[matched_jets_inds]

        matched_dict_counts.update({genobj : matched_jets.counts})
        matched_dict_vars.update({genobj : {
            'pt' : matched_jets.pt.flatten(),
            'eta' : matched_jets.eta.flatten(),
            'phi' : matched_jets.phi.flatten(),
            'mass' : matched_jets.mass.flatten(),
            #'hadronFlav' : matched_jets.hadronFlav.flatten(),
            'jetIdx' : matched_jets_inds.flatten(), # index of jet that the gen object is matched to in the event
            
        }})

        # match lepton closest to gen lepton
    lepDRs = leptons.p4.delta_r(gen_hyp['SL']['Lepton'].p4.flatten())
    lepIdxOfMin = lepDRs.argmin()
    passing_inds = lepDRs[lepIdxOfMin] < 0.4
    matched_leps_inds = lepIdxOfMin[passing_inds]
    matched_leps = leptons[matched_leps_inds]

    matched_dict_counts.update({'Lepton' : matched_leps.counts})
    matched_dict_vars.update({'Lepton' : {remove_fast(key) : matched_leps[key].flatten() for key in matched_leps.columns if key != 'p4'}})
    
        ## create matched perm objects
    matched_BHad = JaggedCandidateArray.candidatesfromcounts(
        counts= matched_dict_counts['BHad'],
        **matched_dict_vars['BHad']
    ).pad(1)
    matched_BLep = JaggedCandidateArray.candidatesfromcounts(
        counts= matched_dict_counts['BLep'],
        **matched_dict_vars['BLep']
    ).pad(1)
    matched_WJa = JaggedCandidateArray.candidatesfromcounts(
        counts= matched_dict_counts['WJa'],
        **matched_dict_vars['WJa']
    ).pad(1)
    matched_WJb = JaggedCandidateArray.candidatesfromcounts(
        counts= matched_dict_counts['WJb'],
        **matched_dict_vars['WJb']
    ).pad(1)
    matched_Lep = JaggedCandidateArray.candidatesfromcounts(
        counts= matched_dict_counts['Lepton'],
        **matched_dict_vars['Lepton']
    ).pad(1)

        # solve for neutrino
    nu_array = np.zeros((jets.counts.size, 4), dtype='float64')
    for idx, blep in enumerate(matched_BLep):
        if blep.size < 1: continue
        if blep.size > 1: raise ValueError("More than one matched blep. Investigate")

            # must specify same dtype for all arrays
        blep_inputs = np.array([blep.p4.x[0], blep.p4.y[0], blep.p4.z[0], blep.p4.energy[0]], dtype='float64')
        lep_inputs = np.array([matched_Lep[idx].p4.x[0], matched_Lep[idx].p4.y[0], matched_Lep[idx].p4.z[0], matched_Lep[idx].p4.energy[0]], dtype='float64')
        met_inputs = np.array([met[idx].p4.x[0], met[idx].p4.y[0]], dtype='float64')
        nu = nu_array[idx]
        pynusolver.run_nu_solver(lep_inputs, blep_inputs, met_inputs, nu) # input order for nusolver is (lepton, jet, met, nu), returns Nu(px, py, pz, chi2)

    valid_nu = ~((nu_array[:, 3] > 1e20) | (nu_array[:, 3] == 0)) # events that have a solution and matched blep
    matched_Nu = JaggedCandidateArray.candidatesfromcounts(
        counts = valid_nu.astype(int),
        px = nu_array[:, 0][valid_nu],
        py = nu_array[:, 1][valid_nu],
        pz = nu_array[:, 2][valid_nu],
        mass = np.zeros(valid_nu.sum()),
        chi2 = nu_array[:, 3][valid_nu],
    ).pad(1)

    matched_perm = make_perm_table(bhad=matched_BHad, blep=matched_BLep, wja=matched_WJa, wjb=matched_WJb, lepton=matched_Lep, met=met, nu=matched_Nu)
    #set_trace()

    return matched_perm
