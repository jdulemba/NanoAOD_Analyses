from pdb import set_trace
import awkward
from coffea.analysis_objects import JaggedCandidateArray
import numpy as np


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
    for genobj in ['Had_B', 'Lep_B', 'First_Had', 'Second_Had']:
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
            'hadronFlav' : matched_jets.hadronFlav.flatten(),
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
    
    #set_trace()

        ## create matched perm objects
    matched_BHad = JaggedCandidateArray.candidatesfromcounts(
        counts= matched_dict_counts['Had_B'],
        **matched_dict_vars['Had_B']
    )
    matched_BLep = JaggedCandidateArray.candidatesfromcounts(
        counts= matched_dict_counts['Lep_B'],
        **matched_dict_vars['Lep_B']
    )
    matched_WJa = JaggedCandidateArray.candidatesfromcounts(
        counts= matched_dict_counts['First_Had'],
        **matched_dict_vars['First_Had']
    )
    matched_WJb = JaggedCandidateArray.candidatesfromcounts(
        counts= matched_dict_counts['Second_Had'],
        **matched_dict_vars['Second_Had']
    )
    matched_Lep = JaggedCandidateArray.candidatesfromcounts(
        counts= matched_dict_counts['Lepton'],
        **matched_dict_vars['Lepton']
    )

    matched_perm = awkward.Table(
        BHad  = matched_BHad,
        BLep  = matched_BLep,
        WJa   = matched_WJa,
        WJb   = matched_WJb,
        Lepton= matched_Lep,
        MET   = met,
    )

    return matched_perm
