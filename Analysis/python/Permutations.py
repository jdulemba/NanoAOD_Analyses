import numpy as np
from pdb import set_trace
import awkward
from coffea.analysis_objects import JaggedCandidateArray

def make_perm_table(bhad, blep, wja, wjb, lepton, met, nu):
    '''
    Inputs:
        Jets, leptons, neutrinos, jet assignment object ordering, array of disrminiant probabilities (Total, mass discriminant, neutrino discrminant), and associated event weights
    Returns:
        awkward Table containing
            1: Jet objects (BLeps, BHads, WJas, WJbs)
            2: Leptons, Neutrinos
            3: Calculated probabilities (Total -> Prob, mass discriminant -> MassDiscr, neutrino discrminant -> NuDiscr)
            4: Number of jets in events (njets)
            5: Event weights (evt_wts)
    '''
    
        ## these attributes are based on those from URTTbar/interface/Permutation.h https://gitlab.cern.ch/jdulemba/URTTbar/-/blob/htt_analysis_2016legacydata_9410/interface/Permutation.h
    isWLepComplete = (lepton[~(lepton == None)].counts == 1) & (nu[~(nu == None)].counts == 1)
    isTLepComplete = isWLepComplete & (blep[~(blep == None)].counts == 1)

    wlep_p4 = lepton[isWLepComplete].p4 + nu[isWLepComplete].p4
    WLep = JaggedCandidateArray.candidatesfromcounts(
        counts = isWLepComplete.astype(int),
        pt     = wlep_p4.pt.flatten(),
        eta    = wlep_p4.eta.flatten(),
        phi    = wlep_p4.phi.flatten(),
        mass   = wlep_p4.mass.flatten(),
        charge = lepton[isWLepComplete].charge.flatten(),
    )#.pad(1) ?

    tlep_p4 = WLep[isTLepComplete].p4 + blep[isTLepComplete].p4
    TLep = JaggedCandidateArray.candidatesfromcounts(
        counts = isTLepComplete.astype(int),
        pt     = tlep_p4.pt.flatten(),
        eta    = tlep_p4.eta.flatten(),
        phi    = tlep_p4.phi.flatten(),
        mass   = tlep_p4.mass.flatten(),
    )

    ## Event categories
        # fill [None] events with values for easier comparisons
    bhad_jetIdx = bhad.jetIdx.fillna(-999)
    blep_jetIdx = blep.jetIdx.fillna(-999)
    wja_jetIdx  = wja.jetIdx.fillna(-999)
    wjb_jetIdx  = wjb.jetIdx.fillna(-999)
        # merged jets event categories
            # only bhad and blep merged
    Merged_BHadBLep = (bhad_jetIdx >= 0) & (bhad_jetIdx == blep_jetIdx) & (bhad_jetIdx != wja_jetIdx) & (bhad_jetIdx != wjb_jetIdx)
            # only bhad and wja merged
    Merged_BHadWJa = (bhad_jetIdx >= 0) & (bhad_jetIdx == wja_jetIdx) & (bhad_jetIdx != blep_jetIdx) & (bhad_jetIdx != wjb_jetIdx)
            # only bhad and wjb merged
    Merged_BHadWJb = (bhad_jetIdx >= 0) & (bhad_jetIdx == wjb_jetIdx) & (bhad_jetIdx != blep_jetIdx) & (bhad_jetIdx != wja_jetIdx)
            # only blep and wja merged
    Merged_BLepWJa = (blep_jetIdx >= 0) & (blep_jetIdx == wja_jetIdx) & (blep_jetIdx != bhad_jetIdx) & (blep_jetIdx != wjb_jetIdx)
            # only blep and wjb merged
    Merged_BLepWJb = (blep_jetIdx >= 0) & (blep_jetIdx == wjb_jetIdx) & (blep_jetIdx != bhad_jetIdx) & (blep_jetIdx != wja_jetIdx)
            # only wja and wjb merged
    Merged_WJets = (wja_jetIdx >= 0) & (wja_jetIdx == wjb_jetIdx) & (wja_jetIdx != bhad_jetIdx) & (wja_jetIdx != blep_jetIdx)
            # only two jets merged
    Merged_Event = (Merged_BHadBLep) | (Merged_BHadWJa) | (Merged_BHadWJb) | (Merged_BLepWJa) | (Merged_BLepWJb) | (Merged_WJets)
            # only bhad, blep, wja merged
    Merged_BHadBLepWJa = (bhad_jetIdx >= 0) & (bhad_jetIdx == blep_jetIdx) & (bhad_jetIdx == wja_jetIdx) & (bhad_jetIdx != wjb_jetIdx)
            # only bhad, blep, wjb merged
    Merged_BHadBLepWJb = (bhad_jetIdx >= 0) & (bhad_jetIdx == blep_jetIdx) & (bhad_jetIdx == wjb_jetIdx) & (bhad_jetIdx != wja_jetIdx)
            # only bhad, wja, wjb merged
    Merged_BHadWJaWJb = (bhad_jetIdx >= 0) & (bhad_jetIdx == wja_jetIdx) & (bhad_jetIdx == wjb_jetIdx) & (bhad_jetIdx != blep_jetIdx)
            # only blep, wja, wjb merged
    Merged_BLepWJaWJb = (blep_jetIdx >= 0) & (blep_jetIdx == wja_jetIdx) & (blep_jetIdx == wjb_jetIdx) & (blep_jetIdx != bhad_jetIdx)
        #

        # lost jet event categories
            # only bhad is lost, other jets exist and are resolved
    Lost_BHad = (bhad_jetIdx < 0) & (blep_jetIdx >= 0) & (wja_jetIdx >= 0) & (wjb_jetIdx >= 0) & (blep_jetIdx != wja_jetIdx) & (blep_jetIdx != wjb_jetIdx) & (wja_jetIdx != wjb_jetIdx)
            # only blep is lost, other jets exist and are resolved
    Lost_BLep = (blep_jetIdx < 0) & (bhad_jetIdx >= 0) & (wja_jetIdx >= 0) & (wjb_jetIdx >= 0) & (bhad_jetIdx != wja_jetIdx) & (bhad_jetIdx != wjb_jetIdx) & (wja_jetIdx != wjb_jetIdx)
            # only wja is lost, other jets exist and are resolved
    Lost_WJa = (wja_jetIdx < 0) & (bhad_jetIdx >= 0) & (blep_jetIdx >= 0) & (wjb_jetIdx >= 0) & (bhad_jetIdx != blep_jetIdx) & (bhad_jetIdx != wjb_jetIdx) & (blep_jetIdx != wjb_jetIdx)
            # only wjb is lost, other jets exist and are resolved
    Lost_WJb = (wjb_jetIdx < 0) & (bhad_jetIdx >= 0) & (blep_jetIdx >= 0) & (wja_jetIdx >= 0) & (bhad_jetIdx != blep_jetIdx) & (bhad_jetIdx != wja_jetIdx) & (blep_jetIdx != wja_jetIdx)
            # only one jet is lost, others exist and are resolved
    Lost_Event = (Lost_BHad) | (Lost_BLep) | (Lost_WJa) | (Lost_WJb)
    ##

    n_perm_matches = awkward.JaggedArray.fromcounts(np.ones(bhad.size, dtype=int), bhad[~(bhad == None)].counts + blep[~(blep == None)].counts + wja[~(wja == None)].counts + wjb[~(wjb == None)].counts)
    isEmpty = (bhad_jetIdx < 0) & (blep_jetIdx < 0) & (wja_jetIdx < 0) & (wjb_jetIdx < 0)

        # find number of unique matches
    jetIdx_stack = np.stack((blep_jetIdx.flatten(), bhad_jetIdx.flatten(), wja_jetIdx.flatten(), wjb_jetIdx.flatten()), axis=1)
    unique_matches = awkward.JaggedArray.fromcounts(np.ones(bhad.size, dtype=int), np.array([len(list(set([ind for ind in inds if ind >= 0]))) for inds in jetIdx_stack.tolist()]))


        # create WHad
            # fill [None] entries with zeros
    wja_dict = {
        'px' : wja.p4.x.fillna(0.).flatten(),
        'py' : wja.p4.y.fillna(0.).flatten(),
        'pz' : wja.p4.z.fillna(0.).flatten(),
        'energy' : wja.p4.energy.fillna(0.).flatten(),
    }
    wja_p4 = JaggedCandidateArray.candidatesfromcounts(
        counts = wja.counts,
        **wja_dict
    )
    wjb_dict = {
        'px' : wjb.p4.x.fillna(0.).flatten(),
        'py' : wjb.p4.y.fillna(0.).flatten(),
        'pz' : wjb.p4.z.fillna(0.).flatten(),
        'energy' : wjb.p4.energy.fillna(0.).flatten(),
    }
    wjb_p4 = JaggedCandidateArray.candidatesfromcounts(
        counts = wjb.counts,
        **wjb_dict
    )

            ## init variables
    whad_px = np.full_like(bhad.counts, np.nan, dtype=np.double) 
    whad_py = np.full_like(bhad.counts, np.nan, dtype=np.double) 
    whad_pz = np.full_like(bhad.counts, np.nan, dtype=np.double) 
    whad_E  = np.full_like(bhad.counts, np.nan, dtype=np.double) 

            # inds where only WJb p4 is used as WHad
    use_wjb_inds = (Merged_BHadWJa | Merged_BLepWJa | Merged_WJets | Lost_WJa).flatten()
    whad_px[use_wjb_inds] = wjb[use_wjb_inds].p4.x.flatten()
    whad_py[use_wjb_inds] = wjb[use_wjb_inds].p4.y.flatten()
    whad_pz[use_wjb_inds] = wjb[use_wjb_inds].p4.z.flatten()
    whad_E[use_wjb_inds]  = wjb[use_wjb_inds].p4.energy.flatten()

            # inds where only WJa p4 is used as WHad
    use_wja_inds = (Merged_BHadWJb | Merged_BLepWJb | Lost_WJb).flatten()
    whad_px[use_wja_inds] = wja[use_wja_inds].p4.x.flatten()
    whad_py[use_wja_inds] = wja[use_wja_inds].p4.y.flatten()
    whad_pz[use_wja_inds] = wja[use_wja_inds].p4.z.flatten()
    whad_E[use_wja_inds]  = wja[use_wja_inds].p4.energy.flatten()

            # inds where combined p4 from WJa and WJb is used as WHad (all other inds)
    use_comb_inds = ~(use_wjb_inds | use_wja_inds)
    comb_wjets_p4 = wja_p4[use_comb_inds].p4 + wjb_p4[use_comb_inds].p4
    whad_px[use_comb_inds] = comb_wjets_p4.x.flatten()
    whad_px[whad_px == 0.] = np.nan
    whad_py[use_comb_inds] = comb_wjets_p4.y.flatten()
    whad_py[whad_py == 0.] = np.nan
    whad_pz[use_comb_inds] = comb_wjets_p4.z.flatten()
    whad_pz[whad_pz == 0.] = np.nan
    whad_E[use_comb_inds]  = comb_wjets_p4.energy.flatten()
    whad_E[whad_E == 0.]   = np.nan

    WHad = JaggedCandidateArray.candidatesfromcounts(
        counts = (~np.isnan(whad_px)).astype(int),
        px     = whad_px[(~np.isnan(whad_px))].flatten(),
        py     = whad_py[(~np.isnan(whad_py))].flatten(),
        pz     = whad_pz[(~np.isnan(whad_pz))].flatten(),
        energy = whad_E[(~np.isnan(whad_E))].flatten(),
    )

    isWHadComplete = (WHad.counts == 1)
    isTHadComplete = (isWHadComplete) & (bhad[~(bhad == None)].counts == 1)

    thad_p4 = WHad[isTHadComplete].p4 + bhad[isTHadComplete].p4
    THad = JaggedCandidateArray.candidatesfromcounts(
        counts = isTHadComplete.astype(int),
        pt     = thad_p4.pt.flatten(),
        eta    = thad_p4.eta.flatten(),
        phi    = thad_p4.phi.flatten(),
        mass   = thad_p4.mass.flatten(),
    )

    isComplete = isTHadComplete & isTLepComplete
    ttbar_p4 = THad[isComplete].p4 + TLep[isComplete].p4
    TTbar = JaggedCandidateArray.candidatesfromcounts(
        counts = isComplete.astype(int),
        pt     = ttbar_p4.pt.flatten(),
        eta    = ttbar_p4.eta.flatten(),
        phi    = ttbar_p4.phi.flatten(),
        mass   = ttbar_p4.mass.flatten(),
    )

        ## Combine everything into a single table, all objects are JaggedArrays
    permutations = awkward.Table(
        BHad = bhad[~(bhad == None)],
        BLep = blep[~(blep == None)],
        WJa = wja[~(wja == None)],
        WJb = wjb[~(wjb == None)],
        Lepton = lepton[~(lepton == None)],
        MET = met[~(met == None)],
        Nu = nu[~(nu == None)],
        WLep = WLep,
        TLep = TLep,
        WHad = WHad,
        THad = THad,
        TTbar = TTbar,
        n_perm_matches = n_perm_matches,
        isEmpty = isEmpty,
        unique_matches = unique_matches,
        Merged_BHadBLep = Merged_BHadBLep,
        Merged_BHadWJa = Merged_BHadWJa,
        Merged_BHadWJb = Merged_BHadWJb,
        Merged_BLepWJa = Merged_BLepWJa,
        Merged_BLepWJb = Merged_BLepWJb,
        Merged_WJets = Merged_WJets,
        Merged_Event = Merged_Event,
        Lost_BHad = Lost_BHad,
        Lost_BLep = Lost_BLep,
        Lost_WJa = Lost_WJa,
        Lost_WJb = Lost_WJb,
        Lost_Event = Lost_Event,
    )

    #set_trace()
    
    return permutations


