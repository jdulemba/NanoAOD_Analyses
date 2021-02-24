import numpy as np
from pdb import set_trace
import awkward as ak

def make_perm_table(bhad, blep, wja, wjb, lepton, met, nu):
    '''
    Inputs:
        Jets, leptons, neutrinos, jet assignment object ordering, array of disrminiant probabilities (Total, mass discriminant, neutrino discrminant), and associated event weights
    Returns:
        awkward Table containing
            1: Jet objects (BLeps, BHads, WJas, WJbs)
            2: Leptons, Neutrinos
            3: Calculated probabilities (Total -> Prob, mass discriminant -> MassDiscr, neutrino discrminant -> NuDiscr)
    '''

        ## these attributes are based on those from URTTbar/interface/Permutation.h https://gitlab.cern.ch/jdulemba/URTTbar/-/blob/htt_analysis_2016legacydata_9410/interface/Permutation.h
    isWLepComplete = (ak.num(lepton.pt) == 1) & (ak.num(nu.pt) == 1)
    isTLepComplete = isWLepComplete & (ak.num(blep.pt) == 1)

    WLep = ak.Array({
        'pt'    : ak.fill_none( (ak.mask(lepton, isWLepComplete)+ak.mask(nu, isWLepComplete)).pt, []),
        'eta'   : ak.fill_none( (ak.mask(lepton, isWLepComplete)+ak.mask(nu, isWLepComplete)).eta, []),
        'phi'   : ak.fill_none( (ak.mask(lepton, isWLepComplete)+ak.mask(nu, isWLepComplete)).phi, []),
        'mass'  : ak.fill_none( (ak.mask(lepton, isWLepComplete)+ak.mask(nu, isWLepComplete)).mass, []),
        'charge': ak.fill_none( (ak.mask(lepton, isWLepComplete)).charge, []),
    }, with_name="PtEtaPhiMLorentzVector")
    TLep = ak.Array({
        'pt'    : ak.fill_none( (ak.mask(WLep, isTLepComplete)+ak.mask(blep, isTLepComplete)).pt, []),
        'eta'   : ak.fill_none( (ak.mask(WLep, isTLepComplete)+ak.mask(blep, isTLepComplete)).eta, []),
        'phi'   : ak.fill_none( (ak.mask(WLep, isTLepComplete)+ak.mask(blep, isTLepComplete)).phi, []),
        'mass'  : ak.fill_none( (ak.mask(WLep, isTLepComplete)+ak.mask(blep, isTLepComplete)).mass, []),
    }, with_name="PtEtaPhiMLorentzVector")


    ## Event categories
        # fill empty [] events with values for easier comparisons
    bhad_jetIdx = ak.fill_none(ak.pad_none(bhad, 1).jetIdx, -999)
    blep_jetIdx = ak.fill_none(ak.pad_none(blep, 1).jetIdx, -999)
    wja_jetIdx = ak.fill_none(ak.pad_none(wja, 1).jetIdx, -999)
    wjb_jetIdx = ak.fill_none(ak.pad_none(wjb, 1).jetIdx, -999)

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

    n_perm_matches = ak.unflatten(ak.num(bhad.pt)+ak.num(blep.pt)+ak.num(wja.pt)+ak.num(wjb.pt), np.ones(len(bhad.pt), dtype=int))
    #isEmpty = (bhad_jetIdx < 0) & (blep_jetIdx < 0) & (wja_jetIdx < 0) & (wjb_jetIdx < 0)

        # find number of unique matches
    jetIdx_stack = np.stack( (ak.to_numpy(ak.flatten(bhad_jetIdx)), ak.to_numpy(ak.flatten(blep_jetIdx)), ak.to_numpy(ak.flatten(wja_jetIdx)), ak.to_numpy(ak.flatten(wjb_jetIdx)) ), axis=1)
    unique_matches = ak.unflatten(np.array([len(list(set([ind for ind in inds if ind >= 0]))) for inds in jetIdx_stack.tolist()]), np.ones(len(bhad.pt), dtype=int))


    # create WHad
        # init variables
    whad_pt  = np.zeros(len(bhad.pt))
    whad_eta = np.zeros(len(bhad.pt))
    whad_phi = np.zeros(len(bhad.pt))
    whad_mass= np.zeros(len(bhad.pt))

            # inds where only WJb p4 is used as WHad
    use_wjb_inds = ak.flatten(Merged_BHadWJa | Merged_BLepWJa | Merged_WJets | Lost_WJa)
    whad_pt[use_wjb_inds]  = ak.to_numpy(ak.flatten(ak.fill_none(ak.pad_none(wjb.pt, 1), np.nan)[use_wjb_inds]))
    whad_eta[use_wjb_inds] = ak.to_numpy(ak.flatten(ak.fill_none(ak.pad_none(wjb.eta, 1), np.nan)[use_wjb_inds]))
    whad_phi[use_wjb_inds] = ak.to_numpy(ak.flatten(ak.fill_none(ak.pad_none(wjb.phi, 1), np.nan)[use_wjb_inds]))
    whad_mass[use_wjb_inds]= ak.to_numpy(ak.flatten(ak.fill_none(ak.pad_none(wjb.mass, 1), np.nan)[use_wjb_inds]))

            # inds where only WJa p4 is used as WHad
    use_wja_inds = ak.flatten(Merged_BHadWJb | Merged_BLepWJb | Lost_WJb)
    whad_pt[use_wja_inds]  = ak.to_numpy(ak.flatten(ak.fill_none(ak.pad_none(wja.pt, 1), np.nan)[use_wja_inds]))
    whad_eta[use_wja_inds] = ak.to_numpy(ak.flatten(ak.fill_none(ak.pad_none(wja.eta, 1), np.nan)[use_wja_inds]))
    whad_phi[use_wja_inds] = ak.to_numpy(ak.flatten(ak.fill_none(ak.pad_none(wja.phi, 1), np.nan)[use_wja_inds]))
    whad_mass[use_wja_inds]= ak.to_numpy(ak.flatten(ak.fill_none(ak.pad_none(wja.mass, 1), np.nan)[use_wja_inds]))

            # inds where combined p4 from WJa and WJb is used as WHad (all other inds)
    use_comb_inds = ~(use_wjb_inds | use_wja_inds)
    whad_pt[use_comb_inds]  = ak.to_numpy(ak.flatten( (ak.fill_none(ak.pad_none( wja, 1), 0) + ak.fill_none(ak.pad_none( wjb, 1), 0)).pt[use_comb_inds] ))
    whad_eta[use_comb_inds] = ak.to_numpy(ak.flatten( (ak.fill_none(ak.pad_none( wja, 1), 0) + ak.fill_none(ak.pad_none( wjb, 1), 0)).eta[use_comb_inds] ))
    whad_phi[use_comb_inds] = ak.to_numpy(ak.flatten( (ak.fill_none(ak.pad_none( wja, 1), 0) + ak.fill_none(ak.pad_none( wjb, 1), 0)).phi[use_comb_inds] ))
    whad_mass[use_comb_inds]= ak.to_numpy(ak.flatten( (ak.fill_none(ak.pad_none( wja, 1), 0) + ak.fill_none(ak.pad_none( wjb, 1), 0)).mass[use_comb_inds] ))
    whad_pt[whad_pt == 0.] = np.nan
    whad_eta[whad_eta == 0.] = np.nan
    whad_phi[whad_phi == 0.] = np.nan
    whad_mass[whad_mass == 0.] = np.nan
    
    WHad = ak.Array({
        'pt'    : ak.unflatten(whad_pt[~np.isnan(whad_pt)], (~np.isnan(whad_pt)).astype(int)),
        'eta'   : ak.unflatten(whad_eta[~np.isnan(whad_eta)], (~np.isnan(whad_eta)).astype(int)),
        'phi'   : ak.unflatten(whad_phi[~np.isnan(whad_phi)], (~np.isnan(whad_phi)).astype(int)),
        'mass'  : ak.unflatten(whad_mass[~np.isnan(whad_mass)], (~np.isnan(whad_mass)).astype(int)),
        'charge': -1*ak.fill_none(ak.mask(WLep, ~np.isnan(whad_mass)).charge, []), # opposite charge as WLep for events that exist
    }, with_name="PtEtaPhiMLorentzVector")


    isWHadComplete = (ak.num(WHad.pt) == 1)
    isTHadComplete = (isWHadComplete) & (ak.num(bhad.pt) == 1)

    #set_trace()

    # create THad
    THad = ak.Array({
        'pt'    : ak.fill_none( (ak.mask(WHad, isTHadComplete)+ak.mask(bhad, isTHadComplete)).pt, []),
        'eta'   : ak.fill_none( (ak.mask(WHad, isTHadComplete)+ak.mask(bhad, isTHadComplete)).eta, []),
        'phi'   : ak.fill_none( (ak.mask(WHad, isTHadComplete)+ak.mask(bhad, isTHadComplete)).phi, []),
        'mass'  : ak.fill_none( (ak.mask(WHad, isTHadComplete)+ak.mask(bhad, isTHadComplete)).mass, []),
    }, with_name="PtEtaPhiMLorentzVector")

    # create TTbar
    isComplete = isTHadComplete & isTLepComplete
    TTbar = ak.Array({
        'pt'    : ak.fill_none( (ak.mask(THad, isComplete)+ak.mask(TLep, isComplete)).pt, []),
        'eta'   : ak.fill_none( (ak.mask(THad, isComplete)+ak.mask(TLep, isComplete)).eta, []),
        'phi'   : ak.fill_none( (ak.mask(THad, isComplete)+ak.mask(TLep, isComplete)).phi, []),
        'mass'  : ak.fill_none( (ak.mask(THad, isComplete)+ak.mask(TLep, isComplete)).mass, []),
    }, with_name="PtEtaPhiMLorentzVector")

    
    #set_trace()    
        ## Combine everything into a single table, all objects are JaggedArrays
    permutations = ak.zip({
        "BHad" : bhad,
        "BLep" : blep,
        "WJa" : wja,
        "WJb" : wjb,
        "Lepton" : lepton,
        "MET" : met,
        "Nu" : nu,
        "WLep" : WLep,
        "TLep" : TLep,
        "WHad" : WHad,
        "THad" : THad,
        "TTbar" : TTbar,
        "n_perm_matches" : n_perm_matches,
        #"isEmpty" : isEmpty,
        "unique_matches" : unique_matches,
        "Merged_BHadBLep" : Merged_BHadBLep,
        "Merged_BHadWJa" : Merged_BHadWJa,
        "Merged_BHadWJb" : Merged_BHadWJb,
        "Merged_BLepWJa" : Merged_BLepWJa,
        "Merged_BLepWJb" : Merged_BLepWJb,
        "Merged_WJets" : Merged_WJets,
        "Merged_Event" : Merged_Event,
        "Lost_BHad" : Lost_BHad,
        "Lost_BLep" : Lost_BLep,
        "Lost_WJa" : Lost_WJa,
        "Lost_WJb" : Lost_WJb,
        "Lost_Event" : Lost_Event,
    }, depth_limit=1)

    #set_trace()
    
    return permutations



def compare_matched_best_perms(mp, bp, njets, bp_mask=None):
#def compare_matched_best_perms(mp, bp, njets, mp_mask=None, bp_mask=None):
    '''
    Compare object assignments across two permutations.
    Inputs: matched perm, best perm, njets category, matched perm mask, best perm mask
    '''
    set_trace()
    if mp['BLep'].size != bp['BLep'][bp_mask].size:
    #if mp['BLep'][mp_mask].size != bp['BLep'][bp_mask].size:
        raise ValueError("Permutations must have the same size in order to be compared!")
    if not ((njets == '3Jets') or (njets == '4PJets')):
        raise ValueError("Only 3Jets or 4+ jets categorizations are supported")

    #set_trace()
    #mp_blep_idx = mp['BLep'][mp_mask].pad(1).jetIdx.fillna(-999).flatten()
    #mp_bhad_idx = mp['BHad'][mp_mask].pad(1).jetIdx.fillna(-999).flatten()
    #mp_wja_idx = mp['WJa'][mp_mask].pad(1).jetIdx.fillna(-999).flatten()
    #mp_wjb_idx = mp['WJb'][mp_mask].pad(1).jetIdx.fillna(-999).flatten()
    #mp_lep_pt = mp['Lepton'][mp_mask].pad(1).pt.fillna(-999).flatten()
    mp_blep_idx = mp['BLep'].pad(1).jetIdx.fillna(-999).flatten()
    mp_bhad_idx = mp['BHad'].pad(1).jetIdx.fillna(-999).flatten()
    mp_wja_idx = mp['WJa'].pad(1).jetIdx.fillna(-999).flatten()
    mp_wjb_idx = mp['WJb'].pad(1).jetIdx.fillna(-999).flatten()
    mp_lep_pt = mp['Lepton'].pad(1).pt.fillna(-999).flatten()

    bp_blep_idx = bp['BLep'][bp_mask].pad(1).jetIdx.fillna(-999).flatten()
    bp_bhad_idx = bp['BHad'][bp_mask].pad(1).jetIdx.fillna(-999).flatten()
    bp_wja_idx = bp['WJa'][bp_mask].pad(1).jetIdx.fillna(-999).flatten()
    bp_wjb_idx = bp['WJb'][bp_mask].pad(1).jetIdx.fillna(-999).flatten()
    bp_lep_pt = bp['Lepton'][bp_mask].pad(1).pt.fillna(-999).flatten()

        # index comparisons
    same_blep = (mp_blep_idx == bp_blep_idx) & (mp_blep_idx >= 0)
    same_bhad = (mp_bhad_idx == bp_bhad_idx) & (mp_bhad_idx >= 0)
    same_wja = (mp_wja_idx == bp_wja_idx) & (mp_wja_idx >= 0)
    same_wjb = (mp_wjb_idx == bp_wjb_idx) & (mp_wjb_idx >= 0)

    same_bs = same_blep & same_bhad

    cats = np.zeros(mp['BLep'].size) # 0 == '' (no gen matching), 1 == 'right', 2 == 'matchable', 3 == 'unmatchable', 4 == 'noslep'
    #cats = np.zeros(mp['BLep'][mp_mask].size) # 0 == '' (no gen matching), 1 == 'right', 2 == 'matchable', 3 == 'unmatchable', 4 == 'noslep'
    if njets == '3Jets':
        valid_evts = ((mp['TTbar'].counts > 0) & ((mp['unique_matches'] >= 3).flatten()))

            # merged events
        merged_evts = valid_evts & mp['Merged_Event'].flatten()
        correct_merged = merged_evts & (mp['Merged_BHadWJa'] | mp['Merged_BHadWJb'] | mp['Merged_WJets']).flatten()
        wrong_merged = merged_evts & ~(mp['Merged_BHadWJa'] | mp['Merged_BHadWJb'] | mp['Merged_WJets']).flatten()

            # lost events
        lost_evts = valid_evts & mp['Lost_Event'].flatten()
        correct_lost = lost_evts & ((mp['Lost_WJa'] | mp['Lost_WJb'])).flatten()
        wrong_lost = lost_evts & ~((mp['Lost_WJa'] | mp['Lost_WJb'])).flatten()

        right_matching = same_bs & (((bp_wja_idx == mp_wja_idx) | (bp_wja_idx == mp_wjb_idx)) & (bp_wja_idx >= 0)) # bs are matched correctly, bp wjet is matched to one of the matched perm wjets

        # event categorization
            # unmatchable events
        unmatchable_evts = (~valid_evts | wrong_merged | wrong_lost)
            # right events
        right_perm_evts = (correct_lost & right_matching) | (correct_merged & right_matching) # matched perm is correct event type and right object matching
            # matchable events
        matchable_evts = (correct_lost & ~right_matching) | (correct_merged & ~right_matching) # matched perm is correct event type but wrong object matching

    else:
        valid_evts = ((mp['TTbar'].counts > 0) & ((mp['unique_matches'] == 4).flatten()))
        isWHadCorrect = ((bp_wja_idx == mp_wja_idx) & (bp_wjb_idx == mp_wjb_idx)) | ((bp_wja_idx == mp_wjb_idx) & (bp_wjb_idx == mp_wja_idx))
        isTHadCorrect = same_bhad & isWHadCorrect
        isTLepCorrect = same_blep & (bp_lep_pt == mp_lep_pt)
        isCorrect = isTLepCorrect & isTHadCorrect

        # event categorization
            # unmatchable events
        unmatchable_evts = ~valid_evts
            # right events
        right_perm_evts = (isCorrect & valid_evts)
            #matchable events
        matchable_evts = (~isCorrect & valid_evts)


        # check that there's no overlap in event categories
    all_cats_as_int = unmatchable_evts.astype(int)+right_perm_evts.astype(int)+matchable_evts.astype(int)
    if (all_cats_as_int != 1).any():
        raise ValueError("Events %s have overlapping categories!\nEvents %s have no categories!" % (np.where(all_cats_as_int > 1), np.where(all_cats_as_int == 0)))

        # set values in category array
    cats[right_perm_evts] = 1 # matched perm is correct event type and right object matching
    cats[matchable_evts] = 2 # matched perm is correct event type but wrong object matching
    cats[unmatchable_evts] = 3 # matched perm isn't valid or has wrong event type

    return cats


