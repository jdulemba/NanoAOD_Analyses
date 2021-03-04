from pdb import set_trace
import awkward as ak
import python.Filters_and_Triggers as Filters_and_Triggers
import python.IDJet as IDJet
import python.IDMuon as IDMuon
import python.IDElectron as IDElectron
import numpy as np
from coffea.analysis_tools import PackedSelection

def select_jets(jets, muons, electrons, year, cutflow=None):

        ## pt and eta cuts
    pass_pt_eta_cuts = IDJet.make_pt_eta_cuts(jets)
    if cutflow is not None: cutflow['jets pass pT and eta cuts'] += ak.sum(pass_pt_eta_cuts)
    
    #    ## HEM issue
    #if year == '2018':
    #    hem_region = IDJet.HEM_15_16_issue(jets)
    
        ## tightID
    # check https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD for definitions
    if year == '2016':
        #jetId = 7 # pass loose, tight, tightLepVeto ID
        jetId = 3 # pass loose and tight, but not tightLepVeto ID
    else:
        #jetId = 4 # pass tight and tightLepVeto ID
        jetId = 2 # pass tight but not tightLepVeto ID
    jet_ID = (jets['jetId'] >= jetId) # pass at least tight
    if cutflow is not None: cutflow['jets pass ID'] += ak.sum(jet_ID)
    #jet_ID = (jets.Id >= jetId) # pass at least tight
    #if cutflow is not None: cutflow['jets pass ID'] += jet_ID.sum().sum()
    
        ## remove jets that don't pass ID and pt/eta cuts
    jets = jets[(jet_ID & pass_pt_eta_cuts)]

        ## clean jets wrt veto+loose+tight el and mu
    jets_ak = ak.with_name(jets[["pt", "eta", "phi", "mass"]],"PtEtaPhiMLorentzVector")
    muons_ak = ak.with_name(muons[(muons['TIGHTMU'] | muons['LOOSEMU'] | muons['VETOMU']) == True][["pt", "eta", "phi", "mass"]],"PtEtaPhiMLorentzVector")
    electrons_ak = ak.with_name(electrons[(electrons['TIGHTEL'] | electrons['LOOSEEL'] | electrons['VETOEL']) == True][["pt", "eta", "phi", "mass"]],"PtEtaPhiMLorentzVector")

            # all mu
    jets_akc, muons_akc = ak.unzip(ak.cartesian([jets_ak, muons_ak], nested=True))
    clean_j_Mu_mask = ak.all((jets_akc.delta_r(muons_akc) >= 0.4), axis=2)
            # all el
    jets_akc, electrons_akc = ak.unzip(ak.cartesian([jets_ak, electrons_ak], nested=True))
    clean_j_El_mask = ak.all((jets_akc.delta_r(electrons_akc) >= 0.4), axis=2)
    clean_jets_mask = (clean_j_Mu_mask & clean_j_El_mask)
    jets = jets[clean_jets_mask]

        ## leading jet pt cut
    leadpt_cut = IDJet.make_leadjet_pt_cut(jets)
    if cutflow is not None: cutflow['jets pass lead jet pT cut'] += ak.sum(leadpt_cut)
    
        ## 3 or more jets
    njet_restriction = 3
    njets_cuts = (ak.num(jets) >= njet_restriction)
    if cutflow is not None: cutflow['nEvts with %s+ clean jets passing ID and kin selection' % njet_restriction] += ak.sum(njets_cuts)
    
    passing_jets = (leadpt_cut & njets_cuts)
    if cutflow is not None: cutflow['nEvts with jets passing selection'] += ak.sum(passing_jets)

    if cutflow is not None:
        return jets, passing_jets, cutflow
    else:
        return jets, passing_jets


def select_leptons(events, year, noIso=False, cutflow=None):

    evt_sel = PackedSelection()

    ### trigger selection
    evt_sel.add('mu_triggers', Filters_and_Triggers.get_triggers(HLT=events['HLT'], leptype='Muon', year=year, noIso=noIso))
    evt_sel.add('el_triggers', Filters_and_Triggers.get_triggers(HLT=events['HLT'], leptype='Electron', year=year, noIso=noIso))

    ### filter selection
    evt_sel.add('pass_filters', Filters_and_Triggers.get_filters(Flags=events['Flag'], year=year))

    ### lepton selection
    events['Muon'] = IDMuon.process_muons(events['Muon'], year)
    events['Electron'] = IDElectron.process_electrons(events['Electron'], year)
    evt_sel.add('single_lep', np.logical_xor(( (ak.sum(events['Muon']['TIGHTMU'], axis=1) + ak.sum(events['Muon']['LOOSEMU'], axis=1)) == 1 ), ( (ak.sum(events['Electron']['TIGHTEL'], axis=1) + ak.sum(events['Electron']['LOOSEEL'], axis=1)) == 1 ))) # only single LOOSE or TIGHT el or mu (not counting vetos)
    evt_sel.add('single_looseMu', ak.sum(events['Muon']['LOOSEMU'], axis=1) == 1)
    evt_sel.add('single_tightMu', ak.sum(events['Muon']['TIGHTMU'], axis=1) == 1)
    evt_sel.add('single_looseEl', ak.sum(events['Electron']['LOOSEEL'], axis=1) == 1)
    evt_sel.add('single_tightEl', ak.sum(events['Electron']['TIGHTEL'], axis=1) == 1)
    evt_sel.add('zero_vetoMu', ak.sum(events['Muon']['VETOMU'], axis=1) == 0)
    evt_sel.add('zero_vetoEl', ak.sum(events['Electron']['VETOEL'], axis=1) == 0)
    evt_sel.add('tightMu_pass_vetoMu', (ak.sum(((events['Muon']['TIGHTMU']*1 + events['Muon']['VETOMU']*1) > 0), axis=1) == 1)) # only veto mu (if it exists) is tight
    evt_sel.add('looseMu_pass_vetoMu', (ak.sum(((events['Muon']['LOOSEMU']*1 + events['Muon']['VETOMU']*1) > 0), axis=1) == 1)) # only veto mu (if it exists) is loose
    evt_sel.add('tightEl_pass_vetoEl', (ak.sum(((events['Electron']['TIGHTEL']*1 + events['Electron']['VETOEL']*1) > 0), axis=1) == 1)) # only veto el (if it exists) is tight
    evt_sel.add('looseEl_pass_vetoEl', (ak.sum(((events['Electron']['LOOSEEL']*1 + events['Electron']['VETOEL']*1) > 0), axis=1) == 1)) # only veto el (if it exists) is loose

        # leptons pass loose/tight requirements and triggers (mu only pass mu triggers, el only pass el triggers)
    evt_sel.add('tightMu_pass', evt_sel.require(single_lep=True, single_tightMu=True, zero_vetoEl=True, tightMu_pass_vetoMu=True, mu_triggers=True) & ~evt_sel.require(el_triggers=True))
    evt_sel.add('looseMu_pass', evt_sel.require(single_lep=True, single_looseMu=True, zero_vetoEl=True, looseMu_pass_vetoMu=True, mu_triggers=True) & ~evt_sel.require(el_triggers=True))
    evt_sel.add('tightEl_pass', evt_sel.require(single_lep=True, single_tightEl=True, zero_vetoMu=True, tightEl_pass_vetoEl=True, el_triggers=True) & ~evt_sel.require(mu_triggers=True))
    evt_sel.add('looseEl_pass', evt_sel.require(single_lep=True, single_looseEl=True, zero_vetoMu=True, looseEl_pass_vetoEl=True, el_triggers=True) & ~evt_sel.require(mu_triggers=True))

    evt_sel.add('passing_mu', evt_sel.require(tightMu_pass=True) | evt_sel.require(looseMu_pass=True))
    evt_sel.add('passing_el', evt_sel.require(tightEl_pass=True) | evt_sel.require(looseEl_pass=True))
    evt_sel.add('passing_lep', evt_sel.require(passing_mu=True) | evt_sel.require(passing_el=True))

    evt_sel.add('lep_and_filter_pass', evt_sel.require(passing_lep=True, pass_filters=True))
    if cutflow is not None: cutflow['lep_and_filter_pass'] += evt_sel.require(lep_and_filter_pass=True).sum()

    return evt_sel.require(lep_and_filter_pass=True)

def jets_selection(events, year, corrections, cutflow=None, shift=None):
        ## build corrected jets and MET
    events['Jet'], events['MET'] = IDJet.process_jets(events, year, corrections['JetCor'])
            # get jet selection for systematic shift (can only support one at a time)
    if shift == 'JES_UP':
        jets_to_use = events['Jet']['JES_jes']['up']
    elif shift == 'JES_DW':
        jets_to_use = events['Jet']['JES_jes']['down']
    elif shift == 'JER_UP':
        jets_to_use = events['Jet']['JER']['up']
    elif shift == 'JER_DW':
        jets_to_use = events['Jet']['JER']['down']
    elif shift == 'MET_UP':
        raise ValueError("MET systematic shifts not supported yet!")
    elif shift == 'MET_DW':
        raise ValueError("MET systematic shifts not supported yet!")
    else:
        jets_to_use = events['Jet']

        # evaluate selection on jets
    if cutflow is not None:
        new_jets, passing_jets, cutflow = select_jets(jets_to_use, events['Muon'], events['Electron'], year, cutflow)
    else:
        new_jets, passing_jets = select_jets(jets_to_use, events['Muon'], events['Electron'], year)

        ## substitute jets for ones that pass requirements from select_jets
    events['SelectedJets'] = new_jets

    return passing_jets


def select(events, year, corrections, noIso=False, cutflow=None, shift=None):

    evt_sel = PackedSelection()

    ### trigger selection
    evt_sel.add('mu_triggers', Filters_and_Triggers.get_triggers(HLT=events['HLT'], leptype='Muon', year=year, noIso=noIso))
    evt_sel.add('el_triggers', Filters_and_Triggers.get_triggers(HLT=events['HLT'], leptype='Electron', year=year, noIso=noIso))

    ### filter selection
    evt_sel.add('pass_filters', Filters_and_Triggers.get_filters(Flags=events['Flag'], year=year))

    ### lepton selection
    events['Muon'] = IDMuon.process_muons(events['Muon'], year)
    events['Electron'] = IDElectron.process_electrons(events['Electron'], year)
    evt_sel.add('single_lep', np.logical_xor(( (ak.sum(events['Muon']['TIGHTMU'], axis=1) + ak.sum(events['Muon']['LOOSEMU'], axis=1)) == 1 ), ( (ak.sum(events['Electron']['TIGHTEL'], axis=1) + ak.sum(events['Electron']['LOOSEEL'], axis=1)) == 1 ))) # only single LOOSE or TIGHT el or mu (not counting vetos)
    evt_sel.add('single_looseMu', ak.sum(events['Muon']['LOOSEMU'], axis=1) == 1)
    evt_sel.add('single_tightMu', ak.sum(events['Muon']['TIGHTMU'], axis=1) == 1)
    evt_sel.add('single_looseEl', ak.sum(events['Electron']['LOOSEEL'], axis=1) == 1)
    evt_sel.add('single_tightEl', ak.sum(events['Electron']['TIGHTEL'], axis=1) == 1)
    evt_sel.add('zero_vetoMu', ak.sum(events['Muon']['VETOMU'], axis=1) == 0)
    evt_sel.add('zero_vetoEl', ak.sum(events['Electron']['VETOEL'], axis=1) == 0)
    evt_sel.add('tightMu_pass_vetoMu', (ak.sum(((events['Muon']['TIGHTMU']*1 + events['Muon']['VETOMU']*1) > 0), axis=1) == 1)) # only veto mu (if it exists) is tight
    evt_sel.add('looseMu_pass_vetoMu', (ak.sum(((events['Muon']['LOOSEMU']*1 + events['Muon']['VETOMU']*1) > 0), axis=1) == 1)) # only veto mu (if it exists) is loose
    evt_sel.add('tightEl_pass_vetoEl', (ak.sum(((events['Electron']['TIGHTEL']*1 + events['Electron']['VETOEL']*1) > 0), axis=1) == 1)) # only veto el (if it exists) is tight
    evt_sel.add('looseEl_pass_vetoEl', (ak.sum(((events['Electron']['LOOSEEL']*1 + events['Electron']['VETOEL']*1) > 0), axis=1) == 1)) # only veto el (if it exists) is loose

        # leptons pass loose/tight requirements and triggers (mu only pass mu triggers, el only pass el triggers)
    evt_sel.add('tightMu_pass', evt_sel.require(single_lep=True, single_tightMu=True, zero_vetoEl=True, tightMu_pass_vetoMu=True, mu_triggers=True) & ~evt_sel.require(el_triggers=True))
    evt_sel.add('looseMu_pass', evt_sel.require(single_lep=True, single_looseMu=True, zero_vetoEl=True, looseMu_pass_vetoMu=True, mu_triggers=True) & ~evt_sel.require(el_triggers=True))
    evt_sel.add('tightEl_pass', evt_sel.require(single_lep=True, single_tightEl=True, zero_vetoMu=True, tightEl_pass_vetoEl=True, el_triggers=True) & ~evt_sel.require(mu_triggers=True))
    evt_sel.add('looseEl_pass', evt_sel.require(single_lep=True, single_looseEl=True, zero_vetoMu=True, looseEl_pass_vetoEl=True, el_triggers=True) & ~evt_sel.require(mu_triggers=True))

    evt_sel.add('passing_mu', evt_sel.require(tightMu_pass=True) | evt_sel.require(looseMu_pass=True))
    evt_sel.add('passing_el', evt_sel.require(tightEl_pass=True) | evt_sel.require(looseEl_pass=True))
    evt_sel.add('passing_lep', evt_sel.require(passing_mu=True) | evt_sel.require(passing_el=True))

    evt_sel.add('lep_and_filter_pass', evt_sel.require(passing_lep=True, pass_filters=True))
    if cutflow is not None: cutflow['lep_and_filter_pass'] += evt_sel.require(lep_and_filter_pass=True).sum()

        ## build corrected jets and MET
    events['Jet'], events['MET'] = IDJet.process_jets(events, year, corrections['JetCor'])
            # get jet selection for systematic shift (can only support one at a time)
    if shift == 'JES_UP':
        jets_to_use = events['Jet']['JES_jes']['up']
    elif shift == 'JES_DW':
        jets_to_use = events['Jet']['JES_jes']['down']
    elif shift == 'JER_UP':
        jets_to_use = events['Jet']['JER']['up']
    elif shift == 'JER_DW':
        jets_to_use = events['Jet']['JER']['down']
    elif shift == 'MET_UP':
        raise ValueError("MET systematic shifts not supported yet!")
    elif shift == 'MET_DW':
        raise ValueError("MET systematic shifts not supported yet!")
    else:
        jets_to_use = events['Jet']

        # evaluate selection on jets
    if cutflow is not None:
        new_jets, passing_jets, cutflow = select_jets(jets_to_use, events['Muon'], events['Electron'], year, cutflow)
    else:
        new_jets, passing_jets = select_jets(jets_to_use, events['Muon'], events['Electron'], year)

        ## substitute jets for ones that pass requirements from select_jets
    events['Jet'] = new_jets

    passing_evts = evt_sel.require(lep_and_filter_pass=True) & passing_jets

    return passing_evts
