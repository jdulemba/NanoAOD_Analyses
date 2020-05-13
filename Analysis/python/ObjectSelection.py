from pdb import set_trace
import coffea.processor.dataframe
import awkward
import python.Filters_and_Triggers as Filters_and_Triggers
import python.IDJet as IDJet
import python.IDMuon as IDMuon
import python.IDElectron as IDElectron
import python.IDMet as IDMet
import numpy as np
import coffea.processor as processor

def select_jets(jets, muons, electrons, year, accumulator=None):
    
    if isinstance(jets, awkward.array.base.AwkwardArray):
            ## pt and eta cuts
        pass_pt_eta_cuts = IDJet.make_pt_eta_cuts(jets)
        if accumulator: accumulator['cutflow']['jets pass pT and eta cuts'] += pass_pt_eta_cuts.sum().sum()

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
        jet_ID = (jets.Id >= jetId) # pass at least tight
        #jet_ID = (jets.Id == jetId)
        if accumulator: accumulator['cutflow']['jets pass ID'] += jet_ID.sum().sum()

            ## remove jets that don't pass ID and pt/eta cuts
        jets = jets[(jet_ID & pass_pt_eta_cuts)]

            ## clean jets wrt loose+tight el and mu
        clean_j_tightMu = ~(jets.match(muons[muons['TIGHTMU'] == True], deltaRCut=0.4))
        clean_j_looseMu = ~(jets.match(muons[muons['LOOSEMU'] == True], deltaRCut=0.4))
        clean_j_tightEl = ~(jets.match(electrons[electrons['TIGHTEL'] == True], deltaRCut=0.4))
        clean_j_looseEl = ~(jets.match(electrons[electrons['LOOSEEL'] == True], deltaRCut=0.4))
            ## clean jets wrt veto el and mu
        clean_j_vetoMu = ~(jets.match(muons[muons['VETOMU'] == True], deltaRCut=0.4))
        clean_j_vetoEl = ~(jets.match(electrons[electrons['VETOEL'] == True], deltaRCut=0.4))
        jets = jets[(clean_j_looseEl & clean_j_tightEl & clean_j_looseMu & clean_j_tightMu) & (clean_j_vetoMu & clean_j_vetoEl)]

            ## leading jet pt cut
        leadpt_cut = IDJet.make_leadjet_pt_cut(jets)
        if accumulator: accumulator['cutflow']['jets pass lead jet pT cut'] += leadpt_cut.sum()

            ## 3 or more jets
        njet_restriction = 3
        njets_cuts = (jets.counts >= njet_restriction)
        #njets_cuts = (jets.counts == njet_restriction)
        if accumulator: accumulator['cutflow']['nEvts with %s+ clean jets passing ID and kin selection' % njet_restriction] += njets_cuts.sum()

        passing_jets = (leadpt_cut & njets_cuts)

    else:
        raise ValueError("Only AwkwardArrays are supported")

    if accumulator:
        return jets, passing_jets, accumulator
    else:
        return jets, passing_jets


def select(df, year, corrections, accumulator=None, shift=None):

    if not isinstance(df, coffea.processor.dataframe.LazyDataFrame):
        raise IOError("This function only works for LazyDataFrame objects")
    evt_sel = processor.PackedSelection()

    ### trigger selection
    evt_sel.add('mu_triggers', Filters_and_Triggers.get_triggers(df=df, leptype='Muon', year=year))
    evt_sel.add('el_triggers', Filters_and_Triggers.get_triggers(df=df, leptype='Electron', year=year))

    ### filter selection
    evt_sel.add('pass_filters', Filters_and_Triggers.get_filters(df=df, year=year))

    ### lepton selection
    df['Muon'] = IDMuon.process_muons(df, year)
    df['Electron'] = IDElectron.process_electrons(df, year)
    evt_sel.add('single_lep', np.logical_xor(( (df['Muon']['TIGHTMU'].sum() + df['Muon']['LOOSEMU'].sum()) == 1 ), ( (df['Electron']['TIGHTEL'].sum() + df['Electron']['LOOSEEL'].sum()) == 1 ))) # only single LOOSE or TIGHT el or mu (not counting vetos)
    evt_sel.add('single_looseMu', df['Muon']['LOOSEMU'].sum() == 1)
    evt_sel.add('single_tightMu', df['Muon']['TIGHTMU'].sum() == 1)
    evt_sel.add('single_looseEl', df['Electron']['LOOSEEL'].sum() == 1)
    evt_sel.add('single_tightEl', df['Electron']['TIGHTEL'].sum() == 1)
    evt_sel.add('zero_vetoMu', df['Muon']['VETOMU'].sum() == 0)
    evt_sel.add('zero_vetoEl', df['Electron']['VETOEL'].sum() == 0)
    evt_sel.add('tightMu_pass_vetoMu', (((df['Muon']['TIGHTMU']*1 + df['Muon']['VETOMU']*1) > 0).sum() == 1)) # only veto mu (if it exists) is tight
    evt_sel.add('looseMu_pass_vetoMu', (((df['Muon']['LOOSEMU']*1 + df['Muon']['VETOMU']*1) > 0).sum() == 1)) # only veto mu (if it exists) is loose
    evt_sel.add('tightEl_pass_vetoEl', (((df['Electron']['TIGHTEL']*1 + df['Electron']['VETOEL']*1) > 0).sum() == 1)) # only veto el (if it exists) is tight
    evt_sel.add('looseEl_pass_vetoEl', (((df['Electron']['LOOSEEL']*1 + df['Electron']['VETOEL']*1) > 0).sum() == 1)) # only veto el (if it exists) is loose

        # leptons pass loose/tight requirements and triggers (mu only pass mu triggers, el only pass el triggers)
    evt_sel.add('tightMu_pass', evt_sel.require(single_lep=True, single_tightMu=True, zero_vetoEl=True, tightMu_pass_vetoMu=True, mu_triggers=True) & ~evt_sel.require(el_triggers=True))
    evt_sel.add('looseMu_pass', evt_sel.require(single_lep=True, single_looseMu=True, zero_vetoEl=True, looseMu_pass_vetoMu=True, mu_triggers=True) & ~evt_sel.require(el_triggers=True))
    evt_sel.add('tightEl_pass', evt_sel.require(single_lep=True, single_tightEl=True, zero_vetoMu=True, tightEl_pass_vetoEl=True, el_triggers=True) & ~evt_sel.require(mu_triggers=True))
    evt_sel.add('looseEl_pass', evt_sel.require(single_lep=True, single_looseEl=True, zero_vetoMu=True, looseEl_pass_vetoEl=True, el_triggers=True) & ~evt_sel.require(mu_triggers=True))

    evt_sel.add('passing_mu', evt_sel.require(tightMu_pass=True) | evt_sel.require(looseMu_pass=True))
    evt_sel.add('passing_el', evt_sel.require(tightEl_pass=True) | evt_sel.require(looseEl_pass=True))
    evt_sel.add('passing_lep', evt_sel.require(passing_mu=True) | evt_sel.require(passing_el=True))

    #set_trace()
    evt_sel.add('lep_and_filter_pass', evt_sel.require(passing_lep=True, pass_filters=True))

        ## SetMET
    df['MET'] = IDMet.process_met(df)

    ### jets selection
    df['Jet'] = IDJet.process_jets(df, year, corrections['JetCor']) # initialize jets
    if accumulator:
        new_jets, passing_jets, accumulator = select_jets(df['Jet'], df['Muon'], df['Electron'], year, accumulator)
    else:
        new_jets, passing_jets = select_jets(df['Jet'], df['Muon'], df['Electron'], year)
    
        ## substitute jets for ones that pass requirements from select_jets
    df['Jet'] = new_jets

    passing_evts = evt_sel.require(lep_and_filter_pass=True) & passing_jets

    return passing_evts
