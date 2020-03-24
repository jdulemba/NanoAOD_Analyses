from pdb import set_trace
import coffea.nanoaod.nanoevents
import coffea.processor.dataframe
import awkward
import python.Filters_and_Triggers as Filters_and_Triggers
import python.IDJet as IDJet
import python.IDMuon as IDMuon
import python.IDElectron as IDElectron
import python.IDMet as IDMet


def select_muons(muons, accumulator=None):
    
    if isinstance(muons, awkward.array.base.AwkwardArray):
        if accumulator: accumulator['cutflow']['before mu sel'] += muons.size

            ## tight muons
        tight_muons_mask = (muons.tightId)
        if accumulator: accumulator['cutflow']['n tight muons'] += tight_muons_mask.any().sum()

            ## single muon
        single_muon_mask = (muons.counts == 1)
        if accumulator: accumulator['cutflow']['1 muon'] += single_muon_mask.sum()

            ## single tight muon
        one_mu_tight_mask = (single_muon_mask & tight_muons_mask).any()
        if accumulator: accumulator['cutflow']['1 muon, tightID'] += one_mu_tight_mask.sum()

            ## pass all muon criteria
        passing_mus = (one_mu_tight_mask)

    else:
        raise ValueError("Only AwkwardArrays are supported")

    if accumulator:
        return passing_mus, accumulator
    else:
        return passing_mus


def select_jets(jets, year, accumulator=None):
    
    if isinstance(jets, awkward.array.base.AwkwardArray):

            ## pt and eta cuts
        pass_pt_eta_cuts = IDJet.make_pt_eta_cuts(jets)
        if accumulator: accumulator['cutflow']['jets pass pT and eta cuts'] += pass_pt_eta_cuts.sum().sum()

            ## tightID
        # check https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD for definitions
        if year == '2016':
            jetId = 3
        else:
            jetId = 2
        jet_tightID = (jets.Id >= jetId) # 3 for passing loose+tight for 2016, other years are different https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD
        if accumulator: accumulator['cutflow']['jets tightID'] += jet_tightID.sum().sum()

            ## remove jets that don't pass tightID and pt/eta cuts
        jets = jets[(jet_tightID & pass_pt_eta_cuts)]

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


def select_electrons(electrons, accumulator=None):

    if isinstance(electrons, awkward.array.base.AwkwardArray):

        if accumulator: accumulator['cutflow']['before el sel'] += electrons.size

        #set_trace()
            ## tight electrons
        tight_electrons_mask = (electrons.tightID)
        if accumulator: accumulator['cutflow']['n tight electrons'] += tight_electrons_mask.any().sum()

            ## single electron
        single_electron_mask = (electrons.counts == 1)
        if accumulator: accumulator['cutflow']['1 electron'] += single_electron_mask.sum()

            ## single tight electron
        one_el_tight_mask = (single_electron_mask & tight_electrons_mask).any()
        if accumulator: accumulator['cutflow']['1 electron, tightID'] += one_el_tight_mask.sum()

            ## pass all electron criteria
        passing_els = (one_el_tight_mask)

    else:
        raise ValueError("Only AwkwardArrays from NanoEvents are supported right now")

    if accumulator:
        return passing_els, accumulator
    else:
        return passing_els


def select(df, leptype, year, accumulator=None, shift=None):

    #set_trace()
    if leptype != 'Muon' and leptype != 'Electron':
        raise IOError("Only events analyzing muons OR electrons supported right now")

    if not isinstance(df, coffea.processor.dataframe.LazyDataFrame):
        raise IOError("This function only works for LazyDataFrame objects")

        ## get triggers
    if accumulator:
        pass_triggers, accumulator = Filters_and_Triggers.get_triggers(df=df, leptype=leptype, year=year, accumulator=accumulator)
    else:
        pass_triggers = Filters_and_Triggers.get_triggers(df=df, leptype=leptype, year=year)

        ## get filters
    if accumulator:
        pass_filters, accumulator = Filters_and_Triggers.get_filters(df=df, year=year, accumulator=accumulator)
    else:
        pass_filters = Filters_and_Triggers.get_filters(df=df, year=year)
    #set_trace()

    ### lepton selection
    if leptype == 'Muon':
        df['Muon'] = IDMuon.process_muons(df)
        if accumulator:
            passing_leps, accumulator = select_muons(df['Muon'], accumulator)
        else:
            passing_leps = select_muons(df['Muon'])

    else:
    #elif leptype == 'Electron':
        df['Electron'] = IDElectron.process_electrons(df)
        if accumulator:
            passing_leps, accumulator = select_electrons(df['Electron'], accumulator)
        else:
            passing_leps = select_electrons(df['Electron'])

    if accumulator:
        accumulator['cutflow']['passing %s' % leptype] += (pass_triggers & pass_filters & passing_leps).sum()

    ### jets selection
    df['Jet'] = IDJet.process_jets(df, year) # initialize jets
        ## clean jets 
    df['Jet'] = df['Jet'][(~df['Jet'].match(df[leptype], deltaRCut=0.4))] ## get only clean jets based on DeltaR(jet, lepton) = 0.4
    if accumulator:
        new_jets, passing_jets, accumulator = select_jets(df['Jet'], year, accumulator)
    else:
        new_jets, passing_jets = select_jets(df['Jet'], year)
    
        ## substitute jets for ones that pass requirements from select_jets
    df['Jet'] = new_jets

    #set_trace()
    passing_evts = pass_triggers & pass_filters & passing_jets & passing_leps

        ## SetMET
    df['MET'] = IDMet.SetMET(df)

    return passing_evts
