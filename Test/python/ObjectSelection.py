from pdb import set_trace
import coffea.nanoaod.nanoevents
import coffea.processor.dataframe
import awkward

def select_muons(muons, accumulator=None):
    import python.IDMuon as IDMuon
    
    muons = IDMuon.make_muon_ids(muons)
    #set_trace()

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
        raise ValueError("Only AwkwardArrays from NanoEvents are supported right now")

    if accumulator:
        return passing_mus, accumulator
    else:
        return passing_mus


def select_jets(jets, accumulator=None):
    import python.IDJet as IDJet
    
    ## hardcoded for now
    btagger = 'DEEPCSV'
    btag_wps = ['MEDIUM', 'MEDIUM']
    
    #set_trace()
    if isinstance(jets, awkward.array.base.AwkwardArray):
        jets = IDJet.add_btag_disc(jets, btagger=btagger, tightb=btag_wps[0], looseb=btag_wps[1])

            ## only 4 jets
        four_jets = (jets.counts == 4)
        if accumulator: accumulator['cutflow']['4 jets'] += four_jets.sum()

            #btag reqs
        if len(list(set(btag_wps))) == 1:
            btag_pass = (jets[btagger+btag_wps[0]]).sum() >= 2
            if accumulator: accumulator['cutflow']['>=2 jets pass %s' % btagger+btag_wps[0]] += btag_pass.sum()
        else:
            raise IOError("Only 1 unique btag working point supported now")

        passing_jets = (four_jets) & (btag_pass)
        if accumulator: accumulator['cutflow']['pass btag + nJets'] += passing_jets.sum()

    else:
        raise ValueError("Only AwkwardArrays from NanoEvents are supported right now")

    if accumulator:
        return passing_jets, accumulator
    else:
        return passing_jets


def select_electrons(electrons, accumulator=None):
    import python.IDElectron as IDElectron

    electrons = IDElectron.make_etaSC(electrons)
    electrons = IDElectron.make_electron_ids(electrons)

    if isinstance(electrons, awkward.array.base.AwkwardArray):
        if accumulator: accumulator['cutflow']['before el sel'] += electrons.size

        #set_trace()
            ## tight electrons
        tight_electrons_mask = (electrons.mvaFall17V2noIso_WP80)
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


def get_triggers(triggers, leptype, accumulator=None):
    ## event triggers to be used found here: https://twiki.cern.ch/twiki/bin/view/CMS/TopTriggerYear2016 or 2017, 2018...
    if isinstance(triggers, awkward.array.base.AwkwardArray):
        if leptype == 'Muon':
            trigger = (triggers.IsoMu24) | (triggers.IsoTkMu24)
        elif leptype == 'Electron':
            trigger = (triggers.Ele27_WPTight_Gsf)
        else:
            raise ValueError("Only events analyzing muons OR electrons supported right now")

    else:
        raise ValueError("Only AwkwardArrays from NanoEvents are supported right now")

    if accumulator:
        accumulator['cutflow']['nEvts pass %s triggers' % leptype] += trigger.sum()
        return trigger, accumulator
    else:
        return trigger


def get_filters(filters, accumulator=None):
    ## Supported filters found here: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
    if isinstance(filters, awkward.array.base.AwkwardArray):
        pass_filters = (filters.goodVertices) & (filters.globalSuperTightHalo2016Filter) & (filters.HBHENoiseFilter) & (filters.HBHENoiseIsoFilter) & (filters.EcalDeadCellTriggerPrimitiveFilter) & (filters.BadPFMuonFilter)

    else:
        raise ValueError("Only AwkwardArrays from NanoEvents are supported right now")

    if accumulator:
        accumulator['cutflow']['pass filters'] += pass_filters.sum()
        return pass_filters, accumulator
    else:
        return pass_filters


def select(df, leptype, accumulator=None, shift=None):

    #set_trace()
    if leptype != 'Muon' and leptype != 'Electron':
        raise IOError("Only events analyzing muons OR electrons supported right now")

    if isinstance(df, coffea.nanoaod.nanoevents.NanoEvents):
            ## get triggers
        if accumulator:
            pass_triggers, accumulator = get_triggers(df['HLT'], leptype, accumulator)
        else:
            pass_triggers = get_triggers(df['HLT'], leptype)

            ## get filters
        if accumulator:
            pass_filters, accumulator = get_filters(df['Flag'], accumulator)
        else:
            pass_filters = get_filters(df['Flag'])

        ### lepton selection
        if leptype == 'Muon':
            if accumulator:
                passing_leps, accumulator = select_muons(df['Muon'], accumulator)
            else:
                passing_leps = select_muons(df['Muon'])

        else:
        #elif leptype == 'Electron':
            if accumulator:
                passing_leps, accumulator = select_electrons(df['Electron'], accumulator)
            else:
                passing_leps = select_electrons(df['Electron'])

        if accumulator:
            accumulator['cutflow']['passing %s' % leptype] += (pass_triggers & passing_leps).sum()

        ### jets selection
        if accumulator:
            passing_jets, accumulator = select_jets(df['Jet'], accumulator)
        else:
            passing_jets = select_jets(df['Jet'])

    passing_evts = passing_jets & passing_leps & pass_triggers & pass_filters

    return passing_evts
