from pdb import set_trace
import coffea.nanoaod.nanoevents
import coffea.processor.dataframe

def select_muons(df, accumulator=None):
    if isinstance(df, coffea.nanoaod.nanoevents.NanoEvents):
        if accumulator: accumulator['cutflow']['before mu sel'] += df['Muon'].size

            ## tight muons
        tight_muons_mask = (df['Muon'].tightId)
        if accumulator: accumulator['cutflow']['n tight muons'] += tight_muons_mask.any().sum()

            ## single muon
        single_muon_mask = (df['Muon'].counts == 1)
        if accumulator: accumulator['cutflow']['1 muon'] += single_muon_mask.sum()

            ## single tight muon
        one_mu_tight_mask = (single_muon_mask & tight_muons_mask).any()
        if accumulator: accumulator['cutflow']['1 muon, tightID'] += one_mu_tight_mask.sum()

            ## pass all muon criteria
        passing_mus = (one_mu_tight_mask)

    #elif isinstance(df, coffea.processor.dataframe.LazyDataFrame):
    #    from python.IDMuon import process_muons as proc_mus
    #    muons = proc_mus(df) # get muons

    #        # pass triggers
    #    trig_muons = triggers.mu_triggers(df)
    #    muons = muons[trig_muons]

    #        # single muon    
    #    onemuon = (muons.counts == 1)
    #    muons = muons[onemuon]

    #        # tight muon    
    #    tight_mu = (muons.tightId > 0)
    #    muons = muons[tight_mu]

    else:
        raise ValueError("Only NanoEvents and LazyDataFrame formats supported right now")

    #set_trace()
    if accumulator:
        return passing_mus, accumulator
    else:
        return passing_mus


def select_jets(df, accumulator=None):
    import python.IDJet as IDJet
    
    ## hardcoded for now
    btagger = 'DEEPCSV'
    btag_wps = ['MEDIUM', 'MEDIUM']
    
    #set_trace()
    if isinstance(df, coffea.nanoaod.nanoevents.NanoEvents):
        df = IDJet.add_btag_disc(df, btagger=btagger, tightb=btag_wps[0], looseb=btag_wps[1])

            ## only 4 jets
        four_jets = (df['Jet'].counts == 4)
        if accumulator: accumulator['cutflow']['4 jets'] += four_jets.sum()

            #btag reqs
        if len(list(set(btag_wps))) == 1:
            btag_pass = (df['Jet'][btagger+btag_wps[0]]).sum() >= 2
            if accumulator: accumulator['cutflow']['>=2 jets pass %s' % btagger+btag_wps[0]] += btag_pass.sum()
        else:
            raise IOError("Only 1 unique btag working point supported now")

        passing_jets = (four_jets) & (btag_pass)
        if accumulator: accumulator['cutflow']['pass btag + nJets'] += passing_jets.sum()

    #elif isinstance(df, coffea.processor.dataframe.LazyDataFrame):
    #    from python.IDJet import process_jets as proc_jets
    #    jets = proc_jets(df) # get jets

    else:
        raise ValueError("Only NanoEvents and LazyDataFrame formats supported right now")

    if accumulator:
        return passing_jets, accumulator
    else:
        return passing_jets


def select_electrons(df, accumulator=None):
    if isinstance(df, coffea.nanoaod.nanoevents.NanoEvents):
        if accumulator: accumulator['cutflow']['before el sel'] += df['Electron'].size

        set_trace()
            ## tight electrons
        tight_electrons_mask = (df['Electron'].tightId)
        if accumulator: accumulator['cutflow']['n tight electrons'] += tight_electrons_mask.any().sum()

            ## single electron
        single_electron_mask = (df['Electron'].counts == 1)
        if accumulator: accumulator['cutflow']['1 electron'] += single_electron_mask.sum()

            ## single tight electron
        one_el_tight_mask = (single_electron_mask & tight_electrons_mask).any()
        if accumulator: accumulator['cutflow']['1 electron, tightID'] += one_el_tight_mask.sum()

            ## pass all electron criteria
        passing_els = (one_el_tight_mask) & (el_trig_mask)

#    #elif isinstance(df, coffea.processor.dataframe.LazyDataFrame):
#    #    from python.IDElectron import process_electrons as proc_els
#
    else:
        raise ValueError("Only NanoEvents and LazyDataFrame formats supported right now")

    if accumulator:
        return passing_els, accumulator
    else:
        return passing_els


def get_triggers(df, leptype, accumulator=None):
    ## event triggers to be used found here: https://twiki.cern.ch/twiki/bin/view/CMS/TopTriggerYear2016 or 2017, 2018...
    if isinstance(df, coffea.nanoaod.nanoevents.NanoEvents):
        if leptype == 'Muon':
            trigger = (df['HLT']['IsoMu24']) | (df['HLT']['IsoTkMu24'])
        elif leptype == 'Electron':
            trigger = (df['HLT']['Ele27_WPTight_Gsf'])
        else:
            raise ValueError("Only events analyzing muons OR electrons supported right now")
    
    elif isinstance(df, coffea.processor.dataframe.LazyDataFrame):
        if leptype == 'Muon':
            trigger = (df['HLT_IsoMu24'] > 0) | (df['HLT_IsoTkMu24'] > 0)
        elif leptype == 'Electron':
            trigger = df['HLT_Ele27_WPTight_Gsf']
        else:
            raise ValueError("Only events analyzing muons OR electrons supported right now")

    else:
        raise ValueError("Only NanoEvents and LazyDataFrame formats supported right now")

    if accumulator:
        accumulator['cutflow']['nEvts pass %s triggers' % leptype] += trigger.sum()
        return trigger, accumulator
    else:
        return trigger


def get_filters(df, accumulator=None):
    ## Supported filters found here: https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
    if isinstance(df, coffea.nanoaod.nanoevents.NanoEvents):
        pass_filters = (df['Flag']['goodVertices']) & (df['Flag']['globalSuperTightHalo2016Filter']) & (df['Flag']['HBHENoiseFilter']) & (df['Flag']['HBHENoiseIsoFilter']) & (df['Flag']['EcalDeadCellTriggerPrimitiveFilter']) & (df['Flag']['BadPFMuonFilter'])
    #elif isinstance(df, coffea.processor.dataframe.LazyDataFrame):
        #pass_filters = (df['Flag']['goodVertices']) & (df['Flag']['globalSuperTightHalo2016Filter']) & (df['Flag']['HBHENoiseFilter']) & (df['Flag']['HBHENoiseIsoFilter']) & (df['Flag']['EcalDeadCellTriggerPrimitiveFilter']) & (df['Flag']['BadPFMuonFilter'])
    else:
        raise ValueError("Only NanoEvents and LazyDataFrame formats supported right now")

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
            pass_triggers, accumulator = get_triggers(df, leptype, accumulator)
        else:
            pass_triggers = get_triggers(df, leptype)

            ## get filters
        if accumulator:
            pass_filters, accumulator = get_filters(df, accumulator)
        else:
            pass_filters = get_filters(df)

        ### lepton selection
        if leptype == 'Muon':
            if accumulator:
                passing_leps, accumulator = select_muons(df, accumulator)
            else:
                passing_leps = select_muons(df)

        else:
        #elif leptype == 'Electron':
            if accumulator:
                passing_leps, accumulator = select_electrons(df, accumulator)
            else:
                passing_leps = select_electrons(df)

        if accumulator:
            accumulator['cutflow']['passing %s' % leptype] += (pass_triggers & passing_leps).sum()

        ### jets selection
        if accumulator:
            passing_jets, accumulator = select_jets(df, accumulator)
        else:
            passing_jets = select_jets(df)

    passing_evts = passing_jets & passing_leps & pass_triggers & pass_filters

    return passing_evts
