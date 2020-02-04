from pdb import set_trace
import python.triggers as triggers
import coffea.processor as processor
import coffea.nanoaod.nanoevents
import coffea.processor.dataframe

def select_muons(df, accumulator=None):
    if isinstance(df, coffea.nanoaod.nanoevents.NanoEvents):
        if accumulator: accumulator['cutflow']['before mu sel'] += df['Muon'].size

            ## pass muon trigger
        mu_trig_mask = triggers.mu_triggers(df)
        if accumulator: accumulator['cutflow']['nEvts pass mu triggers'] += mu_trig_mask.sum()

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
        passing_mus = (one_mu_tight_mask) & (mu_trig_mask)
        if accumulator: accumulator['cutflow']['passing muons'] += passing_mus.sum()

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

            ## pass electron trigger
        el_trig_mask = triggers.el_triggers(df)
        if accumulator: accumulator['cutflow']['nEvts pass el triggers'] += el_trig_mask.sum()

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
