from pdb import set_trace
import coffea.nanoaod.nanoevents
import coffea.processor.dataframe
import awkward

def select_muons(muons, accumulator=None):
    import python.IDMuon as IDMuon
    
    if isinstance(muons, awkward.array.base.AwkwardArray):
        muons = IDMuon.make_muon_ids(muons)
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
    
    if isinstance(jets, awkward.array.base.AwkwardArray):
        jets = IDJet.build_jets(jets)

            ## kinematic cuts (pt, leadpt, eta)
        kin_cut_jets = jets.kin_cuts
        if accumulator: accumulator['cutflow']['jets pass kin cuts'] += kin_cut_jets.sum().sum()

            ## only 4 jets
        four_jets = (jets.counts == 4)
        if accumulator: accumulator['cutflow']['nEvts with 4 jets'] += four_jets.sum()

            #btag reqs
        btag_wps = [col for col in jets.columns if 'BTAG_' in col]
        if len(list(set(btag_wps))) == 1:
            btag_pass = (jets[btag_wps[0]]).sum() >= 2
            if accumulator: accumulator['cutflow']['nEvts >=2 jets pass %s' % btag_wps[0]] += btag_pass.sum()
        else:
            raise IOError("Only 1 unique btag working point supported now")

        passing_jets = (four_jets & btag_pass & kin_cut_jets).any()
        if accumulator: accumulator['cutflow']['nEvts pass btag + nJets + kin cuts'] += passing_jets.sum()

    else:
        raise ValueError("Only AwkwardArrays from NanoEvents are supported right now")

    if accumulator:
        return passing_jets, accumulator
    else:
        return passing_jets


def select_electrons(electrons, accumulator=None):
    import python.IDElectron as IDElectron

    if isinstance(electrons, awkward.array.base.AwkwardArray):
        electrons = IDElectron.build_electrons(electrons)
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

    #set_trace()
    passing_evts = passing_jets & passing_leps & pass_triggers & pass_filters

    return passing_evts
