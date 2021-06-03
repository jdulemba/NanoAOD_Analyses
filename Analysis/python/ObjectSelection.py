from pdb import set_trace
import awkward as ak
import python.Filters_and_Triggers as Filters_and_Triggers
import python.IDJet as IDJet
import python.IDMuon as IDMuon
import python.IDElectron as IDElectron
import numpy as np
from coffea.analysis_tools import PackedSelection
from functools import partial
import operator
from coffea.jetmet_tools.CorrectedJetsFactory import awkward_rewrap
from coffea.jetmet_tools.CorrectedJetsFactory import rewrap_recordarray

def hem1516_corr(jets, MET, corr, lazy_cache): # HEM region that had issues in 2018 data
    # all values based on recommendation from https://hypernews.cern.ch/HyperNews/CMS/get/JetMET/2000.html
    jets_scale = ak.to_numpy(ak.ones_like(ak.flatten(jets['pt'])))
    jets_20pc_mask = ak.to_numpy(ak.flatten((-1.57 < jets['phi']) & (jets['phi'] < -0.87) & (-2.5 < jets['eta']) & (jets['eta'] < -1.3)))
    jets_35pc_mask = ak.to_numpy(ak.flatten((-1.57 < jets['phi']) & (jets['phi'] < -0.87) & (-3.0 < jets['eta']) & (jets['eta'] < -2.5)))
    jets_scale[jets_20pc_mask] = jets_scale[jets_20pc_mask]*0.8
    jets_scale[jets_35pc_mask] = jets_scale[jets_35pc_mask]*0.65

        ## rescale jets based on HEM recommendation, using form from CorrectedJetsFactory
    out_jets = ak.flatten(jets)
    wrap = partial(awkward_rewrap, like_what=jets, gfunc=rewrap_recordarray)
    scalar_form = ak.without_parameters(
        out_jets.energy
    ).layout.form
    # Scale jet pt
    init_pt_corr = partial(
        ak.virtual,
        operator.mul,
        args=(ak.Array(jets_scale), out_jets["pt"]),
        cache=lazy_cache,
    )
    out_jets["pt"] = init_pt_corr(length=len(out_jets), form=scalar_form)
    # Scale jet mass
    init_mass_corr = partial(
        ak.virtual,
        operator.mul,
        args=(ak.Array(jets_scale), out_jets["mass"]),
        cache=lazy_cache,
    )
    out_jets["mass"] = init_mass_corr(length=len(out_jets), form=scalar_form)
    # remake jets with updated values
    new_jets = wrap(out_jets)

        ## rescale MET based on updated jet values, using form from CorrectedMETFactory
    #set_trace()
    new_met = corr['MC']['METFactory'].build(MET, new_jets, lazy_cache=lazy_cache)

    return new_jets, new_met


def remove_HEM_objs(obj, isData=None):
    in_hem_region = (-1.57 < obj['phi']) & (obj['phi'] < -0.87) & (-3.0 < obj['eta']) & (obj['eta'] < -1.3)
    remove_objs = ak.to_numpy(ak.flatten(in_hem_region))
    if isData is not None:
        runs_broadcast = ak.broadcast_arrays(isData, obj.pt)[0] # make runs same shape as pt
        failing_runs = runs_broadcast >= 319077 # runs corresponding to HEM failure
        remove_hem_objs = (in_hem_region & failing_runs)
    else:
            # remove ~60% of MC in HEM region
        np.random.seed(10)
        remove_mc = np.random.choice(2, ak.sum(in_hem_region), p=[0.4, 0.6]).astype(bool) # 60% of MC events in hem region don't pass on average
        remove_objs[remove_objs == True] = remove_mc
        remove_hem_objs = ak.unflatten(remove_objs, ak.num(obj))

    return obj[~remove_hem_objs]


def select_jets(jets, muons, electrons, year, cutflow=None):

    #set_trace()
        ## pt and eta cuts
    pass_pt_eta_cuts = IDJet.make_pt_eta_cuts(jets)
    if cutflow is not None: cutflow['jets pass pT and eta cuts'] += ak.sum(pass_pt_eta_cuts)
    
        ## tightID
    # check https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVUL for definitions
    #if year == '2016':
    #    #jetId = 7 # pass loose, tight, tightLepVeto ID
    #    jetId = 3 # pass loose and tight, but not tightLepVeto ID
    #else:
    #    #jetId = 4 # pass tight and tightLepVeto ID
    #    jetId = 2 # pass tight but not tightLepVeto ID
    #jet_ID = (jets['jetId'] >= jetId) # pass at least tight
    jet_ID = jets.isTight # pass at least tight
    if cutflow is not None: cutflow['jets pass ID'] += ak.sum(jet_ID)

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
    if ak.sum(clean_jets_mask) > 0: jets = jets[clean_jets_mask] 

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


def select_leptons(events, year, noIso=False, cutflow=None, hem_15_16=False):

    evt_sel = PackedSelection()

    ### trigger selection
    evt_sel.add('mu_triggers', Filters_and_Triggers.get_triggers(HLT=events['HLT'], leptype='Muon', year=year, noIso=noIso))
    evt_sel.add('el_triggers', Filters_and_Triggers.get_triggers(HLT=events['HLT'], leptype='Electron', year=year, noIso=noIso))

    ### filter selection
    evt_sel.add('pass_filters', Filters_and_Triggers.get_filters(Flags=events['Flag'], year=year, isMC=not events.metadata['dataset'].startswith('data_Single')))

    ### lepton selection
    events['Muon'] = IDMuon.process_muons(events['Muon'], year)
    events['Electron'] = IDElectron.process_electrons(events['Electron'], year)
    if (year == '2018') and (hem_15_16):
        events['Electron'] = remove_HEM_objs(obj=events['Electron'], isData=events.run if events.metadata['dataset'].startswith('data_Single') else None)

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

def jets_selection(events, year, cutflow=None, shift=None, hem_15_16=False):
#def jets_selection(events, year, cutflow=None, shift=None, hem_15_16=None):
        # get jet selection for systematic shift (can only support one at a time)
    if shift == 'JES_UP':
        jets_to_use = events['Jet']['JES_jes']['up']
        met_to_use = events['MET']['JES_jes']['up']
    elif shift == 'JES_DW':
        jets_to_use = events['Jet']['JES_jes']['down']
        met_to_use = events['MET']['JES_jes']['down']
    elif shift == 'JER_UP':
        jets_to_use = events['Jet']['JER']['up']
        met_to_use = events['MET']['JER']['up']
    elif shift == 'JER_DW':
        jets_to_use = events['Jet']['JER']['down']
        met_to_use = events['MET']['JER']['down']
    elif shift == 'MET_UP':
        jets_to_use = events['Jet']
        met_to_use = events['MET']['MET_UnclusteredEnergy']['up']
    elif shift == 'MET_DW':
        jets_to_use = events['Jet']
        met_to_use = events['MET']['MET_UnclusteredEnergy']['down']
    else:
        jets_to_use = events['Jet']
        met_to_use = events['MET']

    if (year == '2018') and (hem_15_16):
        #set_trace()
        jets_to_use = remove_HEM_objs(obj=jets_to_use, isData=events.run if events.metadata['dataset'].startswith('data_Single') else None)
        #jets_to_use, met_to_use = hem1516_corr(jets=jets_to_use, MET=met_to_use, corr=hem_15_16, lazy_cache=events.caches[0])
        #set_trace()

        # evaluate selection on jets
    if cutflow is not None:
        new_jets, passing_jets, cutflow = select_jets(jets_to_use, events['Muon'], events['Electron'], year, cutflow)
    else:
        new_jets, passing_jets = select_jets(jets_to_use, events['Muon'], events['Electron'], year)

        ## substitute jets for ones that pass requirements from select_jets
    events['SelectedJets'] = new_jets
    events['SelectedMET'] = met_to_use

    return passing_jets

