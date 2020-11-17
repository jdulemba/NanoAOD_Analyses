#!/usr/bin/env python

from coffea import hist
from coffea.util import save, load
import coffea.processor as processor
from pdb import set_trace
import os, sys
import numpy as np
import Utilities.prettyjson as prettyjson
#import Utilities.make_variables as make_vars
import python.GenParticleSelector as genpsel
import re
import awkward

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'ttdecay_fractions'

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('fset', type=str, help='Fileset dictionary (in string form) to be used for the processor')
parser.add_argument('year', choices=['2016', '2017', '2018'], help='Specify which year to run over')
parser.add_argument('outfname', type=str, help='Specify output filename, including directory and file extension')
parser.add_argument('--signal', type=str, help='Specify output filename, including directory and file extension')
parser.add_argument('--debug', action='store_true', help='Uses iterative_executor for debugging purposes, otherwise futures_excutor will be used (faster)')
parser.add_argument('--evt_sys', type=str, help='Specify event systematics to run. Default is (NONE,NONE) and all opts are capitalized through run_analyzer')
parser.add_argument('--rewt_sys', type=str, help='Specify reweighting systematics to run. Default is (NONE,NONE) and all opts are capitalized through run_analyzer')
parser.add_argument('--only_sys', type=int, help='Only run specified systematics and not nominal weights (nosys)')

args = parser.parse_args()

# convert input string of fileset dictionary to actual dictionary
fdict = (args.fset).replace("\'", "\"")
fileset = prettyjson.loads(fdict)

if len(fileset.keys()) > 1:
    raise ValueError("Only one topology run at a time in order to determine which corrections and systematics to run")
samplename = list(fileset.keys())[0]

# get dataset classification, used for corrections/systematics
    ## specify ttJets samples

is_ttJets_ = samplename in ['ttJetsSL', 'ttJetsHad', 'ttJetsDiLep']
if not is_ttJets_:
    raise ValueError("This should only be run on nominal ttbar events!")

    ## get list of signal samples to run over
# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class ttdecay_fractions(processor.ProcessorABC):
    def __init__(self):

            ## make dictionary of hists
        histo_dict = {}
        for sample in fileset.keys():
            histo_dict['%s_SL' % sample] = processor.defaultdict_accumulator(int)
            histo_dict['%s_Had' % sample] = processor.defaultdict_accumulator(int)
            histo_dict['%s_DiLep' % sample] = processor.defaultdict_accumulator(int)

        self._accumulator = processor.dict_accumulator(histo_dict)

    
    @property
    def accumulator(self):
        return self._accumulator


    def process(self, df):
        np.random.seed(10) # sets seed so values from random distributions are reproducible (JER corrections)
        output = self.accumulator.identity()

        genparts = genpsel.process_genParts(df)

        #set_trace()
            # get tops, defined as last copy
        is_last_copy = genparts.statusFlags >> 13 & 1 == 1
        gen_tops = genparts[(is_last_copy) & (genparts.pdgId == 6)]
        gen_tbars = genparts[(is_last_copy) & (genparts.pdgId == -6)]


            # get direct top decay products (will be first copy)
        is_first_copy = genparts.statusFlags >> 12 & 1 == 1
        is_hard_process = genparts.statusFlags >> 7 & 1 == 1
        hard_gps = genparts[is_first_copy & is_hard_process]
        abspdg = abs(hard_gps.pdgId)
        sgn = np.sign(hard_gps.pdgId)

        gen_bs = hard_gps[(hard_gps.pdgId == 5) & (hard_gps.mompdgId == 6)]
        gen_bbars = hard_gps[(hard_gps.pdgId == -5) & (hard_gps.mompdgId == -6)]
        gen_wplus = hard_gps[(hard_gps.pdgId == 24) & (hard_gps.mompdgId == 6)]
        gen_wminus = hard_gps[(hard_gps.pdgId == -24) & (hard_gps.mompdgId == -6)]

        gen_wpartons_up = hard_gps[(np.mod(hard_gps.pdgId, 2) == 0) & (abspdg < 6) & (hard_gps.mompdgId == sgn * 24)]
        gen_wpartons_dw = hard_gps[(np.mod(hard_gps.pdgId, 2) == 1) & (abspdg < 6) & (hard_gps.mompdgId == sgn * -24)]

        gen_els = hard_gps[(abspdg == 11) & (hard_gps.mompdgId == sgn * -24)]
        gen_nu_els = hard_gps[(abspdg == 12) & (hard_gps.mompdgId == sgn * 24)]
        gen_mus = hard_gps[(abspdg == 13) & (hard_gps.mompdgId == sgn * -24)]
        gen_nu_mus = hard_gps[(abspdg == 14) & (hard_gps.mompdgId == sgn * 24)]
        gen_taus = hard_gps[(abspdg == 15) & (hard_gps.mompdgId == sgn * -24)]
        gen_nu_taus = hard_gps[(abspdg == 16) & (hard_gps.mompdgId == sgn * 24)]

                # get direct tau decay products from hard processes (subset of gen_taus events above)
        isDirectHardProcessTauDecayProduct = genparts.statusFlags >> 10 & 1 == 1
        tau_decay_prods = genparts[isDirectHardProcessTauDecayProduct]

        tau_charged_kaons = tau_decay_prods[np.abs(tau_decay_prods.pdgId) == 321]
        tau_charged_pions = tau_decay_prods[np.abs(tau_decay_prods.pdgId) == 211]
        tau_neutral_pions = tau_decay_prods[np.abs(tau_decay_prods.pdgId) == 111]
        tau_tau_nu = tau_decay_prods[np.abs(tau_decay_prods.pdgId) == 16]
        tau_electron = tau_decay_prods[np.abs(tau_decay_prods.pdgId) == 11]
        tau_electron_nu = tau_decay_prods[np.abs(tau_decay_prods.pdgId) == 12]
        tau_muon = tau_decay_prods[(np.abs(tau_decay_prods.pdgId) == 13)]
        tau_muon_nu = tau_decay_prods[np.abs(tau_decay_prods.pdgId) == 14]


            # set charge
        gen_tops.charge = gen_tops.charge.ones_like()*(2./3.)
        gen_tbars.charge = gen_tbars.charge.ones_like()*(-2./3.)
        gen_bs.charge = gen_bs.charge.ones_like()*(-1./3.)
        gen_bbars.charge = gen_bbars.charge.ones_like()*(1./3.)
        gen_wplus.charge = gen_wplus.charge.ones_like()
        gen_wminus.charge = gen_wminus.charge.ones_like()*(-1.)

        gen_wpartons_up.charge = sgn[(np.mod(hard_gps.pdgId, 2) == 0) & (abspdg < 6) & (hard_gps.mompdgId == sgn * 24)]*(2./3.)
        gen_wpartons_dw.charge = sgn[(np.mod(hard_gps.pdgId, 2) == 1) & (abspdg < 6) & (hard_gps.mompdgId == sgn * -24)]*(-1./3.)

        gen_els.charge = sgn[(abspdg == 11) & (hard_gps.mompdgId == sgn * -24)]* -1
        gen_nu_els.charge = gen_nu_els.charge.zeros_like()
        gen_mus.charge = sgn[(abspdg == 13) & (hard_gps.mompdgId == sgn * -24)]* -1
        gen_nu_els.charge = gen_nu_mus.charge.zeros_like()
        gen_taus.charge = sgn[(abspdg == 15) & (hard_gps.mompdgId == sgn * -24)]* -1
        gen_nu_taus.charge = gen_nu_mus.charge.zeros_like()

            # same number of electrons and electron nus
        if not (gen_els.counts == gen_nu_els.counts).all(): raise ValueError("Different number of electrons and electron neutrinos in events")
            # same number of muons and muon nus
        if not (gen_mus.counts == gen_nu_mus.counts).all(): raise ValueError("Different number of muons and muon neutrinos in events")
            # same number of taus and tau nus
        if not (gen_taus.counts == gen_nu_taus.counts).all(): raise ValueError("Different number of taus and tau neutrinos in events")

        # event classification
            # semilep evts
        el_jets_evts = (gen_els.counts == 1) & (gen_wpartons_up.counts == 1) & (gen_wpartons_dw.counts == 1)
        mu_jets_evts = (gen_mus.counts == 1) & (gen_wpartons_up.counts == 1) & (gen_wpartons_dw.counts == 1)
        tau_jets_evts = (gen_taus.counts == 1) & (gen_wpartons_up.counts == 1) & (gen_wpartons_dw.counts == 1)
        semilep_evts = (el_jets_evts | mu_jets_evts | tau_jets_evts)
                # tau decays
        semilep_tau_to_tauNu_el_elNu = (tau_electron[tau_jets_evts].counts == 1) & (tau_electron_nu[tau_jets_evts].counts == 1) & (tau_tau_nu[tau_jets_evts].counts == 1)
        semilep_tau_to_tauNu_mu_muNu = (tau_muon[tau_jets_evts].counts == 1) & (tau_muon_nu[tau_jets_evts].counts == 1) & (tau_tau_nu[tau_jets_evts].counts == 1)
        semilep_tau_leptonic_decay = semilep_tau_to_tauNu_el_elNu | semilep_tau_to_tauNu_mu_muNu
        semilep_tau_hadronic_decay = ~semilep_tau_leptonic_decay

        output['%s_SL' % df.dataset]['e'] += el_jets_evts.sum()
        output['%s_SL' % df.dataset]['mu'] += mu_jets_evts.sum()
        output['%s_SL' % df.dataset]['tau->l'] += semilep_tau_leptonic_decay.sum()
        output['%s_SL' % df.dataset]['tau->h'] += semilep_tau_hadronic_decay.sum()
        output['%s_SL' % df.dataset]['Total'] += semilep_evts.sum()


            # all hadronic evts
        #set_trace()
        hadronic_evts = ((gen_els.counts + gen_mus.counts + gen_taus.counts) == 0) & (gen_wpartons_up.counts == 2) & (gen_wpartons_dw.counts == 2)
        output['%s_Had' % df.dataset]['Total'] += hadronic_evts.sum()

            # dilep evts
        el_el_evts = (gen_els.counts == 2) & (gen_wpartons_up.counts == 0) & (gen_wpartons_dw.counts == 0)
        el_mu_evts = (gen_els.counts == 1) & (gen_mus.counts == 1) & (gen_wpartons_up.counts == 0) & (gen_wpartons_dw.counts == 0)
        el_tau_evts = (gen_els.counts == 1) & (gen_taus.counts == 1) & (gen_wpartons_up.counts == 0) & (gen_wpartons_dw.counts == 0)
        mu_mu_evts = (gen_mus.counts == 2) & (gen_wpartons_up.counts == 0) & (gen_wpartons_dw.counts == 0)
        mu_tau_evts = (gen_mus.counts == 1) & (gen_taus.counts == 1) & (gen_wpartons_up.counts == 0) & (gen_wpartons_dw.counts == 0)
        tau_tau_evts = (gen_taus.counts == 2) & (gen_wpartons_up.counts == 0) & (gen_wpartons_dw.counts == 0)
        dilep_evts = (el_el_evts | el_mu_evts | el_tau_evts | mu_mu_evts | mu_tau_evts | tau_tau_evts)
                # tau decays
                    # tau tau
                        # tau tau -> ll
        dilep_TauTau_to_ElEl = (tau_electron[tau_tau_evts].counts == 2) & (tau_electron_nu[tau_tau_evts].counts == 2) & (tau_tau_nu[tau_tau_evts].counts == 2)
        dilep_TauTau_to_MuMu = (tau_muon[tau_tau_evts].counts == 2) & (tau_muon_nu[tau_tau_evts].counts == 2) & (tau_tau_nu[tau_tau_evts].counts == 2)
        dilep_TauTau_to_ElMu = (tau_electron[tau_tau_evts].counts == 1) & (tau_electron_nu[tau_tau_evts].counts == 1) & (tau_muon[tau_tau_evts].counts == 1) & (tau_muon_nu[tau_tau_evts].counts == 1) & (tau_tau_nu[tau_tau_evts].counts == 2)
        dilep_TauTau_ll_decay = (dilep_TauTau_to_ElEl | dilep_TauTau_to_MuMu | dilep_TauTau_to_ElMu)
                        # tau tau -> hh
        dilep_TauTau_hh_decay = ((tau_electron[tau_tau_evts].counts + tau_electron_nu[tau_tau_evts].counts + tau_muon[tau_tau_evts].counts + tau_muon_nu[tau_tau_evts].counts) == 0) & (tau_tau_nu[tau_tau_evts].counts == 2)
                        # tau tau -> lh
        dilep_TauTau_lh_decay = ~(dilep_TauTau_ll_decay | dilep_TauTau_hh_decay)

                # e tau
        dilep_el_tau_to_El = (tau_electron[el_tau_evts].counts == 1) & (tau_electron_nu[el_tau_evts].counts == 1) & (tau_tau_nu[el_tau_evts].counts == 1)
        dilep_el_tau_to_Mu = (tau_muon[el_tau_evts].counts == 1) & (tau_muon_nu[el_tau_evts].counts == 1) & (tau_tau_nu[el_tau_evts].counts == 1)
        dilep_el_tau_leptonic_decay = (dilep_el_tau_to_El | dilep_el_tau_to_Mu)
        dilep_el_tau_hadronic_decay =  ~dilep_el_tau_leptonic_decay

                # mu tau
        dilep_mu_tau_to_El = (tau_electron[mu_tau_evts].counts == 1) & (tau_electron_nu[mu_tau_evts].counts == 1) & (tau_tau_nu[mu_tau_evts].counts == 1)
        dilep_mu_tau_to_Mu = (tau_muon[mu_tau_evts].counts == 1) & (tau_muon_nu[mu_tau_evts].counts == 1) & (tau_tau_nu[mu_tau_evts].counts == 1)
        dilep_mu_tau_leptonic_decay = (dilep_mu_tau_to_El | dilep_mu_tau_to_Mu)
        dilep_mu_tau_hadronic_decay =  ~dilep_mu_tau_leptonic_decay

        output['%s_DiLep' % df.dataset]['e e'] += el_el_evts.sum()
        output['%s_DiLep' % df.dataset]['e mu'] += el_mu_evts.sum()
        output['%s_DiLep' % df.dataset]['e tau->l'] += dilep_el_tau_leptonic_decay.sum()
        output['%s_DiLep' % df.dataset]['e tau->h'] += dilep_el_tau_hadronic_decay.sum()
        output['%s_DiLep' % df.dataset]['mu mu'] += mu_mu_evts.sum()
        output['%s_DiLep' % df.dataset]['mu tau->l'] += dilep_mu_tau_leptonic_decay.sum()
        output['%s_DiLep' % df.dataset]['mu tau->h'] += dilep_mu_tau_hadronic_decay.sum()
        output['%s_DiLep' % df.dataset]['tau tau->ll'] += dilep_TauTau_ll_decay.sum()
        output['%s_DiLep' % df.dataset]['tau tau->lh'] += dilep_TauTau_lh_decay.sum()
        output['%s_DiLep' % df.dataset]['tau tau->hh'] += dilep_TauTau_hh_decay.sum()
        output['%s_DiLep' % df.dataset]['Total'] += dilep_evts.sum()


        #set_trace()

        return output



    def postprocess(self, accumulator):
        return accumulator

proc_executor = processor.iterative_executor if args.debug else processor.futures_executor

#set_trace()
output = processor.run_uproot_job(fileset,
    treename='Events',
    processor_instance=ttdecay_fractions(),
    executor=proc_executor,
    executor_args={
        'workers': 8,
        'flatten' : True,
        'compression': 5,
    },
    chunksize=10000 if args.debug else 100000,
    #chunksize=10000 if args.debug else 50000,
    #chunksize=50000,
)

save(output, args.outfname)
print('%s has been written' % args.outfname)
