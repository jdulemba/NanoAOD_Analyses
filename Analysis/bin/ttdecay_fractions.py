#!/usr/bin/env python

from coffea import hist
from coffea.util import save, load
import coffea.processor as processor
from pdb import set_trace
import os, sys
import numpy as np
import Utilities.prettyjson as prettyjson
import Utilities.make_variables as make_vars
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

            ## make binning for hists
                ## binning made to be compatible with hists provided by DESY group for validation
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.pt_axis = hist.Bin("pt", "p_{T} [GeV]", 200, 0, 1000)
        self.eta_axis = hist.Bin("eta", r"$\eta$", 200, -5., 5.)
        self.rapidity_axis = hist.Bin("rap", "rapidity", 200, -5., 5.)
        self.phi_axis = hist.Bin("phi", r"$\phi$", 160, -4, 4)
        self.energy_axis = hist.Bin("energy", "E [GeV]", 200, 0, 1000)
        self.mtop_axis = hist.Bin("mtop", "m(top) [GeV]", 100, 150, 200)
        self.mtt_axis = hist.Bin("mtt", "m($t\overline{t}$) [GeV]", 170, 250, 1500)
        self.ctstar_axis = hist.Bin("ctstar", "cos($\\theta^{*}$)", 24, -1., 1.)
        #self.mtt_axis = hist.Bin("mtt", "m($t\overline{t}$) [GeV]", 360, 200, 2000)
        #self.ctstar_axis = hist.Bin("ctstar", "cos($\\theta^{*}$)", 40, -1., 1.)

            ## make dictionary of hists
        histo_dict = {}
                ## make selection plots
        selection_hists = self.make_selection_hists()
        histo_dict.update(selection_hists)        

        self._accumulator = processor.dict_accumulator(histo_dict)

    
    @property
    def accumulator(self):
        return self._accumulator



    def make_selection_hists(self):
        histo_dict = {}
        histo_dict['mtt']      = hist.Hist("Events", self.dataset_axis, self.mtt_axis)
        histo_dict['top_mass']    = hist.Hist("Events", self.dataset_axis, self.mtop_axis)
        histo_dict['antitop_mass']    = hist.Hist("Events", self.dataset_axis, self.mtop_axis)
        histo_dict['top_pt']  = hist.Hist("Events", self.dataset_axis, self.pt_axis)
        histo_dict['antitop_pt']  = hist.Hist("Events", self.dataset_axis, self.pt_axis)
        histo_dict['ttbar_pt']    = hist.Hist("Events", self.dataset_axis, self.pt_axis)
        histo_dict['top_eta'] = hist.Hist("Events", self.dataset_axis, self.eta_axis)
        histo_dict['antitop_eta'] = hist.Hist("Events", self.dataset_axis, self.eta_axis)
        histo_dict['ttbar_eta']   = hist.Hist("Events", self.dataset_axis, self.eta_axis)
        histo_dict['top_rapidity'] = hist.Hist("Events", self.dataset_axis, self.rapidity_axis)
        histo_dict['antitop_rapidity'] = hist.Hist("Events", self.dataset_axis, self.rapidity_axis)
        histo_dict['ttbar_rapidity']   = hist.Hist("Events", self.dataset_axis, self.rapidity_axis)

        histo_dict['top_ctstar']     = hist.Hist("Events", self.dataset_axis, self.ctstar_axis)
        histo_dict['antitop_ctstar']     = hist.Hist("Events", self.dataset_axis, self.ctstar_axis)
        #histo_dict['tlep_ctstar_abs'] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.ctstar_abs_axis)

        #histo_dict['mtt_vs_tlep_ctstar_abs'] = hist.Hist("Events", self.dataset_axis, self.sys_axis, self.jetmult_axis, self.leptype_axis, self.btag_axis, self.lepcat_axis, self.mtt_axis, self.ctstar_abs_axis)
        histo_dict['top_rap_vs_ctstar'] = hist.Hist("Events", self.dataset_axis, self.ctstar_axis, self.rapidity_axis)

        return histo_dict

    def process(self, df):
        np.random.seed(10) # sets seed so values from random distributions are reproducible (JER corrections)
        output = self.accumulator.identity()

        genparts = genpsel.process_genParts(df)
        #genparts = df['GenParts']
        #GenTTbar = genpsel.select(df, systype='FINAL', mode='NORMAL')

        #set_trace()
        is_first_copy = genparts.statusFlags >> 12 & 1 == 1
        is_hard_process = genparts.statusFlags >> 7 & 1 == 1
        hard_gps = genparts[is_first_copy & is_hard_process]
        abspdg = abs(hard_gps.pdgId)
        sgn = np.sign(hard_gps.pdgId)

        gen_top = hard_gps[(hard_gps.pdgId == 6)]
        gen_tbar = hard_gps[(hard_gps.pdgId == -6)]
        gen_b = hard_gps[(hard_gps.pdgId == 5) & (hard_gps.mompdgId == 6)]
        gen_bbar = hard_gps[(hard_gps.pdgId == -5) & (hard_gps.mompdgId == -6)]
        gen_wplus = hard_gps[(hard_gps.pdgId == 24) & (hard_gps.mompdgId == 6)]
        gen_wminus = hard_gps[(hard_gps.pdgId == -24) & (hard_gps.mompdgId == -6)]

        gen_wpartons_up = hard_gps[(np.mod(hard_gps.pdgId, 2) == 0) & (abspdg < 6) & (hard_gps.mompdgId == sgn * 24)]
        gen_wpartons_dw = hard_gps[(np.mod(hard_gps.pdgId, 2) == 1) & (abspdg < 6) & (hard_gps.mompdgId == sgn * -24)]

        gen_charged_leptons = hard_gps[((abspdg == 11) | (abspdg == 13)) & (hard_gps.mompdgId == sgn * -24)]
        gen_neutral_leptons = hard_gps[((abspdg == 12) | (abspdg == 14)) & (hard_gps.mompdgId == sgn * 24)]
        gen_taus = hard_gps[(abspdg == 15) & (hard_gps.mompdgId == sgn * -24)]

            # set charge
        gen_top.charge = gen_top.charge.ones_like()*(2./3.)
        gen_tbar.charge = gen_tbar.charge.ones_like()*(-2./3.)
        gen_b.charge = gen_b.charge.ones_like()*(-1./3.)
        gen_bbar.charge = gen_bbar.charge.ones_like()*(1./3.)
        gen_wplus.charge = gen_wplus.charge.ones_like()
        gen_wminus.charge = gen_wminus.charge.ones_like()*(-1.)

        gen_wpartons_up.charge = sgn[(np.mod(hard_gps.pdgId, 2) == 0) & (abspdg < 6) & (hard_gps.mompdgId == sgn * 24)]*(2./3.)
        gen_wpartons_dw.charge = sgn[(np.mod(hard_gps.pdgId, 2) == 1) & (abspdg < 6) & (hard_gps.mompdgId == sgn * -24)]*(-1./3.)

        gen_charged_leptons.charge = sgn[((abspdg == 11) | (abspdg == 13)) & (hard_gps.mompdgId == sgn * -24)]* -1
        gen_neutral_leptons.charge = gen_neutral_leptons.charge.zeros_like()
        gen_taus.charge = sgn[(abspdg == 15) & (hard_gps.mompdgId == sgn * -24)]* -1

        GenObjs = genpsel.select_normal(df, 24)
        set_trace()
        selection = processor.PackedSelection()
        selection.add('isLastCopy', ((genparts.statusFlags & 8192) == 8192).flatten() )
        selection.add('top_quarks_pdgId', (np.abs(genparts.pdgId) == 6).flatten())
        #selection.add('Wbosons_pdgId', (np.abs(genparts.pdgId) == 24).flatten())
        #selection.add('top_quarks', selection.require(isLastCopy=True, top_quarks_pdgId=True) )
        #selection.add('Wbosons', ( ((genparts[genparts.momIdx].statusFlags & 8192) == 8192) & ( np.abs(genparts[genparts.momIdx].pdgId) == 6 ) & (np.abs(genparts.pdgId) == 24) ).flatten())
        selection.add('W_comes_from_top', ( ((genparts[genparts.momIdx].statusFlags & 8192) == 8192) & ( np.abs(genparts[genparts.momIdx].pdgId) == 6 ) & (np.abs(genparts.pdgId) == 24) ).flatten())
        selection.add('comes_from_W', ( (np.abs(genparts[genparts.momIdx].pdgId) == 24) & (np.abs(genparts.pdgId) != 24) ).flatten())
        selection.add('electron_pdgId', (np.abs(genparts.pdgId) == 11).flatten())
        selection.add('el_nu_pdgId', (np.abs(genparts.pdgId) == 12).flatten())
        selection.add('muon_pdgId', (np.abs(genparts.pdgId) == 13).flatten())
        selection.add('mu_nu_pdgId', (np.abs(genparts.pdgId) == 14).flatten())
        selection.add('tau_pdgId', (np.abs(genparts.pdgId) == 15).flatten())
        selection.add('tau_nu_pdgId', (np.abs(genparts.pdgId) == 16).flatten())
        selection.add('comes_from_tau', ( (np.abs(genparts[genparts.momIdx].pdgId) == 15) & (np.abs(genparts.pdgId) != 15) ).flatten())
        tau_decays = processor.PackedSelection()
        tau_decays.add('chargedPion_neutralPion_tauNu', ( (np.abs(genparts[genparts.momIdx].pdgId) == 15) & (np.abs(genparts.pdgId) != 15) ).flatten())

        top_sel = awkward.JaggedArray.fromcounts(genparts.counts, selection.require(isLastCopy=True, top_quarks_pdgId=True) )
        top_quarks = genparts[top_sel]
        Wbos_sel = awkward.JaggedArray.fromcounts(genparts.counts, selection.require(W_comes_from_top=True))
        w_bosons = genparts[Wbos_sel]

        electron_sel = awkward.JaggedArray.fromcounts(genparts.counts, selection.require(comes_from_W=True, electron_pdgId=True))
        electrons = genparts[electron_sel]
        el_nu_sel = awkward.JaggedArray.fromcounts(genparts.counts, selection.require(comes_from_W=True, el_nu_pdgId=True))
        el_nus = genparts[el_nu_sel]

        muon_sel = awkward.JaggedArray.fromcounts(genparts.counts, selection.require(comes_from_W=True, muon_pdgId=True))
        muons = genparts[muon_sel]
        mu_nu_sel = awkward.JaggedArray.fromcounts(genparts.counts, selection.require(comes_from_W=True, mu_nu_pdgId=True))
        mu_nus = genparts[mu_nu_sel]

        tau_sel = awkward.JaggedArray.fromcounts(genparts.counts, selection.require(comes_from_W=True, tau_pdgId=True))
        taus = genparts[tau_sel]
        tau_nu_sel = awkward.JaggedArray.fromcounts(genparts.counts, selection.require(comes_from_W=True, tau_nu_pdgId=True))
        tau_nus = genparts[tau_nu_sel]
        tau_decay_prods_sel = awkward.JaggedArray.fromcounts(genparts.counts, selection.require(comes_from_tau=True))
        tau_decay_prods = genparts[tau_decay_prods_sel]

        all_electron_nu_pairs = electrons.cross(el_nus)
        valid_e_nu_pairs = all_electron_nu_pairs[(all_electron_nu_pairs.i0.momIdx == all_electron_nu_pairs.i1.momIdx)]
        all_muon_nu_pairs = muons.cross(mu_nus)
        valid_mu_nu_pairs = all_muon_nu_pairs[(all_muon_nu_pairs.i0.momIdx == all_muon_nu_pairs.i1.momIdx)]
        all_tau_nu_pairs = taus.cross(tau_nus)
        valid_tau_nu_pairs = all_tau_nu_pairs[(all_tau_nu_pairs.i0.momIdx == all_tau_nu_pairs.i1.momIdx)]
        set_trace()
        
            # dict of requirements for each parton
        leptons = {
            'electrons' : {'comes_from_W', 'electron_pdgId'},
            'muons' : {'comes_from_W', 'muon_pdgId'},
            'taus' : {'comes_from_W', 'tau_pdgId'},
            #'top' : {'isLastCopy', 'top_pdgId'},
            #'antitop' : {'isLastCopy', 'antitop_pdgId'},
            #'Wboson' : {'isLastCopy'
            ##'top' : {'quarks', 'top_pdgId'},
            ##'antitop' : {'quarks', 'antitop_pdgId'},
        }

        comes_from_W_sel = awkward.JaggedArray.fromcounts(genparts.counts, selection.require(comes_from_W=True))
        electron_sel = awkward.JaggedArray.fromcounts(genparts.counts, selection.all(*leptons['electrons']))
        muon_sel = awkward.JaggedArray.fromcounts(genparts.counts, selection.all(*leptons['muons']))
        tau_sel = awkward.JaggedArray.fromcounts(genparts.counts, selection.all(*leptons['taus']))

        electrons = genparts[electron_sel]
        muons = genparts[muon_sel]
        taus = genparts[tau_sel]

        top_partons = genparts[top_sel]
        antitop_partons = genparts[antitop_sel]

        gen_mtt = (top_partons.p4+antitop_partons.p4).mass.flatten()
        gen_top_ctstar, gen_antitop_ctstar = make_vars.ctstar_flat(top_partons.p4.flatten(), antitop_partons.p4.flatten())

        #set_trace()
            # fill plots for SM ttbar
        output = self.fill_selection_hists(acc=output, dname=df.dataset, top=top_partons.flatten(), antitop=antitop_partons.flatten(), evt_weights=evt_wts)
        output['top_ctstar'].fill(dataset=df.dataset, ctstar=gen_top_ctstar, weight=evt_wts)
        output['antitop_ctstar'].fill(dataset=df.dataset, ctstar=gen_antitop_ctstar, weight=evt_wts)
        output['top_rap_vs_ctstar'].fill(dataset=df.dataset, ctstar=gen_top_ctstar, rap=top_partons.flatten().p4.rapidity, weight=evt_wts)
        
        ##set_trace()
        #    ### get signal weights and fill hists
        #for signal in self.signals_to_run:
        #    #set_trace()
        #    pos_wts = self.corrections['Signal'][signal]['pos']['Central'](gen_mtt, gen_top_ctstar)
        #    output = self.fill_selection_hists(acc=output, dname='%s_%s_pos' % (df.dataset, signal), top=top_partons.flatten(), antitop=antitop_partons.flatten(), evt_weights=pos_wts*evt_wts)
        #    output['top_ctstar'].fill(dataset='%s_%s_pos' % (df.dataset, signal), ctstar=gen_top_ctstar, weight=pos_wts*evt_wts)
        #    output['antitop_ctstar'].fill(dataset='%s_%s_pos' % (df.dataset, signal), ctstar=gen_antitop_ctstar, weight=pos_wts*evt_wts)
        #    output['top_rap_vs_ctstar'].fill(dataset='%s_%s_pos' % (df.dataset, signal), ctstar=gen_top_ctstar, rap=top_partons.flatten().p4.rapidity, weight=pos_wts*evt_wts)
        #    if 'Int' in signal:
        #        neg_wts = self.corrections['Signal'][signal]['neg']['Central'](gen_mtt, gen_top_ctstar)
        #        output = self.fill_selection_hists(acc=output, dname='%s_%s_neg' % (df.dataset, signal), top=top_partons.flatten(), antitop=antitop_partons.flatten(), evt_weights=neg_wts*evt_wts)
        #        output['top_ctstar'].fill(dataset='%s_%s_neg' % (df.dataset, signal), ctstar=gen_top_ctstar, weight=neg_wts*evt_wts)
        #        output['antitop_ctstar'].fill(dataset='%s_%s_neg' % (df.dataset, signal), ctstar=gen_antitop_ctstar, weight=neg_wts*evt_wts)
        #        output['top_rap_vs_ctstar'].fill(dataset='%s_%s_neg' % (df.dataset, signal), ctstar=gen_top_ctstar, rap=top_partons.flatten().p4.rapidity, weight=neg_wts*evt_wts)

        return output



    def fill_selection_hists(self, acc, dname, top, antitop, evt_weights):
        acc['mtt'].fill(         dataset=dname, mtt=(top.p4 + antitop.p4).mass, weight=evt_weights)
        acc['top_mass'].fill(    dataset=dname, mtop=top.p4.mass, weight=evt_weights)
        acc['top_pt'].fill(      dataset=dname, pt=top.p4.pt, weight=evt_weights)
        acc['top_eta'].fill(     dataset=dname, eta=top.p4.eta, weight=evt_weights)
        acc['antitop_mass'].fill(dataset=dname, mtop=antitop.p4.mass, weight=evt_weights)
        acc['antitop_pt'].fill(  dataset=dname, pt=antitop.p4.pt, weight=evt_weights)
        acc['antitop_eta'].fill( dataset=dname, eta=antitop.p4.eta, weight=evt_weights)
        acc['ttbar_pt'].fill(    dataset=dname, pt=(top.p4 + antitop.p4).pt, weight=evt_weights)
        acc['ttbar_eta'].fill(   dataset=dname, eta=(top.p4 + antitop.p4).eta, weight=evt_weights)
        acc['top_rapidity'].fill(dataset=dname, rap=top.p4.rapidity, weight=evt_weights)
        acc['antitop_rapidity'].fill(dataset=dname, rap=antitop.p4.rapidity, weight=evt_weights)
        acc['ttbar_rapidity'].fill(dataset=dname, rap=(top.p4 + antitop.p4).rapidity, weight=evt_weights)

        #acc['tlep_ctstar_abs'].fill(dataset=dname, ctstar_abs=np.abs(tlep_ctstar), weight=evt_weights)
        #acc['mtt_vs_tlep_ctstar_abs'].fill(dataset=dname, mtt=(top.p4 + antitop.p4).mass, ctstar_abs=np.abs(lep_ctstar), weight=evt_weights)

        return acc        



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
        #'nano': True,
    },
    chunksize=10000 if args.debug else 100000,
    #chunksize=10000 if args.debug else 50000,
    #chunksize=50000,
)

save(output, args.outfname)
print('%s has been written' % args.outfname)
