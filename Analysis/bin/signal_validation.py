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
analyzer = 'signal_validation'

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
if args.year == '2016':
    Nominal_ttJets = ['ttJets_PS']#, 'ttJets']
else:
    Nominal_ttJets = ['ttJetsSL', 'ttJetsHad', 'ttJetsDiLep']

isNominal_ttJets_ = samplename in Nominal_ttJets
if not isNominal_ttJets_:
    raise ValueError("This should only be run on nominal ttbar events!")

signal_corrections = load(os.path.join(proj_dir, 'Corrections', 'signal_weights_SMppTOtt.coffea'))
#signal_corrections = load(os.path.join(proj_dir, 'Corrections', 'signal_weights_SMggTOtt.coffea'))

    ## get list of signal samples to run over
possible_signals = [sig.strip('\n') for sig in open('%s/inputs/signal_opts.txt' % proj_dir, 'r') if not sig.startswith('#')]
signals_to_use = sorted([sig for sig in possible_signals if re.match(args.signal, sig)])
if not signals_to_use:
    raise ValueError("No signal samples to run!")
print('Running signal samples:', *signals_to_use)

corrections = {
    'Signal' : signal_corrections,
}


# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class signal_validation(processor.ProcessorABC):
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

        self.corrections = corrections
        self.signals_to_run = signals_to_use

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
        selection = processor.PackedSelection()
        selection.add('isLastCopyBeforeFSR', (genparts.statusFlags & 8192).flatten() == 8192)
        selection.add('top_pdgId', (genparts.pdgId == 6).flatten())
        selection.add('antitop_pdgId', (genparts.pdgId == -6).flatten())
        #selection.add('quarks', ( (genparts.status > 21) & (genparts.status < 30) & (genparts.momIdx != -1) ).flatten())

        # use ttbar events with indices % 10 == 1, 2 for applying signal weights
        events = df.event
        ttbar_mask = np.stack([((events % 10) == idx) for idx in [1, 2]], axis=1).any(axis=1)

            # dict of requirements for each parton
        partons = {
            'top' : {'isLastCopyBeforeFSR', 'top_pdgId'},
            'antitop' : {'isLastCopyBeforeFSR', 'antitop_pdgId'},
            #'top' : {'quarks', 'top_pdgId'},
            #'antitop' : {'quarks', 'antitop_pdgId'},
        }

        evt_wts = df['LHEWeight_originalXWGTUP'][ttbar_mask]
        top_sel = awkward.JaggedArray.fromcounts(genparts.counts, selection.all(*partons['top']))
        antitop_sel = awkward.JaggedArray.fromcounts(genparts.counts, selection.all(*partons['antitop']))

        top_partons = genparts[top_sel][ttbar_mask]
        antitop_partons = genparts[antitop_sel][ttbar_mask]

        gen_mtt = (top_partons.p4+antitop_partons.p4).mass.flatten()
        gen_top_ctstar, gen_antitop_ctstar = make_vars.ctstar_flat(top_partons.p4.flatten(), antitop_partons.p4.flatten())

        #set_trace()
            # fill plots for SM ttbar
        output = self.fill_selection_hists(acc=output, dname=df.dataset, top=top_partons.flatten(), antitop=antitop_partons.flatten(), evt_weights=evt_wts)
        output['top_ctstar'].fill(dataset=df.dataset, ctstar=gen_top_ctstar, weight=evt_wts)
        output['antitop_ctstar'].fill(dataset=df.dataset, ctstar=gen_antitop_ctstar, weight=evt_wts)
        output['top_rap_vs_ctstar'].fill(dataset=df.dataset, ctstar=gen_top_ctstar, rap=top_partons.flatten().p4.rapidity, weight=evt_wts)
        
        #set_trace()
            ### get signal weights and fill hists
        for signal in self.signals_to_run:
            #set_trace()
            pos_wts = self.corrections['Signal'][signal]['pos']['Central'](gen_mtt, gen_top_ctstar)
            output = self.fill_selection_hists(acc=output, dname='%s_%s_pos' % (df.dataset, signal), top=top_partons.flatten(), antitop=antitop_partons.flatten(), evt_weights=pos_wts*evt_wts)
            output['top_ctstar'].fill(dataset='%s_%s_pos' % (df.dataset, signal), ctstar=gen_top_ctstar, weight=pos_wts*evt_wts)
            output['antitop_ctstar'].fill(dataset='%s_%s_pos' % (df.dataset, signal), ctstar=gen_antitop_ctstar, weight=pos_wts*evt_wts)
            output['top_rap_vs_ctstar'].fill(dataset='%s_%s_pos' % (df.dataset, signal), ctstar=gen_top_ctstar, rap=top_partons.flatten().p4.rapidity, weight=pos_wts*evt_wts)
            if 'Int' in signal:
                neg_wts = self.corrections['Signal'][signal]['neg']['Central'](gen_mtt, gen_top_ctstar)
                output = self.fill_selection_hists(acc=output, dname='%s_%s_neg' % (df.dataset, signal), top=top_partons.flatten(), antitop=antitop_partons.flatten(), evt_weights=neg_wts*evt_wts)
                output['top_ctstar'].fill(dataset='%s_%s_neg' % (df.dataset, signal), ctstar=gen_top_ctstar, weight=neg_wts*evt_wts)
                output['antitop_ctstar'].fill(dataset='%s_%s_neg' % (df.dataset, signal), ctstar=gen_antitop_ctstar, weight=neg_wts*evt_wts)
                output['top_rap_vs_ctstar'].fill(dataset='%s_%s_neg' % (df.dataset, signal), ctstar=gen_top_ctstar, rap=top_partons.flatten().p4.rapidity, weight=neg_wts*evt_wts)

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

output = processor.run_uproot_job(fileset,
    treename='Events',
    processor_instance=signal_validation(),
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
