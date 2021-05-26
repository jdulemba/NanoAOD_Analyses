#!/usr/bin/env python

from coffea import hist, processor
from coffea.nanoevents import NanoAODSchema
import awkward as ak
from coffea.nanoevents.methods import vector
ak.behavior.update(vector.behavior)
from pdb import set_trace
import os
import python.GenParticleSelector as genpsel
from coffea.util import load, save
import numpy as np
import Utilities.make_variables as make_vars
import Utilities.prettyjson as prettyjson

base_jobid = os.environ['base_jobid']

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('fset', type=str, help='Fileset dictionary (in string form) to be used for the processor')
parser.add_argument('year', choices=['2016APV', '2016', '2017', '2018'] if base_jobid == 'ULnanoAOD' else ['2016', '2017', '2018'], help='Specify which year to run over')
parser.add_argument('outfname', type=str, help='Specify output filename, including directory and file extension')
parser.add_argument('--debug', action='store_true', help='Uses iterative_executor for debugging purposes, otherwise futures_excutor will be used (faster)')
args = parser.parse_args()

# convert input string of fileset dictionary to actual dictionary
fdict = (args.fset).replace("\'", "\"")
fileset = prettyjson.loads(fdict)

if len(fileset.keys()) > 1:
    raise ValueError("Only one topology run at a time in order to determine which corrections and systematics to run")
samplename = list(fileset.keys())[0]

   ## specify ttJets samples
ttJets_semilep = ['ttJetsSL'] if base_jobid == 'ULnanoAOD' else ['ttJetsSL', 'ttJets_PS']
isttJets_semilep_ = samplename in ttJets_semilep
if not isttJets_semilep_:
    raise ValueError("Should only be run on ttbar l+jets datasets")

# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class Analyzer(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        #self.pt_axis = hist.Bin("pt", "p_{T} [GeV]", np.array([0., 40., 80., 120., 160., 200., 250., 300., 350., 400., 450., 500., 600., 700., 800., 950., 1100., 1600.])) ## same binning as in MATRIX17
        #self.mtt_axis = hist.Bin("mtt", "mtt [GeV]", np.array([250., 420., 520., 620., 800., 1000., 3500.])) ## same binning as in MATRIX17
        #self.ctstar_axis = hist.Bin("ctstar", "cos theta*", np.array([-1., -0.65, -0.3, 0., 0.3, 0.65, 1.])) ## same binning as in MATRIX17
        self.pt_axis = hist.Bin("pt", "p_{T} [GeV]", 200, 0., 2000.)
        self.mtt_axis = hist.Bin("mtt", "mtt [GeV]", 380, 200., 4000.)
        self.ctstar_axis = hist.Bin("ctstar", "cos theta*", 40, -1., 1.)

            ## make dictionary of hists
        histo_dict = {}
        histo_dict['thad_pt'] = hist.Hist("thad_pt", self.dataset_axis, self.pt_axis)
        histo_dict['mtt_vs_thad_ctstar'] = hist.Hist("mtt_vs_thad_ctstar", self.dataset_axis, self.mtt_axis, self.ctstar_axis)

        self._accumulator = processor.dict_accumulator(histo_dict)
    
    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        output = self.accumulator.identity()

        genWeights = events['genWeight']

            # get e/mu SL gen events
        genpsel.select(events, mode='NORMAL')
        el_mu_evts = (np.abs(events['SL']['Lepton'].pdgId) != 15)
            # get kin vars
        thad_pt = events['SL']['THad'].pt
        thad_ctstar, tlep_ctstar = make_vars.ctstar(events['SL']['THad'], events['SL']['TLep'])
        mtt = events['SL']['TTbar'].mass
            # fill hists
        output['thad_pt'].fill(dataset=samplename, pt=ak.flatten(thad_pt[el_mu_evts], axis=None), weight=genWeights[ak.flatten(el_mu_evts, axis=None)])
        output['mtt_vs_thad_ctstar'].fill(dataset=samplename, mtt=ak.flatten(mtt[el_mu_evts], axis=None), ctstar=ak.flatten(thad_ctstar[el_mu_evts], axis=None), weight=genWeights[ak.flatten(el_mu_evts, axis=None)])

        return output


    def postprocess(self, accumulator):
        return accumulator


proc_executor = processor.iterative_executor if args.debug else processor.futures_executor
proc_exec_args = {"schema": NanoAODSchema} if args.debug else {"schema": NanoAODSchema, "workers": 8}
output = processor.run_uproot_job(
    fileset,
    treename="Events",
    processor_instance=Analyzer(),
    executor=proc_executor,
    executor_args=proc_exec_args,
    chunksize=100000,
    #chunksize=10000 if args.debug else 100000,
)

save(output, args.outfname)
print('%s has been written' % args.outfname)
