#!/usr/bin/env python

import time
tic = time.time()

from coffea import hist, processor
from coffea.nanoevents import NanoAODSchema
from coffea.analysis_tools import PackedSelection
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

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
base_jobid = os.environ['base_jobid']

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('fset', type=str, help='Fileset dictionary (in string form) to be used for the processor')
parser.add_argument('year', choices=['2016', '2017', '2018'] if base_jobid == 'NanoAODv6' else ['2016APV', '2016', '2017', '2018'], help='Specify which year to run over')
parser.add_argument('outfname', type=str, help='Specify output filename, including directory and file extension')
parser.add_argument('opts', type=str, help='Fileset dictionary (in string form) to be used for the processor')
args = parser.parse_args()

# convert input string of fileset dictionary to actual dictionary
fdict = (args.fset).replace("\'", "\"")
fileset = prettyjson.loads(fdict)

if len(fileset.keys()) > 1:
    raise ValueError("Only one topology run at a time in order to determine which corrections and systematics to run")
samplename = list(fileset.keys())[0]

if samplename != "ttJetsSL": raise ValueError("Should only be run on ttbar l+jets datasets")

# convert input string of options dictionary to actual dictionary
odict = (args.opts).replace("\'", "\"")
opts_dict = prettyjson.loads(odict)

    ## set config options passed through argparse
import ast
to_debug = ast.literal_eval(opts_dict.get('debug', 'False'))

#   ## specify ttJets samples
#ttJets_semilep = ['ttJetsSL'] if base_jobid == 'ULnanoAOD' else ['ttJetsSL', 'ttJets_PS']
#isttJets_semilep_ = samplename in ttJets_semilep
#if not isttJets_semilep_:
#    raise ValueError("Should only be run on ttbar l+jets datasets")

nnlo_reweighting = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'NNLO_to_Tune_Orig_Interp_Ratios_%s.coffea' % base_jobid))[args.year]

# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class Analyzer(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.reweighting_axis = hist.Cat("rewt", "Reweighting Type")
        self.pt_axis = hist.Bin("pt", "p_{T} [GeV]",  200, 0, 1000.)
        self.eta_axis = hist.Bin("eta", r"$\eta$", 200, -5., 5.)
        self.phi_axis = hist.Bin("phi", r"$\phi$", 160, -4, 4)
        self.mtop_axis = hist.Bin("mtop", "m(top) [GeV]", 300, 0, 300)
        self.mtt_axis = hist.Bin("mtt", "m($t\overline{t}$) [GeV]", 760, 200, 4000)
        self.ctstar_axis = hist.Bin("ctstar", "cos($\\theta^{*}$)", 200, -1., 1.)
        self.ctstar_abs_axis = hist.Bin("ctstar_abs", "|cos($\\theta^{*}$)|", 200, 0., 1.)
        self.mtt_2D_axis = hist.Bin("mtt_2D", "m($t\overline{t}$) [GeV]", np.array([250., 420., 520., 620., 800., 1000., 3500.])) ## same binning as in MATRIX17
        self.ctstar_2D_axis = hist.Bin("ctstar_2D", "cos($\\theta^{*}$)", np.array([-1., -0.65, -0.3, 0., 0.3, 0.65, 1.])) ## same binning as in MATRIX17

            ## make dictionary of hists
        histo_dict = {}
            # thad
        histo_dict['pt_thad'] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.pt_axis)
        histo_dict['eta_thad'] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.eta_axis)
        histo_dict['phi_thad'] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.phi_axis)
        histo_dict['mass_thad'] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.mtop_axis)
        histo_dict['ctstar_thad'] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.ctstar_axis)
        histo_dict['ctstar_abs_thad'] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.ctstar_abs_axis)

            # tlep
        histo_dict['pt_tlep'] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.pt_axis)
        histo_dict['eta_tlep'] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.eta_axis)
        histo_dict['phi_tlep'] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.phi_axis)
        histo_dict['mass_tlep'] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.mtop_axis)
        histo_dict['ctstar_tlep'] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.ctstar_axis)
        histo_dict['ctstar_abs_tlep'] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.ctstar_abs_axis)

            # ttbar
        histo_dict['pt_tt'] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.pt_axis)
        histo_dict['eta_tt'] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.eta_axis)
        histo_dict['phi_tt'] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.phi_axis)
        histo_dict['mtt'] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.mtt_axis)

        histo_dict['mtt_vs_thad_ctstar'] = hist.Hist("Events", self.dataset_axis, self.reweighting_axis, self.mtt_2D_axis, self.ctstar_2D_axis)

        self.reweighting = nnlo_reweighting

        self._accumulator = processor.dict_accumulator(histo_dict)
        self.sample_name = ''
    
    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):
        output = self.accumulator.identity()

        self.sample_name = events.metadata['dataset']

        genpsel.select(events, mode='NORMAL')
        sl_evts = ak.num(events['SL']) > 0
        el_mu_evts = ak.flatten(np.abs(events['SL']['Lepton'].pdgId != 15), axis=None)

        genWeights = events['genWeight'][sl_evts] # only select semilep evts

        thad_ctstar, tlep_ctstar = make_vars.ctstar(events['SL']['THad'], events['SL']['TLep'])
        thad_ctstar, tlep_ctstar = ak.flatten(thad_ctstar, axis=None), ak.flatten(tlep_ctstar, axis=None)

        thad_pt_weights = self.reweighting['thad_pt'](ak.flatten(events['SL']['THad'].pt, axis=None))
        mtt_vs_thad_ctstar_weights = self.reweighting['mtt_vs_thad_ctstar'](ak.flatten(events['SL']['TTbar'].mass, axis=None), thad_ctstar)
        thad_pt_Interp_weights = self.reweighting['thad_pt_Interp'](ak.flatten(events['SL']['THad'].pt, axis=None))
        mtt_vs_thad_ctstar_Interp_weights = self.reweighting['mtt_vs_thad_ctstar_Interp'](ak.flatten(events['SL']['TTbar'].mass, axis=None), thad_ctstar)

            # fill hists
        for rewt_type in ['Nominal', 'thad_pt', 'mtt_vs_thad_ctstar', 'thad_pt_Interp', 'mtt_vs_thad_ctstar_Interp']:
            #set_trace()
            if rewt_type == 'thad_pt':
                evt_weights = genWeights*thad_pt_weights
            elif rewt_type == 'mtt_vs_thad_ctstar':
                evt_weights = genWeights*mtt_vs_thad_ctstar_weights
            elif rewt_type == 'thad_pt_Interp':
                evt_weights = genWeights*thad_pt_Interp_weights
            elif rewt_type == 'mtt_vs_thad_ctstar_Interp':
                evt_weights = genWeights*mtt_vs_thad_ctstar_Interp_weights
            else:
                evt_weights = genWeights

                    # thad
            output['pt_thad'].fill(dataset=self.sample_name, rewt=rewt_type, pt=ak.flatten(events['SL']['THad'].pt, axis=None)[el_mu_evts], weight=evt_weights[el_mu_evts])
            output['eta_thad'].fill(dataset=self.sample_name, rewt=rewt_type, eta=ak.flatten(events['SL']['THad'].eta, axis=None)[el_mu_evts], weight=evt_weights[el_mu_evts])
            output['phi_thad'].fill(dataset=self.sample_name, rewt=rewt_type, phi=ak.flatten(events['SL']['THad'].phi, axis=None)[el_mu_evts], weight=evt_weights[el_mu_evts])
            output['mass_thad'].fill(dataset=self.sample_name, rewt=rewt_type, mtop=ak.flatten(events['SL']['THad'].mass, axis=None)[el_mu_evts], weight=evt_weights[el_mu_evts])
            output['ctstar_thad'].fill(dataset=self.sample_name, rewt=rewt_type, ctstar=thad_ctstar[el_mu_evts], weight=evt_weights[el_mu_evts])
            output['ctstar_abs_thad'].fill(dataset=self.sample_name, rewt=rewt_type, ctstar_abs=np.abs(thad_ctstar[el_mu_evts]), weight=evt_weights[el_mu_evts])

                    # tlep
            output['pt_tlep'].fill(dataset=self.sample_name, rewt=rewt_type, pt=ak.flatten(events['SL']['TLep'].pt, axis=None)[el_mu_evts], weight=evt_weights[el_mu_evts])
            output['eta_tlep'].fill(dataset=self.sample_name, rewt=rewt_type, eta=ak.flatten(events['SL']['TLep'].eta, axis=None)[el_mu_evts], weight=evt_weights[el_mu_evts])
            output['phi_tlep'].fill(dataset=self.sample_name, rewt=rewt_type, phi=ak.flatten(events['SL']['TLep'].phi, axis=None)[el_mu_evts], weight=evt_weights[el_mu_evts])
            output['mass_tlep'].fill(dataset=self.sample_name, rewt=rewt_type, mtop=ak.flatten(events['SL']['TLep'].mass, axis=None)[el_mu_evts], weight=evt_weights[el_mu_evts])
            output['ctstar_tlep'].fill(dataset=self.sample_name, rewt=rewt_type, ctstar=tlep_ctstar[el_mu_evts], weight=evt_weights[el_mu_evts])
            output['ctstar_abs_tlep'].fill(dataset=self.sample_name, rewt=rewt_type, ctstar_abs=np.abs(tlep_ctstar[el_mu_evts]), weight=evt_weights[el_mu_evts])

                    # ttbar system
            output['pt_tt'].fill(dataset=self.sample_name, rewt=rewt_type, pt=ak.flatten(events['SL']['TTbar'].pt, axis=None)[el_mu_evts], weight=evt_weights[el_mu_evts])
            output['eta_tt'].fill(dataset=self.sample_name, rewt=rewt_type, eta=ak.flatten(events['SL']['TTbar'].eta, axis=None)[el_mu_evts], weight=evt_weights[el_mu_evts])
            output['phi_tt'].fill(dataset=self.sample_name, rewt=rewt_type, phi=ak.flatten(events['SL']['TTbar'].phi, axis=None)[el_mu_evts], weight=evt_weights[el_mu_evts])
            output['mtt'].fill(dataset=self.sample_name, rewt=rewt_type, mtt=ak.flatten(events['SL']['TTbar'].mass, axis=None)[el_mu_evts], weight=evt_weights[el_mu_evts])
            output['mtt_vs_thad_ctstar'].fill(dataset=self.sample_name, rewt=rewt_type, mtt_2D=ak.flatten(events['SL']['TTbar'].mass, axis=None)[el_mu_evts], ctstar_2D=thad_ctstar[el_mu_evts], weight=evt_weights[el_mu_evts])


        return output


    def postprocess(self, accumulator):
        return accumulator


proc_executor = processor.iterative_executor if to_debug else processor.futures_executor
proc_exec_args = {"schema": NanoAODSchema} if to_debug else {"schema": NanoAODSchema, "workers": 8}
output = processor.run_uproot_job(
    fileset,
    treename="Events",
    processor_instance=Analyzer(),
    executor=proc_executor,
    executor_args=proc_exec_args,
    #chunksize=10000 if to_debug else 100000,
    chunksize=100000,
)

save(output, args.outfname)
print('%s has been written' % args.outfname)

toc = time.time()
print("Total time: %.1f" % (toc - tic))
