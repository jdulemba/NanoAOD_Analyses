#!/usr/bin/env python

from coffea import hist
import coffea.processor as processor
from pdb import set_trace
import os
from argparse import ArgumentParser
import coffea.processor.dataframe
import python.GenParticleSelector as genpsel
from coffea.util import load, save
import numpy as np
import Utilities.make_variables as make_vars
import Utilities.prettyjson as prettyjson
import python.MCWeights as MCWeights

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
base_jobid = os.environ['base_jobid']

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

        #self.reweighting = kin_reweighting
        self.reweighting = nnlo_reweighting

        self._accumulator = processor.dict_accumulator(histo_dict)
        self.sample_name = ''
    
    @property
    def accumulator(self):
        return self._accumulator

    def process(self, df):
        output = self.accumulator.identity()

        if not isinstance(df, coffea.processor.dataframe.LazyDataFrame):
            raise IOError("This function only works for LazyDataFrame objects")

        #if args.debug: set_trace()
        events = df.event
        self.sample_name = df.dataset

        GenTTbar = genpsel.select(df, mode='NORMAL')
        sl_evts = GenTTbar['SL']['TTbar'].counts > 0
        el_mu_evts = (np.abs(GenTTbar['SL']['Lepton'].pdgId) != 15).flatten()

        genWeights = (df.genWeight)[sl_evts] # only select semilep evts

        thad_p4, tlep_p4, ttbar_p4 = GenTTbar['SL']['THad'].p4.flatten(), GenTTbar['SL']['TLep'].p4.flatten(), GenTTbar['SL']['TTbar'].p4.flatten()
        thad_ctstar, tlep_ctstar = make_vars.ctstar_flat(thad_p4, tlep_p4)
        #set_trace()

        thad_pt_weights = self.reweighting['thad_pt'](thad_p4.pt)
        mtt_vs_thad_ctstar_weights = self.reweighting['mtt_vs_thad_ctstar'](ttbar_p4.mass, thad_ctstar)
        thad_pt_Interp_weights = self.reweighting['thad_pt_Interp'](thad_p4.pt)
        mtt_vs_thad_ctstar_Interp_weights = self.reweighting['mtt_vs_thad_ctstar_Interp'](ttbar_p4.mass, thad_ctstar)

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
            output['pt_thad'].fill(dataset=self.sample_name, rewt=rewt_type, pt=thad_p4.pt[el_mu_evts], weight=evt_weights[el_mu_evts])
            output['eta_thad'].fill(dataset=self.sample_name, rewt=rewt_type, eta=thad_p4.eta[el_mu_evts], weight=evt_weights[el_mu_evts])
            output['phi_thad'].fill(dataset=self.sample_name, rewt=rewt_type, phi=thad_p4.phi[el_mu_evts], weight=evt_weights[el_mu_evts])
            output['mass_thad'].fill(dataset=self.sample_name, rewt=rewt_type, mtop=thad_p4.mass[el_mu_evts], weight=evt_weights[el_mu_evts])
            output['ctstar_thad'].fill(dataset=self.sample_name, rewt=rewt_type, ctstar=thad_ctstar[el_mu_evts], weight=evt_weights[el_mu_evts])
            output['ctstar_abs_thad'].fill(dataset=self.sample_name, rewt=rewt_type, ctstar_abs=np.abs(thad_ctstar[el_mu_evts]), weight=evt_weights[el_mu_evts])

                    # tlep
            output['pt_tlep'].fill(dataset=self.sample_name, rewt=rewt_type, pt=tlep_p4.pt[el_mu_evts], weight=evt_weights[el_mu_evts])
            output['eta_tlep'].fill(dataset=self.sample_name, rewt=rewt_type, eta=tlep_p4.eta[el_mu_evts], weight=evt_weights[el_mu_evts])
            output['phi_tlep'].fill(dataset=self.sample_name, rewt=rewt_type, phi=tlep_p4.phi[el_mu_evts], weight=evt_weights[el_mu_evts])
            output['mass_tlep'].fill(dataset=self.sample_name, rewt=rewt_type, mtop=tlep_p4.mass[el_mu_evts], weight=evt_weights[el_mu_evts])
            output['ctstar_tlep'].fill(dataset=self.sample_name, rewt=rewt_type, ctstar=tlep_ctstar[el_mu_evts], weight=evt_weights[el_mu_evts])
            output['ctstar_abs_tlep'].fill(dataset=self.sample_name, rewt=rewt_type, ctstar_abs=np.abs(tlep_ctstar[el_mu_evts]), weight=evt_weights[el_mu_evts])

                    # ttbar system
            output['pt_tt'].fill(dataset=self.sample_name, rewt=rewt_type, pt=ttbar_p4.pt[el_mu_evts], weight=evt_weights[el_mu_evts])
            output['eta_tt'].fill(dataset=self.sample_name, rewt=rewt_type, eta=ttbar_p4.eta[el_mu_evts], weight=evt_weights[el_mu_evts])
            output['phi_tt'].fill(dataset=self.sample_name, rewt=rewt_type, phi=ttbar_p4.phi[el_mu_evts], weight=evt_weights[el_mu_evts])
            output['mtt'].fill(dataset=self.sample_name, rewt=rewt_type, mtt=ttbar_p4.mass[el_mu_evts], weight=evt_weights[el_mu_evts])
            output['mtt_vs_thad_ctstar'].fill(dataset=self.sample_name, rewt=rewt_type, mtt_2D=ttbar_p4.mass[el_mu_evts], ctstar_2D=thad_ctstar[el_mu_evts], weight=evt_weights[el_mu_evts])


        return output


    def postprocess(self, accumulator):
        return accumulator

#set_trace()
proc_executor = processor.iterative_executor if args.debug else processor.futures_executor
output = processor.run_uproot_job(fileset,
    treename='Events',
    processor_instance=Analyzer(),
    executor=proc_executor,
    executor_args={
        'workers': 8,
        'flatten' : True,
        'compression' : 5,
    },
    chunksize=10000 if args.debug else 100000,
)

if args.debug: print(output)
#print(output)

save(output, args.outfname)
print('%s has been written' % args.outfname)
