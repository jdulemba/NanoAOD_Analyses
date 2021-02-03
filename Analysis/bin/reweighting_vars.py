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

# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class Meta_Analyzer(processor.ProcessorABC):
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
        self.sample_name = ''
    
    @property
    def accumulator(self):
        return self._accumulator

    def process(self, df):
        output = self.accumulator.identity()

        if not isinstance(df, coffea.processor.dataframe.LazyDataFrame):
            raise IOError("This function only works for LazyDataFrame objects")

        if args.debug: set_trace()
        events = df.event
        self.sample_name = df.dataset


        genWeights = df.genWeight

        # used for top pt reweighting
        GenTTbar = genpsel.select(df, mode='NORMAL')
        sl_evts = GenTTbar['SL']['TTbar'].counts > 0
        el_mu_evts = (np.abs(GenTTbar['SL']['Lepton'].pdgId) != 15).flatten()
        thad_pt = GenTTbar['SL']['THad'].p4.pt.flatten()
        output['thad_pt'].fill(dataset=self.sample_name, pt=thad_pt[el_mu_evts], weight=genWeights[sl_evts][el_mu_evts])

        #set_trace()
        thad_p4, tlep_p4 = GenTTbar['SL']['THad'].p4.flatten(), GenTTbar['SL']['TLep'].p4.flatten()
        thad_ctstar, tlep_ctstar = make_vars.ctstar_flat(thad_p4, tlep_p4)
        mtt = GenTTbar['SL']['TTbar'].p4.mass.flatten()
        output['mtt_vs_thad_ctstar'].fill(dataset=self.sample_name, mtt=mtt[el_mu_evts], ctstar=thad_ctstar[el_mu_evts], weight=genWeights[sl_evts][el_mu_evts])


        return output


    def postprocess(self, accumulator):
        return accumulator

#set_trace()
proc_executor = processor.iterative_executor if args.debug else processor.futures_executor
output = processor.run_uproot_job(fileset,
    treename='Events',
    processor_instance=Meta_Analyzer(),
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
