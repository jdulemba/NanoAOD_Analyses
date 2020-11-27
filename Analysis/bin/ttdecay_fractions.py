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

is_ttJets_ = samplename in ['ttJetsSL', 'ttJetsHad', 'ttJetsDiLep', 'ttJets_PS']
if not is_ttJets_:
    raise ValueError("This should only be run on nominal ttbar events!")

    ## get list of signal samples to run over
# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class ttdecay_fractions(processor.ProcessorABC):
    def __init__(self):


        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.ttdecay_axis = hist.Cat("ttdecay", "tt Decay Type")

            ## make dictionary of hists
        histo_dict = {}
        histo_dict['tt_decay'] = hist.Hist("Events", self.dataset_axis, self.ttdecay_axis)

        self._accumulator = processor.dict_accumulator(histo_dict)

    
    @property
    def accumulator(self):
        return self._accumulator


    def process(self, df):
        np.random.seed(10) # sets seed so values from random distributions are reproducible (JER corrections)
        output = self.accumulator.identity()

        GenTTbar = genpsel.select(df, mode='NORMAL')
            # semileptonic decays
        output['tt_decay'].fill(dataset=df.dataset, ttdecay='SL e', weight=np.ones(1)*((np.abs(GenTTbar['SL']['Lepton'].pdgId) == 11).sum().sum()))
        output['tt_decay'].fill(dataset=df.dataset, ttdecay='SL mu', weight=np.ones(1)*((np.abs(GenTTbar['SL']['Lepton'].pdgId) == 13).sum().sum()))
        output['tt_decay'].fill(dataset=df.dataset, ttdecay='SL tau', weight=np.ones(1)*((np.abs(GenTTbar['SL']['Lepton'].pdgId) == 15).sum().sum()))
        output['tt_decay'].fill(dataset=df.dataset, ttdecay='SL tau->l', weight=np.ones(1)*((GenTTbar['SL']['Lepton'].decaytype == 1).sum().sum()))
        output['tt_decay'].fill(dataset=df.dataset, ttdecay='SL tau->h', weight=np.ones(1)*((GenTTbar['SL']['Lepton'].decaytype == 2).sum().sum()))
        output['tt_decay'].fill(dataset=df.dataset, ttdecay='SL Total', weight=np.ones(1)*(GenTTbar['SL']['Lepton'].counts.sum()))

            # dileptonic decays
        output['tt_decay'].fill(dataset=df.dataset, ttdecay='DL e e', weight=np.ones(1)*(((np.abs(GenTTbar['DL']['First_plus'].pdgId) == 11) & (np.abs(GenTTbar['DL']['First_minus'].pdgId) == 11)).sum().sum()))
        output['tt_decay'].fill(dataset=df.dataset, ttdecay='DL e mu', weight=np.ones(1)*((((np.abs(GenTTbar['DL']['First_plus'].pdgId) == 11) & (np.abs(GenTTbar['DL']['First_minus'].pdgId) == 13)) | ((np.abs(GenTTbar['DL']['First_plus'].pdgId) == 13) & (np.abs(GenTTbar['DL']['First_minus'].pdgId) == 11))).sum().sum()))
        output['tt_decay'].fill(dataset=df.dataset, ttdecay='DL e tau->l', weight=np.ones(1)*((((np.abs(GenTTbar['DL']['First_plus'].pdgId) == 11) & ((np.abs(GenTTbar['DL']['First_minus'].pdgId) == 15) & (GenTTbar['DL']['First_minus'].decaytype == 1))) | ((np.abs(GenTTbar['DL']['First_minus'].pdgId) == 11) & ((np.abs(GenTTbar['DL']['First_plus'].pdgId) == 15) & (GenTTbar['DL']['First_plus'].decaytype == 1)))).sum().sum()))
        output['tt_decay'].fill(dataset=df.dataset, ttdecay='DL e tau->h', weight=np.ones(1)*((((np.abs(GenTTbar['DL']['First_plus'].pdgId) == 11) & ((np.abs(GenTTbar['DL']['First_minus'].pdgId) == 15) & (GenTTbar['DL']['First_minus'].decaytype == 2))) | ((np.abs(GenTTbar['DL']['First_minus'].pdgId) == 11) & ((np.abs(GenTTbar['DL']['First_plus'].pdgId) == 15) & (GenTTbar['DL']['First_plus'].decaytype == 2)))).sum().sum()))
        output['tt_decay'].fill(dataset=df.dataset, ttdecay='DL mu mu', weight=np.ones(1)*(((np.abs(GenTTbar['DL']['First_plus'].pdgId) == 13) & (np.abs(GenTTbar['DL']['First_minus'].pdgId) == 13)).sum().sum()))
        output['tt_decay'].fill(dataset=df.dataset, ttdecay='DL mu tau->l', weight=np.ones(1)*((((np.abs(GenTTbar['DL']['First_plus'].pdgId) == 13) & ((np.abs(GenTTbar['DL']['First_minus'].pdgId) == 15) & (GenTTbar['DL']['First_minus'].decaytype == 1))) | ((np.abs(GenTTbar['DL']['First_minus'].pdgId) == 13) & ((np.abs(GenTTbar['DL']['First_plus'].pdgId) == 15) & (GenTTbar['DL']['First_plus'].decaytype == 1)))).sum().sum()))
        output['tt_decay'].fill(dataset=df.dataset, ttdecay='DL mu tau->h', weight=np.ones(1)*((((np.abs(GenTTbar['DL']['First_plus'].pdgId) == 13) & ((np.abs(GenTTbar['DL']['First_minus'].pdgId) == 15) & (GenTTbar['DL']['First_minus'].decaytype == 2))) | ((np.abs(GenTTbar['DL']['First_minus'].pdgId) == 13) & ((np.abs(GenTTbar['DL']['First_plus'].pdgId) == 15) & (GenTTbar['DL']['First_plus'].decaytype == 2)))).sum().sum()))
        output['tt_decay'].fill(dataset=df.dataset, ttdecay='DL tau tau->ll', weight=np.ones(1)*((((np.abs(GenTTbar['DL']['First_plus'].pdgId) == 15) & (GenTTbar['DL']['First_plus'].decaytype == 1)) & ((np.abs(GenTTbar['DL']['First_minus'].pdgId) == 15) & (GenTTbar['DL']['First_minus'].decaytype == 1))).sum().sum()))
        output['tt_decay'].fill(dataset=df.dataset, ttdecay='DL tau tau->lh', weight=np.ones(1)*(((((np.abs(GenTTbar['DL']['First_plus'].pdgId) == 15) & (GenTTbar['DL']['First_plus'].decaytype == 1)) & ((np.abs(GenTTbar['DL']['First_minus'].pdgId) == 15) & (GenTTbar['DL']['First_minus'].decaytype == 2))) | (((np.abs(GenTTbar['DL']['First_plus'].pdgId) == 15) & (GenTTbar['DL']['First_plus'].decaytype == 2)) & ((np.abs(GenTTbar['DL']['First_minus'].pdgId) == 15) & (GenTTbar['DL']['First_minus'].decaytype == 1)))).sum().sum()))
        output['tt_decay'].fill(dataset=df.dataset, ttdecay='DL tau tau->hh', weight=np.ones(1)*((((np.abs(GenTTbar['DL']['First_plus'].pdgId) == 15) & (GenTTbar['DL']['First_plus'].decaytype == 2)) & ((np.abs(GenTTbar['DL']['First_minus'].pdgId) == 15) & (GenTTbar['DL']['First_minus'].decaytype == 2))).sum().sum()))
        output['tt_decay'].fill(dataset=df.dataset, ttdecay='DL Total', weight=np.ones(1)*(GenTTbar['DL']['First_plus'].counts.sum()))

            # all hadronic decays
        output['tt_decay'].fill(dataset=df.dataset, ttdecay='Had Total', weight=np.ones(1)*(GenTTbar['Had']['First_plus'].counts.sum()))

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
