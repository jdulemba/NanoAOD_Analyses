#!/usr/bin/env python

from coffea import hist
from coffea.util import save, load
import coffea.processor as processor
from pdb import set_trace
import os, sys
import python.ObjectSelection as objsel
import coffea.processor.dataframe
import Utilities.plot_tools as plt_tools
import python.MCWeights as MCWeights
import numpy as np
import Utilities.prettyjson as prettyjson
import coffea.lumi_tools.lumi_tools as lumi_tools
import python.GenParticleSelector as genpsel

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'GenPartons'

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('fset', type=str, help='Fileset dictionary (in string form) to be used for the processor')
parser.add_argument('year', choices=['2016', '2017', '2018'], help='Specify which year to run over')
parser.add_argument('outfname', type=str, help='Specify output filename, including directory and file extension')
parser.add_argument('--debug', action='store_true', help='Uses iterative_executor for debugging purposes, otherwise futures_excutor will be used (faster)')

args = parser.parse_args()

# convert input string of fileset dictionary to actual dictionary
fdict = (args.fset).replace("\'", "\"")
fileset = prettyjson.loads(fdict)

## load corrections for event weights
pu_correction = load('%s/Corrections/%s/MC_PU_Weights.coffea' % (proj_dir, jobid))
lepSF_correction = load('%s/Corrections/leptonSFs.coffea' % proj_dir)
jet_corrections = load('%s/Corrections/JetCorrections.coffea' % proj_dir)[args.year]
corrections = {
    'Pileup' : pu_correction,
    'Prefire' : True,
    'LeptonSF' : lepSF_correction,
    'BTagSF' : False,
    'JetCor' : jet_corrections,
}

if corrections['BTagSF'] == True:
    import python.BTagScaleFactors as btagSF
    threejets_btagSFs = btagSF.create_btag_sf_computer(args.year, '3')
    fourPlusjets_btagSFs = btagSF.create_btag_sf_computer(args.year, '4+')


# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class GenPartons(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.object_axis = hist.Cat("objtype", "Gen Object")
        self.ttdecaymode_axis = hist.Cat("ttdecay", "tt Decay Mode")
        self.wdecaymode_axis = hist.Cat("wdecay", "W Decay Mode")
        self.pt_axis = hist.Bin("pt", "p_{T} [GeV]", 200, 0, 1000)
        self.eta_axis = hist.Bin("eta", r"$\eta$", 60, -3, 3)
        self.phi_axis = hist.Bin("phi", r"$\phi$", 160, -4, 4)
        self.mass_axis = hist.Bin("mass", "mass [GeV]", 1000, 0, 1000)
        self.energy_axis = hist.Bin("energy", "E [GeV]", 200, 0, 1000)

            ## make dictionary of hists
        histo_dict = {}
                ## make jet hists
        gen_hists = self.genpartons_hists()
        histo_dict.update(gen_hists)

        histo_dict['cutflow'] = processor.defaultdict_accumulator(int)

        self._accumulator = processor.dict_accumulator(histo_dict)
        self.sample_name = ''
        self.corrections = corrections
        if args.year == '2016':
            self.Nominal_ttJets = ['ttJets_PS', 'ttJets']
        else:
            self.Nominal_ttJets = ['ttJetsSL', 'ttJetsHad', 'ttJetsDiLep']
        if not self.Nominal_ttJets:
            raise ValueError("This should only be run on ttbar events!")

    
    @property
    def accumulator(self):
        return self._accumulator


    def genpartons_hists(self):
        histo_dict = {}
        histo_dict['pt']    = hist.Hist("Events", self.dataset_axis, self.object_axis, self.ttdecaymode_axis, self.pt_axis)
        histo_dict['eta']   = hist.Hist("Events", self.dataset_axis, self.object_axis, self.ttdecaymode_axis, self.eta_axis)
        histo_dict['phi']   = hist.Hist("Events", self.dataset_axis, self.object_axis, self.ttdecaymode_axis, self.phi_axis)
        histo_dict['mass']  = hist.Hist("Events", self.dataset_axis, self.object_axis, self.ttdecaymode_axis, self.mass_axis)
        histo_dict['energy']= hist.Hist("Events", self.dataset_axis, self.object_axis, self.ttdecaymode_axis, self.energy_axis)

        #histo_dict['pt']    = hist.Hist("Events", self.dataset_axis, self.objtype_axis, self.ttdecaymode_axis, self.wdecaymode_axis, self.pt_axis)
        #histo_dict['eta']   = hist.Hist("Events", self.dataset_axis, self.objtype_axis, self.ttdecaymode_axis, self.wdecaymode_axis, self.eta_axis)
        #histo_dict['phi']   = hist.Hist("Events", self.dataset_axis, self.objtype_axis, self.ttdecaymode_axis, self.wdecaymode_axis, self.phi_axis)
        #histo_dict['mass']  = hist.Hist("Events", self.dataset_axis, self.objtype_axis, self.ttdecaymode_axis, self.wdecaymode_axis, self.mass_axis)
        #histo_dict['energy']= hist.Hist("Events", self.dataset_axis, self.objtype_axis, self.ttdecaymode_axis, self.wdecaymode_axis, self.energy_axis)
        return histo_dict



    def process(self, df):
        np.random.seed(10) # sets seed so values from random distributions are reproducible (JER corrections)
        output = self.accumulator.identity()

        if not isinstance(df, coffea.processor.dataframe.LazyDataFrame):
            raise IOError("This function only works for LazyDataFrame objects")

        #if args.debug: set_trace()
        self.sample_name = df.dataset

            ## make event weights
                # data or MC distinction made internally
        evt_weights = MCWeights.get_event_weights(df, year=args.year, corrections=self.corrections)


        GenTTbar = genpsel.select(df, systype='FINAL', mode='NORMAL')
        #GenTTbar = genpsel.select(df, systype='', mode='NORMAL')
        #genpsel.select(df, mode='LHE')

        #set_trace()
        for ttdecay in GenTTbar.columns:
            evt_weights_to_use = evt_weights.weight()[GenTTbar[ttdecay]['TTbar'].counts > 0]
            for gen_obj in GenTTbar[ttdecay].columns:
                output = self.fill_genp_hists(accumulator=output, genp_type=gen_obj, ttdecaymode=ttdecay, obj=GenTTbar[ttdecay][gen_obj], evt_weights=evt_weights_to_use)
                #output = self.fill_genp_hists(accumulator=output, genp_type=gen_obj, ttdecaymode=ttdecay, wdecaymode=, obj=GenTTbar[ttdecay][gen_obj], evt_weights=evt_weights_to_use)
                

        return output

    def fill_genp_hists(self, accumulator, genp_type, ttdecaymode, obj, evt_weights):
    #def fill_genp_hists(self, accumulator, genp_type, ttdecaymode, wdecaymode, obj, evt_weights):
        #set_trace()
        accumulator['pt'].fill(    dataset=self.sample_name, objtype=genp_type, ttdecay=ttdecaymode, pt=obj.p4.pt.flatten(), weight=evt_weights)
        accumulator['eta'].fill(   dataset=self.sample_name, objtype=genp_type, ttdecay=ttdecaymode, eta=obj.p4.eta.flatten(), weight=evt_weights)
        accumulator['phi'].fill(   dataset=self.sample_name, objtype=genp_type, ttdecay=ttdecaymode, phi=obj.p4.phi.flatten(), weight=evt_weights)
        accumulator['mass'].fill(  dataset=self.sample_name, objtype=genp_type, ttdecay=ttdecaymode, mass=obj.p4.mass.flatten(), weight=evt_weights)
        accumulator['energy'].fill(dataset=self.sample_name, objtype=genp_type, ttdecay=ttdecaymode, energy=obj.p4.E.flatten(), weight=evt_weights)

        return accumulator        


    def postprocess(self, accumulator):
        return accumulator

proc_executor = processor.iterative_executor if args.debug else processor.futures_executor

output = processor.run_uproot_job(fileset,
    treename='Events',
    processor_instance=GenPartons(),
    executor=proc_executor,
    executor_args={
        'workers': 8,
        'flatten' : True,
        'compression': 5,
    },
    chunksize=10000,
    #chunksize=500000,
)


if args.debug:
    print(output)
#set_trace()
#print(output['cutflow'])

save(output, args.outfname)
print('%s has been written' % args.outfname)
