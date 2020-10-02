from coffea import hist
import coffea.processor as processor
from pdb import set_trace
import os
from argparse import ArgumentParser
import coffea.processor.dataframe
import python.GenParticleSelector as genpsel
from Utilities.make_variables import ctstar as ctstar
from coffea.util import load, save
import numpy as np
import coffea.lumi_tools.lumi_tools as lumi_tools
import Utilities.prettyjson as prettyjson
import python.MCWeights as MCWeights

parser = ArgumentParser()
parser.add_argument('fset', type=str, help='Fileset dictionary (in string form) to be used for the processor')
parser.add_argument('year', choices=['2016', '2017', '2018'], help='Specify which year to run over')
parser.add_argument('outfname', type=str, help='Specify output filename, including directory and file extension')
parser.add_argument('--debug', action='store_true', help='Uses iterative_executor for debugging purposes, otherwise futures_excutor will be used (faster)')

args = parser.parse_args()

# convert input string of fileset dictionary to actual dictionary
fdict = (args.fset).replace("\'", "\"")
fileset = prettyjson.loads(fdict)

proj_dir = os.environ['PROJECT_DIR']

# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class Meta_Analyzer(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.evtIdx_axis = hist.Cat("evtIdx", "Event Index % 10")
        self.pu_axis = hist.Bin("pu", "nTrueInt", 200, 0, 200)
        self.mtt_axis = hist.Bin("mtt", r"$m_{t\\bart}$ [GeV]", 340, 300., 2000.)
        self.ctstar_axis = hist.Bin("ctstar", r"cos($\theta^{*}$)", 200, -1., 1.)

            ## make dictionary of hists
        histo_dict = {}
        histo_dict['mtt_topctstar'] = hist.Hist("mtt_vs_ctstar", self.dataset_axis, self.evtIdx_axis, self.mtt_axis, self.ctstar_axis)
        histo_dict['PUDistribution'] = hist.Hist("PUDistribution", self.dataset_axis, self.pu_axis)

        #set_trace()        
            ## construct dictionary of dictionaries to hold meta info for each sample
        #for sample in samples:
        for sample in fileset.keys():
            histo_dict[sample] = processor.defaultdict_accumulator(int)
            histo_dict['%s_runs_to_lumis' % sample] = processor.value_accumulator(list)

        self._accumulator = processor.dict_accumulator(histo_dict)
        self.sample_name = ''
        if args.year == '2016':
            self.Nominal_ttJets = ['ttJets_PS', 'ttJets']
        else:
            self.Nominal_ttJets = ['ttJetsSL', 'ttJetsHad', 'ttJetsDiLep']
    
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

        if self.sample_name.startswith('data_Single'):
            runs = df.run
            lumis = df.luminosityBlock
            Golden_Json_LumiMask = lumi_tools.LumiMask('%s/inputs/data/LumiMasks/%s_GoldenJson.txt' % (proj_dir, args.year))
            LumiMask = Golden_Json_LumiMask.__call__(runs, lumis) ## returns array of valid events

            output[self.sample_name]['nEvents'] += events[LumiMask].size
            output[self.sample_name]['nWeightedEvts'] += events[LumiMask].size

            if events[LumiMask].size > 0:
                valid_runs_lumis = np.unique(np.stack((runs[LumiMask], lumis[LumiMask]), axis=1), axis=0) ## make 2D array of uniqe valid [[run, lumi], [run, lumi]...] pairs
                    # make dictionary of valid runs: sorted list of unique lumisections for each valid run
                lumi_map = {str(valid_run):sorted(list(set(valid_runs_lumis[:, 1][valid_runs_lumis[:, 0] == valid_run]))) for valid_run in list(set(valid_runs_lumis[:, 0]))}

                output['%s_runs_to_lumis' % self.sample_name].add(list(lumi_map.items()))

        else:

            output['PUDistribution'].fill(dataset=self.sample_name, pu=df.Pileup_nTrueInt)

            output[self.sample_name]['nEvents'] += events.size

            genWeights = df.genWeight
            output[self.sample_name]['nWeightedEvts'] += (genWeights != 0).sum()
            output[self.sample_name]['sumGenWeights'] += genWeights.sum()

                ## create mtt vs cos theta* dists for nominal ttJets
            if self.sample_name in self.Nominal_ttJets:
                genParts = genpsel.process_genParts(df)

                #set_trace()
                    # same requirements as in GenParticleSelector
                tops = genParts[((genParts.status > 21) & (genParts.status < 30) & (genParts.momIdx != -1) & (genParts.pdgId == 6))]
                antitops = genParts[((genParts.status > 21) & (genParts.status < 30) & (genParts.momIdx != -1) & (genParts.pdgId == -6))]

                mtt = (tops+antitops).p4.mass.flatten()
                top_ctstar, tbar_ctstar = ctstar(tops.p4, antitops.p4)

                for idx in range(10):
                    output['mtt_topctstar'].fill(dataset=self.sample_name, evtIdx='%s' % idx, mtt=mtt[(events % 10) == idx], ctstar=top_ctstar[(events % 10) == idx])
                    output[self.sample_name]['nWeightedEvts_%s' % idx] += (genWeights[(events % 10) == idx] != 0).sum()

            if 'LHEPdfWeight' in df.columns:
                MCWeights.get_pdf_weights(df) # add pdf weights to df
                LHEpdfWeights = df.PDFWeights
                pdfwts = np.array(LHEpdfWeights) # convert to numpy so summing of each weight is easier
                    ## get sum of each pdf weight over all events
                sumLHEpdfWeights = pdfwts.sum()
                output[self.sample_name]['sumLHEpdfWeights'] += sumLHEpdfWeights


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
