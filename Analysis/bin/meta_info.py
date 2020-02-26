from coffea import hist
import coffea.processor as processor
from pdb import set_trace
import os
from argparse import ArgumentParser
import coffea.processor.dataframe
#import itertools
#import python.MCWeights as MCWeights
import python.Partons as Partons
from Utilities.make_variables import ctstar as ctstar
from coffea.util import load, save
import Utilities.plot_tools as plt_tools
import numpy as np

parser = ArgumentParser()
parser.add_argument('frange', type=str, help='Specify start:stop indices for files')
parser.add_argument('--year', choices=['2016', '2017', '2018'], default=2016, help='Specify which year to run over')
parser.add_argument('--debug', action='store_true', help='Uses iterative_executor for debugging purposes, otherwise futures_excutor will be used (faster)')

args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'meta_info'

    ## get file that has names for all datasets to use
fpath = '/'.join([proj_dir, 'inputs', jobid, '%s_inputs.txt' % analyzer])
if not os.path.isfile(fpath):
    raise IOError("File with samples %s_inputs.txt not found" % analyzer)

txt_file = open(fpath, 'r')
samples = [sample.strip('\n') for sample in txt_file if not sample.startswith('#')]
if not samples:
    raise IOError("No samples found as inputs")

#if ('ttJets_PS' in samples and 'ttJets' in samples):
#    raise IOError("ttJets_PS and ttJets samples can't be run at the same time because double counting for the mtt_cth distributions will occur.")

    ## add files to fileset
fileset = {}
for sample in samples:
    if sample.startswith('data_Single'):
        raise IOError("Meta Info should only be run over simulation")

    spath = '/'.join([proj_dir, 'inputs', jobid, '%s.txt' % sample])
    if not os.path.isfile(spath):
        raise IOError("Sample file %s.txt not found" % sample)

    sfiles = open(spath, 'r')
    files_to_use = [fname.strip('\n') for fname in sfiles]

    #set_trace()
    if ':' in args.frange:
        file_start, file_stop = int((args.frange).split(':')[0]), int((args.frange).split(':')[1])
    else:
        file_start = 0
        file_stop = len(files_to_use)-1 if (args.frange).lower() == 'all' else int(args.frange)

    if file_start >= 0 and file_stop <= len(files_to_use)-1:
        files_to_use = files_to_use[file_start:file_stop]
    else:
        raise IOError("The number of root files available for the %s sample is %i. args.frange must be less than or equal to this." % (sample, len(files_to_use) ) )

    fileset[sample] = files_to_use


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
        for sample in samples:
            histo_dict[sample] = processor.defaultdict_accumulator(int)

        self._accumulator = processor.dict_accumulator(histo_dict)
        self.sample_name = ''
    
    @property
    def accumulator(self):
        return self._accumulator


    def process(self, df):
        output = self.accumulator.identity()

        if not isinstance(df, coffea.processor.dataframe.LazyDataFrame):
            raise IOError("This function only works for LazyDataFrame objects")

        events = df.event
        self.sample_name = df.dataset

        output['PUDistribution'].fill(dataset=self.sample_name, pu=df.Pileup_nTrueInt)

        if args.debug: set_trace()
        output[self.sample_name]['nEvents'] += events.size

        genWeights = df.genWeight
        output[self.sample_name]['nWeightedEvts'] += (genWeights != 0).sum()
        output[self.sample_name]['sumGenWeights'] += genWeights.sum()

            ## check if there's the same number of pdf weights in every event
        if not np.equal(df.nLHEPdfWeight, df.nLHEPdfWeight[0]).all():
            raise IOError("Events don't have the same number of LHE PDF weights!")
        LHEpdfWeights = df.LHEPdfWeight
            ## reshape because it's just a single fucking array instead of array of weights per event
        LHEpdfWeights = LHEpdfWeights.reshape((df.nLHEPdfWeight[0], events.size))
            ## get sum of each pdf weight over all events
        sumLHEpdfWeights = LHEpdfWeights.sum(axis=1)
        output[self.sample_name]['sumLHEpdfWeights'] += sumLHEpdfWeights


            ## create mtt vs cos theta* dists for nominal ttJets
        if (self.sample_name == 'ttJets_PS' or self.sample_name == 'ttJets'):
            genParts = Partons.process_genParts(df)

                ## pick gen particles whose mother index == 0
            gps = genParts[(genParts.momIdx == 0)]
            tops = gps[(gps.pdgId == 6)]
            antitops = gps[(gps.pdgId == -6)]

            mtt = (tops+antitops).p4.mass.flatten()
            top_ctstar, tbar_ctstar = ctstar(tops.p4, antitops.p4)
            #set_trace()

            for idx in range(10):
                output['mtt_topctstar'].fill(dataset=self.sample_name, evtIdx='%s' % idx, mtt=mtt[(events % 10) == idx], ctstar=top_ctstar[(events % 10) == idx])

        #lumiBlocks = df.luminosityBlock
        #runs = df.run
        #set_trace()
        #nom_LHEweight = df.LHEWeight_originalXWGTUP
        #eff_lumi = expected_evts/xsection

        return output


    def postprocess(self, accumulator):
        return accumulator


proc_executor = processor.iterative_executor if args.debug else processor.futures_executor
output = processor.run_uproot_job(fileset,
    treename='Events',
    processor_instance=Meta_Analyzer(),
    executor=proc_executor,
    executor_args={'workers': 4, 'flatten' : True},
    #chunksize=500000,
)

if args.debug: print(output)

    ## save output to coffea pkl file
if (args.frange).lower() == 'all':
    outdir = '/'.join([proj_dir, 'results', jobid, analyzer])
    cfname = '%s/%s.coffea' % (outdir, 'test')
    #cfname = '%s/%s.coffea' % (outdir, args.sample)
else:
    if ':' in args.frange:
        outdir = '/'.join([proj_dir, 'results', jobid, analyzer])
        cfname = '%s/%sto%s.coffea' % (outdir, file_start, file_stop)
        #cfname = '%s/%s_%sto%s.coffea' % (outdir, args.sample, file_start, file_stop)
    else:
        outdir = proj_dir
        cfname = '%s/%s.test.coffea' % (outdir, analyzer)
        #cfname = '%s/%s.test.%s.coffea' % (outdir, args.sample, analyzer)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

save(output, cfname)
print('%s has been written' % cfname)

