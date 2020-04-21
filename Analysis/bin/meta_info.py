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

parser = ArgumentParser()
parser.add_argument('frange', type=str, help='Specify start:stop indices for files')
parser.add_argument('year', choices=['2016', '2017', '2018'], help='Specify which year to run over')
parser.add_argument('--sample', type=str, help='Use specific sample')
parser.add_argument('--fname', type=str, help='Specify output filename')
parser.add_argument('--debug', action='store_true', help='Uses iterative_executor for debugging purposes, otherwise futures_excutor will be used (faster)')

args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'meta_info'

    ## get samples to use
indir = '/'.join([proj_dir, 'inputs', '%s_%s' % (args.year, jobid)])
if args.sample:
        ## sample specified
    if not os.path.isfile('%s/%s.txt' % (indir, args.sample)):
        raise IOError("File with samples %s.txt not found" % args.sample)

    samples = [args.sample]
    if args.fname:
        print("  --- Sample name %s will be overridden by fname %s ---  \n" % (args.sample, args.fname))

else:
        ## get file that has names for all datasets to use
    fpath = '/'.join([indir, '%s_inputs.txt' % analyzer])
    if not os.path.isfile(fpath):
        raise IOError("File with samples %s_inputs.txt not found" % analyzer)
    
    txt_file = open(fpath, 'r')
    samples = [sample.strip('\n') for sample in txt_file if not sample.startswith('#')]
    if not samples:
        raise IOError("No samples found as inputs")

    ## add files to fileset
fileset = {}
for sample in samples:

    spath = '/'.join([indir, '%s.txt' % sample])
    if not os.path.isfile(spath):
        raise IOError("Sample file %s.txt not found" % sample)

    sfiles = open(spath, 'r')
    files_to_use = [fname.strip('\n') for fname in sfiles if not fname.startswith('#')]

    if ':' in args.frange:
        file_start, file_stop = int((args.frange).split(':')[0]), int((args.frange).split(':')[1])
    else:
        file_start = 0
        file_stop = len(files_to_use) if (args.frange).lower() == 'all' else int(args.frange)

    if file_start >= 0 and file_stop <= len(files_to_use):
        files_to_use = files_to_use[file_start:file_stop]
    else:
        raise IOError("The number of root files available for the %s sample is %i. args.frange must be less than or equal to this." % (sample, len(files_to_use) ) )

    #set_trace()
    fileset[sample] = files_to_use

# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class Meta_Analyzer(processor.ProcessorABC):
    def __init__(self, columns=[]):

        #if args.debug: set_trace()
        ## only get columns that are used
        self._columns = columns

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

    @property
    def columns(self):
        return self._columns

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

                    ## pick gen particles whose mother index == 0
                gps = genParts[(genParts.momIdx == 0)]
                tops = gps[(gps.pdgId == 6)]
                antitops = gps[(gps.pdgId == -6)]

                mtt = (tops+antitops).p4.mass.flatten()
                top_ctstar, tbar_ctstar = ctstar(tops.p4, antitops.p4)

                for idx in range(10):
                    output['mtt_topctstar'].fill(dataset=self.sample_name, evtIdx='%s' % idx, mtt=mtt[(events % 10) == idx], ctstar=top_ctstar[(events % 10) == idx])
                    output[self.sample_name]['nWeightedEvts_%s' % idx] += (genWeights[(events % 10) == idx] != 0).sum()

            if 'LHEPdfWeight' in df.columns:
                    ## check if there's the same number of pdf weights in every event
                pdf_wt_list = list(set(df.nLHEPdfWeight))
                if len(pdf_wt_list) > 1:
                    print(pdf_wt_list)
                    return output
                #set_trace()
                LHEpdfWeights = df.LHEPdfWeight
                    ## reshape because it's just a single fucking array instead of array of weights per event
                LHEpdfWeights = LHEpdfWeights.reshape((df.nLHEPdfWeight[0], events.size))
                    ## get sum of each pdf weight over all events
                sumLHEpdfWeights = LHEpdfWeights.sum(axis=1)
                output[self.sample_name]['sumLHEpdfWeights'] += sumLHEpdfWeights


        return output


    def postprocess(self, accumulator):
        return accumulator

# define columns to use (not working currently)
columns = [
    'event', 'dataset', 'Pileup_nTrueInt',
    'genWeight', 'nLHEPdfWeight', 'LHEPdfWeight',
    'nGenPart', 'GenPart_pt', 'GenPart_eta', 'GenPart_phi', 'GenPart_mass', 'GenPart_genPartIdxMother', 'GenPart_pdgId'
]

#set_trace()
proc_executor = processor.iterative_executor if args.debug else processor.futures_executor
output = processor.run_uproot_job(fileset,
    treename='Events',
    processor_instance=Meta_Analyzer(columns=columns),
    executor=proc_executor,
    executor_args={
        'workers': 8,
        #'workers': 4,
        'flatten' : True,
        'processor_compression' : 5,
    },
    chunksize=50000,
    #chunksize=500000,
)

if args.debug: print(output)
#print(output)

    ## save output to coffea pkl file
if (args.frange).lower() == 'all':
    outdir = '/'.join([proj_dir, 'results', '%s_%s' % (args.year, jobid), analyzer])
    if args.fname:
        cfname = '%s/%s.coffea' % (outdir, args.fname)
    elif args.sample:
        cfname = '%s/%s.coffea' % (outdir, args.sample)
    else:
        cfname = '%s/test_%s.coffea' % (outdir, analyzer)
else:
    if ':' in args.frange:
        outdir = '/'.join([proj_dir, 'results', '%s_%s' % (args.year, jobid), analyzer])
        if args.fname:
            cfname = '%s/%s_%sto%s.coffea' % (outdir, args.fname, file_start, file_stop)
        elif args.sample:
            cfname = '%s/%s_%sto%s.coffea' % (outdir, args.sample, file_start, file_stop)
        else:
            cfname = '%s/test_%sto%s.coffea' % (outdir, file_start, file_stop)
    else:
        outdir = proj_dir
        if args.fname:
            cfname = '%s/%s_%s.test.coffea' % (outdir, args.fname, analyzer)
        elif args.sample:
            cfname = '%s/%s_%s.test.coffea' % (outdir, args.sample, analyzer)
        else:
            cfname = '%s/%s.test.coffea' % (outdir, analyzer)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

save(output, cfname)
print('%s has been written' % cfname)

