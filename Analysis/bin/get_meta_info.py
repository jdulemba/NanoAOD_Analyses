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

parser = ArgumentParser()
parser.add_argument('sample', default='ttJets', help='Samples to run over')
parser.add_argument('frange', type=str, help='Specify start:stop indices for files')
parser.add_argument('--year', choices=['2016', '2017', '2018'], default=2016, help='Specify which year to run over')
parser.add_argument('--debug', action='store_true', help='Uses iterative_executor for debugging purposes, otherwise futures_excutor will be used (faster)')
parser.add_argument('--routput', action='store_true', help='Output (1D) histograms to root file. Only valid during debugging.')

args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'get_meta_info'

if args.sample.startswith('data_Single'):
    raise IOError("Meta Info should only be run over simulation")

isNominalTTbar = True if (args.sample == 'ttJets_PS' or args.sample == 'ttJets') else False

sample = '/'.join([proj_dir, 'inputs', '%s_Testing' % args.year, '%s.txt' % args.sample])
#sample = '/'.join([proj_dir, 'inputs', jobid, '%s.txt' % args.sample])

if not os.path.isfile(sample):
    raise IOError("Sample file %s.txt not found" % args.sample)

infiles = open(sample, 'r')
input_files = [fname.strip('\n') for fname in infiles]

if ':' in args.frange:
    file_start, file_stop = int((args.frange).split(':')[0]), int((args.frange).split(':')[1])
else:
    file_start = 0
    file_stop = len(input_files)-1 if (args.frange).lower() == 'all' else int(args.frange)

if file_start >= 0 and file_stop <= len(input_files)-1:
    input_files = input_files[file_start:file_stop]
else:
    raise IOError("The number of root files available for the %s sample is %i. args.frange must be less than or equal to this." % (args.sample, len(input_files)-1 ) )

fileset = {
    args.sample : input_files
}

# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class Meta_Analyzer(processor.ProcessorABC):
    def __init__(self):

        #    ## make binning for hists
        self.evtIdx_axis = hist.Cat("evtIdx", "Event Index % 10")
        self.pu_axis = hist.Bin("pu", "nTrueInt", 200, 0, 200)
        self.mtt_axis = hist.Bin("mtt", "m_{tt}", 340, 300., 2000.)
        self.ctstar_axis = hist.Bin("ctstar", r"cos($\theta^{*}$)", 200, -1., 1.)

        #set_trace()        
            ## make dictionary of hists
        histo_dict = {}
        if isNominalTTbar:
            for idx in range(10):
                histo_dict['mtt_idx%i' % idx] = hist.Hist("mtt_idx%i" % idx, self.evtIdx_axis, self.mtt_axis)
                histo_dict['top_ctstar_idx%i' % idx] = hist.Hist("ctstar_idx%i" % idx, self.evtIdx_axis, self.ctstar_axis)
                histo_dict['mtt_topctstar_idx%i' % idx] = hist.Hist("mtt_vs_ctstar_idx%i" % idx, self.mtt_axis, self.ctstar_axis)

        histo_dict['PUDistribution'] = hist.Hist("PUDistribution", self.pu_axis)
        histo_dict['MetaInfo'] = processor.defaultdict_accumulator(int)

        self._accumulator = processor.dict_accumulator(histo_dict)
    
    @property
    def accumulator(self):
        return self._accumulator


    def process(self, df):
        output = self.accumulator.identity()

        if not isinstance(df, coffea.processor.dataframe.LazyDataFrame):
            raise IOError("This function only works for LazyDataFrame objects")


        events = df.event
        if isNominalTTbar:
            genParts = Partons.process_genParts(df)

                ## pick gen particles whose mother index == 0
            gps = genParts[(genParts.momIdx == 0)]
            tops = gps[(gps.pdgId == 6)]
            antitops = gps[(gps.pdgId == -6)]

            mtt = (tops+antitops).p4.mass.flatten()
            top_ctstar, tbar_ctstar = ctstar(tops.p4, antitops.p4)
            #set_trace()

            for idx in range(10):
                output['mtt_idx%i' % idx].fill(evtIdx='Idx %% 10 == %s' % idx, mtt=mtt[(events % 10) == idx])
                output['top_ctstar_idx%i' % idx].fill(evtIdx='Idx %% 10 == %s' % idx, ctstar=top_ctstar[(events % 10) == idx])
                output['mtt_topctstar_idx%i' % idx].fill(mtt=mtt[(events % 10) == idx], ctstar=top_ctstar[(events % 10) == idx])

        #lumiBlocks = df.luminosityBlock
        #runs = df.run
        #set_trace()
        genWeights = df.genWeight
        output['MetaInfo']['sumWeights'] += genWeights.sum()
        #nom_LHEweight = df.LHEWeight_originalXWGTUP
        #eff_lumi = expected_evts/xsection


        output['PUDistribution'].fill(pu=df.Pileup_nTrueInt)

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
    cfname = '%s/%s.coffea' % (outdir, args.sample)
    if (args.debug and args.routput):
        rfname = '%s/%s.root' % (outdir, args.sample)
else:
    if ':' in args.frange:
        outdir = '/'.join([proj_dir, 'results', jobid, analyzer])
        cfname = '%s/%s_%sto%s.coffea' % (outdir, args.sample, file_start, file_stop)
        if (args.debug and args.routput):
            rfname = '%s/%s_%sto%s.root' % (outdir, args.sample, file_start, file_stop)
    else:
        outdir = proj_dir
        cfname = '%s/%s.test.%s.coffea' % (outdir, args.sample, analyzer)
        if (args.debug and args.routput):
            rfname = '%s/%s.test.%s.root' % (outdir, args.sample, analyzer)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

save(output, cfname)
print('%s has been written' % cfname)

if (args.debug and args.routput):
        ## write hists to root file
    import uproot
    fout = uproot.recreate(rfname) if os.path.isfile(rfname) else uproot.create(rfname)
    histos = [key for key in output.keys() if key != 'MetaInfo']
    #set_trace()
    for histo in histos:
        if output[histo].dense_dim() == 1:
            fout[histo] = hist.export1d(output[histo])
            print(histo, ' written to file')
        #elif output[histo].dense_dim() == 2:
        #    #continue
        #    fout[histo] = hist.export2d(output[histo])
        else:
            print('Only 1D histograms supproted for writing, %s skipped' % histo)
    fout.close()
    
    print('%s has been written' % rfname)
