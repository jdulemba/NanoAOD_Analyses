from coffea import hist
import coffea.processor as processor
from pdb import set_trace
import os
from argparse import ArgumentParser
import coffea.processor.dataframe
#import itertools
#import python.MCWeights as MCWeights
import python.Partons as Partons

parser = ArgumentParser()
parser.add_argument('sample', default='ttJets', help='Samples to run over')
parser.add_argument('--nfiles', default=-1, type=int, help='Specify the first number of files in the txt to run over. -1 means all')
parser.add_argument('--year', choices=['2016', '2017', '2018'], default=2016, help='Specify which year to run over')
parser.add_argument('--debug', action='store_true', help='Uses iterative_executor for debugging purposes, otherwise futures_excutor will be used (faster)')

args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'get_meta_info'

sample = '/'.join([proj_dir, 'inputs', jobid, '%s.txt' % args.sample])
if not os.path.isfile(sample):
    raise IOError("Sample file %s.txt not found" % args.sample)

infiles = open(sample, 'r')
input_files = [fname.strip('\n') for fname in infiles]
if args.nfiles <= len(input_files):
    input_files = input_files[0:args.nfiles]
else:
    raise IOError("The number of root files available for the %s sample is %i. args.nfiles must be less than or equal to this." % (args.sample, len(input_files) ) )

fileset = {
    args.sample : input_files
}

# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class Test_Analyzer(processor.ProcessorABC):
    def __init__(self):

        #    ## make binning for hists
        self.pu_axis = hist.Bin("pu", "nTrueInt", 200, 0, 200)

        #set_trace()        
            ## make dictionary of hists
        histo_dict = {}
        histo_dict['PUDistribution'] = hist.Hist("PUDistribution", self.pu_axis)
        histo_dict['cutflow'] = processor.defaultdict_accumulator(int)

        self._accumulator = processor.dict_accumulator(histo_dict)
    
    @property
    def accumulator(self):
        return self._accumulator


    def process(self, df):
        output = self.accumulator.identity()

        if not isinstance(df, coffea.processor.dataframe.LazyDataFrame):
            raise IOError("This function only works for LazyDataFrame objects")

        #set_trace()

        ##df['GenPartons'] = Partons.process_genParts(df)
        #df['LHEPartons'] = Partons.process_lheParts(df)
        #lumiBlocks = df.luminosityBlock
        #runs = df.run
        #events = df.event
        #genWeights = df.genWeight

        #nom_LHEweight = df.LHEWeight_originalXWGTUP
        #eff_lumi = expected_evts/xsection


        output['PUDistribution'].fill(pu=df.Pileup_nTrueInt)

        return output


    #def fill_lep_hists(self, accumulator, tdir, obj):
    #    #set_trace()
    #    #accumulator['%s_mass' % tdir].fill(mass=obj.mass.flatten())
    #    accumulator['%s_pt' % tdir].fill(pt=obj.pt.flatten())
    #    accumulator['%s_eta' % tdir].fill(eta=obj.eta.flatten())
    #    accumulator['%s_phi' % tdir].fill(phi=obj.phi.flatten())

    #    return accumulator        

    def postprocess(self, accumulator):
        return accumulator


proc_executor = processor.iterative_executor if args.debug else processor.futures_executor
output = processor.run_uproot_job(fileset,
    treename='Events',
    processor_instance=Test_Analyzer(),
    executor=proc_executor,
    executor_args={'workers': 4, 'flatten' : True},
    #chunksize=500000,
)

#print(output)
#
    ## write hists to root file
if args.nfiles == -1:
    outdir = '/'.join([proj_dir, 'results', jobid])
    rfname = '%s/%s.root' % (outdir, args.sample)
else:
    outdir = proj_dir
    rfname = '%s/%s.test.%s.root' % (outdir, args.sample, analyzer)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

import uproot
fout = uproot.recreate(rfname) if os.path.isfile(rfname) else uproot.create(rfname)
histos = [key for key in output.keys() if key != 'cutflow']
#set_trace()
for histo in histos:
    fout[histo] = hist.export1d(output[histo])
fout.close()

print('%s has been written' % rfname)
