#import matplotlib.pyplot as plt
#plt.switch_backend('agg')
from coffea import hist
#from coffea.analysis_objects import JaggedCandidateArray
import coffea.processor as processor
import uproot
from pdb import set_trace
import os
from argparse import ArgumentParser
#from python.IDMuon import process_muons as proc_mus
#from python.IDElectron import process_electrons as proc_els
#from python.IDJet import process_jets as proc_jets
#import python.triggers as triggers
import python.ObjectSelection as objsel


parser = ArgumentParser()
parser.add_argument('sample', default='ttJets', help='Samples to run over')
parser.add_argument('--nfiles', default=-1, type=int, help='Specify the first number of files in the txt to run over. -1 means all')
args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'test_analyzer'

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

            ## make binning for hists
        mass_axis = hist.Bin("mass", "m [GeV]", 100, 0, 5)
        pt_axis = hist.Bin("pt", "p_{T} [GeV]", 200, 0, 1000)
        eta_axis = hist.Bin("eta", r"$\eta$", 200, -5, 5)
        phi_axis = hist.Bin("phi", r"$\phi$", 160, -4, 4)
        
            ## make dictionary of hists
        #dirs = ['muons']
        #dirs = ['tight_muons', 'one_tight_mu']
        dirs = ['muons', 'jets']
        #dirs = ['electrons']
        #dirs = ['muons', 'electrons']
        histo_dict = {}
        #for obj_dir in dirs:
        #    histo_dict['%s_mass' % obj_dir] = hist.Hist("Counts", mass_axis)
        #    histo_dict['%s_pt' % obj_dir] = hist.Hist("Counts", pt_axis)
        #    histo_dict['%s_eta' % obj_dir] = hist.Hist("Counts", eta_axis)
        #    histo_dict['%s_phi' % obj_dir] = hist.Hist("Counts", phi_axis)
        histo_dict['cutflow'] = processor.defaultdict_accumulator(int)

        self._accumulator = processor.dict_accumulator(histo_dict)
    
    @property
    def accumulator(self):
        return self._accumulator
    
    def process(self, df):
        output = self.accumulator.identity()

        #df['muons'] = proc_mus(df)
        #df['jets'] = proc_jets(df)
        #electrons = proc_els(df)

        ### muons selection
        passing_mus, output = objsel.select_muons(df, output)

        ### jets selection
        passing_jets, output = objsel.select_jets(df, output)

        ### final event selection
        passing_evt = passing_jets & passing_mus
        output['cutflow']['passing jet and mu'] += passing_evt.sum()

        #sel_mus = df['Muon'][(passing_evt)]
        #sel_jets = df['Jet'][(passing_evt)]

        #output = self.fill_hists(output, 'muons', sel_mus)        
        #output = self.fill_hists(output, 'jets', sel_jets)        
        
        #    ## electrons
        #output['cutflow']['all electrons'] += electrons.size
        ##trigelectrons = (electrons.trig > 0)
        ##electrons = electrons[trigelectrons]
        ##output['cutflow']['el trigger'] += trigelectrons.any().sum()

        #oneelectron = (electrons.counts == 1)
        #electrons = electrons[oneelectron]
        #output['cutflow']['1 electron'] += oneelectron.sum()

        #tight_el = (electrons.cutBasedId == 4)
        #electrons = electrons[tight_el]
        #output['cutflow']['el tight id'] += tight_el.any().sum()

        #ipcuts = (electrons.IPCuts > 0)
        #output['cutflow']['e IP cuts'] += ipcuts.any().sum()

        #output = self.fill_hists(output, 'electrons', electrons)        

        return output

    def fill_hists(self, accumulator, hdir, obj):
        accumulator['%s_mass' % hdir].fill(mass=obj.mass.flatten())
        accumulator['%s_pt' % hdir].fill(pt=obj.pt.flatten())
        accumulator['%s_eta' % hdir].fill(eta=obj.eta.flatten())
        accumulator['%s_phi' % hdir].fill(phi=obj.phi.flatten())

        return accumulator        

    def postprocess(self, accumulator):
        return accumulator

output = processor.run_uproot_job(fileset,
    treename='Events',
    processor_instance=Test_Analyzer(),
    executor=processor.iterative_executor,
    #executor=processor.futures_executor,
    #executor_args={'workers': 4, 'flatten' : True},
    executor_args={'workers': 4, 'flatten' : True, 'nano' : True},
    #chunksize=500000,
)

print(output)

    ## write hists to root file
if args.nfiles == -1:
    outdir = '/'.join([proj_dir, 'results', jobid])
    rfname = '%s/%s.root' % (outdir, args.sample)
else:
    outdir = proj_dir
    rfname = '%s/%s.test.%s.root' % (outdir, args.sample, analyzer)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

fout = uproot.recreate(rfname) if os.path.isfile(rfname) else uproot.create(rfname)
histos = [key for key in output.keys() if key != 'cutflow']
#set_trace()
for histo in histos:
    fout[histo] = hist.export1d(output[histo])
fout.close()
