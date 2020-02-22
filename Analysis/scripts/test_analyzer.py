from coffea import hist
import coffea.processor as processor
from pdb import set_trace
import os
from argparse import ArgumentParser
import python.ObjectSelection as objsel
#import Utilities.maskedlazy as maskedlazy
import coffea.processor.dataframe
import itertools
import python.Permutations as Permutations
#import python.MCWeights as MCWeights

parser = ArgumentParser()
parser.add_argument('sample', default='ttJets', help='Samples to run over')
parser.add_argument('--nfiles', default=-1, type=int, help='Specify the first number of files in the txt to run over. -1 means all')
parser.add_argument('--year', choices=['2016', '2017', '2018'], default=2016, help='Specify which year to run over')
parser.add_argument('--debug', action='store_true', help='Uses iterative_executor for debugging purposes, otherwise futures_excutor will be used (faster)')
parser.add_argument('--routput', action='store_true', help='Output (1D) histograms to root file. Only valid during debugging.')

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
        self.mass_axis = hist.Bin("mass", "m [GeV]", 100, 0, 5)
        self.pt_axis = hist.Bin("pt", "p_{T} [GeV]", 200, 0, 1000)
        self.eta_axis = hist.Bin("eta", r"$\eta$", 200, -5, 5)
        self.phi_axis = hist.Bin("phi", r"$\phi$", 160, -4, 4)
        self.njets_axis = hist.Bin("njets", "n_{jets}", 10, 0, 10)

        self.lepton = {
            'Muon': 'MU',
            #'Electron' : 'EL'
        }
        self.leptypes = ['LOOSE', 'TIGHT']
        #obj_dirs = [*self.lepton.keys()]+['Jets']
        directories = itertools.product(self.leptypes, self.lepton.values())

        #set_trace()        
            ## make dictionary of hists
        histo_dict = {}
        for dirid in directories:
            tdir = '%s%s' % dirid
    
                ## make jet hists
            jet_hists = self.make_jet_hists('%s_Jets' % tdir)
            histo_dict.update(jet_hists)
                ## make lepton hists
            lep_hists = self.make_lep_hists('%s_%s' % (tdir, [*self.lepton.keys()][0]))
            histo_dict.update(lep_hists)
        histo_dict['cutflow'] = processor.defaultdict_accumulator(int)

        self._accumulator = processor.dict_accumulator(histo_dict)
    
    @property
    def accumulator(self):
        return self._accumulator


    def make_jet_hists(self, tdir):
        histo_dict = {}
        histo_dict['%s_pt' % tdir] = hist.Hist("Counts", self.pt_axis)
        histo_dict['%s_eta' % tdir] = hist.Hist("Counts", self.eta_axis)
        histo_dict['%s_phi' % tdir] = hist.Hist("Counts", self.phi_axis)
        histo_dict['%s_njets' % tdir] = hist.Hist("Counts", self.njets_axis)

        return histo_dict

    
    def make_lep_hists(self, tdir):
        histo_dict = {}
        histo_dict['%s_pt' % tdir] = hist.Hist("Counts", self.pt_axis)
        histo_dict['%s_eta' % tdir] = hist.Hist("Counts", self.eta_axis)
        histo_dict['%s_phi' % tdir] = hist.Hist("Counts", self.phi_axis)

        return histo_dict

    def process(self, df):
        output = self.accumulator.identity()

        if not isinstance(df, coffea.processor.dataframe.LazyDataFrame):
            raise IOError("This function only works for LazyDataFrame objects")

        #set_trace()

        lep_to_use = [*self.lepton.keys()][0]
        presel_evts = objsel.select(df, leptype=lep_to_use, accumulator=output)
        ##df = maskedlazy.MaskedLazyDataFrame()
            
        output['cutflow']['nEvts passing jet and %s preselection' % lep_to_use] += presel_evts.sum()

            #btag reqs
        btag_wps = [col for col in df.Jet.columns if 'BTAG_' in col]
        if len(list(set(btag_wps))) == 1:
            btag_pass = (df.Jet[btag_wps[0]]).sum() >= 2
            if output: output['cutflow']['nEvts >=2 jets pass %s' % btag_wps[0]] += btag_pass.sum()
        else:
            raise IOError("Only 1 unique btag working point supported now")

        #set_trace()
        passing_evts = (presel_evts & btag_pass)
        output['cutflow']['nEvts passing jet and %s selection' % lep_to_use] += passing_evts.sum()

            ## get selected leptons, jets, and MET corresponding to passing events
        sel_leps = df[lep_to_use][(passing_evts)]
        sel_jets = df['Jet'][(passing_evts)]
        sel_met  = df['MET'][(passing_evts)]

            ## only one lepton categorized as tight/loose
        tight_leps = sel_leps['TIGHT%s' % self.lepton[lep_to_use]].flatten()
        loose_leps = sel_leps['LOOSE%s' % self.lepton[lep_to_use]].flatten()

        #set_trace()
        #pref_weights = MCWeights.prefire_weight(df, mask=passing_evts) ## get nominal prefire weight for passing events
        make_perms = Permutations.make_permutations(jets=sel_jets[tight_leps], leptons=sel_leps[tight_leps], MET=sel_met[tight_leps])

            ## fill hists for tight leptons
        output = self.fill_jet_hists(output, 'TIGHT%s_Jets' % self.lepton[lep_to_use], sel_jets[tight_leps])        
        output = self.fill_lep_hists(output, 'TIGHT%s_%s' % (self.lepton[lep_to_use], lep_to_use), sel_leps[tight_leps])        

            ## fill hists for loose leptons
        output = self.fill_jet_hists(output, 'LOOSE%s_Jets' % self.lepton[lep_to_use], sel_jets[loose_leps])        
        output = self.fill_lep_hists(output, 'LOOSE%s_%s' % (self.lepton[lep_to_use], lep_to_use), sel_leps[loose_leps])        

        return output

    def fill_jet_hists(self, accumulator, tdir, obj):
        #accumulator['%s_mass' % tdir].fill(mass=obj.mass.flatten())
        accumulator['%s_pt' % tdir].fill(pt=obj.pt.flatten())
        accumulator['%s_eta' % tdir].fill(eta=obj.eta.flatten())
        accumulator['%s_phi' % tdir].fill(phi=obj.phi.flatten())
        accumulator['%s_njets' % tdir].fill(njets=obj.counts)

        return accumulator        

    def fill_lep_hists(self, accumulator, tdir, obj):
        #set_trace()
        #accumulator['%s_mass' % tdir].fill(mass=obj.mass.flatten())
        accumulator['%s_pt' % tdir].fill(pt=obj.pt.flatten())
        accumulator['%s_eta' % tdir].fill(eta=obj.eta.flatten())
        accumulator['%s_phi' % tdir].fill(phi=obj.phi.flatten())

        return accumulator        

    def postprocess(self, accumulator):
        return accumulator

#set_trace()
proc_executor = processor.iterative_executor if args.debug else processor.futures_executor

#output = processor.run_spark_job(fileset,
output = processor.run_uproot_job(fileset,
    treename='Events',
    processor_instance=Test_Analyzer(),
    #executor=processor.spark_executor,
    #executor=processor.dask_executor,
    executor=proc_executor,
    #executor_args={'workers': 1, 'flatten' : True},
    executor_args={'workers': 4, 'flatten' : True},
    #chunksize=500000,
)

    ## save output to coffea pkl file
if args.nfiles == -1:
    outdir = '/'.join([proj_dir, 'results', jobid])
    cfname = '%s/%s.coffea' % (outdir, args.sample)
else:
    outdir = proj_dir
    cfname = '%s/%s.test.%s.coffea' % (outdir, args.sample, analyzer)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

save(output, cfname)

if args.debug: print(output)

if (args.debug and args.routput):
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
