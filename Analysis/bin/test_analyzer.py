from coffea import hist
from coffea.util import save, load
import coffea.processor as processor
from pdb import set_trace
import os
from argparse import ArgumentParser
import python.ObjectSelection as objsel
#import Utilities.maskedlazy as maskedlazy
import coffea.processor.dataframe
import itertools
#import python.Permutations as Permutations
#import python.MCWeights as MCWeights

parser = ArgumentParser()
parser.add_argument('sample', default='ttJets', help='Samples to run over')
parser.add_argument('frange', type=str, help='Specify start:stop indices for files')
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

#set_trace()
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
class Test_Analyzer(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
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
        self.jetmults = ['3Jets', '4PJets']
        directories = itertools.product(self.jetmults, self.leptypes, self.lepton.values())

        #set_trace()        
            ## make dictionary of hists
        histo_dict = {}
        for dirid in directories:
            tdir = '%s_%s%s' % dirid
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
        histo_dict['%s_pt' % tdir]    = hist.Hist("Counts", self.dataset_axis, self.pt_axis)
        histo_dict['%s_eta' % tdir]   = hist.Hist("Counts", self.dataset_axis, self.eta_axis)
        histo_dict['%s_phi' % tdir]   = hist.Hist("Counts", self.dataset_axis, self.phi_axis)
        histo_dict['%s_njets' % tdir] = hist.Hist("Counts", self.dataset_axis, self.njets_axis)

        return histo_dict

    
    def make_lep_hists(self, tdir):
        histo_dict = {}
        histo_dict['%s_pt' % tdir]  = hist.Hist("Counts", self.dataset_axis, self.pt_axis)
        histo_dict['%s_eta' % tdir] = hist.Hist("Counts", self.dataset_axis, self.eta_axis)
        histo_dict['%s_phi' % tdir] = hist.Hist("Counts", self.dataset_axis, self.phi_axis)

        return histo_dict

    def process(self, df):
        output = self.accumulator.identity()

        if not isinstance(df, coffea.processor.dataframe.LazyDataFrame):
            raise IOError("This function only works for LazyDataFrame objects")

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
                ## get clean jets
        clean_jets = sel_jets[~sel_jets.match(sel_leps, deltaRCut=0.4)] ## make sure no jets are within deltaR=0.4 of lepton
        three_jets_events = (clean_jets.counts == 3)
        fourPlus_jets_events = (clean_jets.counts > 3)

            ## only one lepton categorized as tight/loose
        tight_leps = sel_leps['TIGHT%s' % self.lepton[lep_to_use]].flatten()
        loose_leps = sel_leps['LOOSE%s' % self.lepton[lep_to_use]].flatten()

        #pref_weights = MCWeights.prefire_weight(df, mask=passing_evts) ## get nominal prefire weight for passing events
        #make_perms = Permutations.make_permutations(jets=clean_jets[tight_leps], leptons=sel_leps[tight_leps], MET=sel_met[tight_leps])
        #set_trace()

            ## fill hists for tight leptons
                ## 3 jets
        output = self.fill_jet_hists(output, '3Jets_TIGHT%s_Jets' % self.lepton[lep_to_use], clean_jets[(tight_leps & three_jets_events)])
        output = self.fill_lep_hists(output, '3Jets_TIGHT%s_%s' % (self.lepton[lep_to_use], lep_to_use), sel_leps[(tight_leps & three_jets_events)])
                ## 4+ jets
        output = self.fill_jet_hists(output, '4PJets_TIGHT%s_Jets' % self.lepton[lep_to_use], clean_jets[(tight_leps & fourPlus_jets_events)])
        output = self.fill_lep_hists(output, '4PJets_TIGHT%s_%s' % (self.lepton[lep_to_use], lep_to_use), sel_leps[(tight_leps & fourPlus_jets_events)])

            ## fill hists for loose leptons
                ## 3 jets
        output = self.fill_jet_hists(output, '3Jets_LOOSE%s_Jets' % self.lepton[lep_to_use], clean_jets[(loose_leps & three_jets_events)])
        output = self.fill_lep_hists(output, '3Jets_LOOSE%s_%s' % (self.lepton[lep_to_use], lep_to_use), sel_leps[(loose_leps & three_jets_events)])
                ## 4+ jets
        output = self.fill_jet_hists(output, '4PJets_LOOSE%s_Jets' % self.lepton[lep_to_use], clean_jets[(loose_leps & fourPlus_jets_events)])
        output = self.fill_lep_hists(output, '4PJets_LOOSE%s_%s' % (self.lepton[lep_to_use], lep_to_use), sel_leps[(loose_leps & fourPlus_jets_events)])

        return output

    def fill_jet_hists(self, accumulator, tdir, obj):
        accumulator['%s_pt' % tdir].fill(dataset=args.sample, pt=obj.pt.flatten())
        accumulator['%s_eta' % tdir].fill(dataset=args.sample, eta=obj.eta.flatten())
        accumulator['%s_phi' % tdir].fill(dataset=args.sample, phi=obj.phi.flatten())
        accumulator['%s_njets' % tdir].fill(dataset=args.sample, njets=obj.counts)

        return accumulator        

    def fill_lep_hists(self, accumulator, tdir, obj):
        #set_trace()
        accumulator['%s_pt' % tdir].fill(dataset=args.sample, pt=obj.pt.flatten())
        accumulator['%s_eta' % tdir].fill(dataset=args.sample, eta=obj.eta.flatten())
        accumulator['%s_phi' % tdir].fill(dataset=args.sample, phi=obj.phi.flatten())

        return accumulator        

    def postprocess(self, accumulator):
        return accumulator

proc_executor = processor.iterative_executor if args.debug else processor.futures_executor

#output = processor.run_spark_job(fileset,
output = processor.run_uproot_job(fileset,
    treename='Events',
    processor_instance=Test_Analyzer(),
    #executor=processor.spark_executor,
    #executor=processor.dask_executor,
    executor=proc_executor,
    executor_args={'workers': 6, 'flatten' : True},
    #chunksize=500000,
)

if args.debug:
    print(output)
#set_trace()
print(output['cutflow'])

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
    histos = [key for key in output.keys() if key != 'cutflow']
    #set_trace()
    for histo in histos:
        fout[histo] = hist.export1d(output[histo])
    fout.close()
    
    print('%s has been written' % rfname)
