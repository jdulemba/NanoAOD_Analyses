from coffea import hist
from coffea.util import save, load
import coffea.processor as processor
from pdb import set_trace
import os
from argparse import ArgumentParser
import python.ObjectSelection as objsel
#import Utilities.maskedlazy as maskedlazy
import coffea.processor.dataframe
from coffea.arrays import Initialize
import itertools
import Utilities.plot_tools as plt_tools
#import python.LeptonSF as lepSF
#import python.BTagScaleFactors as btagSF
#import python.Permutations as Permutations
import python.MCWeights as MCWeights
import numpy as np
import Utilities.prettyjson as prettyjson

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'test_analyzer'

parser = ArgumentParser()
parser.add_argument('frange', type=str, help='Specify start:stop indices for files')
parser.add_argument('year', choices=['2016', '2017', '2018'], help='Specify which year to run over')
parser.add_argument('--debug', action='store_true', help='Uses iterative_executor for debugging purposes, otherwise futures_excutor will be used (faster)')

args = parser.parse_args()


    ## get file that has names for all datasets to use
fpath = '/'.join([proj_dir, 'inputs', '%s_%s' % (args.year, jobid), '%s_inputs.txt' % analyzer])
if not os.path.isfile(fpath):
    raise IOError("File with samples %s_inputs.txt not found" % analyzer)

txt_file = open(fpath, 'r')
samples = [sample.strip('\n') for sample in txt_file if not sample.startswith('#')]
if not samples:
    raise IOError("No samples found as inputs")

#    ## load json files for data lumi and cross sections -- testing
#data_lumi =  prettyjson.loads(open('%s/inputs/lumis.json' % proj_dir).read())[args.year]
#samples_file = prettyjson.loads(open('%s/inputs/samples_%s.json' % (proj_dir, args.year)).read())
#plt_weights = {}

    ## add files to fileset and get plotting weights
fileset = {}
for sample in samples:
    txtpath = '/'.join([proj_dir, 'inputs', '%s_%s' % (args.year, jobid), '%s.txt' % sample])
    if not os.path.isfile(txtpath):
        raise IOError("Sample file %s.txt not found" % sample)

    #metafile = '%s/inputs/%s/%s.meta.json' % (proj_dir, jobid, sample)
    #if os.path.isfile(metafile):
    #    nWeightedEvts = prettyjson.loads(open(metafile).read())["nWeightedEvts"]
    #    tmp_sample_name = 'ttJets' if sample == 'TEST_ttJets' else sample ## for testing
    #    xsec = [info['xsection'] for info in samples_file if info['name'] == tmp_sample_name ][0] ## for testing
    #    #xsec = [info['xsection'] for info in samples_file if info['name'] == sample ][0]
    #    plt_weights[sample] = data_lumi/(nWeightedEvts/xsec)
    
    txtfiles = open(txtpath, 'r')
    files_to_use = [fname.strip('\n') for fname in txtfiles]
    
    #set_trace()
    if ':' in args.frange:
        file_start, file_stop = int((args.frange).split(':')[0]), int((args.frange).split(':')[1])
    else:
        file_start = 0
        file_stop = len(files_to_use) if (args.frange).lower() == 'all' else int(args.frange)
    
    if file_start >= 0 and file_stop <= len(files_to_use):
        files_to_use = files_to_use[file_start:file_stop]
    else:
        raise IOError("The number of root files available for the %s sample is %i. args.frange must be less than or equal to this." % (sample, len(files_to_use) ) )

        ## replace sample name with group name ( [WZ]Jets -> EWK for instance)
    group_name = plt_tools.get_group(sample)
    if group_name in fileset.keys():
        for fname in files_to_use:
            fileset[group_name].append(fname)
    else:
        fileset[group_name] = files_to_use

#leptonSFs = lepSF.LeptonSF()
##set_trace()
#threejets_btagSFs = btagSF.create_btag_sf_computer('3')
#fourPlusjets_btagSFs = btagSF.create_btag_sf_computer('4+')

# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class Test_Analyzer(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.mass_axis = hist.Bin("mass", "m [GeV]", 100, 0, 5)
        self.pt_axis = hist.Bin("pt", "p_{T} [GeV]", 200, 0, 1000)
        self.eta_axis = hist.Bin("eta", r"$\eta$", 200, -5, 5)
        self.phi_axis = hist.Bin("phi", r"$\phi$", 160, -4, 4)
        self.energy_axis = hist.Bin("energy", "E [GeV]", 200, 0, 1000)
        self.njets_axis = hist.Bin("njets", "n_{jets}", 15, 0, 15)
        self.discr_axis = hist.Bin("discr", "", 160, 0., 40.)

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
                ## make best perm hists
            bp_hists = self.make_best_perm_hists(tdir)
            histo_dict.update(bp_hists)
        histo_dict['cutflow'] = processor.defaultdict_accumulator(int)

        self._accumulator = processor.dict_accumulator(histo_dict)
        self.sample_name = ''
        #self.plt_weights = plt_weights

        #    # get lepton SF info
        #self.leptonSFs = leptonSFs

    
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

    def make_best_perm_hists(self, tdir):
        histo_dict = {}

        histo_dict['%s_BestPerm_njets' % tdir] = hist.Hist("Counts", self.dataset_axis, self.njets_axis)

        histo_dict['%s_BestPerm_Prob' % tdir] = hist.Hist("Counts", self.dataset_axis, self.discr_axis)
        histo_dict['%s_BestPerm_MassDiscr' % tdir] = hist.Hist("Counts", self.dataset_axis, self.discr_axis)
        histo_dict['%s_BestPerm_NuDiscr' % tdir] = hist.Hist("Counts", self.dataset_axis, self.discr_axis)

        histo_dict['%s_BestPerm_BLep_pt' % tdir]    = hist.Hist("Counts", self.dataset_axis, self.pt_axis)
        histo_dict['%s_BestPerm_BLep_eta' % tdir]   = hist.Hist("Counts", self.dataset_axis, self.eta_axis)
        histo_dict['%s_BestPerm_BLep_phi' % tdir]   = hist.Hist("Counts", self.dataset_axis, self.phi_axis)
        histo_dict['%s_BestPerm_BLep_E' % tdir]   = hist.Hist("Counts", self.dataset_axis, self.energy_axis)

        histo_dict['%s_BestPerm_BHad_pt' % tdir]    = hist.Hist("Counts", self.dataset_axis, self.pt_axis)
        histo_dict['%s_BestPerm_BHad_eta' % tdir]   = hist.Hist("Counts", self.dataset_axis, self.eta_axis)
        histo_dict['%s_BestPerm_BHad_phi' % tdir]   = hist.Hist("Counts", self.dataset_axis, self.phi_axis)
        histo_dict['%s_BestPerm_BHad_E' % tdir]   = hist.Hist("Counts", self.dataset_axis, self.energy_axis)

        histo_dict['%s_BestPerm_WJa_pt' % tdir]    = hist.Hist("Counts", self.dataset_axis, self.pt_axis)
        histo_dict['%s_BestPerm_WJa_eta' % tdir]   = hist.Hist("Counts", self.dataset_axis, self.eta_axis)
        histo_dict['%s_BestPerm_WJa_phi' % tdir]   = hist.Hist("Counts", self.dataset_axis, self.phi_axis)
        histo_dict['%s_BestPerm_WJa_E' % tdir]   = hist.Hist("Counts", self.dataset_axis, self.energy_axis)

        if '4PJets' in tdir:
            histo_dict['%s_BestPerm_WJb_pt' % tdir]    = hist.Hist("Counts", self.dataset_axis, self.pt_axis)
            histo_dict['%s_BestPerm_WJb_eta' % tdir]   = hist.Hist("Counts", self.dataset_axis, self.eta_axis)
            histo_dict['%s_BestPerm_WJb_phi' % tdir]   = hist.Hist("Counts", self.dataset_axis, self.phi_axis)
            histo_dict['%s_BestPerm_WJb_E' % tdir]   = hist.Hist("Counts", self.dataset_axis, self.energy_axis)

        return histo_dict


    def process(self, df):
        output = self.accumulator.identity()

        if not isinstance(df, coffea.processor.dataframe.LazyDataFrame):
            raise IOError("This function only works for LazyDataFrame objects")
        #set_trace()
        self.sample_name = df.dataset
        isData = self.sample_name.startswith('data')

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

        evt_weights = MCWeights.evt_weight(df, mask=passing_evts)
        #evt_weights *= self.plt_weights[self.sample_name] # include plotting weight to event weights

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

        #    ## apply lepton SFs to MC (only applicable to tight leptons)
        #if not isData:
        #    lep_weights = leptonSFs.get_sf_(lepton='%ss' % lep_to_use, pt_array=sel_leps.pt.flatten(), eta_array=sel_leps.eta.flatten())
        #    lep_weights[~tight_leps] = 1.
        #    evt_weights *= lep_weights

        #    ## apply btagging SFs to MC
        #if not isData:
        #    btag_weights = np.ones(clean_jets.size)
        #    #set_trace()
        #        ## get per-jet weights for all systematic variations + central value
        #    threeJ_wts = threejets_btagSFs.get_scale_factor(jets=clean_jets[three_jets_events], passing_cut=btag_wps[0])
        #    fourPJ_wts = fourPlusjets_btagSFs.get_scale_factor(jets=clean_jets[fourPlus_jets_events], passing_cut=btag_wps[0])
        #        ## calculate per-event SFs for central value
        #    btag_weights[three_jets_events] = threeJ_wts['central'].prod()
        #    btag_weights[fourPlus_jets_events] = fourPJ_wts['central'].prod()
        #    evt_weights *= btag_weights


            ## find best permutations and create bp column
        #best_perms = Permutations.find_best_permutations(jets=clean_jets, leptons=sel_leps, MET=sel_met, evt_weights=evt_weights)
        #set_trace()

            ## fill hists for tight leptons
                ## 3 jets
        output = self.fill_jet_hists(output, '3Jets_TIGHT%s_Jets' % self.lepton[lep_to_use], clean_jets[(tight_leps & three_jets_events)], evt_weights[(tight_leps & three_jets_events)])
        output = self.fill_lep_hists(output, '3Jets_TIGHT%s_%s' % (self.lepton[lep_to_use], lep_to_use), sel_leps[(tight_leps & three_jets_events)], evt_weights[(tight_leps & three_jets_events)])
        #output = self.fill_best_perm_hists(output, '3Jets_TIGHT%s' % self.lepton[lep_to_use], best_perms[(best_perms.Leptons['TIGHT%s' % self.lepton[lep_to_use]].flatten()) & (best_perms.njets == 3)])

                ## 4+ jets
        output = self.fill_jet_hists(output, '4PJets_TIGHT%s_Jets' % self.lepton[lep_to_use], clean_jets[(tight_leps & fourPlus_jets_events)], evt_weights[(tight_leps & fourPlus_jets_events)])
        output = self.fill_lep_hists(output, '4PJets_TIGHT%s_%s' % (self.lepton[lep_to_use], lep_to_use), sel_leps[(tight_leps & fourPlus_jets_events)], evt_weights[(tight_leps & fourPlus_jets_events)])
        #output = self.fill_best_perm_hists(output, '4PJets_TIGHT%s' % self.lepton[lep_to_use], best_perms[(best_perms.Leptons['TIGHT%s' % self.lepton[lep_to_use]].flatten()) & (best_perms.njets > 3)])

            ## fill hists for loose leptons
                ## 3 jets
        output = self.fill_jet_hists(output, '3Jets_LOOSE%s_Jets' % self.lepton[lep_to_use], clean_jets[(loose_leps & three_jets_events)], evt_weights[(loose_leps & three_jets_events)])
        output = self.fill_lep_hists(output, '3Jets_LOOSE%s_%s' % (self.lepton[lep_to_use], lep_to_use), sel_leps[(loose_leps & three_jets_events)], evt_weights[(loose_leps & three_jets_events)])
        #output = self.fill_best_perm_hists(output, '3Jets_LOOSE%s' % self.lepton[lep_to_use], best_perms[(best_perms.Leptons['LOOSE%s' % self.lepton[lep_to_use]].flatten()) & (best_perms.njets == 3)])
                ## 4+ jets
        output = self.fill_jet_hists(output, '4PJets_LOOSE%s_Jets' % self.lepton[lep_to_use], clean_jets[(loose_leps & fourPlus_jets_events)], evt_weights[(loose_leps & fourPlus_jets_events)])
        output = self.fill_lep_hists(output, '4PJets_LOOSE%s_%s' % (self.lepton[lep_to_use], lep_to_use), sel_leps[(loose_leps & fourPlus_jets_events)], evt_weights[(loose_leps & fourPlus_jets_events)])
        #output = self.fill_best_perm_hists(output, '4PJets_LOOSE%s' % self.lepton[lep_to_use], best_perms[(best_perms.Leptons['LOOSE%s' % self.lepton[lep_to_use]].flatten()) & (best_perms.njets > 3)])

        return output

    def fill_jet_hists(self, accumulator, tdir, obj, evt_weights):
        #set_trace()
        accumulator['%s_pt' % tdir].fill(dataset=self.sample_name, pt=obj.pt.flatten(), weight=np.repeat(evt_weights, obj.counts))
        accumulator['%s_eta' % tdir].fill(dataset=self.sample_name, eta=obj.eta.flatten(), weight=np.repeat(evt_weights, obj.counts))
        accumulator['%s_phi' % tdir].fill(dataset=self.sample_name, phi=obj.phi.flatten(), weight=np.repeat(evt_weights, obj.counts))
        accumulator['%s_njets' % tdir].fill(dataset=self.sample_name, njets=obj.counts, weight=evt_weights)

        return accumulator        

    def fill_lep_hists(self, accumulator, tdir, obj, evt_weights):
        accumulator['%s_pt' % tdir].fill(dataset=self.sample_name, pt=obj.pt.flatten(), weight=evt_weights)
        accumulator['%s_eta' % tdir].fill(dataset=self.sample_name, eta=obj.eta.flatten(), weight=evt_weights)
        accumulator['%s_phi' % tdir].fill(dataset=self.sample_name, phi=obj.phi.flatten(), weight=evt_weights)

        return accumulator        

    def fill_best_perm_hists(self, accumulator, tdir, table):
        #set_trace()
        accumulator['%s_BestPerm_njets' % tdir].fill(dataset=self.sample_name, njets=table.njets, weight=table.evt_wts)

        accumulator['%s_BestPerm_Prob' % tdir].fill(dataset=self.sample_name, discr=table.Prob, weight=table.evt_wts)
        accumulator['%s_BestPerm_MassDiscr' % tdir].fill(dataset=self.sample_name, discr=table.MassDiscr, weight=table.evt_wts)
        accumulator['%s_BestPerm_NuDiscr' % tdir].fill(dataset=self.sample_name, discr=table.NuDiscr, weight=table.evt_wts)

        accumulator['%s_BestPerm_BLep_pt' % tdir].fill(dataset=self.sample_name,  pt  = table.BLeps.pt.flatten(), weight=table.evt_wts)
        accumulator['%s_BestPerm_BLep_eta' % tdir].fill(dataset=self.sample_name, eta = table.BLeps.eta.flatten(), weight=table.evt_wts)
        accumulator['%s_BestPerm_BLep_phi' % tdir].fill(dataset=self.sample_name, phi = table.BLeps.phi.flatten(), weight=table.evt_wts)
        accumulator['%s_BestPerm_BLep_E' % tdir].fill(dataset=self.sample_name,energy = table.BLeps.E.flatten(), weight=table.evt_wts)

        accumulator['%s_BestPerm_BHad_pt' % tdir].fill(dataset=self.sample_name,  pt  = table.BHads.pt.flatten(), weight=table.evt_wts)
        accumulator['%s_BestPerm_BHad_eta' % tdir].fill(dataset=self.sample_name, eta = table.BHads.eta.flatten(), weight=table.evt_wts)
        accumulator['%s_BestPerm_BHad_phi' % tdir].fill(dataset=self.sample_name, phi = table.BHads.phi.flatten(), weight=table.evt_wts)
        accumulator['%s_BestPerm_BHad_E' % tdir].fill(dataset=self.sample_name,energy = table.BHads.E.flatten(), weight=table.evt_wts)

        accumulator['%s_BestPerm_WJa_pt' % tdir].fill(dataset=self.sample_name,  pt  = table.WJas.pt.flatten(), weight=table.evt_wts)
        accumulator['%s_BestPerm_WJa_eta' % tdir].fill(dataset=self.sample_name, eta = table.WJas.eta.flatten(), weight=table.evt_wts)
        accumulator['%s_BestPerm_WJa_phi' % tdir].fill(dataset=self.sample_name, phi = table.WJas.phi.flatten(), weight=table.evt_wts)
        accumulator['%s_BestPerm_WJa_E' % tdir].fill(dataset=self.sample_name,energy = table.WJas.E.flatten(), weight=table.evt_wts)

        if '4PJets' in tdir:
            accumulator['%s_BestPerm_WJb_pt' % tdir].fill(dataset=self.sample_name,  pt  = table.WJbs.pt.flatten(), weight=table.evt_wts)
            accumulator['%s_BestPerm_WJb_eta' % tdir].fill(dataset=self.sample_name, eta = table.WJbs.eta.flatten(), weight=table.evt_wts)
            accumulator['%s_BestPerm_WJb_phi' % tdir].fill(dataset=self.sample_name, phi = table.WJbs.phi.flatten(), weight=table.evt_wts)
            accumulator['%s_BestPerm_WJb_E' % tdir].fill(dataset=self.sample_name,energy = table.WJbs.E.flatten(), weight=table.evt_wts)

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
    executor_args={
        'workers': 8,
        'flatten' : True,
        'compression': 5,
    },
    chunksize=10000,
    #chunksize=500000,
)

#output['data_lumi'] = data_lumi

#if args.debug:
#    print(output)
#set_trace()
print(output['cutflow'])

    ## save output to coffea pkl file
if (args.frange).lower() == 'all':
    outdir = '/'.join([proj_dir, 'results', '%s_%s' % (args.year, jobid), analyzer])
    cfname = '%s/%s.coffea' % (outdir, 'test')
    #cfname = '%s/%s.coffea' % (outdir, args.sample)

else:
    if ':' in args.frange:
        outdir = '/'.join([proj_dir, 'results', '%s_%s' % (args.year, jobid), analyzer])
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

