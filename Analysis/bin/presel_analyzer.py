#!/usr/bin/env python

from coffea import hist
from coffea.util import save, load
import coffea.processor as processor
from pdb import set_trace
import os, sys
import python.ObjectSelection as objsel
import coffea.processor.dataframe
import Utilities.plot_tools as plt_tools
import python.LeptonSF as lepSF
import python.BTagScaleFactors as btagSF
import python.MCWeights as MCWeights
import numpy as np
import Utilities.prettyjson as prettyjson
import coffea.lumi_tools.lumi_tools as lumi_tools

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'presel_analyzer'

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('frange', type=str, help='Specify start:stop indices for files')
parser.add_argument('year', choices=['2016', '2017', '2018'], help='Specify which year to run over')
parser.add_argument('lepton', choices=['Electron', 'Muon'], help='Choose which lepton to select')
parser.add_argument('--sample', type=str, help='Use specific sample')
parser.add_argument('--outfname', type=str, help='Specify output filename, including directory and file extension')
parser.add_argument('--debug', action='store_true', help='Uses iterative_executor for debugging purposes, otherwise futures_excutor will be used (faster)')

args = parser.parse_args()

    ## get samples to use
indir = '/'.join([proj_dir, 'inputs', '%s_%s' % (args.year, jobid)])
if args.sample:
        ## sample specified
    if not os.path.isfile('%s/%s.txt' % (indir, args.sample)):
        raise IOError("File with samples %s.txt not found" % args.sample)

    samples = [args.sample]
    if args.outfname:
        print("  --- Sample name %s will be overridden by output fname %s ---  \n" % (args.sample, args.outfname))

else:
        ## get file that has names for all datasets to use
    fpath = '/'.join([proj_dir, 'inputs', '%s_%s' % (args.year, jobid), '%s_inputs.txt' % analyzer])
    if not os.path.isfile(fpath):
        raise IOError("File with samples %s_inputs.txt not found" % analyzer)
    
    txt_file = open(fpath, 'r')
    samples = [sample.strip('\n') for sample in txt_file if not sample.startswith('#')]
    if not samples:
        raise IOError("No samples found as inputs")

    ## add files to fileset and get plotting weights
fileset = {}
for sample in samples:
    txtpath = '/'.join([proj_dir, 'inputs', '%s_%s' % (args.year, jobid), '%s.txt' % sample])
    if not os.path.isfile(txtpath):
        raise IOError("Sample file %s.txt not found" % sample)

    txtfiles = open(txtpath, 'r')
    files_to_use = [fname.strip('\n') for fname in txtfiles]
    file_inds = [idx for idx, val in enumerate(files_to_use)]

    if ':' in args.frange:
        file_start, file_stop = int((args.frange).split(':')[0]), int((args.frange).split(':')[1])
    else:
        file_start = 0
        file_stop = len(files_to_use) if (args.frange).lower() == 'all' else int(args.frange)
    
    if file_start >= 0 and file_stop < len(files_to_use):
        files_to_use = files_to_use[file_start:file_stop+1]
        inds_to_use = file_inds[file_start:file_stop+1]
    else:
        raise IOError("The number of root files available for the %s sample is %i. args.frange must be less than this." % (sample, len(files_to_use) ) )

    print(inds_to_use)
    #set_trace()
    fileset[sample] = files_to_use

print(fileset)
#sys.exit()
#set_trace()
## load corrections for event weights
pu_correction = load('%s/Corrections/MC_PU_Weights.coffea' % proj_dir)
lumi_correction = load('%s/Corrections/MC_LumiWeights.coffea' % proj_dir)
corrections = {
    'Pileup' : pu_correction,
    'Lumi' : lumi_correction,
    'Prefire' : True,
    'LeptonSF' : False,
    'BTagSF' : False,
}

if corrections['LeptonSF'] == True:
    leptonSFs = lepSF.LeptonSF()
if corrections['BTagSF'] == True:
    threejets_btagSFs = btagSF.create_btag_sf_computer('3')
    fourPlusjets_btagSFs = btagSF.create_btag_sf_computer('4+')


# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class Presel_Analyzer(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.jetmult_axis = hist.Cat("jmult", "nJets")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")
        self.mass_axis = hist.Bin("mass", "m [GeV]", 100, 0, 5)
        self.pt_axis = hist.Bin("pt", "p_{T} [GeV]", 200, 0, 1000)
        self.eta_axis = hist.Bin("eta", r"$\eta$", 200, -5, 5)
        self.phi_axis = hist.Bin("phi", r"$\phi$", 160, -4, 4)
        self.energy_axis = hist.Bin("energy", "E [GeV]", 200, 0, 1000)
        self.njets_axis = hist.Bin("njets", "n_{jets}", 15, 0, 15)
        self.lepIso_axis = hist.Bin("iso", "pfRelIso", 100, 0., 1.)

        #self.lepton = {
        #    'Muon': 'MU',
        #    #'Electron' : 'EL'
        #}
        self.leptypes = ['LOOSEMU', 'TIGHTMU'] if args.lepton == 'Muon' else ['LOOSEEL', 'TIGHTEL']

        #set_trace()        
            ## make dictionary of hists
        histo_dict = {}
                ## make jet hists
        jet_hists = self.make_jet_hists()
        histo_dict.update(jet_hists)
                ## make lepton hists
        lep_hists = self.make_lep_hists()
        histo_dict.update(lep_hists)        

        histo_dict['cutflow'] = processor.defaultdict_accumulator(int)

        self._accumulator = processor.dict_accumulator(histo_dict)
        self.sample_name = ''
        self.corrections = corrections

        #    # get lepton SF info
        #self.leptonSFs = leptonSFs

    
    @property
    def accumulator(self):
        return self._accumulator


    def make_jet_hists(self):
        histo_dict = {}
        histo_dict['Jets_pt']    = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.pt_axis)
        histo_dict['Jets_eta']   = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.eta_axis)
        histo_dict['Jets_phi']   = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.phi_axis)
        histo_dict['Jets_energy']= hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.energy_axis)
        histo_dict['Jets_njets'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.njets_axis)

        return histo_dict

    
    def make_lep_hists(self):
        histo_dict = {}
        histo_dict['Lep_pt']    = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.pt_axis)
        histo_dict['Lep_eta']   = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.eta_axis)
        histo_dict['Lep_phi']   = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.phi_axis)
        histo_dict['Lep_energy']= hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.energy_axis)
        histo_dict['Lep_iso']   = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepIso_axis)

        return histo_dict



    def process(self, df):
        output = self.accumulator.identity()

        if not isinstance(df, coffea.processor.dataframe.LazyDataFrame):
            raise IOError("This function only works for LazyDataFrame objects")

        #if args.debug: set_trace()
        self.sample_name = df.dataset
        lep_to_use = args.lepton
        #lep_to_use = [*self.lepton.keys()][0]

            ## initialize selections and regions
        selection = processor.PackedSelection()
        regions = {}
        #regions['objsel'] = {'objselection'}
        regions['3Jets'] = {'objselection', '3jets'}
        regions['4PJets'] = {'objselection', '4pjets'}

        isData = self.sample_name.startswith('data_Single')
        if isData:
            runs = df.run
            lumis = df.luminosityBlock
            Golden_Json_LumiMask = lumi_tools.LumiMask('%s/inputs/data/LumiMasks/%s_GoldenJson.txt' % (proj_dir, args.year))
            LumiMask = Golden_Json_LumiMask.__call__(runs, lumis) ## returns array of valid events
            #set_trace()
            selection.add('lumimask', LumiMask)
            for region in regions.keys():
                regions[region].update({'lumimask'})
   
 
            ## make event weights
                # data or MC distinction made internally
        evt_weights = MCWeights.get_event_weights(df, year=args.year, lepton='%ss' % lep_to_use, corrections=self.corrections)

        objsel_evts = objsel.select(df, leptype=lep_to_use, accumulator=output)
        output['cutflow']['nEvts passing jet and %s objection' % lep_to_use] += objsel_evts.sum()
        selection.add('objselection', objsel_evts)

        
        selection.add('3jets', df['Jet'].counts == 3)
        selection.add('4pjets', df['Jet'].counts > 3)

        #set_trace()
        ## fill hists for each region
        for region in regions.keys():
            #set_trace()
            cut = selection.all(*regions[region])
            output = self.fill_jet_hists(output, region, df['Jet'][cut], evt_weights.weight()[cut])
            output = self.fill_lep_hists(output, region, df[lep_to_use][cut], evt_weights.weight()[cut])


        ##    ## apply lepton SFs to MC (only applicable to tight leptons)
        ##if not isData:
        ##    lep_weights = leptonSFs.get_sf_(lepton='%ss' % lep_to_use, pt_array=sel_leps.pt.flatten(), eta_array=sel_leps.eta.flatten())
        ##    lep_weights[~tight_leps] = 1.
        ##    evt_weights *= lep_weights

        ##    ## apply btagging SFs to MC
        ##if not isData:
        ##    btag_weights = np.ones(clean_jets.size)
        ##    #set_trace()
        ##        ## get per-jet weights for all systematic variations + central value
        ##    threeJ_wts = threejets_btagSFs.get_scale_factor(jets=clean_jets[three_jets_events], passing_cut=btag_wps[0])
        ##    fourPJ_wts = fourPlusjets_btagSFs.get_scale_factor(jets=clean_jets[fourPlus_jets_events], passing_cut=btag_wps[0])
        ##        ## calculate per-event SFs for central value
        ##    btag_weights[three_jets_events] = threeJ_wts['central'].prod()
        ##    btag_weights[fourPlus_jets_events] = fourPJ_wts['central'].prod()
        ##    evt_weights *= btag_weights



        return output

    def fill_jet_hists(self, accumulator, jetmult, obj, evt_weights):
        #set_trace()
        accumulator['Jets_pt'].fill(    dataset=self.sample_name, jmult=jetmult, leptype=args.lepton, pt=obj.pt.flatten(), weight=np.repeat(evt_weights, obj.counts))
        accumulator['Jets_eta'].fill(   dataset=self.sample_name, jmult=jetmult, leptype=args.lepton, eta=obj.eta.flatten(), weight=np.repeat(evt_weights, obj.counts))
        accumulator['Jets_phi'].fill(   dataset=self.sample_name, jmult=jetmult, leptype=args.lepton, phi=obj.phi.flatten(), weight=np.repeat(evt_weights, obj.counts))
        accumulator['Jets_energy'].fill(dataset=self.sample_name, jmult=jetmult, leptype=args.lepton, energy=obj.p4.E.flatten(), weight=np.repeat(evt_weights, obj.counts))
        accumulator['Jets_njets'].fill( dataset=self.sample_name, jmult=jetmult, leptype=args.lepton, njets=obj.counts, weight=evt_weights)

        return accumulator        

    def fill_lep_hists(self, accumulator, jetmult, obj, evt_weights):
        #set_trace()
        accumulator['Lep_pt'].fill(    dataset=self.sample_name, jmult=jetmult, leptype=args.lepton, pt=obj.pt.flatten(), weight=evt_weights)
        accumulator['Lep_eta'].fill(   dataset=self.sample_name, jmult=jetmult, leptype=args.lepton, eta=obj.eta.flatten(), weight=evt_weights)
        accumulator['Lep_phi'].fill(   dataset=self.sample_name, jmult=jetmult, leptype=args.lepton, phi=obj.phi.flatten(), weight=evt_weights)
        accumulator['Lep_energy'].fill(dataset=self.sample_name, jmult=jetmult, leptype=args.lepton, energy=obj.p4.E.flatten(), weight=evt_weights)
        accumulator['Lep_iso'].fill(   dataset=self.sample_name, jmult=jetmult, leptype=args.lepton, iso=obj.pfRelIso.flatten(), weight=evt_weights)

        return accumulator        


    def postprocess(self, accumulator):
        return accumulator

proc_executor = processor.iterative_executor if args.debug else processor.futures_executor

output = processor.run_uproot_job(fileset,
    treename='Events',
    processor_instance=Presel_Analyzer(),
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

    ## save output to coffea pkl file
if (args.frange).lower() == 'all':
    outdir = '/'.join([proj_dir, 'results', '%s_%s' % (args.year, jobid), analyzer, args.lepton])
    if args.outfname:
        cfname = args.outfname
    elif args.sample:
        cfname = '%s/%s.coffea' % (outdir, args.sample)
    else:
        cfname = '%s/test_%s.coffea' % (outdir, analyzer)
else:
    if ':' in args.frange:
        outdir = '/'.join([proj_dir, 'results', '%s_%s' % (args.year, jobid), analyzer, args.lepton])
        if args.outfname:
            cfname = args.outfname
        elif args.sample:
            cfname = '%s/%s_%sto%s.coffea' % (outdir, args.sample, file_start, file_stop)
        else:
            cfname = '%s/test_%sto%s.coffea' % (outdir, file_start, file_stop)
    else:
        outdir = proj_dir
        if args.outfname:
            cfname = args.outfname
        elif args.sample:
            cfname = '%s/%s_%s.test.coffea' % (outdir, args.sample, analyzer)
        else:
            cfname = '%s/%s.test.coffea' % (outdir, analyzer)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

save(output, cfname)
print('%s has been written' % cfname)

