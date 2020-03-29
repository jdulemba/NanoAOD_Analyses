#!/usr/bin/env python

from coffea import hist
from coffea.util import save, load
import coffea.processor as processor
from pdb import set_trace
import os, sys
import python.ObjectSelection as objsel
import coffea.processor.dataframe
import Utilities.plot_tools as plt_tools
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
parser.add_argument('fset', type=str, help='Fileset dictionary (in string form) to be used for the processor')
parser.add_argument('year', choices=['2016', '2017', '2018'], help='Specify which year to run over')
parser.add_argument('outfname', type=str, help='Specify output filename, including directory and file extension')
parser.add_argument('--debug', action='store_true', help='Uses iterative_executor for debugging purposes, otherwise futures_excutor will be used (faster)')

args = parser.parse_args()

# convert input string of fileset dictionary to actual dictionary
fdict = (args.fset).replace("\'", "\"")
fileset = prettyjson.loads(fdict)

## load corrections for event weights
pu_correction = load('%s/Corrections/%s/MC_PU_Weights.coffea' % (proj_dir, jobid))
lepSF_correction = load('%s/Corrections/leptonSFs.coffea' % proj_dir)
corrections = {
    'Pileup' : pu_correction,
    'Prefire' : True,
    'LeptonSF' : lepSF_correction,
    'BTagSF' : False,
}

if corrections['BTagSF'] == True:
    threejets_btagSFs = btagSF.create_btag_sf_computer(args.year, '3')
    fourPlusjets_btagSFs = btagSF.create_btag_sf_computer(args.year, '4+')


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
        histo_dict['Jets_LeadJet_pt']    = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.pt_axis)
        histo_dict['Jets_LeadJet_eta']   = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.eta_axis)
        histo_dict['Jets_LeadJet_phi']   = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.phi_axis)
        histo_dict['Jets_LeadJet_energy']= hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.energy_axis)

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

            ## make event weights
                # data or MC distinction made internally
        evt_weights = MCWeights.get_event_weights(df, year=args.year, corrections=self.corrections)

            ## initialize selections and regions
        selection = processor.PackedSelection()
        regions = {
            'Muon' : {
                '3Jets'  : {'objselection_mu', 'mu_jets_3', 'loose_or_tight_mu'},
                '4PJets' : {'objselection_mu', 'mu_jets_4+', 'loose_or_tight_mu'},
            },
            'Electron' : {
                '3Jets'  : {'objselection_el', 'el_jets_3', 'loose_or_tight_el'},
                '4PJets' : {'objselection_el', 'el_jets_4+', 'loose_or_tight_el'},
            },
        }

        isData = self.sample_name.startswith('data_Single')
        if isData:
            isSE_Data = self.sample_name.startswith('data_SingleElectron')
            isSM_Data = self.sample_name.startswith('data_SingleMuon')
            runs = df.run
            lumis = df.luminosityBlock
            Golden_Json_LumiMask = lumi_tools.LumiMask('%s/inputs/data/LumiMasks/%s_GoldenJson.txt' % (proj_dir, args.year))
            LumiMask = Golden_Json_LumiMask.__call__(runs, lumis) ## returns array of valid events
            selection.add('lumimask', LumiMask)
   
                ## object selection and add different selections
            if isSM_Data:
                del regions['Electron']
                        ## muons
                objsel_evts_mu = objsel.select(df, leptype='Muon', year=args.year, accumulator=output)
                output['cutflow']['nEvts passing jet and muon obj selection'] += objsel_evts_mu.sum()
                selection.add('mu_jets_3', df['Jet_Muon'].counts == 3)
                selection.add('mu_jets_4+', df['Jet_Muon'].counts > 3)
                selection.add('objselection_mu', objsel_evts_mu)
                selection.add('tight_mu', df['Muon']['TIGHTMU'].sum() == 1) # one muon passing TIGHT criteria
                #selection.add('loose_mu', df['Muon']['LOOSEMU'].sum() == 1) # one muon passing LOOSE criteria
                selection.add('loose_or_tight_mu', (df['Muon']['LOOSEMU'] | df['Muon']['TIGHTMU']).sum() == 1) # one muon passing LOOSE or TIGHT criteria
            if isSE_Data:
                del regions['Muon']
                        ## electrons
                objsel_evts_el = objsel.select(df, leptype='Electron', year=args.year, accumulator=output)
                output['cutflow']['nEvts passing jet and electron obj selection'] += objsel_evts_el.sum()
                selection.add('el_jets_3', df['Jet_Electron'].counts == 3)
                selection.add('el_jets_4+', df['Jet_Electron'].counts > 3)
                selection.add('objselection_el', objsel_evts_el)
                selection.add('tight_el', df['Electron']['TIGHTEL'].sum() == 1) # one electron passing TIGHT criteria
                #selection.add('loose_el', df['Electron']['LOOSEEL'].sum() == 1) # one electron passing LOOSE criteria
                selection.add('loose_or_tight_el', (df['Electron']['LOOSEEL'] | df['Electron']['TIGHTEL']).sum() == 1) # one electron passing LOOSE or TIGHT criteria

            for lepton in regions.keys():
                for jmult in regions[lepton].keys():
                    regions[lepton][jmult].update({'lumimask'})

        if not isData:
                ## object selection
            objsel_evts_mu = objsel.select(df, leptype='Muon', year=args.year, accumulator=output)
            output['cutflow']['nEvts passing jet and muon obj selection'] += objsel_evts_mu.sum()
            objsel_evts_el = objsel.select(df, leptype='Electron', year=args.year, accumulator=output)
            output['cutflow']['nEvts passing jet and electron obj selection'] += objsel_evts_el.sum()

                ## add different selections
                    ## muons
            selection.add('mu_jets_3', df['Jet_Muon'].counts == 3)
            selection.add('mu_jets_4+', df['Jet_Muon'].counts > 3)
            selection.add('objselection_mu', objsel_evts_mu)
            selection.add('tight_mu', df['Muon']['TIGHTMU'].sum() == 1) # one muon passing TIGHT criteria
            #selection.add('loose_mu', df['Muon']['LOOSEMU'].sum() == 1) # one muon passing LOOSE criteria
            selection.add('loose_or_tight_mu', (df['Muon']['LOOSEMU'] | df['Muon']['TIGHTMU']).sum() == 1) # one muon passing LOOSE or TIGHT criteria
                    ## electrons
            selection.add('el_jets_3', df['Jet_Electron'].counts == 3)
            selection.add('el_jets_4+', df['Jet_Electron'].counts > 3)
            selection.add('objselection_el', objsel_evts_el)
            selection.add('tight_el', df['Electron']['TIGHTEL'].sum() == 1) # one electron passing TIGHT criteria
            #selection.add('loose_el', df['Electron']['LOOSEEL'].sum() == 1) # one electron passing LOOSE criteria
            selection.add('loose_or_tight_el', (df['Electron']['LOOSEEL'] | df['Electron']['TIGHTEL']).sum() == 1) # one electron passing LOOSE or TIGHT criteria

            ## apply lepton SFs to MC (only applicable to tight leptons)
            if 'LeptonSF' in corrections.keys():
                tight_mu_cut = selection.require(objselection_mu=True, tight_mu=True) # find events passing muon object selection with one tight muon
                evt_weights._weights['Muon_SF'][tight_mu_cut] = MCWeights.get_lepton_sf(year=args.year, lepton='Muons', corrections=lepSF_correction,
                    pt=df['Muon'][tight_mu_cut].pt.flatten(), eta=df['Muon'][tight_mu_cut].eta.flatten())
                tight_el_cut = selection.require(objselection_el=True, tight_el=True) # find events passing electron object selection with one tight electron
                evt_weights._weights['Electron_SF'][tight_el_cut] = MCWeights.get_lepton_sf(year=args.year, lepton='Electrons', corrections=lepSF_correction,
                    pt=df['Electron'][tight_el_cut].pt.flatten(), eta=df['Electron'][tight_el_cut].etaSC.flatten())

        #set_trace()
        ## fill hists for each region
        for lepton in regions.keys():
            for jmult in regions[lepton].keys():
                cut = selection.all(*regions[lepton][jmult])
                #set_trace()

                evt_weights_to_use = evt_weights.weight()
                if not isData:
                    ## apply lepton SFs to MC (only applicable to tight leptons)
                    if 'LeptonSF' in corrections.keys():
                        evt_weights_to_use = evt_weights.partial_weight(exclude=['Electron_SF']) if lepton == 'Muon' else evt_weights.partial_weight(exclude=['Muon_SF']) # exclude SF from other lepton
                output = self.fill_jet_hists(accumulator=output, jetmult=jmult, leptype=lepton, obj=df['Jet_%s' % lepton][cut], evt_weights=evt_weights_to_use[cut])
                output = self.fill_lep_hists(accumulator=output, jetmult=jmult, leptype=lepton, obj=df[lepton][cut],            evt_weights=evt_weights_to_use[cut])


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

    def fill_jet_hists(self, accumulator, jetmult, leptype, obj, evt_weights):
        #set_trace()
        accumulator['Jets_pt'].fill(    dataset=self.sample_name, jmult=jetmult, leptype=leptype, pt=obj.pt.flatten(), weight=np.repeat(evt_weights, obj.counts))
        accumulator['Jets_eta'].fill(   dataset=self.sample_name, jmult=jetmult, leptype=leptype, eta=obj.eta.flatten(), weight=np.repeat(evt_weights, obj.counts))
        accumulator['Jets_phi'].fill(   dataset=self.sample_name, jmult=jetmult, leptype=leptype, phi=obj.phi.flatten(), weight=np.repeat(evt_weights, obj.counts))
        accumulator['Jets_energy'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, energy=obj.p4.E.flatten(), weight=np.repeat(evt_weights, obj.counts))
        accumulator['Jets_njets'].fill( dataset=self.sample_name, jmult=jetmult, leptype=leptype, njets=obj.counts, weight=evt_weights)
        accumulator['Jets_LeadJet_pt'].fill(    dataset=self.sample_name, jmult=jetmult, leptype=leptype, pt=obj.pt.max(), weight=evt_weights)
        accumulator['Jets_LeadJet_eta'].fill(   dataset=self.sample_name, jmult=jetmult, leptype=leptype, eta=obj.eta[:, 0], weight=evt_weights)
        accumulator['Jets_LeadJet_phi'].fill(   dataset=self.sample_name, jmult=jetmult, leptype=leptype, phi=obj.phi[:, 0], weight=evt_weights)
        accumulator['Jets_LeadJet_energy'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, energy=obj.p4.E[:, 0], weight=evt_weights)

        return accumulator        

    def fill_lep_hists(self, accumulator, jetmult, leptype, obj, evt_weights):
        #set_trace()
        accumulator['Lep_pt'].fill(    dataset=self.sample_name, jmult=jetmult, leptype=leptype, pt=obj.pt.flatten(), weight=evt_weights)
        accumulator['Lep_eta'].fill(   dataset=self.sample_name, jmult=jetmult, leptype=leptype, eta=obj.eta.flatten(), weight=evt_weights)
        accumulator['Lep_phi'].fill(   dataset=self.sample_name, jmult=jetmult, leptype=leptype, phi=obj.phi.flatten(), weight=evt_weights)
        accumulator['Lep_energy'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, energy=obj.p4.E.flatten(), weight=evt_weights)
        accumulator['Lep_iso'].fill(   dataset=self.sample_name, jmult=jetmult, leptype=leptype, iso=obj.pfRelIso.flatten(), weight=evt_weights)

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

save(output, args.outfname)
print('%s has been written' % args.outfname)
