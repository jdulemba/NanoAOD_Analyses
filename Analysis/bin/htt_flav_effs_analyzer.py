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
import fnmatch

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'both_leps_analyzer'

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

allowed_samples = ['ttJets*', 'WJets', 'ZJets', 'singlet*']
for fname in fileset.keys():
    if not any([fnmatch.fnmatch(fname, sample) for sample in allowed_samples]):
        raise IOError("Sample %s not valid for finding btag efficiencies. Only %s are allowed" % (fname, allowed_samples))

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

jet_pars = prettyjson.loads(open('%s/cfg_files/cfg_pars_%s.json' % (proj_dir, jobid)).read())['Jets']
wps_to_use = list(set([jet_pars['permutations']['tightb'], jet_pars['permutations']['looseb']]))
if not( len(wps_to_use) == 1):
    raise IOError("Only 1 unique btag working point supported now")


# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class Both_Leps_Analyzer(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.btagger_axis = hist.Cat("btagger", "B-Tag WP")
        self.jetmult_axis = hist.Cat("jmult", "nJets")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")
        self.hflav_axis = hist.Cat("hFlav", "Hadron Flavour")
        self.pt_axis = hist.Bin("pt", "p_{T} [GeV]", 200, 0, 1000)
        self.eta_axis = hist.Bin("eta", r"$\eta$", 200, -3, 3)

            ## make dictionary of hists
        histo_dict = {}
                ## make jet hists
        jet_hists = self.make_jet_hists()
        histo_dict.update(jet_hists)

        histo_dict['cutflow'] = processor.defaultdict_accumulator(int)

        self._accumulator = processor.dict_accumulator(histo_dict)
        self.sample_name = ''
        self.corrections = corrections

    
    @property
    def accumulator(self):
        return self._accumulator


    def make_jet_hists(self):
        histo_dict = {}
        histo_dict['Jets_pt_all']     = hist.Hist("Events", self.btagger_axis, self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.hflav_axis, self.pt_axis)
        histo_dict['Jets_eta_all']    = hist.Hist("Events", self.btagger_axis, self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.hflav_axis, self.eta_axis)
        histo_dict['Jets_pt_eta_all'] = hist.Hist("Events", self.btagger_axis, self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.hflav_axis, self.pt_axis, self.eta_axis)
        histo_dict['Jets_pt_pass']    = hist.Hist("Events", self.btagger_axis, self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.hflav_axis, self.pt_axis)
        histo_dict['Jets_eta_pass']   = hist.Hist("Events", self.btagger_axis, self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.hflav_axis, self.eta_axis)
        histo_dict['Jets_pt_eta_pass']= hist.Hist("Events", self.btagger_axis, self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.hflav_axis, self.pt_axis, self.eta_axis)

        return histo_dict



    def process(self, df):
        output = self.accumulator.identity()

        if not isinstance(df, coffea.processor.dataframe.LazyDataFrame):
            raise IOError("This function only works for LazyDataFrame objects")

        #if args.debug: set_trace()
        self.sample_name = df.dataset
        isTTbar = self.sample_name.startswith('ttJets')

            ## make event weights
                # data or MC distinction made internally
        evt_weights = MCWeights.get_event_weights(df, year=args.year, corrections=self.corrections)

            ## initialize selections and regions
        selection = processor.PackedSelection()
        regions = {
            'Muon' : {
                '3Jets'  : {
                    'bjet' : {'objselection_mu', 'mu_jets_3', 'loose_or_tight_mu', 'mu_jets_bjet'},
                    'cjet' : {'objselection_mu', 'mu_jets_3', 'loose_or_tight_mu', 'mu_jets_cjet'},
                    'ljet' : {'objselection_mu', 'mu_jets_3', 'loose_or_tight_mu', 'mu_jets_ljet'},
                },
            },
            'Electron' : {
                '3Jets'  : {
                    'bjet' : {'objselection_el', 'el_jets_3', 'loose_or_tight_el', 'el_jets_bjet'},
                    'cjet' : {'objselection_el', 'el_jets_3', 'loose_or_tight_el', 'el_jets_cjet'},
                    'ljet' : {'objselection_el', 'el_jets_3', 'loose_or_tight_el', 'el_jets_ljet'},
                },
            },
        }

            ## object selection
        objsel_evts_mu = objsel.select(df, leptype='Muon', year=args.year, accumulator=output)
        output['cutflow']['nEvts passing jet and muon obj selection'] += objsel_evts_mu.sum()
        objsel_evts_el = objsel.select(df, leptype='Electron', year=args.year, accumulator=output)
        output['cutflow']['nEvts passing jet and electron obj selection'] += objsel_evts_el.sum()

            ## add different selections
                ## muons
        selection.add('objselection_mu', objsel_evts_mu)
        selection.add('mu_jets_3', df['Jet_Muon'].counts == 3)
            ## hadronFlavour definitions found here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools
        selection.add('mu_jets_bjet', (df['Jet_Muon']['hadronFlav'] == 5).sum() > 0) # jets come from bjets (hadronFlav == 5)
        selection.add('mu_jets_cjet', (df['Jet_Muon']['hadronFlav'] == 4).sum() > 0) # jets come from bjets (hadronFlav == 4)
        selection.add('mu_jets_ljet', (df['Jet_Muon']['hadronFlav'] == 0).sum() > 0) # jets come from light jets (hadronFlav == 0)
        selection.add('tight_mu', df['Muon']['TIGHTMU'].sum() == 1) # one muon passing TIGHT criteria
        selection.add('loose_or_tight_mu', (df['Muon']['LOOSEMU'] | df['Muon']['TIGHTMU']).sum() == 1) # one muon passing LOOSE or TIGHT criteria

                ## electrons
        selection.add('objselection_el', objsel_evts_el)
        selection.add('el_jets_3', df['Jet_Electron'].counts == 3)
        selection.add('el_jets_bjet', (df['Jet_Electron']['hadronFlav'] == 5).sum() > 0) # jets come from bjets (hadronFlav == 5)
        selection.add('el_jets_cjet', (df['Jet_Electron']['hadronFlav'] == 4).sum() > 0) # jets come from bjets (hadronFlav == 4)
        selection.add('el_jets_ljet', (df['Jet_Electron']['hadronFlav'] == 0).sum() > 0) # jets come from light jets (hadronFlav == 0)
        selection.add('tight_el', df['Electron']['TIGHTEL'].sum() == 1) # one electron passing TIGHT criteria
        selection.add('loose_or_tight_el', (df['Electron']['LOOSEEL'] | df['Electron']['TIGHTEL']).sum() == 1) # one electron passing LOOSE or TIGHT criteria

        ## apply lepton SFs to MC (only applicable to tight leptons)
        if 'LeptonSF' in corrections.keys():
            tight_mu_cut = selection.require(objselection_mu=True, tight_mu=True) # find events passing muon object selection with one tight muon
            evt_weights._weights['Muon_SF'][tight_mu_cut] = MCWeights.get_lepton_sf(year=args.year, lepton='Muons', corrections=lepSF_correction,
                pt=df['Muon'][tight_mu_cut].pt.flatten(), eta=df['Muon'][tight_mu_cut].eta.flatten())
            tight_el_cut = selection.require(objselection_el=True, tight_el=True) # find events passing electron object selection with one tight electron
            evt_weights._weights['Electron_SF'][tight_el_cut] = MCWeights.get_lepton_sf(year=args.year, lepton='Electrons', corrections=lepSF_correction,
                pt=df['Electron'][tight_el_cut].pt.flatten(), eta=df['Electron'][tight_el_cut].etaSC.flatten())

        if isTTbar:
            ## add 4+ jets categories for ttbar events
            regions['Muon'].update({
                '4PJets' : {
                    'bjet' : {'objselection_mu', 'mu_jets_4+', 'loose_or_tight_mu', 'mu_jets_bjet'},
                    'cjet' : {'objselection_mu', 'mu_jets_4+', 'loose_or_tight_mu', 'mu_jets_cjet'},
                    'ljet' : {'objselection_mu', 'mu_jets_4+', 'loose_or_tight_mu', 'mu_jets_ljet'},
                }
            })
            regions['Electron'].update({
                '4PJets' : {
                    'bjet' : {'objselection_el', 'el_jets_4+', 'loose_or_tight_el', 'el_jets_bjet'},
                    'cjet' : {'objselection_el', 'el_jets_4+', 'loose_or_tight_el', 'el_jets_cjet'},
                    'ljet' : {'objselection_el', 'el_jets_4+', 'loose_or_tight_el', 'el_jets_ljet'},
                }
            })
            selection.add('mu_jets_4+', df['Jet_Muon'].counts > 3)
            selection.add('el_jets_4+', df['Jet_Electron'].counts > 3)


        btag_wps = [wp for wp in df['Jet'].columns if wps_to_use[0] in wp]
        #set_trace()
        ## fill hists for each region
        for btag_wp in btag_wps:
            for lepton in regions.keys():
                for jmult in regions[lepton].keys():
                    for hflav in regions[lepton][jmult].keys():
                        cut = selection.all(*regions[lepton][jmult][hflav])
                        #set_trace()
    
                        evt_weights_to_use = evt_weights.weight()
                        ## apply lepton SFs to MC (only applicable to tight leptons)
                        if 'LeptonSF' in corrections.keys():
                            evt_weights_to_use = evt_weights.partial_weight(exclude=['Electron_SF']) if lepton == 'Muon' else evt_weights.partial_weight(exclude=['Muon_SF']) # exclude SF from other lepton
                        output = self.fill_jet_hists_all( accumulator=output, btag_wp=btag_wp, jetmult=jmult, leptype=lepton, hadFlav=hflav, obj=df['Jet_%s' % lepton][cut], evt_weights=evt_weights_to_use[cut])
                        output = self.fill_jet_hists_pass(accumulator=output, btag_wp=btag_wp, jetmult=jmult, leptype=lepton, hadFlav=hflav, obj=df['Jet_%s' % lepton][cut][(df['Jet_%s' % lepton][cut][btag_wp])], evt_weights=evt_weights_to_use[cut])

        return output


    def fill_jet_hists_all(self, accumulator, btag_wp, jetmult, leptype, hadFlav, obj, evt_weights):
        #set_trace()
        accumulator['Jets_pt_all'].fill(    btagger=btag_wp, dataset=self.sample_name, jmult=jetmult, leptype=leptype, hFlav=hadFlav, pt=obj.pt.flatten(), weight=np.repeat(evt_weights, obj.counts))
        accumulator['Jets_eta_all'].fill(   btagger=btag_wp, dataset=self.sample_name, jmult=jetmult, leptype=leptype, hFlav=hadFlav, eta=obj.eta.flatten(), weight=np.repeat(evt_weights, obj.counts))
        accumulator['Jets_pt_eta_all'].fill(btagger=btag_wp, dataset=self.sample_name, jmult=jetmult, leptype=leptype, hFlav=hadFlav, pt=obj.pt.flatten(), eta=obj.eta.flatten(), weight=np.repeat(evt_weights, obj.counts))
        return accumulator        

    def fill_jet_hists_pass(self, accumulator, btag_wp, jetmult, leptype, hadFlav, obj, evt_weights):
        #set_trace()
        accumulator['Jets_pt_pass'].fill(    btagger=btag_wp, dataset=self.sample_name, jmult=jetmult, leptype=leptype, hFlav=hadFlav, pt=obj.pt.flatten(), weight=np.repeat(evt_weights, obj.counts))
        accumulator['Jets_eta_pass'].fill(   btagger=btag_wp, dataset=self.sample_name, jmult=jetmult, leptype=leptype, hFlav=hadFlav, eta=obj.eta.flatten(), weight=np.repeat(evt_weights, obj.counts))
        accumulator['Jets_pt_eta_pass'].fill(btagger=btag_wp, dataset=self.sample_name, jmult=jetmult, leptype=leptype, hFlav=hadFlav, pt=obj.pt.flatten(), eta=obj.eta.flatten(), weight=np.repeat(evt_weights, obj.counts))
        return accumulator        


    def postprocess(self, accumulator):
        return accumulator

proc_executor = processor.iterative_executor if args.debug else processor.futures_executor

output = processor.run_uproot_job(fileset,
    treename='Events',
    processor_instance=Both_Leps_Analyzer(),
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
