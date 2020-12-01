#!/usr/bin/env python

from coffea import hist
from coffea.util import save, load
import coffea.processor as processor
from pdb import set_trace
import os, sys
import python.ObjectSelection as objsel
import python.MCWeights as MCWeights
import numpy as np
import Utilities.prettyjson as prettyjson
import coffea.lumi_tools.lumi_tools as lumi_tools
import fnmatch

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'htt_flav_effs'

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('fset', type=str, help='Fileset dictionary (in string form) to be used for the processor')
parser.add_argument('year', choices=['2016APV', '2016', '2017', '2018'], help='Specify which year to run over')
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
pu_correction = load(os.path.join(proj_dir, 'Corrections', jobid, 'MC_PU_Weights.coffea'))
lepSF_correction = load(os.path.join(proj_dir, 'Corrections', jobid, 'leptonSFs.coffea'))
jet_corrections = load(os.path.join(proj_dir, 'Corrections', jobid, 'JetCorrections.coffea'))[args.year]
corrections = {
    'Pileup' : pu_correction,
    'Prefire' : True,
    'LeptonSF' : lepSF_correction,
    'JetCor' : jet_corrections,
    'BTagSF' : False,
}

jet_pars = prettyjson.loads(open(os.path.join(proj_dir, 'cfg_files', 'cfg_pars_%s.json' % jobid)).read())['Jets']
wps_to_use = list(set([jet_pars['permutations']['tightb'], jet_pars['permutations']['looseb']]))
if not( len(wps_to_use) == 1):
    raise IOError("Only 1 unique btag working point supported now")


# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class Htt_Flav_Effs(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.btagger_axis = hist.Cat("btagger", "B-Tag WP")
        self.jetmult_axis = hist.Cat("jmult", "nJets")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")
        self.hflav_axis = hist.Cat("hFlav", "Hadron Flavour")
        self.pt_axis = hist.Bin("pt", "p_{T} [GeV]", 200, 0, 1000)
        self.eta_axis = hist.Bin("eta", r"$\eta$", 120, -3, 3)

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
        #histo_dict['Jets_pt_all']     = hist.Hist("Events", self.btagger_axis, self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.hflav_axis, self.pt_axis)
        #histo_dict['Jets_eta_all']    = hist.Hist("Events", self.btagger_axis, self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.hflav_axis, self.eta_axis)
        histo_dict['Jets_pt_eta_all'] = hist.Hist("Events", self.btagger_axis, self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.hflav_axis, self.pt_axis, self.eta_axis)
        #histo_dict['Jets_pt_pass']    = hist.Hist("Events", self.btagger_axis, self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.hflav_axis, self.pt_axis)
        #histo_dict['Jets_eta_pass']   = hist.Hist("Events", self.btagger_axis, self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.hflav_axis, self.eta_axis)
        histo_dict['Jets_pt_eta_pass']= hist.Hist("Events", self.btagger_axis, self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.hflav_axis, self.pt_axis, self.eta_axis)

        return histo_dict



    def process(self, df):
        output = self.accumulator.identity()

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
                    'objselection', 'jets_3', 'loose_or_tight_MU'
                },
            },
            'Electron' : {
                '3Jets'  : {
                    'objselection', 'jets_3', 'loose_or_tight_EL'
                },
            },
        }

            ## object selection
        objsel_evts = objsel.select(df, year=args.year, corrections=self.corrections, cutflow=output['cutflow'])
        output['cutflow']['nEvts passing jet and muon obj selection'] += objsel_evts.sum()
        selection.add('objselection', objsel_evts)

        ## add different selections
        selection.add('jets_3', df['Jet'].counts == 3)

                ## muons
        selection.add('tight_MU', df['Muon']['TIGHTMU'].sum() == 1) # one muon passing TIGHT criteria
        selection.add('loose_or_tight_MU', (df['Muon']['LOOSEMU'] | df['Muon']['TIGHTMU']).sum() == 1) # one muon passing LOOSE or TIGHT criteria

                ## electrons
        selection.add('tight_EL', df['Electron']['TIGHTEL'].sum() == 1) # one electron passing TIGHT criteria
        selection.add('loose_or_tight_EL', (df['Electron']['LOOSEEL'] | df['Electron']['TIGHTEL']).sum() == 1) # one electron passing LOOSE or TIGHT criteria

        ## apply lepton SFs to MC (only applicable to tight leptons)
        if 'LeptonSF' in corrections.keys():
            tight_mu_cut = selection.require(objselection=True, tight_MU=True) # find events passing muon object selection with one tight muon
            tight_muons = df['Muon'][tight_mu_cut][(df['Muon'][tight_mu_cut]['TIGHTMU'] == True)]
            muSFs_dict =  MCWeights.get_lepton_sf(year=args.year, lepton='Muons', corrections=lepSF_correction,
                pt=tight_muons.pt.flatten(), eta=tight_muons.eta.flatten())
            mu_reco_cen = np.ones(df.size)
            mu_reco_err = np.zeros(df.size)
            mu_trig_cen = np.ones(df.size)
            mu_trig_err = np.zeros(df.size)
            mu_reco_cen[tight_mu_cut] = muSFs_dict['RECO_CEN']
            mu_reco_err[tight_mu_cut] = muSFs_dict['RECO_ERR']
            mu_trig_cen[tight_mu_cut] = muSFs_dict['TRIG_CEN']
            mu_trig_err[tight_mu_cut] = muSFs_dict['TRIG_ERR']
            evt_weights.add('Muon_RECO', mu_reco_cen, mu_reco_err, mu_reco_err, shift=True)
            evt_weights.add('Muon_TRIG', mu_trig_cen, mu_trig_err, mu_trig_err, shift=True)

            tight_el_cut = selection.require(objselection=True, tight_EL=True) # find events passing electron object selection with one tight electron
            tight_electrons = df['Electron'][tight_el_cut][(df['Electron'][tight_el_cut]['TIGHTEL'] == True)]
            elSFs_dict = MCWeights.get_lepton_sf(year=args.year, lepton='Electrons', corrections=lepSF_correction,
                pt=tight_electrons.pt.flatten(), eta=tight_electrons.etaSC.flatten())
            el_reco_cen = np.ones(df.size)
            el_reco_err = np.zeros(df.size)
            el_trig_cen = np.ones(df.size)
            el_trig_err = np.zeros(df.size)
            el_reco_cen[tight_el_cut] = elSFs_dict['RECO_CEN']
            el_reco_err[tight_el_cut] = elSFs_dict['RECO_ERR']
            el_trig_cen[tight_el_cut] = elSFs_dict['TRIG_CEN']
            el_trig_err[tight_el_cut] = elSFs_dict['TRIG_ERR']
            evt_weights.add('Electron_RECO', el_reco_cen, el_reco_err, el_reco_err, shift=True)
            evt_weights.add('Electron_TRIG', el_trig_cen, el_trig_err, el_trig_err, shift=True)

        if isTTbar:
            ## add 4+ jets categories for ttbar events
            selection.add('jets_4+', df['Jet'].counts > 3)
            regions['Muon'].update({
                '4PJets' : {'objselection', 'jets_4+', 'loose_or_tight_MU'}
            })
            regions['Electron'].update({
                '4PJets' : {'objselection', 'jets_4+', 'loose_or_tight_EL'}
            })


        btag_wps = [wp for wp in df['Jet'].columns if wps_to_use[0] in wp]
        #set_trace()
        ## fill hists for each region
        for btag_wp in btag_wps:
            for lepton in regions.keys():
                for jmult in regions[lepton].keys():
                    cut = selection.all(*regions[lepton][jmult])
 
                    evt_weights_to_use = evt_weights.weight()
                    ## apply lepton SFs to MC (only applicable to tight leptons)
                    if 'LeptonSF' in corrections.keys():
                        evt_weights_to_use = evt_weights.partial_weight(exclude=['Electron_SF']) if lepton == 'Muon' else evt_weights.partial_weight(exclude=['Muon_SF']) # exclude SF from other lepton
                    #set_trace()
                    jets = df['Jet'][cut]
                        ## hadronFlavour definitions found here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools
                    bjets = jets[(jets['hadronFlav'] == 5)]
                    cjets = jets[(jets['hadronFlav'] == 4)]
                    ljets = jets[(jets['hadronFlav'] == 0)]

                    output = self.fill_jet_hists_all( acc=output, btag_wp=btag_wp, jetmult=jmult, leptype=lepton, hadFlav='bjet', obj=bjets, evt_weights=evt_weights_to_use[cut])
                    output = self.fill_jet_hists_pass(acc=output, btag_wp=btag_wp, jetmult=jmult, leptype=lepton, hadFlav='bjet', obj=bjets[(bjets[btag_wp])], evt_weights=evt_weights_to_use[cut])
                    output = self.fill_jet_hists_all( acc=output, btag_wp=btag_wp, jetmult=jmult, leptype=lepton, hadFlav='cjet', obj=cjets, evt_weights=evt_weights_to_use[cut])
                    output = self.fill_jet_hists_pass(acc=output, btag_wp=btag_wp, jetmult=jmult, leptype=lepton, hadFlav='cjet', obj=cjets[(cjets[btag_wp])], evt_weights=evt_weights_to_use[cut])
                    output = self.fill_jet_hists_all( acc=output, btag_wp=btag_wp, jetmult=jmult, leptype=lepton, hadFlav='ljet', obj=ljets, evt_weights=evt_weights_to_use[cut])
                    output = self.fill_jet_hists_pass(acc=output, btag_wp=btag_wp, jetmult=jmult, leptype=lepton, hadFlav='ljet', obj=ljets[(ljets[btag_wp])], evt_weights=evt_weights_to_use[cut])

        return output


    def fill_jet_hists_all(self, acc, btag_wp, jetmult, leptype, hadFlav, obj, evt_weights):
        #set_trace()
        #acc['Jets_pt_all'].fill(    btagger=btag_wp, dataset=self.sample_name, jmult=jetmult, leptype=leptype, hFlav=hadFlav, pt=obj.pt.flatten(), weight=(obj.pt.ones_like()*evt_weights).flatten())
        #acc['Jets_eta_all'].fill(   btagger=btag_wp, dataset=self.sample_name, jmult=jetmult, leptype=leptype, hFlav=hadFlav, eta=obj.eta.flatten(), weight=(obj.pt.ones_like()*evt_weights).flatten())
        acc['Jets_pt_eta_all'].fill(btagger=btag_wp, dataset=self.sample_name, jmult=jetmult, leptype=leptype, hFlav=hadFlav, pt=obj.pt.flatten(), eta=obj.eta.flatten(), weight=(obj.pt.ones_like()*evt_weights).flatten())
        return acc        

    def fill_jet_hists_pass(self, acc, btag_wp, jetmult, leptype, hadFlav, obj, evt_weights):
        #set_trace()
        #acc['Jets_pt_pass'].fill(    btagger=btag_wp, dataset=self.sample_name, jmult=jetmult, leptype=leptype, hFlav=hadFlav, pt=obj.pt.flatten(), weight=(obj.pt.ones_like()*evt_weights).flatten())
        #acc['Jets_eta_pass'].fill(   btagger=btag_wp, dataset=self.sample_name, jmult=jetmult, leptype=leptype, hFlav=hadFlav, eta=obj.eta.flatten(), weight=(obj.pt.ones_like()*evt_weights).flatten())
        acc['Jets_pt_eta_pass'].fill(btagger=btag_wp, dataset=self.sample_name, jmult=jetmult, leptype=leptype, hFlav=hadFlav, pt=obj.pt.flatten(), eta=obj.eta.flatten(), weight=(obj.pt.ones_like()*evt_weights).flatten())
        return acc        


    def postprocess(self, accumulator):
        return accumulator

proc_executor = processor.iterative_executor if args.debug else processor.futures_executor

output = processor.run_uproot_job(fileset,
    treename='Events',
    processor_instance=Htt_Flav_Effs(),
    executor=proc_executor,
    executor_args={
        'workers': 8,
        'flatten' : True,
        'compression': 5,
    },
    chunksize=10000 if args.debug else 100000,
    #chunksize=10000 if args.debug else 50000,
    #chunksize=10000,
)


if args.debug:
    print(output)
#set_trace()
#print(output['cutflow'])

save(output, args.outfname)
print('%s has been written' % args.outfname)
