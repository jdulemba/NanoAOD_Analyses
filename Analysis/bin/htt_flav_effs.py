#!/usr/bin/env python

from coffea import hist, processor
from coffea.nanoevents import NanoAODSchema
from pdb import set_trace
import os
from coffea.util import save, load
import Utilities.prettyjson as prettyjson
import numpy as np
import python.ObjectSelection as objsel
import python.MCWeights as MCWeights
import fnmatch
from coffea.analysis_tools import PackedSelection
import awkward as ak
import python.IDJet as IDJet
import python.GenParticleSelector as genpsel

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'htt_flav_effs'
base_jobid = os.environ['base_jobid']

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('fset', type=str, help='Fileset dictionary (in string form) to be used for the processor')
parser.add_argument('year', choices=['2016APV', '2016', '2017', '2018'] if base_jobid == 'ULnanoAOD' else ['2016', '2017', '2018'], help='Specify which year to run over')
parser.add_argument('outfname', type=str, help='Specify output filename, including directory and file extension')
parser.add_argument('opts', type=str, help='Fileset dictionary (in string form) to be used for the processor')

args = parser.parse_args()

# convert input string of fileset dictionary to actual dictionary
fdict = (args.fset).replace("\'", "\"")
fileset = prettyjson.loads(fdict)

# convert input string of options dictionary to actual dictionary
odict = (args.opts).replace("\'", "\"")
opts_dict = prettyjson.loads(odict)

    ## set config options passed through argparse
import ast
to_debug = ast.literal_eval(opts_dict.get('debug', 'False'))

allowed_samples = ['ttJets*', 'WJets', 'ZJets', 'singlet*']
for fname in fileset.keys():
    if not any([fnmatch.fnmatch(fname, sample) for sample in allowed_samples]):
        raise IOError("Sample %s not valid for finding btag efficiencies. Only %s are allowed" % (fname, allowed_samples))

## load corrections for event weights
pu_correction = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'MC_PU_Weights.coffea'))[args.year]
lepSF_correction = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'leptonSFs.coffea'))[args.year]
jet_corrections = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'JetMETCorrections.coffea'))[args.year]
nnlo_var = 'mtt_vs_thad_ctstar_Interp'
nnlo_reweighting = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'NNLO_to_Tune_Orig_Interp_Ratios_%s.coffea' % base_jobid))[args.year]
corrections = {
    'Pileup' : pu_correction,
    'Prefire' : True,
    'LeptonSF' : lepSF_correction,
    'JetCor' : jet_corrections,
    'BTagSF' : False,
    'NNLO_Rewt' : {'Var' : nnlo_var, 'Correction' : nnlo_reweighting[nnlo_var]},
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



    def process(self, events):
        output = self.accumulator.identity()

        self.sample_name = events.metadata['dataset']
        isTTbar = self.sample_name.startswith('ttJets')

            ## make event weights
                # data or MC distinction made internally
        mu_evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections)
        el_evt_weights = MCWeights.get_event_weights(events, year=args.year, corrections=self.corrections)

            ## initialize selections and regions
        selection = PackedSelection()
        regions = {
            'Muon' : {
                '3Jets'  : {
                    'lep_and_filter_pass', 'passing_jets', 'jets_3', 'loose_or_tight_MU'
                },
            },
            'Electron' : {
                '3Jets'  : {
                    'lep_and_filter_pass', 'passing_jets', 'jets_3', 'loose_or_tight_EL'
                },
            },
        }

            # get all passing leptons
        lep_and_filter_pass = objsel.select_leptons(events, year=args.year, cutflow=output['cutflow'])
        selection.add('lep_and_filter_pass', lep_and_filter_pass)
                ## muons
        selection.add('tight_MU', ak.sum(events['Muon']['TIGHTMU'], axis=1) == 1) # one muon passing TIGHT criteria
        selection.add('loose_or_tight_MU', ak.sum(events['Muon']['LOOSEMU'] | events['Muon']['TIGHTMU'], axis=1) == 1) # one muon passing LOOSE or TIGHT criteria
                ## electrons
        selection.add('tight_EL', ak.sum(events['Electron']['TIGHTEL'], axis=1) == 1) # one muon passing TIGHT criteria
        selection.add('loose_or_tight_EL', ak.sum(events['Electron']['LOOSEEL'] | events['Electron']['TIGHTEL'], axis=1) == 1) # one muon passing LOOSE or TIGHT criteria

            ## build corrected jets and MET
        events['Jet'], events['MET'] = IDJet.process_jets(events, args.year, self.corrections['JetCor'])

            # jet selection
        passing_jets = objsel.jets_selection(events, year=args.year, cutflow=output['cutflow'])
        selection.add('passing_jets', passing_jets)
        output['cutflow']['nEvts passing jet and lepton obj selection'] += ak.sum(passing_jets & lep_and_filter_pass)
        selection.add('jets_3', ak.num(events['SelectedJets']) == 3)

        ## apply lepton SFs to MC (only applicable to tight leptons)
        if 'LeptonSF' in self.corrections.keys():
            tight_mu_cut = selection.require(tight_MU=True) # find events passing muon object selection with one tight muon
            tight_muons = events['Muon'][tight_mu_cut][(events['Muon'][tight_mu_cut]['TIGHTMU'] == True)]
            muSFs_dict =  MCWeights.get_lepton_sf(year=args.year, lepton='Muons', corrections=self.corrections['LeptonSF'],
                pt=ak.flatten(tight_muons['pt']), eta=ak.flatten(tight_muons['eta']))
            mu_reco_cen = np.ones(len(events))
            mu_reco_err = np.zeros(len(events))
            mu_trig_cen = np.ones(len(events))
            mu_trig_err = np.zeros(len(events))
            mu_reco_cen[tight_mu_cut] = muSFs_dict['RECO_CEN']
            mu_reco_err[tight_mu_cut] = muSFs_dict['RECO_ERR']
            mu_trig_cen[tight_mu_cut] = muSFs_dict['TRIG_CEN']
            mu_trig_err[tight_mu_cut] = muSFs_dict['TRIG_ERR']
            mu_evt_weights.add('RECO', mu_reco_cen, mu_reco_err, mu_reco_err, shift=True)
            mu_evt_weights.add('TRIG', mu_trig_cen, mu_trig_err, mu_trig_err, shift=True)

            tight_el_cut = selection.require(tight_EL=True) # find events passing electron object selection with one tight electron
            tight_electrons = events['Electron'][tight_el_cut][(events['Electron'][tight_el_cut]['TIGHTEL'] == True)]
            elSFs_dict = MCWeights.get_lepton_sf(year=args.year, lepton='Electrons', corrections=self.corrections['LeptonSF'],
                pt=ak.flatten(tight_electrons['pt']), eta=ak.flatten(tight_electrons['etaSC']))
            el_reco_cen = np.ones(len(events))
            el_reco_err = np.zeros(len(events))
            el_trig_cen = np.ones(len(events))
            el_trig_err = np.zeros(len(events))
            el_reco_cen[tight_el_cut] = elSFs_dict['RECO_CEN']
            el_reco_err[tight_el_cut] = elSFs_dict['RECO_ERR']
            el_trig_cen[tight_el_cut] = elSFs_dict['TRIG_CEN']
            el_trig_err[tight_el_cut] = elSFs_dict['TRIG_ERR']
            el_evt_weights.add('RECO', el_reco_cen, el_reco_err, el_reco_err, shift=True)
            el_evt_weights.add('TRIG', el_trig_cen, el_trig_err, el_trig_err, shift=True)

        if isTTbar:
            ## add 4+ jets categories for ttbar events
            selection.add('jets_4+', ak.num(events['SelectedJets']) > 3)
            regions['Muon'].update({
                '4PJets' : {'lep_and_filter_pass', 'passing_jets', 'jets_4+', 'loose_or_tight_MU'}
            })
            regions['Electron'].update({
                '4PJets' : {'lep_and_filter_pass', 'passing_jets', 'jets_4+', 'loose_or_tight_EL'}
            })
            if 'NNLO_Rewt' in self.corrections.keys():
                    # find gen level particles for ttbar system
                genpsel.select(events, mode='NORMAL')
                nnlo_wts = MCWeights.get_nnlo_weights(self.corrections['NNLO_Rewt'], events)
                mu_evt_weights.add('%s_reweighting' % self.corrections['NNLO_Rewt']['Var'], nnlo_wts)
                el_evt_weights.add('%s_reweighting' % self.corrections['NNLO_Rewt']['Var'], nnlo_wts)


        btag_wps = [wp for wp in events['SelectedJets'].fields if wps_to_use[0] in wp]
        #set_trace()
        ## fill hists for each region
        for btag_wp in btag_wps:
            for lepton in regions.keys():
                evt_weights = mu_evt_weights if lepton == 'Muon' else el_evt_weights
                for jmult in regions[lepton].keys():
                    cut = selection.all(*regions[lepton][jmult])
 
                    evt_weights_to_use = evt_weights.weight()
                    jets = events['SelectedJets'][cut]
                        ## hadronFlavour definitions found here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideBTagMCTools
                    bjets = jets[(jets['hadronFlavour'] == 5)]
                    cjets = jets[(jets['hadronFlavour'] == 4)]
                    ljets = jets[(jets['hadronFlavour'] == 0)]

                    output = self.fill_jet_hists_all( acc=output, btag_wp=btag_wp, jetmult=jmult, leptype=lepton, hadFlav='bjet', obj=bjets, evt_weights=evt_weights_to_use[cut])
                    output = self.fill_jet_hists_pass(acc=output, btag_wp=btag_wp, jetmult=jmult, leptype=lepton, hadFlav='bjet', obj=bjets[(bjets[btag_wp])], evt_weights=evt_weights_to_use[cut])
                    output = self.fill_jet_hists_all( acc=output, btag_wp=btag_wp, jetmult=jmult, leptype=lepton, hadFlav='cjet', obj=cjets, evt_weights=evt_weights_to_use[cut])
                    output = self.fill_jet_hists_pass(acc=output, btag_wp=btag_wp, jetmult=jmult, leptype=lepton, hadFlav='cjet', obj=cjets[(cjets[btag_wp])], evt_weights=evt_weights_to_use[cut])
                    output = self.fill_jet_hists_all( acc=output, btag_wp=btag_wp, jetmult=jmult, leptype=lepton, hadFlav='ljet', obj=ljets, evt_weights=evt_weights_to_use[cut])
                    output = self.fill_jet_hists_pass(acc=output, btag_wp=btag_wp, jetmult=jmult, leptype=lepton, hadFlav='ljet', obj=ljets[(ljets[btag_wp])], evt_weights=evt_weights_to_use[cut])

        return output


    def fill_jet_hists_all(self, acc, btag_wp, jetmult, leptype, hadFlav, obj, evt_weights):
        #set_trace()
        #acc['Jets_pt_all'].fill(    btagger=btag_wp, dataset=self.sample_name, jmult=jetmult, leptype=leptype, hFlav=hadFlav, pt=ak.flatten(obj['pt']), weight=ak.flatten((ak.ones_like(obj['pt'])*evt_weights)))
        #acc['Jets_eta_all'].fill(   btagger=btag_wp, dataset=self.sample_name, jmult=jetmult, leptype=leptype, hFlav=hadFlav, eta=ak.flatten(obj['eta']), weight=ak.flatten((ak.ones_like(obj['pt'])*evt_weights)))
        acc['Jets_pt_eta_all'].fill(btagger=btag_wp, dataset=self.sample_name, jmult=jetmult, leptype=leptype, hFlav=hadFlav, pt=ak.flatten(obj['pt']), eta=ak.flatten(obj['eta']), weight=ak.flatten((ak.ones_like(obj['pt'])*evt_weights)))
        return acc        

    def fill_jet_hists_pass(self, acc, btag_wp, jetmult, leptype, hadFlav, obj, evt_weights):
        #set_trace()
        #acc['Jets_pt_pass'].fill(    btagger=btag_wp, dataset=self.sample_name, jmult=jetmult, leptype=leptype, hFlav=hadFlav, pt=ak.flatten(obj['pt']), weight=ak.flatten((ak.ones_like(obj['pt'])*evt_weights)))
        #acc['Jets_eta_pass'].fill(   btagger=btag_wp, dataset=self.sample_name, jmult=jetmult, leptype=leptype, hFlav=hadFlav, eta=ak.flatten(obj['eta']), weight=ak.flatten((ak.ones_like(obj['pt'])*evt_weights)))
        acc['Jets_pt_eta_pass'].fill(btagger=btag_wp, dataset=self.sample_name, jmult=jetmult, leptype=leptype, hFlav=hadFlav, pt=ak.flatten(obj['pt']), eta=ak.flatten(obj['eta']), weight=ak.flatten((ak.ones_like(obj['pt'])*evt_weights)))
        return acc        


    def postprocess(self, accumulator):
        return accumulator


proc_executor = processor.iterative_executor if to_debug else processor.futures_executor
proc_exec_args = {"schema": NanoAODSchema} if to_debug else {"schema": NanoAODSchema, "workers": 8}
output = processor.run_uproot_job(
    fileset,
    treename="Events",
    processor_instance=Htt_Flav_Effs(),
    executor=proc_executor,
    executor_args=proc_exec_args,
    chunksize=100000,
)

save(output, args.outfname)
print('%s has been written' % args.outfname)
