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
import Utilities.make_variables as make_vars

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
base_jobid = jobid.split('_')[0]
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
pu_correction = load(os.path.join(proj_dir, 'Corrections', jobid, 'MC_PU_Weights.coffea'))
lepSF_correction = load(os.path.join(proj_dir, 'Corrections', jobid, 'leptonSFs.coffea'))
jet_corrections = load(os.path.join(proj_dir, 'Corrections', jobid, 'JetCorrections.coffea'))[args.year]
corrections = {
    #'Pileup' : pu_correction,
    'Prefire' : True,
    'LeptonSF' : lepSF_correction,
    'BTagSF' : False,
    'JetCor' : jet_corrections,
    #'JetCor' : None,
}


# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class presel_analyzer(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.jetmult_axis = hist.Cat("jmult", "nJets")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")
        self.pt_axis = hist.Bin("pt", "p_{T} [GeV]", 200, 0, 1000)
        self.eta_axis = hist.Bin("eta", r"$\eta$", 60, -3, 3)
        self.phi_axis = hist.Bin("phi", r"$\phi$", 160, -4, 4)
        self.energy_axis = hist.Bin("energy", "E [GeV]", 200, 0, 1000)
        self.njets_axis = hist.Bin("njets", "n_{jets}", 20, 0, 20)
        self.lepIso_axis = hist.Bin("iso", "pfRelIso", 100, 0., 1.)
        self.mt_axis = hist.Bin("mt", "M_{T}", 200, 0., 1000.)

            ## make dictionary of hists
        histo_dict = {}
                ## make hists
        hists = self.make_hists()
        histo_dict.update(hists)

        histo_dict['cutflow'] = processor.defaultdict_accumulator(int)

        self._accumulator = processor.dict_accumulator(histo_dict)
        self.sample_name = ''
        self.corrections = corrections
        self.isData = True

    
    @property
    def accumulator(self):
        return self._accumulator


    def make_hists(self):
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

        histo_dict['MT']        = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.mt_axis)

        histo_dict['Lep_pt']    = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.pt_axis)
        histo_dict['Lep_eta']   = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.eta_axis)
        histo_dict['Lep_etaSC'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.eta_axis)
        histo_dict['Lep_phi']   = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.phi_axis)
        histo_dict['Lep_energy']= hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.energy_axis)
        histo_dict['Lep_iso']   = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.lepIso_axis)

        return histo_dict



    def process(self, df):
        np.random.seed(10) # sets seed so values from random distributions are reproducible (JER corrections)
        output = self.accumulator.identity()

        self.sample_name = df.dataset

            ## make event weights
                # data or MC distinction made internally
                    # data or MC distinction made internally
        mu_evt_weights = MCWeights.get_event_weights(df, year=args.year, corrections=self.corrections)
        el_evt_weights = MCWeights.get_event_weights(df, year=args.year, corrections=self.corrections)

            ## initialize selections and regions
        selection = processor.PackedSelection()
        regions = {
            'Muon' : {
                '3Jets'  : {'objselection', 'jets_3', 'loose_or_tight_MU'},
                '4PJets' : {'objselection', 'jets_4+', 'loose_or_tight_MU'},
            },
            'Electron' : {
                '3Jets'  : {'objselection', 'jets_3', 'loose_or_tight_EL'},
                '4PJets' : {'objselection', 'jets_4+', 'loose_or_tight_EL'},
            },
        }

            ## object selection
        objsel_evts = objsel.select(df, year=args.year, corrections=self.corrections, cutflow=output['cutflow'])
        output['cutflow']['nEvts passing jet and lepton obj selection'] += objsel_evts.sum()
        selection.add('jets_3', df['Jet'].counts == 3)
        selection.add('jets_4+', df['Jet'].counts > 3)
        selection.add('objselection', objsel_evts)

        self.isData = self.sample_name.startswith('data_Single')
        if self.isData:
            isSE_Data = self.sample_name.startswith('data_SingleElectron')
            isSM_Data = self.sample_name.startswith('data_SingleMuon')
            runs = df.run
            lumis = df.luminosityBlock
            Golden_Json_LumiMask = lumi_tools.LumiMask(os.path.join(proj_dir, 'inputs', 'data', base_jobid, 'LumiMasks', '%s_GoldenJson_%s.txt' % (args.year, base_jobid)))
            LumiMask = Golden_Json_LumiMask.__call__(runs, lumis) ## returns array of valid events
            selection.add('lumimask', LumiMask)
   
                ## object selection and add different selections
            if isSM_Data:
                del regions['Electron']
                        ## muons
                selection.add('tight_MU', df['Muon']['TIGHTMU'].sum() == 1) # one muon passing TIGHT criteria
                #selection.add('loose_MU', df['Muon']['LOOSEMU'].sum() == 1) # one muon passing LOOSE criteria
                selection.add('loose_or_tight_MU', (df['Muon']['LOOSEMU'] | df['Muon']['TIGHTMU']).sum() == 1) # one muon passing LOOSE or TIGHT criteria
            if isSE_Data:
                del regions['Muon']
                        ## electrons
                selection.add('tight_EL', df['Electron']['TIGHTEL'].sum() == 1) # one electron passing TIGHT criteria
                #selection.add('loose_EL', df['Electron']['LOOSEEL'].sum() == 1) # one electron passing LOOSE criteria
                selection.add('loose_or_tight_EL', (df['Electron']['LOOSEEL'] | df['Electron']['TIGHTEL']).sum() == 1) # one electron passing LOOSE or TIGHT criteria

            for lepton in regions.keys():
                for jmult in regions[lepton].keys():
                    regions[lepton][jmult].update({'lumimask'})

        if not self.isData:
                ## add different selections
                    ## muons
            selection.add('tight_MU', df['Muon']['TIGHTMU'].sum() == 1) # one muon passing TIGHT criteria
            #selection.add('loose_MU', df['Muon']['LOOSEMU'].sum() == 1) # one muon passing LOOSE criteria
            selection.add('loose_or_tight_MU', (df['Muon']['LOOSEMU'] | df['Muon']['TIGHTMU']).sum() == 1) # one muon passing LOOSE or TIGHT criteria
                    ## electrons
            selection.add('tight_EL', df['Electron']['TIGHTEL'].sum() == 1) # one electron passing TIGHT criteria
            #selection.add('loose_EL', df['Electron']['LOOSEEL'].sum() == 1) # one electron passing LOOSE criteria
            selection.add('loose_or_tight_EL', (df['Electron']['LOOSEEL'] | df['Electron']['TIGHTEL']).sum() == 1) # one electron passing LOOSE or TIGHT criteria

            #set_trace()
            ### apply lepton SFs to MC (only applicable to tight leptons)
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
                mu_evt_weights.add('Lep_RECO', mu_reco_cen, mu_reco_err, mu_reco_err, shift=True)
                mu_evt_weights.add('Lep_TRIG', mu_trig_cen, mu_trig_err, mu_trig_err, shift=True)

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
                el_evt_weights.add('Lep_RECO', el_reco_cen, el_reco_err, el_reco_err, shift=True)
                el_evt_weights.add('Lep_TRIG', el_trig_cen, el_trig_err, el_trig_err, shift=True)



        #set_trace()
        ## fill hists for each region
        for lepton in regions.keys():
            evt_weights = mu_evt_weights if lepton == 'Muon' else el_evt_weights
            for jmult in regions[lepton].keys():
                cut = selection.all(*regions[lepton][jmult])
                #set_trace()

                if cut.sum() > 0:
                    leptype = 'MU' if lepton == 'Muon' else 'EL'
                    if 'loose_or_tight_%s' % leptype in regions[lepton][jmult]:
                        lep_mask = ((df[lepton][cut]['TIGHT%s' % leptype] == True) | (df[lepton][cut]['LOOSE%s' % leptype] == True))
                    elif 'tight_%s' % leptype in regions[lepton][jmult]:
                        lep_mask = (df[lepton][cut]['TIGHT%s' % leptype] == True)
                    elif 'loose_%s' % leptype in regions[lepton][jmult]:
                        lep_mask = (df[lepton][cut]['LOOSE%s' % leptype] == True)
                    else:
                        raise ValueError("Not sure what lepton type to choose for event")


                    evt_weights_to_use = evt_weights.weight()
                    if not self.isData:
                        ## apply lepton SFs to MC (only applicable to tight leptons)
                        if 'LeptonSF' in corrections.keys():
                            evt_weights_to_use = evt_weights.partial_weight(exclude=['Electron_SF']) if lepton == 'Muon' else evt_weights.partial_weight(exclude=['Muon_SF']) # exclude SF from other lepton
                    output = self.fill_jet_hists(acc=output, jetmult=jmult, leptype=lepton, obj=df['Jet'][cut], evt_weights=evt_weights_to_use[cut])
                    output = self.fill_lep_hists(acc=output, jetmult=jmult, leptype=lepton, obj=df[lepton][cut][lep_mask],evt_weights=evt_weights_to_use[cut])

                            ## get MT
                    MT = make_vars.MT(df[lepton][cut][lep_mask], df['MET'][cut])
                    output['MT'].fill(dataset=self.sample_name, jmult=jmult, leptype=lepton, mt=MT.flatten(), weight=(MT.ones_like()*evt_weights_to_use[cut]).flatten())

        return output


    def fill_jet_hists(self, acc, jetmult, leptype, obj, evt_weights):
        #set_trace()
        acc['Jets_pt'].fill(    dataset=self.sample_name, jmult=jetmult, leptype=leptype, pt=obj.pt.flatten(), weight=(obj.pt.ones_like()*evt_weights).flatten())
        acc['Jets_eta'].fill(   dataset=self.sample_name, jmult=jetmult, leptype=leptype, eta=obj.eta.flatten(), weight=(obj.pt.ones_like()*evt_weights).flatten())
        acc['Jets_phi'].fill(   dataset=self.sample_name, jmult=jetmult, leptype=leptype, phi=obj.phi.flatten(), weight=(obj.pt.ones_like()*evt_weights).flatten())
        acc['Jets_energy'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, energy=obj.p4.E.flatten(), weight=(obj.pt.ones_like()*evt_weights).flatten())
        acc['Jets_njets'].fill( dataset=self.sample_name, jmult=jetmult, leptype=leptype, njets=obj.counts, weight=evt_weights)
        acc['Jets_LeadJet_pt'].fill(    dataset=self.sample_name, jmult=jetmult, leptype=leptype, pt=obj.pt.max(), weight=evt_weights)
        acc['Jets_LeadJet_eta'].fill(   dataset=self.sample_name, jmult=jetmult, leptype=leptype, eta=obj.eta[:, 0], weight=evt_weights)
        acc['Jets_LeadJet_phi'].fill(   dataset=self.sample_name, jmult=jetmult, leptype=leptype, phi=obj.phi[:, 0], weight=evt_weights)
        acc['Jets_LeadJet_energy'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, energy=obj.p4.E[:, 0], weight=evt_weights)

        return acc        

    def fill_lep_hists(self, acc, jetmult, leptype, obj, evt_weights):
        #set_trace()
        acc['Lep_pt'].fill(    dataset=self.sample_name, jmult=jetmult, leptype=leptype, pt=obj.pt.flatten(), weight=evt_weights)
        acc['Lep_eta'].fill(   dataset=self.sample_name, jmult=jetmult, leptype=leptype, eta=obj.eta.flatten(), weight=evt_weights)
        if leptype == 'Electron': acc['Lep_etaSC'].fill(   dataset=self.sample_name, jmult=jetmult, leptype=leptype, eta=obj.etaSC.flatten(), weight=evt_weights)
        acc['Lep_phi'].fill(   dataset=self.sample_name, jmult=jetmult, leptype=leptype, phi=obj.phi.flatten(), weight=evt_weights)
        acc['Lep_energy'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, energy=obj.p4.E.flatten(), weight=evt_weights)
        acc['Lep_iso'].fill(   dataset=self.sample_name, jmult=jetmult, leptype=leptype, iso=obj.pfRelIso.flatten(), weight=evt_weights)

        return acc        


    def postprocess(self, accumulator):
        return accumulator



proc_executor = processor.iterative_executor if args.debug else processor.futures_executor
output = processor.run_uproot_job(fileset,
    treename='Events',
    processor_instance=presel_analyzer(),
    executor=proc_executor,
    executor_args={
        'workers': 8,
        'flatten' : True,
        'compression': 5,
    },
    chunksize=10000 if args.debug else 100000,
    #chunksize=10000 if args.debug else 50000,
    #chunksize=50000,
)


if args.debug:
    print(output)

save(output, args.outfname)
print('%s has been written' % args.outfname)
