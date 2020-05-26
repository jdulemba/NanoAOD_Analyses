#!/usr/bin/env python

from coffea import hist
from coffea.util import save, load
import coffea.processor as processor
from pdb import set_trace
import os, sys
import python.ObjectSelection as objsel
#import coffea.processor.dataframe
import Utilities.plot_tools as plt_tools
import python.MCWeights as MCWeights
import numpy as np
import Utilities.prettyjson as prettyjson
import python.GenParticleSelector as genpsel
import python.TTGenMatcher as ttmatcher
import python.TTPermutator as ttpermutator
import Utilities.make_variables as make_vars

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'ttbar_alpha_reco'

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

    ## run on correct samples
if args.year == '2016':
    Nominal_ttJets = ['ttJets_PS', 'ttJets']
else:
    Nominal_ttJets = ['ttJetsSL', 'ttJetsHad', 'ttJetsDiLep']
isTTbar = np.array([(key in Nominal_ttJets) for key in fileset.keys()]).all()
if not isTTbar:
    raise ValueError("This should only be run on nominal ttbar events!")

init_btag = ~(np.array([key.startswith('data') for key in fileset.keys()]).all())

## load corrections for event weights
pu_correction = load('%s/Corrections/%s/MC_PU_Weights.coffea' % (proj_dir, jobid))
lepSF_correction = load('%s/Corrections/leptonSFs.coffea' % proj_dir)
jet_corrections = load('%s/Corrections/JetCorrections.coffea' % proj_dir)[args.year]
corrections = {
    'Pileup' : pu_correction,
    'Prefire' : True,
    'LeptonSF' : lepSF_correction,
    'BTagSF' : init_btag,
    'JetCor' : jet_corrections,
}

    ## parameters for b-tagging
jet_pars = prettyjson.loads(open('%s/cfg_files/cfg_pars_%s.json' % (proj_dir, jobid)).read())['Jets']
btagger = jet_pars['btagger']
wps_to_use = list(set([jet_pars['permutations']['tightb'],jet_pars['permutations']['looseb']]))
if not( len(wps_to_use) == 1):
    raise ValueError("Only 1 unique btag working point supported now")
btag_wp = btagger+wps_to_use[0]

if corrections['BTagSF'] == True:
    sf_file = '%s/Corrections/%s/%s' % (proj_dir, jobid, jet_pars['btagging']['btagSF_file'])
    if not os.path.isfile(sf_file):
        raise ValueError("BTag SF file %s doesn't exist" % sf_file)

    btag_sfs = load(sf_file)
    threeJets = btag_sfs[args.year][btagger]['3Jets'][wps_to_use[0]]
    fourPJets = btag_sfs[args.year][btagger]['4PJets'][wps_to_use[0]]
    corrections.update({'BTag_Constructors' : {'3Jets' : threeJets, '4PJets' : fourPJets} })


# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class ttbar_alpha_reco(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")
        self.eventtype_axis = hist.Cat("evttype", "Event Type (Lost or Merged)")
        self.objtype_axis = hist.Cat("objtype", "Object Type")
        self.pt_axis     = hist.Bin("pt", "p_{T} [GeV]", 200, 0, 1000)
        self.eta_axis    = hist.Bin("eta", r"$\eta$", 60, -3, 3)
        self.phi_axis    = hist.Bin("phi", r"$\phi$", 160, -4, 4)
        self.mass_axis   = hist.Bin("mass", "m [GeV]", 1000, 0., 2000.)
        self.energy_axis = hist.Bin("energy", "E [GeV]", 200, 0, 1000)
        self.reso_pt_axis     = hist.Bin("reso_pt", "p_{T} [GeV]", 200, -200, 200)
        self.reso_eta_axis    = hist.Bin("reso_eta", r"$\eta$", 100, -5, 5)
        self.reso_phi_axis    = hist.Bin("reso_phi", r"$\phi$", 100, -5, 5)
        self.reso_mass_axis   = hist.Bin("reso_mass", "m [GeV]", 1000, -200., 200.)
        self.reso_energy_axis = hist.Bin("reso_energy", "E [GeV]", 1000, -200, 200)

            ## make dictionary of hists
        histo_dict = {}
            ## make jet hists
                # Gen
        gen_hists = self.make_gen_hists()
        histo_dict.update(gen_hists)
                # Reco
        reco_hists = self.make_reco_hists()
        histo_dict.update(reco_hists)
                # Reso
        reso_hists = self.make_reso_hists()
        histo_dict.update(reso_hists)

        histo_dict['cutflow'] = processor.defaultdict_accumulator(int)

        self._accumulator = processor.dict_accumulator(histo_dict)
        self.sample_name = ''
        self.corrections = corrections

    
    @property
    def accumulator(self):
        return self._accumulator

    def make_gen_hists(self):
        histo_dict = {}
        histo_dict['Gen_pt']     = hist.Hist("Events", self.dataset_axis, self.leptype_axis, self.eventtype_axis, self.objtype_axis, self.pt_axis)
        histo_dict['Gen_eta']    = hist.Hist("Events", self.dataset_axis, self.leptype_axis, self.eventtype_axis, self.objtype_axis, self.eta_axis)
        histo_dict['Gen_phi']    = hist.Hist("Events", self.dataset_axis, self.leptype_axis, self.eventtype_axis, self.objtype_axis, self.phi_axis)
        histo_dict['Gen_mass']   = hist.Hist("Events", self.dataset_axis, self.leptype_axis, self.eventtype_axis, self.objtype_axis, self.mass_axis)
        histo_dict['Gen_energy'] = hist.Hist("Events", self.dataset_axis, self.leptype_axis, self.eventtype_axis, self.objtype_axis, self.energy_axis)

        return histo_dict

    def make_reco_hists(self):
        histo_dict = {}
        histo_dict['Reco_pt']     = hist.Hist("Events", self.dataset_axis, self.leptype_axis, self.eventtype_axis, self.objtype_axis, self.pt_axis)
        histo_dict['Reco_eta']    = hist.Hist("Events", self.dataset_axis, self.leptype_axis, self.eventtype_axis, self.objtype_axis, self.eta_axis)
        histo_dict['Reco_phi']    = hist.Hist("Events", self.dataset_axis, self.leptype_axis, self.eventtype_axis, self.objtype_axis, self.phi_axis)
        histo_dict['Reco_mass']   = hist.Hist("Events", self.dataset_axis, self.leptype_axis, self.eventtype_axis, self.objtype_axis, self.mass_axis)
        histo_dict['Reco_energy'] = hist.Hist("Events", self.dataset_axis, self.leptype_axis, self.eventtype_axis, self.objtype_axis, self.energy_axis)

        return histo_dict

    def make_reso_hists(self):
        histo_dict = {}
        histo_dict['Reso_pt']     = hist.Hist("Events", self.dataset_axis, self.leptype_axis, self.eventtype_axis, self.objtype_axis, self.reso_pt_axis)
        histo_dict['Reso_eta']    = hist.Hist("Events", self.dataset_axis, self.leptype_axis, self.eventtype_axis, self.objtype_axis, self.reso_eta_axis)
        histo_dict['Reso_phi']    = hist.Hist("Events", self.dataset_axis, self.leptype_axis, self.eventtype_axis, self.objtype_axis, self.reso_phi_axis)
        histo_dict['Reso_mass']   = hist.Hist("Events", self.dataset_axis, self.leptype_axis, self.eventtype_axis, self.objtype_axis, self.reso_mass_axis)
        histo_dict['Reso_energy'] = hist.Hist("Events", self.dataset_axis, self.leptype_axis, self.eventtype_axis, self.objtype_axis, self.reso_energy_axis)

        return histo_dict

    
    def process(self, df):
        np.random.seed(10) # sets seed so values from random distributions are reproducible (JER corrections)
        output = self.accumulator.identity()

        self.sample_name = df.dataset

            ## make event weights
                # data or MC distinction made internally
        evt_weights = MCWeights.get_event_weights(df, year=args.year, corrections=self.corrections)

            ## initialize selections and regions
        selection = processor.PackedSelection()
        regions = {
            'Muon' : {
                'objselection', 'jets_3' , 'tight_MU', 'btag_pass', 'semilep'
            },
            'Electron' : {
                'objselection', 'jets_3' , 'tight_EL', 'btag_pass', 'semilep'
            },
        }

            ## object selection
        objsel_evts = objsel.select(df, year=args.year, corrections=self.corrections, accumulator=output)
        output['cutflow']['nEvts passing jet and lepton obj selection'] += objsel_evts.sum()
        selection.add('jets_3', df['Jet'].counts == 3)
        selection.add('objselection', objsel_evts)
        selection.add('btag_pass', df['Jet'][btag_wp].sum() >= 2)

            ## add different selections
                ## muons
        selection.add('tight_MU', df['Muon']['TIGHTMU'].sum() == 1) # one muon passing TIGHT criteria
        #selection.add('loose_or_tight_MU', (df['Muon']['LOOSEMU'] | df['Muon']['TIGHTMU']).sum() == 1) # one muon passing LOOSE or TIGHT criteria

                ## electrons
        selection.add('tight_EL', df['Electron']['TIGHTEL'].sum() == 1) # one electron passing TIGHT criteria
        #selection.add('loose_or_tight_EL', (df['Electron']['LOOSEEL'] | df['Electron']['TIGHTEL']).sum() == 1) # one electron passing LOOSE or TIGHT criteria

        #set_trace()
        ### apply lepton SFs to MC (only applicable to tight leptons)
        if 'LeptonSF' in corrections.keys():
            tight_mu_cut = selection.require(objselection=True, tight_MU=True) # find events passing muon object selection with one tight muon
            tight_muons = df['Muon'][tight_mu_cut][(df['Muon'][tight_mu_cut]['TIGHTMU'] == True)]
            evt_weights._weights['Muon_SF'][tight_mu_cut] = MCWeights.get_lepton_sf(year=args.year, lepton='Muons', corrections=lepSF_correction,
                pt=tight_muons.pt.flatten(), eta=tight_muons.eta.flatten())
            tight_el_cut = selection.require(objselection=True, tight_EL=True) # find events passing electron object selection with one tight electron
            tight_electrons = df['Electron'][tight_el_cut][(df['Electron'][tight_el_cut]['TIGHTEL'] == True)]
            evt_weights._weights['Electron_SF'][tight_el_cut] = MCWeights.get_lepton_sf(year=args.year, lepton='Electrons', corrections=lepSF_correction,
                pt=tight_electrons.pt.flatten(), eta=tight_electrons.etaSC.flatten())

        if corrections['BTagSF'] == True:
            #set_trace()
            threeJets_cut = selection.require(objselection=True, jets_3=True)
            threeJets_btagwts = self.corrections['BTag_Constructors']['3Jets'].get_scale_factor(jets=df['Jet'][threeJets_cut], passing_cut=btag_wp)
            evt_weights._weights['Btag_SF'][threeJets_cut] = threeJets_btagwts['central'].prod()


            # find gen level particles for ttbar system
        #set_trace()
        genp_mode = 'NORMAL'
        GenTTbar = genpsel.select(df, systype='FINAL', mode=genp_mode)
        #genpsel.select(df, mode='LHE')
        #set_trace()
        selection.add('semilep', GenTTbar['SL']['TTbar'].counts > 0)

        ## fill hists for each region
        for lepton in regions.keys():
            cut = selection.all(*regions[lepton])

            #set_trace()
            if cut.sum() > 0:
                leptype = 'MU' if lepton == 'Muon' else 'EL'
                if 'loose_or_tight_%s' % leptype in regions[lepton]:
                    lep_mask = ((df[lepton][cut]['TIGHT%s' % leptype]) | (df[lepton][cut]['LOOSE%s' % leptype]))
                elif 'tight_%s' % leptype in regions[lepton]:
                    lep_mask = (df[lepton][cut]['TIGHT%s' % leptype])
                elif 'loose_%s' % leptype in regions[lepton]:
                    lep_mask = (df[lepton][cut]['LOOSE%s' % leptype])
                else:
                    raise ValueError("Not sure what lepton type to choose for event")

                matched_perm = ttmatcher.best_match(gen_hyp=GenTTbar[cut], jets=df['Jet'][cut], leptons=df[lepton][cut][lep_mask], met=df['MET'][cut])

                    ## create MT regions
                MT = make_vars.MT(df[lepton][cut][lep_mask], df['MET'][cut])
                MTHigh = (MT >= 40.).flatten()

                    # find best permutations
                best_perms = ttpermutator.find_best_permutations(jets=df['Jet'][cut], leptons=df[lepton][cut][lep_mask], MET=df['MET'][cut], btagWP=btag_wp)
                valid_perms = best_perms['TTbar'].counts > 0

                evt_weights_to_use = evt_weights.weight()
                ## apply lepton SFs to MC (only applicable to tight leptons)
                if 'LeptonSF' in corrections.keys():
                    evt_weights_to_use = evt_weights.partial_weight(exclude=['Electron_SF']) if lepton == 'Muon' else evt_weights.partial_weight(exclude=['Muon_SF']) # exclude SF from other lepton

                output = self.make_3j_categories(accumulator=output, leptype=lepton, genttbar=GenTTbar[cut][MTHigh], matched_perm=matched_perm[MTHigh], evt_weights=evt_weights_to_use[cut][MTHigh])

        return output


    def make_3j_categories(self, accumulator, leptype, genttbar, matched_perm, evt_weights):
        #unmatchable_evts = ~valid_evts
    
        valid_evts = (matched_perm['TTbar'].counts > 0) & ((matched_perm['unique_matches'] >= 3).flatten())

        correct_mp_merged = valid_evts & (matched_perm['Merged_Event'] & (matched_perm['Merged_BHadWJa'] | matched_perm['Merged_BHadWJb'] | matched_perm['Merged_WJets'])).flatten()
        correct_mp_lost   = valid_evts & (matched_perm['Lost_Event'] & (matched_perm['Lost_WJa'] | matched_perm['Lost_WJb'])).flatten()

        #set_trace()

        for evt_type in ['Lost', 'Merged']:
            mask = correct_mp_lost if evt_type == 'Lost' else correct_mp_merged
            for obj in ['TTbar', 'THad', 'TLep', 'BHad', 'BLep', 'WHad', 'WLep', 'Lepton']:
                accumulator['Gen_pt'    ].fill(dataset=self.sample_name, leptype=leptype, evttype=evt_type, objtype=obj, pt=np.array(genttbar['SL'][obj][mask].p4.pt.flatten()), weight=evt_weights[mask])
                accumulator['Gen_eta'   ].fill(dataset=self.sample_name, leptype=leptype, evttype=evt_type, objtype=obj, eta=np.array(genttbar['SL'][obj][mask].p4.eta.flatten()), weight=evt_weights[mask])
                accumulator['Gen_phi'   ].fill(dataset=self.sample_name, leptype=leptype, evttype=evt_type, objtype=obj, phi=np.array(genttbar['SL'][obj][mask].p4.phi.flatten()), weight=evt_weights[mask])
                accumulator['Gen_mass'  ].fill(dataset=self.sample_name, leptype=leptype, evttype=evt_type, objtype=obj, mass=np.array(genttbar['SL'][obj][mask].p4.mass.flatten()), weight=evt_weights[mask])
                accumulator['Gen_energy'].fill(dataset=self.sample_name, leptype=leptype, evttype=evt_type, objtype=obj, energy=np.array(genttbar['SL'][obj][mask].p4.energy.flatten()), weight=evt_weights[mask])

                accumulator['Reco_pt'    ].fill(dataset=self.sample_name, leptype=leptype, evttype=evt_type, objtype=obj, pt=np.array(matched_perm[obj][mask].p4.pt.flatten()), weight=evt_weights[mask])
                accumulator['Reco_eta'   ].fill(dataset=self.sample_name, leptype=leptype, evttype=evt_type, objtype=obj, eta=np.array(matched_perm[obj][mask].p4.eta.flatten()), weight=evt_weights[mask])
                accumulator['Reco_phi'   ].fill(dataset=self.sample_name, leptype=leptype, evttype=evt_type, objtype=obj, phi=np.array(matched_perm[obj][mask].p4.phi.flatten()), weight=evt_weights[mask])
                accumulator['Reco_mass'  ].fill(dataset=self.sample_name, leptype=leptype, evttype=evt_type, objtype=obj, mass=np.array(matched_perm[obj][mask].p4.mass.flatten()), weight=evt_weights[mask])
                accumulator['Reco_energy'].fill(dataset=self.sample_name, leptype=leptype, evttype=evt_type, objtype=obj, energy=np.array(matched_perm[obj][mask].p4.energy.flatten()), weight=evt_weights[mask])

                accumulator['Reso_pt'    ].fill(dataset=self.sample_name, leptype=leptype, evttype=evt_type, objtype=obj, reso_pt=np.array(genttbar['SL'][obj][mask].p4.pt.flatten())-np.array(matched_perm[obj][mask].p4.pt.flatten()), weight=evt_weights[mask])
                accumulator['Reso_eta'   ].fill(dataset=self.sample_name, leptype=leptype, evttype=evt_type, objtype=obj, reso_eta=np.array(genttbar['SL'][obj][mask].p4.eta.flatten())-np.array(matched_perm[obj][mask].p4.eta.flatten()), weight=evt_weights[mask])
                accumulator['Reso_phi'   ].fill(dataset=self.sample_name, leptype=leptype, evttype=evt_type, objtype=obj, reso_phi=np.array(genttbar['SL'][obj][mask].p4.phi.flatten())-np.array(matched_perm[obj][mask].p4.phi.flatten()), weight=evt_weights[mask])
                accumulator['Reso_mass'  ].fill(dataset=self.sample_name, leptype=leptype, evttype=evt_type, objtype=obj, reso_mass=np.array(genttbar['SL'][obj][mask].p4.mass.flatten())-np.array(matched_perm[obj][mask].p4.mass.flatten()), weight=evt_weights[mask])
                accumulator['Reso_energy'].fill(dataset=self.sample_name, leptype=leptype, evttype=evt_type, objtype=obj, reso_energy=np.array(genttbar['SL'][obj][mask].p4.energy.flatten())-np.array(matched_perm[obj][mask].p4.energy.flatten()), weight=evt_weights[mask])

        return accumulator

    def postprocess(self, accumulator):
        return accumulator

proc_executor = processor.iterative_executor if args.debug else processor.futures_executor

output = processor.run_uproot_job(fileset,
    treename='Events',
    processor_instance=ttbar_alpha_reco(),
    executor=proc_executor,
    executor_args={
        'workers': 8,
        'flatten' : True,
        'compression': 5,
    },
    chunksize=10000 if args.debug else 50000,
    #chunksize=50000,
)


if args.debug:
    print(output)
#set_trace()
#print(output['cutflow'])

save(output, args.outfname)
print('%s has been written' % args.outfname)
