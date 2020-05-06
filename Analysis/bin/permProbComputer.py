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

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'permProbComputer'

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
class permProbComputer(processor.ProcessorABC):
    def __init__(self):

            ## make binning for hists
        self.dataset_axis = hist.Cat("dataset", "Event Process")
        self.jetmult_axis = hist.Cat("jmult", "nJets")
        self.leptype_axis = hist.Cat("leptype", "Lepton Type")
        self.tMass_axis = hist.Bin("topmass", "m(t_{had}) [GeV]", 500, 0., 500.)
        self.wMass_axis = hist.Bin("wmass", "m(W_{had}) [GeV]", 500, 0., 500.)
        self.nu_chi2_axis = hist.Bin("nu_chi2", r"$\nu$ $\chi^{2}$", 1000, 0., 1000.)
        self.nu_dist_axis = hist.Bin("nu_dist", r"$\nu$ dist", 150, 0., 150.)

            ## make dictionary of hists
        histo_dict = {}
                ## make jet hists
        perm_4pj_hists = self.make_4pj_hists()
        histo_dict.update(perm_4pj_hists)

        histo_dict['cutflow'] = processor.defaultdict_accumulator(int)

        self._accumulator = processor.dict_accumulator(histo_dict)
        self.sample_name = ''
        self.corrections = corrections

    
    @property
    def accumulator(self):
        return self._accumulator


    def make_4pj_hists(self):
        histo_dict = {}
        histo_dict['mWHad_vs_mTHad'] = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.tMass_axis, self.wMass_axis)
        histo_dict['nusolver_chi2']  = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.nu_chi2_axis)
        histo_dict['nusolver_dist']  = hist.Hist("Events", self.dataset_axis, self.jetmult_axis, self.leptype_axis, self.nu_dist_axis)

        return histo_dict

    


    def process(self, df):
        np.random.seed(10) # sets seed so values from random distributions are reproducible (JER corrections)
        output = self.accumulator.identity()

        #if not isinstance(df, coffea.processor.dataframe.LazyDataFrame):
        #    raise ValueError("This function only works for LazyDataFrame objects")

        #if args.debug: set_trace()
        self.sample_name = df.dataset

            ## make event weights
                # data or MC distinction made internally
        evt_weights = MCWeights.get_event_weights(df, year=args.year, corrections=self.corrections)

            ## initialize selections and regions
        selection = processor.PackedSelection()
        regions = {
            'Muon' : {
                #'3Jets'  : {'objselection', 'jets_3' , 'loose_or_tight_MU', 'btag_pass', 'semilep'},
                '4PJets' : {'objselection', 'jets_4p', 'loose_or_tight_MU', 'btag_pass', 'semilep'},
            },
            'Electron' : {
                #'3Jets'  : {'objselection', 'jets_3' , 'loose_or_tight_EL', 'btag_pass', 'semilep'},
                '4PJets' : {'objselection', 'jets_4p', 'loose_or_tight_EL', 'btag_pass', 'semilep'},
            },
        }

            ## object selection
        objsel_evts = objsel.select(df, year=args.year, corrections=self.corrections, accumulator=output)
        output['cutflow']['nEvts passing jet and lepton obj selection'] += objsel_evts.sum()
        selection.add('jets_3', df['Jet'].counts == 3)
        selection.add('jets_4p', df['Jet'].counts > 3)
        selection.add('objselection', objsel_evts)
        selection.add('btag_pass', df['Jet'][btag_wp].sum() >= 2)

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
            fourplusJets_cut = selection.require(objselection=True, jets_4p=True)
            fourplusJets_btagwts = self.corrections['BTag_Constructors']['4PJets'].get_scale_factor(jets=df['Jet'][fourplusJets_cut], passing_cut=btag_wp)
            evt_weights._weights['Btag_SF'][fourplusJets_cut] = fourplusJets_btagwts['central'].prod()


            # find gen level particles for ttbar system
        #set_trace()
        genp_mode = 'NORMAL'
        GenTTbar = genpsel.select(df, systype='FINAL', mode=genp_mode)
        #genpsel.select(df, mode='LHE')
        #set_trace()
        selection.add('semilep', GenTTbar['SL']['TTbar'].counts > 0)

        ## fill hists for each region
        for lepton in regions.keys():
            for jmult in regions[lepton].keys():
                cut = selection.all(*regions[lepton][jmult])

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

                    matched_perm = ttmatcher.best_match(gen_hyp=GenTTbar[cut], jets=df['Jet'][cut], leptons=df[lepton][cut][lep_mask], met=df['MET'][cut])
                    #set_trace()

                    evt_weights_to_use = evt_weights.weight()
                    ## apply lepton SFs to MC (only applicable to tight leptons)
                    if 'LeptonSF' in corrections.keys():
                        evt_weights_to_use = evt_weights.partial_weight(exclude=['Electron_SF']) if lepton == 'Muon' else evt_weights.partial_weight(exclude=['Muon_SF']) # exclude SF from other lepton
                    output = self.fill_4pj_hists(accumulator=output, jetmult=jmult, leptype=lepton, perm=matched_perm, evt_weights=evt_weights_to_use[cut])

        return output

    def fill_4pj_hists(self, accumulator, jetmult, leptype, perm, evt_weights):
        #set_trace()
        valid_evts = (perm['TTbar'].counts > 0) & ((perm['unique_matches'] == 4).flatten())
        accumulator['mWHad_vs_mTHad'].fill(dataset=self.sample_name, jmult=jetmult, leptype=leptype, topmass=np.array(perm['THad'][valid_evts].mass.flatten()), wmass=np.array(perm['WHad'][valid_evts].mass.flatten()), weight=np.array((perm['Nu'][valid_evts].pt.ones_like()*evt_weights[valid_evts]).flatten()))
        accumulator['nusolver_chi2'].fill( dataset=self.sample_name, jmult=jetmult, leptype=leptype, nu_chi2=np.array(perm['Nu'][valid_evts].chi2.flatten()), weight=np.array((perm['Nu'][valid_evts].pt.ones_like()*evt_weights[valid_evts]).flatten()))
        accumulator['nusolver_dist'].fill( dataset=self.sample_name, jmult=jetmult, leptype=leptype, nu_dist=np.array(np.sqrt(perm['Nu'][valid_evts].chi2.flatten())), weight=np.array((perm['Nu'][valid_evts].pt.ones_like()*evt_weights[valid_evts]).flatten()))

        return accumulator        


    def postprocess(self, accumulator):
        return accumulator

proc_executor = processor.iterative_executor if args.debug else processor.futures_executor

output = processor.run_uproot_job(fileset,
    treename='Events',
    processor_instance=permProbComputer(),
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
