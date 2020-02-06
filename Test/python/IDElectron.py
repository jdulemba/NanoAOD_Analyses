from coffea.analysis_objects import JaggedCandidateArray
import numpy as np
from pdb import set_trace
import Utilities.prettyjson as prettyjson
import os

def process_electrons(dataframe):

    electrons = JaggedCandidateArray.candidatesfromcounts(
        dataframe['nElectron'],
        pt=dataframe['Electron_pt'],
        eta=dataframe['Electron_eta'],
        phi=dataframe['Electron_phi'],
        mass=dataframe['Electron_mass'],
        charge=dataframe['Electron_charge'],
        cutBasedId=dataframe['Electron_cutBased'],
        dxy=dataframe['Electron_dxy'],
        dz=dataframe['Electron_dz'],
        deltaEtaSC=dataframe['Electron_deltaEtaSC'],
        pfRelIsoAll=dataframe['Electron_pfRelIso03_all']
        #trig=dataframe['HLT_Ele27_WPTight_Gsf'],
    )

        ## add attributes that must be computed
    eta_sc = np.add(electrons.deltaEtaSC, electrons.eta)
    electrons.add_attributes(
        etaSC = eta_sc
    )

    ecal_gap = (np.abs(electrons.etaSC) <= 1.4442) | (np.abs(electrons.etaSC) >= 1.5660)
    ip_cuts = ((np.abs(electrons.etaSC) < 1.479) & (np.abs(electrons.dxy) < 0.05) & (np.abs(electrons.dz) < 0.10)) | ((np.abs(electrons.etaSC) >= 1.479) & (np.abs(electrons.dxy) < 0.10) & (np.abs(electrons.dz) < 0.20))
    electrons.add_attributes(
        IPCuts = ip_cuts,
        ecalgap = ecal_gap,
    )


    return electrons

    ## various functions
def make_electron_id_functs(electrons):
    # makes etaSC
    electrons['etaSC'] = electrons.deltaEtaSC + electrons.eta

    # make IP cuts
    ipcuts = ( (np.abs(electrons.etaSC) < 1.479) & (np.abs(electrons.dxy) < 0.05) & (np.abs(electrons.dz) < 0.10) ) | ( (np.abs(electrons.etaSC) >= 1.479) & (np.abs(electrons.dxy) < 0.10) & (np.abs(electrons.dz) < 0.20) )
    electrons['IPCuts'] = ipcuts
    
    return electrons


    ## Electron ID types
#def fail(electrons):
#    return electrons
#
#def veto_15(electrons):
#
#def loose_15(electrons):
#
#def medium_15(electrons):
#
#def tight_15(electrons):
#
def tight_15_NoECAL_Gap(electrons):
    ID = (electrons.mvaFall17V2noIso_WP80)
    Iso = (electrons.pfRelIso03_all < 0.15) #???
    ecalgap = (np.abs(electrons.etaSC) <= 1.4442) | (np.abs(electrons.etaSC) >= 1.5660)
    ipcuts = (electrons.IPCuts)

    return (ID & Iso & ecalgap & ipcuts)

def fakes(electrons):
    ID = (electrons.mvaFall17V2noIso_WP80)
    Iso = (electrons.pfRelIso03_all >= 0.15) #???
    ecalgap = (np.abs(electrons.etaSC) <= 1.4442) | (np.abs(electrons.etaSC) >= 1.5660)
    ipcuts = (electrons.IPCuts)

    return (ID & Iso & ecalgap & ipcuts)


def make_electron_ids(electrons):

    el_pars = prettyjson.loads(open('%s/cfg_files/cfg_pars.json' % os.environ['PROJECT_DIR']).read())['Electrons']

    id_names = {
        #'FAIL' : fail,
        #'VETO_15' : veto_15,
        #'LOOSE_15' : loose_15,
        #'MEDIUM_15' : medium_15,
        #'TIGHT_15' : tight_15,
        'TIGHT_15_NoECAL_Gap' : tight_15_NoECAL_Gap,
        #'NOTVETO_15' : notveto_15,
        'FAKES' : fakes
    }

    #if el_pars['VETOEL']['id'] not in id_names.keys():
    #    raise IOError("veto Electron ID name not valid")
    if el_pars['LOOSEEL']['id'] not in id_names.keys():
        raise IOError("loose Electron ID name not valid")
    if el_pars['TIGHTEL']['id'] not in id_names.keys():
        raise IOError("tight Electron ID name not valid")

    #set_trace()
    for elID in el_pars.keys():
        pt_cut = (electrons.pt >= el_pars[elID]['ptmin'])
        etaSC_cut = (np.abs(electrons.etaSC) <= el_pars[elID]['etascmax'])
        pass_id = id_names[el_pars[elID]['id']](electrons)
        electrons[elID] = (pass_id) & (pt_cut) & (etaSC_cut)

    return electrons


def build_electrons(electrons):
    electrons = make_electron_id_functs(electrons)
    electrons = make_electron_ids(electrons)

    return electrons
