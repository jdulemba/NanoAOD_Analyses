import numpy as np
from pdb import set_trace
import Utilities.prettyjson as prettyjson
import os

def process_electrons(df, year):
    from coffea.analysis_objects import JaggedCandidateArray
    electrons = JaggedCandidateArray.candidatesfromcounts(
        df['nElectron'],
        pt=df['Electron_pt'],
        eta=df['Electron_eta'],
        phi=df['Electron_phi'],
        mass=df['Electron_mass'],
        charge=df['Electron_charge'],
        dxy=df['Electron_dxy'],
        dz=df['Electron_dz'],
        tightID=df['Electron_mvaFall17V2noIso_WP80'],
        vetoID=df['Electron_mvaFall17V2Iso_WPL'],
        #vetoID=df['Electron_mvaFall17V2noIso_WPL'],
        deltaEtaSC=df['Electron_deltaEtaSC'],
        pfRelIso=df['Electron_pfRelIso03_all'],
        cutBasedID=df['Electron_cutBased'],
        bitmap = df['Electron_vidNestedWPBitmap'],
    )

        # makes etaSC
    electrons['etaSC'] = electrons.deltaEtaSC + electrons.eta
    electrons['ECAL_GAP'] = (np.abs(electrons.etaSC) <= 1.4442) | (np.abs(electrons.etaSC) >= 1.5660)

    electrons['IPCuts'] = ((np.abs(electrons.etaSC) < 1.479) & (np.abs(electrons.dxy) < 0.05) & (np.abs(electrons.dz) < 0.10)) | ((np.abs(electrons.etaSC) >= 1.479) & (np.abs(electrons.dxy) < 0.10) & (np.abs(electrons.dz) < 0.20))

    #set_trace()
        ## make electron ID/Iso categories
    electrons = make_electron_ids(electrons, year)

    return electrons


    ## Electron ID types
'''
* Isolation values are initially taken from the EXO-19-016 paper (Section 4.1)
'''


#def fail(electrons):
#    return electrons
#
def veto_15(electrons):
    ID = (electrons.cutBasedID == 1) # 0 is Fail, 1 Veto, 2 Loose, 3 Medium, 4 Tight
    #ID = (electrons.vetoID)
    ##Iso = (electrons.pfRelIso < 0.25) # based on muon loose Iso def
    return ID

def tight_15_NoECAL_Gap(electrons):
    ID = (electrons.cutBasedID == 4) # 0 is Fail, 1 Veto, 2 Loose, 3 Medium, 4 Tight
    #ID = (electrons.tightID)
    #Iso = (electrons.pfRelIso < 0.15) #*
    ecalgap = (electrons.ECAL_GAP)
    ipcuts = (electrons.IPCuts)

    return (ID & ecalgap & ipcuts)
    #return (ID & Iso & ecalgap & ipcuts)

def fakes(electrons):
    '''
    For ID values:
        613566692 corresponds to Iso of 011 (pass medium and loose but not tight) part of binary but 100 for other cuts
        613566564 corresponds to 001 for Iso but 100 for other cuts
        613566628 corresponds to 010 for Iso but 100 for other cuts
        passing tight cutBased ID corresponds to 613566756 == 100 100 100 100 100 100 100 100 100 100 in binary
    '''
    ID = ((electrons.bitmap == 613566692) | (electrons.bitmap == 613566564) | (electrons.bitmap == 613566628)) 
    #ID = (electrons.tightID)
    #Iso = (electrons.pfRelIso >= 0.15) #*
    ecalgap = (electrons.ECAL_GAP)
    ipcuts = (electrons.IPCuts)

    return (ID & ecalgap & ipcuts)
    #return (ID & Iso & ecalgap & ipcuts)


def make_electron_ids(electrons, year):

    el_pars = prettyjson.loads(open('%s/cfg_files/cfg_pars_%s.json' % (os.environ['PROJECT_DIR'], os.environ['jobid'])).read())['Electrons'][year]

    id_names = {
        #'FAIL' : fail,
        'VETO_15' : veto_15,
        #'LOOSE_15' : loose_15,
        #'MEDIUM_15' : medium_15,
        #'TIGHT_15' : tight_15,
        'TIGHT_15_NoECAL_Gap' : tight_15_NoECAL_Gap,
        #'NOTVETO_15' : notveto_15,
        'FAKES' : fakes
    }

    if el_pars['VETOEL']['id'] not in id_names.keys():
        raise IOError("veto Electron ID name not valid")
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


