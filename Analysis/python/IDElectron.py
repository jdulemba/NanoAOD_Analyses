import numpy as np
from pdb import set_trace
import Utilities.prettyjson as prettyjson
import os

def process_electrons(electrons, year):
        # makes etaSC
    electrons['etaSC'] = electrons['deltaEtaSC'] + electrons['eta']
    electrons['ECAL_GAP'] = (np.abs(electrons['etaSC']) <= 1.4442) | (np.abs(electrons['etaSC']) >= 1.5660)
    electrons['IPCuts'] = ((np.abs(electrons['etaSC']) < 1.479) & (np.abs(electrons['dxy']) < 0.05) & (np.abs(electrons['dz']) < 0.10)) | ((np.abs(electrons['etaSC']) >= 1.479) & (np.abs(electrons['dxy']) < 0.10) & (np.abs(electrons['dz']) < 0.20))

        ## make electron ID/Iso categories
    electrons = make_electron_ids(electrons, year)

    return electrons


    ## Electron ID types
'''
* Isolation values are initially taken from the EXO-19-016 paper (Section 4.1)
'''


def veto_15(electrons):
    ID = (electrons['cutBased'] == 1) # 0 is Fail, 1 Veto, 2 Loose, 3 Medium, 4 Tight
    return ID

def tight_15_NoECAL_Gap(electrons):
    ID = (electrons['cutBased'] == 4) # 0 is Fail, 1 Veto, 2 Loose, 3 Medium, 4 Tight
    ecalgap = (electrons['ECAL_GAP'])
    ipcuts = (electrons['IPCuts'])
    return (ID & ecalgap & ipcuts)

def fakes(electrons):
    '''
    For ID values:
        passing tight cutBased ID corresponds to 613566756 == 100 100 100 100 100 100 100 100 100 100 in binary
        I think the 3rd group of bits corresponds to isolation
        609372452 corresponds to 010 
        607275300 corresponds to 001
        605178148 corresponds to 000
        611469604 corresponds to 011 but has entries that also pass the isolation criteria so it's not included
        inverted Iso bitvals = [609372452, 607275300, 605178148]#, 611469604]
        
        invertedIso_bitvals = [613566692, 613566564, 613566628, 613566500]
    '''
    ID = ((electrons['vidNestedWPBitmap'] == 609372452) | (electrons['vidNestedWPBitmap'] == 607275300) | (electrons['vidNestedWPBitmap'] == 605178148)) #| (electrons.bitmap == 611469604))
    ecalgap = (electrons['ECAL_GAP'])
    ipcuts = (electrons['IPCuts'])
    return (ID & ecalgap & ipcuts)


def make_electron_ids(electrons, year):

    el_pars = prettyjson.loads(open(os.path.join(os.environ['PROJECT_DIR'], 'cfg_files', 'cfg_pars_%s.json' % os.environ['jobid'])).read())['Electrons'][year]

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

    for elID in el_pars.keys():
        pt_cut = (electrons['pt'] >= el_pars[elID]['ptmin'])
        #etaSC_cut = (np.abs(electrons['etaSC']) <= el_pars[elID]['etascmax'])
        eta_cut = (np.abs(electrons['eta']) <= el_pars[elID]['etamax'])
        pass_id = id_names[el_pars[elID]['id']](electrons)
        electrons[elID] = (pass_id) & (pt_cut) & (eta_cut)
        #electrons[elID] = (pass_id) & (pt_cut) & (etaSC_cut)

    return electrons


