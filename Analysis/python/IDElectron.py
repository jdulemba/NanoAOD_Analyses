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
        mvaTightID=df['Electron_mvaFall17V2noIso_WP80'],
        mvaVetoID=df['Electron_mvaFall17V2noIso_WPL'],
        #mvaVetoID=df['Electron_mvaFall17V2noIso_WPL'],
        deltaEtaSC=df['Electron_deltaEtaSC'],
        pfRelIso=df['Electron_pfRelIso03_all'],
        minipfRelIso=df['Electron_miniPFRelIso_all'],
        cutBasedID=df['Electron_cutBased'],
        bitmap=df['Electron_vidNestedWPBitmap'],
        sieie=df['Electron_sieie'],
        hoe=df['Electron_hoe'],
        eInvMpInv=df['Electron_eInvMinusPInv'],
        missingHits=df['Electron_lostHits'],
        convVeto=df['Electron_convVeto'],
    )

        # makes etaSC
    electrons['etaSC'] = electrons.deltaEtaSC + electrons.eta
    electrons['ECAL_GAP'] = (np.abs(electrons.etaSC) <= 1.4442) | (np.abs(electrons.etaSC) >= 1.5660)

    electrons['IPCuts'] = ((np.abs(electrons.etaSC) < 1.479) & (np.abs(electrons.dxy) < 0.05) & (np.abs(electrons.dz) < 0.10)) | ((np.abs(electrons.etaSC) >= 1.479) & (np.abs(electrons.dxy) < 0.10) & (np.abs(electrons.dz) < 0.20))


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
    #ID = (electrons.mvaVetoID)
    ##Iso = (electrons.pfRelIso < 0.25) # based on muon loose Iso def
    return ID

def tight_15_NoECAL_Gap(electrons):
    ID = (electrons.cutBasedID == 4) # 0 is Fail, 1 Veto, 2 Loose, 3 Medium, 4 Tight
    #ID = (electrons.mvaTightID)
    #Iso = (electrons.pfRelIso < 0.15) #*
    ecalgap = (electrons.ECAL_GAP)
    ipcuts = (electrons.IPCuts)

    return (ID & ecalgap & ipcuts)
    #return (ID & Iso & ecalgap & ipcuts)

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
    #    # testing electron iso inversion
    #basic_cut = (electrons.pt > 30) & (np.abs(electrons.eta) < 2.4) & (electrons.ECAL_GAP) & (electrons.IPCuts)
    #els = electrons[basic_cut]

    #barrel = (np.abs(els.etaSC) <= 1.479)
    #endcap = (np.abs(els.etaSC) > 1.479)
    #bitmap = els.bitmap

    #rho = els.pt.ones_like()*df['fixedGridRhoFastjetAll']
    #tightID_barrel = barrel & (els.sieie < 0.0104) & (els.hoe < 0.026+1.15/els.etaSC+0.0324*rho/els.etaSC) & (np.abs(els.eInvMpInv) < 0.159) & (els.missingHits <= 1) & els.convVeto
    #tightID_endcap = endcap & (els.sieie < 0.0353) & (els.hoe < 0.0188+2.06/els.etaSC+0.183*rho/els.etaSC) & (np.abs(els.eInvMpInv) < 0.0197) & (els.missingHits <= 1) & els.convVeto
    #tightID_els = tightID_barrel | tightID_endcap

    #barrel_els = els[tightID_barrel]
    #endcap_els = els[tightID_endcap]
    #tightIso_barrel = (barrel_els.pfRelIso < (0.0287 + 0.506/barrel_els.pt))
    #tightIso_endcap = (endcap_els.pfRelIso < (0.0445 + 0.963/endcap_els.pt))
    #tightIso_bitvals = set(bitmap[tightID_barrel][tightIso_barrel].flatten().tolist()+bitmap[tightID_endcap][tightIso_endcap].flatten().tolist())
    ##tightIso_bitvals = list(set(bitmap[tightID_barrel][tightIso_barrel].flatten().tolist()+bitmap[tightID_endcap][tightIso_endcap].flatten().tolist()))
    #invTightIso_barrel = (barrel_els.pfRelIso >= (0.0287 + 0.506/barrel_els.pt))
    #invTightIso_endcap = (endcap_els.pfRelIso >= (0.0445 + 0.963/endcap_els.pt))
    #invTightIso_bitvals = set(bitmap[tightID_barrel][invTightIso_barrel].flatten().tolist()+bitmap[tightID_endcap][invTightIso_endcap].flatten().tolist())
    ##invTightIso_bitvals = list(set(bitmap[tightID_barrel][invTightIso_barrel].flatten().tolist()+bitmap[tightID_endcap][invTightIso_endcap].flatten().tolist()))

    #only_inv_bits = list(invTightIso_bitvals - tightIso_bitvals)
    #import re
    #binary_nums = ['{0:0b}'.format(num) for num in only_inv_bits]
    #split_binary_nums = [re.findall('...', num) for num in binary_nums]
    #diff_inds = [[idx for idx, val in enumerate(re.findall('...', num)) if val != '100'] for num in binary_nums]
    #flat_list = [item for sublist in diff_inds for item in sublist]
    #flat_set = set(flat_list)
    #set_trace()
    #invest_barrel_els = els[(bitmap == 611469604) & (np.abs(els.etaSC) <= 1.479)]
    #invest_endcap_els = els[(bitmap == 611469604) & (np.abs(els.etaSC) > 1.479)]

    #invIso_barrel = invest_barrel_els.pfRelIso < 0.0287 + 0.506/invest_barrel_els.pt
    #invIso_endcap = invest_endcap_els.pfRelIso < 0.0445 + 0.963/invest_endcap_els.pt
        ##

    ID = ((electrons.bitmap == 609372452) | (electrons.bitmap == 607275300) | (electrons.bitmap == 605178148)) #| (electrons.bitmap == 611469604))
    #ID = (electrons.mvaTightID)
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
        #etaSC_cut = (np.abs(electrons.etaSC) <= el_pars[elID]['etascmax'])
        eta_cut = (np.abs(electrons.eta) <= el_pars[elID]['etamax'])
        pass_id = id_names[el_pars[elID]['id']](electrons)
        electrons[elID] = (pass_id) & (pt_cut) & (eta_cut)
        #electrons[elID] = (pass_id) & (pt_cut) & (etaSC_cut)

    return electrons


