from coffea.analysis_objects import JaggedCandidateArray
import numpy as np
from pdb import set_trace

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

    #TIGHT_15_NoECAL_Gap = 

    return electrons

