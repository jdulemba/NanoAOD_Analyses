from coffea.analysis_objects import JaggedCandidateArray
from pdb import set_trace

def process_jets(dataframe):
    jets = JaggedCandidateArray.candidatesfromcounts(
        dataframe['nJet'],
        pt=dataframe['Jet_pt'],
        eta=dataframe['Jet_eta'],
        phi=dataframe['Jet_phi'],
        mass=dataframe['Jet_mass'],
        DeepCSVb=dataframe['Jet_btagDeepB'],
        DeepJetb=dataframe['Jet_btagDeepFlavB'],
        Id=dataframe['Jet_jetId'],
        cleanmask=dataframe['Jet_cleanmask'],
    )

    #tight_id = (jets.Id == 3) # 3 for passing loose+tight for 2016, other years are different https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookNanoAOD

        ## add btag wps
    deepcsvL = (jets.DeepCSVb > 0.2217)
    deepcsvM = (jets.DeepCSVb > 0.6321)
    deepcsvT = (jets.DeepCSVb > 0.8953)

    deepjetL = (jets.DeepJetb > 0.0614)
    deepjetM = (jets.DeepJetb > 0.3093)
    deepjetT = (jets.DeepJetb > 0.7221)

    jets.add_attributes(
        #tightID = tight_id,
        deepcsvL = deepcsvL,
        deepcsvM = deepcsvM,
        deepcsvT = deepcsvT,
        deepjetL = deepjetL,
        deepjetM = deepjetM,
        deepjetT = deepjetT,
    )

    return jets

