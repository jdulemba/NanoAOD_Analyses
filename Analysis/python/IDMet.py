import numpy as np
#from pdb import set_trace

def process_met(df):
    from coffea.analysis_objects import JaggedCandidateArray
    met = JaggedCandidateArray.candidatesfromcounts(
        counts=np.ones(df.size, dtype=int),
        pt=df['MET_pt'],
        eta=0,
        phi=df['MET_phi'],
        mass=0,
        MetUnclustEnUpDeltaX=df['MET_MetUnclustEnUpDeltaX'],
        MetUnclustEnUpDeltaY=df['MET_MetUnclustEnUpDeltaY'],
        covXX=df['MET_covXX'],
        covXY=df['MET_covXY'],
        covYY=df['MET_covYY'],
        sumEt=df['MET_sumEt'],
    )

    return met

