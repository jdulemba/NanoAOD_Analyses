from coffea.arrays import Initialize

def SetMET(df):
    met = Initialize({
        'pt' : df['MET_pt'],
        'eta' : 0,
        'phi' : df['MET_phi'],
        'mass' : 0
    })

    return met
