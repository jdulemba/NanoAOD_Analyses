from pdb import set_trace
import coffea.processor.dataframe
import numpy as np

def get_prefire_weights(df, shift=None, mask=None):
    if not isinstance(df, coffea.processor.dataframe.LazyDataFrame):
        raise IOError("This function only works for LazyDataFrame objects")

    #set_trace()
    if shift == 'PREFIRE_UP':
        weights = df['L1PreFiringWeight_Up'][(mask)]
    elif shift == 'PREFIRE_DW':
        weights = df['L1PreFiringWeight_Dn'][(mask)]
    else:
        weights = df['L1PreFiringWeight_Nom'][(mask)]

    return weights


def get_gen_weights(df, shift=None, mask=None):
    'LHEScaleWeight definitions can be found here: https://cms-nanoaod-integration.web.cern.ch/integration/master/mc94X_doc.html#LHE\nPSWeight definitions can be found here: https://cms-nanoaod-integration.web.cern.ch/integration/master/mc94X_doc.html#PSWeight'

    gen_weight = df.genWeight[(mask)]

        ## LHEScale Weight Variations
    if shift == 'FACTOR_DW':
        weight = df['LHEScaleWeight'][3][(mask)] # (muR=1, muF=0.5)
    elif shift == 'FACTOR_UP':
        weight = df['LHEScaleWeight'][5][(mask)] # (muR=1, muF=2)
    if shift == 'RENORM_DW':
        weight = df['LHEScaleWeight'][1][(mask)] # (muR=0.5, muF=1)
    elif shift == 'RENORM_UP':
        weight = df['LHEScaleWeight'][7][(mask)] # (muR=2, muF=1)
    elif shift == 'RENFACTOR_DW':
        weight = df['LHEScaleWeight'][0][(mask)] # (muR=0.5, muF=0.5)
    elif shift == 'RENFACTOR_UP':
        weight = df['LHEScaleWeight'][8][(mask)] # (muR=2, muF=2)
    elif shift == 'RENORM_UP_FACTOR_DW':
        weight = df['LHEScaleWeight'][6][(mask)] # (muR=2, muF=0.5)
    elif shift == 'RENORM_DW_FACTOR_UP':
        weight = df['LHEScaleWeight'][2][(mask)] # (muR=0.5, muF=2)

        ## PS Weight variations
    elif shift == 'ISR_DW':
        weight = df['PSWeight'][0][(mask)] # (ISR=0.5, FSR=1)
    elif shift == 'ISR_UP':
        weight = df['PSWeight'][2][(mask)] # (ISR=2, FSR=1)
    if shift == 'FSR_DW':
        weight = df['PSWeight'][1][(mask)] # (ISR=1, FSR=0.5)
    elif shift == 'FSR_UP':
        weight = df['PSWeight'][3][(mask)] # (ISR=1, FSR=2)
    else: #nominal
        weight = 1.

    return gen_weight*weight        


def get_pu_weights(df, shift=None, mask=None):
    'Pileup info defined here: https://cms-nanoaod-integration.web.cern.ch/integration/master/mc94X_doc.html#Pileup'

    weight = df['Pileup_nTrueInt']

    return weight


def get_toppt_weights(pt1=np.array([-1.]), pt2=np.array([-1.]), shift=None):
    if pt1.any() < 0. or pt2.any() < 0.:
        raise IOError("top pt inputs not valid")

    if shift == 'TOPPT_UP_NOM':
        nu_1, nu_2 = 1., 0.
    elif shift == 'TOPPT_DW_NOM':
        nu_1, nu_2 = -1., 0.
    elif shift == 'TOPPT_NOM_UP':
        nu_1, nu_2 = 0., 1.
    elif shift == 'TOPPT_NOM_DW':
        nu_1, nu_2 = 0., -1.
    else:
        nu_1, nu_2 = 0., 0.

    # top pT reweighting formula based on section 5.5 of AN-16-272
    p0 =  6.15025*np.power(10., -2.) + 3.243*np.power(10., -2.)*nu_1 - 4.353*np.power(10., -7.)*nu_2
    p1 = -5.17833*np.power(10., -4.) - 1.404*np.power(10., -4.)*nu_1 - 1.005*np.power(10., -4.)*nu_2

    weight = exp(p0+p1*( (pt1+pt2)/2 ))
    return weight


def evt_weight(df, shift=None, mask=None, use_weight=True):
    #set_trace()
    if use_weight:
        gen_weights = get_gen_weights(df, shift, mask)
    else:
        gen_weights = np.ones(len(mask)) if mask else np.ones(df.size)

    pref_weights = get_prefire_weights(df, shift, mask)
    #pu_weights = get_pu_weights(df, shift, mask)
    evt_weights = pref_weights*gen_weights#*pu_weight
    return evt_weights
