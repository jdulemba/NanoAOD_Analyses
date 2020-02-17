from pdb import set_trace
import coffea.processor.dataframe
import numpy as np

def prefire_weight(df, shift=None, mask=None):
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

'''
def gen_weight(df, shift=None, mask=None):

def pu_weight(df, shift=None, mask=None):


def evt_weight(df, shift=None, mask=None):
'''
def toppt_weight(pt1=np.array([-1.]), pt2=np.array([-1.]), shift=None):
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
