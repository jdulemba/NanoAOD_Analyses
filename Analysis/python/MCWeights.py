from pdb import set_trace
import coffea.processor.dataframe
import numpy as np
import coffea.processor as processor

def get_event_weights(df, year: str, corrections, BTagSFs = []):
    weights = processor.Weights(df.size, storeIndividual=True)# store individual variations

        ## only apply to MC
    if not df.dataset.startswith('data_Single'):
            ## Prefire Corrections
        if (year != '2018') and (corrections['Prefire'] == True):
            weights.add('prefire_weight',
                df['L1PreFiringWeight_Nom'],
                df['L1PreFiringWeight_Up'],
                df['L1PreFiringWeight_Dn']
            )

            ## Generator Weights (normalize them)
        genWeights = np.ones(df.genWeight.size)
        genWeights[df.genWeight < 0] = -1.
        genWeights[df.genWeight == 0] = 0.
        weights.add('genweight', genWeights)

            ## Initialize Lepton Scale Factors
        if 'LeptonSF' in corrections.keys():
            weights.add('Muon_SF',
                np.ones(df.genWeight.size),
                np.ones(df.genWeight.size),
                np.ones(df.genWeight.size),
                shift=True # makes up/down variations relative to nominal
            )
            weights.add('Electron_SF',
                np.ones(df.genWeight.size),
                np.ones(df.genWeight.size),
                np.ones(df.genWeight.size),
                shift=True # makes up/down variations relative to nominal
            )
    
            ## Pileup Reweighting
        if 'Pileup' in corrections.keys():
            weights.add('pileup_weight',
                corrections['Pileup'][year][df['dataset']]['central'](df['Pileup_nTrueInt']),
                corrections['Pileup'][year][df['dataset']]['up'](df['Pileup_nTrueInt']),
                corrections['Pileup'][year][df['dataset']]['down'](df['Pileup_nTrueInt'])
            )
    
            ## BTag SFs
        if corrections['BTagSF'] == True:
            if BTagSFs:
                for btagSF in BTagSFs:
                    weights.add(btagSF, np.ones(df.genWeight.size))
            else:
                weights.add('Btag_SF',
                    np.ones(df.genWeight.size),
                    #np.ones(df.genWeight.size),
                    #np.ones(df.genWeight.size),
                    #shift=True # makes up/down variations relative to nominal
                )
    ## Need to add at some point
            ## LHEScale Weight Variations
            ## PS Weight variations

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


def get_lepton_sf(year: str, lepton: str, corrections, pt: np.ndarray, eta: np.ndarray, shift='Central'):
    if not (shift == 'Central' or shift == 'Error'):
        raise ValueError('Shift value %s not defined' % shift)
    sf_dict = corrections[year][lepton]
    eta_ranges = sf_dict['eta_ranges']
    lepSFs = np.ones(pt.size)

    for idx, eta_range in enumerate(eta_ranges):
        mask = (eta >= eta_range[0]) & (eta < eta_range[1]) # find inds that are within given eta range
        if not mask.any(): continue # no values fall within eta range
        #set_trace()
        recoSF = sf_dict['Reco_ID'][shift]['eta_bin%i' % idx]
        trigSF = sf_dict['Trig'][shift]['eta_bin%i' % idx]
        lepSFs[mask] = recoSF(pt[mask])*trigSF(pt[mask])

    return lepSFs
