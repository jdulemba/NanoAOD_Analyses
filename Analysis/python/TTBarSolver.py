'''
This file computes the probabilities from mass and neutrino solver distributions (and their sums)
for different jet multiplicities (4+ jets, 3 jets lost, 3 jets merged)
'''

from coffea.lookup_tools.root_converters import convert_histo_root_file
from coffea.lookup_tools.dense_lookup import dense_lookup
import Utilities.prettyjson as prettyjson
from pdb import set_trace
import os
import numpy as np
from numba import njit

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
cfg_pars = prettyjson.loads(open('%s/cfg_files/cfg_pars.json' % proj_dir).read())['ttsolver']

year = '2016'
filename = cfg_pars['filename']
dirname = cfg_pars['dirname']
probs = convert_histo_root_file('%s/inputs/%s_%s/INPUTS/%s' % (proj_dir, year, jobid, filename))

    ## create arrays for binning and values separately for each dist because njit can't handle constant dictionaries currently
if cfg_pars['USEMASS']:
    USEMASS = True
    WTmass_right = dense_lookup(*probs[('%s/mWhad_vs_mtophad_right' % dirname, 'dense_lookup')])
    WTmass_right_binning = WTmass_right._axes
    WTmass_right_values = WTmass_right._values

if cfg_pars['USE3JMERGED']:
    USE3JMERGED = True
    Mass_3J_Merged_right = dense_lookup(*probs[('%s/3J_mbpjet_vs_maxmjet_merged_right' % dirname, 'dense_lookup')])
    Mass_3J_Merged_right_binning = Mass_3J_Merged_right._axes
    Mass_3J_Merged_right_values = Mass_3J_Merged_right._values

if cfg_pars['USE3JLOST']:
    USE3JLOST = True
    Mass_3J_Lost_right = dense_lookup(*probs[('%s/3J_mbpjet_lost_right' % dirname, 'dense_lookup')])
    Mass_3J_Lost_right_binning = Mass_3J_Lost_right._axes
    Mass_3J_Lost_right_values = Mass_3J_Lost_right._values

if cfg_pars['USENS']:
    USENS = True
    NS_4PJ_right = dense_lookup(*probs[('%s/nusolver_chi2_right' % dirname, 'dense_lookup')])
    NS_4PJ_right_binning = NS_4PJ_right._axes
    NS_4PJ_right_values = NS_4PJ_right._values

    # merged 3 jet vars
    NS_3J_Merged_right = dense_lookup(*probs[('%s/3J_nusolver_chi2_merged_right' % dirname, 'dense_lookup')])
    NS_3J_Merged_right_binning = NS_3J_Merged_right._axes
    NS_3J_Merged_right_values = NS_3J_Merged_right._values

    # lost 3 jet vars
    NS_3J_Lost_right = dense_lookup(*probs[('%s/3J_nusolver_chi2_lost_right' % dirname, 'dense_lookup')])
    NS_3J_Lost_right_binning = NS_3J_Lost_right._axes
    NS_3J_Lost_right_values = NS_3J_Lost_right._values


@njit()
def solve_4PJ(mthad, mwhad, nschi):
    'Inputs: m(thad), m(whad), and the NS solver chi2 values (mthad, mwhad, nschi)\nOutputs: tuple of the NS, mass, and combined disriminants (NuDiscr, MassDiscr, Prob)'
    ## initialize values used to build up likelihood
    result = 0.
    nstest = np.inf
    masstest = np.inf

        ## evaluate probabilities
    if USENS:
            ## find binning to use for when sqrt(nschi) is in under/overflow
        if np.sqrt(nschi) > NS_4PJ_right_binning[-1]:
            bin_to_use = int(NS_4PJ_right_binning[-1] - 1)
        elif np.sqrt(nschi) < NS_4PJ_right_binning[0]:
            bin_to_use = 0
        else:
            bin_to_use = np.argmax(NS_4PJ_right_binning > np.sqrt(nschi)) - 1 ## -1 since binning is edges, not center

        nsdisval = NS_4PJ_right_values[bin_to_use]/np.sum(NS_4PJ_right_values) # normalizes dist based on integral
        if nsdisval > 1.0E-10:
            nstest = -1*np.log(nsdisval)
        result += nstest

    if USEMASS:
            ## find x binning to use for when mwhad is in under/overflow
        if mwhad > WTmass_right_binning[0][-1]:
            xbin_to_use = int(WTmass_right_binning[0][-1] - 1)
        elif mwhad < WTmass_right_binning[0][0]:
            xbin_to_use = 0
        else:
            xbin_to_use = np.argmax(WTmass_right_binning[0] > mwhad) - 1

            ## find y binning to use for when mwhad is in under/overflow
        if mthad > WTmass_right_binning[1][-1]:
            ybin_to_use = int(WTmass_right_binning[1][-1] - 1)
        elif mthad < WTmass_right_binning[1][0]:
            ybin_to_use = 0
        else:
            ybin_to_use = np.argmax(WTmass_right_binning[1] > mthad) - 1
        
        massdisval = WTmass_right_values[xbin_to_use][ybin_to_use]/np.sum(WTmass_right_values) # normalizes dist based on integral
        if massdisval > 1.0E-10:
            masstest = -1*np.log(massdisval)
        result += masstest

    return nstest, masstest, result

@njit()
def solve_3J_lost(mbpjet, nschi):
    'Inputs: m(b+jet) and the NS solver chi2 values (mbpjet, nschi)\nOutputs: tuple of the NS, mass, and combined disriminants (NuDiscr, MassDiscr, Prob)'
        ## initialize values used to build up likelihood
    result = 0.
    nstest = np.inf
    masstest = np.inf

        ## evaluate probabilities
    if USENS:
            ## find binning to use for when sqrt(nschi) is in under/overflow
        if np.sqrt(nschi) > NS_3J_Lost_right_binning[-1]:
            bin_to_use = int(NS_3J_Lost_right_binning[-1] - 1)
        elif np.sqrt(nschi) < NS_3J_Lost_right_binning[0]:
            bin_to_use = 0
        else:
            bin_to_use = np.argmax(NS_3J_Lost_right_binning > np.sqrt(nschi)) - 1 ## -1 since binning is edges, not center

        nsdisval = NS_3J_Lost_right_values[bin_to_use]/np.sum(NS_3J_Lost_right_values) # normalizes dist based on integral
        if nsdisval > 1.0E-10:
            nstest = -1*np.log(nsdisval)
        result += nstest

    if USE3JLOST:
        if mbpjet > Mass_3J_Lost_right_binning[-1]:
            bin_to_use = int(Mass_3J_Lost_right_binning[-1] - 1)
        elif mbpjet < Mass_3J_Lost_right_binning[0]:
            bin_to_use = 0
        else:
            bin_to_use = np.argmax(Mass_3J_Lost_right_binning > mbpjet) - 1 ## -1 since binning is edges, not center

        massdisval = Mass_3J_Lost_right_values[bin_to_use]/np.sum(Mass_3J_Lost_right_values) # normalizes dist based on integral
        if massdisval > 1.0E-10:
            masstest = -1*np.log(massdisval)
        result += masstest

    return nstest, masstest, result


@njit()
def solve_3J_merged(maxmjet, mbpjet, nschi):
    'Inputs: max(jet mass), m(b+jet), and the NS solver chi2 values (maxmjet, mbpjet, nschi)\nOutputs: tuple of the NS, mass, and combined disriminants (NuDiscr, MassDiscr, Prob)'
        ## initialize values used to build up likelihood
    result = 0.
    nstest = np.inf
    masstest = np.inf

        ## evaluate probabilities
    if USENS:
            ## find binning to use for when sqrt(nschi) is in under/overflow
        if np.sqrt(nschi) > NS_3J_Merged_right_binning[-1]:
            bin_to_use = int(NS_3J_Merged_right_binning[-1] - 1)
        elif np.sqrt(nschi) < NS_3J_Merged_right_binning[0]:
            bin_to_use = 0
        else:
            bin_to_use = np.argmax(NS_3J_Merged_right_binning > np.sqrt(nschi)) - 1 ## -1 since binning is edges, not center

        nsdisval = NS_3J_Merged_right_values[bin_to_use]/np.sum(NS_3J_Merged_right_values) # normalizes dist based on integral
        if nsdisval > 1.0E-10:
            nstest = -1*np.log(nsdisval)
        result += nstest

    if USE3JMERGED:
            ## find x binning to use for when maxmjet is in under/overflow
        if maxmjet > Mass_3J_Merged_right_binning[0][-1]:
            xbin_to_use = int(Mass_3J_Merged_right_binning[0][-1] - 1)
        elif maxmjet < Mass_3J_Merged_right_binning[0][0]:
            xbin_to_use = 0
        else:
            xbin_to_use = np.argmax(Mass_3J_Merged_right_binning[0] > maxmjet) - 1

            ## find y binning to use for when maxmjet is in under/overflow
        if mbpjet > Mass_3J_Merged_right_binning[1][-1]:
            ybin_to_use = int(Mass_3J_Merged_right_binning[1][-1] - 1)
        elif mbpjet < Mass_3J_Merged_right_binning[1][0]:
            ybin_to_use = 0
        else:
            ybin_to_use = np.argmax(Mass_3J_Merged_right_binning[1] > mbpjet) - 1
        
        massdisval = Mass_3J_Merged_right_values[xbin_to_use][ybin_to_use]/np.sum(Mass_3J_Merged_right_values) # normalizes dist based on integral
        if massdisval > 1.0E-10:
            masstest = -1*np.log(massdisval)
        result += masstest

    return nstest, masstest, result

