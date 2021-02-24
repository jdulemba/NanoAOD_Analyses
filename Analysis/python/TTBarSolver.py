'''
This file computes the probabilities from mass and neutrino solver distributions (and their sums)
for different jet multiplicities (4+ jets, 3 jets lost, 3 jets merged)
'''

import Utilities.prettyjson as prettyjson
from pdb import set_trace
import os
import numpy as np
from numba import njit
from coffea.util import load

class TTSolver(object):

    def __init__(self, year):
        print('TTBarSolver:', year)
        proj_dir = os.environ['PROJECT_DIR']
        jobid = os.environ['jobid']
        base_jobid = os.environ['base_jobid']
        cfg_pars = prettyjson.loads(open(os.path.join(proj_dir, 'cfg_files', 'cfg_pars_%s.json' % jobid)).read())['ttsolver']
        
        probs = load(os.path.join(proj_dir, 'Corrections', base_jobid, cfg_pars['filename']))[year]

            ## create arrays for binning and values separately for each dist because njit can't handle constant dictionaries currently
        self.USEMASS = cfg_pars['USEMASS']
        WTmass_right = probs['4PJets']['mWHad_vs_mTHad']
        self.WTmass_right_binning = WTmass_right._axes
        self.WTmass_right_values  = WTmass_right._values
        
        self.USE3JMERGED = cfg_pars['USE3JMERGED']
        Mass_3J_Merged_right = probs['3Jets']['Merged_mTHadProxy_vs_maxmjet']
        self.Mass_3J_Merged_right_binning = Mass_3J_Merged_right._axes
        self.Mass_3J_Merged_right_values  = Mass_3J_Merged_right._values
        
        self.USE3JLOST = cfg_pars['USE3JLOST']
        Mass_3J_Lost_right = probs['3Jets']['Lost_mTHadProxy']
        self.Mass_3J_Lost_right_binning = Mass_3J_Lost_right._axes
        self.Mass_3J_Lost_right_values  = Mass_3J_Lost_right._values
        
        self.USENS = cfg_pars['USENS']
        NS_4PJ_right = probs['4PJets']['nusolver_dist']
        #NS_4PJ_right = probs['4PJets']['nusolver_chi2']
        self.NS_4PJ_right_binning = NS_4PJ_right._axes
        self.NS_4PJ_right_values  = NS_4PJ_right._values
        
        # merged 3 jet vars
        NS_3J_Merged_right = probs['3Jets']['Merged_nusolver_dist']
        #NS_3J_Merged_right = probs['3Jets']['Merged_nusolver_chi2']
        self.NS_3J_Merged_right_binning = NS_3J_Merged_right._axes
        self.NS_3J_Merged_right_values  = NS_3J_Merged_right._values
        
        # lost 3 jet vars
        NS_3J_Lost_right = probs['3Jets']['Lost_nusolver_dist']
        #NS_3J_Lost_right = probs['3Jets']['Lost_nusolver_chi2']
        self.NS_3J_Lost_right_binning = NS_3J_Lost_right._axes
        self.NS_3J_Lost_right_values  = NS_3J_Lost_right._values

    
    def solve_4PJ(self, mthad, mwhad, nschi):
        'Inputs: m(thad), m(whad), and the NS solver chi2 values (mthad, mwhad, nschi)\nOutputs: numpy array of the combined, mass, and NS disriminants (Prob, MassDiscr, NuDiscr)'
        ns_pars = (self.USENS, self.NS_4PJ_right_binning, self.NS_4PJ_right_values)
        mass_pars = (self.USEMASS, self.WTmass_right_binning, self.WTmass_right_values)
        return self._solve_4PJ(mthad, mwhad, nschi, ns_pars, mass_pars)
 
    @staticmethod
    @njit()
    def _solve_4PJ(mthad, mwhad, nschi, ns_pars, mass_pars):
        'Inputs: m(thad), m(whad), and the NS solver chi2 values (mthad, mwhad, nschi)\nOutputs: numpy array of the combined, mass, and NS disriminants (Prob, MassDiscr, NuDiscr)'
        ## initialize values used to build up likelihood
        result = 0.
        nstest = np.inf
        masstest = np.inf
    
            ## evaluate probabilities
        if ns_pars[0]:
            if nschi < 1e12:
                nsbinning = ns_pars[1] 
                nsvalues  = ns_pars[2] 
                    ## find binning to use for when sqrt(nschi) is in under/overflow
                if np.sqrt(nschi) > nsbinning[-1]:
                    bin_to_use = int(nsbinning[-1] - 1)
                elif np.sqrt(nschi) < nsbinning[0]:
                    bin_to_use = 0
                else:
                    bin_to_use = np.argmax(nsbinning > np.sqrt(nschi)) - 1 ## -1 since binning is edges, not center
    
                nsdisval = nsvalues[bin_to_use]/np.sum(nsvalues) # normalizes dist based on integral
                if nsdisval > 1.0E-10:
                    nstest = -1*np.log(nsdisval)
            result += nstest
    
        if mass_pars[0]:
            mbinning = mass_pars[1]
            mvalues  = mass_pars[2]
                ## find x binning to use for when mwhad is in under/overflow
            if mthad > mbinning[0][-1]:
                xbin_to_use = int(mbinning[0][-1] - 1)
            elif mthad < mbinning[0][0]:
                xbin_to_use = 0
            else:
                xbin_to_use = np.argmax(mbinning[0] > mthad) - 1
    
                ## find y binning to use for when mwhad is in under/overflow
            if mwhad > mbinning[1][-1]:
                ybin_to_use = int(mbinning[1][-1] - 1)
            elif mwhad < mbinning[1][0]:
                ybin_to_use = 0
            else:
                ybin_to_use = np.argmax(mbinning[1] > mwhad) - 1
            
            massdisval = mvalues[xbin_to_use][ybin_to_use]/np.sum(mvalues) # normalizes dist based on integral
            if massdisval > 1.0E-10:
                masstest = -1*np.log(massdisval)
            result += masstest
    
        return np.array([result, masstest, nstest])

    def solve_3J_lost(self, mbpjet, nschi):
        'Inputs: m(b+jet) and the NS solver chi2 values (mbpjet, nschi)\nOutputs: numpy array of the combined, mass, and NS disriminants (Prob, MassDiscr, NuDiscr)'
        ns_pars = (self.USENS, self.NS_3J_Lost_right_binning, self.NS_3J_Lost_right_values)
        mass_pars = (self.USE3JLOST, self.Mass_3J_Lost_right_binning, self.Mass_3J_Lost_right_values)
        return self._solve_3J_lost(mbpjet, nschi, ns_pars, mass_pars)
        
    
    @staticmethod
    @njit()
    def _solve_3J_lost(mbpjet, nschi, ns_pars, mass_pars):
        'Inputs: m(b+jet) and the NS solver chi2 values (mbpjet, nschi)\nOutputs: numpy array of the combined, mass, and NS disriminants (Prob, MassDiscr, NuDiscr)'
            ## initialize values used to build up likelihood
        result = 0.
        nstest = np.inf
        masstest = np.inf
    
            ## evaluate probabilities
        if ns_pars[0]:
            if nschi < 1e12:
                nsbinning = ns_pars[1] 
                nsvalues  = ns_pars[2] 
                    ## find binning to use for when sqrt(nschi) is in under/overflow
                if np.sqrt(nschi) > nsbinning[-1]:
                    bin_to_use = int(nsbinning[-1] - 1)
                elif np.sqrt(nschi) < nsbinning[0]:
                    bin_to_use = 0
                else:
                    bin_to_use = np.argmax(nsbinning > np.sqrt(nschi)) - 1 ## -1 since binning is edges, not center
    
                nsdisval = nsvalues[bin_to_use]/np.sum(nsvalues) # normalizes dist based on integral
                if nsdisval > 1.0E-10:
                    nstest = -1*np.log(nsdisval)
            result += nstest
    
        if mass_pars[0]:
            mbinning = mass_pars[1] 
            mvalues  = mass_pars[2] 
            if mbpjet > mbinning[-1]:
                bin_to_use = int(mbinning[-1] - 1)
            elif mbpjet < mbinning[0]:
                bin_to_use = 0
            else:
                bin_to_use = np.argmax(mbinning > mbpjet) - 1 ## -1 since binning is edges, not center
    
            massdisval = mvalues[bin_to_use]/np.sum(mvalues) # normalizes dist based on integral
            if massdisval > 1.0E-10:
                masstest = -1*np.log(massdisval)
            result += masstest
    
        return np.array([result, masstest, nstest])
    
    def solve_3J_merged(self, maxmjet, mbpjet, nschi):
        'Inputs: max(jet mass), m(b+jet), and the NS solver chi2 values (maxmjet, mbpjet, nschi)\nOutputs: numpy array of the combined, mass, and NS disriminants (Prob, MassDiscr, NuDiscr)'
        ns_pars = (self.USENS, self.NS_3J_Merged_right_binning, self.NS_3J_Merged_right_values)
        mass_pars = (self.USE3JLOST, self.Mass_3J_Merged_right_binning, self.Mass_3J_Merged_right_values)
        return self._solve_3J_merged(maxmjet, mbpjet, nschi, ns_pars, mass_pars)
       
 
    @staticmethod
    @njit()
    def _solve_3J_merged(maxmjet, mbpjet, nschi, ns_pars, mass_pars):
        'Inputs: max(jet mass), m(b+jet), and the NS solver chi2 values (maxmjet, mbpjet, nschi)\nOutputs: numpy array of the combined, mass, and NS disriminants (Prob, MassDiscr, NuDiscr)'
            ## initialize values used to build up likelihood
        result = 0.
        nstest = np.inf
        masstest = np.inf
    
            ## evaluate probabilities
        if ns_pars[0]:
            if nschi < 1e12:
                nsbinning = ns_pars[1] 
                nsvalues  = ns_pars[2] 
                    ## find binning to use for when sqrt(nschi) is in under/overflow
                if np.sqrt(nschi) > nsbinning[-1]:
                    bin_to_use = int(nsbinning[-1] - 1)
                elif np.sqrt(nschi) < nsbinning[0]:
                    bin_to_use = 0
                else:
                    bin_to_use = np.argmax(nsbinning > np.sqrt(nschi)) - 1 ## -1 since binning is edges, not center
    
                nsdisval = nsvalues[bin_to_use]/np.sum(nsvalues) # normalizes dist based on integral
                if nsdisval > 1.0E-10:
                    nstest = -1*np.log(nsdisval)
            result += nstest
    
        if mass_pars[0]:
            mbinning = mass_pars[1] 
            mvalues  = mass_pars[2] 
                ## find x binning to use for when maxmjet is in under/overflow
            if maxmjet > mbinning[0][-1]:
                xbin_to_use = int(mbinning[0][-1] - 1)
            elif maxmjet < mbinning[0][0]:
                xbin_to_use = 0
            else:
                xbin_to_use = np.argmax(mbinning[0] > maxmjet) - 1
    
                ## find y binning to use for when maxmjet is in under/overflow
            if mbpjet > mbinning[1][-1]:
                ybin_to_use = int(mbinning[1][-1] - 1)
            elif mbpjet < mbinning[1][0]:
                ybin_to_use = 0
            else:
                ybin_to_use = np.argmax(mbinning[1] > mbpjet) - 1
            
            massdisval = mvalues[xbin_to_use][ybin_to_use]/np.sum(mvalues) # normalizes dist based on integral
            if massdisval > 1.0E-10:
                masstest = -1*np.log(massdisval)
            result += masstest
    
        return np.array([result, masstest, nstest])

