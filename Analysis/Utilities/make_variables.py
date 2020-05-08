from pdb import set_trace
#from numba import njit, objmode
import numpy as np
#import awkward

def ctstar(top_p4, tbar_p4, debug=False):

    if debug:
        ## for testing
        tops = top_p4
        tbars = tbar_p4
    else:
        # convert 4vecs to cartesian
        tops = top_p4._to_cartesian().flatten()
        tbars = tbar_p4._to_cartesian().flatten()

    ttbar = (tops+tbars)

    # find top quarks in rest frame
    tt_boost = ttbar.boostp3
    rf_tops = tops.boost(tt_boost*-1)
    rf_tbars = tbars.boost(tt_boost*-1)

    # find ctstar
    top_ctstar = ttbar.p3.unit.dot(rf_tops.p3.unit)
    tbar_ctstar = ttbar.p3.unit.dot(rf_tbars.p3.unit)

    return top_ctstar, tbar_ctstar

def MT(leptons, met, debug=False):
    if debug: set_trace()

    met_pt = (leptons.pt.ones_like())*(met.pt.flatten())
    met_px = (leptons.pt.ones_like())*(met.p4.x.flatten())
    met_py = (leptons.pt.ones_like())*(met.p4.y.flatten())

    return np.sqrt( np.square(leptons.pt + met_pt) - np.square(leptons.p4.x + met_px) - np.square(leptons.p4.y + met_py) )
