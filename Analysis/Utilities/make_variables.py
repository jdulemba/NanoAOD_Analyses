from pdb import set_trace
import numpy as np
import awkward as ak

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
        # broadcast met into same shape as leptons
    met_pt = (ak.ones_like(leptons.pt))*(met.pt)
    met_px = (ak.ones_like(leptons.pt))*(met.px)
    met_py = (ak.ones_like(leptons.pt))*(met.py)

    return np.sqrt( np.square(leptons.pt + met_pt) - np.square(leptons.px + met_px) - np.square(leptons.py + met_py) )

def ctstar_flat(top_p4, tbar_p4):
    # convert 4vecs to cartesian
    tops = top_p4._to_cartesian()
    tbars = tbar_p4._to_cartesian()

    ttbar = (tops+tbars)

    # find top quarks in rest frame
    tt_boost = ttbar.boostp3
    rf_tops = tops.boost(tt_boost*-1)
    rf_tbars = tbars.boost(tt_boost*-1)

    # find ctstar
    top_ctstar = ttbar.p3.unit.dot(rf_tops.p3.unit)
    tbar_ctstar = ttbar.p3.unit.dot(rf_tbars.p3.unit)

    return top_ctstar, tbar_ctstar
