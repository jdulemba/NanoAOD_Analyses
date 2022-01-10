from pdb import set_trace
import numpy as np
import awkward as ak

def ctstar(top1, top2):
    ttbar = (top1+top2)

    # find top quarks in rest frame
    tt_boost = ttbar.boostvec
    rf_top1, rf_top2 = top1.boost(tt_boost*-1), top2.boost(tt_boost*-1)

    # find ctstar
    top1_ctstar, top2_ctstar = ttbar.unit.dot(rf_top1.unit), ttbar.unit.dot(rf_top2.unit)

    return top1_ctstar, top2_ctstar


def cpTP(top1, top2):
    """
    cpTP is obtained by boosting into the zero momentum frame of the top pair system, and then taking the pz/|p| of the top quark in this frame
    """
    #set_trace()
    ttbar = (top1+top2)

    # find top quarks in rest frame
    tt_boost = ttbar.boostvec
    rf_top1, rf_top2 = top1.boost(tt_boost*-1), top2.boost(tt_boost*-1)

    # find cpTP
    top1_cpTP, top2_cpTP = rf_top1.pz/abs(rf_top1.p), rf_top2.pz/abs(rf_top2.p)

    return top1_cpTP, top2_cpTP


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
