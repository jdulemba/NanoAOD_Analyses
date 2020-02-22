from pdb import set_trace

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

