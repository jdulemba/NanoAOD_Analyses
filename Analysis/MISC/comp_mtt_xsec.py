#! /bin/env python

"""
This script calculates the ttbar cross section for different values of mtop.
Based on 13TeV values from NNPDF2.3 central values from this twiki: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/TtbarNNLO
"""

from pdb import set_trace

sigma_mref = 843.483
mref = 172.5
a1 = -0.745047
a2 = 0.127417


brs = {"ll" : 0.105, "jj" : 0.457, "lj" : 0.438}
def sigma_mtt(mt):
    sigma_mtt = sigma_mref*(mref/mt)**4*(1+a1*((mt-mref)/mref) + a2*((mt-mref)/mref)**2)
    return sigma_mtt

for mt in [166.5, 169.5, 171.5, 173.5, 175.5, 178.5]:
    mtt_xsec = sigma_mtt(mt)

    #set_trace()
    print(f"Total xsec for {mt} = {mtt_xsec}")
    for decay, br in brs.items():
        print("\t%s %f" % (decay, mtt_xsec*br))
    print("\n")
