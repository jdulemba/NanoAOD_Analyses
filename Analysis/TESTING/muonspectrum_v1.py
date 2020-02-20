from pdb import set_trace
#import time
import uproot
import awkward
import uproot_methods

import matplotlib.pyplot as plt
plt.switch_backend('agg')
from coffea import hist

#tstart = time.time()

files = [
    'data/Run2012B_DoubleMuParked.root',
    'data/Run2012C_DoubleMuParked.root',
]

masshist = hist.Hist("Counts", hist.Bin("mass", r"$m_{\mu\mu}$ [GeV]", 30000, 0.25, 300))

branches = ['nMuon', 'Muon_pt', 'Muon_eta', 'Muon_phi', 'Muon_mass', 'Muon_charge']
for chunk in uproot.iterate(files, 'Events', branches=branches, entrysteps=500000, namedecode='ascii'):
    p4 = uproot_methods.TLorentzVectorArray.from_ptetaphim(
        chunk.pop('Muon_pt'),
        chunk.pop('Muon_eta'),
        chunk.pop('Muon_phi'),
        chunk.pop('Muon_mass'),
    )
    muons = awkward.JaggedArray.zip(p4=p4, charge=chunk['Muon_charge'])

    twomuons = (muons.counts == 2)
    opposite_charge = (muons['charge'].prod() == -1)
    dimuons = muons[twomuons & opposite_charge].distincts()
    dimuon_mass = (dimuons.i0['p4'] + dimuons.i1['p4']).mass
    masshist.fill(mass=dimuon_mass.flatten())
    
#elapsed = time.time() - tstart

#set_trace()
fig = plt.figure()
ax = hist.plot1d(masshist)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(0.1, 1e6)

fig.savefig('muonspectrum_v1')

#print("Events/s:", masshist.values()[()].sum()/elapsed)
