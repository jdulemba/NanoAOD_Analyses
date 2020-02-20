import uproot
import awkward
from coffea import hist
from coffea.analysis_objects import JaggedCandidateArray
from coffea.processor import defaultdict_accumulator
from pdb import set_trace

import matplotlib.pyplot as plt
plt.switch_backend('agg')

files = [
    'data/Run2012B_DoubleMuParked.root',
    'data/Run2012C_DoubleMuParked.root',
]

masshist = hist.Hist("Counts", hist.Bin("mass", r"$m_{\mu\mu}$ [GeV]", 30000, 0.25, 300))
cutflow = defaultdict_accumulator(lambda: 0)

branches = ['nMuon', 'Muon_pt', 'Muon_eta', 'Muon_phi', 'Muon_mass', 'Muon_charge']
for chunk in uproot.iterate(files, 'Events', branches=branches, entrysteps=500000, namedecode='ascii'):
    muons = JaggedCandidateArray.candidatesfromcounts(chunk['nMuon'],
                                            pt=chunk['Muon_pt'].content,
                                            eta=chunk['Muon_eta'].content,
                                            phi=chunk['Muon_phi'].content,
                                            mass=chunk['Muon_mass'].content,
                                            charge=chunk['Muon_charge'].content,
                                           )
    
    cutflow['all events'] += muons.size
    twomuons = (muons.counts == 2)
    cutflow['two muons'] += twomuons.sum()
    opposite_charge = twomuons & (muons['charge'].prod() == -1)
    cutflow['opposite charge'] += opposite_charge.sum()
    dimuons = muons[opposite_charge].distincts()
    masshist.fill(mass=dimuons.mass.flatten())

print(dict(cutflow))

fig = plt.figure()
ax = hist.plot1d(masshist)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(0.1, 1e6)
fig.savefig('muonspectrum_v2')
