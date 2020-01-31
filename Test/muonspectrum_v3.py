import matplotlib.pyplot as plt
plt.switch_backend('agg')
from coffea import hist
from coffea.analysis_objects import JaggedCandidateArray
import coffea.processor as processor

# Look at ProcessorABC documentation to see the expected methods and what they are supposed to do
class DimuonProcessor(processor.ProcessorABC):
    def __init__(self):
        dataset_axis = hist.Cat("dataset", "Primary dataset")
        mass_axis = hist.Bin("mass", r"$m_{\mu\mu}$ [GeV]", 30000, 0.25, 300)
        
        self._accumulator = processor.dict_accumulator({
            'mass': hist.Hist("Counts", dataset_axis, mass_axis),
            'cutflow': processor.defaultdict_accumulator(int),
        })
    
    @property
    def accumulator(self):
        return self._accumulator
    
    def process(self, df):
        output = self.accumulator.identity()
        
        dataset = df['dataset']
        muons = JaggedCandidateArray.candidatesfromcounts(
            df['nMuon'],
            pt=df['Muon_pt'].content,
            eta=df['Muon_eta'].content,
            phi=df['Muon_phi'].content,
            mass=df['Muon_mass'].content,
            charge=df['Muon_charge'].content,
            )
        
        output['cutflow']['all events'] += muons.size
        
        twomuons = (muons.counts == 2)
        output['cutflow']['two muons'] += twomuons.sum()
        
        opposite_charge = twomuons & (muons['charge'].prod() == -1)
        output['cutflow']['opposite charge'] += opposite_charge.sum()
        
        dimuons = muons[opposite_charge].distincts()
        output['mass'].fill(dataset=dataset, mass=dimuons.mass.flatten())
        
        return output

    def postprocess(self, accumulator):
        return accumulator

fileset = {
    'DoubleMuon': [
        'data/Run2012B_DoubleMuParked.root',
        'data/Run2012C_DoubleMuParked.root',
    ]
}

output = processor.run_uproot_job(fileset,
                                  treename='Events',
                                  processor_instance=DimuonProcessor(),
                                  executor=processor.futures_executor,
                                  executor_args={'workers': 4},
        )

print(output)

fig = plt.figure()
ax = hist.plot1d(output['mass'], overlay='dataset')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(0.1, 1e6)
fig.savefig('muonspectrum_v3')
