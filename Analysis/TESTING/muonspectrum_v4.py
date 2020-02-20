from coffea import hist
from coffea.analysis_objects import JaggedCandidateArray
import coffea.processor as processor
from awkward import JaggedArray
import numpy as np
from pdb import set_trace

import matplotlib.pyplot as plt
plt.switch_backend('agg')

# Look at ProcessorABC to see the expected methods and what they are supposed to do
class FancyDimuonProcessor(processor.ProcessorABC):
    def __init__(self):
        dataset_axis = hist.Cat("dataset", "Primary dataset")
        mass_axis = hist.Bin("mass", r"$m_{\mu\mu}$ [GeV]", 600, 0.25, 300)
        pt_axis = hist.Bin("pt", r"$p_{T,\mu}$ [GeV]", 3000, 0.25, 300)
        
        self._accumulator = processor.dict_accumulator({
            'mass': hist.Hist("Counts", dataset_axis, mass_axis),
            'mass_near': hist.Hist("Counts", dataset_axis, mass_axis),
            'mass_far': hist.Hist("Counts", dataset_axis, mass_axis),
            'pt_lead': hist.Hist("Counts", dataset_axis, pt_axis),
            'pt_trail': hist.Hist("Counts", dataset_axis, pt_axis),
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
            pt=df['Muon_pt'],
            eta=df['Muon_eta'],
            phi=df['Muon_phi'],
            mass=df['Muon_mass'],
            charge=df['Muon_charge'],
            softId=df['Muon_softId'],
            tightId=df['Muon_tightId']
            )        
        
        output['cutflow']['all events'] += muons.size
        
        soft_id = (muons.softId > 0)
        muons = muons[soft_id]
        output['cutflow']['soft id'] += soft_id.any().sum()
        
        twomuons = (muons.counts >= 2)
        output['cutflow']['two muons'] += twomuons.sum()
        
        dimuons = muons[twomuons].distincts()
        
        twodimuons = (dimuons.counts >= 2)
        output['cutflow']['>= two dimuons'] += twodimuons.sum()
        dimuons = dimuons[twodimuons]
        
        opposite_charge = (dimuons.i0['charge'] * dimuons.i1['charge'] == -1)
        
        dimuons = dimuons[opposite_charge]
        output['cutflow']['opposite charge'] += opposite_charge.any().sum()
        
        mass_20GeV = (dimuons.mass > 35)
        dimuons = dimuons[mass_20GeV]
        
        exactlytwodimuons = (dimuons.counts == 2)
        output['cutflow']['== two dimuons'] += exactlytwodimuons.sum()
        dimuons = dimuons[exactlytwodimuons].compact()
        
        leading_mu = (dimuons.i0.pt.content > dimuons.i1.pt.content)
        pt_lead = JaggedArray.fromoffsets(dimuons.offsets, np.where(leading_mu, 
                                                                    dimuons.i0.pt.content, dimuons.i1.pt.content))
        pt_trail = JaggedArray.fromoffsets(dimuons.offsets, np.where(~leading_mu, 
                                                                     dimuons.i0.pt.content, dimuons.i1.pt.content))
        
        near_z = np.abs(dimuons.mass - 91.118).argmin()
        far_z = np.abs(dimuons.mass - 91.118).argmax()
        
        output['mass'].fill(dataset=dataset,
                            mass=dimuons.p4.sum().mass)
        output['mass_near'].fill(dataset=dataset, 
                                 mass=dimuons.mass[near_z].flatten())
        output['mass_far'].fill(dataset=dataset, 
                                mass=dimuons.mass[far_z].flatten())
        output['pt_lead'].fill(dataset=dataset,
                               pt=pt_lead.flatten())
        output['pt_trail'].fill(dataset=dataset,
                                pt=pt_trail.flatten())
        return output

    def postprocess(self, accumulator):
        return accumulator


fileset = {
    'DoubleMuon': [
        'data/Run2012B_DoubleMuParked.root',
        'data/Run2012C_DoubleMuParked.root',
    ],
    'ZZ to 4mu': [
        'data/ZZTo4mu.root'
    ]
}

output = processor.run_uproot_job(fileset,
                                  treename='Events',
                                  processor_instance=FancyDimuonProcessor(),
                                  executor=processor.futures_executor,
                                  executor_args={'workers': 6, 'flatten': True},
                                  chunksize=500000,
         )

print(output)

#set_trace()
fig = plt.figure()
ax = hist.plot1d(output['mass'], overlay='dataset')
ax.set_xlim(70,150)
ax.set_ylim(0, 3000)
fig.savefig('muonspectrum_v4_mass')
plt.close()

fig = plt.figure()
ax = hist.plot1d(output['mass_near'], overlay='dataset')
#ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_xlim(60,120)
ax.set_ylim(0.1, 7500)
fig.savefig('muonspectrum_v4_mass_near')
plt.close()

fig = plt.figure()
ax = hist.plot1d(output['mass_far'], overlay='dataset')
#ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_ylim(0.1, 8000)
fig.savefig('muonspectrum_v4_mass_far')
plt.close()

fig = plt.figure()
ax = hist.plot1d(output['pt_lead'], overlay='dataset')
#ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(0.1, 5e3)
fig.savefig('muonspectrum_v4_pt_lead')
plt.close()

fig = plt.figure()
ax = hist.plot1d(output['pt_trail'], overlay='dataset')
#ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim(0.1, 2e4)
fig.savefig('muonspectrum_v4_pt_trail')
plt.close()


