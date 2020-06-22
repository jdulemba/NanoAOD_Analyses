from coffea.util import load, save
from pdb import set_trace
import os
import Utilities.prettyjson as prettyjson

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']

outdir = '%s/Corrections/%s' % (proj_dir, jobid)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

data_lumi =  prettyjson.loads(open('%s/inputs/lumis_data.json' % proj_dir).read()) # file with integrated luminosity for all three years

lumi_weights = {
    '2016' : {
        'Electrons' : {},
        'Muons' :{}
    },
    '2017' : {
        'Electrons' : {},
        'Muons' :{}
    },
    '2018' : {
        'Electrons' : {},
        'Muons' :{}
    },
}

proj_dir = os.environ['PROJECT_DIR']
signal_fname = '%s/signal_scripts/results/signal_effLumi_inds12.coffea' % proj_dir
signalExists = os.path.isfile(signal_fname)
if signalExists:
    signal = load(signal_fname)

# for each year, read nWeightedEvts from all meta.json files
for year in ['2016', '2017', '2018']:
    if year == '2016':
        Nominal_ttJets = ['ttJets_PS', 'ttJets']
    else:
        Nominal_ttJets = ['ttJetsSL', 'ttJetsHad', 'ttJetsDiLep']
    xsec_file = prettyjson.loads(open('%s/inputs/samples_%s.json' % (proj_dir, year)).read()) # file with cross sections
    samples = sorted([fname.split('.')[0] for fname in os.listdir('%s/inputs/%s_%s/' % (proj_dir, year, jobid)) if fname.endswith('.meta.json')])
    for sample in samples:
        if sample in Nominal_ttJets:
            meta_json = prettyjson.loads(open('%s/inputs/%s_%s/%s.meta.json' % (proj_dir, year, jobid, sample)).read())
            nWeightedEvts = sum([meta_json["nWeightedEvts_%i" % idx] for idx in range(3, 10)])
        else:
            nWeightedEvts = prettyjson.loads(open('%s/inputs/%s_%s/%s.meta.json' % (proj_dir, year, jobid, sample)).read())["nWeightedEvts"]
        xsec = [info['xsection'] for info in xsec_file if info['name'] == sample ][0]
        for lep in ['Electrons', 'Muons']:
            lumi_weights[year][lep][sample] = data_lumi[year][lep]/(nWeightedEvts/xsec)

    if signalExists:
        for sig in signal[year].keys():
            effLumi = signal[year][sig]
            for lep in ['Electrons', 'Muons']:
                lumi_weights[year][lep][sig] = data_lumi[year][lep]/effLumi
            

    print("%s calculated" % year)

    # save files
mcweights_name = '%s/MC_LumiWeights_IgnoreSigEvts.coffea' % outdir
save(lumi_weights, mcweights_name)
print('\n', mcweights_name, 'written')
