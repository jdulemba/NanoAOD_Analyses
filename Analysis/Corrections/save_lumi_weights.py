from coffea.util import save
from pdb import set_trace
import os
import Utilities.prettyjson as prettyjson

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']

outdir = '%s/Corrections' % proj_dir
if not os.path.isdir(outdir):
    os.makedirs(outdir)

data_lumi =  prettyjson.loads(open('%s/inputs/data_lumis.json' % proj_dir).read()) # file with integrated luminosity for all three years

lumi_weights = {
    '2016' : {},
    '2017' : {},
    '2018' : {},
}

# for each year, read nWeightedEvts from all meta.json files
for year in ['2016', '2017', '2018']:
    xsec_file = prettyjson.loads(open('%s/inputs/samples_%s.json' % (proj_dir, year)).read()) # file with cross sections
    samples = sorted([fname.split('.')[0] for fname in os.listdir('%s/inputs/%s_%s/' % (proj_dir, year, jobid)) if fname.endswith('.meta.json')])
    for sample in samples:
        nWeightedEvts = prettyjson.loads(open('%s/inputs/%s_%s/%s.meta.json' % (proj_dir, year, jobid, sample)).read())["nWeightedEvts"]
        xsec = [info['xsection'] for info in xsec_file if info['name'] == sample ][0]
        lumi_weights[year][sample] = data_lumi[year]/(nWeightedEvts/xsec)

    # save files
mcweights_name = '%s/MC_LumiWeights.coffea' % outdir
save(lumi_weights, mcweights_name)
print('\n', mcweights_name, 'written')
