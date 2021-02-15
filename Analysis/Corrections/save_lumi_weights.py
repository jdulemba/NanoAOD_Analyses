from coffea.util import load, save
from pdb import set_trace
import os
from fnmatch import fnmatch
import Utilities.prettyjson as prettyjson

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
base_jobid = os.environ['base_jobid']

outdir = os.path.join(proj_dir, 'Corrections', base_jobid)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

data_lumi = prettyjson.loads(open(os.path.join(proj_dir,'inputs', '%s_lumis_data.json' % base_jobid)).read()) # file with integrated luminosity for all three years

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

years_to_run = ['2017', '2018'] if base_jobid == 'ULnanoAOD' else ['2016', '2017', '2018']
# for each year, read sumGenWeights from all meta.json files
for year in years_to_run:
    print(year)
    xsec_file = prettyjson.loads(open(os.path.join(proj_dir, 'inputs', 'samples_%s_%s.json' % (year, base_jobid))).read()) # file with cross sections
    datasets = list(filter(lambda x: fnmatch(x['name'], '*'), xsec_file))
    for dataset in datasets:
        sample = dataset['name']
        if sample.startswith('data_Single'): continue
        if dataset['DBSName'] == 'NOT PRESENT':
            print("Dataset %s not present, will be skipped" % sample)
            continue
        if not os.path.isfile(os.path.join(proj_dir, 'inputs', '%s_%s' % (year, base_jobid), '%s.meta.json' % sample)):
            print("No meta.json file found for dataset %s, skipping" % sample)
            continue
        print('    %s' % sample)
        meta_json = prettyjson.loads(open(os.path.join(proj_dir, 'inputs', '%s_%s' % (year, base_jobid), '%s.meta.json' % sample)).read())
        sumGenWeights = meta_json["sumGenWeights"]
        xsec = dataset['xsection']
        for lep in ['Electrons', 'Muons']:
            lumi_weights[year][lep][sample] = data_lumi[year][lep]/(sumGenWeights/xsec)

    print("%s calculated" % year)

    # save files
mcweights_name = os.path.join(outdir, 'MC_LumiWeights.coffea')
save(lumi_weights, mcweights_name)
print('\n', mcweights_name, 'written')
