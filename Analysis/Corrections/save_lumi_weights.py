from coffea.util import load, save
from pdb import set_trace
import os
from fnmatch import fnmatch
import Utilities.prettyjson as prettyjson

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']

outdir = os.path.join(proj_dir, 'Corrections', jobid)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

data_lumi = prettyjson.loads(open(os.path.join(proj_dir,'inputs', 'lumis_data.json')).read()) # file with integrated luminosity for all three years

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

# for each year, read sumGenWeights from all meta.json files
for year in ['2016', '2017', '2018']:
    print(year)
    #set_trace()
    xsec_file = prettyjson.loads(open(os.path.join(proj_dir, 'inputs', 'samples_%s.json' % year)).read()) # file with cross sections
    #xsec_file = prettyjson.loads(open(os.path.join(proj_dir, 'inputs', 'samples_%s_%s.json' % (year, jobid))).read()) # file with cross sections
    datasets = list(filter(lambda x: fnmatch(x['name'], '*'), xsec_file))
    for dataset in datasets:
        sample = dataset['name']
        if sample.startswith('data_Single'): continue
        if dataset['DBSName'] == 'NOT PRESENT':
            print("Dataset %s not present, will be skipped" % sample)
            continue
        if not os.path.isfile(os.path.join(proj_dir, 'inputs', '%s_%s' % (year, jobid), '%s.meta.json' % sample)):
            print("No meta.json file found for dataset %s, skipping" % sample)
            continue
        print('    %s' % sample)
        meta_json = prettyjson.loads(open(os.path.join(proj_dir, 'inputs', '%s_%s' % (year, jobid), '%s.meta.json' % sample)).read())
        sumGenWeights = meta_json["sumGenWeights"]
        xsec = dataset['xsection']
        for lep in ['Electrons', 'Muons']:
            lumi_weights[year][lep][sample] = data_lumi[year][lep]/(sumGenWeights/xsec)

    print("%s calculated" % year)

#set_trace()
    # save files
mcweights_name = os.path.join(outdir, 'MC_LumiWeights.coffea')
save(lumi_weights, mcweights_name)
print('\n', mcweights_name, 'written')
