from coffea.util import load, save
from pdb import set_trace
import os
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
    if year == '2016':
        Nominal_ttJets = ['ttJets_PS']#, 'ttJets']
    else:
        Nominal_ttJets = ['ttJetsSL', 'ttJetsHad', 'ttJetsDiLep']
    xsec_file = prettyjson.loads(open(os.path.join(proj_dir, 'inputs', 'samples_%s.json' % year)).read()) # file with cross sections
    samples = sorted([fname.split('.')[0] for fname in os.listdir(os.path.join(proj_dir, 'inputs', '%s_%s' % (year, jobid))) if fname.endswith('.meta.json')])
    for sample in samples:
        print('    %s' % sample)
        if sample in Nominal_ttJets:
            #set_trace()
            meta_json = prettyjson.loads(open(os.path.join(proj_dir, 'inputs', '%s_%s' % (year, jobid), '%s.meta.json' % sample)).read())
            sumGenWeights_nominal = sum([meta_json["sumGenWeights_%i" % idx] for idx in range(3, 10)]) # get evtIdx 3-9
            sumGenWeights_signal = sum([meta_json["sumGenWeights_%i" % idx] for idx in [1, 2]])
        else:
            sumGenWeights = prettyjson.loads(open(os.path.join(proj_dir, 'inputs', '%s_%s' % (year, jobid), '%s.meta.json' % sample)).read())["sumGenWeights"]
        xsec = [info['xsection'] for info in xsec_file if info['name'] == sample ][0]
        for lep in ['Electrons', 'Muons']:
            if sample in Nominal_ttJets:
                lumi_weights[year][lep][sample] = data_lumi[year][lep]/(sumGenWeights_nominal/xsec)
                lumi_weights[year][lep]['signal'] = data_lumi[year][lep]/(sumGenWeights_signal/xsec)
            else:
                lumi_weights[year][lep][sample] = data_lumi[year][lep]/(sumGenWeights/xsec)

    print("%s calculated" % year)

#set_trace()
    # save files
mcweights_name = os.path.join(outdir, 'MC_LumiWeights_IgnoreSigEvts.coffea')
save(lumi_weights, mcweights_name)
print('\n', mcweights_name, 'written')
