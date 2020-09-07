import os
from pdb import set_trace
import Utilities.prettyjson as prettyjson

from argparse import ArgumentParser
parser = ArgumentParser()
#parser.add_argument('sample', default='data', nargs='?', choices=['MC', 'data', 'all'], help='specify which samples to print')
args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']

years = ['2016', '2017', '2018']
for year in years:
    samples_list = prettyjson.loads(open('%s/inputs/samples_%s.json' % (proj_dir, year)).read())
    
    samples_dict = {}
    for sample in samples_list:
        name = sample.pop('name')
        samples_dict[name] = sample

    lumiMasks = {}    
    data_output = "Dataset & Integrated luminosity ($\\fbinv$)\\\\ \n\hline \n"
    
    for sample in samples_dict.keys():
        if not 'data' in sample: continue
        dbs_name = samples_dict[sample]['DBSName'] # get full DBS name
        dbs_name = '/'.join(dbs_name.split('/')[:-1]) # remove 'NANOAOD' part of dbs name
        data_output += "{SAMPLE} & {LUMI:.1f} \\\\ \n".format(SAMPLE=dbs_name, LUMI=0.0)

            # add lumimask used for datset
        lumiMasks[sample] = samples_dict[sample]['lumimask']
    
        # add lumimask to output
    if len(sorted(set(lumiMasks.values()))) == 1:
        data_output += "\n\n%s\n\n" % sorted(set(lumiMasks.values()))[0]

    print(data_output)
    data_datasets_lumi = open('%s/AN_Scripts/%s_data_datasets.txt' % (proj_dir, year), 'w')
    data_datasets_lumi.write(data_output)
    data_datasets_lumi.close()
