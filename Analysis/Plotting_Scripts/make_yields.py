from pdb import set_trace
import os
import Utilities.prettyjson as prettyjson
from collections import Counter

from argparse import ArgumentParser
parser = ArgumentParser()

parser.add_argument('analyzer', help='Which analyzer yields to use')
parser.add_argument('year', choices=['2016', '2017', '2018'], help='What year is the ntuple from.')

args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']

leptons = ['Muon', 'Electron']
jmults = ['3Jets', '4PJets']

input_dir = '/'.join([proj_dir, 'plots', '%s_%s' % (args.year, jobid), args.analyzer])

yields_dict = {
    'Electron' : {
        '3Jets' : {},
        '4PJets' : {},
    },
    'Muon' : {
        '3Jets' : {},
        '4PJets' : {},
    },
    '3Jets' : {},
    '4PJets' : {},
}

for jmult in jmults:
    for lep in leptons:
        dtc = '/'.join([input_dir, lep, jmult]) # dir to check
        json_fname = ['%s/%s' % (dtc, fname) for fname in os.listdir(dtc) if fname.endswith('.json')][0]
        if not os.path.isfile(json_fname):
            raise IOError("File %s does not exist" % json_fname)

        yields_dict[lep][jmult] = prettyjson.loads(open(json_fname).read())

for jmult in jmults:
    sum_yields = Counter(yields_dict['Electron'][jmult])+Counter(yields_dict['Muon'][jmult])
    sum_yields['data/SIM'] = round(sum_yields['data']/sum_yields['SIM'], 3)
    yields_dict[jmult] = sum_yields 

with open('%s/yields_compilation.json' % input_dir, 'w') as out:
    out.write(prettyjson.dumps(yields_dict))

