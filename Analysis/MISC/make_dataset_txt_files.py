#! /bin/env python

"""
This script searches DBS for the datasets from the input json file and dumps the associated rootfiles into txt files

Created by Joseph Dulemba
30 October 2020
"""
from pdb import set_trace
import os, sys
import subprocess
from fnmatch import fnmatch
from pdb import set_trace
import Utilities.prettyjson as prettyjson
import Utilities.das as das

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('json', help='json file containing the samples definition')
parser.add_argument('--sample', help='Sample to run on, POSIX regex allowed')
parser.add_argument('--test', action='store_true', help="Output txt file named 'tmp.txt'")
parser.add_argument('--options', help='command-line arguments'
                    ' to be passed to the configuration', default="")

args = parser.parse_args()

jobid = os.environ['jobid']
proj_dir = os.environ['PROJECT_DIR']

if not os.path.isfile(args.json):
   raise ValueError('file %s does not exist' % args.json)

outdir = os.path.join(proj_dir, 'inputs', '_'.join(os.path.basename(args.json).split('.')[0].split('_')[1:])) # get name of json file except for 'samples_'
if not os.path.isdir(outdir): os.path.makedirs(outdir)
#set_trace()

all_samples = prettyjson.loads(open(args.json).read())
samples_to_run = list(filter(
   lambda x: fnmatch(x['name'], args.sample if args.sample else '*'),
     all_samples
))
if not len(samples_to_run):
    raise RuntimeError('Could not find any sample matching the pattern')

for sample in samples_to_run:
    #set_trace()
    
    if 'DBSName' in sample:
        txtname = os.path.join(outdir, '%s_tmp.txt' % sample['name']) if args.test else os.path.join(outdir, '%s.txt' % sample['name'])

        flist = das.query('file dataset=%s' % sample['DBSName'])#, True)
        flist = [fname.replace('/store', 'root://cmsxrootd.fnal.gov//store') for fname in flist] # add xrootd redirector
        fnames = '\n'.join(flist)

        txt_out = open(txtname, 'w')
        txt_out.write(fnames)
        txt_out.close()
        print("%s written" % txtname)
    else:
        raise ValueError("No DBSName found for sample: %s" % sample['name'])
    
