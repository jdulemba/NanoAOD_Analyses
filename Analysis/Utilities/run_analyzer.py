#!/usr/bin/env python

import os, subprocess, time
from pdb import set_trace
import Run_Jobs.tools as tools
import Utilities.prettyjson as prettyjson
import fnmatch

from argparse import ArgumentParser
parser = ArgumentParser('submit analyzer to the batch queues')
parser.add_argument('analyzer', help='Analyzer to use.')
parser.add_argument('frange', type=str, help='Specify start:stop indices for files, inclusive. 0:1 means valid files in indices 0 and 1 will be used')
parser.add_argument('year', choices=['2016', '2017', '2018'], help='Specify which year to run over')
parser.add_argument('lepton', choices=['Electron', 'Muon'], help='Choose which lepton to select')
parser.add_argument('--sample', type=str, help='Use specific sample')
parser.add_argument('--outfname', type=str, help='Specify output filename, including directory and file extension')
parser.add_argument('--debug', action='store_true', help='Uses iterative_executor for debugging purposes, otherwise futures_excutor will be used (faster)')
args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer=args.analyzer

    ## get samples to use
indir = '/'.join([proj_dir, 'inputs', '%s_%s' % (args.year, jobid)])
samples_to_use = tools.get_sample_list(indir=indir, sample=args.sample) if args.sample else tools.get_sample_list(indir=indir, text_file='analyzer_inputs.txt')
samples_to_use = [sample.split('/')[-1] for sample in samples_to_use if not (sample.split('/')[-1].startswith('data') and args.lepton not in sample.split('/')[-1])] # get rid of data samples that don't correspond to lepton chosen

fileset = {}
for sample in samples_to_use:
    if not os.path.isfile(sample):
        raise IOError("Sample file %s.txt not found" % sample)

    sample_name = sample.split('/')[-1].split('.')[0]
    sfiles = open(sample, 'r')
    valid_inds_fnames = [(idx, fname.strip('\n')) for idx, fname in enumerate(sfiles) if not fname.startswith('#')]
    valid_inds = [x[0] for x in valid_inds_fnames]
    valid_files = [x[1] for x in valid_inds_fnames]

    if ':' in args.frange:
        file_start, file_stop = int((args.frange).split(':')[0]), int((args.frange).split(':')[1])
    else:
        file_start = 0
        file_stop = len(valid_files) if (args.frange).lower() == 'all' else int(args.frange)

    if file_start >= 0 and file_stop < len(valid_files):
        valid_files = valid_files[file_start:file_stop+1]
        inds_to_use = valid_inds[file_start:file_stop+1]
    else:
        raise IOError("The number of root files available for the %s sample is %i. args.frange must be less than this." % (sample, len(valid_files) ) )

    print(inds_to_use)
    fileset[sample_name] = valid_files

print(fileset)

    ## save output to coffea pkl file
if (args.frange).lower() == 'all':
    outdir = '/'.join([proj_dir, 'results', '%s_%s' % (args.year, jobid), analyzer, args.lepton])
    if args.outfname:
        cfname = args.outfname
    elif args.sample:
        cfname = '%s/%s.coffea' % (outdir, args.sample)
    else:
        cfname = '%s/test_%s.coffea' % (outdir, analyzer)
else:
    if ':' in args.frange:
        outdir = '/'.join([proj_dir, 'results', '%s_%s' % (args.year, jobid), analyzer, args.lepton])
        if args.outfname:
            cfname = args.outfname
        elif args.sample:
            cfname = '%s/%s_%sto%s.coffea' % (outdir, args.sample, file_start, file_stop)
        else:
            cfname = '%s/test_%sto%s.coffea' % (outdir, file_start, file_stop)
    else:
        outdir = proj_dir
        if args.outfname:
            cfname = args.outfname
        elif args.sample:
            cfname = '%s/%s_%s.test.coffea' % (outdir, args.sample, analyzer)
        else:
            cfname = '%s/%s.test.coffea' % (outdir, analyzer)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

    ## get string for passed options
opts = ""
#if args.sample: opts += " --sample=%s" % args.sample
if args.debug: opts += " --debug"


run_cmd = """python {PROJDIR}/bin/{ANALYZER}.py "{FSET}" {YEAR} {LEPTON} {OUTFNAME} {OPTS}""".format(
        PROJDIR=proj_dir,
        ANALYZER=analyzer,
        FSET=fileset,
        YEAR=args.year,
        LEPTON=args.lepton,
        OUTFNAME=cfname,
        OPTS=opts
)
os.system(run_cmd)
