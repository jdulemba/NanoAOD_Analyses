import os, subprocess, time
from pdb import set_trace
import tools
import Utilities.prettyjson as prettyjson
import fnmatch

from argparse import ArgumentParser
parser = ArgumentParser('submit analyzer to the batch queues')
parser.add_argument('analyzer', help='Analyzer to use.')
parser.add_argument('jobdir', help='Directory name to be created in nobackup area.')
parser.add_argument('year', choices=['2016', '2017', '2018'], help='Specify which year to run over')
parser.add_argument('lepton', choices=['Electron', 'Muon'], help='Choose which lepton to select')
parser.add_argument('--sample', type=str, help='Use specific sample')
args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
nano_dir = os.environ['NANODIR']
jobdir = args.jobdir
analyzer=args.analyzer

def create_batch_job():
    batch_job="""#!/bin/bash

echo "source {NANODIR}/environment.sh"
source {NANODIR}/environment.sh

echo "source {PROJECTDIR}/environment.sh"
source {PROJECTDIR}/environment.sh

EXE="$@"
echo "Executing python {PROJECTDIR}/bin/{ANALYZER}.py " $EXE
python {PROJECTDIR}/bin/presel_analyzer.py $EXE
""".format(NANODIR=nano_dir, PROJECTDIR=proj_dir, ANALYZER=analyzer)

    return batch_job

def base_condor_jdl():
    condorfile = """universe = vanilla
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Executable = {BATCHDIR}/batch_job.sh
+MaxRuntime = 21600
requirements = (OpSysAndVer =?= "CentOS7")

""".format(BATCHDIR=batch_dir)

    return condorfile

def add_condor_jobs(idx, frange, sample):
    condorfile = """
Output = con_{IDX}.stdout
Error = con_{IDX}.stderr
Log = con_{IDX}.log
Arguments = {FRANGE} {YEAR} {LEPTON} --sample={SAMPLE}
Queue
""".format(IDX=idx, FRANGE=frange, YEAR=args.year, LEPTON=args.lepton, SAMPLE=sample)
    return condorfile

    ## get samples to use
indir = '/'.join([proj_dir, 'inputs', '%s_%s' % (args.year, jobid)])
samples_to_use = tools.get_sample_list(indir=indir, sample=args.sample) if args.sample else tools.get_sample_list(indir=indir, text_file='%s_inputs.txt' % analyzer)
for sample in samples_to_use:
    if not os.path.isfile(sample):
        raise IOError("Sample file %s.txt not found" % sample)

    sample_name = sample.split('/')[-1].split('.')[0]
        ## make batch_job.sh file
    batch_dir = '%s/%s/%s' % (proj_dir, jobdir, sample_name)
    if not os.path.isdir(batch_dir): os.makedirs(batch_dir)
    batch_cmd = create_batch_job()
    batch_conf = open(os.path.join(batch_dir, 'batch_job.sh'), 'w')
    batch_conf.write(batch_cmd)
    batch_conf.close()
    
        ## make condor.jdl file
    condor_cmd = base_condor_jdl()

        # add output, error, log, arguments for each job splitting
    sfiles = open(sample, 'r')
    file_inds = [idx for idx, fname in enumerate([fname.strip('\n') for fname in sfiles if not fname.startswith('#')])]
    splitting = tools.get_file_splitting(sample.split('/')[-1].split('.')[0])
    file_chunks = list(tools.get_file_range(file_inds, splitting))
    for idx, chunk in enumerate(file_chunks):
        condor_cmd += add_condor_jobs(idx, chunk, sample.split('/')[-1].split('.')[0])

    condor_conf = open(os.path.join(batch_dir, 'condor.jdl'), 'w')
    condor_conf.write(condor_cmd)
    condor_conf.close()
