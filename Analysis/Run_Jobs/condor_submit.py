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
parser.add_argument('proxy_file', help='name of x509 file in afs private space to use')
parser.add_argument('--sample', type=str, help='Use specific sample')
args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
nano_dir = os.environ['NANODIR']
jobdir = args.jobdir
analyzer=args.analyzer
proxy_path = '/afs/cern.ch/work/j/jdulemba/private/%s' % args.proxy_file

def create_batch_job():
    batch_job="""#!/bin/bash

export X509_USER_PROXY=$1
#voms-proxy-info -all
#voms-proxy-info -all -file $1

echo "source {NANODIR}/environment.sh"
source {NANODIR}/environment.sh

echo "source {PROJECTDIR}/environment.sh"
source {PROJECTDIR}/environment.sh

EXE="${{@:2}}"
echo "Executing python {PROJECTDIR}/Utilities/run_analyzer.py " $EXE
python {PROJECTDIR}/Utilities/run_analyzer.py $EXE
""".format(NANODIR=nano_dir, PROJECTDIR=proj_dir, ANALYZER=analyzer)

    return batch_job

def base_condor_jdl():
    condorfile = """universe = vanilla
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Executable = {BATCHDIR}/batch_job.sh
+MaxRuntime = 21600
requirements = (OpSysAndVer =?= "CentOS7")
Proxy_path = {PROXYPATH}

""".format(BATCHDIR=batch_dir, PROXYPATH=proxy_path)

    return condorfile

def add_condor_jobs(idx, frange, sample):
    condorfile = """
Output = con_{IDX}.stdout
Error = con_{IDX}.stderr
Log = con_{IDX}.log
Arguments = $(Proxy_path) {ANALYZER} {FRANGE} {YEAR} {LEPTON} --sample={SAMPLE} --outfname={BATCHDIR}/{SAMPLE}_out_{IDX}.coffea
Queue
""".format(IDX=idx, ANALYZER=analyzer, FRANGE=frange, YEAR=args.year, LEPTON=args.lepton, SAMPLE=sample, BATCHDIR=batch_dir)
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

    #set_trace()
    # submit job
    orig_dir = os.getcwd()
    print('\nSubmitting jobs for %s' % sample_name)
    os.system('cd ' + batch_dir + ' && condor_submit condor.jdl')

    os.system('cd ' + orig_dir)

os.system('python %s/Utilities/track_jobs.py %s' % (proj_dir, jobdir))