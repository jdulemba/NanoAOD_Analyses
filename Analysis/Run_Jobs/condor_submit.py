import os, subprocess, time
from pdb import set_trace
import tools
import Utilities.prettyjson as prettyjson
import fnmatch
from copy import deepcopy

import argparse
class ParseKwargs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, dict())
        for value in values:
            key, value = value.split('=')
            getattr(namespace, self.dest)[key] = value

parser = argparse.ArgumentParser('submit analyzer to the batch queues')
parser.add_argument('analyzer', help='Analyzer to use.')
parser.add_argument('jobdir', help='Directory name to be created in nobackup area.')
parser.add_argument('year', choices=['2016APV', '2016', '2017', '2018'], help='Specify which year to run over')
parser.add_argument('--opts', nargs='*', action=ParseKwargs, help='Options to pass to analyzers.')
#parser.add_argument('--sample', type=str, help='Use specific sample')
parser.add_argument('--submit', action='store_true', help='Submit jobs')
#parser.add_argument('--signal', type=str, default='.*', help='Signal sample to use, regex')
#parser.add_argument('--evt_sys', type=str, default='NONE', help='Specify event systematics to run, will be capitalized. Default is NONE')
#parser.add_argument('--rewt_sys', type=str, default='NONE', help='Specify reweighting systematics to run, will be capitalized. Default is NONE')
#parser.add_argument('--only_sys', type=int, default=0, help='Only run specified systematics and not nominal weights (nosys)')
args = parser.parse_args()

#set_trace()
# define dictionary of options to pass
opts_dict = {} if args.opts is None else args.opts
opts_dict['apply_hem'] = opts_dict.get('apply_hem', 'True')
opts_dict['signal'] = opts_dict.get('signal', ".*")
opts_dict['evt_sys'] = opts_dict.get('evt_sys', 'NONE')
opts_dict['rewt_sys'] = opts_dict.get('rewt_sys', 'NONE')
opts_dict['only_sys'] = opts_dict.get('only_sys', 'False')
opts_dict['debug'] = opts_dict.get('debug', 'False')

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
base_jobid = os.environ['base_jobid']
nano_dir = os.environ['NANODIR']
jobdir = args.jobdir
analyzer=args.analyzer
proxy_path = '/afs/cern.ch/work/j/jdulemba/private/x509up_u81826'

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
echo "Executing python {PROJECTDIR}/Run_Jobs/run_analyzer.py " $EXE
python {PROJECTDIR}/Run_Jobs/run_analyzer.py $EXE
""".format(NANODIR=nano_dir, PROJECTDIR=proj_dir, ANALYZER=analyzer)

    return batch_job

def base_condor_jdl():
    condorfile = """universe = vanilla
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Executable = {BATCHDIR}/batch_job.sh
+MaxRuntime = 10800
requirements = (OpSysAndVer =?= "CentOS7")
Proxy_path = {PROXYPATH}
""".format(BATCHDIR=batch_dir, PROXYPATH=proxy_path)
    return condorfile

def add_signal_condor_jobs(idx, frange, sample, signal):
    sig_outname = 'AHtoTT' if signal == parser.get_default('signal') else signal
    condorfile = """
Output = con_{IDX}.stdout
Error = con_{IDX}.stderr
Log = con_{IDX}.log
Arguments = $(Proxy_path) {ANALYZER} {FRANGE} {YEAR} --sample={SAMPLE} --signal={SIGNAL} --evt_sys={EVTSYS} --rewt_sys={REWTSYS} --only_sys={ONLYSYS} --outfname={BATCHDIR}/{SAMPLE}_{SIGOUTNAME}_out_{IDX}.coffea
Queue
""".format(IDX=idx, ANALYZER=analyzer, FRANGE=frange, YEAR=args.year, SIGNAL=signal, SAMPLE=sample, EVTSYS=args.evt_sys, REWTSYS=args.rewt_sys, ONLYSYS=args.only_sys, BATCHDIR=batch_dir, SIGOUTNAME=sig_outname)
    return condorfile

def add_condor_jobs(idx, frange, opts):
    condorfile = """
Output = con_{IDX}.stdout
Error = con_{IDX}.stderr
Log = con_{IDX}.log
Arguments = $(Proxy_path) {ANALYZER} {FRANGE} {YEAR} --opts {OPTS}
Queue
""".format(IDX=idx, ANALYZER=analyzer, FRANGE=frange, YEAR=args.year, OPTS=opts)
    return condorfile

    ## get samples to use
indir = os.path.join(proj_dir, 'inputs', '%s_%s' % (args.year, base_jobid))
#set_trace()
samples_to_use = tools.get_sample_list(indir=indir, sample=opts_dict['sample']) if 'sample' in opts_dict.keys() else tools.get_sample_list(indir=indir, text_file='analyzer_inputs.txt')
for sample in samples_to_use:
    if not os.path.isfile(sample):
        raise IOError("Sample file %s.txt not found" % sample)

    sample_name = sample.split('/')[-1].split('.')[0]

        # add output, error, log, arguments for each job splitting
    sfiles = open(sample, 'r')
    file_inds = [idx for idx, fname in enumerate([fname.strip('\n') for fname in sfiles if not fname.startswith('#')])]
    splitting = tools.get_file_splitting('AtoTT') if analyzer == 'htt_signal_reweight' else tools.get_file_splitting(sample.split('/')[-1].split('.')[0])
    file_chunks = list(tools.get_file_range(file_inds, splitting))

    if (analyzer == 'htt_signal_reweight') or (analyzer == 'signal_validation'):
            ## make batch_job.sh file
        batch_dir = '%s/%s/%s_%s' % (proj_dir, jobdir, sample_name, args.signal) if not args.signal == parser.get_default('signal') else '%s/%s/%s_AHtoTT' % (proj_dir, jobdir, sample_name)
        if not os.path.isdir(batch_dir): os.makedirs(batch_dir)
        batch_cmd = create_batch_job()
        batch_conf = open(os.path.join(batch_dir, 'batch_job.sh'), 'w')
        batch_conf.write(batch_cmd)
        batch_conf.close()
        
            ## make condor.jdl file
        condor_cmd = base_condor_jdl()

        for idx, chunk in enumerate(file_chunks):
            condor_cmd += add_signal_condor_jobs(idx, chunk, sample.split('/')[-1].split('.')[0], args.signal)

        condor_conf = open(os.path.join(batch_dir, 'condor.jdl'), 'w')
        condor_conf.write(condor_cmd)
        condor_conf.close()

        #set_trace()
        # submit job
        if args.submit:
            orig_dir = os.getcwd()
            print('\nSubmitting jobs for %s' % args.signal) if not args.signal == parser.get_default('signal') else print('\nSubmitting jobs for all signal')
            os.system('cd ' + batch_dir + ' && condor_submit condor.jdl')
            os.system('cd ' + orig_dir)

    else:
            ## make batch_job.sh file
        batch_dir = os.path.join(proj_dir, jobdir, sample_name)
        if not os.path.isdir(batch_dir): os.makedirs(batch_dir)
        batch_cmd = create_batch_job()
        batch_conf = open(os.path.join(batch_dir, 'batch_job.sh'), 'w')
        batch_conf.write(batch_cmd)
        batch_conf.close()
        
            ## make condor.jdl file
        condor_cmd = base_condor_jdl()
        for idx, chunk in enumerate(file_chunks):
                # make list of options to pass to analyzers
            tmp_opts_dict = deepcopy(opts_dict)
            tmp_opts_dict['sample'] = sample_name
            tmp_opts_dict['outfname'] = os.path.join(batch_dir, '%s_out_%s.coffea' % (sample_name, idx))
            opts_list = ["%s=%s" % (key, val) for key, val in tmp_opts_dict.items()]

            condor_cmd += add_condor_jobs(idx, chunk, " ".join(opts_list))

        condor_conf = open(os.path.join(batch_dir, 'condor.jdl'), 'w')
        condor_conf.write(condor_cmd)
        condor_conf.close()

        # submit job
        if args.submit:
            orig_dir = os.getcwd()
            print('\nSubmitting jobs for %s' % sample_name)
            os.system('cd ' + batch_dir + ' && condor_submit condor.jdl')
            os.system('cd ' + orig_dir)

