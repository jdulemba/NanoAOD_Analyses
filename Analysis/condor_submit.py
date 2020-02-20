import os, subprocess, time
from pdb import set_trace
from argparse import ArgumentParser

parser = ArgumentParser('submit analyzer to the batch queues')
parser.add_argument('jobdir', help='Directory name to be created in nobackup area.')

args = parser.parse_args()

user = os.environ['USER']
eos_base_dir = '/eos/uscms/store/user/%s' % user
noback_base_dir = '/uscms_data/d3/%s' % user

jobdir = args.jobdir

def zip_release(fname):

    timestamp = str(time.time()).split('.')[0]
    gz_name = '_'.join([fname, timestamp])
    fname = '%s/%s' % (eos_base_dir, gz_name)
    #eos_dir = '%s/%s' % (eos_base_dir, fname)
    if os.path.isdir(fname):
        raise IOError('%s already exists! Change name' % fname)

    #set_trace()
    zip_cmd = "time tar --exclude-caches-all --exclude-vcs --exclude-from=../.gitignore -zcf {FNAME}.tar.gz -C {BASE} miniconda3 ".format(
    #zip_cmd = "time tar --exclude-caches-all --exclude-vcs --exclude-from=../.gitignore -zcf {FNAME}.tar.gz -C {BASE}/miniconda3 NanoAOD_Analyses ".format(
        FNAME=fname, BASE=noback_base_dir)

    return zip_cmd, gz_name

def create_batch_job():
    batch_job="""#!/bin/bash
WORKINGDIR=$PWD
echo "WORKINGDIR: "$WORKINGDIR

xrdcp -s root://cmseos.fnal.gov//store/user/{USER}/{TAR_FNAME}.tar.gz .

tar -zxf {TAR_FNAME}.tar.gz
rm {TAR_FNAME}.tar.gz

ls -lht
export PATH="miniconda3/bin:$PATH"

echo "cd $WORKINGDIR/miniconda3/NanoAOD_Analyses/Test"
cd $WORKINGDIR/miniconda3/NanoAOD_Analyses/Test
ls -lht

echo "source environment.sh"
source environment.sh

python helloworld.py
#EXE=$1
#echo "Executing " $EXE
#./runMadGraph.sh $EXE
#
##ls -lht
#
#xrdcp -r {JOBDIR} root://cmseos.fnal.gov//store/user/{USER}/{TAR_FNAME}
#xrdcp -s {JOBDIR}*.txt root://cmseos.fnal.gov//store/user/{USER}/{TAR_FNAME}
""".format( USER=user, TAR_FNAME=zip_fname, JOBDIR=jobdir)

    return batch_job

def create_condor_jdl():
    condorfile = """universe = vanilla
Executable = batch_job.sh
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT

Output = con_0.stdout
Error = con_0.stderr
Log = con_0.log
Arguments = scripts/test_analyzer.py --nfiles=1
request_memory = 5000
Queue
"""

    return condorfile    



zip_dir_cmd, zip_fname = zip_release('BATCH_NanoAOD_%s' % jobdir)

print(zip_dir_cmd)
#set_trace()
batch_dir = '%s/%s' % (noback_base_dir, zip_fname)
if not os.path.isdir(batch_dir): os.mkdir(batch_dir)
#condor_cmd = """python Condor_Sub/make_condor_file.py {BATCHDIR} {CFG}""".format(BATCHDIR=batch_dir, CFG='%s/configs/%s' % (tt_dir, jobdir))
batch_cmd = create_batch_job()
batch_conf = open(os.path.join(batch_dir, 'batch_job.sh'), 'w')
batch_conf.write(batch_cmd)
batch_conf.close()

condor_cmd = create_condor_jdl()
condor_conf = open(os.path.join(batch_dir, 'condor.jdl'), 'w')
condor_conf.write(condor_cmd)
condor_conf.close()

