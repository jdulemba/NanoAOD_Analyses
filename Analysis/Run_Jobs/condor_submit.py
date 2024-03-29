import os
import time, datetime
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
            key, value = value.split("=")
            getattr(namespace, self.dest)[key] = value

parser = argparse.ArgumentParser("submit analyzer to the batch queues")
parser.add_argument("analyzer", help="Analyzer to use.")
parser.add_argument("jobdir", help="Directory name to be created in nobackup area.")
parser.add_argument("year", choices=["2016preVFP", "2016postVFP", "2017", "2018"] if os.environ["base_jobid"] == "ULnanoAODv9" else ["2016APV", "2016", "2017", "2018"], help="Specify which year to run over")
parser.add_argument("--opts", nargs="*", action=ParseKwargs, help="Options to pass to analyzers.")
parser.add_argument("--submit", action="store_true", help="Submit jobs")
args = parser.parse_args()

# define dictionary of options to pass
opts_dict = {} if args.opts is None else args.opts
opts_dict["apply_hem"] = opts_dict.get("apply_hem", "True")
opts_dict["evt_sys"] = opts_dict.get("evt_sys", "NONE")
opts_dict["rewt_sys"] = opts_dict.get("rewt_sys", "NONE")
opts_dict["only_sys"] = opts_dict.get("only_sys", "False")
opts_dict["debug"] = opts_dict.get("debug", "False")
opts_dict["jobid"] = opts_dict.get("jobid", os.environ["jobid"])
opts_dict["allowed_masses"] = opts_dict.get("allowed_masses", "All")
opts_dict["allowed_widths"] = opts_dict.get("allowed_widths", "All")
sample = opts_dict.get("sample", None)
if sample: opts_dict.pop("sample")
opts_dict["isCondor"] = "True"

proj_dir = os.environ["PROJECT_DIR"]
jobid, base_jobid = opts_dict["jobid"], os.environ["base_jobid"]
#jobid, base_jobid = opts_dict["jobid"], opts_dict["base_jobid"]
#jobid = os.environ["jobid"]
#base_jobid = os.environ["base_jobid"]
analyzer=args.analyzer
proxy_path = "/afs/cern.ch/work/j/jdulemba/private/x509up_u81826"

    # get jobdir
year, month, day = time.localtime().tm_year, time.localtime().tm_mon, time.localtime().tm_mday
dtime = datetime.datetime(year, month, day)
dtime.strftime("%d%B%Y")
jobdir = "_".join([args.jobdir, dtime.strftime("%d%B%Y"), args.year, jobid])
jobdir = f"BATCH_{jobdir}" if not jobdir.startswith("BATCH") else jobdir
#set_trace()
eos_dir = os.path.join("/eos", "user", os.environ["USER"][0], os.environ["USER"], "NanoAOD_Analyses")
print(f"{os.path.join(eos_dir, jobdir)} written")


def create_batch_job():
    batch_job="""#!/bin/bash

export X509_USER_PROXY=$1

EXE="${{@:2}}"
echo "Executing python Run_Jobs/run_analyzer.py " $EXE from within singularity

singularity exec --bind /afs/cern.ch/work/j/jdulemba/private --bind {PROJECTDIR}:/scratch  --home $PWD:/srv   /cvmfs/unpacked.cern.ch/registry.hub.docker.com/coffeateam/coffea-base:latest   /bin/bash -c "source /scratch/environment.sh {JOBID} && python /scratch/Run_Jobs/run_analyzer.py $EXE"
#singularity exec --bind /afs/cern.ch/work/j/jdulemba/private --bind {PROJECTDIR}:/scratch  --home $PWD:/srv   /cvmfs/unpacked.cern.ch/registry.hub.docker.com/coffeateam/coffea-base:latest   /bin/bash -c "source /scratch/environment.sh && python /scratch/Run_Jobs/run_analyzer.py $EXE"

outfname_arg="${{@: -1}}"
IFS='='
read -a strarr <<< "$outfname_arg"
outfname=${{strarr[1]}}

echo "Copying " $oufname "to {EOSDIR}"
xrdcp $outfname root://eosuser.cern.ch/{EOSDIR}
rm $outfname
""".format(PROJECTDIR=proj_dir, BATCHDIR=os.path.join(proj_dir, jobdir, sample_name), EOSDIR=eos_batch_dir, JOBID=jobid)

    return batch_job

def base_condor_jdl():
#Should_Transfer_Files = YES
#WhenToTransferOutput = ON_EXIT
    condorfile = """universe = vanilla
Executable = {BATCHDIR}/batch_job.sh
+MaxRuntime = 21600
Proxy_path = {PROXYPATH}
Requirements = HasSingularity
""".format(BATCHDIR=os.path.join(proj_dir, jobdir, sample_name), PROXYPATH=proxy_path)
    return condorfile


def add_condor_jobs(idx, frange, sample, opts):
    condorfile = """
Output = con_{IDX}.stdout
Error = con_{IDX}.stderr
Log = con_{IDX}.log
Arguments = $(Proxy_path) {ANALYZER} {FRANGE} {YEAR} {SAMPLE} --opts {OPTS}
Queue
""".format(IDX=idx, ANALYZER=analyzer, FRANGE=frange, YEAR=args.year, SAMPLE=sample, OPTS=opts)
    return condorfile

    ## get samples to use
indir = os.path.join(proj_dir, "inputs", f"{base_jobid}_{args.year}" if base_jobid == "ULnanoAODv9" else f"{args.year}_{base_jobid}")
samples_to_use = tools.get_sample_list(indir=indir, sample=sample) if sample else tools.get_sample_list(indir=indir, text_file="analyzer_inputs.txt")
for sample in samples_to_use:
    if not os.path.isfile(sample):
        raise IOError(f"Sample file {sample}.txt not found")

    sample_name = sample.split("/")[-1].split(".")[0]

        # add output, error, log, arguments for each job splitting
    sfiles = open(sample, "r")
    file_inds = [idx for idx, fname in enumerate([fname.strip("\n") for fname in sfiles if not fname.startswith("#")])]
    splitting = tools.get_file_splitting(sample=sample_name, analyzer=analyzer)
    file_chunks = list(tools.get_file_range(file_inds, splitting))

        # set dir paths in eos space and project dir
    batch_dir = os.path.join(proj_dir, jobdir, sample_name)
    if not os.path.isdir(batch_dir): os.makedirs(batch_dir)
    eos_batch_dir = os.path.join(eos_dir, jobdir, sample_name)
    if not os.path.isdir(eos_batch_dir): os.makedirs(eos_batch_dir)

        ## make condor.jdl file
    condor_cmd = base_condor_jdl()
    for idx, chunk in enumerate(file_chunks):
            # make list of options to pass to analyzers
        tmp_opts_dict = deepcopy(opts_dict)
        tmp_opts_dict["outfname"] = f"{jobdir}_{sample_name}_out_{idx}.coffea"
        opts_list = ["%s=%s" % (key, val) for key, val in tmp_opts_dict.items()]

        condor_cmd += add_condor_jobs(idx, chunk, sample_name, " ".join(opts_list))

    condor_conf = open(os.path.join(batch_dir, "condor.jdl"), "w")
    condor_conf.write(condor_cmd)
    condor_conf.close()

    #set_trace()
        ## make batch_job.sh file
    batch_cmd = create_batch_job()
    batch_conf = open(os.path.join(batch_dir, "batch_job.sh"), "w")
    batch_conf.write(batch_cmd)
    batch_conf.close()
    
    # submit job
    if args.submit:
        orig_dir = os.getcwd()
        print(f"\nSubmitting jobs for {sample_name}")
        os.system("cd " + batch_dir + " && condor_submit condor.jdl")
        os.system("cd " + orig_dir)
