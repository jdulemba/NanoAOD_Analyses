import os
import time, datetime
from pdb import set_trace
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
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="Specify which year to run over")
parser.add_argument("templates_to_run", type=str, help="Choose which type of templates to run, multiple options can be input as ':' separated strings.")
parser.add_argument("--opts", nargs="*", action=ParseKwargs, help="Options to pass to analyzers.")
parser.add_argument("--MEopts", nargs="*", action=ParseKwargs, help="Options to pass for ME reweight sig.")
parser.add_argument("--submit", action="store_true", help="Submit jobs")
args = parser.parse_args()

proxy_path = "/afs/cern.ch/work/j/jdulemba/private/x509up_u81826"

def create_batch_job():
    batch_job="""#!/bin/bash

EXE="$@"
#EXE="${{@:2}}"
echo "Executing python Templates/$EXE" from within singularity

#singularity exec --bind /afs/cern.ch/work/j/jdulemba/private --bind {PROJECTDIR}:/scratch  --home $PWD:/srv -B /eos/user/j/jdulemba/NanoAOD_Analyses  /cvmfs/unpacked.cern.ch/registry.hub.docker.com/coffeateam/coffea-base:latest   /bin/bash -c "source /scratch/environment.sh && python /scratch/Templates/$EXE"
singularity exec --bind /afs/cern.ch/work/j/jdulemba/private --bind {PROJECTDIR}:/scratch  --home $PWD:/srv   /cvmfs/unpacked.cern.ch/registry.hub.docker.com/coffeateam/coffea-base:latest   /bin/bash -c "source /scratch/environment.sh && python /scratch/Templates/$EXE"
""".format(PROJECTDIR=proj_dir)

    return batch_job

def condor_jdl(opts, MEopts):
    condorfile = """universe = vanilla
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Executable = {BATCHDIR}/batch_job.sh
+MaxRuntime = 43200
Proxy_path = {PROXYPATH}
Requirements = HasSingularity

Output = con_0.stdout
Error = con_0.stderr
Log = con_0.log
Arguments = $(Proxy_path) {ANALYZER}.py {YEAR} {TEMPLATES_TO_RUN} {OPTS} {ME_OPTS}
Queue
""".format(BATCHDIR=batch_dir, ANALYZER=analyzer, YEAR=args.year, TEMPLATES_TO_RUN=templates_to_run, OPTS=opts, ME_OPTS=MEopts, PROXYPATH=proxy_path)
    return condorfile


def rfile_create_batch_job():
    batch_job="""#!/bin/bash

export X509_USER_PROXY=$1
EXE="${{@:2}}"
echo "Executing python $PROJECT_DIR/Templates/" $EXE from within singularity

singularity exec --bind /afs/cern.ch/work/j/jdulemba/private --bind {PROJECTDIR}:/scratch  --home $PWD:/srv -B /eos/user/j/jdulemba/NanoAOD_Analyses  /cvmfs/unpacked.cern.ch/registry.hub.docker.com/coffeateam/coffea-base:latest   /bin/bash -c "source /scratch/environment.sh && python /scratch/Templates/$EXE"
#singularity exec --bind /afs/cern.ch/work/j/jdulemba/private --bind {PROJECTDIR}:/scratch  --home $PWD:/srv   /cvmfs/unpacked.cern.ch/registry.hub.docker.com/coffeateam/coffea-base:latest   /bin/bash -c "source /scratch/environment.sh && python /scratch/Templates/$EXE"

""".format(PROJECTDIR=proj_dir)

    return batch_job

def rfile_condor_jdl(opts):
    condorfile = """universe = vanilla
Should_Transfer_Files = YES
WhenToTransferOutput = ON_EXIT
Executable = {BATCHDIR}/batch_job.sh
+MaxRuntime = 43200
Proxy_path = {PROXYPATH}
Requirements = HasSingularity

Output = con_0.stdout
Error = con_0.stderr
Log = con_0.log
Arguments = $(Proxy_path) {ANALYZER}.py {YEAR} {SIG_TYPE} {OPTS}
Queue
""".format(BATCHDIR=batch_dir, ANALYZER=analyzer, YEAR=args.year, SIG_TYPE=args.opts["sig_type"], OPTS=opts, PROXYPATH=proxy_path)
    return condorfile


if __name__ == "__main__":
    allowed_template_options = ["bkg", "sig", "MEreweight_sig", "PDF"]
    templates_to_run = [template for template in (args.templates_to_run).split(":") if template in allowed_template_options]
    templates_to_not_run = [template for template in (args.templates_to_run).split(":") if template not in allowed_template_options]
    if templates_to_not_run:
        print(f"{templates_to_not_run} are not valid options for making templates, will be skipped")
    templates_to_run = ":".join(templates_to_run)

    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]
    base_jobid = os.environ["base_jobid"]
    analyzer = args.analyzer
    
    # define dictionary of options to pass
    opts_dict = {} if args.opts is None else args.opts
    MEopts_dict = {} if args.MEopts is None else args.MEopts

    #set_trace()
    if analyzer == "convert_Otto_signal_to_combine":
        opts_dict["dir"] = opts_dict.get("dir")
    else:
        opts_dict["maskData"] = opts_dict.get("maskData")
    
        # get jobdir
    year, month, day = time.localtime().tm_year, time.localtime().tm_mon, time.localtime().tm_mday
    dtime = datetime.datetime(year, month, day)
    dtime.strftime("%d%B%Y")
    jobdir = "_".join([args.jobdir, dtime.strftime("%d%B%Y"), args.year, jobid])
    jobdir = f"BATCH_{jobdir}" if not jobdir.startswith("BATCH") else jobdir

        ## make batch_job.sh file
    batch_dir = os.path.join(proj_dir, jobdir)
    if not os.path.isdir(batch_dir): os.makedirs(batch_dir)
    print(f"{batch_dir} written")
    #set_trace()
    batch_cmd = rfile_create_batch_job() if ((analyzer == "format_final_templates") or (analyzer == "convert_Otto_signal_to_combine")) else create_batch_job()
    batch_conf = open(os.path.join(batch_dir, "batch_job.sh"), "w")
    batch_conf.write(batch_cmd)
    batch_conf.close()

    #set_trace()    
        # make list of options to pass to analyzers
    tmp_opts_dict = deepcopy(opts_dict)
    tmp_MEopts_dict = deepcopy(MEopts_dict)
    if analyzer == "convert_Otto_signal_to_combine":
        opts_list = []
        for key, val in tmp_opts_dict.items():
            if key == "sig_type": continue
            if key == "dir": opts_list.append(f"--{key}={val}")
            else: opts_list.append(f"--{key}")
    else:
        tmp_opts_dict["outfname"] = f"{jobdir}_{analyzer}_out_0.coffea"
        opts_list = []
        MEopts_list = ["--MEreweight_opts"]
        for key, val in tmp_opts_dict.items():
            if val == "True": opts_list.append(f"--{key}")
        for key, val in tmp_MEopts_dict.items():
            MEopts_list.append(f"{key}={val}")
            

        ## make condor.jdl file
    condor_cmd = rfile_condor_jdl(" ".join(opts_list)) if ((analyzer == "format_final_templates") or (analyzer == "convert_Otto_signal_to_combine")) else condor_jdl(" ".join(opts_list), " ".join(MEopts_list))
    #condor_cmd = rfile_condor_jdl(" ".join(opts_list)) if ((analyzer == "format_final_templates") or (analyzer == "convert_Otto_signal_to_combine")) else condor_jdl(" ".join(opts_list))
    
    condor_conf = open(os.path.join(batch_dir, "condor.jdl"), "w")
    condor_conf.write(condor_cmd)
    condor_conf.close()
    
    # submit job
    if args.submit:
        orig_dir = os.getcwd()
        print(f"\nSubmitting jobs for {analyzer}")
        os.system("cd " + batch_dir + " && condor_submit condor.jdl")
        os.system("cd " + orig_dir)
