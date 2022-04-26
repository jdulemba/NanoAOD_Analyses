import os, subprocess, time
from pdb import set_trace
import Utilities.plot_tools as plot_tools

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("jobdir", help="specify name of directory where jobs were submitted")
parser.add_argument("--resub", action="store_false", help="Resubmit failed jobs, default is True.")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
user = os.environ["USER"]

#print(args.resub)
#set_trace()
## once jobs are finished, merge all output coffea files for each sample into one
condor_jobdir = args.jobdir if (args.jobdir).startswith(os.path.join("afs", "cern.ch", "work", "j", user)) else os.path.join(proj_dir, args.jobdir)
condor_jobdir = condor_jobdir[:-1] if condor_jobdir.endswith("/") else condor_jobdir
if not os.path.isdir(condor_jobdir):
    raise IOError(f"Directory: {condor_jobdir} does not exist")

basedir = os.path.basename(condor_jobdir)
eos_dir = os.path.join("/eos", "user", user[0], user, "NanoAOD_Analyses")
eos_jobdir = os.path.join(eos_dir, basedir)
if os.path.isfile(os.path.join(eos_jobdir, f"{args.jobdir}_TOT.coffea")):
    print(f"Combined output file {eos_jobdir}/{args.jobdir}_TOT.coffea already exists")
    import sys; sys.exit()


orig_dir = os.getcwd()
samples = [dirname for dirname in os.listdir(condor_jobdir) if os.path.isdir(os.path.join(condor_jobdir, dirname))]

def check_correctness(sample, dump_rescue=False):
    #set_trace()
    isCorrect = True
    if os.path.isfile(f"{eos_jobdir}/{sample}/{sample}_TOT.coffea"):  ## TOT output file already created
        return isCorrect, 0, True
    jdl = open(f"{condor_jobdir}/{sample}/condor.jdl").read()
    blocks = jdl.split("\n\n")
    header = blocks[0]
    block_map = {}
    #set_trace()
    for block in blocks[1:]:
            ## save output filename to check if it exists
        fname = block.split("Arguments = ")[1].split(" ")[-1].split("=")[-1].split("\n")[0]
        block_map[fname] = block

    npass = 0
    fails = []
    for fname, arguments in block_map.items():
        if not os.path.isfile(os.path.join(eos_jobdir, sample, fname)):
            fails.append(fname)
        else:
            npass += 1

    print('''Run Summary for %s:
Successful jobs: %d
Failed jobs: %d
''' % (sample, npass, len(fails)))

    if fails:
        isCorrect = False
        if dump_rescue:
            print("dumping rescue job")
            with open(f"{condor_jobdir}/{sample}/condor.rescue.jdl", "w") as rescue:
                rescue.write(header)
                rescue.write("\n\n")
                rescue.write(
                    "\n\n".join([block_map[key] for key in fails])
                    )

    return isCorrect, len(fails), False

tot_files = {}

escape = False
finished_samples = []

while not escape:
    stdout= subprocess.Popen("condor_q $USER -nobatch", stdout=subprocess.PIPE, shell=True).stdout.read()
    njobs = 0
    samples_status = {}

        # find files that are already finished
    if os.path.isfile(os.path.join(condor_jobdir, "to_merge.txt")):
        to_merge_txt = open(os.path.join(condor_jobdir, "to_merge.txt")).read()
        finished_samples = to_merge_txt.split("\n")

        ## this loop is just to check how many jobs are still running for each sample
    for sample in samples:
        if sample in finished_samples: continue
            ## checks how many lines correspond to "outfname=sample_out"
        sample_njobs = len([line for line in stdout.split(b"\n") if ((f"outfname={basedir}_{sample}_out").encode() in line)])
        if sample_njobs == 0:
            isCorrect, nfails, alreadyComplete = check_correctness(sample=sample, dump_rescue=args.resub)
            if isCorrect == True:
                status = "MERGED" if alreadyComplete else "COMPLETE" # jobs are already merged or (completed but still need to be merged)
                finished_samples.append(sample)
            else:
                status = "RESUB" # jobs need to be resubmitted
            samples_status[sample] = (status, nfails)
        else:
                ## get lines corresponding to sample
            job_lines = [line.decode() for line in stdout.split(b"\n") if (f"outfname={basedir}_{sample}_out").encode() in line]
                ## get job statuses
            job_statuses = ["".join(line).split()[5] for line in job_lines] # 5 is hardcoded
                ## check how many jobs haven"t been removed
            njobs_running = len([status for status in job_statuses if not (status == "X")])
            if njobs_running == 0:
                isCorrect, nfails, alreadyComplete = check_correctness(sample=sample, dump_rescue=args.resub)
                if isCorrect == True:
                    status = "MERGED" if alreadyComplete else "COMPLETE" # jobs are already merged or (completed but still need to be merged)
                    finished_samples.append(sample)
                else:
                    status = "RESUB" # jobs need to be resubmitted
                samples_status[sample] = (status, nfails)
            else:
                samples_status[sample] = ("RUNNING", njobs_running) # jobs running
                njobs += njobs_running

        # write finished samples to txt file
    if finished_samples:
        to_merge_txt = open(os.path.join(condor_jobdir, "to_merge.txt"),"w")
        to_merge_txt.write("\n".join(finished_samples))
        to_merge_txt.close()
        print(f"{os.path.join(condor_jobdir, 'to_merge.txt')} written\n")
        #set_trace()

        ## this loop is responsible for merging or resubmitting jobs based on status
    for sample, (status, sample_njobs) in samples_status.items():
            ## if status is "RESUB", resubmit jobs to condor
        if status == "RESUB":
            #set_trace()
            print(f"\t{sample} needs to be resubmitted")
            os.system(f"cd {condor_jobdir}/{sample} && condor_submit condor.rescue.jdl && cd {orig_dir}")
            samples_status[sample] = ("RUNNING", sample_njobs)
            njobs += sample_njobs

            ## if status is "RUNNING" do nothing, jobs are already running
        elif status == "RUNNING":
            print("\t%s has %i jobs running" % (sample, sample_njobs) )

        ## recheck all jobs if there are none running
    if njobs == 0:
        failed_samples = []
        finished_samples = [] # reinitialize finished_samples to check if they're actually finished
        for sample in samples:
            isCorrect, nfails, alreadyComplete = check_correctness(sample=sample, dump_rescue=args.resub)
            if isCorrect == False: 
                failed_samples.append(sample)
                os.system(f"cd {condor_jobdir}/{sample} && condor_submit condor.rescue.jdl && cd {orig_dir}")
            else:
                finished_samples.append(sample)
        if len(failed_samples) == 0:
            print("All jobs completed")
            escape = True
        else:
            time.sleep(30)
    else:
        print("%i jobs are still running, checking again in 30 seconds\n" % njobs )
        time.sleep(30)

## merge all TOT files from each sample into one
if sorted(set(tot_files.keys())) == sorted(set(samples)):
    print("All jobs completed. Run merge_jobs.py inside singularity to merge all jobs")
