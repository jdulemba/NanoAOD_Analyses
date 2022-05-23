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

def check_correctness(sample):
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
        fname = [key for key in block.split("\n") if "Arguments" in key][0].split("=")[-1]
        block_map[fname] = block

    #set_trace()
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

    ##set_trace()
    #if fails:
    #    isCorrect = False
    #    #if dump_rescue:
    #    print("dumping rescue job")
    #    with open(f"{condor_jobdir}/{sample}/condor.rescue.jdl", "w") as rescue:
    #        rescue.write(header)
    #        rescue.write("\n\n")
    #        rescue.write(
    #            "\n\n".join([block_map[key] for key in fails])
    #            )

    rewrite_condor = False
    if fails:
        isCorrect = False
        print("dumping rescue job")
        #set_trace()
        #jobs_to_submit = []
        ##if "+MaxRuntime = 10800" in header:
        ##    rewrite_condor = True
        ##    header = header.replace("+MaxRuntime = 10800", "+MaxRuntime = 21600")

        ##jobs_to_remove_from_fails = []
        #for key in fails:
        #    #set_trace()
        #    #if "isCondor=True" not in block_map[key]:
        #    #    #set_trace()
        #    #    jobs_to_submit.append(block_map[key].replace("--opts", "--opts isCondor=True"))

        #    #if ("evt_sys=*") and ("rewt_sys=*") in block_map[key]:
        #    #    rewrite_condor = True
        #    #    condor_job_num = block_map[key].split(".log")[0].split("Log = con")[-1]
        #    #        # split condor job into different reweighting systs
        #    #    jobs_to_submit.append(block_map[key].replace("rewt_sys=*", "rewt_sys=LEP*").replace("evt_sys=*", "evt_sys=NONE"))
        #    #    jobs_to_submit.append(block_map[key].replace("rewt_sys=*", "rewt_sys=[!LEP]*").replace("evt_sys=*", "evt_sys=NONE").replace("only_sys=False", "only_sys=True").replace(f"{condor_job_num}.", f"{condor_job_num}_0."))
        #    #        # split condor job into different event systs
        #    #    jobs_to_submit.append(block_map[key].replace("rewt_sys=*", "rewt_sys=NONE").replace("evt_sys=*", "evt_sys=JER_*").replace("only_sys=False", "only_sys=True").replace(f"{condor_job_num}.", f"{condor_job_num}_1."))
        #    #    jobs_to_submit.append(block_map[key].replace("rewt_sys=*", "rewt_sys=NONE").replace("evt_sys=*", "evt_sys=MET_*").replace("only_sys=False", "only_sys=True").replace(f"{condor_job_num}.", f"{condor_job_num}_2."))
        #    #    jobs_to_submit.append(block_map[key].replace("rewt_sys=*", "rewt_sys=NONE").replace("evt_sys=*", "evt_sys=JES_*").replace("only_sys=False", "only_sys=True").replace(f"{condor_job_num}.", f"{condor_job_num}_3."))
        #    #    block_map[key.replace(f"{condor_job_num}.", f"{condor_job_num}_0.")] = block_map[key].replace("rewt_sys=*", "rewt_sys=[!LEP]*").replace("evt_sys=*", "evt_sys=NONE").replace("only_sys=False", "only_sys=True").replace(f"{condor_job_num}.", f"{condor_job_num}_0.")
        #    #    block_map[key.replace(f"{condor_job_num}.", f"{condor_job_num}_1.")] = block_map[key].replace("rewt_sys=*", "rewt_sys=NONE").replace("evt_sys=*", "evt_sys=JER_*").replace("only_sys=False", "only_sys=True").replace(f"{condor_job_num}.", f"{condor_job_num}_1.")
        #    #    block_map[key.replace(f"{condor_job_num}.", f"{condor_job_num}_2.")] = block_map[key].replace("rewt_sys=*", "rewt_sys=NONE").replace("evt_sys=*", "evt_sys=MET_*").replace("only_sys=False", "only_sys=True").replace(f"{condor_job_num}.", f"{condor_job_num}_2.")
        #    #    block_map[key.replace(f"{condor_job_num}.", f"{condor_job_num}_3.")] = block_map[key].replace("rewt_sys=*", "rewt_sys=NONE").replace("evt_sys=*", "evt_sys=JES_*").replace("only_sys=False", "only_sys=True").replace(f"{condor_job_num}.", f"{condor_job_num}_3.")
        #    #    block_map[key] = block_map[key].replace("rewt_sys=*", "rewt_sys=LEP*").replace("evt_sys=*", "evt_sys=NONE")

        #    ##if "evt_sys=NONE" in block_map[key]:
        #    #    #set_trace()
        #    #    rewrite_condor = True
        #    #    del block_map[key]
        #    #    jobs_to_remove_from_fails.append(key)
        #    #    print("Removing job from condor.jdl file")
        #    #    continue

        #    #if "rewt_sys=*" in block_map[key]:
        #    #    rewrite_condor = True
        #    #    jobs_to_submit.append(block_map[key].replace("rewt_sys=*", "rewt_sys=LEP*"))
        #    #    jobs_to_submit.append(block_map[key].replace("rewt_sys=*", "rewt_sys=[!LEP]*").replace("only_sys=False", "only_sys=True").replace("0_0", "0_0p1"))
        #    #    #set_trace()
        #    #    block_map[key.replace("0_0", "0_0p1")] = block_map[key].replace("rewt_sys=*", "rewt_sys=[!LEP]*").replace("only_sys=False", "only_sys=True").replace("0_0", "0_0p1")
        #    #    block_map[key] = block_map[key].replace("rewt_sys=*", "rewt_sys=LEP*")
        #    else:
        #        jobs_to_submit.append(block_map[key])

        #set_trace()
        #[fails.remove(key) for key in jobs_to_remove_from_fails]
        with open(f"{condor_jobdir}/{sample}/condor.rescue.jdl", "w") as rescue:
            rescue.write(header)
            rescue.write("\n\n")
            rescue.write(
                "\n\n".join([block_map[key] for key in fails])
                )
        #if jobs_to_submit:
        #    with open(f"{condor_jobdir}/{sample}/condor.rescue.jdl", "w") as rescue:
        #        rescue.write(header)
        #        rescue.write("\n\n")
        #        rescue.write(
        #            "\n\n".join([block_map[key] for key in fails])
        #            #"\n\n".join(jobs_to_submit)
        #            )

    if rewrite_condor:
        print(f"Rewriting {condor_jobdir}/{sample}/condor.jdl")
        with open(f"{condor_jobdir}/{sample}/condor.jdl", "w") as updated_jdl:
            updated_jdl.write(header)
            updated_jdl.write("\n\n")
            updated_jdl.write(
                "\n\n".join(sorted(block_map.values()))
                )

    #set_trace()        
    if (len(fails) == 0) and (not isCorrect):
        #set_trace()
        isCorrect = True

    return isCorrect, len(fails), False

tot_files = {}

escape = False
finished_samples = []

job_status_meanings = {
    "I" : "Idle",
    "R" : "Running",
    "X" : "Removed",
    "H" : "Held",
}

def get_statuses(job_statuses):
    status_dict = {}
    for key, val in job_status_meanings.items():
        status_dict[val] = len([stat for stat in job_statuses if (stat == key)])
    job_status_list = [f"{key}={val}" for key, val in status_dict.items()]
    return job_status_list


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
        #set_trace()
        if sample_njobs == 0:
            isCorrect, nfails, alreadyComplete = check_correctness(sample=sample)
            if isCorrect == True:
                status = "MERGED" if alreadyComplete else "COMPLETE" # jobs are already merged or (completed but still need to be merged)
                finished_samples.append(sample)
            else:
                status = "RESUB" # jobs need to be resubmitted
            samples_status[sample] = [f"{status}={nfails}"]
        else:
                ## get lines corresponding to sample
            job_lines = [line.decode() for line in stdout.split(b"\n") if (f"outfname={basedir}_{sample}_out").encode() in line]
                ## get job statuses
            job_statuses = ["".join(line).split()[5] for line in job_lines] # 5 is hardcoded
                ## check how many jobs haven"t been removed
            njobs_running = len([status for status in job_statuses if not (status == "X")])
            if njobs_running == 0:
                isCorrect, nfails, alreadyComplete = check_correctness(sample=sample)
                if isCorrect == True:
                    status = "MERGED" if alreadyComplete else "COMPLETE" # jobs are already merged or (completed but still need to be merged)
                    finished_samples.append(sample)
                else:
                    status = "RESUB" # jobs need to be resubmitted
                samples_status[sample] = [f"{status}={nfails}"]
            else:
                job_status_list = get_statuses(job_statuses)
                samples_status[sample] = job_status_list # jobs running
                njobs += njobs_running

    #set_trace()
        # write finished samples to txt file
    if finished_samples:
        #set_trace()
        to_merge_txt = open(os.path.join(condor_jobdir, "to_merge.txt"),"w")
        to_merge_txt.write("\n".join(finished_samples))
        to_merge_txt.close()
        print(f"{os.path.join(condor_jobdir, 'to_merge.txt')} written\n")
        #set_trace()

        ## this loop is responsible for merging or resubmitting jobs based on status
    #set_trace()
    for sample, status_list in samples_status.items():
        if len(status_list) == 1:
            #set_trace()
            status, sample_njobs = status_list[0].split("=")
            if (status == "RESUB") and (args.resub):
                print(f"\t{sample} needs to be resubmitted")
                os.system(f"cd {condor_jobdir}/{sample} && condor_submit condor.rescue.jdl && cd {orig_dir}")
                samples_status[sample] = [f"Idle={sample_njobs}"]
                njobs += int(sample_njobs)

        else:
            print(f"{sample} job status: {status_list}")

        ## recheck all jobs if there are none running
    if njobs == 0:
        failed_samples = []
        finished_samples = [] # reinitialize finished_samples to check if they're actually finished
        for sample in samples:
            isCorrect, nfails, alreadyComplete = check_correctness(sample=sample)
            if isCorrect == False: 
                failed_samples.append(sample)
                if args.resub: os.system(f"cd {condor_jobdir}/{sample} && condor_submit condor.rescue.jdl && cd {orig_dir}")
            else:
                finished_samples.append(sample)
        if len(failed_samples) == 0:
            print("All jobs completed")
            escape = True
        else:
            time.sleep(30)
    else:
        print(f"{njobs} jobs are submitted, checking again in 30 seconds\n")
        time.sleep(30)

## merge all TOT files from each sample into one
if sorted(set(tot_files.keys())) == sorted(set(samples)):
    print("All jobs completed. Run merge_jobs.py inside singularity to merge all jobs")
