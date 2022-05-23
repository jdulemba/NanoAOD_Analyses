import os
from pdb import set_trace

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("jobdir", help="specify name of directory where jobs were submitted")
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
    print(f"Combined output file {jobdir}/{args.jobdir}_TOT.coffea already exists")
    import sys; sys.exit()

    # find files that are already finished
finished_jobs_fname = os.path.join(condor_jobdir, "to_merge.txt")
if not os.path.isfile(finished_jobs_fname): raise ValueError(f"No file {finished_jobs_fname} with finished jobs found.")
to_merge_txt = open(finished_jobs_fname).read()
finished_samples = to_merge_txt.split("\n")

for sample in finished_samples:
    input_files = [fname for fname in os.listdir(os.path.join(eos_jobdir, sample))]
    if f"{sample}_TOT.coffea" in input_files:
        files_to_remove = [os.path.join(eos_jobdir, sample, fname) for fname in input_files if fname != f"{sample}_TOT.coffea"]
        if files_to_remove:
            cmd = f"rm {' '.join(files_to_remove)}"
            print(f"{sample} Executing:\n\t{cmd}\n")
            os.system(cmd)
        else:
            print(f"{sample}_TOT.coffea file already the only one left.\n")
    else:
        print(f"{sample}:\n\t{sample}_TOT.coffea file not found\n")
