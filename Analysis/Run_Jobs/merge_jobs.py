import os, subprocess, time
from pdb import set_trace
import Utilities.plot_tools as plot_tools

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("jobdir", help="specify name of directory where jobs were submitted")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
user = os.environ["USER"]

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

orig_dir = os.getcwd()
all_possible_samples = sorted([dirname for dirname in os.listdir(condor_jobdir) if os.path.isdir(os.path.join(condor_jobdir, dirname))])

tot_files = {}
def merge_files(samples):
    #set_trace()
    for sample in samples:
        print(f"Merging {sample}")
        merged_fname = os.path.join(eos_jobdir, sample, f"{sample}_TOT.coffea")
        output_files = [f"{eos_jobdir}/{sample}/{fname}" for fname in os.listdir(f"{eos_jobdir}/{sample}") if fname.endswith(".coffea")]
            # don't re-merge files if this has already been done
        if merged_fname in output_files:
            tot_files[sample] = merged_fname
            print(f"{merged_fname} already exists")
            continue
        else:
            if len(output_files) > 1:
                ## merge files
                output_acc = plot_tools.add_coffea_files(output_files)
                plot_tools.save_accumulator(output_acc, merged_fname)
            else:
                ## rename file 
                os.system(f"mv {output_files[0]} {merged_fname}")
                print(f"{merged_fname} written")
            tot_files[sample] = merged_fname

    # open to_merge.txt
if not os.path.isfile(os.path.join(condor_jobdir, "to_merge.txt")):
    raise ValueError(f"No file found to merge in {condor_jobdir}")

#set_trace()
txt_file = open(os.path.join(condor_jobdir, "to_merge.txt"), "r")
samples_to_merge = [sample.strip("\n") for sample in txt_file]
merge_files(samples_to_merge)

## merge all TOT files from each sample into one
if sorted(set(tot_files.keys())) == sorted(set(all_possible_samples)):
    #set_trace()
    tot_acc = plot_tools.add_coffea_files(list(tot_files.values()))
    tot_outname = "%s_TOT.coffea" % args.jobdir.strip("/")
    plot_tools.save_accumulator(tot_acc, f"{eos_jobdir}/{tot_outname}")
