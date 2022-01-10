import os, subprocess, time
from pdb import set_trace
import Utilities.plot_tools as plot_tools

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("jobdir", help="specify name of directory where jobs were submitted")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
user = os.environ["USER"]

## once jobs are finished, merge all output coffea files for each sample into one
jobdir = args.jobdir if (args.jobdir).startswith(os.path.join("afs", "cern.ch", "work", "j", user)) else os.path.join(proj_dir, args.jobdir)
jobdir = jobdir[:-1] if jobdir.endswith("/") else jobdir
if not os.path.isdir(jobdir):
    raise IOError("Directory: %s does not exist" % jobdir)



if os.path.isfile(os.path.join(jobdir, f"{args.jobdir}_TOT.coffea")):
    print(f"Combined output file {jobdir}/{args.jobdir}_TOT.coffea already exists")
    import sys; sys.exit()

orig_dir = os.getcwd()
all_possible_samples = sorted([dirname for dirname in os.listdir(jobdir) if os.path.isdir(os.path.join(jobdir, dirname))])

tot_files = {}
def merge_files(jobdir, samples):
    #set_trace()
    for sample in samples:
        merged_fname = os.path.join(jobdir, sample, f"{sample}_TOT.coffea")
        output_files = [f"{jobdir}/{sample}/{fname}" for fname in os.listdir(f"{jobdir}/{sample}") if fname.endswith(".coffea")]
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
if not os.path.isfile(os.path.join(jobdir, "to_merge.txt")):
    raise ValueError(f"No file found to merge in {jobdir}")

txt_file = open(os.path.join(jobdir, "to_merge.txt"), "r")
samples = [sample.strip("\n") for sample in txt_file]
#set_trace()
for sample in samples:
    merge_files(jobdir, [sample])

## merge all TOT files from each sample into one
if sorted(set(tot_files.keys())) == sorted(set(all_possible_samples)):
    tot_acc = plot_tools.add_coffea_files(list(tot_files.values()))
    tot_outname = "%s_TOT.coffea" % args.jobdir.strip("/")
    plot_tools.save_accumulator(tot_acc, f"{jobdir}/{tot_outname}")

