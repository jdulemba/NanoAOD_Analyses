import os, subprocess, time
from pdb import set_trace
import Utilities.plot_tools as plot_tools

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('jobdir', help='specify name of directory where jobs were submitted')
args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
user = os.environ['USER']

#set_trace()
escape = False
start = time.time()
totjobs = -1

while not escape:
    stdout= subprocess.Popen("condor_q $USER", stdout=subprocess.PIPE, shell=True).stdout.read()
    #set_trace()
    njobs = int(stdout.split(b'Total for query:')[-1].split(b'jobs;')[0])

    eta = 'UNKNOWN'
    if totjobs != -1:
        elapsed = time.time() - start
        completed = totjobs - njobs
        if completed != 0:
            eta_time = njobs*(elapsed/completed)
            m, s = divmod(eta_time, 60)
            h, m = divmod(m, 60)
            eta = "%d:%02d:%02d" % (h, m, s)
    else:
        totjobs = njobs

    if njobs == 0:
        print('Jobs completed: runtime = %i' % elapsed)
        escape = True

    else:
        print("%i jobs are still running, checking again in 30 seconds. ETA: %s" % (njobs, eta))
        time.sleep(30)


## once jobs are finished, merge all output coffea files for each sample into one
jobdir = args.jobdir if (args.jobdir).startswith('/afs/cern.ch/work/j/%s/' % user) else '%s/%s' % (proj_dir, args.jobdir)
if not os.path.isdir(jobdir):
    raise IOError('Directory: %s does not exist' % jobdir)

set_trace()
samples = [dirname for dirname in os.listdir(jobdir) if os.path.isdir('%s/%s' % (jobdir, dirname))]

tot_files = []
for sample in samples:
    output_files = ['%s/%s/%s' % (jobdir, sample, fname) for fname in os.listdir('%s/%s' % (jobdir, sample)) if fname.endswith('.coffea')]
    if len(output_files) > 1:
        ## merge files
        output_acc = plot_tools.add_coffea_files(output_files)
        plot_tools.save_accumulator(output_acc, '%s/%s/%s_TOT.coffea' % (jobdir, sample, sample))
    else:
        ## mv file to
        os.system('mv %s %s/%s/%s_TOT.coffea' % (output_files[0], jobdir, sample, sample))
    tot_files.append('%s/%s/%s_TOT.coffea' % (jobdir, sample, sample))

## merge all TOT files from each sample into one
tot_acc = plot_tools.add_coffea_files(tot_files)
plot_tools.save_accumulator(tot_acc, '%s/total_output.coffea' % jobdir)
