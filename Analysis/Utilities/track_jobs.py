import os, subprocess, time
from pdb import set_trace
import Utilities.plot_tools as plot_tools

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('jobdir', help='specify name of directory where jobs were submitted')
args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
user = os.environ['USER']

## once jobs are finished, merge all output coffea files for each sample into one
jobdir = args.jobdir if (args.jobdir).startswith('/afs/cern.ch/work/j/%s/' % user) else '%s/%s' % (proj_dir, args.jobdir)
if not os.path.isdir(jobdir):
    raise IOError('Directory: %s does not exist' % jobdir)

if os.path.isfile('%s/%s_TOT.coffea' % (jobdir, args.jobdir)):
    print('Combined output file %s/%s_TOT.coffea already exists' % (jobdir, args.jobdir))
    import sys; sys.exit()

samples = [dirname for dirname in os.listdir(jobdir) if os.path.isdir('%s/%s' % (jobdir, dirname))]

def check_correctness(jobdir, dump_rescue = False):
    isCorrect = True
    sample = jobdir.split('/')[-1]
    #set_trace()
    if os.path.isfile('%s/%s_TOT.coffea' % (jobdir, sample)):  ## TOT output file already created
        #print('%s/%s_TOT.coffea already exists' % (jobdir, sample))
        return isCorrect
    jdl = open('%s/condor.jdl' % jobdir).read()
    blocks = jdl.split('\n\n')
    header = blocks[0]
    block_map = {}
    for block in blocks[1:]:
            ## save output filename to check if it exists
        fname = block.split('Arguments = ')[1].split(' ')[-1].split('=')[-1].split('\n')[0]
        block_map[fname] = block

    npass = 0
    fails = []
    for fname, arguments in block_map.items():
        if not os.path.isfile(fname):
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
            print('dumping rescue job')
            with open('%s/condor.rescue.jdl' % jobdir, 'w') as rescue:
                rescue.write(header)
                rescue.write('\n\n')
                rescue.write(
                    '\n\n'.join([block_map[key] for key in fails])
                    )

    return isCorrect

tot_files = {}
def merge_files(jobdir, samples):
    #set_trace()
    for sample in samples:
        merged_fname = '%s/%s/%s_TOT.coffea' % (jobdir, sample, sample)
        output_files = ['%s/%s/%s' % (jobdir, sample, fname) for fname in os.listdir('%s/%s' % (jobdir, sample)) if fname.endswith('.coffea')]
            # don't re-merge files if this has already been done
        if merged_fname in output_files:
            tot_files[sample] = merged_fname
            print('%s already exists' % merged_fname)
            continue
        else:
            if len(output_files) > 1:
                ## merge files
                output_acc = plot_tools.add_coffea_files(output_files)
                plot_tools.save_accumulator(output_acc, merged_fname)
            else:
                ## rename file 
                os.system('mv %s %s' % (output_files[0], merged_fname))
                print('%s written' % merged_fname)
            tot_files[sample] = merged_fname

#set_trace()
escape = False
start = time.time()
elapsed = 0
totjobs = -1
finished_samples = []

while not escape:
    stdout= subprocess.Popen("condor_q $USER -nobatch", stdout=subprocess.PIPE, shell=True).stdout.read()
    njobs = 0
    for sample in samples:
        if sample in finished_samples: continue
            ## checks how many lines correspond to jobdir/sample 
        sample_njobs = len([line for line in stdout.split(b'\n') if ('%s/%s' % (jobdir, sample)).encode() in line])
        if sample_njobs == 0:
            finished_samples.append(sample)
        else:
            njobs += sample_njobs

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
            ## check if all jobs have run correctly
        failed_samples = []
        finished_samples = [] # reinitialize finished_samples to check if they're actually finished
        for sample in samples:
            isCorrect = check_correctness('%s/%s' % (jobdir, sample), dump_rescue=True)
            if isCorrect == False: 
                failed_samples.append(sample)
                os.system('condor_submit %s/%s/condor.rescue.jdl' % (jobdir, sample))
            else:
                finished_samples.append(sample)
        if finished_samples:
            merge_files(jobdir, finished_samples)
        if len(failed_samples) == 0:
            print('All jobs completed: runtime = %i' % elapsed)
            escape = True
        else:
            time.sleep(30)

    else:
        print("%i jobs are still running, checking again in 30 seconds. ETA: %s" % (njobs, eta))
        time.sleep(30)

## merge all TOT files from each sample into one
if list(set(tot_files.keys())) == list(set(samples)):
    tot_acc = plot_tools.add_coffea_files(list(tot_files.values()))
    plot_tools.save_accumulator(tot_acc, '%s/%s_TOT.coffea' % (jobdir, args.jobdir))

