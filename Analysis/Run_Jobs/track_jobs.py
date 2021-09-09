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
jobdir = args.jobdir if (args.jobdir).startswith(os.path.join('afs', 'cern.ch', 'work', 'j', user)) else os.path.join(proj_dir, args.jobdir)
jobdir = jobdir[:-1] if jobdir.endswith('/') else jobdir
if not os.path.isdir(jobdir):
    raise IOError('Directory: %s does not exist' % jobdir)

if os.path.isfile(os.path.join(jobdir, '%s_TOT.coffea' % args.jobdir)):
    print('Combined output file %s/%s_TOT.coffea already exists' % (jobdir, args.jobdir))
    import sys; sys.exit()

orig_dir = os.getcwd()
samples = [dirname for dirname in os.listdir(jobdir) if os.path.isdir(os.path.join(jobdir, dirname))]

def check_correctness(jobdir, dump_rescue = False):
    isCorrect = True
    sample = jobdir.split('/')[-1]
    if os.path.isfile('%s/%s_TOT.coffea' % (jobdir, sample)):  ## TOT output file already created
        #print('%s/%s_TOT.coffea already exists' % (jobdir, sample))
        return isCorrect, 0, True
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

    return isCorrect, len(fails), False

tot_files = {}
def merge_files(jobdir, samples):
    #set_trace()
    for sample in samples:
        merged_fname = os.path.join(jobdir, sample, '%s_TOT.coffea' % sample)
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

escape = False
finished_samples = []

while not escape:
    stdout= subprocess.Popen("condor_q $USER -nobatch", stdout=subprocess.PIPE, shell=True).stdout.read()
    njobs = 0
    samples_status = {}

        ## this loop is just to check how many jobs are still running for each sample
    for sample in samples:
        if sample in finished_samples: continue
            ## checks how many lines correspond to jobdir/sample 
        sample_njobs = len([line for line in stdout.split(b'\n') if (os.path.join(jobdir, '%s/' % sample)).encode() in line])
        if sample_njobs == 0:
            isCorrect, nfails, alreadyComplete = check_correctness(os.path.join(jobdir, sample), dump_rescue=True)
            if isCorrect == True:
                status = 'MERGED' if alreadyComplete else 'COMPLETE' # jobs are already merged or (completed but still need to be merged)
            else:
                status = 'RESUB' # jobs need to be resubmitted
            samples_status[sample] = (status, nfails)
        else:
                ## get lines corresponding to sample
            job_lines = [line.decode() for line in stdout.split(b'\n') if (os.path.join(jobdir, '%s/' % sample)).encode() in line]
                ## get job statuses
            job_statuses = [''.join(line).split()[5] for line in job_lines] # 5 is hardcoded
                ## check how many jobs haven't been removed
            njobs_running = len([status for status in job_statuses if not (status == 'X')])
            if njobs_running == 0:
                isCorrect, nfails, alreadyComplete = check_correctness(os.path.join(jobdir, sample), dump_rescue=True)
                if isCorrect == True:
                    status = 'MERGED' if alreadyComplete else 'COMPLETE' # jobs are already merged or (completed but still need to be merged)
                else:
                    status = 'RESUB' # jobs need to be resubmitted
                samples_status[sample] = (status, nfails)
            else:
                samples_status[sample] = ('RUNNING', njobs_running) # jobs running
                njobs += njobs_running

        ## this loop is responsible for merging or resubmitting jobs based on status
    for sample, (status, sample_njobs) in samples_status.items():
            ## if status is 'RESUB', resubmit jobs to condor
        if status == 'RESUB':
            #set_trace()
            print(f"\t{sample} needs to be resubmitted")
            os.system('cd %s/%s && condor_submit condor.rescue.jdl && cd %s' % (jobdir, sample, orig_dir))
            samples_status[sample] = ('RUNNING', sample_njobs)
            njobs += sample_njobs

            ## if status is 'RUNNING' do nothing, jobs are already running
        elif status == 'RUNNING':
            print('\t%s has %i jobs running' % (sample, sample_njobs) )

            ## if status is 'MERGED' or 'COMPLETE', run merge_jobs and add files to finished_samples list so they're not check anymore
        else:
            merge_files(jobdir, [sample])
            finished_samples.append(sample)

        ## recheck all jobs if there are none running
    if njobs == 0:
        failed_samples = []
        finished_samples = [] # reinitialize finished_samples to check if they're actually finished
        for sample in samples:
            isCorrect, nfails, alreadyComplete = check_correctness('%s/%s' % (jobdir, sample), dump_rescue=True)
            if isCorrect == False: 
                failed_samples.append(sample)
                os.system('cd %s/%s && condor_submit condor.rescue.jdl && cd %s' % (jobdir, sample, orig_dir))
            else:
                finished_samples.append(sample)
        if finished_samples:
            merge_files(jobdir, finished_samples)
        if len(failed_samples) == 0:
            print('All jobs completed')
            escape = True
        else:
            time.sleep(30)
    else:
        print("%i jobs are still running, checking again in 30 seconds\n" % njobs )
        time.sleep(30)

## merge all TOT files from each sample into one
if list(set(tot_files.keys())) == list(set(samples)):
    tot_acc = plot_tools.add_coffea_files(list(tot_files.values()))
    tot_outname = '%s_TOT.coffea' % args.jobdir.strip('/')
    plot_tools.save_accumulator(tot_acc, '%s/%s' % (jobdir, tot_outname))

