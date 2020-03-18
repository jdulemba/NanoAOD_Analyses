import os, subprocess, time
from pdb import set_trace

#from argparse import ArgumentParser
#parser = ArgumentParser()
#parser.add_argument('indir', help='specify name of input directory')
#
#args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
user = os.environ['USER']

##set_trace()
#jobdir = args.indir if (args.indir).startswith('/afs/cern.ch/work/j/%s/' % user) else '%s/%s' % (proj_dir, args.indir)
#if not os.path.isdir(jobdir):
#    raise IOError('Directory: %s does not exist' % jobdir)

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
        print('Jobs completed')
        escape = True

    else:
        print("%i jobs are still running, checking again in 30 seconds. ETA: %s" % (njobs, eta))
        time.sleep(30)

