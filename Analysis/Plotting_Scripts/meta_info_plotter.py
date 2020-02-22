from coffea.hist import plot
# matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from coffea.util import load
from pdb import set_trace
import os
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('sample', default='ttJets', help='Samples to run over')
parser.add_argument('--debug', action='store_true', help='Determines where input file is.')

args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'get_meta_info'

#set_trace()
fname = '%s/%s.test.%s.coffea' % (proj_dir, args.sample, analyzer) if args.debug else '%s/%s.coffea' % ('/'.join([proj_dir, 'results', jobid]), args.sample)
outdir = '/'.join([proj_dir, 'results'])
if not os.path.isdir(outdir):
    os.makedirs(outdir)

hists = load(fname)

for hname in hists.keys():
    histo = hists[hname]

    fig = plt.figure()
    if histo.dim() == 1:
        ax = plot.plot1d(histo)
        plt.xlabel('$%s$' % histo.axes()[0].label)
        plt.ylabel('Events')
    elif histo.dim() == 2:
        set_trace()
        ax = plot.plot2d(histo, 'mtt')
        plt.xlabel('$%s$' % histo.axes()[0].label)
        plt.ylabel(histo.axes()[1].label)


    figname = '%s/%s_%s.png' % (outdir, args.sample, hname)
    fig.savefig(figname)
    print('%s written' % figname)
    plt.close()
