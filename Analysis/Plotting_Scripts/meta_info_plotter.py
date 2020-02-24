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

if not args.debug:
    raise IOError("only debug supported for plotting at the moment")

fname = '%s/%s.test.%s.coffea' % (proj_dir, args.sample, analyzer) if args.debug else '%s/%s.coffea' % ('/'.join([proj_dir, 'results', jobid]), args.sample)
outdir = '/'.join([proj_dir, 'plots', jobid, analyzer, 'Test']) if args.debug else '/'.join([proj_dir, 'plots', jobid, analyzer])

if not os.path.isdir(outdir):
    os.makedirs(outdir)

hists = load(fname)
#set_trace()

variables = {
    'mtt' : '$m_{t\\bart}$',
    'ctstar' : 'cos($\\theta_{t}^{*}$)'
}

#set_trace()
for hname in hists.keys():
    if hname == 'MetaInfo': continue
    histo = hists[hname]

    fig = plt.figure()
    #set_trace()
    if histo.dense_dim() == 1:
        #continue
        ax = plot.plot1d(histo)
        if 'mtt' in hname:
            plt.xlabel(variables['mtt'])
        elif 'ctstar' in hname:
            plt.xlabel(variables['ctstar'])
        else:
            plt.xlabel('$%s$' % histo.axes()[0].label)
        plt.ylabel('Events')
    elif histo.dense_dim() == 2:
        continue
        set_trace()
        ax = plot.plot2d(histo, xaxis='mtt')
        plt.xlabel(variables['mtt'])
        plt.ylabel(variables['ctstar'])


    figname = '%s/%s_%s.png' % (outdir, args.sample, hname)
    fig.savefig(figname)
    print('%s written' % figname)
    plt.close()
