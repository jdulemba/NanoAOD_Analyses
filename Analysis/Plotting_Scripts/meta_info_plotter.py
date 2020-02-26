from coffea.hist import plot
import coffea
# matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from coffea.util import load
from pdb import set_trace
import os
from argparse import ArgumentParser

parser = ArgumentParser()
#parser.add_argument('sample', default='ttJets', help='Samples to run over')
parser.add_argument('--testing', action='store_true', help='Determines where input file is.')

args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'meta_info'

if not args.testing:
    raise IOError("only testing supported for plotting at the moment")

fname = '%s/%s.test.coffea' % (proj_dir, analyzer) if args.testing else '%s/%s.coffea' % ('/'.join([proj_dir, 'results', jobid]), 'test')
#fname = '%s/%s.test.%s.coffea' % (proj_dir, args.sample, analyzer) if args.testing else '%s/%s.coffea' % ('/'.join([proj_dir, 'results', jobid]), args.sample)
outdir = '/'.join([proj_dir, 'plots', jobid, analyzer, 'Test']) if args.testing else '/'.join([proj_dir, 'plots', jobid, analyzer])

if not os.path.isdir(outdir):
    os.makedirs(outdir)

hists = load(fname)
#set_trace()

variables = {
    'mtt' : '$m_{t\\bart}$ [GeV]',
    'ctstar' : 'cos($\\theta_{t}^{*}$)'
}

for hname in hists.keys():
    if not isinstance(hists[hname], coffea.hist.hist_tools.Hist): continue
    #if hname == 'MetaInfo': continue
    histo = hists[hname]
    set_trace()

    if histo.dense_dim() == 1:
        fig = plt.figure()
        ax = plot.plot1d(histo)
        if 'mtt' in hname:
            plt.xlabel(variables['mtt'])
        elif 'ctstar' in hname:
            plt.xlabel(variables['ctstar'])
        else:
            plt.xlabel('$%s$' % histo.axes()[0].label)
        plt.ylabel('Events')
        figname = '%s/%s.png' % (outdir, hname)
        fig.savefig(figname)
        print('%s written' % figname)
        plt.close()

    elif histo.dense_dim() == 2:
        xvar, yvar = histo.axes()[-2].name, histo.axes()[-1].name

        ## make plots for different ttJets samples
        for sample in histo.axes()[0]._sorted:
            

            ## plot x projection
        fig = plt.figure()
        x_proj_histo = histo.sum(yvar)
        x_ax = plot.plot1d(x_proj_histo)
        x_ax.set_xlabel(variables[xvar])
        x_ax.set_ylabel('Events')
        xfigname = '%s/%s.png' % (outdir, xvar)
        fig.savefig(xfigname)
        print('%s written' % xfigname)
        plt.close()

            ## plot y projection
        fig = plt.figure()
        y_proj_histo = histo.sum(xvar)
        y_ax = plot.plot1d(y_proj_histo)
        y_ax.set_xlabel(variables[yvar])
        y_ax.set_ylabel('Events')
        yfigname = '%s/%s.png' % (outdir, yvar)
        fig.savefig(yfigname)
        print('%s written' % yfigname)
        plt.close()


