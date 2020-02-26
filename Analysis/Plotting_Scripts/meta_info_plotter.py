from coffea.hist import plot
import coffea
# matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from coffea.util import load
from pdb import set_trace
import os
import Utilities.prettyjson as prettyjson
from argparse import ArgumentParser

parser = ArgumentParser()
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
        ## plot histograms
    if isinstance(hists[hname], coffea.hist.hist_tools.Hist):
        histo = hists[hname]
        #set_trace()

        if histo.dense_dim() == 1:
            fig = plt.figure()
            ax = plot.plot1d(histo)
            if 'mtt' in hname:
                plt.xlabel(variables['mtt'])
            elif 'ctstar' in hname:
                plt.xlabel(variables['ctstar'])
            else:
                plt.xlabel('$%s$' % histo.axes()[-1].label)
            plt.ylabel('Events')
            figname = '%s/%s.png' % (outdir, hname)
            fig.savefig(figname)
            print('%s written' % figname)
            plt.close()

        elif histo.dense_dim() == 2:
            xvar, yvar = histo.axes()[-2].name, histo.axes()[-1].name

            ## make plots for different ttJets samples
            for sample in histo.axes()[0]._sorted:
                sample_histo = histo[sample].project(histo.axes()[1].name, xvar, yvar)

                    ## plot x projection
                fig = plt.figure()
                x_proj_histo = sample_histo.sum(yvar)
                x_ax = plot.plot1d(x_proj_histo)
                x_ax.set_xlabel(variables[xvar])
                x_ax.set_ylabel('Events')
                xfigname = '%s/%s_%s.png' % (outdir, sample, xvar)
                fig.savefig(xfigname)
                print('%s written' % xfigname)
                plt.close()

                    ## plot y projection
                fig = plt.figure()
                y_proj_histo = sample_histo.sum(xvar)
                y_ax = plot.plot1d(y_proj_histo)
                y_ax.set_xlabel(variables[yvar])
                y_ax.set_ylabel('Events')
                yfigname = '%s/%s_%s.png' % (outdir, sample, yvar)
                fig.savefig(yfigname)
                print('%s written' % yfigname)
                plt.close()

        ## write other meta info to meta.json files
    else:
        if args.testing: continue
        meta_dict = {}
        for key, val in hists[hname].items():
            if isinstance(val, (int, float, list)):
                meta_dict[key] = val
            else:
                meta_dict[key] = val.tolist()
        set_trace()
        with open('%s/inputs/%s/%s.meta.json' % (proj_dir, jobid, hname), 'w') as out:
            out.write(prettyjson.dumps(meta_dict))
            print('%s/inputs/%s/%s.meta.json written' % (proj_dir, jobid, hname))


