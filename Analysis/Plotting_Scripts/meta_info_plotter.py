from coffea.hist import plot
import coffea
# matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from coffea.util import load, save
from pdb import set_trace
import os
import Utilities.prettyjson as prettyjson
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--sample', type=str, help='Input sample to use.')
parser.add_argument('--testing', action='store_true', help='Determines where input file is.')
parser.add_argument('--save_pu', action='store_true', help='Save PU distributions for each sample to a single file.')

args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'meta_info'

variables = {
    'mtt' : '$m_{t\\bart}$ [GeV]',
    'ctstar' : 'cos($\\theta_{t}^{*}$)'
}

input_dir = proj_dir if args.testing else '/'.join([proj_dir, 'results', jobid, analyzer])
f_ext = '%s.test.coffea' % analyzer if args.testing else '.coffea'
outdir = '/'.join([proj_dir, 'plots', jobid, analyzer, 'Test']) if args.testing else '/'.join([proj_dir, 'plots', jobid, analyzer])

if args.testing:
    fnames = ['%s/%s_%s' % (input_dir, args.sample, f_ext)] if args.sample else ['%s/%s' % (input_dir, f_ext)]
else:
    fnames = ['%s/%s%s' % (input_dir, args.sample, f_ext)] if args.sample else ['%s/%s' % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)]

#set_trace()

if not os.path.isdir(outdir):
    os.makedirs(outdir)

if args.save_pu:
    from coffea.lookup_tools.dense_lookup import dense_lookup
    pu_dists = {}

for fname in fnames:
    if not os.path.isfile(fname):
        raise IOError("%s not found" % fname)
    hists = load(fname)
    #set_trace()
    
    for hname in hists.keys():
            ## plot histograms
        if isinstance(hists[hname], coffea.hist.hist_tools.Hist):
            histo = hists[hname]

                ## save pileup array to dict
            if hname == 'PUDistribution' and args.save_pu:
                #set_trace()
                pu_dists[histo.axes()[0]._sorted[0]] = dense_lookup([val for val in histo.values().values()][0], histo.axes()[-1].edges())
    
            if histo.dense_dim() == 1:
                ## make plot for separate samples
                for sample in histo.axes()[0]._sorted:
                    sample_histo = histo[sample]
    
                    fig = plt.figure()
                    ax = plot.plot1d(sample_histo)
                    if 'mtt' in hname:
                        plt.xlabel(variables['mtt'])
                    elif 'ctstar' in hname:
                        plt.xlabel(variables['ctstar'])
                    else:
                        plt.xlabel('$%s$' % sample_histo.axes()[-1].label)
                    plt.ylabel('Events')
                    figname = '%s/%s_%s.png' % (outdir, sample, hname)
                    fig.savefig(figname)
                    print('%s written' % figname)
                    plt.close()
    
    
            elif histo.dense_dim() == 2:
                xvar, yvar = histo.axes()[-2].name, histo.axes()[-1].name
    
                ## make plots for different ttJets samples
                for sample in histo.axes()[0]._sorted:
                    if not (sample == 'ttJets' or sample == 'ttJets_PS'): continue
                    sample_histo = histo[sample].project(histo.axes()[1].name, xvar, yvar)

                    #set_trace()    
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
            #set_trace()
            #if args.testing: set_trace()
            meta_dict = {}
            for key, val in hists[hname].items():
                if isinstance(val, (int, float, list)):
                    meta_dict[key] = val
                else:
                    meta_dict[key] = val.tolist()
            with open('%s/inputs/%s/%s.meta.json' % (proj_dir, jobid, hname), 'w') as out:
                out.write(prettyjson.dumps(meta_dict))
                print('%s/inputs/%s/%s.meta.json written' % (proj_dir, jobid, hname))
    

if args.save_pu:
    #set_trace()
    pu_name = '%s/Corrections/MC_PUDistributions.coffea' % proj_dir
    save(pu_dists, pu_name)
    print('\n', pu_name, 'written')
