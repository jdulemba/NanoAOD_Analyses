from coffea.hist import plot
# matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from coffea.util import load
from pdb import set_trace
import os
import styles
import Utilities.plot_tools as plt_tools
import Utilities.prettyjson as prettyjson
import re
from coffea import hist

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('year', choices=['2016', '2017', '2018'], help='What year is the ntuple from.')
parser.add_argument('--testing', action='store_true', help='Determines where input file is.')

args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'presel_analyzer'

if not args.testing:
    raise IOError("only testing supported for plotting at the moment")

fname = '%s/%s.test.coffea' % (proj_dir, analyzer) if args.testing else '%s/%s.coffea' % ('/'.join([proj_dir, 'results', jobid]), 'test')
#fname = '%s/%s.test.%s.coffea' % (proj_dir, args.sample, analyzer) if args.testing else '%s/%s.coffea' % ('/'.join([proj_dir, 'results', jobid]), args.sample)
outdir = '/'.join([proj_dir, 'plots', jobid, analyzer, 'Test']) if args.testing else '/'.join([proj_dir, 'plots', jobid, analyzer])

if not os.path.isdir(outdir):
    os.makedirs(outdir)

data_lumi_year = prettyjson.loads(open('%s/inputs/lumis_data.json' % proj_dir).read())[args.year]

## make data and mc categories for data/MC plotting
mc_samples = re.compile('(?!data*)')
data_samples = re.compile('(data*)')


hists = load(fname)
#set_trace()

variables = {
    'pt' : '$p_{T}$',
    'eta' : '$\\eta$',
    'phi' : '$\\phi$',
    'energy' : 'E',
    'njets' : '$n_{jets}$',
}


jet_mults = {
    '3Jets' : '3 jets',
    '4PJets' : '4+ jets'
}

objtypes = {
    'Jets' : 'jets',
    'Lep' :  {
        'Muon' : '$\\mu$',
        'Electron' : '$e$',
    }
}

    ## get plotting colors/settings
hstyles = styles.styles
stack_fill_opts = {'alpha': 0.8, 'edgecolor':(0,0,0,.5)}
stack_error_opts = {'edgecolor':(0,0,0,.5)}
error_opts = {
    'label':'Stat. Unc.',
    'hatch':'///',
    'facecolor':'none',
    'edgecolor':(0,0,0,.5),
    'linewidth': 0
}
data_err_opts = {
    #'linestyle':'none',
    'marker': '.',
    'markersize': 10.,
    'color':'k',
    'elinewidth': 1,
}

#set_trace()
process = hist.Cat("process", "Process", sorting='placement')
process_cat = "dataset"
for hname in hists.keys():
    if hname == 'cutflow': continue
    hists[hname] = hists[hname].group(process_cat, process, plt_tools.hardcoded_groups)
    

for hname in hists.keys():
    if hname == 'cutflow': continue
    histo = hists[hname]
    #set_trace()

    if histo.dense_dim() == 1:
        ## hists should have 3 category axes (dataset, jet multiplicity, lepton type) followed by variable
        for jmult in histo.axes()[1]._sorted:
            for lep in histo.axes()[2]._sorted:
                fig, (ax, rax) = plt.subplots(2, 1, figsize=(7,7), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                fig.subplots_adjust(hspace=.07)
                hslice = histo[:, jmult, lep].integrate('jmult').integrate('leptype')

                    ## plot MC and data
                plot.plot1d(hslice[mc_samples],
                    overlay=hslice.axes()[0].name,
                    ax=ax,
                    clear=False,
                    stack=True,
                    line_opts=None,
                    fill_opts=stack_fill_opts,
                    error_opts=stack_error_opts
                )
                plot.plot1d(hslice[data_samples],
                    overlay=hslice.axes()[0].name,
                    ax=ax,
                    clear=False,
                    error_opts=data_err_opts
                )
                ax.autoscale(axis='x', tight=True)
                ax.set_ylim(0, None)
                ax.set_xlabel(None)

                    ## set legend and corresponding colors
                handles, labels = ax.get_legend_handles_labels()
                for idx, sample in enumerate(labels):
                    if sample == 'data' or sample == 'Observed': continue
                    #print(sample)
                    #set_trace()
                    facecolor, legname = plt_tools.get_styles(sample, hstyles)
                    handles[idx].set_facecolor(facecolor)
                    labels[idx] = legname
                # call plt.legend() with the new values
                #set_trace()
                ax.legend(handles,labels)

                    ## plot data/MC ratio
                plot.plotratio(hslice[data_samples].sum(hslice.axes()[0].name), hslice[mc_samples].sum(hslice.axes()[0].name), 
                    ax=rax,
                    error_opts=data_err_opts, 
                    denom_fill_opts={},
                    guide_opts={},
                    unc='num'
                )
                rax.set_ylabel('data/MC')
                rax.set_ylim(0,2)


                    ## set axes labels and titles
                obj, kvar = hname.split('_')
                if kvar in variables.keys():
                    if kvar == 'njets':
                        xtitle = variables[kvar]
                    elif kvar == 'pt' or kvar == 'energy':
                        xtitle = '%s(%s) [GeV]' % (variables[kvar], objtypes[obj]) if obj == 'Jets' else '%s(%s) [GeV]' % (variables[kvar], objtypes[obj][lep])
                    else:
                        xtitle = '%s(%s)' % (variables[kvar], objtypes[obj]) if obj == 'Jets' else '%s(%s)' % (variables[kvar], objtypes[obj][lep])
                    plt.xlabel(xtitle)
                else:
                    plt.xlabel('$%s$' % histo.axes()[-1].label)

                cms_blurb = plt.text(
                    0., 1., r"CMS Preliminary",
                    fontsize=12, 
                    horizontalalignment='left', 
                    verticalalignment='bottom', 
                    transform=ax.transAxes,
                    style='italic'
                )
                lumi_blurb = plt.text(
                    1., 1., r"(13 TeV %.2f fb$^{-1}$, %s + %s)" % (data_lumi_year['%ss' % lep]/1000., jet_mults[jmult], objtypes['Lep'][lep]),
                    fontsize=12, 
                    horizontalalignment='right', 
                    verticalalignment='bottom', 
                    transform=ax.transAxes
                )

                pltdir = outdir if args.testing else '/'.join([outdir, jmult, lep_type])
                if not os.path.isdir(pltdir):
                    os.makedirs(pltdir)
                figname = '%s/%s.png' % (pltdir, '_'.join([jmult, lep, hname]))
                #figname = '%s/%s_%s.png' % (pltdir, args.sample, hname)

                fig.savefig(figname)
                print('%s written' % figname)
                plt.close()
