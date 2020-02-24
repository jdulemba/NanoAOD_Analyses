from coffea.hist import plot
# matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from coffea.util import load
from pdb import set_trace
import os
from argparse import ArgumentParser
import styles
import fnmatch


parser = ArgumentParser()
parser.add_argument('sample', default='ttJets', help='Samples to run over')
parser.add_argument('--year', default='2016', help='What year is the ntuple from.')
parser.add_argument('--debug', action='store_true', help='Determines where input file is.')

args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'test_analyzer'

if not args.debug:
    raise IOError("only debug supported for plotting at the moment")

fname = '%s/%s.test.%s.coffea' % (proj_dir, args.sample, analyzer) if args.debug else '%s/%s.coffea' % ('/'.join([proj_dir, 'results', jobid]), args.sample)
outdir = '/'.join([proj_dir, 'plots', jobid, analyzer, 'Test']) if args.debug else '/'.join([proj_dir, 'plots', jobid, analyzer])

if not os.path.isdir(outdir):
    os.makedirs(outdir)

hists = load(fname)

variables = {
    'pt' : '$p_{T}$',
    'eta' : '$\\eta$',
    'phi' : '$\\phi$',
    'njets' : '$n_{jets}$',
}

jet_mults = {
    '3Jets' : '3 jets',
    '4PJets' : '4+ jets'
}

lep_types = {
    'TIGHTMU' : 'tight $\\mu$',
    'LOOSEMU' : 'loose $\\mu$',
    'TIGHTEL' : 'tight $e$',
    'LOOSEEL' : 'loose $e$',
}

objtypes = {
    'Jets' : 'jets',
    'Muon' : '$\\mu$',
    'Electron' : '$e$'
}

    ## get plotting colors/settings
hstyles = styles.styles
stack_fill_opts = {'alpha': 0.8, 'edgecolor':(0,0,0,.5)}
#stack_error_opts = {'label':'Stat. Unc.', 'hatch':'///', 'facecolor':'none', 'edgecolor':(0,0,0,.5), 'linewidth': 0}

def get_styles(sample, styles):
    best_pattern = ''
    for pattern, style_dict in styles.items():
        if fnmatch.fnmatch(sample, pattern):
            if len(pattern) > len(best_pattern):
                best_pattern = pattern
    if best_pattern:
        return styles[best_pattern]['facecolor'], styles[best_pattern]['name']
    else:
        return None


for hname in hists.keys():
    if hname == 'cutflow': continue
    histo = hists[hname]

    jmult, lep_type, obj_type, kvar = hname.split('_')

    fig, ax = plt.subplots(1, 1, figsize=(7,7))
    fig.subplots_adjust(hspace=.07)
    if histo.dense_dim() == 1:
        #set_trace()
        plot.plot1d(histo, ax=ax, stack=True, fill_opts=stack_fill_opts)
        if kvar in variables.keys():
            xtitle = '%s(%s) [GeV]' % (variables[kvar], objtypes[obj_type]) if kvar == 'pt' else '%s(%s)' % (variables[kvar], objtypes[obj_type])
            plt.xlabel(xtitle)
        else:
            plt.xlabel('$%s$' % histo.axes()[-1].label)
        plt.ylabel('Events')

    ax.autoscale(axis='x', tight=True)

        ## set legend and corresponding colors
    handles, labels = plt.gca().get_legend_handles_labels()
    for idx, sample in enumerate(labels):
        facecolor, legname = get_styles(sample, hstyles)
        handles[idx].set_facecolor(facecolor)
        labels[idx] = legname
    # call plt.legend() with the new values
    plt.legend(handles,labels)

    #cms_blurb = plt.text(1., 1., r"1 fb$^{-1}$ (13 TeV %s)",
    cms_blurb = plt.text(
        1., 1., r"CMS Preliminary %s,  %s (13 TeV %s)" % (jet_mults[jmult], lep_types[lep_type], args.year),
        fontsize=12, 
        horizontalalignment='right', 
        verticalalignment='bottom', 
        transform=ax.transAxes
    )


    pltdir = outdir if args.debug else '/'.join([outdir, jmult, lep_type])
    if not os.path.isdir(pltdir):
        os.makedirs(pltdir)
    figname = '%s/%s_%s.png' % (pltdir, args.sample, hname)

    fig.savefig(figname)
    print('%s written' % figname)
    plt.close()
