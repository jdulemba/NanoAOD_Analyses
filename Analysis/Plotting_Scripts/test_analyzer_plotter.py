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
import Utilities.plot_tools as plt_tools


parser = ArgumentParser()
parser.add_argument('--year', default='2016', help='What year is the ntuple from.')
parser.add_argument('--testing', action='store_true', help='Determines where input file is.')

args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'test_analyzer'

if not args.testing:
    raise IOError("only testing supported for plotting at the moment")

fname = '%s/%s.test.coffea' % (proj_dir, analyzer) if args.testing else '%s/%s.coffea' % ('/'.join([proj_dir, 'results', jobid]), 'test')
#fname = '%s/%s.test.%s.coffea' % (proj_dir, args.sample, analyzer) if args.testing else '%s/%s.coffea' % ('/'.join([proj_dir, 'results', jobid]), args.sample)
outdir = '/'.join([proj_dir, 'plots', jobid, analyzer, 'Test']) if args.testing else '/'.join([proj_dir, 'plots', jobid, analyzer])

if not os.path.isdir(outdir):
    os.makedirs(outdir)

hists = load(fname)
data_lumi = hists['data_lumi']
#set_trace()

variables = {
    'pt' : '$p_{T}$',
    'eta' : '$\\eta$',
    'phi' : '$\\phi$',
    'E' : 'E',
    'njets' : '$n_{jets}$',
    'Prob' : '$\\lambda_{C}$',
    'MassDiscr' : '$\\lambda_{M}$',
    'NuDiscr' : '$\\lambda_{NS}$',
}

best_perm_vars = {
    'BLep' : '$b_{l}$',
    'BHad' : '$b_{h}$',
    'WJa'  : '$leading W_{jet}$',
    'WJb'  : '$subleading W_{jet}$',
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
    'Electron' : '$e$',
    'BLep' : '$b_{l}$',
    'BHad' : '$b_{h}$',
    'WJa'  : '$leading\ W_{jet}$',
    'WJb'  : '$subleading\ W_{jet}$',
}

    ## get plotting colors/settings
hstyles = styles.styles
stack_fill_opts = {'alpha': 0.8, 'edgecolor':(0,0,0,.5)}
stack_error_opts = {'edgecolor':(0,0,0,.5)}


for hname in hists.keys():
    if hname == 'cutflow': continue
    if 'BestPerm' not in hname: continue
    histo = hists[hname]
    #set_trace()

    if 'BestPerm' in hname:
        if 'njets' in hname or 'Prob' in hname or 'Discr' in hname:
            jmult, lep_type, perm_type, kvar = hname.split('_')
        else:
            jmult, lep_type, perm_type, obj_type, kvar = hname.split('_')
    else:
        jmult, lep_type, obj_type, kvar = hname.split('_')

    fig, ax = plt.subplots(1, 1, figsize=(7,7))
    fig.subplots_adjust(hspace=.07)
    if histo.dense_dim() == 1:
        #set_trace()
        plot.plot1d(histo, ax=ax, stack=True,
            fill_opts=stack_fill_opts,
            error_opts=stack_error_opts
        )
        if kvar in variables.keys():
            if kvar == 'njets' or kvar == 'Prob' or 'Discr' in kvar:
                xtitle = variables[kvar]
            elif kvar == 'pt' or kvar == 'E':
                xtitle = '%s(%s) [GeV]' % (variables[kvar], objtypes[obj_type])
            else:
                xtitle = '%s(%s)' % (variables[kvar], objtypes[obj_type])
            plt.xlabel(xtitle)
        else:
            plt.xlabel('$%s$' % histo.axes()[-1].label)
        plt.ylabel('Events')

    ax.autoscale(axis='x', tight=True)

        ## set legend and corresponding colors
    handles, labels = plt.gca().get_legend_handles_labels()
    for idx, sample in enumerate(labels):
        #print(sample)
        #set_trace()
        facecolor, legname = plt_tools.get_styles(sample, hstyles)
        handles[idx].set_facecolor(facecolor)
        labels[idx] = legname
    # call plt.legend() with the new values
    plt.legend(handles,labels)

    cms_blurb = plt.text(
        0., 1., r"CMS Preliminary",
        fontsize=12, 
        horizontalalignment='left', 
        verticalalignment='bottom', 
        transform=ax.transAxes,
        style='italic'
    )
    lumi_blurb = plt.text(
        1., 1., r"(13 TeV %.2f fb$^{-1}$, %s + %s)" % (data_lumi, jet_mults[jmult], lep_types[lep_type]),
        fontsize=12, 
        horizontalalignment='right', 
        verticalalignment='bottom', 
        transform=ax.transAxes
    )

    pltdir = outdir if args.testing else '/'.join([outdir, jmult, lep_type])
    if not os.path.isdir(pltdir):
        os.makedirs(pltdir)
    figname = '%s/%s.png' % (pltdir, hname)
    #figname = '%s/%s_%s.png' % (pltdir, args.sample, hname)

    fig.savefig(figname)
    print('%s written' % figname)
    plt.close()
