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
parser.add_argument('--sample', type=str, help='Input sample to use.')
parser.add_argument('--testing', action='store_true', help='Determines where input file is.')
parser.add_argument('--use_combined', action='store_true', help='Used file that has multiple datasets.')

args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'permProbComputer'

input_dir = proj_dir if args.testing else '/'.join([proj_dir, 'results', '%s_%s' % (args.year, jobid), analyzer])
f_ext = '%s.test.coffea' % analyzer if args.testing else 'TOT.coffea'
outdir = '/'.join([proj_dir, 'plots', '%s_%s' % (args.year, jobid), analyzer, 'Test']) if args.testing else '/'.join([proj_dir, 'plots', '%s_%s' % (args.year, jobid), analyzer])
if not os.path.isdir(outdir):
    os.makedirs(outdir)

if args.testing:
    if args.use_combined:
        fnames = ['%s/%s' % (input_dir, f_ext)]
    else:
        fnames = ['%s/%s' % (input_dir, fname) for fname in os.listdir(input_dir) if fname == '%s_%s' % (args.sample, f_ext)] if args.sample else ['%s/%s' % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)]
        fnames.remove('%s/%s' % (input_dir, f_ext))
else:
    fnames = ['%s/%s%s' % (input_dir, args.sample, f_ext)] if args.sample else ['%s/%s' % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)]
fnames = sorted(fnames)

#set_trace()
hdict = plt_tools.add_coffea_files(fnames) if len(fnames) > 1 else load(fnames[0])

jet_mults = {
    #'3Jets' : '3 jets',
    '4PJets' : '4+ jets',
    '4Jets' : '4 jets',
    '5Jets' : '5 jets',
    '5PJets' : '5+ jets'
}

lep_cats = {
    'LoT' : 'loose or tight $e/\mu$',
    'Tight' : 'tight $e/\mu$',
}

variables_1d = {
    #'Common' : {
        'nusolver_chi2' : ('$\\chi_{\\nu}^{2}$', 5, (0., 1000.)),
        'nusolver_dist' : ('$D_{\\nu, min}$', 1, (0., 150.)),
    #}
}

variables_2d = {
    #'4PJets' : {
        'mWHad_vs_mTHad' : ('m($W_{had}$) [GeV]', 'm($t_{had}$) [GeV]', 1, (0., 500.), 1, (0., 500.))
    #},
}


    ## get plotting colors/settings
hstyles = styles.styles
stack_fill_opts = {'alpha': 0.8, 'edgecolor':(0,0,0,.5)}
stack_error_opts = {'edgecolor':(0,0,0,.5)}

    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open('%s/inputs/lumis_data.json' % proj_dir).read())[args.year]
lumi_to_use = (data_lumi_year['Muons']+data_lumi_year['Electrons'])/2000.
lumi_correction = load('%s/Corrections/%s/MC_LumiWeights.coffea' % (proj_dir, jobid))

## make data and mc categories for data/MC plotting
mc_samples = re.compile('(?!data*)')
data_samples = re.compile('(data*)')

## make groups based on process
process = hist.Cat("process", "Process", sorting='placement')
process_cat = "dataset"
process_groups = plt_tools.make_dataset_groups('Muon', args.year) # works with MC
#set_trace()
for hname in hdict.keys():
    if hname == 'cutflow': continue
    hdict[hname] = hdict[hname].group(process_cat, process, process_groups)
    


    ## make plots
for hname in hdict.keys():
    if (hname not in variables_1d.keys()) and (hname not in variables_2d.keys()): continue

    histo = hdict[hname]
        ## rescale hist by lumi for muons and electrons separately and then combine
    h_mu = histo[:, :, 'Muon', :].integrate('leptype')
    h_mu.scale(lumi_correction[args.year]['Muons'], axis='process')
    h_el = histo[:, :, 'Electron', :].integrate('leptype')
    h_el.scale(lumi_correction[args.year]['Electrons'], axis='process')
    h_tot = h_mu+h_el
    h_tot = h_tot.integrate('process')

    #set_trace()

    if h_tot.dense_dim() == 1:
        xtitle, rebinning, x_lims = variables_1d[hname]

        ## hists should have 3 category axes (dataset, jet multiplicity, lepton type) followed by variable
        for jmult in h_tot.axis('jmult')._sorted:
            for lepcat in h_tot.axis('lepcat')._sorted:
                pltdir = outdir if args.testing else '/'.join([outdir, jmult, lepcat])
                if not os.path.isdir(pltdir):
                    os.makedirs(pltdir)

                fig, ax = plt.subplots(1, 1, figsize=(7,7))
                fig.subplots_adjust(hspace=.07)
                hslice = h_tot[jmult, lepcat].integrate('jmult').integrate('lepcat')

                if rebinning != 1:
                    xaxis_name = hslice.dense_axes()[0].name
                    hslice = hslice.rebin(xaxis_name, rebinning)

                    ## plot MC and data
                plot.plot1d(hslice,
                    #overlay=hslice.axes()[0].name,
                    ax=ax,
                    clear=False,
                    error_opts=hstyles['data_err_opts']
                )
                ax.autoscale(axis='x', tight=True)
                ax.set_ylim(0, None)
                ax.set_xlabel(None)
                ax.set_xlim(x_lims)

                #    ## set legend and corresponding colors
                #handles, labels = ax.get_legend_handles_labels()
                #for idx, sample in enumerate(labels):
                #    if sample == 'data' or sample == 'Observed': continue
                #    facecolor, legname = plt_tools.get_styles(sample, hstyles)
                #    handles[idx].set_facecolor(facecolor)
                #    labels[idx] = legname
                ## call ax.legend() with the new values
                #ax.legend(handles,labels, loc='upper right')
                ax.get_legend().remove()
                #set_trace()

                    # add lep category
                ax.text(
                    0.02, 0.95, lep_cats[lepcat],
                    fontsize=10, 
                    horizontalalignment='left', 
                    verticalalignment='bottom', 
                    transform=ax.transAxes
                )
                    ## set axes labels and titles
                plt.xlabel(xtitle)
                cms_blurb = plt.text(
                    0., 1., r"CMS Preliminary",
                    fontsize=12, 
                    horizontalalignment='left', 
                    verticalalignment='bottom', 
                    transform=ax.transAxes,
                    style='italic'
                )
                lumi_blurb = plt.text(
                        1., 1., r"(13 TeV %.2f fb$^{-1}$, %s)" % (lumi_to_use, jet_mults[jmult]),
                        #1., 1., r"(13 TeV %.2f fb$^{-1}$, %s, %s)" % (lumi_to_use, lep_cats[lepcat], jet_mults[jmult]),
                    fontsize=12, 
                    horizontalalignment='right', 
                    verticalalignment='bottom', 
                    transform=ax.transAxes
                )

                figname = '%s/%s.png' % (pltdir, '_'.join([jmult, lepcat, hname]))
                fig.savefig(figname)
                print('%s written' % figname)
                plt.close()
                #set_trace()

    elif h_tot.dense_dim() == 2:
        xtitle, ytitle, x_rebinning, x_lims, y_rebinning, y_lims = variables_2d[hname]

        #set_trace()
        ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable(s)
        for jmult in h_tot.axis('jmult')._sorted:
            for lepcat in h_tot.axis('lepcat')._sorted:
                pltdir = outdir if args.testing else '/'.join([outdir, jmult, lepcat])
                if not os.path.isdir(pltdir):
                    os.makedirs(pltdir)

                hslice = h_tot[jmult, lepcat].integrate('jmult').integrate('lepcat')

                if x_rebinning != 1:
                    xaxis_name = hslice.dense_axes()[0].name
                    hslice = hslice.rebin(xaxis_name, x_rebinning)
                if y_rebinning != 1:
                    yaxis_name = hslice.dense_axes()[1].name
                    hslice = hslice.rebin(yaxis_name, y_rebinning)

                #set_trace()
                    # make 1D projection along dense axes
                for dax in range(2):
                    fig, ax = plt.subplots(1, 1, figsize=(7,7))
                    fig.subplots_adjust(hspace=.07)

                    #set_trace()
                    hproj = hslice[:, :].integrate(hslice.dense_axes()[dax].name)
                        ## plot MC and data
                    plot.plot1d(hproj,
                        ax=ax,
                        clear=False,
                        error_opts=hstyles['data_err_opts']
                    )
                    ax.autoscale(axis='x', tight=True)
                    ax.set_ylim(0, None)
                    ax.set_xlabel(None)
                    ax.set_xlim(y_lims) if dax == 0 else ax.set_xlim(x_lims)

                    #    ## set legend and corresponding colors
                    #handles, labels = ax.get_legend_handles_labels()
                    #for idx, sample in enumerate(labels):
                    #    if sample == 'data' or sample == 'Observed': continue
                    #    facecolor, legname = plt_tools.get_styles(sample, hstyles)
                    #    handles[idx].set_facecolor(facecolor)
                    #    labels[idx] = legname
                    ## call ax.legend() with the new values
                    #ax.legend(handles,labels, loc='upper right')
                    ax.get_legend().remove()

                    #set_trace()
                        # add lep category
                    ax.text(
                        0.02, 0.95, lep_cats[lepcat],
                        fontsize=10, 
                        horizontalalignment='left', 
                        verticalalignment='bottom', 
                        transform=ax.transAxes
                    )

                        ## set axes labels and titles
                    plt.xlabel(ytitle) if dax == 0 else plt.xlabel(xtitle)
                    cms_blurb = plt.text(
                        0., 1., r"CMS Preliminary",
                        fontsize=12, 
                        horizontalalignment='left', 
                        verticalalignment='bottom', 
                        transform=ax.transAxes,
                        style='italic'
                    )
                    lumi_blurb = plt.text(
                        1., 1., r"(13 TeV %.2f fb$^{-1}$, %s)" % (lumi_to_use, jet_mults[jmult]),
                        fontsize=12, 
                        horizontalalignment='right', 
                        verticalalignment='bottom', 
                        transform=ax.transAxes
                    )

                    #figname = 'test.png'
                    figname = '%s/%s.png' % (pltdir, '_'.join([jmult, lepcat, hname.split('_vs_')[1] if dax == 0 else hname.split('_vs_')[0]]))
                    fig.savefig(figname)
                    print('%s written' % figname)
                    plt.close()
                    #set_trace()

                    # make 2D plots for each topology
                fig, ax = plt.subplots(1, 1, figsize=(7,7))
                fig.subplots_adjust(hspace=.07)

                #set_trace()
                    ## plot 
                plot.plot2d(hslice,
                    xaxis=hslice.axes()[0].name,
                    ax=ax,
                    patch_opts={'cmap' : 'OrRd'},
                    clear=True,
                )
                ax.autoscale(axis='x', tight=True)
                ax.set_ylim(y_lims)
                ax.set_xlim(x_lims)


                    # add lep category
                ax.text(
                    0.02, 0.95, lep_cats[lepcat],
                    fontsize=10, 
                    horizontalalignment='left', 
                    verticalalignment='bottom', 
                    transform=ax.transAxes
                )
                    ## set axes labels and titles
                plt.xlabel(xtitle)
                plt.ylabel(ytitle)
                cms_blurb = plt.text(
                    0., 1., r"CMS Preliminary",
                    fontsize=12, 
                    horizontalalignment='left', 
                    verticalalignment='bottom', 
                    transform=ax.transAxes,
                    style='italic'
                )
                lumi_blurb = plt.text(
                    1., 1., r"(13 TeV %.2f fb$^{-1}$, %s)" % (lumi_to_use, jet_mults[jmult]),
                    #1., 1., r"(13 TeV %.2f fb$^{-1}$, %s, %s)" % (lumi_to_use, lep_cats[lepcat], jet_mults[jmult]),
                    fontsize=12, 
                    horizontalalignment='right', 
                    verticalalignment='bottom', 
                    transform=ax.transAxes
                )

                figname = '%s/%s.png' % (pltdir, '_'.join([jmult, lepcat, hname]))
                fig.savefig(figname)
                print('%s written' % figname)
                plt.close()
                #set_trace()

