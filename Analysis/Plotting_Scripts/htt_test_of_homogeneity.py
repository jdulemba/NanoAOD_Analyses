from coffea.hist import plot
# matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend('agg')
from coffea.util import load
from pdb import set_trace
import os
import styles
import Utilities.plot_tools as plt_tools
import Utilities.prettyjson as prettyjson
import re
from coffea import hist
import numpy as np
from equal_split import partition_list

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('year', choices=['2016', '2017', '2018'], help='What year is the ntuple from.')
parser.add_argument('lepton', choices=['Electron', 'Muon'], help='Choose which lepton to make plots for')

args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'htt_test_of_homogeneity'

input_dir = '/'.join([proj_dir, 'results', '%s_%s' % (args.year, jobid), analyzer])
f_ext = 'TOT.coffea'
outdir = '/'.join([proj_dir, 'plots', '%s_%s' % (args.year, jobid), analyzer])
if not os.path.isdir(outdir):
    os.makedirs(outdir)

fnames = sorted(['%s/%s' % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)])

#set_trace()
hdict = plt_tools.add_coffea_files(fnames) if len(fnames) > 1 else load(fnames[0])

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

lep_cats = {
    'Tight' : 'tight %s' % objtypes['Lep'][args.lepton],
    'Loose' : 'loose %s' % objtypes['Lep'][args.lepton],
}

variables = {
    'mtt_vs_iso' : ('m($t\\bar{t}$) [GeV]', 'pfRelIso, %s' % objtypes['Lep'][args.lepton], 2, (200., 2000.), 1, (0., 2.) if args.lepton == 'Muon' else (0., 1.5), True),
    'tlep_ctstar_vs_iso' : ('cos($\\theta^{*}_{t_{l}}$)', 'pfRelIso, %s' % objtypes['Lep'][args.lepton], 2, (-1., 1.), 1, (0., 2.) if args.lepton == 'Muon' else (0., 1.5), True),
}


    ## get plotting colors/settings
hstyles = styles.styles
stack_fill_opts = {'alpha': 0.8, 'edgecolor':(0,0,0,.5)}
stack_error_opts = {'edgecolor':(0,0,0,.5)}

    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open('%s/inputs/lumis_data.json' % proj_dir).read())[args.year]
lumi_correction = load('%s/Corrections/%s/MC_LumiWeights_IgnoreSigEvts.coffea' % (proj_dir, jobid))
for hname in hdict.keys():
    if hname == 'cutflow': continue
    hdict[hname].scale(lumi_correction[args.year]['%ss' % args.lepton], axis='dataset')


## make data and mc categories for data/MC plotting
mc_samples = re.compile('(?!data*)')
data_samples = re.compile('(data*)')

## make groups based on process
process = hist.Cat("process", "Process", sorting='placement')
process_cat = "dataset"
process_groups = plt_tools.make_dataset_groups(args.lepton, args.year)
#set_trace()
for hname in hdict.keys():
    if hname == 'cutflow': continue
    hdict[hname] = hdict[hname].group(process_cat, process, process_groups)
    


    ## make plots
for hname in variables.keys():
    if hname not in hdict.keys():
        raise ValueError("%s not found in file" % hname)
    #set_trace()
    histo = hdict[hname][:, :, args.lepton, :, :].integrate('leptype') # process, jmult, lepton, lepcat, btagregion

    if histo.dense_dim() == 2:
        xtitle, ytitle, x_rebinning, x_lims, y_rebinning, y_lims, withData = variables[hname]
        xvar, yvar = hname.split('_vs_')

        #set_trace()
        ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
        for jmult in list(set([key[1] for key in histo.values().keys()])):
            for btagregion in list(set([key[3] for key in histo.values().keys()])):
                for lepcat in list(set([key[2] for key in histo.values().keys()])):
                    pltdir = '/'.join([outdir, args.lepton, jmult, lepcat, btagregion])
                    if not os.path.isdir(pltdir):
                        os.makedirs(pltdir)

                    hslice = histo[:, jmult, lepcat, btagregion].integrate('jmult').integrate('lepcat').integrate('btag')

                    if x_rebinning != 1:
                        xaxis_name = hslice.dense_axes()[0].name
                        hslice = hslice.rebin(xaxis_name, x_rebinning)
                    if y_rebinning != 1:
                        yaxis_name = hslice.dense_axes()[1].name
                        hslice = hslice.rebin(yaxis_name, y_rebinning)


                        # make 1D projection along dense axes
                    for dax in range(2):
                        fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                        fig.subplots_adjust(hspace=.07)

                        #set_trace()
                        hproj = hslice.integrate(hslice.dense_axes()[dax].name)
                            ## plot comparison of perm categories
                        plot.plot1d(hproj[mc_samples],
                            overlay=hproj.axes()[0].name,
                            ax=ax,
                            clear=False,
                            stack=True,
                            line_opts=None,
                            fill_opts=stack_fill_opts,
                            error_opts=stack_error_opts,
                            order=['QCD', 'singlet', 'EWK', 'ttJets'],
                        )
                        if withData:
                            plot.plot1d(hproj[data_samples],
                                overlay=hproj.axes()[0].name,
                                ax=ax,
                                clear=False,
                                error_opts=hstyles['data_err_opts']
                            )
                        ax.autoscale(axis='x', tight=True)
                        ax.set_ylim(0, None)
                        ax.set_ylabel(None)
                        ax.set_xlabel(None)
                        if hproj.dense_axes()[0].name == 'iso':
                            mc_vals = np.sum(np.array(list(hproj[mc_samples].values().values())), axis=0)
                            data_vals = np.sum(np.array(list(hproj[data_samples].values().values())), axis=0)
                            max_binval = np.maximum(np.max(mc_vals), np.max(data_vals))
                            min_binval = np.minimum(np.min(mc_vals), np.min(data_vals))
                            ax.set_yscale('log')
                            ax.set_ylim(np.maximum(1e-1, min_binval), 10**(np.ceil(np.log10(max_binval))+1))
                            if lepcat == 'Tight':
                                ax.set_xlim((0., 0.2))
                        else:
                            ax.set_xlim(y_lims) if dax == 0 else ax.set_xlim(x_lims)

                            ## set legend and corresponding colors
                        handles, labels = ax.get_legend_handles_labels()
                        for idx, sample in enumerate(labels):
                            if sample == 'data' or sample == 'Observed': continue
                            facecolor, legname = plt_tools.get_styles(sample, hstyles)
                            handles[idx].set_facecolor(facecolor)
                            labels[idx] = legname
                        # call ax.legend() with the new values
                        ax.legend(handles,labels, loc='upper right')
                        #set_trace()

                        if withData:
                                ## plot data/MC ratio
                            plot.plotratio(hproj[data_samples].sum(hproj.axes()[0].name), hproj[mc_samples].sum(hproj.axes()[0].name),
                                ax=rax,
                                error_opts=hstyles['data_err_opts'],
                                denom_fill_opts={},
                                guide_opts={},
                                unc='num'
                            )
                            rax.set_ylabel('data/MC')
                            rax.set_ylim(0.5, 1.5)
                            if lepcat == 'Tight' and hproj.dense_axes()[0].name == 'iso':
                                rax.set_xlim((0., 0.2))
                            else:
                                rax.set_xlim(y_lims) if dax == 0 else rax.set_xlim(x_lims)

                        #set_trace()
                            # add lep category
                        ax.text(
                            0.02, 0.91, "%s, %s" % (lep_cats[lepcat], jet_mults[jmult]),
                            fontsize=18,
                            horizontalalignment='left',
                            verticalalignment='bottom',
                            transform=ax.transAxes
                        )

                            ## set axes labels and titles
                        plt.xlabel(ytitle) if dax == 0 else plt.xlabel(xtitle)
                        ax = hep.cms.cmslabel(ax=ax, data=withData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))

                        #set_trace()
                        #figname = 'test.png'
                        figname = '%s/%s.png' % (pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, yvar if dax == 0 else xvar]))
                        fig.savefig(figname, bbox_inches='tight')
                        print('%s written' % figname)
                        plt.close()
                        #set_trace()


