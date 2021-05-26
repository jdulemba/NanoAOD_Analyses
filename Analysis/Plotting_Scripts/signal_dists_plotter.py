# matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend('agg')
from matplotlib import rcParams
rcParams['font.size'] = 18
rcParams["savefig.format"] = 'png'
rcParams["savefig.bbox"] = 'tight'
from coffea.util import load#, save
from pdb import set_trace
import os
#import Utilities.prettyjson as prettyjson
import numpy as np
import Plotter as Plotter

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('plot', choices=['sig_inputs', 'weights', 'all'], help='Plot raw signal input distributions, signal weights for each year and signal topology, or all of the above.')
args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'Signal_Dists'

outdir = '/'.join([proj_dir, 'plots', analyzer])
if not os.path.isdir(outdir):
    os.makedirs(outdir)

widthTOname = lambda width : str(width).replace('.', 'p')
nameTOwidth = lambda width : str(width).replace('p', '.')

def get_title(boson, mass, width, sampletype):
    mass_title = 'm($%s$)=%s GeV' % (boson[0], mass.split('M')[-1])
    width_title = '$\\Gamma$/m($%s$)=%s%%' % (boson[0], nameTOwidth(width.split('W')[-1]))
    samp_title = sampletype[0]

    return '%s, %s, %s' % (mass_title, width_title, samp_title)

mtt_title, mtt_lims = 'm($t\\bar{t}$) [GeV]', (300., 2000.)
ct_title, ct_lims = 'cos($\\theta^{*}_{t}$)', (-1., 1.)

if (args.plot == 'sig_inputs') or (args.plot == 'all'):
    signals = load('%s/signal_scripts/inputs/signal/signal_dists.coffea' % proj_dir)

    plt_kwargs = {
        'cmap_label' : 'Events',
        'xtitle' : mtt_title,
        'ytitle' : ct_title,
    } 

    for sig in signals.keys():
        boson, mass, width, sampletype = sig.split('_')
        label = get_title(boson, mass, width, sampletype) 
    
        pltdir = '/'.join([outdir, boson])
        if not os.path.isdir(pltdir):
            os.makedirs(pltdir)
    
        central_lookup = signals[sig]['Central']
        values = central_lookup._values
        mtt_bins = central_lookup._axes[0]
        ctstar_bins = central_lookup._axes[1]
    
        ## make 2D plots
        fig_2d, ax_2d = plt.subplots()
        fig_2d.subplots_adjust(hspace=.07)
        ax_2d = Plotter.plot_2d_norm('', ax=ax_2d, xlimits=mtt_lims, ylimits=ct_lims,
                values=values, xbins=mtt_bins, ybins=ctstar_bins, mask=(values == 0.), **plt_kwargs
            )
        ax_2d = hep.label.lumitext(text=label, ax=ax_2d)
        fig_2d_name = '%s/%s_mtt_vs_ctstar' % (pltdir, sig)
        fig_2d.savefig(fig_2d_name)
        print('%s written' % fig_2d_name)
        plt.close(fig_2d)
    
        ## make 1D plots
            # ctstar
        fig_ct, ax_ct = plt.subplots()
        fig_ct.subplots_adjust(hspace=.07)
        costh_values = np.sum(values, axis=0)
    
        ax_ct = Plotter.plot_1D(costh_values, ctstar_bins, xlimits=ct_lims, xlabel=ct_title, ax=ax_ct, histtype='step')
        ax_ct = hep.label.lumitext(text=label, ax=ax_ct)
        fig_ct_name = '%s/%s_ctstar' % (pltdir, sig)
        fig_ct.savefig(fig_ct_name)
        print('%s written' % fig_ct_name)
        plt.close(fig_ct)
    
            #mtt 
        fig_mtt, ax_mtt = plt.subplots()
        fig_mtt.subplots_adjust(hspace=.07)
        mtt_values = np.sum(values, axis=1)
    
        ax_mtt = Plotter.plot_1D(mtt_values, mtt_bins, xlimits=mtt_lims, xlabel=mtt_title, ax=ax_mtt, histtype='step')
        ax_mtt = hep.label.lumitext(text=label, ax=ax_mtt)
        fig_mtt_name = '%s/%s_mtt' % (pltdir, sig)
        fig_mtt.savefig(fig_mtt_name)
        print('%s written' % fig_mtt_name)
        plt.close(fig_mtt)


    ## plot weights for each year and distribution
if (args.plot == 'weights') or (args.plot == 'all'):
    weights = load('%s/Corrections/%s/signal_weights_%s_inds0.coffea' % (proj_dir, jobid, jobid))

    wt_kwargs = {
        'cmap_label' : 'w=$\Phi$/$t\\bar{t}$',
        'xtitle' : mtt_title,
        'ytitle' : ct_title,
    }

    for year in ['2016', '2017', '2018']:
        for sig in weights[year].keys():
            boson, mass, width, sampletype = sig.split('_')
            label = get_title(boson, mass, width, sampletype) 
    
            pltdir = '/'.join([outdir, 'Weights', year, boson])
            if not os.path.isdir(pltdir):
                os.makedirs(pltdir)
    
            sig_weight = weights[year][sig]
            values = sig_weight._values
            mtt_bins = sig_weight._axes[0]
            ctstar_bins = sig_weight._axes[1]
    
            ## make 2D plots
            fig, ax = plt.subplots()
            fig.subplots_adjust(hspace=.07)
            ax = Plotter.plot_2d_norm('', ax=ax, xlimits=mtt_lims, ylimits=ct_lims,
                    values=values, xbins=mtt_bins, ybins=ctstar_bins, mask=(values == 0.), **wt_kwargs
                )
            ax = hep.label.lumitext(text=label, ax=ax)
            figname = '%s/%s_mtt_vs_ctstar_Weights_%s' % (pltdir, sig, year)
            fig.savefig(figname)
            print('%s written' % figname)
            plt.close(fig)
    
