# matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend('agg')
from matplotlib import rcParams
rcParams['font.size'] = 18
rcParams["savefig.format"] = 'png'
rcParams["savefig.bbox"] = 'tight'
from coffea.util import load, save
from pdb import set_trace
import numpy as np
import os, re
import Utilities.plot_tools as plot_tools
from coffea.lookup_tools.dense_lookup import dense_lookup
import Utilities.prettyjson as prettyjson
from Plotting_Scripts import Plotter

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('--year', type=str, help='Choose year(s) to run')
parser.add_argument('--plot_inputs', action='store_true', help='Plot input distributions for weights')

args = parser.parse_args()

years_to_run = [args.year] if args.year else ['2016', '2017', '2018']

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'Signal_Dists'

tt_samples_dict = {
    '2016' : {
        'ttJets_PS' : 0., # 0. placeholder for xsection weight
    },
    '2017' : {
        'ttJetsSL' : 0.,
        'ttJetsHad': 0.,
        'ttJetsDiLep': 0.,
    },
    '2018' : {
        'ttJetsSL' : 0.,
        'ttJetsHad': 0.,
        'ttJetsDiLep': 0.,
    },
}

derivation_inds = re.compile('(0)')
application_inds = re.compile(r'\b(?:%s)\b' % '|'.join(['1', '2']))

weights_dict = {
    '2016' : {},
    '2017' : {},
    '2018' : {},
}
lumi_dict = {
    '2016' : {},
    '2017' : {},
    '2018' : {},
}

widthTOname = lambda width : str(width).replace('.', 'p')
nameTOwidth = lambda width : str(width).replace('p', '.')

def get_title(boson, mass, width, sampletype):
    mass_title = 'm($%s$)=%s GeV' % (boson[0], mass.split('M')[-1])
    width_title = '$\\Gamma$/m($%s$)=%s%%' % (boson[0], nameTOwidth(width.split('W')[-1]))
    samp_title = sampletype[0]

    return '%s, %s, %s' % (mass_title, width_title, samp_title)


for year in years_to_run:
        ## get branching ratio for each tt decay
    samples_json =  prettyjson.loads(open('%s/inputs/samples_%s.json' % (proj_dir, year)).read())
    ttweight_dict = {tt : [info['xsection'] for info in samples_json if info['name'] == tt ][0] for tt in tt_samples_dict[year].keys()}
    for tt in tt_samples_dict[year].keys():
        tt_samples_dict[year][tt] = ttweight_dict[tt]/np.sum(np.array(list(ttweight_dict.values())))
    ## get ttJets distributions
    tt_fnames = ['/'.join([proj_dir, 'signal_scripts', 'inputs', 'ttJets_inputs', '%s_%s_meta_info.coffea' % (year, fname)]) for fname in tt_samples_dict[year].keys()]

        ## weight hists by branching ratio and then combine
    tt_dict = plot_tools.add_coffea_files(tt_fnames)
    tt_histo = tt_dict['mtt_topctstar']
    tt_histo.scale(tt_samples_dict[year], axis='dataset')
    tt_histo = tt_histo.integrate('dataset')
    tt_lumi_histo = tt_histo[application_inds].integrate('evtIdx')
    tt_histo = tt_histo[derivation_inds].integrate('evtIdx')
    tt_lookup = dense_lookup(*(tt_histo.values().values()), (tt_histo.dense_axes()[0].edges(), tt_histo.dense_axes()[1].edges()))
    norm_tt_vals = tt_lookup._values/np.sum(tt_histo.values(overflow='all')[()]) # normalize array
        ## plot input ttbar distribution
    if args.plot_inputs:
        pltdir = '/'.join([proj_dir, 'plots', 'Signal_Dists', 'Weights', 'inputs'])
        if not os.path.isdir(pltdir):
            os.makedirs(pltdir)

        fig, ax = plt.subplots()
        fig.subplots_adjust(hspace=.07)
        ax = Plotter.plot_2d_norm('', ax=ax, xlimits=(300., 2000.), ylimits=(-1., 1.),
                xbins=tt_histo.axis('mtt').edges(), ybins=tt_histo.axis('ctstar').edges(),
                values=norm_tt_vals, mask=(norm_tt_vals == 0.), **{'cmap_label' : '$t\\bar{t}$ Probability', 'xtitle' : 'm($t\\bar{t}$) [GeV]', 'ytitle' : 'cos($\\theta^{*}_{t}$)'}
            )
        ax = hep.cms.cmslabel(ax=ax, data=False, paper=False, year=year)
        figname = '%s/Norm_ttbar_mtt_vs_ctstar_inputs_%s' % (pltdir, year)
        fig.savefig(figname)
        print('%s written' % figname)
        plt.close(fig)
    
    ## get signal_dists
    signal_dict = load('%s/signal_scripts/inputs/signal/signal_dists.coffea' % proj_dir)
    for sig in signal_dict.keys():
        boson, mass, width, sampletype = sig.split('_')
        label = get_title(boson, mass, width, sampletype)

        xsec = [info['xsection'] for info in samples_json if info['name'] == sig ][0]
        sig_dist = signal_dict[sig]['Central']
        norm_sig_vals = sig_dist._values/(np.sum(sig_dist._values[sig_dist._values > 0]) + np.abs(np.sum(sig_dist._values[sig_dist._values < 0])) ) # normalize based on number of nonzero events
            ## plot input signal distributions
        if args.plot_inputs:
            if year != '2016': continue # only need to make input signal plots once
            pltdir = '/'.join([proj_dir, 'plots', 'Signal_Dists', 'Weights', 'inputs'])
            if not os.path.isdir(pltdir):
                os.makedirs(pltdir)

            fig, ax = plt.subplots()
            fig.subplots_adjust(hspace=.07)
            ax = Plotter.plot_2d_norm('', ax=ax, xlimits=(300., 2000.), ylimits=(-1., 1.),
                    xbins=sig_dist._axes[0], ybins=sig_dist._axes[1],
                    values=norm_sig_vals, mask=(norm_sig_vals == 0.), **{'cmap_label' : 'Probability', 'xtitle' : 'm($t\\bar{t}$) [GeV]', 'ytitle' : 'cos($\\theta^{*}_{t}$)'}
                )
            hep.label.lumitext(text=label, ax=ax)
            #ax = hep.cms.cmslabel(ax=ax, data=False, paper=False, year=year)
            figname = '%s/Norm_%s_mtt_vs_ctstar_inputs' % (pltdir, sig)
            fig.savefig(figname)
            print('%s written' % figname)
            plt.close(fig)
    
            # make phi/tt weights
        phi_over_tt_weights = np.divide(norm_sig_vals, norm_tt_vals)
        phi_over_tt_weights[norm_tt_vals == 0.] = 0. # set all bins where tt binvals == 0. to 0.

        ## save as lookup table
        weights_lookup = dense_lookup((phi_over_tt_weights), sig_dist._axes)
        weights_dict[year][sig] = weights_lookup

            ## get effective lumi
        lumi = np.sum(np.multiply(phi_over_tt_weights, tt_lumi_histo.values()[()]))/xsec
        lumi_dict[year][sig] = lumi

if len(years_to_run) == 3:
    sigweights_outfname = '%s/Corrections/%s/signal_weights_%s_inds0.coffea' % (proj_dir, jobid, jobid)
    save(weights_dict, sigweights_outfname)
    print('%s written' % sigweights_outfname)

    sig_effLumi_fname = '%s/signal_scripts/results/signal_effLumi_inds12.coffea' % proj_dir
    save(lumi_dict, sig_effLumi_fname)
    print('%s written' % sig_effLumi_fname)
