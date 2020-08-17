# matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend('agg')
from matplotlib import rcParams
rcParams['font.size'] = 20
rcParams["savefig.format"] = 'png'
rcParams["savefig.bbox"] = 'tight'
from coffea.util import load, save
from pdb import set_trace
import numpy as np
import os, re
import Utilities.plot_tools as plot_tools
from coffea.lookup_tools.dense_lookup import dense_lookup
import Utilities.prettyjson as prettyjson
from Utilities import Plotter

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('--year', type=str, help='Choose year(s) to run')
parser.add_argument('--plot_inputs', action='store_true', help='Plot ttJets distributions for calculating effective lumi')
parser.add_argument('--plot_weights', action='store_true', help='Plot gen weight distributions')

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
    print(year)
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

    if args.plot_inputs:
        pltdir = '/'.join([proj_dir, 'plots', 'Signal_Dists', 'Weights', 'ttJets'])
        if not os.path.isdir(pltdir):
            os.makedirs(pltdir)
        fig, ax = plt.subplots()
        fig.subplots_adjust(hspace=.07)
        ax = Plotter.plot_2d_norm('', ax=ax, xlimits=(300., 2000.), ylimits=(-1., 1.),
                xbins=tt_lumi_histo.dense_axes()[0].edges(), ybins=tt_lumi_histo.dense_axes()[1].edges(),
                values=tt_lumi_histo.values()[()], mask=(tt_lumi_histo.values()[()] == 0.), **{'cmap_label' : 'SM $t\\bar{t}$', 'xtitle' : 'm($t\\bar{t}$) [GeV]', 'ytitle' : 'cos($\\theta^{*}_{t}$)'}
            )
        ax = hep.cms.cmslabel(ax=ax, data=False, paper=False, year=year)
        figname = '%s/ttbar_mtt_vs_ctstar_appInds_%s' % (pltdir, year)
        fig.savefig(figname)
        print('%s written' % figname)
        plt.close(fig)

    
    ## get gen weights dists
    weights_dict = load('%s/signal_scripts/inputs/signal/signal_genweights.coffea' % proj_dir)
    for sig in weights_dict.keys():
        pltdir = '/'.join([proj_dir, 'plots', 'Signal_Dists', 'Weights', sig])
        if not os.path.isdir(pltdir):
            os.makedirs(pltdir)
    
        boson, mass, width, sampletype = sig.split('_')
        label = get_title(boson, mass, width, sampletype)
        xsec = abs([info['xsection'] for info in samples_json if info['name'] == sig ][0])
    
        if sampletype == 'Int':
            genwt_pos = weights_dict[sig]['pos']['Central']
            genwt_neg = weights_dict[sig]['neg']['Central']

                ## get effective lumi
            lumi_pos = np.sum(np.multiply(genwt_pos._values, tt_lumi_histo.values()[()]))/xsec
            lumi_neg = np.sum(np.multiply(genwt_neg._values, tt_lumi_histo.values()[()]))/xsec
            lumi_dict[year][sig] = {'pos' : lumi_pos, 'neg' : lumi_neg*-1}

            if args.plot_weights and year == '2016':
                pos_fig, pos_ax = plt.subplots()
                pos_fig.subplots_adjust(hspace=.07)
                pos_ax = Plotter.plot_2d_norm('', ax=pos_ax, xlimits=(300., 2000.), ylimits=(-1., 1.),
                        xbins=genwt_pos._axes[0], ybins=genwt_pos._axes[1],
                        values=genwt_pos._values, mask=(genwt_pos._values == 0.), **{'cmap_label' : 'w=$\Phi$/SM $t\\bar{t}$ Gen-Level', 'xtitle' : 'm($t\\bar{t}$) [GeV]', 'ytitle' : 'cos($\\theta^{*}_{t}$)'}
                    )
                hep.label.lumitext(text=label+', pos', ax=pos_ax)
                pos_figname = '%s/%s_mtt_vs_ctstar_GenWeights_pos' % (pltdir, sig)
                pos_fig.savefig(pos_figname)
                print('%s written' % pos_figname)
                plt.close(pos_fig)
    
                neg_fig, neg_ax = plt.subplots()
                neg_fig.subplots_adjust(hspace=.07)
                neg_ax = Plotter.plot_2d_norm('', ax=neg_ax, xlimits=(300., 2000.), ylimits=(-1., 1.),
                        xbins=genwt_neg._axes[0], ybins=genwt_neg._axes[1],
                        values=genwt_neg._values, mask=(genwt_neg._values == 0.), **{'cmap_label' : 'w=$\Phi$/SM $t\\bar{t}$ Gen-Level', 'xtitle' : 'm($t\\bar{t}$) [GeV]', 'ytitle' : 'cos($\\theta^{*}_{t}$)'}
                    )
                hep.label.lumitext(text=label+', neg', ax=neg_ax)
                neg_figname = '%s/%s_mtt_vs_ctstar_GenWeights_neg' % (pltdir, sig)
                neg_fig.savefig(neg_figname)
                print('%s written' % neg_figname)
                plt.close(neg_fig)
    
        else:
            genwt_pos = weights_dict[sig]['pos']['Central']
    
                ## get effective lumi
            lumi_pos = np.sum(np.multiply(genwt_pos._values, tt_lumi_histo.values()[()]))/xsec
            lumi_dict[year][sig] = {'pos' : lumi_pos}

            if args.plot_weights and year == '2016':
                fig, ax = plt.subplots()
                fig.subplots_adjust(hspace=.07)
                ax = Plotter.plot_2d_norm('', ax=ax, xlimits=(300., 2000.), ylimits=(-1., 1.),
                        xbins=genwt_pos._axes[0], ybins=genwt_pos._axes[1],
                        values=genwt_pos._values, mask=(genwt_pos._values == 0.), **{'cmap_label' : 'w=$\Phi$/SM $t\\bar{t}$ Gen-Level', 'xtitle' : 'm($t\\bar{t}$) [GeV]', 'ytitle' : 'cos($\\theta^{*}_{t}$)'}
                    )
                hep.label.lumitext(text=label+', pos', ax=ax)
                figname = '%s/%s_mtt_vs_ctstar_GenWeights_pos' % (pltdir, sig)
                fig.savefig(figname)
                print('%s written' % figname)
                plt.close(fig)


    ## save weights
if len(years_to_run) == 3:
    sig_effLumi_fname = '%s/signal_scripts/results/signal_effLumi_inds12.coffea' % proj_dir
    save(lumi_dict, sig_effLumi_fname)
    print('%s written' % sig_effLumi_fname)
