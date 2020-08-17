# matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend('agg')
from matplotlib import rcParams
rcParams['font.size'] = 20
rcParams["savefig.format"] = 'png'
rcParams["savefig.bbox"] = 'tight'
from coffea.util import load
from pdb import set_trace
import os
import Utilities.plot_tools as plt_tools
import Utilities.prettyjson as prettyjson
#import re
from coffea import hist
import numpy as np
#import fnmatch
import Utilities.Plotter as Plotter
import Utilities.systematics as systematics
from coffea.hist import plot

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('year', choices=['2016', '2017', '2018'], help='What year is the ntuple from.')
parser.add_argument('lepton', choices=['Electron', 'Muon'], help='Choose which lepton to make plots for')
parser.add_argument('--ind', action='store_true', help='Plot individual sys variations')
parser.add_argument('--comp', action='store_true', help='Plot sys variations with nominal mc')
parser.add_argument('--ratio', action='store_true', help='Make sys variation/nominal mc ratio plots')

args = parser.parse_args()

if (not args.comp) and (not args.ratio) and (not args.ind):
    raise ValueError("No plots chosen to make")

sys_to_name = systematics.sys_to_name[args.year]

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'htt_btag_iso_cut'

input_dir = '/'.join([proj_dir, 'plots', '%s_%s' % (args.year, jobid), analyzer])

lep_dname = '%s/%s/%s_Post_QCD_Est_dict.coffea' % (input_dir, args.lepton, args.lepton)
if os.path.isfile(lep_dname):
    hdict = load(lep_dname)
else: raise ValueError("%s doesn't exist" % lep_dname)

outdir = '/'.join([proj_dir, 'plots', '%s_%s' % (args.year, jobid), 'plot_htt_sys', args.lepton])
if not os.path.isdir(outdir):
    os.makedirs(outdir)


jet_mults = {
    '3Jets' : '3 jets',
    '4PJets' : '4+ jets'
}

leptypes = {
    'Muon' : '$\\mu$',
    'Electron' : '$e$',
}

objtypes = {
    'Jets' : 'jets',
    'Lep' :  {
        'Muon' : '$\\mu$',
        'Electron' : '$e$',
    }
}

variables = {
        ## (xaxis title, xaxis limits, mask data)
    'mtt_vs_tlep_ctstar_abs' : ('m($t\\bar{t}$) $\otimes$ |cos($\\theta^{*}_{t_{l}}$)|', None, True),
    'mtt' : ('m($t\\bar{t}$) [GeV]', (200., 2000.), True),
    'tlep_ctstar_abs' : ('|cos($\\theta^{*}_{t_{l}}$)|', (0., 1.), True),
    'Jets_njets' : ('$n_{jets}$', (0, 15), False),
    'mthad' : ('m($t_{h}$) [GeV]', (0., 300.), False),
    'pt_thad' : ('$p_{T}$($t_{h}$) [GeV]', (0., 500.), False),
    'pt_tlep' : ('$p_{T}$($t_{l}$) [GeV]', (0., 500.), False),
    'pt_tt' : ('$p_{T}$($t\\bar{t}$) [GeV]', (0., 500.), False),
    'eta_thad' : ('$\\eta$($t_{h}$)', (-4., 4.), False),
    'eta_tlep' : ('$\\eta$($t_{l}$)', (-4., 4.), False),
    'eta_tt' : ('$\\eta$($t\\bar{t}$)', (-4., 4.), False),
    'tlep_ctstar' : ('cos($\\theta^{*}_{t_{l}}$)', (-1., 1.), False),
    'full_disc' : ('$\\lambda_{C}$', (5, 25.), False),
    'mass_disc' : ('$\\lambda_{M}$', (0, 20.), False),
    'ns_disc' : ('$\\lambda_{NS}$', (3., 10.), False),
    'MT' : ('$M_{T}$ [GeV]', (0., 300.), False),
    'MET_pt' : ('$p_{T}$(MET) [GeV]', (0., 300.), False),
    'MET_phi' : ('$\phi$(MET)', (-3.2, 3.2), False),
    'Jets_pt' : ('$p_{T}$(jets) [GeV]', (0., 300.), False),
    'Jets_eta' : ('$\\eta$(jets)', (-2.6, 2.6), False),
    'Jets_LeadJet_pt' : ('$p_{T}$(leading jet) [GeV]', (0., 300.), False),
    'Jets_LeadJet_eta' : ('$\\eta$(leading jet)', (-2.6, 2.6), False),
    'Lep_pt' : ('$p_{T}$(%s) [GeV]' % objtypes['Lep'][args.lepton], (0., 300.), False),
    'Lep_eta' : ('$\\eta$(%s)' % objtypes['Lep'][args.lepton], (-2.6, 2.6), False),
    'Lep_iso' : ('pfRelIso, %s' % objtypes['Lep'][args.lepton], (0., 0.15) if args.lepton == 'Muon' else (0., 0.1), False),
}


    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open('%s/inputs/lumis_data.json' % proj_dir).read())[args.year]

    # group all ttJets into one
process = hist.Cat("process", "Process", sorting='placement')
process_cat = 'process'
process_groups = {'ttJets' : ['ttJets_right', 'ttJets_matchable', 'ttJets_unmatchable', 'ttJets_other'], 'EWK' : ['EWK'], 'singlet' : ['singlet'], 'QCD' : ['QCD'], 'data' : ['data']}

    ## compare up/down variations to nominal for each systematic
#for hdict in hdicts:
#set_trace()
#lep = sorted(set([key[0] for key in hdict.keys()]))[0]
jmults = sorted(set([key[1] for key in hdict.keys()]))
systs = sorted(set([key[2] for key in hdict.keys()]))
#hnames = sorted(set([key[3] for key in hdict.keys()]))

systypes = sorted(set(['_'.join(sys.split('_')[:-1]) for sys in systs if (not sys == 'nosys') and (not ('UP' and 'DW') in sys)]))
if ('RENORM_UP_FACTOR_DW' and 'RENORM_DW_FACTOR_UP') in systs:
    systypes += ['RENFACTOR_DIFF']
for jmult in jmults:
    for hname in variables.keys():#hnames:
        xtitle, x_lims, maskData = variables[hname]
        
        nom_histo = hdict[(args.lepton, jmult, 'nosys', hname)]
        nom_histo = nom_histo.group(process_cat, process, process_groups)
        if args.ratio:
            nom_mc = nom_histo[Plotter.mc_samples].integrate('process')

        if hname == 'mtt_vs_tlep_ctstar_abs':
            x_lims = (0, nom_histo.dense_axes()[0].centers().size)

        vlines = [x_lims[1]*ybin/5 for ybin in range(1, 5)] if hname == 'mtt_vs_tlep_ctstar_abs' else None

        if args.ind:
            #set_trace()
            for sys in systs:
                print(jmult, hname, sys)
                pltdir = '/'.join([outdir, jmult, 'Individual', sys])
                if not os.path.isdir(pltdir):
                    os.makedirs(pltdir)

                opts = {
                    'legend_title' : sys,
                    'maskData' : maskData,
                }

                sys_histo = hdict[(args.lepton, jmult, sys, hname)]
                sys_histo = sys_histo.group(process_cat, process, process_groups)

                fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                fig.subplots_adjust(hspace=.07)

                if hname == 'Jets_njets':
                    yields_txt, yields_json = Plotter.get_samples_yield_and_frac(sys_histo, data_lumi_year['%ss' % args.lepton]/1000., sys=sys)
                    frac_name = '%s_%s_yields_and_fracs_QCD_Est' % (sys, '_'.join([jmult, args.lepton]))
                    plt_tools.print_table(yields_txt, filename='%s/%s.txt' % (pltdir, frac_name), print_output=True)
                    print('%s/%s.txt written' % (pltdir, frac_name))
                    with open('%s/%s.json' % (pltdir, frac_name), 'w') as out:
                        out.write(prettyjson.dumps(yields_json))
                
                ax, rax = Plotter.plot_stack1d(ax, rax, sys_histo, xlabel=xtitle, xlimits=x_lims, **opts)
                
                    # add lepton/jet multiplicity label
                ax.text(
                    0.02, 0.92, "%s, %s" % (leptypes[args.lepton], jet_mults[jmult]),
                    fontsize=rcParams['font.size']*0.9, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
                )
                        ## draw vertical lines for distinguishing different ctstar bins
                if vlines is not None:
                    for vline in vlines:
                        ax.axvline(vline, color='k', linestyle='--')
                        if rax is not None: rax.axvline(vline, color='k', linestyle='--')
                hep.cms.cmslabel(ax=ax, data=not maskData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
                
                #set_trace()
                figname = '%s/%s' % (pltdir, '_'.join([jmult, args.lepton, sys, hname]))
                fig.savefig(figname)
                print('%s written' % figname)
                plt.close()


        if (args.comp) or (args.ratio):
            for sys in systypes:
                pltdir = '/'.join([outdir, jmult, 'Variations_UP_DW', sys])
                if not os.path.isdir(pltdir):
                    os.makedirs(pltdir)
        
                if sys == 'RENFACTOR_DIFF':
                    up_histo = hdict[(args.lepton, jmult, 'RENORM_UP_FACTOR_DW', hname)]
                    dw_histo = hdict[(args.lepton, jmult, 'RENORM_DW_FACTOR_UP', hname)]
        
                else:
                    up_histo = hdict[(args.lepton, jmult, '%s_UP' % sys, hname)]
                    dw_histo = hdict[(args.lepton, jmult, '%s_DW' % sys, hname)]
                up_histo = up_histo.group(process_cat, process, process_groups)
                dw_histo = dw_histo.group(process_cat, process, process_groups)
        
                opts = {
                    'legend_title' : sys,
                    #'logy' : True,
                    'maskData' : maskData,
                }

                if args.comp:
                        ## make data/MC plots for nominal and variations        
                    fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                    fig.subplots_adjust(hspace=.07)

                        ## plot MC for variations
                    ax = hep.plot.histplot(up_histo[Plotter.mc_samples].integrate('process').values()[()], nom_histo.dense_axes()[0].edges(), ax=ax, histtype='step', **Plotter.hstyles['Up'])# up MC
                    ax = hep.plot.histplot(dw_histo[Plotter.mc_samples].integrate('process').values()[()], nom_histo.dense_axes()[0].edges(), ax=ax, histtype='step', **Plotter.hstyles['Down'])# down MC

                        # plot nominal data and MC
                    ax, rax = Plotter.plot_stack1d(ax, rax, nom_histo, xlabel=xtitle, xlimits=x_lims, **opts)

                    if not maskData:
                            ## add data/MC for variations
                        rax_ylabel = rax.get_ylabel()
                        plot.plotratio(
                            nom_histo[Plotter.data_samples].integrate('process'), up_histo[Plotter.mc_samples].integrate('process'), error_opts=Plotter.hstyles['Up'],
                            unc='num', clear=False, ax=rax, guide_opts={},
                        )
                        plot.plotratio(
                            nom_histo[Plotter.data_samples].integrate('process'), dw_histo[Plotter.mc_samples].integrate('process'), error_opts=Plotter.hstyles['Down'],
                            unc='num', clear=False, ax=rax, guide_opts={},
                        )
                        rax.set_ylabel(rax_ylabel)
                        rax.set_xlabel(xtitle)
                    #set_trace()
                    #ax, rax = Plotter.comp_variations(ax, rax, nom=nom_histo, up=up_histo, dw=dw_histo, xlabel=xtitle, xlimits=x_lims, ylabel='Events', **opts)
        
                        # add lepton/jet multiplicity label
                    ax.text(
                        0.02, 0.92, "%s, %s" % (leptypes[args.lepton], jet_mults[jmult]),
                        fontsize=rcParams['font.size']*0.9, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
                    )
                        ## draw vertical lines for distinguishing different ctstar bins
                    if vlines is not None:
                        for vline in vlines:
                            ax.axvline(vline, color='k', linestyle='--')
                            if rax is not None: rax.axvline(vline, color='k', linestyle='--')
                    hep.cms.cmslabel(ax=ax, data=maskData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
                    
                    #set_trace()
                    figname = '%s/%s' % (pltdir, '_'.join([jmult, args.lepton, sys, hname]))
                    fig.savefig(figname)
                    print('%s written' % figname)
                    plt.close()


                ## make variation/nominal ratio plots
                if args.ratio:
                    up_mc = up_histo[Plotter.mc_samples].integrate('process')
                    dw_mc = dw_histo[Plotter.mc_samples].integrate('process')

                    fig, ax = plt.subplots()
                    fig.subplots_adjust(hspace=.07)

                    #set_trace()
                        ## find first and last valid bins for up variation
                    up_first_valid_bin, up_last_valid_bin = np.where(~np.isnan(up_mc.values()[()]/nom_mc.values()[()]))[0][0], np.where(~np.isnan(up_mc.values()[()]/nom_mc.values()[()]))[0][-1]+1
                    up_masked_vals = np.ma.masked_where(np.isnan((up_mc.values()[()]/nom_mc.values()[()])[up_first_valid_bin:up_last_valid_bin]), (up_mc.values()[()]/nom_mc.values()[()])[up_first_valid_bin:up_last_valid_bin])
                        ## find first and last valid bins for down variation
                    dw_first_valid_bin, dw_last_valid_bin = np.where(~np.isnan(dw_mc.values()[()]/nom_mc.values()[()]))[0][0], np.where(~np.isnan(dw_mc.values()[()]/nom_mc.values()[()]))[0][-1]+1
                    dw_masked_vals = np.ma.masked_where(np.isnan((dw_mc.values()[()]/nom_mc.values()[()])[dw_first_valid_bin:dw_last_valid_bin]), (dw_mc.values()[()]/nom_mc.values()[()])[dw_first_valid_bin:dw_last_valid_bin])
                        ## plot up/down variations
                            ## bins = edges     vals = np.r_[vals, vals[-1]]
                    ax.step(nom_mc.dense_axes()[0].edges()[up_first_valid_bin:up_last_valid_bin+1], np.r_[up_masked_vals, up_masked_vals[-1]], where='post', **Plotter.hstyles['Up'])
                    ax.step(nom_mc.dense_axes()[0].edges()[dw_first_valid_bin:dw_last_valid_bin+1], np.r_[dw_masked_vals, dw_masked_vals[-1]], where='post', **Plotter.hstyles['Down'])

                    ax.axhline(1, **{'linestyle': '--', 'color': (0, 0, 0, 0.5), 'linewidth': 1})
                    ax.legend(loc='upper right', title=sys)
                    ax.autoscale(axis='x', tight=True)
                    ax.set_xlim(x_lims)
                    ax.set_xlabel(xtitle)
                    ax.set_ylabel('Ratio to Nominal')
                        # find y limits using buffer of 1.5 times max deviation from variations
                    max_rat_val = max(max(abs((up_mc.values()[()]/nom_mc.values()[()])[np.isfinite(up_mc.values()[()]/nom_mc.values()[()])]-1)), max(abs((dw_mc.values()[()]/nom_mc.values()[()])[np.isfinite(dw_mc.values()[()]/nom_mc.values()[()])]-1)))
                    max_buffer = round(max_rat_val*1.5, 2)
                    ax.set_ylim(1-max_buffer, 1+max_buffer)
        
                        # add lepton/jet multiplicity label
                    ax.text(
                        0.02, 0.94, "%s, %s" % (leptypes[args.lepton], jet_mults[jmult]),
                        fontsize=rcParams['font.size']*0.9, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
                    )
                        ## draw vertical lines for distinguishing different ctstar bins
                    if vlines is not None:
                        for vline in vlines:
                            ax.axvline(vline, color='k', linestyle='--')
                    hep.cms.cmslabel(ax=ax, data=False, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
                    
                    #set_trace()
                    figname = '%s/%s' % (pltdir, '_'.join([jmult, args.lepton, sys, hname, 'MC_Ratio_to_Nominal']))
                    fig.savefig(figname)
                    print('%s written' % figname)
                    plt.close()


