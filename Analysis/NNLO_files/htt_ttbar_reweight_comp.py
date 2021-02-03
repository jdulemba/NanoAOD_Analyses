#!/usr/bin/env python

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
import re
from coffea import hist
import numpy as np
import fnmatch
import Utilities.Plotter as Plotter

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('year', choices=['2016', '2017', '2018'], help='What year is the ntuple from.')
parser.add_argument('lepton', choices=['Electron', 'Muon'], help='Choose which lepton to make plots for')
args = parser.parse_args()


proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
base_jobid = os.environ['base_jobid']
analyzer = 'htt_btag_iso_cut'

rewt_style_dict = {
    'jpt30_ljpt50_MT40_cutBasedEl' : ('Nominal', 'k'),
    'pt_thad_reweighting' : ('$W_{Orig}$($p_{T}$($t_{h}$))', '#e42a2c'), ## red
    'mtt_thadctstar_reweighting' : ('$W_{Orig}$($m_{t\\bar{t}}$, cos($\\theta^{*}_{t_{h}}$))', '#377eb8'), ## blue
    ##'thad_pt_Interp_reweighting' : ('$W_{Int}$($p_{T}$($t_{h}$))', '#4daf4a'), ## green
    'mtt_thadctstar_Interp_reweighting' : ('$W_{Int}$($m_{t\\bar{t}}$, cos($\\theta^{*}_{t_{h}}$))', '#ff7f00'), ## orange
}

nominal_jobid = 'jpt30_ljpt50_MT40_cutBasedEl'

input_dirs = {key: os.path.join(proj_dir, 'results', '%s_%s_%s' % (args.year, base_jobid, key), analyzer) for key in rewt_style_dict.keys()}
f_ext = 'TOT.coffea'
outdir = os.path.join(proj_dir, 'plots', base_jobid, 'htt_ttbar_reweight_comp', args.year)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

fnames_dict = {rewt : sorted([os.path.join(dir, fname) for fname in os.listdir(dir) if fname.endswith(f_ext)]) for rewt, dir in input_dirs.items()}

#set_trace()
hdicts = {rewt : plt_tools.add_coffea_files(fnames) if len(fnames) > 1 else load(fnames[0]) for rewt, fnames in fnames_dict.items()}

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

btag_cats = {
    'btagFail' : '0 btags',
    'btagPass' : '$n_{btags} \geq$ 2',
}

lep_cats = {
    'Tight' : 'tight %s' % objtypes['Lep'][args.lepton],
    'Loose' : 'loose %s' % objtypes['Lep'][args.lepton],
}

linearize_binning = (
    np.array([300.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 625.0, 650.0, 675.0, 700.0, 730.0, 760.0, 800.0, 850.0, 900.0, 1000.0, 1200.0]),
    #np.array([300.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 620.0, 650., 700.0, 750.0, 800.0, 850.0, 900.0, 2000.0]),
    np.array([0.0, 0.4, 0.6, 0.75, 0.9, 1.0])
)

variables = {
    #'mtt_vs_tlep_ctstar_abs' : ('m($t\\bar{t}$)', '|cos($\\theta^{*}_{t_{l}}$)|', linearize_binning[0], linearize_binning[1], (linearize_binning[0][0], linearize_binning[0][-1]), (linearize_binning[1][0], linearize_binning[1][-1]), True),
    'Jets_njets' : ('$n_{jets}$', 1, (0, 15)),
    ####'mtt' : ('m($t\\bar{t}$) [GeV]', linearize_binning[0], (200., 2000.), True),
    ####'tlep_ctstar_abs' : ('|cos($\\theta^{*}_{t_{l}}$)|', linearize_binning[1], (0., 1.), True),
    'mtt' : ('m($t\\bar{t}$) [GeV]', 4, (200., 2000.)),
    'tlep_ctstar_abs' : ('|cos($\\theta^{*}_{t_{l}}$)|', 1, (0., 1.)),
    'mthad' : ('m($t_{h}$) [GeV]', 2, (0., 300.)),
    'pt_thad' : ('$p_{T}$($t_{h}$) [GeV]', 2, (0., 500.)),
    'pt_tlep' : ('$p_{T}$($t_{l}$) [GeV]', 2, (0., 500.)),
    'pt_tt' : ('$p_{T}$($t\\bar{t}$) [GeV]', 2, (0., 500.)),
    'eta_thad' : ('$\\eta$($t_{h}$)', 2, (-4., 4.)),
    'eta_tlep' : ('$\\eta$($t_{l}$)', 2, (-4., 4.)),
    'eta_tt' : ('$\\eta$($t\\bar{t}$)', 2, (-4., 4.)),
    'tlep_ctstar' : ('cos($\\theta^{*}_{t_{l}}$)', 2, (-1., 1.)),
    'full_disc' : ('$\\lambda_{C}$', 2, (5, 25.)),
    'mass_disc' : ('$\\lambda_{M}$', 2, (0, 20.)),
    'ns_disc' : ('$\\lambda_{NS}$', 2, (3., 10.)),
    'Jets_pt' : ('$p_{T}$(jets) [GeV]', 2, (0., 300.)),
    'Jets_eta' : ('$\\eta$(jets)', 2, (-2.6, 2.6)),
    'Jets_LeadJet_pt' : ('$p_{T}$(leading jet) [GeV]', 2, (0., 300.)),
    'Jets_LeadJet_eta' : ('$\\eta$(leading jet)', 2, (-2.6, 2.6)),
    'Lep_pt' : ('$p_{T}$(%s) [GeV]' % objtypes['Lep'][args.lepton], 2, (0., 300.)),
    'Lep_eta' : ('$\\eta$(%s)' % objtypes['Lep'][args.lepton], 2, (-2.6, 2.6)),
    'Lep_iso' : ('pfRelIso, %s' % objtypes['Lep'][args.lepton], 1, (0., 1.)),
    'MT' : ('$M_{T}$ [GeV]', 1, (0., 300.)),
    'MET_pt' : ('$p_{T}$(MET) [GeV]', 1, (0., 300.)),
    'MET_phi' : ('$\phi$(MET)', 1, (-3.2, 3.2)),
}



    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, 'inputs', 'lumis_data.json')).read())[args.year]
lumi_correction = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'MC_LumiWeights_Test.coffea'))[args.year]['%ss' % args.lepton]
#lumi_correction = load(os.path.join(proj_dir, 'Corrections', jobid, 'MC_LumiWeights_allTTJets.coffea'))[args.year]['%ss' % args.lepton]
        # scale ttJets events, split by reconstruction type, by normal ttJets lumi correction
ttJets_permcats = ['*right', '*matchable', '*unmatchable', '*sl_tau', '*other']
names_list = [[dataset for dataset in sorted(set([key[0] for key in hdict[sorted(variables.keys())[0]].values().keys()]))] for hdict in hdicts.values()]
names = sorted(set(sum(names_list, []))) # get dataset names in hists
ttJets_cats = [name for name in names if any([fnmatch.fnmatch(name, cat) for cat in ttJets_permcats])] # gets ttJets(_PS)_other, ...
if len(ttJets_cats) > 0:
    for tt_cat in ttJets_cats:
        ttJets_lumi_topo = '_'.join(tt_cat.split('_')[:-2]) if 'sl_tau' in tt_cat else '_'.join(tt_cat.split('_')[:-1]) # gets ttJets[SL, Had, DiLep] or ttJets_PS
        ttJets_eff_lumi = lumi_correction[ttJets_lumi_topo]
        lumi_correction.update({tt_cat: ttJets_eff_lumi})

## scale by lumi and make groups based on process
process = hist.Cat("process", "Process", sorting='placement')
process_cat = "dataset"
process_groups = plt_tools.make_dataset_groups(args.lepton, args.year, samples=names)
for hdict in hdicts.values():
    for hname in hdict.keys():
        if 'cutflow' in hname: continue
        #set_trace()
        hdict[hname].scale(lumi_correction, axis='dataset')
        hdict[hname] = hdict[hname].group(process_cat, process, process_groups)
        hdict[hname] = hdict[hname].group('process', hist.Cat("process", "Process", sorting='placement'), {'ttJets' : ['ttJets_right', 'ttJets_matchable', 'ttJets_unmatchable', 'ttJets_sl_tau', 'ttJets_other']}) # only keeps ttbar


    ## make plots
for hname in variables.keys():
    if not all([hname in hdict.keys() for hdict in hdicts.values()]): # check that hist is in all hdicts
        raise ValueError("%s not found in file" % hname)
    histos = {rewt: hdict[hname][:, 'nosys', :, args.lepton, :, :].integrate('sys').integrate('leptype') for rewt, hdict in hdicts.items()} # process, sys, jmult, leptype, btag, lepcat

    if (list(histos.values())[0]).dense_dim() == 1:
        xtitle, rebinning, x_lims = variables[hname]
        xaxis_name = (list(histos.values())[0]).dense_axes()[0].name
        if rebinning != 1:
            histos = {rewt: histo.rebin(xaxis_name, rebinning) for rewt, histo in histos.items()}

        #set_trace()
        # hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
        #for jmult in ['4PJets']:
        #    for lepcat in ['Tight']:
        #        for btagregion in ['btagPass']:
        for jmult in sorted(set([key[1] for key in (list(histos.values())[0]).values().keys()])):
            for lepcat in sorted(set([key[3] for key in (list(histos.values())[0]).values().keys()])):
                for btagregion in sorted(set([key[2] for key in (list(histos.values())[0]).values().keys()])):
                    pltdir = '/'.join([outdir, args.lepton, jmult, lepcat, btagregion])
                    if not os.path.isdir(pltdir):
                        os.makedirs(pltdir)

                    print(', '.join([jmult, lepcat, btagregion, hname])) 
                    fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                    fig.subplots_adjust(hspace=.07)

                    if hname == 'Lep_iso':
                        if args.lepton == 'Muon':
                            x_lims = (0., 0.15) if lepcat == 'Tight' else (0.15, 1.)
                        if args.lepton == 'Electron':
                            x_lims = (0., 0.1) if lepcat == 'Tight' else (0., 0.5)

                    if hname == 'Jets_njets':
                        rows = [("Lumi: %s fb^-1" % format(data_lumi_year['%ss' % args.lepton]/1000., '.1f'), "%s, %s" % (args.lepton, jmult), "%s, %s" % (lepcat, btagregion))]
                        rows += [("Reweighting", "Yield", "Yield/Nominal")]

                    #set_trace()
                    hslices = {rewt : histo['ttJets', jmult, btagregion, lepcat].integrate('jmult').integrate('lepcat').integrate('btag').integrate('process') for rewt, histo in histos.items()}
                    logy_min, logy_max = 1e-1, 0.
                    for rewt, hslice in hslices.items():
                        vals, bins = hslice.values()[()], hslice.axis(xaxis_name).edges()
                        ax = Plotter.plot_1D(vals, bins, xlimits=x_lims, ytitle='$t\\bart$ Events', ax=ax, histtype='step', label=rewt_style_dict[rewt][0], color=rewt_style_dict[rewt][1])

                        if rewt != nominal_jobid:
                            ratio_vals, ratio_bins = Plotter.get_ratio_arrays(num_vals=vals, denom_vals=hslices[nominal_jobid].values()[()], input_bins=bins)
                            rax.step(ratio_bins, ratio_vals, where='post', **{'linestyle' : '-', 'color' : rewt_style_dict[rewt][1]})

                        if np.min(vals[vals > 0.]) < logy_min: logy_min = np.min(vals[vals > 0.])
                        if np.max(vals) > logy_max: logy_max = np.max(vals)

                        if hname == 'Jets_njets':
                            rows += [(rewt.split('_reweighting')[0], format(hslice.values(overflow='all')[()].sum(), '.2f'), format(hslice.values(overflow='all')[()].sum()/(hslices[nominal_jobid].values(overflow='all')[()].sum()), '.6f') )]
                            
                    if hname == 'Jets_njets':
                        yields_fname = os.path.join(pltdir, '%s_ttJets_reweighting_yields.txt' % '_'.join([jmult, args.lepton, lepcat, btagregion]))
                        plt_tools.print_table(rows, filename=yields_fname, print_output=True, header_line=1)
                        print('\n%s written' % yields_fname)

                    ax.set_yscale('log')
                    ax.autoscale(axis='x', tight=True)
                    rax.autoscale(axis='x', tight=True)
                    ax.set_ylim(logy_min, 15.*logy_max)
                    rax.axhline(1, **{'linestyle': '--', 'color': (0, 0, 0, 0.5), 'linewidth': 1})
                    rax.set_ylabel('W(var)/Nominal')
                    rax.set_xlabel(xtitle)
                    ax.legend(loc='upper right')
                        # add lepton/jet multiplicity label
                    #set_trace()
                    ax.text(
                        0.02, 0.88, "%s, %s\n%s" % (lep_cats[lepcat], jet_mults[jmult], btag_cats[btagregion]),
                        fontsize=rcParams['font.size']*0.75, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
                    )
                    hep.cms.cmslabel(ax=ax, data=False, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))

                    #set_trace()
                    figname = os.path.join(pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, hname, 'ttJets']))
                    fig.savefig(figname)
                    print('%s written' % figname)
                    plt.close()

    
    if (list(histos.values())[0]).dense_dim() == 2:
    #if histo.dense_dim() == 2:
        xtitle, ytitle, xrebinning, yrebinning, x_lims, y_lims, withData = variables[hname]
        set_trace()

        xaxis_name = histo.dense_axes()[0].name
        yaxis_name = histo.dense_axes()[1].name
            ## rebin x axis
        if isinstance(xrebinning, np.ndarray):
            new_xbins = hist.Bin(xaxis_name, xaxis_name, xrebinning)
        elif isinstance(xrebinning, float) or isinstance(xrebinning, int):
            new_xbins = xrebinning
        histo = histo.rebin(xaxis_name, new_xbins)
            ## rebin y axis
        if isinstance(yrebinning, np.ndarray):
            new_ybins = hist.Bin(yaxis_name, yaxis_name, yrebinning)
        elif isinstance(yrebinning, float) or isinstance(yrebinning, int):
            new_ybins = yrebinning
        histo = histo.rebin(yaxis_name, new_ybins)

        ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
        for jmult in sorted(set([key[1] for key in histo.values().keys()])):
            #for lepcat in ['Tight']:
            #    for btagregion in ['btagPass']:
            for lepcat in sorted(set([key[3] for key in histo.values().keys()])):
                for btagregion in sorted(set([key[2] for key in histo.values().keys()])):
                    pltdir = '/'.join([outdir, args.lepton, jmult, lepcat, btagregion])
                    if not os.path.isdir(pltdir):
                        os.makedirs(pltdir)

                    if (lepcat, btagregion) == ('Tight', 'btagPass'):
                        withData = False
                    else:
                        withData = variables[hname][-1]
                    print(', '.join([jmult, lepcat, btagregion, hname])) 
                    hslice = histo[:, jmult, btagregion, lepcat].integrate('jmult').integrate('lepcat').integrate('btag')

                    mc_opts = {
                        'maskData' : not withData
                    #    'mcorder' : ['QCD', 'EWK', 'singlet', 'ttJets'] if not ttJets_cats else ['QCD', 'EWK', 'singlet', 'ttJets_other', 'ttJets_unmatchable', 'ttJets_matchable', 'ttJets_right']
                    }

                        # make 1D projection along dense axes
                    for dax in range(2):
                        if dax == 0:
                            xlabel = ytitle
                            xlimits = y_lims
                            figname = '%s/%s' % (pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, yaxis_name, 'proj']))
                        else:
                            xlabel = '%s [GeV]' % xtitle
                            xlimits = x_lims
                            figname = '%s/%s' % (pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, xaxis_name, 'proj']))

                        fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                        fig.subplots_adjust(hspace=.07)

                        #set_trace()
                        hproj = hslice.integrate(hslice.dense_axes()[dax].name)

                        ax, rax = Plotter.plot_stack1d(ax, rax, hproj, xlabel=xlabel, xlimits=xlimits, **mc_opts)

                            # add lepton/jet multiplicity label
                        #set_trace()
                        ax.text(
                            0.02, 0.88, "%s, %s\n%s" % (lep_cats[lepcat], jet_mults[jmult], btag_cats[btagregion]),
                            fontsize=rcParams['font.size']*0.75, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
                        )
                        ax = hep.cms.cmslabel(ax=ax, data=withData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))

                        #set_trace()
                        fig.savefig(figname)
                        print('%s written' % figname)
                        plt.close()


                        ## make 1D plots of mtt for each ctstar bin
                    for ybin in range(len(hslice.dense_axes()[1].edges())-1):
                        fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                        fig.subplots_adjust(hspace=.07)

                        hproj = hslice[:, :, hslice.dense_axes()[1].edges()[ybin]:hslice.dense_axes()[1].edges()[ybin+1]].integrate(hslice.dense_axes()[1].name)

                        ax, rax = Plotter.plot_stack1d(ax, rax, hproj, xlabel='%s [GeV]' % xtitle, xlimits=x_lims, **mc_opts)

                            # add lepton/jet multiplicity label, add ctstar range label
                        binlabel = '%s $\leq$ %s < %s' % (hslice.dense_axes()[1].edges()[ybin], ytitle, hslice.dense_axes()[1].edges()[ybin+1])
                        ax.text(
                            0.02, 0.88, "%s, %s\t%s\n%s" % (lep_cats[lepcat], jet_mults[jmult], binlabel, btag_cats[btagregion]),
                            fontsize=rcParams['font.size']*0.75, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
                        )
                        ax = hep.cms.cmslabel(ax=ax, data=withData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))

                        #set_trace()
                        bintitle = '%sctstar%s' % (hslice.dense_axes()[1].edges()[ybin], hslice.dense_axes()[1].edges()[ybin+1])
                        bintitle = bintitle.replace('.', 'p')
                        figname = '%s/%s' % (pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, bintitle, 'mtt']))
                        fig.savefig(figname)
                        print('%s written' % figname)
                        plt.close()


                        # plot linearized view 
                    fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                    fig.subplots_adjust(hspace=.07)

                    #set_trace()
                    hline = Plotter.linearize_hist(hslice)
                    #set_trace()
                    
                    ax, rax = Plotter.plot_stack1d(ax, rax, hline, xlabel='%s $\otimes$ %s' % (xtitle, ytitle), xlimits=(0, len(hline.axis(hline.dense_axes()[0].name).edges())-1), **mc_opts)

                    #set_trace()
                        # draw vertical lines separating ctstar bins
                    bin_sep_lines = [hslice.values()[('EWK',)].shape[0]*ybin for ybin in range(1, hslice.values()[('EWK',)].shape[1])]
                    for binline in bin_sep_lines:
                        ax.axvline(binline, color='k', linestyle='--')
                        if rax is not None: rax.axvline(binline, color='k', linestyle='--')
                        # add lepton/jet multiplicity label
                    ax.text(
                        0.02, 0.88, "%s, %s\n%s" % (lep_cats[lepcat], jet_mults[jmult], btag_cats[btagregion]),
                        fontsize=rcParams['font.size']*0.75, 
                        horizontalalignment='left', 
                        verticalalignment='bottom', 
                        transform=ax.transAxes
                    )
                    ax = hep.cms.cmslabel(ax=ax, data=withData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))

                    #set_trace()
                    figname = '%s/%s' % (pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, hname]))
                    fig.savefig(figname)
                    print('%s written' % figname)
                    plt.close()



 
