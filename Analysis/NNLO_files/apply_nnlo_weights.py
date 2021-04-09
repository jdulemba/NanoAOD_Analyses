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

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
base_jobid = os.environ['base_jobid']
analyzer = 'apply_nnlo_weights'

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('year', choices=['2016APV', '2016', '2017', '2018'] if base_jobid == 'ULnanoAOD' else ['2016', '2017', '2018'], help='Specify which year to run over')
parser.add_argument('lepton', choices=['Electron', 'Muon'], help='Choose which lepton to make plots for')
args = parser.parse_args()


rewt_style_dict = {
    'Nominal' : ('Nominal', 'k'),
    'thad_pt' : ('$W_{Orig}$($p_{T}$($t_{h}$))', '#e42a2c'), ## red
    'mtt_vs_thad_ctstar' : ('$W_{Orig}$($m_{t\\bar{t}}$, cos($\\theta^{*}_{t_{h}}$))', '#377eb8'), ## blue
    'thad_pt_Interp' : ('$W_{Int}$($p_{T}$($t_{h}$))', '#4daf4a'), ## green
    'mtt_vs_thad_ctstar_Interp' : ('$W_{Int}$($m_{t\\bar{t}}$, cos($\\theta^{*}_{t_{h}}$))', '#ff7f00'), ## orange
}

input_dir = os.path.join(proj_dir, 'results', '%s_%s' % (args.year, jobid), analyzer)
f_ext = 'TOT.coffea'
outdir = os.path.join(proj_dir, 'plots', '%s_%s' % (args.year, jobid), analyzer)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

fnames = [os.path.join(input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)]
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

btag_cats = {
    'btagFail' : '0 btags',
    'btagPass' : '$n_{btags} \geq$ 2',
}

lep_cats = {
    'Tight' : 'tight %s' % objtypes['Lep'][args.lepton],
    'Loose' : 'loose %s' % objtypes['Lep'][args.lepton],
}

evt_type_cats = {
    'RECO' : 'detector level',
    'GEN' : 'parton level',
    'RESO' : '',
}

perm_styles = {
    'ttJetsSL_right' : {
        'name' : "$t\\bart$ SL correct",
    },
    'ttJetsSL_matchable' : {
        'name' : "$t\\bart$ SL matchable",
    },
    'ttJetsSL_unmatchable' : {
        'name' : "$t\\bart$ SL unmatchable",
    },
    'ttJetsSL_sl_tau' : {
        'name' : "$t\\bart$ SL $\\tau$",
    },
    'ttJetsDiLep' : {
        'name' : "$t\\bart$ DL",
    },
    'ttJetsHad' : {
        'name' : "$t\\bart$ FH",
    },
}

linearize_binning = (
    np.array([300.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 625.0, 650.0, 675.0, 700.0, 730.0, 760.0, 800.0, 850.0, 900.0, 1000.0, 1200.0]),
    #np.array([300.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 620.0, 650., 700.0, 750.0, 800.0, 850.0, 900.0, 2000.0]),
    np.array([0.0, 0.4, 0.6, 0.75, 0.9, 1.0])
)

variables = {
    'mtt' : ('m($t\\bar{t}$) [GeV]', 4, (200., 2000.)),
    'mthad' : ('m($t_{h}$) [GeV]', 2, (0., 300.)),
    'pt_thad' : ('$p_{T}$($t_{h}$) [GeV]', 2, (0., 500.)),
    'pt_tlep' : ('$p_{T}$($t_{l}$) [GeV]', 2, (0., 500.)),
    'pt_tt' : ('$p_{T}$($t\\bar{t}$) [GeV]', 2, (0., 500.)),
    'eta_thad' : ('$\\eta$($t_{h}$)', 2, (-4., 4.)),
    'eta_tlep' : ('$\\eta$($t_{l}$)', 2, (-4., 4.)),
    'eta_tt' : ('$\\eta$($t\\bar{t}$)', 2, (-4., 4.)),
    'tlep_ctstar' : ('cos($\\theta^{*}_{t_{l}}$)', 2, (-1., 1.)),
    'tlep_ctstar_abs' : ('|cos($\\theta^{*}_{t_{l}}$)|', 1, (0., 1.)),
    ##'mtt_vs_tlep_ctstar_abs' : ('m($t\\bar{t}$)', '|cos($\\theta^{*}_{t_{l}}$)|', linearize_binning[0], linearize_binning[1], (linearize_binning[0][0], linearize_binning[0][-1]), (linearize_binning[1][0], linearize_binning[1][-1]), True),
}
reso_variables = {
    'Reso_mtt' : ('m($t\\bar{t}$) Resolution (Gen-Reco) [GeV]', 1, (-200., 200.)),
    'Reso_mthad' : ('m($t_{h}$) Resolution (Gen-Reco) [GeV]', 2, (0., 300.)),
    'Reso_pt_thad' : ('$p_{T}$($t_{h}$) Resolution (Gen-Reco) [GeV]', 2, (0., 500.)),
    'Reso_pt_tlep' : ('$p_{T}$($t_{l}$) Resolution (Gen-Reco) [GeV]', 2, (0., 500.)),
    'Reso_pt_tt' : ('$p_{T}$($t\\bar{t}$) Resolution (Gen-Reco) [GeV]', 2, (0., 500.)),
    'Reso_eta_thad' : ('$\\eta$($t_{h}$) Resolution (Gen-Reco)', 2, (-4., 4.)),
    'Reso_eta_tlep' : ('$\\eta$($t_{l}$) Resolution (Gen-Reco)', 2, (-4., 4.)),
    'Reso_eta_tt' : ('$\\eta$($t\\bar{t}$) Resolution (Gen-Reco)', 2, (-4., 4.)),
    'Reso_tlep_ctstar' : ('cos($\\theta^{*}_{t_{l}}$) Resolution (Gen-Reco)', 2, (-1., 1.)),
    'Reso_tlep_ctstar_abs' : ('|cos($\\theta^{*}_{t_{l}}$)| Resolution (Gen-Reco)', 1, (0., 1.)),
    ###'Reso_mtt_vs_tlep_ctstar_abs' : ('m($t\\bar{t}$)', '|cos($\\theta^{*}_{t_{l}}$)|', linearize_binning[0], linearize_binning[1], (linearize_binning[0][0], linearize_binning[0][-1]), (linearize_binning[1][0], linearize_binning[1][-1]), True),
}


    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, 'inputs', '%s_lumis_data.json' % base_jobid)).read())[args.year]
lumi_correction = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'MC_LumiWeights.coffea'))[args.year]['%ss' % args.lepton]
        # scale ttJets events, split by reconstruction type, by normal ttJets lumi correction
ttJets_permcats = ['*right', '*matchable', '*unmatchable', '*sl_tau', '*other']
names = [dataset for dataset in sorted(set([key[0] for key in hdict['mtt'].values().keys()]))] # 'mtt' hardcoded because it has all ttJets event cats
ttJets_cats = [name for name in names if any([fnmatch.fnmatch(name, cat) for cat in ttJets_permcats])] # gets ttJets(_PS)_other, ...
if len(ttJets_cats) > 0:
    for tt_cat in ttJets_cats:
        ttJets_lumi_topo = '_'.join(tt_cat.split('_')[:-2]) if 'sl_tau' in tt_cat else '_'.join(tt_cat.split('_')[:-1]) # gets ttJets[SL, Had, DiLep] or ttJets_PS
        ttJets_eff_lumi = lumi_correction[ttJets_lumi_topo]
        lumi_correction.update({tt_cat: ttJets_eff_lumi})

## scale by lumi and make groups based on process
process = hist.Cat("process", "Process", sorting='placement')
process_cat = "dataset"
process_groups = {'ttJetsDiLep' : ['ttJetsDiLep_other'], 'ttJetsHad' : ['ttJetsHad_other'], 'ttJetsSL_matchable' : ['ttJetsSL_matchable'], 'ttJetsSL_right' : ['ttJetsSL_right'], 'ttJetsSL_sl_tau' : ['ttJetsSL_sl_tau'], 'ttJetsSL_unmatchable' : ['ttJetsSL_unmatchable']}
for hname in hdict.keys():
    if 'cutflow' in hname: continue
    #set_trace()
    hdict[hname].scale(lumi_correction, axis='dataset')
    hdict[hname] = hdict[hname].group(process_cat, process, process_groups)


    ## make plots
#for plt_type in ['RECO']:
#for plt_type in ['GEN']:
for plt_type in ['RECO', 'GEN', 'RESO']:
    vars_dict = reso_variables if plt_type == 'RESO' else variables
    for hname in vars_dict.keys():
        if not hname in hdict.keys(): # check that hist is in hdict
            raise ValueError("%s not found in file" % hname)
        #set_trace()
        
        histo = hdict[hname][:, :, :, args.lepton, :, :].integrate('leptype') if plt_type == 'RESO' else hdict[hname][:, :, plt_type, :, args.lepton, :, :].integrate('evt_type').integrate('leptype') # process, evt_type, jmult, leptype, btag, lepcat
    
        if histo.dense_dim() == 1:
            xtitle, rebinning, x_lims = vars_dict[hname]
            xaxis_name = histo.dense_axes()[0].name
            if rebinning != 1:
                histo = histo.rebin(xaxis_name, rebinning)
    
            #set_trace()
            # hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
            #for jmult in ['4PJets']:
            #    for lepcat in ['Tight']:
            #        for btagregion in ['btagPass']:
            for jmult in sorted(set([key[2] for key in histo.values().keys()])):
                for lepcat in sorted(set([key[4] for key in histo.values().keys()])):
                    for btagregion in sorted(set([key[3] for key in histo.values().keys()])):
                        pltdir = os.path.join(outdir, args.lepton, jmult, lepcat, btagregion, plt_type)
                        if not os.path.isdir(pltdir):
                            os.makedirs(pltdir)
    
                        print(', '.join([jmult, lepcat, btagregion, hname])) 
                        #set_trace()

                            # inclusive ttbar categories
                        allTT_hslice = histo[:, :, jmult, btagregion, lepcat].integrate('jmult').integrate('lepcat').integrate('btag').integrate('process') # use SL+DL+Had events

                        fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                        fig.subplots_adjust(hspace=.07)
    
                        for rewt in sorted(set([key[0] for key in allTT_hslice.values().keys()])):
                            vals, bins = (allTT_hslice[rewt].integrate('rewt')).values()[()], allTT_hslice.axis(xaxis_name).edges()
                            ax = Plotter.plot_1D(vals, bins, xlimits=x_lims, ax=ax, histtype='step', label=rewt_style_dict[rewt][0], color=rewt_style_dict[rewt][1])
    
                            if (rewt != 'Nominal'):
                                ratio_vals, ratio_bins = Plotter.get_ratio_arrays(num_vals=vals, denom_vals=(allTT_hslice['Nominal'].integrate('rewt')).values()[()], input_bins=bins)
                                rax.step(ratio_bins, ratio_vals, where='post', **{'linestyle' : '-', 'color' : rewt_style_dict[rewt][1]})
                                
                        ax.autoscale(axis='x', tight=True)
                        rax.autoscale(axis='x', tight=True)
                        ax.set_ylim(None)
                        rax.axhline(1, **{'linestyle': '--', 'color': (0, 0, 0, 0.5), 'linewidth': 1})
                        rax.set_ylabel('W(var)/Nominal')
                        rax.set_xlabel(xtitle)
                        ax.legend(title='All $t\\bart$', loc='upper right')
                            # add lepton/jet multiplicity label
                        #set_trace()
                        ax.text(
                            0.02, 0.88, "%s, %s\n%s\t\t%s" % (lep_cats[lepcat], jet_mults[jmult], btag_cats[btagregion], evt_type_cats[plt_type]),
                            fontsize=rcParams['font.size']*0.75, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
                        )
                        hep.cms.label(ax=ax, data=False, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
    
                        #set_trace()
                        figname = os.path.join(pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, 'ttJetsAll', plt_type, hname]))
                        fig.savefig(figname)
                        print('%s written' % figname)
                        plt.close()

    
                            # individual ttbar categories
                        for ttcat in process_groups.keys():
                            #set_trace()
                            if ttcat not in sorted(set([key[0] for key in histo.values().keys()])): continue
                            TTcat_hslice = histo[ttcat, :, jmult, btagregion, lepcat].integrate('jmult').integrate('lepcat').integrate('btag').integrate('process')

                            fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                            fig.subplots_adjust(hspace=.07)
    
                            for rewt in sorted(set([key[0] for key in TTcat_hslice.values().keys()])):
                                vals, bins = (TTcat_hslice[rewt].integrate('rewt')).values()[()], TTcat_hslice.axis(xaxis_name).edges()
                                ax = Plotter.plot_1D(vals, bins, xlimits=x_lims, ax=ax, histtype='step', label=rewt_style_dict[rewt][0], color=rewt_style_dict[rewt][1])
    
                                if (rewt != 'Nominal'):
                                    ratio_vals, ratio_bins = Plotter.get_ratio_arrays(num_vals=vals, denom_vals=(TTcat_hslice['Nominal'].integrate('rewt')).values()[()], input_bins=bins)
                                    rax.step(ratio_bins, ratio_vals, where='post', **{'linestyle' : '-', 'color' : rewt_style_dict[rewt][1]})
                                    
                            ax.autoscale(axis='x', tight=True)
                            rax.autoscale(axis='x', tight=True)
                            ax.set_ylim(None)
                            rax.axhline(1, **{'linestyle': '--', 'color': (0, 0, 0, 0.5), 'linewidth': 1})
                            rax.set_ylabel('W(var)/Nominal')
                            rax.set_xlabel(xtitle)
                            ax.legend(title=perm_styles[ttcat]['name'], loc='upper right')
                                # add lepton/jet multiplicity label
                            #set_trace()
                            ax.text(
                                0.02, 0.88, "%s, %s\n%s\t\t%s" % (lep_cats[lepcat], jet_mults[jmult], btag_cats[btagregion], evt_type_cats[plt_type]),
                                fontsize=rcParams['font.size']*0.75, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
                            )
                            hep.cms.label(ax=ax, data=False, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
    
                            #set_trace()
                            figname = os.path.join(pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, ttcat, plt_type, hname]))
                            fig.savefig(figname)
                            print('%s written' % figname)
                            plt.close()
    
        
        
        if histo.dense_dim() == 2:
            set_trace()
            xtitle, ytitle, xrebinning, yrebinning, x_lims, y_lims, withData = vars_dict[hname]
    
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
                            ax = hep.cms.label(ax=ax, data=withData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
    
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
                            ax = hep.cms.label(ax=ax, data=withData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
    
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
                        ax = hep.cms.label(ax=ax, data=withData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
    
                        #set_trace()
                        figname = '%s/%s' % (pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, hname]))
                        fig.savefig(figname)
                        print('%s written' % figname)
                        plt.close()
    
    
