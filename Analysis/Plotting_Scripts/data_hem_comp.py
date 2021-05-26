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
parser.add_argument('year', choices=['2018'], help='What year is the ntuple from.')
parser.add_argument('lepton', choices=['Electron', 'Muon'], help='Choose which lepton to make plots for')
args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
base_jobid = os.environ['base_jobid']
analyzer = 'data_hem_comp'

input_dir = os.path.join(proj_dir, 'results', '%s_%s' % (args.year, jobid), analyzer)
f_ext = 'TOT.coffea'
outdir = os.path.join(proj_dir, 'plots', '%s_%s' % (args.year, jobid), analyzer)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

fnames = sorted(['%s/%s' % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)])
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

# binning and bin labels for mtt x ctstar
linearize_binning = (
    np.array([300.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 625.0, 650.0, 675.0, 700.0, 730.0, 760.0, 800.0, 850.0, 900.0, 1000.0, 1200.0]),
    #np.array([300.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 620.0, 650., 700.0, 750.0, 800.0, 850.0, 900.0, 2000.0]),
    np.array([0.0, 0.4, 0.6, 0.75, 0.9, 1.0])
)
#set_trace()
ctstar_binlabels = ['%s $\leq$ |cos($\\theta^{*}_{t_{l}}$)| $\leq$ %s' % (linearize_binning[1][bin], linearize_binning[1][bin+1]) for bin in range(len(linearize_binning[1])-1)]
ctstar_bin_locs = np.linspace((len(linearize_binning[0])-1)/2, (len(linearize_binning[0])-1)*(len(linearize_binning[1])-1) - (len(linearize_binning[0])-1)/2, len(linearize_binning[1])-1)
mtt_binlabels = ['%s $\leq$ m($t\\bar{t}$) $\leq$ %s' % (linearize_binning[0][bin], linearize_binning[0][bin+1]) for bin in range(len(linearize_binning[0])-1)]*len(ctstar_binlabels)

# binning and bin labels for phi x eta
phi_eta_binning = (
    np.array([-3.2, -2.4, -1.6, -0.85, 0., 0.85, 1.6, 2.4, 3.2]), # phi binning
    np.array([-2.5, -1.3, -0.7, 0., 0.7, 1.3, 2.5]) # eta binning
)
eta_binlabels = ['%s $\leq$ $\\eta$ $\leq$ %s' % (phi_eta_binning[1][bin], phi_eta_binning[1][bin+1]) for bin in range(len(phi_eta_binning[1])-1)]
eta_bin_locs = np.linspace((len(phi_eta_binning[0])-1)/2, (len(phi_eta_binning[0])-1)*(len(phi_eta_binning[1])-1) - (len(phi_eta_binning[0])-1)/2, len(phi_eta_binning[1])-1)
phi_binlabels = ['%s $\leq$ $\\phi$ $\leq$ %s' % (phi_eta_binning[0][bin], phi_eta_binning[0][bin+1]) for bin in range(len(phi_eta_binning[0])-1)]*len(eta_binlabels)


variables = {
    'mtt_vs_tlep_ctstar_abs' : ('m($t\\bar{t}$)', '|cos($\\theta^{*}_{t_{l}}$)|', linearize_binning[0], linearize_binning[1], (linearize_binning[0][0], linearize_binning[0][-1]),
        (linearize_binning[1][0], linearize_binning[1][-1]), True, None, (ctstar_binlabels, ctstar_bin_locs)),
    'Jets_phi_vs_eta' : ('$\\phi$(jets)', '$\\eta$(jets)', phi_eta_binning[0], phi_eta_binning[1], (-3.3, 3.3), (-2.6, 2.6), True, phi_binlabels, (eta_binlabels, eta_bin_locs)),
    'Lep_phi_vs_eta' : ('$\\phi$(%s)' % objtypes['Lep'][args.lepton], '$\\eta$(%s)' % objtypes['Lep'][args.lepton], phi_eta_binning[0], phi_eta_binning[1],
        (-3.3, 3.3), (-2.6, 2.6), True, phi_binlabels, (eta_binlabels, eta_bin_locs)),
    #'Jets_njets' : ('$n_{jets}$', 1, (0, 15), True),
    #'mtt' : ('m($t\\bar{t}$) [GeV]', 4, (200., 2000.), True),
    #'tlep_ctstar_abs' : ('|cos($\\theta^{*}_{t_{l}}$)|', 1, (0., 1.), True),
    #'mthad' : ('m($t_{h}$) [GeV]', 2, (0., 300.), True),
    #'pt_thad' : ('$p_{T}$($t_{h}$) [GeV]', 2, (0., 500.), True),
    #'pt_tlep' : ('$p_{T}$($t_{l}$) [GeV]', 2, (0., 500.), True),
    #'pt_tt' : ('$p_{T}$($t\\bar{t}$) [GeV]', 2, (0., 500.), True),
    #'eta_thad' : ('$\\eta$($t_{h}$)', 2, (-4., 4.), True),
    #'eta_tlep' : ('$\\eta$($t_{l}$)', 2, (-4., 4.), True),
    #'eta_tt' : ('$\\eta$($t\\bar{t}$)', 2, (-4., 4.), True),
    #'tlep_ctstar' : ('cos($\\theta^{*}_{t_{l}}$)', 2, (-1., 1.), True),
    #'full_disc' : ('$\\lambda_{C}$', 2, (5, 25.), True),
    #'mass_disc' : ('$\\lambda_{M}$', 2, (0, 20.), True),
    #'ns_disc' : ('$\\lambda_{NS}$', 2, (3., 10.), True),
    #'ns_dist' : ('$D_{\\nu, min}$', 1, (0., 150.), True),
    #'Jets_pt' : ('$p_{T}$(jets) [GeV]', 2, (0., 300.), True),
    #'Jets_eta' : ('$\\eta$(jets)', 2, (-2.6, 2.6), True),
    #'Jets_phi' : ('$\\phi$(jets)', 2, (-4., 4.), True),
    #'Jets_LeadJet_pt' : ('$p_{T}$(leading jet) [GeV]', 2, (0., 300.), True),
    #'Jets_LeadJet_eta' : ('$\\eta$(leading jet)', 2, (-2.6, 2.6), True),
    #'Lep_pt' : ('$p_{T}$(%s) [GeV]' % objtypes['Lep'][args.lepton], 2, (0., 300.), True),
    #'Lep_eta' : ('$\\eta$(%s)' % objtypes['Lep'][args.lepton], 2, (-2.6, 2.6), True),
    #'Lep_phi' : ('$\\phi$(%s)' % objtypes['Lep'][args.lepton], 2, (-4, 4), True),
    #'Lep_iso' : ('pfRelIso, %s' % objtypes['Lep'][args.lepton], 1, (0., 1.), True),
    #'MT' : ('$M_{T}$ [GeV]', 1, (0., 300.), True),
    #'MET_pt' : ('$p_{T}$(MET) [GeV]', 1, (0., 300.), True),
    #'MET_phi' : ('$\phi$(MET)', 1, (-3.2, 3.2), True),
}


    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, 'inputs', '%s_lumis_data.json' % base_jobid)).read())[args.year]
#lumi_correction = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'MC_LumiWeights.coffea'))[args.year]['%ss' % args.lepton]
        # scale ttJets events, split by reconstruction type, by normal ttJets lumi correction
#ttJets_permcats = ['*right', '*matchable', '*unmatchable', '*sl_tau', '*other']
names = [dataset for dataset in sorted(set([key[0] for key in hdict[sorted(variables.keys())[0]].values().keys()]))] # get dataset names in hists
#ttJets_cats = [name for name in names if any([fnmatch.fnmatch(name, cat) for cat in ttJets_permcats])] # gets ttJets(_PS)_other, ...
#if len(ttJets_cats) > 0:
#    for tt_cat in ttJets_cats:
#        ttJets_lumi_topo = '_'.join(tt_cat.split('_')[:-2]) if 'sl_tau' in tt_cat else '_'.join(tt_cat.split('_')[:-1]) # gets ttJets[SL, Had, DiLep] or ttJets_PS
#        ttJets_eff_lumi = lumi_correction[ttJets_lumi_topo]
#        lumi_correction.update({tt_cat: ttJets_eff_lumi})

## make groups based on process
process = hist.Cat("process", "Process", sorting='placement')
process_cat = "dataset"
process_groups = plt_tools.make_dataset_groups(args.lepton, args.year, samples=names)

    # scale and group hists by process
for hname in hdict.keys():
    if 'cutflow' in hname: continue
    #set_trace()
    #hdict[hname].scale(lumi_correction, axis='dataset') # scale hists
    hdict[hname] = hdict[hname].group(process_cat, process, process_groups) # group by process
    hdict[hname] = hdict[hname][:, :, :, args.lepton].integrate('leptype') # only pick out specified lepton



    ## make plots
for hname in variables.keys():
    if hname not in hdict.keys():
        raise ValueError("%s not found in file" % hname)
    #set_trace()
    histo = hdict[hname].integrate('process')

    if histo.dense_dim() == 1:
        xtitle, rebinning, x_lims, withData = variables[hname]
        xaxis_name = histo.dense_axes()[0].name
        if rebinning != 1:
            histo = histo.rebin(xaxis_name, rebinning)

        ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
        for jmult in sorted(set([key[1] for key in histo.values().keys()])):
            #for lepcat in ['Tight']:
            #    for btagregion in ['btagPass']:
            for lepcat in sorted(set([key[3] for key in histo.values().keys()])):
                for btagregion in sorted(set([key[2] for key in histo.values().keys()])):
                    pltdir = os.path.join(outdir, args.lepton, jmult, lepcat, btagregion)
                    if not os.path.isdir(pltdir):
                        os.makedirs(pltdir)

                    print(', '.join([jmult, lepcat, btagregion, hname])) 
                    hslice = histo[:, jmult, btagregion, lepcat].integrate('jmult').integrate('lepcat').integrate('btag')

                    if hname == 'Lep_iso':
                        if args.lepton == 'Muon':
                            x_lims = (0., 0.15) if lepcat == 'Tight' else (0.15, 1.)
                        if args.lepton == 'Electron':
                            x_lims = (0., 0.1) if lepcat == 'Tight' else (0., 0.5)

                        # plot original yields
                    fig, ax = plt.subplots()
                    fig.subplots_adjust(hspace=.07)

                    Plotter.plot_1D(hslice.values()[('Before',)], hslice.axis(xaxis_name).edges(), xlimits=x_lims, ax=ax, label='Runs $<$ 319077')
                    Plotter.plot_1D(hslice.values()[('After',)], hslice.axis(xaxis_name).edges(), xlimits=x_lims, xlabel=xtitle, color='r', ax=ax, label='Runs $\\geq$ 319077')
                    ax.legend(loc='upper right')

                    if hname == 'Jets_njets':
                        print(jmult)
                        #set_trace() 
                        rows = [("Lumi: %s fb^-1" % format(data_lumi_year['%ss' % args.lepton]/1000., '.1f'), "Yield", "Error", "Frac")]
                        rows += [("Runs < 319077", format(sum(hslice.values(overflow='all')[('Before',)]), '.1f'), format(np.sqrt(sum(hslice.values(overflow='all', sumw2=True)[('Before',)][1])), '.1f'),
                            format(sum(hslice.values(overflow='all')[('Before',)])/sum(hslice.sum('hem').values(overflow='all')[()]), '.3f'))]
                        rows += [("Runs >= 319077", format(sum(hslice.values(overflow='all')[('After',)]), '.1f'), format(np.sqrt(sum(hslice.values(overflow='all', sumw2=True)[('After',)][1])), '.1f'),
                            format(sum(hslice.values(overflow='all')[('After',)])/sum(hslice.sum('hem').values(overflow='all')[()]), '.3f'))]
                        rows += [("Total", format(sum(hslice.sum('hem').values(overflow='all')[()]), '.1f'), format(np.sqrt(sum(hslice.sum('hem').values(overflow='all', sumw2=True)[()][1])), '.1f'), "")]

                        frac_name = '%s_yields_and_fracs.txt' % '_'.join([jmult, args.lepton, lepcat, btagregion])
                        plt_tools.print_table(rows, filename=os.path.join(pltdir, frac_name), print_output=True)
                        print('%s/%s.txt written' % (pltdir, frac_name))

                        # add lepton/jet multiplicity label
                    ax.text(
                        0.02, 0.90, "%s, %s\n%s" % (lep_cats[lepcat], jet_mults[jmult], btag_cats[btagregion]),
                        fontsize=rcParams['font.size']*0.75, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
                    )
                    hep.cms.label(ax=ax, data=True, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))

                    #set_trace()
                    figname = os.path.join(pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, hname]))
                    fig.savefig(figname)
                    print('%s written' % figname)
                    plt.close()
    
                        # plot normalized yields
                    fig_norm, ax_norm = plt.subplots()
                    fig_norm.subplots_adjust(hspace=.07)

                    Plotter.plot_1D(hslice.values()[('Before',)]/sum(hslice.values()[('Before',)]), hslice.axis(xaxis_name).edges(), xlimits=x_lims, ylabel='A.U.', ax=ax_norm, label='Runs $<$ 319077')
                    Plotter.plot_1D(hslice.values()[('After',)]/sum(hslice.values()[('After',)]), hslice.axis(xaxis_name).edges(), xlimits=x_lims, xlabel=xtitle, ylabel='A.U.', color='r', ax=ax_norm, label='Runs $\\geq$ 319077')
                    ax_norm.legend(loc='upper right')

                        # add lepton/jet multiplicity label
                    ax_norm.text(
                        0.02, 0.90, "%s, %s\n%s" % (lep_cats[lepcat], jet_mults[jmult], btag_cats[btagregion]),
                        fontsize=rcParams['font.size']*0.75, horizontalalignment='left', verticalalignment='bottom', transform=ax_norm.transAxes
                    )
                    hep.cms.label(ax=ax_norm, data=True, paper=False, year=args.year)

                    #set_trace()
                    figname_norm = os.path.join(pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, hname, 'Norm']))
                    fig_norm.savefig(figname_norm)
                    print('%s written' % figname_norm)
                    plt.close()
    
    if histo.dense_dim() == 2:
        xtitle, ytitle, xrebinning, yrebinning, x_lims, y_lims, withData, plot_xlabels, plot_ylabels = variables[hname]

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
                    pltdir = os.path.join(outdir, args.lepton, jmult, lepcat, btagregion)
                    if not os.path.isdir(pltdir):
                        os.makedirs(pltdir)

                    if ((lepcat, btagregion) == ('Tight', 'btagPass')) and (hname == 'mtt_vs_tlep_ctstar_abs'):
                        continue
                    #else:
                    #    withData = variables[hname][-1]
                    print(', '.join([jmult, lepcat, btagregion, hname])) 

                    hslice = histo[:, jmult, btagregion, lepcat].integrate('jmult').integrate('lepcat').integrate('btag')

                    # plot linearized view 
                    hline_before = Plotter.linearize_hist(hslice['Before'].integrate('hem'))
                    hline_after = Plotter.linearize_hist(hslice['After'].integrate('hem'))
                    hline_edges = hline_before.axis(hline_before.dense_axes()[0].name).edges()

                        # plot original yields
                    fig, ax = plt.subplots()
                    fig.subplots_adjust(hspace=.07)

                    Plotter.plot_1D(hline_before.values()[()], hline_edges, xlimits=(min(hline_edges), max(hline_edges)), ax=ax, label='Runs $<$ 319077')
                    Plotter.plot_1D(hline_after.values()[()], hline_edges, xlimits=(min(hline_edges), max(hline_edges)), xlabel='%s $\otimes$ %s' % (xtitle, ytitle), color='r', ax=ax, label='Runs $\\geq$ 319077')
                    ax.legend(loc='upper right')

                        # draw vertical lines separating ctstar bins
                    bin_sep_lines = [hslice.values()[('Before',)].shape[0]*ybin for ybin in range(1, hslice.values()[('Before',)].shape[1])]
                    for binline in bin_sep_lines:
                        ax.axvline(binline, color='k', linestyle='--')

                        # plot unrolled x and y labels for each bin
                    ## plot x labels
                    if plot_xlabels is not None:
                        ax.set_xticks(np.arange(len(plot_xlabels)))
                        ax.set_xticklabels(plot_xlabels)
                        ax.tick_params(which='minor', bottom=False, top=False)
                        plt.setp(ax.get_xticklabels(), rotation=90, ha="left", fontsize=rcParams['font.size']*0.5)
                    
                    ## plot y labels
                    if plot_ylabels is not None: # (binlabels, bin_locs)
                        for idx, label in enumerate(plot_ylabels[0]):
                            ax.annotate(label, xy=(plot_ylabels[1][idx], 0), xycoords=('data', 'axes fraction'),
                                xytext=(0, -50 if plot_xlabels is None else -120), textcoords='offset points', va='top', ha='center', fontsize=rcParams['font.size']*0.5, rotation=45)

                        # add lepton/jet multiplicity label
                    ax.text(
                        0.02, 0.90, "%s, %s\n%s" % (lep_cats[lepcat], jet_mults[jmult], btag_cats[btagregion]),
                        fontsize=rcParams['font.size']*0.75, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
                    )
                    hep.cms.label(ax=ax, data=True, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))

                    #set_trace()
                    figname = os.path.join(pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, hname]))
                    fig.savefig(figname)
                    print('%s written' % figname)
                    plt.close()

                        # plot normalized yields
                    fig_norm, ax_norm = plt.subplots()
                    fig_norm.subplots_adjust(hspace=.07)

                    Plotter.plot_1D(hline_before.values()[()]/sum(hline_before.values()[()]), hline_edges, xlimits=(min(hline_edges), max(hline_edges)), ylabel='A.U.', ax=ax_norm, label='Runs $<$ 319077')
                    Plotter.plot_1D(hline_after.values()[()]/sum(hline_after.values()[()]), hline_edges, xlimits=(min(hline_edges), max(hline_edges)),
                        xlabel='%s $\otimes$ %s' % (xtitle, ytitle), ylabel='A.U.', color='r', ax=ax_norm, label='Runs $\\geq$ 319077')
                    ax_norm.legend(loc='upper right')

                        # draw vertical lines separating ctstar bins
                    bin_sep_lines = [hslice.values()[('Before',)].shape[0]*ybin for ybin in range(1, hslice.values()[('Before',)].shape[1])]
                    for binline in bin_sep_lines:
                        ax_norm.axvline(binline, color='k', linestyle='--')

                        # plot unrolled x and y labels for each bin
                    ## plot x labels
                    if plot_xlabels is not None:
                        ax_norm.set_xticks(np.arange(len(plot_xlabels)))
                        ax_norm.set_xticklabels(plot_xlabels)
                        ax_norm.tick_params(which='minor', bottom=False, top=False)
                        plt.setp(ax_norm.get_xticklabels(), rotation=90, ha="left", fontsize=rcParams['font.size']*0.5)
                    
                    ## plot y labels
                    if plot_ylabels is not None: # (binlabels, bin_locs)
                        for idx, label in enumerate(plot_ylabels[0]):
                            ax_norm.annotate(label, xy=(plot_ylabels[1][idx], 0), xycoords=('data', 'axes fraction'),
                                xytext=(0, -50 if plot_xlabels is None else -120), textcoords='offset points', va='top', ha='center', fontsize=rcParams['font.size']*0.5, rotation=45)

                        # add lepton/jet multiplicity label
                    ax_norm.text(
                        0.02, 0.90, "%s, %s\n%s" % (lep_cats[lepcat], jet_mults[jmult], btag_cats[btagregion]),
                        fontsize=rcParams['font.size']*0.75, horizontalalignment='left', verticalalignment='bottom', transform=ax_norm.transAxes
                    )
                    hep.cms.label(ax=ax_norm, data=True, paper=False, year=args.year)

                    figname_norm = os.path.join(pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, hname, 'Norm']))
                    fig_norm.savefig(figname_norm)
                    print('%s written' % figname_norm)
                    plt.close()

