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
parser.add_argument('--nosys', action='store_true', help='Make plots without systematics and no qcd estimation')
parser.add_argument('--qcd_est', action='store_true', help='Estimate qcd contribution')
parser.add_argument('--hem', choices=['not', 'scaled', 'removed'], default='not', help='Make plots when HEM correction is applied or not (default is not applied).')
args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
base_jobid = os.environ['base_jobid']
analyzer = 'hem_invest'

#set_trace()
if args.hem == 'not':
    hem_odir = 'Not_Applied'
    fname_str = 'NotApplied'
if args.hem == 'scaled':
    hem_odir = 'Scaled'
    fname_str = 'SCALED'
if args.hem == 'removed':
    hem_odir = 'Removed'
    fname_str = 'REMOVE'
    

input_dir = os.path.join(proj_dir, 'results', '%s_%s' % (args.year, jobid), analyzer)
f_ext = 'TOT.coffea'
outdir = os.path.join(proj_dir, 'plots', '%s_%s' % (args.year, jobid), analyzer, hem_odir)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

fnames = sorted(['%s/%s' % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)])
fname_to_use = [fname for fname in fnames if fname_str in fname]
hdict = plt_tools.add_coffea_files(fname_to_use) if len(fname_to_use) > 1 else load(fname_to_use[0])

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
    #np.array([-3.2, -2.4, -1.6, -0.85, 0., 0.85, 1.6, 2.4, 3.2]), # phi binning
    np.array([-3.2, -2.4, -1.57, -0.87, 0., 0.87, 1.57, 2.4, 3.2]), # phi
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
    'Jets_njets' : ('$n_{jets}$', 1, (0, 15), True),
    'mtt' : ('m($t\\bar{t}$) [GeV]', 4, (200., 2000.), True),
    'tlep_ctstar_abs' : ('|cos($\\theta^{*}_{t_{l}}$)|', 1, (0., 1.), True),
    'mthad' : ('m($t_{h}$) [GeV]', 2, (0., 300.), True),
    'pt_thad' : ('$p_{T}$($t_{h}$) [GeV]', 2, (0., 500.), True),
    'pt_tlep' : ('$p_{T}$($t_{l}$) [GeV]', 2, (0., 500.), True),
    'pt_tt' : ('$p_{T}$($t\\bar{t}$) [GeV]', 2, (0., 500.), True),
    'eta_thad' : ('$\\eta$($t_{h}$)', 2, (-4., 4.), True),
    'eta_tlep' : ('$\\eta$($t_{l}$)', 2, (-4., 4.), True),
    'eta_tt' : ('$\\eta$($t\\bar{t}$)', 2, (-4., 4.), True),
    'tlep_ctstar' : ('cos($\\theta^{*}_{t_{l}}$)', 2, (-1., 1.), True),
    'full_disc' : ('$\\lambda_{C}$', 2, (5, 25.), True),
    'mass_disc' : ('$\\lambda_{M}$', 2, (0, 20.), True),
    'ns_disc' : ('$\\lambda_{NS}$', 2, (3., 10.), True),
    'ns_dist' : ('$D_{\\nu, min}$', 1, (0., 150.), True),
    'Jets_pt' : ('$p_{T}$(jets) [GeV]', 2, (0., 300.), True),
    'Jets_eta' : ('$\\eta$(jets)', 2, (-2.6, 2.6), True),
    'Jets_phi' : ('$\\phi$(jets)', 2, (-4., 4.), True),
    'Jets_LeadJet_pt' : ('$p_{T}$(leading jet) [GeV]', 2, (0., 300.), True),
    'Jets_LeadJet_eta' : ('$\\eta$(leading jet)', 2, (-2.6, 2.6), True),
    'Lep_pt' : ('$p_{T}$(%s) [GeV]' % objtypes['Lep'][args.lepton], 2, (0., 300.), True),
    'Lep_eta' : ('$\\eta$(%s)' % objtypes['Lep'][args.lepton], 2, (-2.6, 2.6), True),
    'Lep_phi' : ('$\\phi$(%s)' % objtypes['Lep'][args.lepton], 2, (-4, 4), True),
    'Lep_iso' : ('pfRelIso, %s' % objtypes['Lep'][args.lepton], 1, (0., 1.), True),
    'MT' : ('$M_{T}$ [GeV]', 1, (0., 300.), True),
    'MET_pt' : ('$p_{T}$(MET) [GeV]', 1, (0., 300.), True),
    'MET_phi' : ('$\phi$(MET)', 1, (-3.2, 3.2), True),
}


    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, 'inputs', '%s_lumis_data.json' % base_jobid)).read())[args.year]
lumi_correction = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'MC_LumiWeights.coffea'))[args.year]['%ss' % args.lepton]
        # scale ttJets events, split by reconstruction type, by normal ttJets lumi correction
ttJets_permcats = ['*right', '*matchable', '*unmatchable', '*sl_tau', '*other']
names = [dataset for dataset in sorted(set([key[0] for key in hdict[sorted(variables.keys())[0]].values().keys()]))] # get dataset names in hists
ttJets_cats = [name for name in names if any([fnmatch.fnmatch(name, cat) for cat in ttJets_permcats])] # gets ttJets(_PS)_other, ...
if len(ttJets_cats) > 0:
    for tt_cat in ttJets_cats:
        ttJets_lumi_topo = '_'.join(tt_cat.split('_')[:-2]) if 'sl_tau' in tt_cat else '_'.join(tt_cat.split('_')[:-1]) # gets ttJets[SL, Had, DiLep] or ttJets_PS
        ttJets_eff_lumi = lumi_correction[ttJets_lumi_topo]
        lumi_correction.update({tt_cat: ttJets_eff_lumi})

## make groups based on process
process = hist.Cat("process", "Process", sorting='placement')
process_cat = "dataset"
process_groups = plt_tools.make_dataset_groups(args.lepton, args.year, samples=names)

    # scale and group hists by process
for hname in hdict.keys():
    if 'cutflow' in hname: continue
    #set_trace()
    hdict[hname].scale(lumi_correction, axis='dataset') # scale hists
    hdict[hname] = hdict[hname].group(process_cat, process, process_groups) # group by process
    hdict[hname] = hdict[hname][:, :, :, args.lepton].integrate('leptype') # only pick out specified lepton



    ## make nosys plots
if args.nosys:
    for hname in variables.keys():
        if hname not in hdict.keys():
            raise ValueError("%s not found in file" % hname)
        histo = hdict[hname][:, 'nosys', :, :, :].integrate('sys') # process, sys, jmult, btag, lepcat
        #set_trace()
    
        if histo.dense_dim() == 1:
            xtitle, rebinning, x_lims, withData = variables[hname]
            if rebinning != 1:
                xaxis_name = histo.dense_axes()[0].name
                histo = histo.rebin(xaxis_name, rebinning)
    
            ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
            for jmult in sorted(set([key[1] for key in histo.values().keys()])):
                for lepcat in sorted(set([key[3] for key in histo.values().keys()])):
                #for lepcat in ['Tight']:
                    #for btagregion in ['btagPass']:
                    for btagregion in sorted(set([key[2] for key in histo.values().keys()])):
                        pltdir = os.path.join(outdir, args.lepton, jmult, lepcat, btagregion)
                        if not os.path.isdir(pltdir):
                            os.makedirs(pltdir)
    
                        print(', '.join([jmult, lepcat, btagregion, hname])) 
                        fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                        fig.subplots_adjust(hspace=.07)
    
                        hslice = histo[:, jmult, btagregion, lepcat].integrate('jmult').integrate('lepcat').integrate('btag')
    
                        if hname == 'Lep_iso':
                            if args.lepton == 'Muon':
                                x_lims = (0., 0.15) if lepcat == 'Tight' else (0.15, 1.)
                            if args.lepton == 'Electron':
                                x_lims = (0., 0.1) if lepcat == 'Tight' else (0., 0.5)
    
                        mc_opts = {
                            #'mcorder' : ['QCD', 'EWK', 'singlet', 'ttJets'] if not ttJets_cats else ['QCD', 'EWK', 'singlet', 'ttJets_other', 'ttJets_unmatchable', 'ttJets_matchable', 'ttJets_right']
                            #'maskData' : not withData
                        }
    
                        Plotter.plot_stack1d(ax, rax, hslice, xlabel=xtitle, xlimits=x_lims, **mc_opts)
    
                        #set_trace() 
                        if hname == 'Jets_njets':
                            print(jmult)
                            yields_txt, yields_json = Plotter.get_samples_yield_and_frac(hslice, data_lumi_year['%ss' % args.lepton]/1000., promptmc=True)
                            frac_name = '%s_yields_and_fracs' % '_'.join([jmult, args.lepton, lepcat, btagregion])
                            plt_tools.print_table(yields_txt, filename='%s/%s.txt' % (pltdir, frac_name), print_output=True)
                            print('%s/%s.txt written' % (pltdir, frac_name))
                            with open('%s/%s.json' % (pltdir, frac_name), 'w') as out:
                                out.write(prettyjson.dumps(yields_json))
    
                            # add lepton/jet multiplicity label
                        #set_trace()
                        ax.text(
                            0.02, 0.88, "%s, %s\n%s" % (lep_cats[lepcat], jet_mults[jmult], btag_cats[btagregion]),
                            fontsize=rcParams['font.size']*0.75, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
                        )
                        hep.cms.label(ax=ax, data=withData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
    
                        #set_trace()
                        figname = os.path.join(pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, hname]))
                        fig.savefig(figname)
                        print('%s written' % figname)
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
                                figname = os.path.join(pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, hname, yaxis_name, 'proj']))
                            else:
                                xlabel = '%s [GeV]' % xtitle
                                xlimits = x_lims
                                figname = os.path.join(pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, hname, xaxis_name, 'proj']))
    
                            fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                            fig.subplots_adjust(hspace=.07)
    
                            hproj = hslice.integrate(hslice.dense_axes()[dax].name)
    
                            Plotter.plot_stack1d(ax, rax, hproj, xlabel=xlabel, xlimits=xlimits)
    
                                # add lepton/jet multiplicity label
                            ax.text(
                                0.02, 0.88, "%s, %s\n%s" % (lep_cats[lepcat], jet_mults[jmult], btag_cats[btagregion]),
                                fontsize=rcParams['font.size']*0.75, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
                            )
                            hep.cms.label(ax=ax, data=withData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
    
                            #set_trace()
                            fig.savefig(figname)
                            print('%s written' % figname)
                            plt.close()
    
                        #    ## make 1D plots of mtt for each ctstar bin
                        #for ybin in range(len(hslice.dense_axes()[1].edges())-1):
                        #    fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                        #    fig.subplots_adjust(hspace=.07)
    
                        #    hproj = hslice[:, :, hslice.dense_axes()[1].edges()[ybin]:hslice.dense_axes()[1].edges()[ybin+1]].integrate(hslice.dense_axes()[1].name)
    
                        #    Plotter.plot_stack1d(ax, rax, hproj, xlabel='%s [GeV]' % xtitle, xlimits=x_lims, **mc_opts)
    
                        #        # add lepton/jet multiplicity label, add ctstar range label
                        #    binlabel = '%s $\leq$ %s < %s' % (hslice.dense_axes()[1].edges()[ybin], ytitle, hslice.dense_axes()[1].edges()[ybin+1])
                        #    ax.text(
                        #        0.02, 0.88, "%s, %s\t%s\n%s" % (lep_cats[lepcat], jet_mults[jmult], binlabel, btag_cats[btagregion]),
                        #        fontsize=rcParams['font.size']*0.75, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
                        #    )
                        #    hep.cms.label(ax=ax, data=withData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
    
                        #    #set_trace()
                        #    bintitle = '%sctstar%s' % (hslice.dense_axes()[1].edges()[ybin], hslice.dense_axes()[1].edges()[ybin+1])
                        #    bintitle = bintitle.replace('.', 'p')
                        #    figname = os.path.join(pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, bintitle, 'mtt']))
                        #    fig.savefig(figname)
                        #    print('%s written' % figname)
                        #    plt.close()
    
    
                            # plot linearized view 
                        fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                        fig.subplots_adjust(hspace=.07)
    
                        hline = Plotter.linearize_hist(hslice)
                        
                        Plotter.plot_stack1d(ax, rax, hline, xlabel='%s $\otimes$ %s' % (xtitle, ytitle), xlimits=(0, len(hline.axis(hline.dense_axes()[0].name).edges())-1), **mc_opts)
    
                            # draw vertical lines separating ctstar bins
                        bin_sep_lines = [hslice.values()[('EWK',)].shape[0]*ybin for ybin in range(1, hslice.values()[('EWK',)].shape[1])]
                        for binline in bin_sep_lines:
                            ax.axvline(binline, color='k', linestyle='--')
                            if rax is not None: rax.axvline(binline, color='k', linestyle='--')
    
                            # plot unrolled x and y labels for each bin
                        ## plot x labels
                        if plot_xlabels is not None:
                            rax.set_xticks(np.arange(len(plot_xlabels)))
                            rax.set_xticklabels(plot_xlabels)
                            ax.tick_params(which='minor', bottom=False, top=False)
                            rax.tick_params(which='minor', bottom=False, top=False)
                            plt.setp(rax.get_xticklabels(), rotation=90, ha="left", fontsize=rcParams['font.size']*0.5)
                        
                        ## plot y labels
                        if plot_ylabels is not None: # (binlabels, bin_locs)
                            for idx, label in enumerate(plot_ylabels[0]):
                                rax.annotate(label, xy=(plot_ylabels[1][idx], 0), xycoords=('data', 'axes fraction'),
                                    xytext=(0, -50 if plot_xlabels is None else -120), textcoords='offset points', va='top', ha='center', fontsize=rcParams['font.size']*0.5, rotation=45)
    
                            # add lepton/jet multiplicity label
                        ax.text(
                            0.02, 0.88, "%s, %s\n%s" % (lep_cats[lepcat], jet_mults[jmult], btag_cats[btagregion]),
                            fontsize=rcParams['font.size']*0.75, 
                            horizontalalignment='left', 
                            verticalalignment='bottom', 
                            transform=ax.transAxes
                        )
                        hep.cms.label(ax=ax, data=withData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
    
                        #set_trace()
                        figname = os.path.join(pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, hname]))
                        fig.savefig(figname)
                        print('%s written' % figname)
                        plt.close()

# plot with QCD estimation
if args.qcd_est:
    import Utilities.systematics as systematics
    sys_to_name = systematics.sys_to_name[args.year]

    systs_to_run = ['nosys']
    #systs_to_run = sorted(sys_to_name.keys()) if args.save_sys else ['nosys']
    print('\nPlotting distributions for systematics:\n\t', *systs_to_run)

    #set_trace()
        ## save post qcd-est dists for all systematics
    for hname in variables.keys():
        if hname not in hdict.keys():
            raise ValueError("%s not found in file" % hname)
        #set_trace()
        histo = hdict[hname] # process, sys, jmult, btag, lepcat
        #histo = hdict[hname][:, :, :, :, :] # process, sys, jmult, btag, lepcat

        if histo.dense_dim() == 1:
            xtitle, xrebinning, x_lims, withData = variables[hname]
            xaxis_name = histo.dense_axes()[0].name
                ## rebin x axis
            if isinstance(xrebinning, np.ndarray):
                new_xbins = hist.Bin(xaxis_name, xaxis_name, xrebinning)
            elif isinstance(xrebinning, float) or isinstance(xrebinning, int):
                new_xbins = xrebinning
            histo = histo.rebin(xaxis_name, new_xbins)

        if histo.dense_dim() == 2:
            xtitle, ytitle, xrebinning, yrebinning, x_lims, y_lims, withData, plot_xlabels, plot_ylabels = variables[hname]
            #xtitle, ytitle, xrebinning, yrebinning, x_lims, y_lims, withData = variables[hname]

            #withData = False
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
        mc_opts = {
        #    'mcorder' : ['QCD', 'EWK', 'singlet', 'ttJets'] if not ttJets_cats else ['QCD', 'EWK', 'singlet', 'ttJets_other', 'ttJets_unmatchable', 'ttJets_matchable', 'ttJets_right']
            'maskData' : not withData,
        }

        vlines = [(len(xrebinning)-1)*ybin for ybin in range(1, len(yrebinning)-1)] if histo.dense_dim() == 2 else None
        nbins = (len(xrebinning)-1)*(len(yrebinning)-1) if histo.dense_dim() == 2 else None

        if hname == 'Lep_iso':
            x_lims = (0., 0.15) if args.lepton == 'Muon' else (0., 0.1)

        ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
        for jmult in sorted(set([key[2] for key in histo.values().keys()])):
        #for jmult in ['3Jets']:
        #for jmult in ['4PJets']:
            print(jmult)
            if histo.dense_dim() == 1:
                iso_sb    = histo[:, 'nosys', jmult, 'btagPass', 'Loose'].integrate('sys').integrate('jmult').integrate('lepcat').integrate('btag')
                btag_sb   = histo[:, 'nosys', jmult, 'btagFail', 'Tight'].integrate('sys').integrate('jmult').integrate('lepcat').integrate('btag')
                double_sb = histo[:, 'nosys', jmult, 'btagFail', 'Loose'].integrate('sys').integrate('jmult').integrate('lepcat').integrate('btag')
                sig_histo = histo[:, :, jmult, 'btagPass', 'Tight'].integrate('jmult').integrate('lepcat').integrate('btag')
            else:
                iso_sb_histo    = sideband_group[:, 'nosys', jmult, 'btagPass', 'Loose'].integrate('sys').integrate('jmult').integrate('lepcat').integrate('btag')
                btag_sb_histo   = sideband_group[:, 'nosys', jmult, 'btagFail', 'Tight'].integrate('sys').integrate('jmult').integrate('lepcat').integrate('btag')
                double_sb_histo = sideband_group[:, 'nosys', jmult, 'btagFail', 'Loose'].integrate('sys').integrate('jmult').integrate('lepcat').integrate('btag')
                iso_sb    = Plotter.linearize_hist(iso_sb_histo)
                btag_sb   = Plotter.linearize_hist(btag_sb_histo)
                double_sb = Plotter.linearize_hist(double_sb_histo)
                sig_histo = Plotter.linearize_hist(sig_group[:, :, jmult, 'btagPass', 'Tight'].integrate('jmult').integrate('lepcat').integrate('btag')) 

                ## plot systematic variations after qcd estimation is made
            #set_trace()
            for sys in systs_to_run:
                if sys not in histo.axis('sys')._sorted:
                    print('\n\n   Systematic %s not available, skipping\n\n' % sys)
                    continue
                print('QCD est:', jmult, sys, hname)
                #set_trace()

                shape_reg = 'BTAG'
                for norm in ['Sideband']:
                #for norm in ['ABCD', 'Sideband']:
                    qcd_est_histo = Plotter.QCD_Est(sig_reg=sig_histo, iso_sb=iso_sb, btag_sb=btag_sb, double_sb=double_sb, norm_type=norm, shape_region=shape_reg, norm_region=shape_reg if norm=='Sideband' else None, sys=sys)
                    if sys == 'nosys':
                        qcd_name = '%s%s_Norm' % (shape_reg, norm) if norm == 'Sideband' else '%s_Norm' % norm
                        qcd_dir = os.path.join(outdir, args.lepton, jmult, 'QCD_Est', sys_to_name[sys], qcd_name)
                        if not os.path.isdir(qcd_dir):
                            os.makedirs(qcd_dir)

                        fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                        fig.subplots_adjust(hspace=.07)

                        if hname == 'Jets_njets':
                            print('QCD est:', jmult, sys)
                            yields_txt, yields_json = Plotter.get_samples_yield_and_frac(qcd_est_histo, data_lumi_year['%ss' % args.lepton]/1000., sys=sys)
                            frac_name = '%s_yields_and_fracs_QCD_Est_%s' % ('_'.join([sys, jmult, args.lepton]), qcd_name)
                            plt_tools.print_table(yields_txt, filename='%s/%s.txt' % (qcd_dir, frac_name), print_output=True)
                            print('%s/%s.txt written' % (qcd_dir, frac_name))
                            with open('%s/%s.json' % (qcd_dir, frac_name), 'w') as out:
                                out.write(prettyjson.dumps(yields_json))

                        #set_trace()
                        Plotter.plot_stack1d(ax, rax, qcd_est_histo, xlabel=xtitle, xlimits=x_lims, **mc_opts) if histo.dense_dim() == 1\
                            else Plotter.plot_stack1d(ax, rax, qcd_est_histo, xlabel='%s $\otimes$ %s' % (xtitle, ytitle), xlimits=(0, nbins), **mc_opts)

                            # add lepton/jet multiplicity label
                        ax.text(
                            0.02, 0.92, "%s, %s" % (objtypes['Lep'][args.lepton], jet_mults[jmult]),
                            fontsize=rcParams['font.size']*0.9, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
                        )
                                ## draw vertical lines for distinguishing different ctstar bins
                        if vlines is not None:
                            for vline in vlines:
                                ax.axvline(vline, color='k', linestyle='--')
                                if rax is not None: rax.axvline(vline, color='k', linestyle='--')
                        hep.cms.label(ax=ax, data=withData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))

                        #set_trace()
                        figname = '%s/%s_QCD_Est_%s' % (qcd_dir, '_'.join([sys, jmult, args.lepton, hname]), qcd_name)
                        fig.savefig(figname)
                        print('%s written' % figname)
                        plt.close()
