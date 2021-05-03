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
import Utilities.systematics as systematics
from copy import deepcopy

base_jobid = os.environ['base_jobid']

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('year', choices=['2016APV', '2016', '2017', '2018'] if base_jobid == 'ULnanoAOD' else ['2016', '2017', '2018'], help='What year is the ntuple from.')
parser.add_argument('lepton', choices=['Electron', 'Muon'], help='Choose which lepton to make plots for')
parser.add_argument('--nosys', action='store_true', help='Make plots without systematics and no qcd estimation')
parser.add_argument('--qcd_est', action='store_true', help='Estimate qcd contribution')
parser.add_argument('--plot_shapes', action='store_true', help='Make plots of qcd shape, and shape and normalization uncertainties')
parser.add_argument('--save_sys', action='store_true', help='Save dists for all systematic variations')
args = parser.parse_args()

sys_to_name = systematics.sys_to_name[args.year]

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'htt_btag_iso_cut'

input_dir = os.path.join(proj_dir, 'results', '%s_%s' % (args.year, jobid), analyzer)
f_ext = 'TOT.coffea'
outdir = os.path.join(proj_dir, 'plots', '%s_%s' % (args.year, jobid), analyzer)
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

ctstar_binlabels = ['%s $\leq$ |cos($\\theta^{*}_{t_{l}}$)| $\leq$ %s' % (linearize_binning[1][bin], linearize_binning[1][bin+1]) for bin in range(len(linearize_binning[1])-1)]
ctstar_bin_locs = np.linspace((len(linearize_binning[0])-1)/2, (len(linearize_binning[0])-1)*(len(linearize_binning[1])-1) - (len(linearize_binning[0])-1)/2, len(linearize_binning[1])-1)
mtt_binlabels = ['%s $\leq$ m($t\\bar{t}$) $\leq$ %s' % (linearize_binning[0][bin], linearize_binning[0][bin+1]) for bin in range(len(linearize_binning[0])-1)]*len(ctstar_binlabels)

# binning and bin labels for phi x eta
phi_eta_binning = (
    np.array([-3.2, -2.4, -1.57, -0.87, 0., 0.87, 1.57, 2.4, 3.2]), # phi
    np.array([-2.5, -1.3, -0.7, 0., 0.7, 1.3, 2.5]) # eta binning
)
eta_binlabels = ['%s $\leq$ $\\eta$ $\leq$ %s' % (phi_eta_binning[1][bin], phi_eta_binning[1][bin+1]) for bin in range(len(phi_eta_binning[1])-1)]
eta_bin_locs = np.linspace((len(phi_eta_binning[0])-1)/2, (len(phi_eta_binning[0])-1)*(len(phi_eta_binning[1])-1) - (len(phi_eta_binning[0])-1)/2, len(phi_eta_binning[1])-1)
phi_binlabels = ['%s $\leq$ $\\phi$ $\leq$ %s' % (phi_eta_binning[0][bin], phi_eta_binning[0][bin+1]) for bin in range(len(phi_eta_binning[0])-1)]*len(eta_binlabels)


variables = {
    'mtt_vs_tlep_ctstar_abs' : ('m($t\\bar{t}$)', '|cos($\\theta^{*}_{t_{l}}$)|', linearize_binning[0], linearize_binning[1], (linearize_binning[0][0], linearize_binning[0][-1]),
        (linearize_binning[1][0], linearize_binning[1][-1]), True, None, (ctstar_binlabels, ctstar_bin_locs)),
    'Jets_njets' : ('$n_{jets}$', 1, (0, 15), True),
    'Jets_phi_vs_eta' : ('$\\phi$(jets)', '$\\eta$(jets)', phi_eta_binning[0], phi_eta_binning[1], (-3.3, 3.3), (-2.6, 2.6), True, phi_binlabels, (eta_binlabels, eta_bin_locs)),
    'Lep_phi_vs_eta' : ('$\\phi$(%s)' % objtypes['Lep'][args.lepton], '$\\eta$(%s)' % objtypes['Lep'][args.lepton], phi_eta_binning[0], phi_eta_binning[1],
        (-3.3, 3.3), (-2.6, 2.6), True, phi_binlabels, (eta_binlabels, eta_bin_locs)),
    ##'mtt' : ('m($t\\bar{t}$) [GeV]', linearize_binning[0], (200., 2000.), True),
    ##'tlep_ctstar_abs' : ('|cos($\\theta^{*}_{t_{l}}$)|', linearize_binning[1], (0., 1.), True),
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
    #hdict[hname] = hdict[hname].group(process_cat, process, process_groups) # group by process
    hdict[hname] = hdict[hname][:, :, :, args.lepton].integrate('leptype') # only pick out specified lepton



def plot_shape_uncs(numerators=[], denom={}, **opts):
    """
    Inputs:
        numerators is a list of dictionaries containing the histogram and style options
        denom is a dictionary containing the histogram and style options
        opts is a dictionary of plotting options for the entire plot
    """
    if not denom: return

    fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
    fig.subplots_adjust(hspace=.07)

        # get opts
    xlimits = opts.get('xlims')
    xlabel = opts.get('xtitle')
    jmult = opts.get('jmult')
    pltdir = opts.get('pltdir')
    hname = opts.get('hname')
    vlines = opts.get('vlines')

    labels = []
    dmp_denom_hist = denom['histo'].integrate('process').copy()

    # get normalized arrays of data-promptMC
        # denominator
    dmp_denom_sumw, dmp_denom_sumw2 = Plotter.get_qcd_shape(Plotter.data_minus_prompt(denom['histo']))
    for idx in range(len(dmp_denom_sumw)):
        dmp_denom_hist.values(overflow='all', sumw2=True)[()][0][idx] = dmp_denom_sumw[idx]
        dmp_denom_hist.values(overflow='all', sumw2=True)[()][1][idx] = dmp_denom_sumw2[idx]

        # numerators
    for reg in numerators:
        dmp_num_hist = dmp_denom_hist.copy()
        dmp_num_sumw, dmp_num_sumw2 = Plotter.get_qcd_shape(Plotter.data_minus_prompt(reg['histo']))
        for idx in range(len(dmp_num_sumw)):
            dmp_num_hist.values(overflow='all', sumw2=True)[()][0][idx] = dmp_num_sumw[idx]
            dmp_num_hist.values(overflow='all', sumw2=True)[()][1][idx] = dmp_num_sumw2[idx]

            # plot step hist for bin ranges
        hist.plot1d(dmp_num_hist, ax=ax, clear=False,
            line_opts={'color' : reg['color']},
        )
            # plot error bars
        hist.plot1d(dmp_num_hist, ax=ax, clear=False,
            error_opts={'color' : reg['color']},
        )
        labels.append(reg['label'])

            # plot num/denom ratio
        hist.plotratio(
            dmp_num_hist, dmp_denom_hist, error_opts={'marker': '.', 'markersize': 10., 'color': reg['color'], 'elinewidth': 1},
            unc='num', clear=False, ax=rax, guide_opts={}, label='/'.join([reg['label'].split(' ')[0], denom['label'].split(' ')[0]])
        )

    # plot denom
        # plot step hist for bin ranges
    hist.plot1d(dmp_denom_hist, ax=ax, clear=False,
        line_opts={'color' : denom['color']},
    )
        # plot error bars
    hist.plot1d(dmp_denom_hist, ax=ax, clear=False,
        error_opts={'color' : denom['color']},
    )
    labels.append(denom['label'])

        ## set legend and corresponding colors
    handles, def_labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc='upper right', title='data-$MC_{prompt}$')

    ax.autoscale()
    ax.set_ylabel('Probability Density')
    ax.set_ylim(0, 1 if ax.get_ylim()[1] > 1 else None)
    ax.set_xlabel(None)
    ax.set_xlim(xlimits)

    #rax.legend(loc='lower right')
    rax.legend(loc='upper right')
    rax.autoscale()
    rax.set_ylabel('Ratio')
    #set_trace()
    rax.set_ylim(0., 5. if rax.get_ylim()[1] > 5. else None)
    #rax.set_ylim(0., 5.)
    rax.set_xlim(xlimits)
    rax.set_xlabel(xlabel)

    if vlines is not None:
        #set_trace()
            # draw vertical lines separating ctstar bins
        for vline in vlines:
            ax.axvline(vline, color='k', linestyle='--')
            rax.axvline(vline, color='k', linestyle='--')

    ax.text(
        0.02, 0.92, "%s, %s" % (objtypes['Lep'][args.lepton], jet_mults[jmult]),
        fontsize=rcParams['font.size'], horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
    )
    hep.cms.label(ax=ax, data=withData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
    #set_trace()

    if not os.path.isdir(pltdir):
        os.makedirs(pltdir)
    pltname = opts.get('fname')
    figname = os.path.join(pltdir, pltname)
    fig.savefig(figname)
    print('%s written' % figname)
    #set_trace()


def plot_qcd_dmp_shape(btag_sb=None, iso_sb=None, double_sb=None,**opts):
    fig, ax = plt.subplots()
    fig.subplots_adjust(hspace=.07)

        # get opts
    xlimits = opts.get('xlims')
    xlabel = opts.get('xtitle')
    jmult = opts.get('jmult')
    pltdir = opts.get('pltdir')
    hname = opts.get('hname')
    vlines = opts.get('vlines')

    labels = []
    if btag_sb is not None:
        btag_dmp_hist = btag_sb.integrate('process').copy()

            # get normalized arrays of data-promptMC    
        btag_dmp_sumw, btag_dmp_sumw2 = Plotter.get_qcd_shape(Plotter.data_minus_prompt(btag_sb))

        for idx in range(len(btag_dmp_sumw)):
            btag_dmp_hist.values(overflow='all', sumw2=True)[()][0][idx] = btag_dmp_sumw[idx]
            btag_dmp_hist.values(overflow='all', sumw2=True)[()][1][idx] = btag_dmp_sumw2[idx]

        hist.plot1d(btag_dmp_hist, ax=ax, clear=False,
            line_opts={'color' : 'k'},
        )
        hist.plot1d(btag_dmp_hist, ax=ax, clear=False,
            error_opts={'color' : 'k'},
        )
        labels.append('B (BTAG)')
    if iso_sb is not None:
        iso_dmp_hist = iso_sb.integrate('process').copy()

            # get normalized arrays of data-promptMC    
        iso_dmp_sumw, iso_dmp_sumw2 = Plotter.get_qcd_shape(Plotter.data_minus_prompt(iso_sb))

        for idx in range(len(iso_dmp_sumw)):
            iso_dmp_hist.values(overflow='all', sumw2=True)[()][0][idx] = iso_dmp_sumw[idx]
            iso_dmp_hist.values(overflow='all', sumw2=True)[()][1][idx] = iso_dmp_sumw2[idx]

        hist.plot1d(iso_dmp_hist, ax=ax, clear=False,
            line_opts={'color' : 'r'},
        )
        hist.plot1d(iso_dmp_hist, ax=ax, clear=False,
            error_opts={'color' : 'r'},
        )
        labels.append('C (ISO)')
    if double_sb is not None:
        double_dmp_hist = double_sb.integrate('process').copy()

            # get normalized arrays of data-promptMC    
        double_dmp_sumw, double_dmp_sumw2 = Plotter.get_qcd_shape(Plotter.data_minus_prompt(double_sb))

        for idx in range(len(double_dmp_sumw)):
            double_dmp_hist.values(overflow='all', sumw2=True)[()][0][idx] = double_dmp_sumw[idx]
            double_dmp_hist.values(overflow='all', sumw2=True)[()][1][idx] = double_dmp_sumw2[idx]

        hist.plot1d(double_dmp_hist, ax=ax, clear=False,
            line_opts={'color' : 'b'},
        )
        hist.plot1d(double_dmp_hist, ax=ax, clear=False,
            error_opts={'color' : 'b'},
        )
        labels.append('D (BTAG-ISO)')
        
        ## set legend and corresponding colors
    handles, def_labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc='upper right', title='data-$MC_{prompt}$')

    ax.autoscale()
    ax.set_ylabel('Probability Density')
    ax.set_ylim(0, None)
    ax.set_xlabel(xlabel)
    ax.set_xlim(xlimits)

    if vlines is not None:
        #set_trace()
            # draw vertical lines separating ctstar bins
        for vline in vlines:
            ax.axvline(vline, color='k', linestyle='--')

    ax.text(
        0.02, 0.94, "%s, %s" % (objtypes['Lep'][args.lepton], jet_mults[jmult]),
        fontsize=rcParams['font.size'],
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax.transAxes
    )
    hep.cms.label(ax=ax, data=withData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))

    #set_trace()
    if not os.path.isdir(pltdir):
        os.makedirs(pltdir)
    pltname = opts.get('fname')
    figname = os.path.join(pltdir, pltname)
    fig.savefig(figname)
    print('%s written' % figname)


def plot_qcd_mc_shape(sig_reg=None, btag_sb=None, iso_sb=None, double_sb=None, print_ratios=False, **opts):
    fig, ax = plt.subplots()
    fig.subplots_adjust(hspace=.07)

        # get opts
    xlimits = opts.get('xlims')
    xlabel = opts.get('xtitle')
    jmult = opts.get('jmult')
    pltdir = opts.get('pltdir')
    hname = opts.get('hname')
    vlines = opts.get('vlines')

    labels = []
    yields_and_ratios_dict = {}

    if sig_reg is not None:
        sig_qcd_hist = sig_reg[('QCD',)].integrate('process').copy()
        mc_sumw, mc_sumw2 = np.sum(sig_qcd_hist.values(overflow='all', sumw2=True)[()][0]), np.sum(sig_qcd_hist.values(overflow='all', sumw2=True)[()][1])
        yields_and_ratios_dict['A'] = [mc_sumw, mc_sumw2]

            # get normalized arrays of QCD MC    
        sig_mc_sumw, sig_mc_sumw2 = Plotter.get_qcd_shape(sig_qcd_hist)

        for idx in range(len(sig_mc_sumw)):
            sig_qcd_hist.values(overflow='all', sumw2=True)[()][0][idx] = sig_mc_sumw[idx]
            sig_qcd_hist.values(overflow='all', sumw2=True)[()][1][idx] = sig_mc_sumw2[idx]

        hist.plot1d(sig_qcd_hist, ax=ax, clear=False,
            line_opts={'color' : 'g'},
        )
        #hist.plot1d(sig_qcd_hist, ax=ax, clear=False,
        #    error_opts={'color' : 'g'},
        #)
        labels.append('A (Signal)')

    if btag_sb is not None:
        btag_qcd_hist = btag_sb[('QCD',)].integrate('process').copy()
        mc_sumw, mc_sumw2 = np.sum(btag_qcd_hist.values(overflow='all', sumw2=True)[()][0]), np.sum(btag_qcd_hist.values(overflow='all', sumw2=True)[()][1])
        yields_and_ratios_dict['B'] = [mc_sumw, mc_sumw2]

            # get normalized arrays of QCD MC    
        btag_mc_sumw, btag_mc_sumw2 = Plotter.get_qcd_shape(btag_qcd_hist)

        for idx in range(len(btag_mc_sumw)):
            btag_qcd_hist.values(overflow='all', sumw2=True)[()][0][idx] = btag_mc_sumw[idx]
            btag_qcd_hist.values(overflow='all', sumw2=True)[()][1][idx] = btag_mc_sumw2[idx]

        hist.plot1d(btag_qcd_hist, ax=ax, clear=False,
            line_opts={'color' : 'k'},
        )
        #hist.plot1d(btag_qcd_hist, ax=ax, clear=False,
        #    error_opts={'color' : 'k'},
        #)
        labels.append('B (BTAG)')

    if iso_sb is not None:
        iso_qcd_hist = iso_sb[('QCD',)].integrate('process').copy()
        mc_sumw, mc_sumw2 = np.sum(iso_qcd_hist.values(overflow='all', sumw2=True)[()][0]), np.sum(iso_qcd_hist.values(overflow='all', sumw2=True)[()][1])
        yields_and_ratios_dict['C'] = [mc_sumw, mc_sumw2]

            # get normalized arrays of QCD MC    
        iso_mc_sumw, iso_mc_sumw2 = Plotter.get_qcd_shape(iso_qcd_hist)

        for idx in range(len(iso_mc_sumw)):
            iso_qcd_hist.values(overflow='all', sumw2=True)[()][0][idx] = iso_mc_sumw[idx]
            iso_qcd_hist.values(overflow='all', sumw2=True)[()][1][idx] = iso_mc_sumw2[idx]

        hist.plot1d(iso_qcd_hist, ax=ax, clear=False,
            line_opts={'color' : 'r'},
        )
        #hist.plot1d(iso_qcd_hist, ax=ax, clear=False,
        #    error_opts={'color' : 'r'},
        #)
        labels.append('C (ISO)')

    if double_sb is not None:
        double_qcd_hist = double_sb[('QCD',)].integrate('process').copy()
        mc_sumw, mc_sumw2 = np.sum(double_qcd_hist.values(overflow='all', sumw2=True)[()][0]), np.sum(double_qcd_hist.values(overflow='all', sumw2=True)[()][1])
        yields_and_ratios_dict['D'] = [mc_sumw, mc_sumw2]

            # get normalized arrays of QCD MC    
        double_mc_sumw, double_mc_sumw2 = Plotter.get_qcd_shape(double_qcd_hist)

        for idx in range(len(double_mc_sumw)):
            double_qcd_hist.values(overflow='all', sumw2=True)[()][0][idx] = double_mc_sumw[idx]
            double_qcd_hist.values(overflow='all', sumw2=True)[()][1][idx] = double_mc_sumw2[idx]

        hist.plot1d(double_qcd_hist, ax=ax, clear=False,
            line_opts={'color' : 'b'},
        )
        #hist.plot1d(double_qcd_hist, ax=ax, clear=False,
        #    error_opts={'color' : 'b'},
        #)
        labels.append('D (BTAG-ISO)')
        
        ## set legend and corresponding colors
    handles, def_labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc='upper right', title='QCD Simulation')

    ax.autoscale()
    ax.set_ylabel('Probability Density')
    ax.set_ylim(0, None)
    ax.set_xlabel(xlabel)
    ax.set_xlim(xlimits)

    if vlines is not None:
        #set_trace()
            # draw vertical lines separating ctstar bins
        for vline in vlines:
            ax.axvline(vline, color='k', linestyle='--')

    ax.text(
        0.02, 0.94, "%s, %s" % (objtypes['Lep'][args.lepton], jet_mults[jmult]),
        fontsize=rcParams['font.size'],
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax.transAxes
    )
    hep.cms.label(ax=ax, data=withData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))

    # save yields, errors, and ratios of QCD MC for different regions
    if print_ratios:
            # get ratio of A/C
        A_over_C_val, A_over_C_err = Plotter.get_ratio_and_uncertainty(a_val=yields_and_ratios_dict['A'][0], a_err=np.sqrt(yields_and_ratios_dict['A'][1]),
            b_val=yields_and_ratios_dict['C'][0], b_err=np.sqrt(yields_and_ratios_dict['C'][1]))
            # get ratio of B/D
        B_over_D_val, B_over_D_err = Plotter.get_ratio_and_uncertainty(a_val=yields_and_ratios_dict['B'][0], a_err=np.sqrt(yields_and_ratios_dict['B'][1]),
            b_val=yields_and_ratios_dict['D'][0], b_err=np.sqrt(yields_and_ratios_dict['D'][1]))
            # get ratio of A/B
        A_over_B_val, A_over_B_err = Plotter.get_ratio_and_uncertainty(a_val=yields_and_ratios_dict['A'][0], a_err=np.sqrt(yields_and_ratios_dict['A'][1]),
            b_val=yields_and_ratios_dict['B'][0], b_err=np.sqrt(yields_and_ratios_dict['B'][1]))
            # get ratio of C/D
        C_over_D_val, C_over_D_err = Plotter.get_ratio_and_uncertainty(a_val=yields_and_ratios_dict['C'][0], a_err=np.sqrt(yields_and_ratios_dict['C'][1]),
            b_val=yields_and_ratios_dict['D'][0], b_err=np.sqrt(yields_and_ratios_dict['D'][1]))

        rows = [("Lumi: %s fb^-1" % format(round(data_lumi_year['%ss' % args.lepton]/1000., 1), '.1f'), "QCD MC Yield", "Error")]
        rows += [("A", format(yields_and_ratios_dict['A'][0], '.3f'), format(np.sqrt(yields_and_ratios_dict['A'][1]), '.3f'))]
        rows += [("B", format(yields_and_ratios_dict['B'][0], '.3f'), format(np.sqrt(yields_and_ratios_dict['B'][1]), '.3f'))]
        rows += [("C", format(yields_and_ratios_dict['C'][0], '.3f'), format(np.sqrt(yields_and_ratios_dict['C'][1]), '.3f'))]
        rows += [("D", format(yields_and_ratios_dict['D'][0], '.3f'), format(np.sqrt(yields_and_ratios_dict['D'][1]), '.3f'))]
        rows += [("", "", "")]
        rows += [("A/C", format(A_over_C_val, '.3f'), format(A_over_C_err, '.3f'))]
        rows += [("B/D", format(B_over_D_val, '.3f'), format(B_over_D_err, '.3f'))]
        rows += [("A/B", format(A_over_B_val, '.3f'), format(A_over_B_err, '.3f'))]
        rows += [("C/D", format(C_over_D_val, '.3f'), format(C_over_D_err, '.3f'))]
        ratios_txt_name = os.path.join(pltdir, 'QCD_MC_Yield_Ratios_%s_%s.txt' % (args.lepton, jmult))
        plt_tools.print_table(rows, filename=ratios_txt_name, print_output=True)
        print('%s written' % ratios_txt_name)

    if not os.path.isdir(pltdir):
        os.makedirs(pltdir)
    pltname = opts.get('fname')
    figname = os.path.join(pltdir, pltname)
    fig.savefig(figname)
    print('%s written' % figname)


def get_norm_unc(reg_A, reg_B):
    ## relative norm unc = N_dmp,B/N_qcd,B - 1
    N_dmp_B = np.sum(Plotter.data_minus_prompt(reg_B).values(overflow='all')[()])
    N_qcd_B = np.sum(reg_B['QCD*'].integrate('process').values(overflow='all')[()])
    rel_norm_unc = N_dmp_B/N_qcd_B - 1.

    ## relative stat. unc = sum over bins(stat error of QCD,A)/N_qcd,A
    sumw_qcd_A, sumw2_qcd_A = reg_A['QCD*'].integrate('process').values(overflow='all', sumw2=True)[()]
    rel_stat_unc = np.sum(sumw2_qcd_A)/np.sum(sumw_qcd_A)

    rows = [("Relative Norm. Unc:", format(rel_norm_unc, '.3f'))]
    rows += [("Relative Stat. Unc:", format(rel_stat_unc, '.3f'))]
    return rows


    ## make plots
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
            #for jmult in ['3Jets']:
            #for jmult in ['4PJets']:
            #    for lepcat in ['Tight']:
            #        for btagregion in ['btagPass']:
            for jmult in sorted(set([key[1] for key in histo.values().keys()])):
                for lepcat in sorted(set([key[3] for key in histo.values().keys()])):
                    for btagregion in sorted(set([key[2] for key in histo.values().keys()])):
                        pltdir = os.path.join(outdir, args.lepton, jmult, lepcat, btagregion)
                        if not os.path.isdir(pltdir):
                            os.makedirs(pltdir)
   
                        print(', '.join([jmult, lepcat, btagregion, hname])) 
                        fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                        fig.subplots_adjust(hspace=.07)
    
                        #set_trace()
                        tmp_groups = deepcopy(process_groups)
                        if (jobid == 'ULnanoAOD_noMT_noLJpt') and (args.year == '2018') and (args.lepton == 'Electron') and (jmult == '3Jets') and (lepcat == 'Tight') and (btagregion == 'btagPass'): tmp_groups['QCD'].remove('QCD_EM_30to50')
                        if (jobid == 'ULnanoAOD_noMT_noLJpt') and (args.year == '2017') and (args.lepton == 'Muon') and (jmult == '4PJets') and (lepcat == 'Tight') and (btagregion == 'btagPass'): tmp_groups['QCD'].remove('QCD_Mu_80to120')
                        if (jobid == 'ULnanoAOD_noMT_noLJpt') and (args.year == '2017') and (args.lepton == 'Electron') and (jmult == '3Jets') and (lepcat == 'Tight') and (btagregion == 'btagPass'):
                            tmp_groups['QCD'].remove('QCD_Mu_80to120')
                            tmp_groups['QCD'].remove('QCD_EM_50to80')

                        h_group = histo.group(process_cat, process, tmp_groups) # group by process
                        hslice = h_group[:, jmult, btagregion, lepcat].integrate('jmult').integrate('lepcat').integrate('btag')
   
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
            #for jmult in ['3Jets']:
            ##for jmult in ['4PJets']:
            #    for lepcat in ['Tight']:
            #        for btagregion in ['btagPass']:
            for jmult in sorted(set([key[1] for key in histo.values().keys()])):
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

                        ## TEMPORARY
                        tmp_groups = deepcopy(process_groups)
                        if (jobid == 'ULnanoAOD_noMT_noLJpt') and (args.year == '2018') and (args.lepton == 'Electron') and (jmult == '3Jets') and (lepcat == 'Tight') and (btagregion == 'btagPass'): tmp_groups['QCD'].remove('QCD_EM_30to50')
                        if (jobid == 'ULnanoAOD_noMT_noLJpt') and (args.year == '2017') and (args.lepton == 'Muon') and (jmult == '4PJets') and (lepcat == 'Tight') and (btagregion == 'btagPass'): tmp_groups['QCD'].remove('QCD_Mu_80to120')
                        if (jobid == 'ULnanoAOD_noMT_noLJpt') and (args.year == '2017') and (args.lepton == 'Electron') and (jmult == '3Jets') and (lepcat == 'Tight') and (btagregion == 'btagPass'):
                            tmp_groups['QCD'].remove('QCD_Mu_80to120')
                            tmp_groups['QCD'].remove('QCD_EM_50to80')

                        h_group = histo.group(process_cat, process, tmp_groups) # group by process
                        hslice = h_group[:, jmult, btagregion, lepcat].integrate('jmult').integrate('lepcat').integrate('btag')
   
                        mc_opts = {
                            'maskData' : not withData
                        #    'mcorder' : ['QCD', 'EWK', 'singlet', 'ttJets'] if not ttJets_cats else ['QCD', 'EWK', 'singlet', 'ttJets_other', 'ttJets_unmatchable', 'ttJets_matchable', 'ttJets_right']
                        }

                            # make 1D projection along dense axes
                        for dax in range(2):
                            if dax == 0:
                                xlabel = ytitle
                                xlimits = y_lims
                                figname = os.path.join(pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, yaxis_name, 'proj']))
                            else:
                                xlabel = '%s [GeV]' % xtitle
                                xlimits = x_lims
                                figname = os.path.join(pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, xaxis_name, 'proj']))

                            fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                            fig.subplots_adjust(hspace=.07)

                            hproj = hslice.integrate(hslice.dense_axes()[dax].name)
    
                            Plotter.plot_stack1d(ax, rax, hproj, xlabel=xlabel, xlimits=xlimits, **mc_opts)
    
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
    systs_to_run = sorted(sys_to_name.keys()) if args.save_sys else ['nosys']
    print('\nPlotting distributions for systematics:\n\t', *systs_to_run)

    #set_trace()
        ## save post qcd-est dists for all systematics
    if (args.save_sys):
        import itertools
            # keys are (lepton, jmult, sys, hname)
        save_keys = sorted(itertools.product([args.lepton], ['3Jets', '4PJets'], [sys_to_name[sys] for sys in systs_to_run], [name for name in variables.keys()]))
        post_qcd_est = {key:None for key in save_keys}

    for hname in variables.keys():
        if hname not in hdict.keys():
            raise ValueError("%s not found in file" % hname)
        #set_trace()
        histo = hdict[hname][:, :, :, :, :] # process, sys, jmult, btag, lepcat

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

            if (hname == 'mtt_vs_tlep_ctstar_abs'):
                withData = False
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
            sig_groups = deepcopy(process_groups)
                # remove single event with really high weight
            if (jobid == 'ULnanoAOD_noMT_noLJpt') and (args.year == '2018') and (args.lepton == 'Electron') and (jmult == '3Jets'): sig_groups['QCD'].remove('QCD_EM_30to50')
            if (jobid == 'ULnanoAOD_noMT_noLJpt') and (args.year == '2017') and (args.lepton == 'Muon') and (jmult == '4PJets'): sig_groups['QCD'].remove('QCD_Mu_80to120')
            if (jobid == 'ULnanoAOD_noMT_noLJpt') and (args.year == '2017') and (args.lepton == 'Electron') and (jmult == '3Jets'):
                sig_groups['QCD'].remove('QCD_Mu_80to120')
                sig_groups['QCD'].remove('QCD_EM_50to80')

            sig_group = histo.group(process_cat, process, sig_groups) # group by process
            sideband_group = histo.group(process_cat, process, process_groups) # group by process

            if histo.dense_dim() == 1:
                iso_sb    = sideband_group[:, 'nosys', jmult, 'btagPass', 'Loose'].integrate('sys').integrate('jmult').integrate('lepcat').integrate('btag')
                btag_sb   = sideband_group[:, 'nosys', jmult, 'btagFail', 'Tight'].integrate('sys').integrate('jmult').integrate('lepcat').integrate('btag')
                double_sb = sideband_group[:, 'nosys', jmult, 'btagFail', 'Loose'].integrate('sys').integrate('jmult').integrate('lepcat').integrate('btag')
                sig_histo = sig_group[:, :, jmult, 'btagPass', 'Tight'].integrate('jmult').integrate('lepcat').integrate('btag')
            else:
                iso_sb_histo    = sideband_group[:, 'nosys', jmult, 'btagPass', 'Loose'].integrate('sys').integrate('jmult').integrate('lepcat').integrate('btag')
                btag_sb_histo   = sideband_group[:, 'nosys', jmult, 'btagFail', 'Tight'].integrate('sys').integrate('jmult').integrate('lepcat').integrate('btag')
                double_sb_histo = sideband_group[:, 'nosys', jmult, 'btagFail', 'Loose'].integrate('sys').integrate('jmult').integrate('lepcat').integrate('btag')
                nosys_sig_reg_histo = sig_group[:, 'nosys', jmult, 'btagPass', 'Tight'].integrate('sys').integrate('jmult').integrate('lepcat').integrate('btag')
                iso_sb    = Plotter.linearize_hist(iso_sb_histo)
                btag_sb   = Plotter.linearize_hist(btag_sb_histo)
                double_sb = Plotter.linearize_hist(double_sb_histo)
                sig_histo = Plotter.linearize_hist(sig_group[:, :, jmult, 'btagPass', 'Tight'].integrate('jmult').integrate('lepcat').integrate('btag'))

                ## plot shape uncertainties and find normalization uncs
            if args.plot_shapes:
                #set_trace()
                if (hname == 'Jets_njets') or (hname == 'mtt_vs_tlep_ctstar_abs'):
                    unc_dir = os.path.join(outdir, args.lepton, jmult, 'QCD_Est', 'Uncertainties')
                    if not os.path.isdir(unc_dir): os.makedirs(unc_dir)
                    shape_dir = os.path.join(outdir, args.lepton, jmult, 'QCD_Est', 'Shapes')
                    if not os.path.isdir(shape_dir): os.makedirs(shape_dir)

                    if hname == 'Jets_njets':
                        sig = sig_group[:, 'nosys', jmult, 'btagPass', 'Tight'].integrate('sys').integrate('jmult').integrate('lepcat').integrate('btag')
                        norm_uncs = get_norm_unc(reg_A=sig, reg_B=btag_sb)
                        unc_txt_name = 'Norm_Uncs_%s_%s' % (args.lepton, jmult)
                        #set_trace()
                        plt_tools.print_table(norm_uncs, filename='%s/%s.txt' % (unc_dir, unc_txt_name), print_output=True)
                        print('%s/%s.txt written' % (os.path.join(outdir, args.lepton, jmult, 'QCD_Est'), unc_txt_name))

                    if hname == 'mtt_vs_tlep_ctstar_abs':
                        ## plot qcd shapes (MC and data-prompt)
                            ## QCD MC
                        #set_trace()
                                ## plot projection onto mtt
                        plot_qcd_mc_shape(sig_reg=nosys_sig_reg_histo.integrate('ctstar_abs'), btag_sb=btag_sb_histo.integrate('ctstar_abs'),
                            iso_sb=iso_sb_histo.integrate('ctstar_abs'), double_sb=double_sb_histo.integrate('ctstar_abs'), print_ratios=True,
                            **{'xtitle':'%s [GeV]' % xtitle, 'xlims':x_lims, 'pltdir':shape_dir, 'hname': 'mtt', 'jmult':jmult, 'fname': 'QCD_Simulation_Shapes_%s_%s_mtt' % (args.lepton, jmult)})
                                ## plot projection onto ctstar
                        plot_qcd_mc_shape(sig_reg=nosys_sig_reg_histo.integrate('mtt'), btag_sb=btag_sb_histo.integrate('mtt'),
                            iso_sb=iso_sb_histo.integrate('mtt'), double_sb=double_sb_histo.integrate('mtt'),
                            **{'xtitle':ytitle, 'xlims':y_lims, 'pltdir':shape_dir, 'hname': 'tlep_ctstar_abs', 'jmult':jmult, 'fname': 'QCD_Simulation_Shapes_%s_%s_tlep_ctstar_abs' % (args.lepton, jmult)})
                                ## plot linearized 2D mtt ctstar dist
                        plot_qcd_mc_shape(sig_reg=sig_histo[:, 'nosys', :].integrate('sys'), btag_sb=btag_sb, iso_sb=iso_sb, double_sb=double_sb,
                            **{'xtitle':'%s $\otimes$ %s' % (xtitle, ytitle), 'xlims':(0, nbins), 'pltdir':shape_dir, 'hname': hname, 'jmult':jmult,
                            'vlines':vlines, 'fname': 'QCD_Simulation_Shapes_%s_%s_%s' % (args.lepton, jmult, hname)}
                        )
                            ## data-prompt
                                ## plot projection onto mtt
                        plot_qcd_dmp_shape(btag_sb=btag_sb_histo.integrate('ctstar_abs'), iso_sb=iso_sb_histo.integrate('ctstar_abs'), double_sb=double_sb_histo.integrate('ctstar_abs'),
                            **{'xtitle':'%s [GeV]' % xtitle, 'xlims':x_lims, 'pltdir':shape_dir, 'hname': 'mtt', 'jmult':jmult, 'fname': 'QCD_DMP_Shapes_%s_%s_mtt' % (args.lepton, jmult)})
                                ## plot projection onto ctstar
                        plot_qcd_dmp_shape(btag_sb=btag_sb_histo.integrate('mtt'), iso_sb=iso_sb_histo.integrate('mtt'), double_sb=double_sb_histo.integrate('mtt'),
                            **{'xtitle':ytitle, 'xlims':y_lims, 'pltdir':shape_dir, 'hname': 'tlep_ctstar_abs', 'jmult':jmult, 'fname': 'QCD_DMP_Shapes_%s_%s_tlep_ctstar_abs' % (args.lepton, jmult)})
                                ## plot linearized 2D mtt ctstar dist
                        plot_qcd_dmp_shape(btag_sb=btag_sb, iso_sb=iso_sb, double_sb=double_sb,
                            **{'xtitle':'%s $\otimes$ %s' % (xtitle, ytitle), 'xlims':(0, nbins), 'pltdir':shape_dir, 'hname': hname, 'jmult':jmult, 'vlines':vlines, 'fname': 'QCD_DMP_Shapes_%s_%s_%s' % (args.lepton, jmult, hname)})

                        ## plot shape uncertainties
                            ## plot projection onto mtt
                        plot_shape_uncs(numerators=[{'histo' : btag_sb_histo.integrate('ctstar_abs'), 'label' : 'B (BTAG)', 'color' : 'k'}, {'histo' : iso_sb_histo.integrate('ctstar_abs'), 'label' : 'C (ISO)', 'color' : 'r'}],
                            denom={'histo' : double_sb_histo.integrate('ctstar_abs'), 'label' : 'D (BTAG-ISO)', 'color' : 'b'}, 
                            **{'xtitle':'%s [GeV]' % xtitle, 'xlims':x_lims, 'pltdir':unc_dir, 'hname': 'mtt', 'jmult':jmult, 'fname':'Shape_Uncertainties_%s_%s_mtt' % (args.lepton, jmult)})
                            ## plot projection onto ctstar
                        plot_shape_uncs(numerators=[{'histo' : btag_sb_histo.integrate('mtt'), 'label' : 'B (BTAG)', 'color' : 'k'}, {'histo' : iso_sb_histo.integrate('mtt'), 'label' : 'C (ISO)', 'color' : 'r'}],
                            denom={'histo' : double_sb_histo.integrate('mtt'), 'label' : 'D (BTAG-ISO)', 'color' : 'b'}, 
                            **{'xtitle':ytitle, 'xlims':y_lims, 'pltdir':unc_dir, 'hname': 'tlep_ctstar_abs', 'jmult':jmult, 'fname':'Shape_Uncertainties_%s_%s_tlep_ctstar_abs' % (args.lepton, jmult)})
                            ## plot linearized 2D mtt ctstar dist
                        plot_shape_uncs(numerators=[{'histo' : btag_sb, 'label' : 'B (BTAG)', 'color' : 'k'}, {'histo' : iso_sb, 'label' : 'C (ISO)', 'color' : 'r'}],
                            denom={'histo' : double_sb, 'label' : 'D (BTAG-ISO)', 'color' : 'b'}, 
                            **{'xtitle':'%s $\otimes$ %s' % (xtitle, ytitle), 'xlims':(0, nbins), 'pltdir':unc_dir, 'hname': hname, 'jmult':jmult, 'vlines':vlines, 'fname':'Shape_Uncertainties_%s_%s_%s' % (args.lepton, jmult, hname)})


                ## plot systematic variations after qcd estimation is made
            #set_trace()
            for sys in systs_to_run:
                if sys not in histo.axis('sys')._sorted:
                    print('\n\n   Systematic %s not available, skipping\n\n' % sys)
                    continue
                print('QCD est:', jmult, sys, hname)
                #set_trace()

                shape_reg = 'BTAG'
                #for norm in ['Sideband']:
                #for norm in ['ABCD-ML']:
                for norm in ['ABCD']:
                #for norm in ['ABCD', 'Sideband']:
                    #set_trace()
                    qcd_est_histo = Plotter.QCD_Est(sig_reg=sig_histo, iso_sb=iso_sb, btag_sb=btag_sb, double_sb=double_sb, norm_type=norm, shape_region=shape_reg, norm_region=shape_reg if norm=='Sideband' else None, sys=sys)
                    #set_trace()
                        # save post qcd-est dist to dict
                    if norm == 'Sideband':
                        if args.save_sys:
                            post_qcd_est[(args.lepton, jmult, sys_to_name[sys], hname)] = qcd_est_histo

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

                            # plot unrolled x and y labels for each bin
                        if histo.dense_dim() == 2:
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

                        hep.cms.label(ax=ax, data=withData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
    
                        #set_trace()
                        figname = '%s/%s_QCD_Est_%s' % (qcd_dir, '_'.join([sys, jmult, args.lepton, hname]), qcd_name)
                        fig.savefig(figname)
                        print('%s written' % figname)
                        plt.close()


    if args.save_sys:
        from coffea.util import save
        outname = '%s/%s/%s_Post_QCD_Est_dict.coffea' % (outdir, args.lepton, args.lepton)
        save(post_qcd_est, outname)
        #outname = '%s/%s/%s_QCD_Est_mtt_ctstar_dict.coffea' % (outdir, args.lepton, args.lepton)
        #save(mtt_ctstar_qcd_est, outname)
        print('%s written' % outname)

