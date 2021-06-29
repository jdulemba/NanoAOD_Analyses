# matplotlib
import matplotlib
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
#parser.add_argument('--dmp', action='store_true', help='Make plots of data minus prompt MC shapes')
parser.add_argument('--bkg_shapes', action='store_true', help='Make plots comparing EKW+QCD MC and data-driven shapes.')
parser.add_argument('--bkg_est', action='store_true', help='Estimate ewk+qcd contribution')
parser.add_argument('--vary_norm', action='store_true', help='Make plots varying normalization uncertainty for ewk+qcd estimation.')
parser.add_argument('--vary_shape', action='store_true', help='Make plots varying shape for ewk+qcd estimation.')
parser.add_argument('--vary_shape_and_norm', action='store_true', help='Make plots varying shape and norm for ewk+qcd estimation.')
#parser.add_argument('--save_sys', action='store_true', help='Save dists for all systematic variations')
args = parser.parse_args()

sys_to_name = systematics.sys_to_name[args.year]

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'htt_btag_sb_regions'

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


#btag_cats = {
#    '0p05' : '0.00 < max(btag discr) $\leq$ 0.05',
#    '0p10' : '0.05 < max(btag discr) $\leq$ 0.10',
#    '0p15' : '0.10 < max(btag discr) $\leq$ 0.15',
#    '0p20' : '0.15 < max(btag discr) $\leq$ 0.20',
#    '0p25' : '0.20 < max(btag discr) $\leq$ 0.25',
#    '0p30' : '0.25 < max(btag discr) $\leq$ 0.30',
#    '0p35' : '0.30 < max(btag discr) $\leq$ 0.35',
#    '0p40' : '0.35 < max(btag discr) $\leq$ 0.40',
#    'btagPass' : '$n_{btags} \geq$ 2',
#    '0p0to0p2' : '0.0 < max(btag discr) $\leq$ 0.2',
#    '0p1to0p4' : '0.1 < max(btag discr) $\leq$ 0.4',
#    '0p2to0p5' : '0.2 < max(btag discr) $\leq$ 0.5',
#    '0p0to0p15' : '0.0 < max(btag discr) $\leq$ 0.15',
#    '0p15to0p3' : '0.15 < max(btag discr) $\leq$ 0.3',
#    '0p3to0p45' : '0.3 < max(btag discr) $\leq$ 0.45',
#}
btag_cats = {
    'p00p15' : '0.00 < max(btag discr) $\leq$ 0.15',
    'p15p30' : '0.15 < max(btag discr) $\leq$ 0.30',
    'p30p45' : '0.30 < max(btag discr) $\leq$ 0.45',
    'btagPass' : '$n_{btags} \geq$ 2',
}

btag_reg_names_dict = {
    'Signal' : {'reg' : 'btagPass'},
    'Central': {'reg' : 'p15p30', 'label' : 'Cen (0.15-0.3)', 'color' : 'k'},
    'Up'     : {'reg' : 'p30p45', 'label' : 'Up (0.3-0.45)', 'color' : 'r'},
    'Down'   : {'reg' : 'p00p15', 'label' : 'Down (0.0-0.15)', 'color' : 'b'}
    #'Central': {'reg' : '0p15to0p3', 'label' : 'Cen (0.15-0.3)', 'color' : 'k'},
    #'Up'     : {'reg' : '0p3to0p45', 'label' : 'Up (0.3-0.45)', 'color' : 'r'},
    #'Down'   : {'reg' : '0p0to0p15', 'label' : 'Down (0.0-0.15)', 'color' : 'b'}
}

lep_cats = {
    'Tight' : 'tight %s' % objtypes['Lep'][args.lepton],
    'Loose' : 'loose %s' % objtypes['Lep'][args.lepton],
}

linearize_binning = (
    #np.array([300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0]),
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
    'Jets_phi_vs_eta' : ('$\\phi$(jets)', '$\\eta$(jets)', phi_eta_binning[0], phi_eta_binning[1], (-3.3, 3.3), (-2.6, 2.6), True, phi_binlabels, (eta_binlabels, eta_bin_locs)),
    'Lep_phi_vs_eta' : ('$\\phi$(%s)' % objtypes['Lep'][args.lepton], '$\\eta$(%s)' % objtypes['Lep'][args.lepton], phi_eta_binning[0], phi_eta_binning[1],
        (-3.3, 3.3), (-2.6, 2.6), True, phi_binlabels, (eta_binlabels, eta_bin_locs)),
    ##'mtt' : ('m($t\\bar{t}$) [GeV]', linearize_binning[0], (200., 2000.), True),
    ##'tlep_ctstar_abs' : ('|cos($\\theta^{*}_{t_{l}}$)|', linearize_binning[1], (0., 1.), True),
    'Jets_njets' : ('$n_{jets}$', 1, (0, 15), True),
    'mtt' : ('m($t\\bar{t}$) [GeV]', 4, (200., 2000.), True),
    'tlep_ctstar_abs' : ('|cos($\\theta^{*}_{t_{l}}$)|', 1, (0., 1.), True),
    'mthad' : ('m($t_{h}$) [GeV]', 2, (0., 300.), True),
    'mWHad' : ('m($W_{h}$) [GeV]', 2, (0., 300.), True),
    'mWLep' : ('m($W_{l}$) [GeV]', 2, (0., 300.), True),
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
    'Jets_pt' : ('$p_{T}$(jets) [GeV]', 1, (0., 300.), True),
    'Jets_eta' : ('$\\eta$(jets)', 2, (-2.6, 2.6), True),
    'Jets_phi' : ('$\\phi$(jets)', 2, (-4., 4.), True),
    'Jets_LeadJet_pt' : ('$p_{T}$(leading jet) [GeV]', 1, (0., 300.), True),
    'Jets_LeadJet_eta' : ('$\\eta$(leading jet)', 2, (-2.6, 2.6), True),
    'Lep_pt' : ('$p_{T}$(%s) [GeV]' % objtypes['Lep'][args.lepton], 1, (0., 300.), True),
    'Lep_eta' : ('$\\eta$(%s)' % objtypes['Lep'][args.lepton], 2, (-2.6, 2.6), True),
    'Lep_phi' : ('$\\phi$(%s)' % objtypes['Lep'][args.lepton], 2, (-4, 4), True),
    'Lep_iso' : ('pfRelIso, %s' % objtypes['Lep'][args.lepton], 1, (0., 1.), True),
    'MT' : ('$M_{T}$ [GeV]', 1, (0., 300.), True),
    'MET_pt' : ('$p_{T}$(MET) [GeV]', 1, (0., 300.), True),
    'MET_phi' : ('$\phi$(MET)', 1, (-3.2, 3.2), True),
    #'DeepCSV_bDisc' : ('DeepCSV b Disc', 1, (-0.01, 1.), True),
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

## make btagging sideband groups
btag_bin = hist.Cat("btag", "btag", sorting='placement')
btag_cat = "btag"
btag_groups = {
    'btagPass' : ['btagPass'],
    ##'0p0to0p2' : ['0p05', '0p10', '0p15', '0p20'],
    ##'0p1to0p4' : ['0p15', '0p20', '0p25', '0p30', '0p35', '0p40'],
    ##'0p2to0p5' : ['0p25', '0p30', '0p35', '0p40', '0p45', '0p50'],
    #'0p0to0p15' : ['0p05', '0p10', '0p15'],
    #'0p15to0p3' : ['0p20', '0p25', '0p30'],
    #'0p3to0p45' : ['0p35', '0p40', '0p45'],
    'p00p15' : ['p00p15'],
    'p15p30' : ['p15p30'],
    'p30p45' : ['p30p45'],
}

    # scale and group hists by process
for hname in hdict.keys():
    if 'cutflow' in hname: continue
    #set_trace()
    hdict[hname].scale(lumi_correction, axis='dataset') # scale hists
    hdict[hname] = hdict[hname].group(btag_cat, btag_bin, btag_groups) # group by btagging sideband regions
    hdict[hname] = hdict[hname].group(process_cat, process, process_groups) # group by process
    hdict[hname] = hdict[hname][:, :, :, :, args.lepton].integrate('leptype') # only pick out specified lepton



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


def plot_qcd_dmp_shape(cen_sb={}, up_sb={}, dw_sb={}, **opts):
    fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
    fig.subplots_adjust(hspace=.07)

        # get opts
    xlimits = opts.get('xlims')
    xlabel = opts.get('xtitle')
    jmult = opts.get('jmult')
    pltdir = opts.get('pltdir')
    vlines = opts.get('vlines')

    labels = []
    if cen_sb:
        cen_dmp_hist = cen_sb['histo'].integrate('process').copy()

            # get normalized arrays of data-promptMC    
        cen_dmp_sumw, cen_dmp_sumw2 = Plotter.get_qcd_shape(Plotter.data_minus_prompt(cen_sb['histo']))

        for idx in range(len(cen_dmp_sumw)):
            cen_dmp_hist.values(overflow='all', sumw2=True)[()][0][idx] = cen_dmp_sumw[idx]
            cen_dmp_hist.values(overflow='all', sumw2=True)[()][1][idx] = cen_dmp_sumw2[idx]

        hist.plot1d(cen_dmp_hist, ax=ax, clear=False,
            line_opts={'color' : cen_sb['color']},
        )
        hist.plot1d(cen_dmp_hist, ax=ax, clear=False,
            error_opts={'color' : cen_sb['color']},
        )
        labels.append(cen_sb['label'])

    if up_sb:
        up_dmp_hist = up_sb['histo'].integrate('process').copy()

            # get normalized arrays of data-promptMC    
        up_dmp_sumw, up_dmp_sumw2 = Plotter.get_qcd_shape(Plotter.data_minus_prompt(up_sb['histo']))

        for idx in range(len(up_dmp_sumw)):
            up_dmp_hist.values(overflow='all', sumw2=True)[()][0][idx] = up_dmp_sumw[idx]
            up_dmp_hist.values(overflow='all', sumw2=True)[()][1][idx] = up_dmp_sumw2[idx]

        hist.plot1d(up_dmp_hist, ax=ax, clear=False,
            line_opts={'color' : up_sb['color']},
        )
        hist.plot1d(up_dmp_hist, ax=ax, clear=False,
            error_opts={'color' : up_sb['color']},
        )
        labels.append(up_sb['label'])

            # plot num/denom ratio
        hist.plotratio(
            up_dmp_hist, cen_dmp_hist, error_opts={'marker': '.', 'markersize': 10., 'color': up_sb['color'], 'elinewidth': 1},
            unc='num', clear=False, ax=rax, guide_opts={}#, label='/'.join([reg['label'].split(' ')[0], denom['label'].split(' ')[0]])
        )

    if dw_sb:
        dw_dmp_hist = dw_sb['histo'].integrate('process').copy()

            # get normalized arrays of data-promptMC    
        dw_dmp_sumw, dw_dmp_sumw2 = Plotter.get_qcd_shape(Plotter.data_minus_prompt(dw_sb['histo']))

        for idx in range(len(dw_dmp_sumw)):
            dw_dmp_hist.values(overflow='all', sumw2=True)[()][0][idx] = dw_dmp_sumw[idx]
            dw_dmp_hist.values(overflow='all', sumw2=True)[()][1][idx] = dw_dmp_sumw2[idx]

        hist.plot1d(dw_dmp_hist, ax=ax, clear=False,
            line_opts={'color' : dw_sb['color']},
        )
        hist.plot1d(dw_dmp_hist, ax=ax, clear=False,
            error_opts={'color' : dw_sb['color']},
        )
        labels.append(dw_sb['label'])
        
            # plot num/denom ratio
        hist.plotratio(
            dw_dmp_hist, cen_dmp_hist, error_opts={'marker': '.', 'markersize': 10., 'color': dw_sb['color'], 'elinewidth': 1},
            unc='num', clear=False, ax=rax, guide_opts={}#, label='/'.join([reg['label'].split(' ')[0], denom['label'].split(' ')[0]])
        )


    #set_trace()
        ## set legend and corresponding colors
    handles, def_labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc='upper right', title='data-$MC_{prompt}$')

    ax.autoscale()
    ax.set_ylabel('Probability Density')
    ax.set_ylim(0, None)
    ax.set_xlabel(None)
    ax.set_xlim(xlimits)

    rax.autoscale()
    rax.set_ylabel('Sys/Cen')
    rax.set_ylim(max(0.0, rax.get_ylim()[0]), min(1.5, rax.get_ylim()[1]))
    #rax.set_ylim(max(0.0, rax.get_ylim()[0]), min(5.0, rax.get_ylim()[1]))
    rax.set_xlim(xlimits)
    rax.set_xlabel(xlabel)

    if vlines is not None:
        #set_trace()
            # draw vertical lines separating ctstar bins
        for vline in vlines:
            ax.axvline(vline, color='k', linestyle='--')

    ax.text(
        0.02, 0.94, "%s, %s" % (objtypes['Lep'][args.lepton], jet_mults[jmult]),
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


def plot_ewk_qcd_dmtop_shape(cen_sb={}, up_sb={}, dw_sb={}, **opts):
    fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
    fig.subplots_adjust(hspace=.07)

        # get opts
    xlimits = opts.get('xlims')
    xlabel = opts.get('xtitle')
    jmult = opts.get('jmult')
    pltdir = opts.get('pltdir')
    vlines = opts.get('vlines')
    overflow = opts.get("overflow", "none")

    labels = []
    if cen_sb:
        cen_dmp_hist = cen_sb['histo'].integrate('process').copy()

            # get normalized arrays of data-promptMC    
        cen_dmp_sumw, cen_dmp_sumw2 = Plotter.get_qcd_shape(Plotter.data_minus_top(cen_sb['histo']))

        for idx in range(len(cen_dmp_sumw)):
            cen_dmp_hist.values(overflow='all', sumw2=True)[()][0][idx] = cen_dmp_sumw[idx]
            cen_dmp_hist.values(overflow='all', sumw2=True)[()][1][idx] = cen_dmp_sumw2[idx]

        hist.plot1d(cen_dmp_hist, ax=ax, clear=False,
            line_opts={'color' : cen_sb['color']},
            overflow=overflow,
            overlay_overflow=overflow,
        )
        hist.plot1d(cen_dmp_hist, ax=ax, clear=False,
            error_opts={'color' : cen_sb['color']},
            overflow=overflow,
            overlay_overflow=overflow,
        )
        labels.append(cen_sb['label'])

    if up_sb:
        up_dmp_hist = up_sb['histo'].integrate('process').copy()

            # get normalized arrays of data-promptMC    
        up_dmp_sumw, up_dmp_sumw2 = Plotter.get_qcd_shape(Plotter.data_minus_top(up_sb['histo']))

        for idx in range(len(up_dmp_sumw)):
            up_dmp_hist.values(overflow='all', sumw2=True)[()][0][idx] = up_dmp_sumw[idx]
            up_dmp_hist.values(overflow='all', sumw2=True)[()][1][idx] = up_dmp_sumw2[idx]

        hist.plot1d(up_dmp_hist, ax=ax, clear=False,
            line_opts={'color' : up_sb['color']},
            overflow=overflow,
            overlay_overflow=overflow,
        )
        hist.plot1d(up_dmp_hist, ax=ax, clear=False,
            error_opts={'color' : up_sb['color']},
            overflow=overflow,
            overlay_overflow=overflow,
        )
        labels.append(up_sb['label'])

            # plot num/denom ratio
        hist.plotratio(
            up_dmp_hist, cen_dmp_hist, error_opts={'marker': '.', 'markersize': 10., 'color': up_sb['color'], 'elinewidth': 1},
            unc='num', clear=False, ax=rax, guide_opts={},
            overflow=overflow,
            #overlay_overflow=overflow,
        )

    if dw_sb:
        dw_dmp_hist = dw_sb['histo'].integrate('process').copy()

            # get normalized arrays of data-promptMC    
        dw_dmp_sumw, dw_dmp_sumw2 = Plotter.get_qcd_shape(Plotter.data_minus_top(dw_sb['histo']))

        for idx in range(len(dw_dmp_sumw)):
            dw_dmp_hist.values(overflow='all', sumw2=True)[()][0][idx] = dw_dmp_sumw[idx]
            dw_dmp_hist.values(overflow='all', sumw2=True)[()][1][idx] = dw_dmp_sumw2[idx]

        hist.plot1d(dw_dmp_hist, ax=ax, clear=False,
            line_opts={'color' : dw_sb['color']},
            overflow=overflow,
            overlay_overflow=overflow,
        )
        hist.plot1d(dw_dmp_hist, ax=ax, clear=False,
            error_opts={'color' : dw_sb['color']},
            overflow=overflow,
            overlay_overflow=overflow,
        )
        labels.append(dw_sb['label'])
        
            # plot num/denom ratio
        hist.plotratio(
            dw_dmp_hist, cen_dmp_hist, error_opts={'marker': '.', 'markersize': 10., 'color': dw_sb['color'], 'elinewidth': 1},
            unc='num', clear=False, ax=rax, guide_opts={},
            overflow=overflow,
            #overlay_overflow=overflow,
        )


    #set_trace()
        ## set legend and corresponding colors
    handles, def_labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc='upper right', title='data-$MC_{st/t\\bar{t}}$')

    ax.autoscale()
    ax.set_ylabel('Probability Density')
    ax.set_ylim(0, None)
    ax.set_xlabel(None)
    ax.set_xlim(xlimits)

    rax.autoscale()
    rax.set_ylabel('Sys/Cen')
    rax.set_ylim(max(0.0, rax.get_ylim()[0]), min(1.5, rax.get_ylim()[1]))
    #rax.set_ylim(max(0.0, rax.get_ylim()[0]), min(5.0, rax.get_ylim()[1]))
    rax.set_xlim(xlimits)
    rax.set_xlabel(xlabel)

    if vlines is not None:
        #set_trace()
            # draw vertical lines separating ctstar bins
        for vline in vlines:
            ax.axvline(vline, color='k', linestyle='--')
            rax.axvline(vline, color='k', linestyle='--')

    ax.text(
        0.02, 0.94, "%s, %s" % (objtypes['Lep'][args.lepton], jet_mults[jmult]),
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


def plot_ewk_qcd_cont(signal={}, cen_sb={}, up_sb={}, dw_sb={}, **opts):
    if not signal: raise ValueError("You need a distribution from the signal region")

    fig, ax = plt.subplots()
    fig.subplots_adjust(hspace=.07)

        # get opts
    xlimits = opts.get('xlims')
    xlabel = opts.get('xtitle')
    jmult = opts.get('jmult')
    pltdir = opts.get('pltdir')
    vlines = opts.get('vlines')
    cols_to_label = {}

    #set_trace()
    if signal:
            # get EWK+QCD MC
        sig_bkg_hist = signal['histo'][Plotter.bkg_samples]
        sig_norm = sig_bkg_hist.sum(*[ax.name for ax in sig_bkg_hist.axes()]).values()[()]
        hist.plot.plot1d(sig_bkg_hist,
            overlay=sig_bkg_hist.axes()[0].name,
            ax=ax, clear=False, stack=True, line_opts=None,
            fill_opts=Plotter.stack_fill_opts,
            error_opts=Plotter.stack_error_opts,
        )

    if cen_sb:
        cen_data_minus_top = Plotter.data_minus_top(cen_sb['histo'])
            # find scale to get shape of data-top and then scale to sig norm
        sb_to_sig_scale = sig_norm/cen_data_minus_top.sum(*[ax.name for ax in cen_data_minus_top.axes()]).values()[()]
        cen_data_minus_top.scale(sb_to_sig_scale)

        hist.plot.plot1d(cen_data_minus_top, ax=ax, clear=False,
            line_opts={'color' : cen_sb['color']},
        )
        cols_to_label.update({cen_sb['color']: cen_sb['label']})

    if up_sb:
        up_data_minus_top = Plotter.data_minus_top(up_sb['histo'])
            # find scale to get shape of data-top and then scale to sig norm
        sb_to_sig_scale = sig_norm/up_data_minus_top.sum(*[ax.name for ax in up_data_minus_top.axes()]).values()[()]
        up_data_minus_top.scale(sb_to_sig_scale)

        hist.plot.plot1d(up_data_minus_top, ax=ax, clear=False,
            line_opts={'color' : up_sb['color']},
        )
        cols_to_label.update({up_sb['color']: up_sb['label']})

    if dw_sb:
        dw_data_minus_top = Plotter.data_minus_top(dw_sb['histo'])
            # find scale to get shape of data-top and then scale to sig norm
        sb_to_sig_scale = sig_norm/dw_data_minus_top.sum(*[ax.name for ax in dw_data_minus_top.axes()]).values()[()]
        dw_data_minus_top.scale(sb_to_sig_scale)

        hist.plot.plot1d(dw_data_minus_top, ax=ax, clear=False,
            line_opts={'color' : dw_sb['color']},
        )
        cols_to_label.update({dw_sb['color']: dw_sb['label']})


        ## set legend and corresponding colors
    handles, labels = ax.get_legend_handles_labels()
    for idx, sample in enumerate(labels):
        if sample == 'data' or sample == 'Observed': continue
        if isinstance(handles[idx], matplotlib.lines.Line2D):
            labels[idx] = cols_to_label[handles[idx].get_color()]
        else:
            facecolor, legname = plt_tools.get_styles(sample, Plotter.hstyles)
            handles[idx].set_facecolor(facecolor)
            labels[idx] = legname

    handles, def_labels = ax.get_legend_handles_labels()
    ax.legend(handles,labels, loc='upper right')

    ax.autoscale()
    ax.set_ylabel('Events')
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


def plot_bkg_mc_dd_comp(signal={}, cen_sb={}, up_sb={}, dw_sb={}, **opts):
    if not signal: raise ValueError("You need a distribution from the signal region")

    fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
    #fig, ax = plt.subplots()
    fig.subplots_adjust(hspace=.07)

        # get opts
    xlimits = opts.get('xlims')
    xlabel = opts.get('xtitle')
    jmult = opts.get('jmult')
    pltdir = opts.get('pltdir')
    vlines = opts.get('vlines')
    overflow = opts.get("overflow", "none")
    cols_to_label = {}

    #set_trace()
    if signal:
            # get EWK+QCD MC
        sig_bkg_hist = signal['histo'][Plotter.bkg_samples]
        sig_norm = sig_bkg_hist.sum(*[ax.name for ax in sig_bkg_hist.axes()]).values()[()]
        hist.plot.plot1d(sig_bkg_hist,
            overlay=sig_bkg_hist.axes()[0].name,
            ax=ax, clear=False, stack=True, line_opts=None,
            fill_opts=Plotter.stack_fill_opts,
            error_opts=Plotter.stack_error_opts,
            overflow=overflow,
            overlay_overflow=overflow,
        )

    if cen_sb:
        cen_data_minus_top = Plotter.data_minus_top(cen_sb['histo'])
            # find scale to get shape of data-top and then scale to sig norm
        sb_to_sig_scale = sig_norm/cen_data_minus_top.sum(*[ax.name for ax in cen_data_minus_top.axes()]).values()[()]
        cen_data_minus_top.scale(sb_to_sig_scale)

        hist.plot.plot1d(cen_data_minus_top, ax=ax, clear=False,
            line_opts={'color' : cen_sb['color']},
            overflow=overflow,
            overlay_overflow=overflow,
        )
        cols_to_label.update({cen_sb['color']: cen_sb['label']})

    if up_sb:
        up_data_minus_top = Plotter.data_minus_top(up_sb['histo'])
            # find scale to get shape of data-top and then scale to sig norm
        sb_to_sig_scale = sig_norm/up_data_minus_top.sum(*[ax.name for ax in up_data_minus_top.axes()]).values()[()]
        up_data_minus_top.scale(sb_to_sig_scale)

        hist.plot.plot1d(up_data_minus_top, ax=ax, clear=False,
            line_opts={'color' : up_sb['color']},
            overflow=overflow,
            overlay_overflow=overflow,
        )
        cols_to_label.update({up_sb['color']: up_sb['label']})

            # plot num/denom ratio
        hist.plotratio(
            up_data_minus_top, cen_data_minus_top, error_opts={'marker': '.', 'markersize': 10., 'color': up_sb['color'], 'elinewidth': 1},
            unc='num', clear=False, ax=rax, guide_opts={},
            overflow=overflow,
            #overlay_overflow=overflow,
        )

    if dw_sb:
        dw_data_minus_top = Plotter.data_minus_top(dw_sb['histo'])
            # find scale to get shape of data-top and then scale to sig norm
        sb_to_sig_scale = sig_norm/dw_data_minus_top.sum(*[ax.name for ax in dw_data_minus_top.axes()]).values()[()]
        dw_data_minus_top.scale(sb_to_sig_scale)

        hist.plot.plot1d(dw_data_minus_top, ax=ax, clear=False,
            line_opts={'color' : dw_sb['color']},
            overflow=overflow,
            overlay_overflow=overflow,
        )
        cols_to_label.update({dw_sb['color']: dw_sb['label']})

            # plot num/denom ratio
        hist.plotratio(
            dw_data_minus_top, cen_data_minus_top, error_opts={'marker': '.', 'markersize': 10., 'color': dw_sb['color'], 'elinewidth': 1},
            unc='num', clear=False, ax=rax, guide_opts={},
            overflow=overflow,
            #overlay_overflow=overflow,
        )

        ## set legend and corresponding colors
    handles, labels = ax.get_legend_handles_labels()
    for idx, sample in enumerate(labels):
        if sample == 'data' or sample == 'Observed': continue
        if isinstance(handles[idx], matplotlib.lines.Line2D):
            labels[idx] = cols_to_label[handles[idx].get_color()]
        else:
            facecolor, legname = plt_tools.get_styles(sample, Plotter.hstyles)
            handles[idx].set_facecolor(facecolor)
            labels[idx] = legname

    handles, def_labels = ax.get_legend_handles_labels()
    ax.legend(handles,labels, loc='upper right', ncol=2)

    ax.autoscale()
    ax.set_ylabel('Events')
    ax.set_ylim(0, None)
    ax.set_xlabel(None)
    ax.set_xlim(xlimits)

    rax.autoscale()
    rax.set_ylabel('Sys/Cen')
    rax.set_ylim(max(0.0, rax.get_ylim()[0]), min(1.5, rax.get_ylim()[1]))
    #rax.set_ylim(max(0.0, rax.get_ylim()[0]), min(5.0, rax.get_ylim()[1]))
    rax.set_xlim(xlimits)
    rax.set_xlabel(xlabel)

    if vlines is not None:
        #set_trace()
            # draw vertical lines separating ctstar bins
        for vline in vlines:
            ax.axvline(vline, color='k', linestyle='--')
            rax.axvline(vline, color='k', linestyle='--')

    ax.text(
        0.02, 0.94, "%s, %s" % (objtypes['Lep'][args.lepton], jet_mults[jmult]),
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
        #set_trace()
        histo = hdict[hname][:, :, 'nosys', :].integrate('sys') # process, sys, jmult, btag
        #set_trace()
    
        if histo.dense_dim() == 1:
            xtitle, rebinning, x_lims, withData = variables[hname]
            if rebinning != 1:
                xaxis_name = histo.dense_axes()[0].name
                histo = histo.rebin(xaxis_name, rebinning)
    
            ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
            #for jmult in ['3Jets']:
            #for jmult in ['4PJets']:
            for jmult in sorted(set([key[2] for key in histo.values().keys()])):
                #for btagregion in ['btagPass']:
                for btagregion in sorted(set([key[1] for key in histo.values().keys()])):
                    print(*[jmult, btagregion, hname], sep=', ')
                    pltdir = os.path.join(outdir, args.lepton, jmult, btagregion)
                    if not os.path.isdir(pltdir):
                        os.makedirs(pltdir)
   
                    fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                    fig.subplots_adjust(hspace=.07)
    
                    #set_trace()
                    hslice = histo[:, btagregion, jmult].integrate('jmult').integrate('btag')
   
                    if hname == 'Lep_iso':
                        x_lims = (0., 0.15) if (args.lepton == 'Muon') else (0., 0.1)
    
                    mc_opts = {
                        #'mcorder' : ['QCD', 'EWK', 'singlet', 'ttJets'] if not ttJets_cats else ['QCD', 'EWK', 'singlet', 'ttJets_other', 'ttJets_unmatchable', 'ttJets_matchable', 'ttJets_right'],
                        #"maskData" : not withData,
                        "overflow" : "under" if 'DeepCSV_bDisc' in hname else "none",
                    }

                    #set_trace() 
                    Plotter.plot_stack1d(ax, rax, hslice, xlabel=xtitle, xlimits=x_lims, **mc_opts)
    
                    #set_trace() 
                    if hname == 'Jets_njets':
                        print(jmult)
                        yields_txt, yields_json = Plotter.get_samples_yield_and_frac(hslice, data_lumi_year['%ss' % args.lepton]/1000., promptmc=True)
                        frac_name = '%s_yields_and_fracs' % '_'.join([jmult, args.lepton, btagregion])
                        plt_tools.print_table(yields_txt, filename='%s/%s.txt' % (pltdir, frac_name), print_output=True)
                        print('%s/%s.txt written' % (pltdir, frac_name))
                        with open('%s/%s.json' % (pltdir, frac_name), 'w') as out:
                            out.write(prettyjson.dumps(yields_json))
    
                        # add lepton/jet multiplicity label
                    #set_trace()
                    ax.text(
                        0.02, 0.88, "%s, %s\n%s" % (objtypes['Lep'][args.lepton], jet_mults[jmult], btag_cats[btagregion]),
                        fontsize=rcParams['font.size']*0.75, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
                    )
                    hep.cms.label(ax=ax, data=withData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
    
                    #set_trace()
                    figname = os.path.join(pltdir, '_'.join([jmult, args.lepton, btagregion, hname]))
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
            #for jmult in ['4PJets']:
            #   for btagregion in ['btagPass']:
            for jmult in sorted(set([key[2] for key in histo.values().keys()])):
                #for btagregion in ['btagPass']:
                for btagregion in sorted(set([key[1] for key in histo.values().keys()])):
                    print(*[jmult, btagregion, hname], sep=', ')
                    pltdir = os.path.join(outdir, args.lepton, jmult, btagregion)
                    if not os.path.isdir(pltdir):
                        os.makedirs(pltdir)
  
                    if (btagregion == 'btagPass') and (hname == 'mtt_vs_tlep_ctstar_abs'): 
                        withData = False
                    else:
                        withData = variables[hname][-1]

                    hslice = histo[:, btagregion, jmult].integrate('jmult').integrate('btag')
   
                    mc_opts = {
                        #'maskData' : not withData
                        #'mcorder' : ['QCD', 'EWK', 'singlet', 'ttJets'] if not ttJets_cats else ['QCD', 'EWK', 'singlet', 'ttJets_other', 'ttJets_unmatchable', 'ttJets_matchable', 'ttJets_right']
                    }

                    if 'DeepCSV_bDisc' in hname:
                        underflow_hist = hslice[:, :, :0.00].integrate('deepcsv_bdisc', overflow='under')
                        normal_reg_hist = hslice.integrate('deepcsv_bdisc')

                            # plot 1D projection of underflow bin
                        fig_und, (ax_und, rax_und) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                        fig_und.subplots_adjust(hspace=.07)

                        Plotter.plot_stack1d(ax_und, rax_und, underflow_hist, xlabel='%s DeepCSV underflow' % xtitle, xlimits=x_lims, **mc_opts)
    
                            # add lepton/jet multiplicity label
                        ax_und.text(
                            0.02, 0.88, "%s, %s\n%s" % (objtypes['Lep'][args.lepton], jet_mults[jmult], btag_cats[btagregion]),
                            fontsize=rcParams['font.size']*0.75, horizontalalignment='left', verticalalignment='bottom', transform=ax_und.transAxes
                        )
                        hep.cms.label(ax=ax_und, data=withData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
    
                        #set_trace()
                        und_figname = os.path.join(pltdir, '_'.join([jmult, args.lepton, btagregion, hname, 'DeepCSV_underflow', 'proj']))
                        fig_und.savefig(und_figname)
                        print('%s written' % und_figname)
                        plt.close()
 
                            # plot 1D projection of normal bin region
                        fig_norm, (ax_norm, rax_norm) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                        fig_norm.subplots_adjust(hspace=.07)

                        Plotter.plot_stack1d(ax_norm, rax_norm, normal_reg_hist, xlabel='%s DeepCSV normal region' % xtitle, xlimits=x_lims, **mc_opts)
    
                            # add lepton/jet multiplicity label
                        ax_norm.text(
                            0.02, 0.88, "%s, %s\n%s" % (objtypes['Lep'][args.lepton], jet_mults[jmult], btag_cats[btagregion]),
                            fontsize=rcParams['font.size']*0.75, horizontalalignment='left', verticalalignment='bottom', transform=ax_norm.transAxes
                        )
                        hep.cms.label(ax=ax_norm, data=withData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
    
                        #set_trace()
                        norm_figname = os.path.join(pltdir, '_'.join([jmult, args.lepton, btagregion, hname, 'DeepCSV_normal', 'proj']))
                        fig_norm.savefig(norm_figname)
                        print('%s written' % norm_figname)
                        plt.close()
 
                        continue


                        # make 1D projection along dense axes
                    for dax in range(2):
                        if dax == 0:
                            xlabel = ytitle
                            xlimits = y_lims
                            figname = os.path.join(pltdir, '_'.join([jmult, args.lepton, btagregion, yaxis_name, 'proj']))
                        else:
                            xlabel = '%s [GeV]' % xtitle
                            xlimits = x_lims
                            figname = os.path.join(pltdir, '_'.join([jmult, args.lepton, btagregion, xaxis_name, 'proj']))

                        fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                        fig.subplots_adjust(hspace=.07)

                        hproj = hslice.integrate(hslice.dense_axes()[dax].name)
    
                        Plotter.plot_stack1d(ax, rax, hproj, xlabel=xlabel, xlimits=xlimits, **mc_opts)
    
                            # add lepton/jet multiplicity label
                        ax.text(
                            0.02, 0.88, "%s, %s\n%s" % (objtypes['Lep'][args.lepton], jet_mults[jmult], btag_cats[btagregion]),
                            fontsize=rcParams['font.size']*0.75, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
                        )
                        hep.cms.label(ax=ax, data=withData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
    
                        #set_trace()
                        fig.savefig(figname)
                        print('%s written' % figname)
                        plt.close()
 
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
                        0.02, 0.88, "%s, %s\n%s" % (objtypes['Lep'][args.lepton], jet_mults[jmult], btag_cats[btagregion]),
                        fontsize=rcParams['font.size']*0.75, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
                    )
                    hep.cms.label(ax=ax, data=withData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
    
                    #set_trace()
                    figname = os.path.join(pltdir, '_'.join([jmult, args.lepton, btagregion, hname]))
                    fig.savefig(figname)
                    print('%s written' % figname)
                    plt.close()



#if args.dmp:
#    for hname in variables.keys():
#        if hname not in hdict.keys():
#            raise ValueError("%s not found in file" % hname)
#        if (hname != 'mtt_vs_tlep_ctstar_abs'): continue
#        histo = hdict[hname][:, :, 'nosys', :].integrate('sys') # process, sys, jmult, btag
#        #set_trace()
#        
#        xtitle, ytitle, xrebinning, yrebinning, x_lims, y_lims, withData, plot_xlabels, plot_ylabels = variables[hname]
#        vlines = [(len(xrebinning)-1)*ybin for ybin in range(1, len(yrebinning)-1)] if histo.dense_dim() == 2 else None
#        nbins = (len(xrebinning)-1)*(len(yrebinning)-1) if histo.dense_dim() == 2 else None
#    
#        xaxis_name = histo.dense_axes()[0].name
#        yaxis_name = histo.dense_axes()[1].name
#            ## rebin x axis
#        if isinstance(xrebinning, np.ndarray):
#            new_xbins = hist.Bin(xaxis_name, xaxis_name, xrebinning)
#        elif isinstance(xrebinning, float) or isinstance(xrebinning, int):
#            new_xbins = xrebinning
#        histo = histo.rebin(xaxis_name, new_xbins)
#            ## rebin y axis
#        if isinstance(yrebinning, np.ndarray):
#            new_ybins = hist.Bin(yaxis_name, yaxis_name, yrebinning)
#        elif isinstance(yrebinning, float) or isinstance(yrebinning, int):
#            new_ybins = yrebinning
#        histo = histo.rebin(yaxis_name, new_ybins)
#
#        #set_trace()   
#        ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
#        #for jmult in ['3Jets']:
#        #for jmult in ['4PJets']:
#        for jmult in sorted(set([key[2] for key in histo.values().keys()])):
#            shape_dir = os.path.join(outdir, args.lepton, jmult, 'QCD_Est', 'Shapes')
#            if not os.path.isdir(shape_dir): os.makedirs(shape_dir)
#                # get sideband and signal region hists (hardcoded)
#            cen_sb_histo = histo[:, '0p1to0p4', jmult].integrate('jmult').integrate('btag')
#            up_sb_histo = histo[:, '0p2to0p5', jmult].integrate('jmult').integrate('btag')
#            dw_sb_histo = histo[:, '0p0to0p2', jmult].integrate('jmult').integrate('btag')
#            sig_histo = histo[:, 'btagPass', jmult].integrate('jmult').integrate('btag')
#                # linearize sideband and signal region hists (hardcoded)
#            cen_sb_lin = Plotter.linearize_hist(cen_sb_histo)
#            up_sb_lin = Plotter.linearize_hist(up_sb_histo)
#            dw_sb_lin = Plotter.linearize_hist(dw_sb_histo)
#            sig_lin = Plotter.linearize_hist(sig_histo)
#
#               ## data-prompt
#                   ## plot projection onto mtt
#            plot_qcd_dmp_shape(cen_sb={'histo' : cen_sb_histo.integrate('ctstar_abs'), 'label' : 'Cen (0.1-0.4)', 'color' : 'k'},
#                up_sb={'histo' : up_sb_histo.integrate('ctstar_abs'), 'label' : 'Up (0.2-0.5)', 'color' : 'r'},
#                dw_sb={'histo' : dw_sb_histo.integrate('ctstar_abs'), 'label' : 'Down (0.0-0.2)', 'color' : 'b'},
#               **{'xtitle':'%s [GeV]' % xtitle, 'xlims':x_lims, 'pltdir':shape_dir, 'jmult':jmult, 'fname': 'QCD_DMP_Shapes_%s_%s_mtt' % (args.lepton, jmult)})
#                   ## plot projection onto ctstar
#            plot_qcd_dmp_shape(cen_sb={'histo' : cen_sb_histo.integrate('mtt'), 'label' : 'Cen (0.1-0.4)', 'color' : 'k'},
#                up_sb={'histo' : up_sb_histo.integrate('mtt'), 'label' : 'Up (0.2-0.5)', 'color' : 'r'},
#                dw_sb={'histo' : dw_sb_histo.integrate('mtt'), 'label' : 'Down (0.0-0.2)', 'color' : 'b'},
#               **{'xtitle':ytitle, 'xlims':y_lims, 'pltdir':shape_dir, 'hname': 'tlep_ctstar_abs', 'jmult':jmult, 'fname': 'QCD_DMP_Shapes_%s_%s_tlep_ctstar_abs' % (args.lepton, jmult)})
#                   ## plot linearized 2D mtt ctstar dist
#            plot_qcd_dmp_shape(cen_sb={'histo' : cen_sb_lin, 'label' : 'Cen (0.1-0.4)', 'color' : 'k'},
#                up_sb={'histo' : up_sb_lin, 'label' : 'Up (0.2-0.5)', 'color' : 'r'},
#                dw_sb={'histo' : dw_sb_lin, 'label' : 'Down (0.0-0.2)', 'color' : 'b'},
#               **{'xtitle':'%s $\otimes$ %s' % (xtitle, ytitle), 'xlims':(0, nbins), 'pltdir':shape_dir, 'hname': hname, 'jmult':jmult, 'vlines':vlines, 'fname': 'QCD_DMP_Shapes_%s_%s_%s' % (args.lepton, jmult, hname)})
#
#               ## data-(single top/ttbar)
#                   ## plot projection onto mtt
#            plot_ewk_qcd_dmtop_shape(cen_sb={'histo' : cen_sb_histo.integrate('ctstar_abs'), 'label' : 'Cen (0.1-0.4)', 'color' : 'k'},
#                up_sb={'histo' : up_sb_histo.integrate('ctstar_abs'), 'label' : 'Up (0.2-0.5)', 'color' : 'r'},
#                dw_sb={'histo' : dw_sb_histo.integrate('ctstar_abs'), 'label' : 'Down (0.0-0.2)', 'color' : 'b'},
#               **{'xtitle':'%s [GeV]' % xtitle, 'xlims':x_lims, 'pltdir':shape_dir, 'jmult':jmult, 'fname': 'EWK_plus_QCD_DD_Shapes_%s_%s_mtt' % (args.lepton, jmult)})
#                   ## plot projection onto ctstar
#            plot_ewk_qcd_dmtop_shape(cen_sb={'histo' : cen_sb_histo.integrate('mtt'), 'label' : 'Cen (0.1-0.4)', 'color' : 'k'},
#                up_sb={'histo' : up_sb_histo.integrate('mtt'), 'label' : 'Up (0.2-0.5)', 'color' : 'r'},
#                dw_sb={'histo' : dw_sb_histo.integrate('mtt'), 'label' : 'Down (0.0-0.2)', 'color' : 'b'},
#               **{'xtitle':ytitle, 'xlims':y_lims, 'pltdir':shape_dir, 'hname': 'tlep_ctstar_abs', 'jmult':jmult, 'fname': 'EWK_plus_QCD_DD_Shapes_%s_%s_tlep_ctstar_abs' % (args.lepton, jmult)})
#                   ## plot linearized 2D mtt ctstar dist
#            plot_ewk_qcd_dmtop_shape(cen_sb={'histo' : cen_sb_lin, 'label' : 'Cen (0.1-0.4)', 'color' : 'k'},
#                up_sb={'histo' : up_sb_lin, 'label' : 'Up (0.2-0.5)', 'color' : 'r'},
#                dw_sb={'histo' : dw_sb_lin, 'label' : 'Down (0.0-0.2)', 'color' : 'b'},
#               **{'xtitle':'%s $\otimes$ %s' % (xtitle, ytitle), 'xlims':(0, nbins), 'pltdir':shape_dir, 'hname': hname, 'jmult':jmult, 'vlines':vlines, 'fname': 'EWK_plus_QCD_DD_Shapes_%s_%s_%s' % (args.lepton, jmult, hname)})
#
#            # make plots comparing data-driven background (data-ttbar-st) in sidebands to EWK+QCD MC in signal
#                   ## plot projection onto mtt
#            plot_ewk_qcd_cont(cen_sb={'histo' : cen_sb_histo.integrate('ctstar_abs'), 'label' : 'Cen (0.1-0.4)', 'color' : 'k'},
#                up_sb={'histo' : up_sb_histo.integrate('ctstar_abs'), 'label' : 'Up (0.2-0.5)', 'color' : 'r'},
#                dw_sb={'histo' : dw_sb_histo.integrate('ctstar_abs'), 'label' : 'Down (0.0-0.2)', 'color' : 'b'},
#                signal={'histo' : sig_histo.integrate('ctstar_abs')},
#               **{'xtitle':'%s [GeV]' % xtitle, 'xlims':x_lims, 'pltdir':shape_dir, 'jmult':jmult, 'fname': 'EWK_QCD_Comparison_%s_%s_mtt' % (args.lepton, jmult)})
#                   ## plot projection onto ctstar
#            plot_ewk_qcd_cont(cen_sb={'histo' : cen_sb_histo.integrate('mtt'), 'label' : 'Cen (0.1-0.4)', 'color' : 'k'},
#                up_sb={'histo' : up_sb_histo.integrate('mtt'), 'label' : 'Up (0.2-0.5)', 'color' : 'r'},
#                dw_sb={'histo' : dw_sb_histo.integrate('mtt'), 'label' : 'Down (0.0-0.2)', 'color' : 'b'},
#                signal={'histo' : sig_histo.integrate('mtt')},
#               **{'xtitle':ytitle, 'xlims':y_lims, 'pltdir':shape_dir, 'hname': 'tlep_ctstar_abs', 'jmult':jmult, 'fname': 'EWK_QCD_Comparison_%s_%s_tlep_ctstar_abs' % (args.lepton, jmult)})
#                   ## plot linearized 2D mtt ctstar dist
#            plot_ewk_qcd_cont(cen_sb={'histo' : cen_sb_lin, 'label' : 'Cen (0.1-0.4)', 'color' : 'k'},
#                up_sb={'histo' : up_sb_lin, 'label' : 'Up (0.2-0.5)', 'color' : 'r'},
#                dw_sb={'histo' : dw_sb_lin, 'label' : 'Down (0.0-0.2)', 'color' : 'b'},
#                signal={'histo' : sig_lin},
#               **{'xtitle':'%s $\otimes$ %s' % (xtitle, ytitle), 'xlims':(0, nbins), 'pltdir':shape_dir, 'hname': hname, 'jmult':jmult, 'vlines':vlines, 'fname': 'EWK_QCD_Comparison_%s_%s_%s' % (args.lepton, jmult, hname)})


if args.bkg_shapes:
    for hname in variables.keys():
        if hname not in hdict.keys():
            raise ValueError("%s not found in file" % hname)
        histo = hdict[hname][:, :, 'nosys', :].integrate('sys') # process, sys, jmult, btag
        #set_trace()
        
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
                xrebinning = np.array([300.0, 400.0, 500.0, 600.0, 800.0, 1000.0, 1500.0])
                x_lims = (xrebinning[0], xrebinning[-1])

            if 'DeepCSV_bDisc' not in hname:
                vlines = [(len(xrebinning)-1)*ybin for ybin in range(1, len(yrebinning)-1)] if histo.dense_dim() == 2 else None
                nbins = (len(xrebinning)-1)*(len(yrebinning)-1) if histo.dense_dim() == 2 else None
    
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

        #set_trace()   
        if hname == 'Lep_iso':
            x_lims = (0., 0.15) if (args.lepton == 'Muon') else (0., 0.1)
    
        ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
        #for jmult in ['3Jets']:
        #for jmult in ['4PJets']:
        for jmult in sorted(set([key[2] for key in histo.values().keys()])):
            shape_dir = os.path.join(outdir, args.lepton, jmult, 'BKG_Est_orthog', 'Shapes', 'DMT')
            if not os.path.isdir(shape_dir): os.makedirs(shape_dir)
            comp_dir = os.path.join(outdir, args.lepton, jmult, 'BKG_Est_orthog', 'Shapes', 'MC_DD_Comp')
            if not os.path.isdir(comp_dir): os.makedirs(comp_dir)

                # get sideband and signal region hists
            cen_sb_histo = histo[:, btag_reg_names_dict['Central']['reg'], jmult].integrate('jmult').integrate('btag')
            up_sb_histo = histo[:, btag_reg_names_dict['Up']['reg'], jmult].integrate('jmult').integrate('btag')
            dw_sb_histo = histo[:, btag_reg_names_dict['Down']['reg'], jmult].integrate('jmult').integrate('btag')
            sig_histo = histo[:, btag_reg_names_dict['Signal']['reg'], jmult].integrate('jmult').integrate('btag')

            if histo.dense_dim() == 1:
                # make data-(single top/ttbar) shape plots
                plot_ewk_qcd_dmtop_shape(
                    cen_sb={'histo' : cen_sb_histo, 'label' : btag_reg_names_dict['Central']['label'], 'color' : btag_reg_names_dict['Central']['color']},
                    up_sb={'histo' : up_sb_histo, 'label' : btag_reg_names_dict['Up']['label'], 'color' : btag_reg_names_dict['Up']['color']},
                    dw_sb={'histo' : dw_sb_histo, 'label' : btag_reg_names_dict['Down']['label'], 'color' : btag_reg_names_dict['Down']['color']},
                    signal={'histo' : sig_histo},
                   **{'xtitle':xtitle, 'xlims':x_lims, 'pltdir':shape_dir, 'hname': hname, 'jmult':jmult, 
                        'fname': 'BKG_DD_Shapes_%s_%s_%s' % (args.lepton, jmult, hname), "overflow" : "under" if 'DeepCSV_bDisc' in hname else "none"})

                # make plots comparing data-driven background (data-ttbar-st) in sidebands to EWK+QCD MC in signal
                plot_bkg_mc_dd_comp(
                    cen_sb={'histo' : cen_sb_histo, 'label' : btag_reg_names_dict['Central']['label'], 'color' : btag_reg_names_dict['Central']['color']},
                    up_sb={'histo' : up_sb_histo, 'label' : btag_reg_names_dict['Up']['label'], 'color' : btag_reg_names_dict['Up']['color']},
                    dw_sb={'histo' : dw_sb_histo, 'label' : btag_reg_names_dict['Down']['label'], 'color' : btag_reg_names_dict['Down']['color']},
                    signal={'histo' : sig_histo},
                   **{'xtitle':xtitle, 'xlims':x_lims, 'pltdir':comp_dir, 'hname': hname, 'jmult':jmult,
                        'fname': 'BKG_MC_DD_Comparison_%s_%s_%s' % (args.lepton, jmult, hname), "overflow" : "under" if 'DeepCSV_bDisc' in hname else "none"})

            if histo.dense_dim() == 2:
                if 'DeepCSV_bDisc' in hname:
                        ## plot 1D projection of underflow bin
                    plot_ewk_qcd_dmtop_shape(
                        cen_sb={'histo' : cen_sb_histo[:, :, :0.00].integrate('deepcsv_bdisc', overflow='under'), 'label' : btag_reg_names_dict['Central']['label'], 'color' : btag_reg_names_dict['Central']['color']},
                        up_sb={'histo' : up_sb_histo[:, :, :0.00].integrate('deepcsv_bdisc', overflow='under'), 'label' : btag_reg_names_dict['Up']['label'], 'color' : btag_reg_names_dict['Up']['color']},
                        dw_sb={'histo' : dw_sb_histo[:, :, :0.00].integrate('deepcsv_bdisc', overflow='under'), 'label' : btag_reg_names_dict['Down']['label'], 'color' : btag_reg_names_dict['Down']['color']},
                        signal={'histo' : sig_histo[:, :, :0.00].integrate('deepcsv_bdisc', overflow='under')},
                       **{'xtitle':'%s DeepCSV underflow' % xtitle, 'xlims':x_lims, 'pltdir':shape_dir, 'jmult':jmult, 'fname': '_'.join(['BKG_DD_Shapes', args.lepton, jmult, hname, 'DeepCSV_underflow', 'proj'])})

                    plot_bkg_mc_dd_comp(
                        cen_sb={'histo' : cen_sb_histo[:, :, :0.00].integrate('deepcsv_bdisc', overflow='under'), 'label' : btag_reg_names_dict['Central']['label'], 'color' : btag_reg_names_dict['Central']['color']},
                        up_sb={'histo' : up_sb_histo[:, :, :0.00].integrate('deepcsv_bdisc', overflow='under'), 'label' : btag_reg_names_dict['Up']['label'], 'color' : btag_reg_names_dict['Up']['color']},
                        dw_sb={'histo' : dw_sb_histo[:, :, :0.00].integrate('deepcsv_bdisc', overflow='under'), 'label' : btag_reg_names_dict['Down']['label'], 'color' : btag_reg_names_dict['Down']['color']},
                        signal={'histo' : sig_histo[:, :, :0.00].integrate('deepcsv_bdisc', overflow='under')},
                       **{'xtitle':'%s DeepCSV underflow' % xtitle, 'xlims':x_lims, 'pltdir':comp_dir, 'jmult':jmult, 'fname': '_'.join(['BKG_MC_DD_Comparison', args.lepton, jmult, hname, 'DeepCSV_underflow', 'proj'])})

                        ## plot 1D projection of normal bin content
                    plot_ewk_qcd_dmtop_shape(
                        cen_sb={'histo' : cen_sb_histo.integrate('deepcsv_bdisc'), 'label' : btag_reg_names_dict['Central']['label'], 'color' : btag_reg_names_dict['Central']['color']},
                        up_sb={'histo' : up_sb_histo.integrate('deepcsv_bdisc'), 'label' : btag_reg_names_dict['Up']['label'], 'color' : btag_reg_names_dict['Up']['color']},
                        dw_sb={'histo' : dw_sb_histo.integrate('deepcsv_bdisc'), 'label' : btag_reg_names_dict['Down']['label'], 'color' : btag_reg_names_dict['Down']['color']},
                        signal={'histo' : sig_histo.integrate('deepcsv_bdisc')},
                       **{'xtitle':'%s DeepCSV normal region' % xtitle, 'xlims':x_lims, 'pltdir':shape_dir, 'jmult':jmult, 'fname': '_'.join(['BKG_DD_Shapes', args.lepton, jmult, hname, 'DeepCSV_normal', 'proj'])})

                    plot_bkg_mc_dd_comp(
                        cen_sb={'histo' : cen_sb_histo.integrate('deepcsv_bdisc'), 'label' : btag_reg_names_dict['Central']['label'], 'color' : btag_reg_names_dict['Central']['color']},
                        up_sb={'histo' : up_sb_histo.integrate('deepcsv_bdisc'), 'label' : btag_reg_names_dict['Up']['label'], 'color' : btag_reg_names_dict['Up']['color']},
                        dw_sb={'histo' : dw_sb_histo.integrate('deepcsv_bdisc'), 'label' : btag_reg_names_dict['Down']['label'], 'color' : btag_reg_names_dict['Down']['color']},
                        signal={'histo' : sig_histo.integrate('deepcsv_bdisc')},
                       **{'xtitle':'%s DeepCSV normal region' % xtitle, 'xlims':x_lims, 'pltdir':comp_dir, 'jmult':jmult, 'fname': '_'.join(['BKG_MC_DD_Comparison', args.lepton, jmult, hname, 'DeepCSV_normal', 'proj'])})
 
                    continue


                    # linearize sideband and signal region hists (hardcoded)
                cen_sb_lin = Plotter.linearize_hist(cen_sb_histo)
                up_sb_lin = Plotter.linearize_hist(up_sb_histo)
                dw_sb_lin = Plotter.linearize_hist(dw_sb_histo)
                sig_lin = Plotter.linearize_hist(sig_histo)

                # make data-(single top/ttbar) shape plots
                plot_ewk_qcd_dmtop_shape(
                    cen_sb={'histo' : cen_sb_lin, 'label' : btag_reg_names_dict['Central']['label'], 'color' : btag_reg_names_dict['Central']['color']},
                    up_sb={'histo' : up_sb_lin, 'label' : btag_reg_names_dict['Up']['label'], 'color' : btag_reg_names_dict['Up']['color']},
                    dw_sb={'histo' : dw_sb_lin, 'label' : btag_reg_names_dict['Down']['label'], 'color' : btag_reg_names_dict['Down']['color']},
                   **{'xtitle':'%s $\otimes$ %s' % (xtitle, ytitle), 'xlims':(0, nbins), 'pltdir':shape_dir, 'hname': hname, 'jmult':jmult, 'vlines':vlines, 'fname': 'BKG_DD_Shapes_%s_%s_%s' % (args.lepton, jmult, hname)})

                # make 1D projection along dense axes
                for dax in range(2):
                    if dax == 0:
                        ## plot projection onto yaxis
                        plot_ewk_qcd_dmtop_shape(
                            cen_sb={'histo' : cen_sb_histo.integrate(cen_sb_histo.dense_axes()[dax].name), 'label' : btag_reg_names_dict['Central']['label'], 'color' : btag_reg_names_dict['Central']['color']},
                            up_sb={'histo' : up_sb_histo.integrate(up_sb_histo.dense_axes()[dax].name), 'label' : btag_reg_names_dict['Up']['label'], 'color' : btag_reg_names_dict['Up']['color']},
                            dw_sb={'histo' : dw_sb_histo.integrate(dw_sb_histo.dense_axes()[dax].name), 'label' : btag_reg_names_dict['Down']['label'], 'color' : btag_reg_names_dict['Down']['color']},
                            signal={'histo' : sig_histo.integrate(sig_histo.dense_axes()[dax].name)},
                           **{'xtitle':ytitle, 'xlims':y_lims, 'pltdir':shape_dir, 'jmult':jmult, 'fname': 'BKG_DD_Shapes_%s_%s_%s_proj' % (args.lepton, jmult, hname.split('_vs_')[1])})
                    else:
                        ## plot projection onto xaxis
                        plot_ewk_qcd_dmtop_shape(
                            cen_sb={'histo' : cen_sb_histo.integrate(cen_sb_histo.dense_axes()[dax].name), 'label' : btag_reg_names_dict['Central']['label'], 'color' : btag_reg_names_dict['Central']['color']},
                            up_sb={'histo' : up_sb_histo.integrate(up_sb_histo.dense_axes()[dax].name), 'label' : btag_reg_names_dict['Up']['label'], 'color' : btag_reg_names_dict['Up']['color']},
                            dw_sb={'histo' : dw_sb_histo.integrate(dw_sb_histo.dense_axes()[dax].name), 'label' : btag_reg_names_dict['Down']['label'], 'color' : btag_reg_names_dict['Down']['color']},
                            signal={'histo' : sig_histo.integrate(sig_histo.dense_axes()[dax].name)},
                           **{'xtitle':xtitle, 'xlims':x_lims, 'pltdir':shape_dir, 'jmult':jmult, 'fname': 'BKG_DD_Shapes_%s_%s_%s_proj' % (args.lepton, jmult, hname.split('_vs_')[0])})


                # make plots comparing data-driven background (data-ttbar-st) in sidebands to EWK+QCD MC in signal
                plot_bkg_mc_dd_comp(
                    cen_sb={'histo' : cen_sb_lin, 'label' : btag_reg_names_dict['Central']['label'], 'color' : btag_reg_names_dict['Central']['color']},
                    up_sb={'histo' : up_sb_lin, 'label' : btag_reg_names_dict['Up']['label'], 'color' : btag_reg_names_dict['Up']['color']},
                    dw_sb={'histo' : dw_sb_lin, 'label' : btag_reg_names_dict['Down']['label'], 'color' : btag_reg_names_dict['Down']['color']},
                    signal={'histo' : sig_lin},
                   **{'xtitle':'%s $\otimes$ %s' % (xtitle, ytitle), 'xlims':(0, nbins), 'pltdir':comp_dir, 'hname': hname, 'jmult':jmult, 'vlines':vlines, 'fname': 'BKG_MC_DD_Comparison_%s_%s_%s' % (args.lepton, jmult, hname)})

                # make 1D projection along dense axes
                for dax in range(2):
                    if dax == 0:
                        ## plot projection onto yaxis
                        plot_bkg_mc_dd_comp(
                            cen_sb={'histo' : cen_sb_histo.integrate(cen_sb_histo.dense_axes()[dax].name), 'label' : btag_reg_names_dict['Central']['label'], 'color' : btag_reg_names_dict['Central']['color']},
                            up_sb={'histo' : up_sb_histo.integrate(up_sb_histo.dense_axes()[dax].name), 'label' : btag_reg_names_dict['Up']['label'], 'color' : btag_reg_names_dict['Up']['color']},
                            dw_sb={'histo' : dw_sb_histo.integrate(dw_sb_histo.dense_axes()[dax].name), 'label' : btag_reg_names_dict['Down']['label'], 'color' : btag_reg_names_dict['Down']['color']},
                            signal={'histo' : sig_histo.integrate(sig_histo.dense_axes()[dax].name)},
                           **{'xtitle':ytitle, 'xlims':y_lims, 'pltdir':comp_dir, 'jmult':jmult, 'fname': 'BKG_MC_DD_Comparison_%s_%s_%s_proj' % (args.lepton, jmult, hname.split('_vs_')[1])})
                    else:
                        ## plot projection onto xaxis
                        plot_bkg_mc_dd_comp(
                            cen_sb={'histo' : cen_sb_histo.integrate(cen_sb_histo.dense_axes()[dax].name), 'label' : btag_reg_names_dict['Central']['label'], 'color' : btag_reg_names_dict['Central']['color']},
                            up_sb={'histo' : up_sb_histo.integrate(up_sb_histo.dense_axes()[dax].name), 'label' : btag_reg_names_dict['Up']['label'], 'color' : btag_reg_names_dict['Up']['color']},
                            dw_sb={'histo' : dw_sb_histo.integrate(dw_sb_histo.dense_axes()[dax].name), 'label' : btag_reg_names_dict['Down']['label'], 'color' : btag_reg_names_dict['Down']['color']},
                            signal={'histo' : sig_histo.integrate(sig_histo.dense_axes()[dax].name)},
                           **{'xtitle':xtitle, 'xlims':x_lims, 'pltdir':comp_dir, 'jmult':jmult, 'fname': 'BKG_MC_DD_Comparison_%s_%s_%s_proj' % (args.lepton, jmult, hname.split('_vs_')[0])})


 
# plot with EWK+QCD estimation
if args.bkg_est:
    systs_to_run = ['nosys']
    #systs_to_run = sorted(sys_to_name.keys()) if args.save_sys else ['nosys']
    print('\nPlotting distributions for systematics:\n\t', *systs_to_run)

    ## combine EWK and QCD processes into BKG groups
    proc_bin = hist.Cat("process", "process", sorting='placement')
    bkg_groups = {proc:proc for proc in process_groups.keys()}
    del bkg_groups['EWK']
    del bkg_groups['QCD']
    bkg_groups['BKG'] = ['EWK', 'QCD']
    
    #set_trace()
    #    ## save post qcd-est dists for all systematics
    #if (args.save_sys):
    #    import itertools
    #        # keys are (lepton, jmult, sys, hname)
    #    save_keys = sorted(itertools.product([args.lepton], ['3Jets', '4PJets'], [sys_to_name[sys] for sys in systs_to_run], [name for name in variables.keys()]))
    #    post_bkg_est = {key:None for key in save_keys}

    for hname in variables.keys():
        if 'DeepCSV_bDisc_' in hname: continue
        if hname not in hdict.keys():
            raise ValueError("%s not found in file" % hname)
        histo = hdict[hname].group("process", proc_bin, bkg_groups) # group by EWK+QCD into BKG
        #set_trace()

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
            "maskData" : not withData,
            "overflow" : "under" if hname == 'DeepCSV_bDisc' else "none",
        }

        vlines = [(len(xrebinning)-1)*ybin for ybin in range(1, len(yrebinning)-1)] if histo.dense_dim() == 2 else None
        nbins = (len(xrebinning)-1)*(len(yrebinning)-1) if histo.dense_dim() == 2 else None
    
        if hname == 'Lep_iso':
            x_lims = (0., 0.15) if args.lepton == 'Muon' else (0., 0.1)
    

            # loop over different jet multiplicities
        for jmult in sorted(set([key[3] for key in histo.values().keys()])):
        #for jmult in ['3Jets']:
        #for jmult in ['4PJets']:
            #set_trace()
            print(jmult)
                # get sideband and signal region hists
            cen_sb_histo = histo[:, btag_reg_names_dict['Central']['reg'], :, jmult].integrate('jmult').integrate('btag')
            sig_histo = histo[:, btag_reg_names_dict['Signal']['reg'], jmult].integrate('jmult').integrate('btag')

            if histo.dense_dim() == 1:
                cen_sb = histo[:, btag_reg_names_dict['Central']['reg'], 'nosys', jmult].integrate('btag').integrate('sys').integrate('jmult')
                sig_histo = histo[:, btag_reg_names_dict['Signal']['reg'], :, jmult].integrate('btag').integrate('jmult')
            else:
                cen_sb = Plotter.linearize_hist(histo[:, btag_reg_names_dict['Central']['reg'], 'nosys', jmult].integrate('btag').integrate('sys').integrate('jmult'))
                sig_histo = Plotter.linearize_hist(histo[:, btag_reg_names_dict['Signal']['reg'], :, jmult].integrate('btag').integrate('jmult'))

                ## plot systematic variations after bkg estimation is made
            #set_trace()
            for sys in systs_to_run:
                if sys not in histo.axis('sys')._sorted:
                    print('\n\n   Systematic %s not available, skipping\n\n' % sys)
                    continue
                print('BKG est:', jmult, sys, hname)
                #set_trace()

                for norm in ['SigMC']:
                    #set_trace()
                    bkg_est_histo = Plotter.BKG_Est(sig_reg=sig_histo, sb_reg=cen_sb, norm_type=norm, sys=sys, ignore_uncs=True)
                    #set_trace()

                        # save post bkg-est dist to dict
                    #if norm == 'Sideband':
                    #    if args.save_sys:
                    #        post_bkg_est[(args.lepton, jmult, sys_to_name[sys], hname)] = bkg_est_histo

                    if sys == 'nosys':                 
                        bkg_name = '%s%s_Norm' % (shape_reg, norm) if norm == 'Sideband' else '%s_Norm' % norm
                        bkg_dir = os.path.join(outdir, args.lepton, jmult, 'BKG_Est_orthog', bkg_name, sys_to_name[sys])
                        if not os.path.isdir(bkg_dir):
                            os.makedirs(bkg_dir)
    
                        fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                        fig.subplots_adjust(hspace=.07)
   
                        if hname == 'Jets_njets':
                            print('BKG est:', jmult, sys)
                            yields_txt, yields_json = Plotter.get_samples_yield_and_frac(bkg_est_histo, data_lumi_year['%ss' % args.lepton]/1000., sys=sys)
                            frac_name = '%s_yields_and_fracs_BKG_Est_orthog_%s' % ('_'.join([sys, jmult, args.lepton]), bkg_name)
                            plt_tools.print_table(yields_txt, filename='%s/%s.txt' % (bkg_dir, frac_name), print_output=True)
                            print('%s/%s.txt written' % (bkg_dir, frac_name))
                            with open('%s/%s.json' % (bkg_dir, frac_name), 'w') as out:
                                out.write(prettyjson.dumps(yields_json))

                        #set_trace() 
                        Plotter.plot_stack1d(ax, rax, bkg_est_histo, xlabel=xtitle, xlimits=x_lims, **mc_opts) if histo.dense_dim() == 1\
                            else Plotter.plot_stack1d(ax, rax, bkg_est_histo, xlabel='%s $\otimes$ %s' % (xtitle, ytitle), xlimits=(0, nbins), **mc_opts)
    
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
                        figname = os.path.join(bkg_dir, '%s_BKG_Est_orthog_%s' % ('_'.join([sys, jmult, args.lepton, hname]), bkg_name))
                        fig.savefig(figname)
                        print('%s written' % figname)
                        plt.close()


# plots with EWK+QCD estimation
if args.vary_norm:
    systs_to_run = ['nosys']
    print('\nPlotting distributions for systematics:\n\t', *systs_to_run)

    #set_trace()
    ## combine EWK and QCD processes into BKG groups
    proc_bin = hist.Cat("process", "process", sorting='placement')
    bkg_groups = {proc:proc for proc in process_groups.keys()}
    mc_procs = [proc for proc in bkg_groups.keys() if 'data' not in proc]
    for proc in mc_procs:
        del bkg_groups[proc]
    bkg_groups['ttJets'] = ['ttJets_right', 'ttJets_matchable', 'ttJets_unmatchable', 'ttJets_sl_tau', 'ttJets_other']
    bkg_groups['singlet'] = ['singlet']
    bkg_groups['BKG'] = ['EWK', 'QCD']

    for hname in variables.keys():
        if hname not in hdict.keys():
            raise ValueError("%s not found in file" % hname)
        if 'DeepCSV_bDisc_' in hname: continue
        histo = hdict[hname].group("process", proc_bin, bkg_groups) # group by EWK+QCD into BKG, all ttJets together
        #set_trace()

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
            "overflow" : "under" if hname == 'DeepCSV_bDisc' else "none",
        }

        vlines = [(len(xrebinning)-1)*ybin for ybin in range(1, len(yrebinning)-1)] if histo.dense_dim() == 2 else None
        nbins = (len(xrebinning)-1)*(len(yrebinning)-1) if histo.dense_dim() == 2 else None
    
        if hname == 'Lep_iso':
            x_lims = (0., 0.15) if args.lepton == 'Muon' else (0., 0.1)
    

            # loop over different jet multiplicities
        for jmult in sorted(set([key[3] for key in histo.values().keys()])):
        #for jmult in ['3Jets']:
        #for jmult in ['4PJets']:
            #set_trace()
            print(jmult)
                # get sideband and signal region hists (hardcoded)
            cen_sb_histo = histo[:, btag_reg_names_dict['Central']['reg'], :, jmult].integrate('jmult').integrate('btag')
            sig_histo = histo[:, btag_reg_names_dict['Signal']['reg'], jmult].integrate('jmult').integrate('btag')

            if histo.dense_dim() == 1:
                cen_sb = histo[:, btag_reg_names_dict['Central']['reg'], 'nosys', jmult].integrate('btag').integrate('sys').integrate('jmult')
                sig_histo = histo[:, btag_reg_names_dict['Signal']['reg'], :, jmult].integrate('btag').integrate('jmult')
            else:
                cen_sb = Plotter.linearize_hist(histo[:, btag_reg_names_dict['Central']['reg'], 'nosys', jmult].integrate('btag').integrate('sys').integrate('jmult'))
                sig_histo = Plotter.linearize_hist(histo[:, btag_reg_names_dict['Signal']['reg'], :, jmult].integrate('btag').integrate('jmult'))

                ## plot systematic variations after bkg estimation is made
            #set_trace()
            for sys in systs_to_run:
                if sys not in histo.axis('sys')._sorted:
                    print('\n\n   Systematic %s not available, skipping\n\n' % sys)
                    continue
                print('BKG est:', jmult, sys, hname)
                #set_trace()

                for norm in ['SigMC']:
                    #set_trace()
                    bkg_est_histo = Plotter.BKG_Est(sig_reg=sig_histo, sb_reg=cen_sb, norm_type=norm, sys=sys, ignore_uncs=True)
                    #set_trace()

                    # create bkg hist variations with same shape but +100%, -50% yields compared to nominal
                    bkg_shape, _ = Plotter.get_qcd_shape(bkg_est_histo[Plotter.bkg_mask].integrate('process'), overflow='none')
                    bkg_norm = np.sum(bkg_est_histo[Plotter.bkg_mask].integrate('process').values()[()])
                    st_tt_sumw = bkg_est_histo[Plotter.top_samples].integrate('process').values()[()]
                        # +100% bkg yield (plus other MC contribution for plotting)
                    bkg_up_sumw = bkg_shape*(bkg_norm*2) + st_tt_sumw
                        # -50% bkg yield (plus other MC contribution for plotting)
                    bkg_dw_sumw = bkg_shape*(bkg_norm*0.5) + st_tt_sumw

                    if sys == 'nosys':                 
                        bkg_name = '%s%s_Norm' % (shape_reg, norm) if norm == 'Sideband' else '%s_Norm' % norm
                        bkg_dir = os.path.join(outdir, args.lepton, jmult, 'BKG_Est_orthog', bkg_name, sys_to_name[sys], 'VARY_NORM')
                        if not os.path.isdir(bkg_dir):
                            os.makedirs(bkg_dir)
    
                        fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                        fig.subplots_adjust(hspace=.07)
   
                        if hname == 'Jets_njets':
                            print('BKG est:', jmult, sys)
                            yields_txt, yields_json = Plotter.get_samples_yield_and_frac(bkg_est_histo, data_lumi_year['%ss' % args.lepton]/1000., sys=sys)
                            frac_name = '%s_yields_and_fracs_BKG_Est_orthog_%s' % ('_'.join([sys, jmult, args.lepton]), bkg_name)
                            plt_tools.print_table(yields_txt, filename='%s/%s.txt' % (bkg_dir, frac_name), print_output=True)
                            print('%s/%s.txt written' % (bkg_dir, frac_name))
                            with open('%s/%s.json' % (bkg_dir, frac_name), 'w') as out:
                                out.write(prettyjson.dumps(yields_json))

                            # plot bkg variations
                        Plotter.plot_1D(values=bkg_up_sumw, bins=bkg_est_histo.dense_axes()[0].edges(), xlabel=xtitle, xlimits=x_lims, ax=ax, color='r', label='BKG +100%')
                        Plotter.plot_1D(values=bkg_dw_sumw, bins=bkg_est_histo.dense_axes()[0].edges(), xlabel=xtitle, xlimits=x_lims, ax=ax, color='b', label='BKG -50%')

                        #set_trace() 
                        Plotter.plot_stack1d(ax, rax, bkg_est_histo, xlabel=xtitle, xlimits=x_lims, **mc_opts) if histo.dense_dim() == 1\
                            else Plotter.plot_stack1d(ax, rax, bkg_est_histo, xlabel='%s $\otimes$ %s' % (xtitle, ytitle), xlimits=(0, nbins), **mc_opts)

                            # plot bkg variations
                                # get ratio arrays
                        if (hname != 'mtt_vs_tlep_ctstar_abs'):
                            bkg_up_ratio_vals, bkg_up_ratio_bins = Plotter.get_ratio_arrays(num_vals=bkg_est_histo[Plotter.data_samples].integrate('process').values()[()], denom_vals=bkg_up_sumw, input_bins=bkg_est_histo.dense_axes()[0].edges())
                            bkg_dw_ratio_vals, bkg_dw_ratio_bins = Plotter.get_ratio_arrays(num_vals=bkg_est_histo[Plotter.data_samples].integrate('process').values()[()], denom_vals=bkg_dw_sumw, input_bins=bkg_est_histo.dense_axes()[0].edges())
                            rax.step(bkg_up_ratio_bins, bkg_up_ratio_vals, where='post', **{'linestyle' : '-', 'color' : 'r'})
                            rax.step(bkg_dw_ratio_bins, bkg_dw_ratio_vals, where='post', **{'linestyle' : '-', 'color' : 'b'})
    
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
                        figname = os.path.join(bkg_dir, '%s_BKG_Est_orthog_VARY_NORM_%s' % ('_'.join([sys, jmult, args.lepton, hname]), bkg_name))
                        fig.savefig(figname)
                        print('%s written' % figname)
                        plt.close()


if args.vary_shape:
    systs_to_run = ['nosys']
    print('\nPlotting distributions for systematics:\n\t', *systs_to_run)

    #set_trace()
    ## combine EWK and QCD processes into BKG groups
    proc_bin = hist.Cat("process", "process", sorting='placement')
    bkg_groups = {proc:proc for proc in process_groups.keys()}
    mc_procs = [proc for proc in bkg_groups.keys() if 'data' not in proc]
    for proc in mc_procs:
        del bkg_groups[proc]
    bkg_groups['ttJets'] = ['ttJets_right', 'ttJets_matchable', 'ttJets_unmatchable', 'ttJets_sl_tau', 'ttJets_other']
    bkg_groups['singlet'] = ['singlet']
    bkg_groups['BKG'] = ['EWK', 'QCD']

    for hname in variables.keys():
        if hname not in hdict.keys():
            raise ValueError("%s not found in file" % hname)
        if 'DeepCSV_bDisc_' in hname: continue
        histo = hdict[hname].group("process", proc_bin, bkg_groups) # group by EWK+QCD into BKG, all ttJets together
        #set_trace()

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
            "overflow" : "under" if hname == 'DeepCSV_bDisc' else "none",
        }

        vlines = [(len(xrebinning)-1)*ybin for ybin in range(1, len(yrebinning)-1)] if histo.dense_dim() == 2 else None
        nbins = (len(xrebinning)-1)*(len(yrebinning)-1) if histo.dense_dim() == 2 else None
    
        if hname == 'Lep_iso':
            x_lims = (0., 0.15) if args.lepton == 'Muon' else (0., 0.1)
    

            # loop over different jet multiplicities
        for jmult in sorted(set([key[3] for key in histo.values().keys()])):
        #for jmult in ['3Jets']:
        #for jmult in ['4PJets']:
            #set_trace()
            print(jmult)

            if histo.dense_dim() == 1:
                cen_sb = histo[:, btag_reg_names_dict['Central']['reg'], 'nosys', jmult].integrate('btag').integrate('sys').integrate('jmult')
                up_sb = histo[:, btag_reg_names_dict['Up']['reg'], 'nosys', jmult].integrate('btag').integrate('sys').integrate('jmult')
                dw_sb = histo[:, btag_reg_names_dict['Down']['reg'], 'nosys', jmult].integrate('btag').integrate('sys').integrate('jmult')
                sig_histo = histo[:, btag_reg_names_dict['Signal']['reg'], :, jmult].integrate('btag').integrate('jmult')
            else:
                cen_sb = Plotter.linearize_hist(histo[:, btag_reg_names_dict['Central']['reg'], 'nosys', jmult].integrate('btag').integrate('sys').integrate('jmult'))
                up_sb = Plotter.linearize_hist(histo[:, btag_reg_names_dict['Up']['reg'], 'nosys', jmult].integrate('btag').integrate('sys').integrate('jmult'))
                dw_sb = Plotter.linearize_hist(histo[:, btag_reg_names_dict['Down']['reg'], 'nosys', jmult].integrate('btag').integrate('sys').integrate('jmult'))
                sig_histo = Plotter.linearize_hist(histo[:, btag_reg_names_dict['Signal']['reg'], :, jmult].integrate('btag').integrate('jmult'))

                ## plot systematic variations after bkg estimation is made
            #set_trace()
            for sys in systs_to_run:
                if sys not in histo.axis('sys')._sorted:
                    print('\n\n   Systematic %s not available, skipping\n\n' % sys)
                    continue
                print('BKG est:', jmult, sys, hname)
                #set_trace()

                for norm in ['SigMC']:
                    #set_trace()
                    bkg_est_histo = Plotter.BKG_Est(sig_reg=sig_histo, sb_reg=cen_sb, norm_type=norm, sys=sys, ignore_uncs=True)
                    up_bkg_est_histo = Plotter.BKG_Est(sig_reg=sig_histo, sb_reg=up_sb, norm_type=norm, sys=sys, ignore_uncs=True)
                    dw_bkg_est_histo = Plotter.BKG_Est(sig_reg=sig_histo, sb_reg=dw_sb, norm_type=norm, sys=sys, ignore_uncs=True)
                    #set_trace()

                    # create bkg hist shape variations with same yield
                        # bkg up variation yield (plus other MC contribution for plotting)
                    bkg_up_sumw = up_bkg_est_histo[Plotter.mc_samples].integrate('process').values()[()]
                        # bkg dw variation yield (plus other MC contribution for plotting)
                    bkg_dw_sumw = dw_bkg_est_histo[Plotter.mc_samples].integrate('process').values()[()]

                    if sys == 'nosys':                 
                        bkg_name = '%s%s_Norm' % (shape_reg, norm) if norm == 'Sideband' else '%s_Norm' % norm
                        bkg_dir = os.path.join(outdir, args.lepton, jmult, 'BKG_Est_orthog', bkg_name, sys_to_name[sys], 'VARY_SHAPE')
                        if not os.path.isdir(bkg_dir):
                            os.makedirs(bkg_dir)
    
                        fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                        fig.subplots_adjust(hspace=.07)
   
                        if hname == 'Jets_njets':
                            print('BKG est:', jmult, sys)
                            yields_txt, yields_json = Plotter.get_samples_yield_and_frac(bkg_est_histo, data_lumi_year['%ss' % args.lepton]/1000., sys=sys)
                            frac_name = '%s_yields_and_fracs_BKG_Est_orthog_%s' % ('_'.join([sys, jmult, args.lepton]), bkg_name)
                            plt_tools.print_table(yields_txt, filename='%s/%s.txt' % (bkg_dir, frac_name), print_output=True)
                            print('%s/%s.txt written' % (bkg_dir, frac_name))
                            with open('%s/%s.json' % (bkg_dir, frac_name), 'w') as out:
                                out.write(prettyjson.dumps(yields_json))

                            # plot bkg variations
                        Plotter.plot_1D(values=bkg_up_sumw, bins=bkg_est_histo.dense_axes()[0].edges(), xlabel=xtitle, xlimits=x_lims, ax=ax,
                            color=btag_reg_names_dict['Up']['color'], label=btag_reg_names_dict['Up']['label'])
                        Plotter.plot_1D(values=bkg_dw_sumw, bins=bkg_est_histo.dense_axes()[0].edges(), xlabel=xtitle, xlimits=x_lims, ax=ax,
                            color=btag_reg_names_dict['Down']['color'], label=btag_reg_names_dict['Down']['label'])

                        #set_trace() 
                        Plotter.plot_stack1d(ax, rax, bkg_est_histo, xlabel=xtitle, xlimits=x_lims, **mc_opts) if histo.dense_dim() == 1\
                            else Plotter.plot_stack1d(ax, rax, bkg_est_histo, xlabel='%s $\otimes$ %s' % (xtitle, ytitle), xlimits=(0, nbins), **mc_opts)

                            # plot bkg variations
                                # get ratio arrays
                        if (hname != 'mtt_vs_tlep_ctstar_abs'):
                            bkg_up_ratio_vals, bkg_up_ratio_bins = Plotter.get_ratio_arrays(num_vals=bkg_est_histo[Plotter.data_samples].integrate('process').values()[()], denom_vals=bkg_up_sumw, input_bins=bkg_est_histo.dense_axes()[0].edges())
                            bkg_dw_ratio_vals, bkg_dw_ratio_bins = Plotter.get_ratio_arrays(num_vals=bkg_est_histo[Plotter.data_samples].integrate('process').values()[()], denom_vals=bkg_dw_sumw, input_bins=bkg_est_histo.dense_axes()[0].edges())
                            rax.step(bkg_up_ratio_bins, bkg_up_ratio_vals, where='post', **{'linestyle' : '-', 'color' : btag_reg_names_dict['Up']['color']})
                            rax.step(bkg_dw_ratio_bins, bkg_dw_ratio_vals, where='post', **{'linestyle' : '-', 'color' : btag_reg_names_dict['Down']['color']})
    
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
                        figname = os.path.join(bkg_dir, '%s_BKG_Est_orthog_VARY_SHAPE_%s' % ('_'.join([sys, jmult, args.lepton, hname]), bkg_name))
                        fig.savefig(figname)
                        print('%s written' % figname)
                        plt.close()


if args.vary_shape_and_norm:
    systs_to_run = ['nosys']
    print('\nPlotting distributions for systematics:\n\t', *systs_to_run)

    ## combine EWK and QCD processes into BKG groups
    proc_bin = hist.Cat("process", "process", sorting='placement')
    bkg_groups = {proc:proc for proc in process_groups.keys()}
    mc_procs = [proc for proc in bkg_groups.keys() if 'data' not in proc]
    for proc in mc_procs:
        del bkg_groups[proc]
    bkg_groups['ttJets'] = ['ttJets_right', 'ttJets_matchable', 'ttJets_unmatchable', 'ttJets_sl_tau', 'ttJets_other']
    bkg_groups['singlet'] = ['singlet']
    bkg_groups['BKG'] = ['EWK', 'QCD']

    for hname in variables.keys():
        if hname not in hdict.keys():
            raise ValueError("%s not found in file" % hname)
        if 'DeepCSV_bDisc_' in hname: continue
        histo = hdict[hname].group("process", proc_bin, bkg_groups) # group by EWK+QCD into BKG, all ttJets together
        #set_trace()

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
            "overflow" : "under" if hname == 'DeepCSV_bDisc' else "none",
        }

        vlines = [(len(xrebinning)-1)*ybin for ybin in range(1, len(yrebinning)-1)] if histo.dense_dim() == 2 else None
        nbins = (len(xrebinning)-1)*(len(yrebinning)-1) if histo.dense_dim() == 2 else None
    
        if hname == 'Lep_iso':
            x_lims = (0., 0.15) if args.lepton == 'Muon' else (0., 0.1)
    

            # loop over different jet multiplicities
        for jmult in sorted(set([key[3] for key in histo.values().keys()])):
        #for jmult in ['3Jets']:
        #for jmult in ['4PJets']:
            #set_trace()
            print(jmult)

            if histo.dense_dim() == 1:
                cen_sb = histo[:, btag_reg_names_dict['Central']['reg'], 'nosys', jmult].integrate('btag').integrate('sys').integrate('jmult')
                up_sb = histo[:, btag_reg_names_dict['Up']['reg'], 'nosys', jmult].integrate('btag').integrate('sys').integrate('jmult')
                dw_sb = histo[:, btag_reg_names_dict['Down']['reg'], 'nosys', jmult].integrate('btag').integrate('sys').integrate('jmult')
                sig_histo = histo[:, btag_reg_names_dict['Signal']['reg'], :, jmult].integrate('btag').integrate('jmult')
            else:
                cen_sb = Plotter.linearize_hist(histo[:, btag_reg_names_dict['Central']['reg'], 'nosys', jmult].integrate('btag').integrate('sys').integrate('jmult'))
                up_sb = Plotter.linearize_hist(histo[:, btag_reg_names_dict['Up']['reg'], 'nosys', jmult].integrate('btag').integrate('sys').integrate('jmult'))
                dw_sb = Plotter.linearize_hist(histo[:, btag_reg_names_dict['Down']['reg'], 'nosys', jmult].integrate('btag').integrate('sys').integrate('jmult'))
                sig_histo = Plotter.linearize_hist(histo[:, btag_reg_names_dict['Signal']['reg'], :, jmult].integrate('btag').integrate('jmult'))

                ## plot systematic variations after bkg estimation is made
            #set_trace()
            for sys in systs_to_run:
                if sys not in histo.axis('sys')._sorted:
                    print('\n\n   Systematic %s not available, skipping\n\n' % sys)
                    continue
                print('BKG est:', jmult, sys, hname)
                #set_trace()

                for norm in ['SigMC']:
                    #set_trace()
                    bkg_est_histo = Plotter.BKG_Est(sig_reg=sig_histo, sb_reg=cen_sb, norm_type=norm, sys=sys, ignore_uncs=True)
                    up_bkg_est_histo = Plotter.BKG_Est(sig_reg=sig_histo, sb_reg=up_sb, norm_type=norm, sys=sys, ignore_uncs=True)
                    dw_bkg_est_histo = Plotter.BKG_Est(sig_reg=sig_histo, sb_reg=dw_sb, norm_type=norm, sys=sys, ignore_uncs=True)
                    #set_trace()

                    # create bkg hists varying both shape and norm
                    bkg_norm = np.sum(bkg_est_histo[Plotter.bkg_mask].integrate('process').values()[()])
                    st_tt_sumw = bkg_est_histo[Plotter.top_samples].integrate('process').values()[()]
                        # central sb
                    cen_bkg_shape, _ = Plotter.get_qcd_shape(bkg_est_histo[Plotter.bkg_mask].integrate('process'), overflow='none')
                            # +100% bkg yield (plus other MC contribution for plotting) for central sb region
                    cen_bkg_up_sumw = cen_bkg_shape*(bkg_norm*2) + st_tt_sumw
                            # -50% bkg yield (plus other MC contribution for plotting) for central sb region
                    cen_bkg_dw_sumw = cen_bkg_shape*(bkg_norm*0.5) + st_tt_sumw

                        # up sb
                    up_bkg_shape, _ = Plotter.get_qcd_shape(up_bkg_est_histo[Plotter.bkg_mask].integrate('process'), overflow='none')
                            # +100% bkg yield (plus other MC contribution for plotting) for up sb region
                    up_bkg_up_sumw = up_bkg_shape*(bkg_norm*2) + st_tt_sumw
                            # -50% bkg yield (plus other MC contribution for plotting) for up sb region
                    up_bkg_dw_sumw = up_bkg_shape*(bkg_norm*0.5) + st_tt_sumw

                        # down sb
                    dw_bkg_shape, _ = Plotter.get_qcd_shape(dw_bkg_est_histo[Plotter.bkg_mask].integrate('process'), overflow='none')
                            # +100% bkg yield (plus other MC contribution for plotting) for dw sb region
                    dw_bkg_up_sumw = dw_bkg_shape*(bkg_norm*2) + st_tt_sumw
                            # -50% bkg yield (plus other MC contribution for plotting) for dw sb region
                    dw_bkg_dw_sumw = dw_bkg_shape*(bkg_norm*0.5) + st_tt_sumw

                    if sys == 'nosys':                 
                        bkg_name = '%s%s_Norm' % (shape_reg, norm) if norm == 'Sideband' else '%s_Norm' % norm
                        bkg_dir = os.path.join(outdir, args.lepton, jmult, 'BKG_Est_orthog', bkg_name, sys_to_name[sys], 'VARY_SHAPE_AND_NORM')
                        if not os.path.isdir(bkg_dir):
                            os.makedirs(bkg_dir)
    
                        fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                        fig.subplots_adjust(hspace=.07)
   
                        if hname == 'Jets_njets':
                            print('BKG est:', jmult, sys)
                            yields_txt, yields_json = Plotter.get_samples_yield_and_frac(bkg_est_histo, data_lumi_year['%ss' % args.lepton]/1000., sys=sys)
                            frac_name = '%s_yields_and_fracs_BKG_Est_orthog_%s' % ('_'.join([sys, jmult, args.lepton]), bkg_name)
                            plt_tools.print_table(yields_txt, filename='%s/%s.txt' % (bkg_dir, frac_name), print_output=True)
                            print('%s/%s.txt written' % (bkg_dir, frac_name))
                            with open('%s/%s.json' % (bkg_dir, frac_name), 'w') as out:
                                out.write(prettyjson.dumps(yields_json))

                            # plot bkg variations
                                # central region
                        Plotter.plot_1D(values=cen_bkg_up_sumw, bins=bkg_est_histo.dense_axes()[0].edges(), xlabel=xtitle, xlimits=x_lims, ax=ax,
                            color='#e41a1c', label='%s +100%%' % btag_reg_names_dict['Central']['label'])
                        Plotter.plot_1D(values=cen_bkg_dw_sumw, bins=bkg_est_histo.dense_axes()[0].edges(), xlabel=xtitle, xlimits=x_lims, ax=ax,
                            color='#377eb8', label='%s -50%%' % btag_reg_names_dict['Central']['label'])
                                # up region
                        Plotter.plot_1D(values=up_bkg_up_sumw, bins=bkg_est_histo.dense_axes()[0].edges(), xlabel=xtitle, xlimits=x_lims, ax=ax,
                            color='#4daf4a', label='%s +100%%' % btag_reg_names_dict['Up']['label'])
                        Plotter.plot_1D(values=up_bkg_dw_sumw, bins=bkg_est_histo.dense_axes()[0].edges(), xlabel=xtitle, xlimits=x_lims, ax=ax,
                            color='#984ea3', label='%s -50%%' % btag_reg_names_dict['Up']['label'])
                                # down region
                        Plotter.plot_1D(values=dw_bkg_dw_sumw, bins=bkg_est_histo.dense_axes()[0].edges(), xlabel=xtitle, xlimits=x_lims, ax=ax,
                            color='#ff7f00', label='%s +100%%' % btag_reg_names_dict['Down']['label'])
                        Plotter.plot_1D(values=dw_bkg_dw_sumw, bins=bkg_est_histo.dense_axes()[0].edges(), xlabel=xtitle, xlimits=x_lims, ax=ax,
                            color='#ffff33', label='%s -50%%' % btag_reg_names_dict['Down']['label'])

                        #set_trace() 
                        Plotter.plot_stack1d(ax, rax, bkg_est_histo, xlabel=xtitle, xlimits=x_lims, **mc_opts) if histo.dense_dim() == 1\
                            else Plotter.plot_stack1d(ax, rax, bkg_est_histo, xlabel='%s $\otimes$ %s' % (xtitle, ytitle), xlimits=(0, nbins), **mc_opts)

                            # plot bkg variations
                                # get ratio arrays
                        if (hname != 'mtt_vs_tlep_ctstar_abs'):
                                # central sb region
                            cen_bkg_up_ratio_vals, cen_bkg_up_ratio_bins = Plotter.get_ratio_arrays(num_vals=bkg_est_histo[Plotter.data_samples].integrate('process').values()[()], denom_vals=cen_bkg_up_sumw, input_bins=bkg_est_histo.dense_axes()[0].edges())
                            cen_bkg_dw_ratio_vals, cen_bkg_dw_ratio_bins = Plotter.get_ratio_arrays(num_vals=bkg_est_histo[Plotter.data_samples].integrate('process').values()[()], denom_vals=cen_bkg_dw_sumw, input_bins=bkg_est_histo.dense_axes()[0].edges())
                            rax.step(cen_bkg_up_ratio_bins, cen_bkg_up_ratio_vals, where='post', **{'linestyle' : '-', 'color' : '#e41a1c'})
                            rax.step(cen_bkg_dw_ratio_bins, cen_bkg_dw_ratio_vals, where='post', **{'linestyle' : '-', 'color' : '#377eb8'})
                                # up sb region
                            up_bkg_up_ratio_vals, up_bkg_up_ratio_bins = Plotter.get_ratio_arrays(num_vals=bkg_est_histo[Plotter.data_samples].integrate('process').values()[()], denom_vals=up_bkg_up_sumw, input_bins=bkg_est_histo.dense_axes()[0].edges())
                            up_bkg_dw_ratio_vals, up_bkg_dw_ratio_bins = Plotter.get_ratio_arrays(num_vals=bkg_est_histo[Plotter.data_samples].integrate('process').values()[()], denom_vals=up_bkg_dw_sumw, input_bins=bkg_est_histo.dense_axes()[0].edges())
                            rax.step(up_bkg_up_ratio_bins, up_bkg_up_ratio_vals, where='post', **{'linestyle' : '-', 'color' : '#4daf4a'})
                            rax.step(up_bkg_dw_ratio_bins, up_bkg_dw_ratio_vals, where='post', **{'linestyle' : '-', 'color' : '#984ea3'})
                                # down sb region
                            dw_bkg_up_ratio_vals, dw_bkg_up_ratio_bins = Plotter.get_ratio_arrays(num_vals=bkg_est_histo[Plotter.data_samples].integrate('process').values()[()], denom_vals=dw_bkg_up_sumw, input_bins=bkg_est_histo.dense_axes()[0].edges())
                            dw_bkg_dw_ratio_vals, dw_bkg_dw_ratio_bins = Plotter.get_ratio_arrays(num_vals=bkg_est_histo[Plotter.data_samples].integrate('process').values()[()], denom_vals=dw_bkg_dw_sumw, input_bins=bkg_est_histo.dense_axes()[0].edges())
                            rax.step(dw_bkg_up_ratio_bins, dw_bkg_up_ratio_vals, where='post', **{'linestyle' : '-', 'color' : '#ff7f00'})
                            rax.step(dw_bkg_dw_ratio_bins, dw_bkg_dw_ratio_vals, where='post', **{'linestyle' : '-', 'color' : '#ffff33'})
    
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
                        figname = os.path.join(bkg_dir, '%s_BKG_Est_orthog_VARY_SHAPE_AND_NORM_%s' % ('_'.join([sys, jmult, args.lepton, hname]), bkg_name))
                        fig.savefig(figname)
                        print('%s written' % figname)
                        plt.close()


#    if args.save_sys:
#        from coffea.util import save
#        outname = '%s/%s/%s_Post_QCD_Est_dict.coffea' % (outdir, args.lepton, args.lepton)
#        save(post_bkg_est, outname)
#        #outname = '%s/%s/%s_QCD_Est_mtt_ctstar_dict.coffea' % (outdir, args.lepton, args.lepton)
#        #save(mtt_ctstar_bkg_est, outname)
#        print('%s written' % outname)

