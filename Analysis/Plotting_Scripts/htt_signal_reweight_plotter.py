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
#import Utilities.systematics as systematics

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('year', choices=['2016', '2017', '2018'], help='What year is the ntuple from.')
parser.add_argument('lepton', choices=['Electron', 'Muon'], help='Choose which lepton to make plots for')
parser.add_argument('--plot_indiv', action='store_true', help='Make plots for each mass and width point separately')
parser.add_argument('--comp_widths', action='store_true', help='Compare widths for single mass point')
parser.add_argument('--comp_masses', action='store_true', help='Compare masses for single width point')

args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'htt_signal_reweight'

input_dir = '/'.join([proj_dir, 'results', '%s_%s' % (args.year, jobid), analyzer])
f_ext = 'TOT.coffea'
outdir = '/'.join([proj_dir, 'plots', '%s_%s' % (args.year, jobid), analyzer])
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

variables = {
    'mtt' : ('m($t\\bar{t}$) [GeV]', 2, (200., 2000.), False),
    'tlep_ctstar' : ('cos($\\theta^{*}_{t_{l}}$)', 2, (-1., 1.), False),
    'tlep_ctstar_abs' : ('|cos($\\theta^{*}_{t_{l}}$)|', 2, (0., 1.), True),
    'mthad' : ('m($t_{h}$) [GeV]', 2, (0., 300.), False),
    'pt_thad' : ('$p_{T}$($t_{h}$) [GeV]', 2, (0., 500.), True),
    'pt_tlep' : ('$p_{T}$($t_{l}$) [GeV]', 2, (0., 500.), True),
    'pt_tt' : ('$p_{T}$($t\\bar{t}$) [GeV]', 2, (0., 500.), True),
    'eta_thad' : ('$\\eta$($t_{h}$)', 2, (-4., 4.), True),
    'eta_tlep' : ('$\\eta$($t_{l}$)', 2, (-4., 4.), True),
    'eta_tt' : ('$\\eta$($t\\bar{t}$)', 2, (-4., 4.), True),
    'full_disc' : ('$\\lambda_{C}$', 2, (5, 25.), True),
    'mass_disc' : ('$\\lambda_{M}$', 2, (0, 20.), True),
    'ns_disc' : ('$\\lambda_{NS}$', 2, (3., 10.), True),
    'Jets_pt' : ('$p_{T}$(jets) [GeV]', 2, (0., 300.), True),
    'Jets_eta' : ('$\\eta$(jets)', 2, (-2.6, 2.6), True),
    'Jets_njets' : ('$n_{jets}$', 1, (0, 15), True),
    'Jets_LeadJet_pt' : ('$p_{T}$(leading jet) [GeV]', 2, (0., 300.), True),
    'Jets_LeadJet_eta' : ('$\\eta$(leading jet)', 2, (-2.6, 2.6), True),
    'Lep_pt' : ('$p_{T}$(%s) [GeV]' % objtypes['Lep'][args.lepton], 2, (0., 300.), True),
    'Lep_eta' : ('$\\eta$(%s)' % objtypes['Lep'][args.lepton], 2, (-2.6, 2.6), True),
    'Lep_iso' : ('pfRelIso, %s' % objtypes['Lep'][args.lepton], 1, (0., 1.), True),
    'MT' : ('$M_{T}$ [GeV]', 1, (0., 300.), True),
    'MET_pt' : ('$p_{T}$(MET) [GeV]', 1, (0., 300.), True),
    'MET_phi' : ('$\phi$(MET)', 1, (-3.2, 3.2), True),
}

widthTOname = lambda width : str(width).replace('.', 'p')
nameTOwidth = lambda width : str(width).replace('p', '.')

def get_title(boson, mass, width, shape):
    mass_title = 'm($%s$)=%s GeV' % (boson[0], mass.split('M')[-1])
    width_title = '$\\Gamma$/m($%s$)=%s%%' % (boson[0], nameTOwidth(width.split('W')[-1]))
    samp_title = shape[0]

    return '%s, %s, %s' % (mass_title, width_title, samp_title)



    ## get plotting colors/settings
stack_fill_opts = {'alpha': 0.8, 'edgecolor':(0,0,0,.5)}
stack_error_opts = {'edgecolor':(0,0,0,.5)}

    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open('%s/inputs/lumis_data.json' % proj_dir).read())[args.year]
lumi_correction = load('%s/Corrections/%s/MC_LumiWeights_IgnoreSigEvts.coffea' % (proj_dir, jobid))[args.year]['%ss' % args.lepton]

for idx, hname in enumerate(hdict.keys()):
    if 'cutflow' in hname: continue
    hdict[hname].scale(lumi_correction, axis='dataset')
    if idx == 0:
        signals = sorted(set(['_'.join(key[0].split('_')[:-1]) for key in hdict[hname].values().keys()]))
        bosons = sorted(set([key[0].split('_')[0] for key in hdict[hname].values().keys()]))
        masses = sorted(set([key[0].split('_')[1] for key in hdict[hname].values().keys()]))
        widths = sorted(set([key[0].split('_')[2].strip('W') for key in hdict[hname].values().keys()]), key=lambda width: float(nameTOwidth(width)))
        shapes = sorted(set([key[0].split('_')[3] for key in hdict[hname].values().keys()]))

    ## make plots
for hname in variables.keys():
    if hname not in hdict.keys():
        raise ValueError("%s not found in file" % hname)
    #set_trace()
    histo = hdict[hname][:, :, :, args.lepton, 'btagPass', 'Tight'].integrate('leptype').integrate('btag').integrate('lepcat') # dataset, jmult, leptype, btag, lepcat

    #set_trace()
    systs = sorted(set([key[1] for key in histo.values().keys()]))
    systypes = sorted(set(['_'.join(sys.split('_')[:-1]) for sys in systs if (not sys == 'nosys') and (not ('UP' and 'DW') in sys)]))
    if ('RENORM_UP_FACTOR_DW' and 'RENORM_DW_FACTOR_UP') in systs:
        systypes += ['RENFACTOR_DIFF']

    if histo.dense_dim() == 1:
        xtitle, rebinning, x_lims, withData = variables[hname]
        if rebinning != 1:
            xaxis_name = histo.dense_axes()[0].name
            histo = histo.rebin(xaxis_name, rebinning)

        if hname == 'Lep_iso':
            x_lims = (0., 0.15) if args.lepton == 'Muon' else (0., 0.1)
        #set_trace()
        ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
        for jmult in sorted(set([key[2] for key in histo.values().keys()])):
            if (hname == 'mtt') or (hname == 'tlep_ctstar') or (hname == 'tlep_ctstar_abs'):
                if args.comp_widths:
                    #set_trace()
                    for boson in bosons:
                        pltdir = '/'.join([outdir, args.lepton, jmult, lepcat, btagregion, 'Comp_Widths', boson])
                        if not os.path.isdir(pltdir):
                            os.makedirs(pltdir)
                        for mass in masses:
                            mass_title = 'm($%s$)=%s GeV' % (boson[0], mass.split('M')[-1])
                            for shape in shapes:
                                #set_trace()
                                    ## sort samples by increasing width (2.5, 5, 10, 25)
                                samples = sorted([sig for sig in signals if fnmatch.fnmatch(sig, '%s*%s*%s' % (boson, mass, shape))], key=lambda sample: float(nameTOwidth(sample.split('_')[2])[1:]))

                                if shape == 'Int':
                                    #set_trace()
                                    pos_fig, pos_ax = plt.subplots()
                                    pos_fig.subplots_adjust(hspace=.07)
                                    neg_fig, neg_ax = plt.subplots()
                                    neg_fig.subplots_adjust(hspace=.07)

                                    for sig in samples:
                                        boson, mass, width, shape = sig.split('_')
                                        opts = Plotter.styles.styles[width]

                                        pos_histo = histo['%s_pos' % sig, jmult, btagregion, lepcat].integrate('jmult').integrate('lepcat').integrate('btag').integrate('dataset')
                                        Plotter.plot_1D(pos_histo.values()[()], pos_histo.dense_axes()[0].edges(), ax=pos_ax, xlimits=x_lims, xlabel=xtitle, color=opts['color'], label=opts['name'], histtype='step')
                                        neg_histo = histo['%s_neg' % sig, jmult, btagregion, lepcat].integrate('jmult').integrate('lepcat').integrate('btag').integrate('dataset')
                                        Plotter.plot_1D(neg_histo.values()[()], neg_histo.dense_axes()[0].edges(), ax=neg_ax, xlimits=x_lims, xlabel=xtitle, color=opts['color'], label=opts['name'], histtype='step')

                                    pos_ax.legend(loc='upper right', title='%s, Int, w > 0' % mass_title)
                                        # add labels
                                    pos_ax.text(
                                        0.02, 0.90, "%s, %s\n%s" % (lep_cats[lepcat], jet_mults[jmult], btag_cats[btagregion]),
                                        fontsize=rcParams['font.size']*0.75, horizontalalignment='left', verticalalignment='bottom', transform=pos_ax.transAxes
                                    )
                                    pos_ax = hep.cms.cmslabel(ax=pos_ax, data=False, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
                                    #set_trace()
                                    pos_figname = '%s/%s' % (pltdir, '_'.join([boson, mass, 'Int', 'pos', 'WidthComp', jmult, args.lepton, lepcat, btagregion, hname]))
                                    pos_fig.savefig(pos_figname)
                                    print('%s written' % pos_figname)

                                    neg_ax.legend(loc='upper right', title='%s, Int, w < 0' % mass_title)
                                        # add labels
                                    neg_ax.text(
                                        0.02, 0.90, "%s, %s\n%s" % (lep_cats[lepcat], jet_mults[jmult], btag_cats[btagregion]),
                                        fontsize=rcParams['font.size']*0.75, horizontalalignment='left', verticalalignment='bottom', transform=neg_ax.transAxes
                                    )
                                    neg_ax = hep.cms.cmslabel(ax=neg_ax, data=False, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
                                    #set_trace()
                                    neg_figname = '%s/%s' % (pltdir, '_'.join([boson, mass, 'Int', 'neg', 'WidthComp', jmult, args.lepton, lepcat, btagregion, hname]))
                                    neg_fig.savefig(neg_figname)
                                    print('%s written' % neg_figname)

                                    plt.close()
                                else:
                                    fig, ax = plt.subplots()
                                    fig.subplots_adjust(hspace=.07)

                                    for sig in samples:
                                        boson, mass, width, shape = sig.split('_')
                                        opts = Plotter.styles.styles[width]

                                        w_histo = histo['%s_pos' % sig, jmult, btagregion, lepcat].integrate('jmult').integrate('lepcat').integrate('btag').integrate('dataset')
                                        Plotter.plot_1D(w_histo.values()[()], w_histo.dense_axes()[0].edges(), xlimits=x_lims, xlabel=xtitle, color=opts['color'], label=opts['name'], histtype='step')
                                    #set_trace()
                                    ax.legend(loc='upper right', title='%s, Res' % mass_title)

                                        # add labels
                                    ax.text(
                                        0.02, 0.90, "%s, %s\n%s" % (lep_cats[lepcat], jet_mults[jmult], btag_cats[btagregion]),
                                        fontsize=rcParams['font.size']*0.75, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
                                    )
                                    ax = hep.cms.cmslabel(ax=ax, data=False, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))

                                    figname = '%s/%s' % (pltdir, '_'.join([boson, mass, 'Res', 'WidthComp', jmult, args.lepton, lepcat, btagregion, hname]))
                                    fig.savefig(figname)
                                    print('%s written' % figname)
                                    plt.close()

                if args.comp_masses:
                    #set_trace()
                    for boson in bosons:
                        pltdir = '/'.join([outdir, args.lepton, jmult, lepcat, btagregion, 'Comp_Masses', boson])
                        if not os.path.isdir(pltdir):
                            os.makedirs(pltdir)
                        for width in widths:
                            width_title = '$\\Gamma$/m($%s$)=%s%%' % (boson[0], nameTOwidth(width.split('W')[-1]))
                            for shape in shapes:
                                #set_trace()
                                    ## sort samples by increasing mass (400, 500, 600, 750)
                                samples = sorted([sig for sig in signals if fnmatch.fnmatch(sig, '%s*%s*%s' % (boson, width, shape))])

                                if shape == 'Int':
                                    #set_trace()
                                    pos_fig, pos_ax = plt.subplots()
                                    pos_fig.subplots_adjust(hspace=.07)
                                    neg_fig, neg_ax = plt.subplots()
                                    neg_fig.subplots_adjust(hspace=.07)

                                    for sig in samples:
                                        boson, mass, width, shape = sig.split('_')
                                        opts = Plotter.styles.styles[mass]

                                        pos_histo = histo['%s_pos' % sig, jmult, btagregion, lepcat].integrate('jmult').integrate('lepcat').integrate('btag').integrate('dataset')
                                        Plotter.plot_1D(pos_histo.values()[()], pos_histo.dense_axes()[0].edges(), ax=pos_ax, xlimits=x_lims, xlabel=xtitle, color=opts['color'], label=opts['name'], histtype='step')
                                        neg_histo = histo['%s_neg' % sig, jmult, btagregion, lepcat].integrate('jmult').integrate('lepcat').integrate('btag').integrate('dataset')
                                        Plotter.plot_1D(neg_histo.values()[()], neg_histo.dense_axes()[0].edges(), ax=neg_ax, xlimits=x_lims, xlabel=xtitle, color=opts['color'], label=opts['name'], histtype='step')

                                    pos_ax.legend(loc='upper right', title='%s, Int, w > 0' % width_title)
                                        # add labels
                                    pos_ax.text(
                                        0.02, 0.90, "%s, %s\n%s" % (lep_cats[lepcat], jet_mults[jmult], btag_cats[btagregion]),
                                        fontsize=rcParams['font.size']*0.75, horizontalalignment='left', verticalalignment='bottom', transform=pos_ax.transAxes
                                    )
                                    pos_ax = hep.cms.cmslabel(ax=pos_ax, data=False, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
                                    #set_trace()
                                    pos_figname = '%s/%s' % (pltdir, '_'.join([boson, width, 'Int', 'pos', 'MassComp', jmult, args.lepton, lepcat, btagregion, hname]))
                                    pos_fig.savefig(pos_figname)
                                    print('%s written' % pos_figname)

                                    neg_ax.legend(loc='upper right', title='%s, Int, w < 0' % width_title)
                                        # add labels
                                    neg_ax.text(
                                        0.02, 0.90, "%s, %s\n%s" % (lep_cats[lepcat], jet_mults[jmult], btag_cats[btagregion]),
                                        fontsize=rcParams['font.size']*0.75, horizontalalignment='left', verticalalignment='bottom', transform=neg_ax.transAxes
                                    )
                                    neg_ax = hep.cms.cmslabel(ax=neg_ax, data=False, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
                                    #set_trace()
                                    neg_figname = '%s/%s' % (pltdir, '_'.join([boson, width, 'Int', 'neg', 'MassComp', jmult, args.lepton, lepcat, btagregion, hname]))
                                    neg_fig.savefig(neg_figname)
                                    print('%s written' % neg_figname)

                                    plt.close()
                                else:
                                    fig, ax = plt.subplots()
                                    fig.subplots_adjust(hspace=.07)

                                    for sig in samples:
                                        boson, mass, width, shape = sig.split('_')
                                        opts = Plotter.styles.styles[mass]

                                        w_histo = histo['%s_pos' % sig, jmult, btagregion, lepcat].integrate('jmult').integrate('lepcat').integrate('btag').integrate('dataset')
                                        Plotter.plot_1D(w_histo.values()[()], w_histo.dense_axes()[0].edges(), xlimits=x_lims, xlabel=xtitle, color=opts['color'], label=opts['name'], histtype='step')
                                    #set_trace()
                                    ax.legend(loc='upper right', title='%s, Res' % width_title)

                                        # add labels
                                    ax.text(
                                        0.02, 0.90, "%s, %s\n%s" % (lep_cats[lepcat], jet_mults[jmult], btag_cats[btagregion]),
                                        fontsize=rcParams['font.size']*0.75, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
                                    )
                                    ax = hep.cms.cmslabel(ax=ax, data=False, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))

                                    figname = '%s/%s' % (pltdir, '_'.join([boson, width, 'Res', 'MassComp', jmult, args.lepton, lepcat, btagregion, hname]))
                                    fig.savefig(figname)
                                    print('%s written' % figname)
                                    plt.close()

            if args.plot_indiv:
                #set_trace()
                for sig in signals:
                    boson, mass, width, shape = sig.split('_')
                    label = get_title(boson, mass, width, shape)

                    #set_trace()
                    for sys in systs:
                    #for sys in ['nosys']:
                        print(jmult, sig, sys, hname)
                        pltdir = '/'.join([outdir, args.lepton, jmult, 'Individual', sig, sys])
                        if not os.path.isdir(pltdir):
                            os.makedirs(pltdir)

                        fig, ax = plt.subplots()
                        fig.subplots_adjust(hspace=.07)

                        if shape == 'Int':
                            pos_hist = histo['%s_pos' % sig, sys, jmult].integrate('jmult').integrate('sys').integrate('dataset')
                            neg_hist = histo['%s_neg' % sig, sys, jmult].integrate('jmult').integrate('sys').integrate('dataset')
                            int_hist = pos_hist.copy()
                            int_hist.add(neg_hist)
                            #set_trace()

                            Plotter.plot_1D(pos_hist.values()[()], pos_hist.dense_axes()[0].edges(), xlimits=x_lims, xlabel=xtitle, color='r', label='w > 0', histtype='step')
                            Plotter.plot_1D(neg_hist.values()[()], neg_hist.dense_axes()[0].edges(), xlimits=x_lims, xlabel=xtitle, color='b', label='w < 0', histtype='step')
                            Plotter.plot_1D(int_hist.values()[()], int_hist.dense_axes()[0].edges(), xlimits=x_lims, xlabel=xtitle, color='k', label='Sum', histtype='step')
                            ymin, ymax = min(min(pos_hist.values()[()]), min(neg_hist.values()[()]), min(int_hist.values()[()])), max(max(pos_hist.values()[()]), max(neg_hist.values()[()]), max(int_hist.values()[()]))
                            ax.set_ylim(ymin*1.05, ymax*1.2)

                        else:
                            hslice = histo['%s_pos' % sig, sys, jmult].integrate('jmult').integrate('sys').integrate('dataset')

                            Plotter.plot_1D(hslice.values()[()], hslice.dense_axes()[0].edges(), xlimits=x_lims, xlabel=xtitle, label='w > 0', histtype='step')

                        ax.legend(loc='upper right', title=label)

                            # add labels
                        ax.text(
                            0.02, 0.94, "%s, %s" % (objtypes['Lep'][args.lepton], jet_mults[jmult]),
                            fontsize=rcParams['font.size']*0.75, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
                        )
                        ax = hep.cms.cmslabel(ax=ax, data=False, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))

                        #set_trace()
                        figname = '%s/%s' % (pltdir, '_'.join([sig, jmult, args.lepton, sys, hname]))
                        fig.savefig(figname)
                        print('%s written' % figname)
                        plt.close()
