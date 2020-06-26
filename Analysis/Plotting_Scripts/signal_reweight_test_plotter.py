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
parser.add_argument('--qcd_est', action='store_true', help='Estimate qcd contribution')

args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'signal_reweight_test'

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
    'mthad' : ('m($t_{h}$) [GeV]', 2, (0., 300.), False),
    ##'pt_thad' : ('$p_{T}$($t_{h}$) [GeV]', 2, (0., 500.), True),
    ##'pt_tlep' : ('$p_{T}$($t_{l}$) [GeV]', 2, (0., 500.), True),
    ##'pt_tt' : ('$p_{T}$($t\\bar{t}$) [GeV]', 2, (0., 500.), True),
    ##'eta_thad' : ('$\\eta$($t_{h}$)', 2, (-4., 4.), True),
    ##'eta_tlep' : ('$\\eta$($t_{l}$)', 2, (-4., 4.), True),
    ##'eta_tt' : ('$\\eta$($t\\bar{t}$)', 2, (-4., 4.), True),
    'tlep_ctstar' : ('cos($\\theta^{*}_{t_{l}}$)', 2, (-1., 1.), False),
    ##'tlep_ctstar_abs' : ('|cos($\\theta^{*}_{t_{l}}$)|', 2, (0., 1.), True),
    ##'full_disc' : ('$\\lambda_{C}$', 2, (5, 25.), True),
    ##'mass_disc' : ('$\\lambda_{M}$', 2, (0, 20.), True),
    ##'ns_disc' : ('$\\lambda_{NS}$', 2, (3., 10.), True),
    ##'Jets_pt' : ('$p_{T}$(jets) [GeV]', 2, (0., 300.), True),
    ##'Jets_eta' : ('$\\eta$(jets)', 2, (-2.6, 2.6), True),
    ##'Jets_njets' : ('$n_{jets}$', 1, (0, 15), True),
    ##'Jets_LeadJet_pt' : ('$p_{T}$(leading jet) [GeV]', 2, (0., 300.), True),
    ##'Jets_LeadJet_eta' : ('$\\eta$(leading jet)', 2, (-2.6, 2.6), True),
    ##'Lep_pt' : ('$p_{T}$(%s) [GeV]' % objtypes['Lep'][args.lepton], 2, (0., 300.), True),
    ##'Lep_eta' : ('$\\eta$(%s)' % objtypes['Lep'][args.lepton], 2, (-2.6, 2.6), True),
    ##'Lep_iso' : ('pfRelIso, %s' % objtypes['Lep'][args.lepton], 1, (0., 1.), True),
    ##'MT' : ('$M_{T}$ [GeV]', 1, (0., 300.), True),
}

widthTOname = lambda width : str(width).replace('.', 'p')
nameTOwidth = lambda width : str(width).replace('p', '.')

def get_title(boson, mass, width, sampletype):
    mass_title = 'm($%s$)=%s GeV' % (boson[0], mass.split('M')[-1])
    width_title = '$\\Gamma$/m($%s$)=%s%%' % (boson[0], nameTOwidth(width.split('W')[-1]))
    samp_title = sampletype[0]

    return '%s, %s, %s' % (mass_title, width_title, samp_title)



    ## get plotting colors/settings
stack_fill_opts = {'alpha': 0.8, 'edgecolor':(0,0,0,.5)}
stack_error_opts = {'edgecolor':(0,0,0,.5)}

    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open('%s/inputs/lumis_data.json' % proj_dir).read())[args.year]
lumi_correction = load('%s/Corrections/%s/MC_LumiWeights_IgnoreSigEvts.coffea' % (proj_dir, jobid))[args.year]['%ss' % args.lepton]

for hname in hdict.keys():
    if hname == 'cutflow': continue
    hdict[hname].scale(lumi_correction, axis='dataset')


    ## make plots
for hname in variables.keys():
    if hname not in hdict.keys():
        raise ValueError("%s not found in file" % hname)
    #set_trace()
    histo = hdict[hname][:, :, args.lepton, :, :].integrate('leptype') # dataset, jmult, leptype, btag, lepcat

    if histo.dense_dim() == 1:
        xtitle, rebinning, x_lims, withData = variables[hname]
        if rebinning != 1:
            xaxis_name = histo.dense_axes()[0].name
            histo = histo.rebin(xaxis_name, rebinning)

        #set_trace()
        ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
        for jmult in list(set([key[1] for key in histo.values().keys()])):
            for lepcat in list(set([key[3] for key in histo.values().keys()])):
                for btagregion in list(set([key[2] for key in histo.values().keys()])):
                    for dataset in list(set([key[0] for key in histo.values().keys()])):
                        boson, mass, width, sampletype = dataset.split('_')
                        label = get_title(boson, mass, width, sampletype)

                        pltdir = '/'.join([outdir, args.lepton, jmult, lepcat, btagregion, dataset])
                        if not os.path.isdir(pltdir):
                            os.makedirs(pltdir)

                        fig, ax = plt.subplots()
                        fig.subplots_adjust(hspace=.07)

                        #hslice = histo[:, jmult, btagregion, lepcat].integrate('jmult').integrate('lepcat').integrate('btag')
                        hslice = histo[dataset, jmult, btagregion, lepcat].integrate('jmult').integrate('lepcat').integrate('btag').integrate('dataset')

                        ax = Plotter.plot_1D(hslice.values()[()], hslice.dense_axes()[0].edges(), xlimits=x_lims, xlabel=xtitle, ax=ax, histtype='step')
                        #set_trace()
                        upper_ylim = ax.get_ylim()[1]*3 if ax.get_ylim()[1] > 0 else ax.get_ylim()[1]-ax.get_ylim()[1]*0.5
                        ax.set_ylim(ax.get_ylim()[0], upper_ylim)
                        #ax = Plotter.plot_1D(hslice.values()[()], hslice.dense_axes()[0].edges(), xlimits=x_lims, xlabel=xtitle, ax=ax, label=label, histtype='step')
                        #ax.legend(loc='upper right')

                            # add lepton/jet multiplicity label
                        #set_trace()
                        ax.text(0.95, 0.95, label, fontsize=rcParams['font.size']*0.75, 
                            horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes
                        )
                        ax.text(
                            0.02, 0.90, "%s, %s\n%s" % (lep_cats[lepcat], jet_mults[jmult], btag_cats[btagregion]),
                            fontsize=rcParams['font.size']*0.75, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
                        )
                        ax = hep.cms.cmslabel(ax=ax, data=withData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))

                        #set_trace()
                        figname = '%s/%s' % (pltdir, '_'.join([dataset, jmult, args.lepton, lepcat, btagregion, hname]))
                        fig.savefig(figname)
                        print('%s written' % figname)
                        plt.close()


