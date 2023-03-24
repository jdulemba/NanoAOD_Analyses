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

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('year', choices=['2016', '2017', '2018'], help='What year is the ntuple from.')
parser.add_argument('lepton', choices=['Electron', 'Muon'], help='Choose which lepton to make plots for')
args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
base_jobid = os.environ['base_jobid']
analyzer = 'evtWeights'

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

variables = {
    'GenWeights' : ('Generator Weight', 1, (-500., 500.)),
    'TotWeights' : ('Total Event Weight', 1, (-500., 500.)),
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


    ## make plots
for hname in variables.keys():
    if hname not in hdict.keys():
        raise ValueError("%s not found in file" % hname)
    histo = hdict[hname][:, 'nosys', :, args.lepton, :, :].integrate('sys').integrate('leptype') # process, sys, jmult, btag, lepcat
    #set_trace()

    xtitle, rebinning, x_lims = variables[hname]
    xaxis_name = histo.dense_axes()[0].name
    if rebinning != 1:
        histo = histo.rebin(xaxis_name, rebinning)

    ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
    for jmult in sorted(set([key[1] for key in histo.values().keys()])):
        for lepcat in sorted(set([key[3] for key in histo.values().keys()])):
            for btagregion in sorted(set([key[2] for key in histo.values().keys()])):
                pltdir = os.path.join(outdir, args.lepton, jmult, lepcat, btagregion)
                if not os.path.isdir(pltdir):
                    os.makedirs(pltdir)

                for sample in sorted(set([key[0] for key in histo.values().keys()])):
                    print(', '.join([jmult, lepcat, btagregion, hname, sample])) 
                    fig, ax = plt.subplots()
                    fig.subplots_adjust(hspace=.07)

                    hslice = histo[sample, jmult, btagregion, lepcat].integrate('jmult').integrate('lepcat').integrate('btag').integrate('dataset')

                    if hslice.values():
                        Plotter.plot_1D(*hslice.values().values(), histo.axis(xaxis_name).edges(), xlimits=(-5., 5.) if sample.startswith('QCD') else x_lims, xlabel=xtitle,
                            ylabel='Events (Unweighted)', ax=ax, label='%s\n%s' % (sample, format(lumi_correction[sample], '.3f')))
                        ax.legend(loc='upper right')

                            # add lepton/jet multiplicity label
                        ax.text(
                            0.02, 0.88, "%s, %s\n%s" % (lep_cats[lepcat], jet_mults[jmult], btag_cats[btagregion]),
                            fontsize=rcParams['font.size']*0.75, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
                        )
                        hep.cms.label(ax=ax, rlabel=args.year)

                        #set_trace()
                        figname = os.path.join(pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, hname, sample]))
                        fig.savefig(figname)
                        print('%s written' % figname)
                        plt.close()
