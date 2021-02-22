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
#from coffea import hist
#import numpy as np
import Utilities.Plotter as Plotter

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
base_jobid = os.environ['base_jobid']
analyzer = 'genpartons'

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('year', choices=['2017', '2018'] if base_jobid == 'ULnanoAOD' else ['2016', '2017', '2018'], help='What year is the ntuple from.')
args = parser.parse_args()

input_dir = os.path.join(proj_dir, 'results', '%s_%s' % (args.year, jobid), analyzer)
f_ext = 'TOT.coffea'
outdir = os.path.join(proj_dir, 'plots', '%s_%s' % (args.year, jobid), analyzer)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

fnames = sorted(['%s/%s' % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)])
hdict = plt_tools.add_coffea_files(fnames) if len(fnames) > 1 else load(fnames[0])

objects = {
    'SL' : {
        'TTbar' : ('$t\\bar{t}$', (200., 2000.)),
        'THad' : ('$t_{h}$', (150., 200.)),
        'TLep' : ('$t_{l}$', (150., 200.)),
        'BHad' : ('$b_{h}$', (0., 5.)),
        'BLep' : ('$b_{l}$', (0., 5.)),
        'WHad' : ('$W_{h}$', (50., 100.)),
        'WLep' : ('$W_{l}$', (50., 100.)),
        'Lepton' : ('l', (0., 5.)),
    },
    'DL' : {
        'TTbar' : ('$t\\bar{t}$', (200., 2000.)),
        'Top' : ('$t$', (150., 200.)),
        'Tbar' : ('$\\bar{t}$', (150., 200.)),
        'B' : ('$b$', (0., 5.)),
        'Bbar' : ('$\\bar{b}$', (0., 5.)),
        'Wplus' : ('$W^{+}$', (50., 100.)),
        'Wminus' : ('$W^{-}$', (50., 100.)),
        'First_plus' : ('$l^{+}$', (0., 5.)),
        'Second_plus' : ('$\\nu_{l}$', (0., 5.)),
        'First_minus' : ('$l^{-}$', (0., 5.)),
        'Second_minus' : ('$\\bar{\\nu_{l}}$', (0., 5.)),
    },
    'Had' : {
        'TTbar' : ('$t\\bar{t}$', (200., 2000.)),
        'Top' : ('$t$', (150., 200.)),
        'Tbar' : ('$\\bar{t}$', (150., 200.)),
        'B' : ('$b$', (0., 5.)),
        'Bbar' : ('$\\bar{b}$', (0., 5.)),
        'Wplus' : ('$W^{+}$', (50., 100.)),
        'Wminus' : ('$W^{-}$', (50., 100.)),
        'First_plus' : ('leading +charged partons from W', (0., 5.)),
        'Second_plus' : ('subleading +charged partons from W', (0., 5.)),
        'First_minus' : ('leading -charged partons from W', (0., 5.)),
        'Second_minus' : ('subleading -charged partons from W', (0., 5.)),
        'Up_plus' : ('+charged up-type partons from W', (0., 5.)),
        'Down_plus' : ('+charged down-type partons from W', (0., 5.)),
        'Up_minus' : ('-charged up-type partons from W', (0., 5.)),
        'Down_minus' : ('-charged down-type partons from W', (0., 5.)),
    }
}

ttdecay_types = {
    'SL' : '$t\\bar{t} \\rightarrow$ lj',
    'DL' : '$t\\bar{t} \\rightarrow$ ll',
    'Had': '$t\\bar{t} \\rightarrow$ jj',
}

variables = {
    'pt' : ('$p_{T}$(obj) [GeV]', 2, (0., 500.)),
    'eta': ('$\\eta$(obj)', 2, (-2.6, 2.6)),
    'phi': ('$\\phi$(obj)', 2, (-4, 4)),
    'mass': ('m(obj) [GeV]', 1, (0., 300.)),
    'energy': ('E(obj) [GeV]', 2, (0., 1000.)),
}


    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, 'inputs', '%s_lumis_data.json' % base_jobid)).read())[args.year]
        ## temporary
if base_jobid == 'ULnanoAOD':
    lumi_to_use = data_lumi_year['Muons']/1000. if args.year == '2018' else data_lumi_year['Electrons']/1000.
else:
    lumi_to_use = (data_lumi_year['Muons']+data_lumi_year['Electrons'])/2000.
lumi_correction = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'MC_LumiWeights.coffea'))[args.year]

        # scale events by lumi correction
for hname in hdict.keys():
    if hname == 'cutflow': continue
    hdict[hname].scale(lumi_correction['Muons'], axis='dataset')
    hdict[hname] = hdict[hname].integrate('dataset')


    ## make bp plots
for hname in variables.keys():
    if hname not in hdict.keys():
        raise ValueError("%s not found in file" % hname)
    histo = hdict[hname]
    xtitle, rebinning, x_lims = variables[hname]
    xaxis_name = histo.dense_axes()[0].name
    if rebinning != 1:
        histo = histo.rebin(xaxis_name, rebinning)
    #set_trace()
    #for ttbar_type in ['SL']:
    for ttbar_type in sorted(set([key[1] for key in histo.values().keys()])):
        ttdecay_label = ttdecay_types[ttbar_type]
        pltdir = os.path.join(outdir, ttbar_type)
        if not os.path.isdir(pltdir):
            os.makedirs(pltdir)
        for genobj, (objlabel, mass_range) in objects[ttbar_type].items():
            new_xtitle = xtitle.replace('obj', objlabel)
            if hname == 'mass':
                x_lims = mass_range
            tt_histo = histo[genobj, ttbar_type].integrate('objtype').integrate('ttdecay')
    
            fig, ax = plt.subplots()
            fig.subplots_adjust(hspace=.07)

            Plotter.plot_1D(tt_histo.values()[()], tt_histo.axis(xaxis_name).edges(), xlabel=new_xtitle, xlimits=x_lims, ax=ax, histtype='step')
            hep.cms.label(ax=ax, data=False, paper=False, year=args.year, lumi=round(lumi_to_use, 1))

                # add lepton/jet multiplicity label
            ax.text(
                0.95, 0.90, ttdecay_label,
                fontsize=rcParams['font.size'], horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes
            )

            #set_trace()
            figname = os.path.join(pltdir, '_'.join([ttbar_type, genobj, hname]))
            fig.savefig(figname)
            print('%s written' % figname)
            plt.close()


