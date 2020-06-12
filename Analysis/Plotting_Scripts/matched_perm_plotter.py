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
import styles
import Utilities.plot_tools as plt_tools
import Utilities.prettyjson as prettyjson
#import re
from coffea import hist
import numpy as np
import fnmatch
import Plotter as Plotter

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('year', choices=['2016', '2017', '2018'], help='What year is the ntuple from.')
parser.add_argument('lepton', choices=['Electron', 'Muon'], help='Choose which lepton to make plots for')

args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'matched_perm'

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

objects = {
    'TTbar' : ('$t\\bar{t}$', (200., 2000.)),
    'THad' : ('$t_{h}$', (0., 300.)),
    'TLep' : ('$t_{l}$', (150., 200.)),
    'BHad' : ('$b_{h}$', (0., 5.)),
    'BLep' : ('$b_{l}$', (0., 5.)),
    'WHad' : ('$W_{h}$', (0., 150.)),
    'WLep' : ('$W_{l}$', (50., 100.)),
    'Lepton' : (objtypes['Lep'][args.lepton], (0., 5.)),
}

variables = {
    'Gen_pt' : ('$p_{T}$(Gen obj) [GeV]', 1, (0., 500.)),
    'Gen_eta': ('$\\eta$(Gen obj)', 1, (-2.6, 2.6)),
    'Gen_phi': ('$\\phi$(Gen obj)', 1, (-4, 4)),
    'Gen_mass': ('m(Gen obj) [GeV]', 1, (0., 300.)),
    'Gen_energy': ('E(Gen obj) [GeV]', 1, (0., 1000.)),
    'Reco_pt' : ('$p_{T}$(Reco obj) [GeV]', 1, (0., 500.)),
    'Reco_eta': ('$\\eta$(Reco obj)', 1, (-2.6, 2.6)),
    'Reco_phi': ('$\\phi$(Reco obj)', 1, (-4, 4)),
    'Reco_mass': ('m(Reco obj) [GeV]', 1, (0., 300.)),
    'Reco_energy': ('E(Reco obj) [GeV]', 1, (0., 1000.)),
    'Reso_pt' : ('$p_{T}$(Gen-Reco obj) [GeV]', 1, (-200., 200.)),
    'Reso_eta': ('$\\eta$(Gen-Reco obj)', 1, (-2.6, 2.6)),
    'Reso_phi': ('$\\phi$(Gen-Reco obj)', 1, (-4, 4)),
    'Reso_mass': ('m(Gen-Reco obj) [GeV]', 1, (-200., 200.)),
    'Reso_energy': ('E(Gen-Reco obj) [GeV]', 1, (-200., 200.)),
}


    ## get plotting colors/settings
hstyles = styles.styles
stack_fill_opts = {'alpha': 0.8, 'edgecolor':(0,0,0,.5)}
stack_error_opts = {'edgecolor':(0,0,0,.5)}

    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open('%s/inputs/lumis_data.json' % proj_dir).read())[args.year]
lumi_correction = load('%s/Corrections/%s/MC_LumiWeights_IgnoreSigEvts.coffea' % (proj_dir, jobid))[args.year]['%ss' % args.lepton]
        # scale ttJets events, split by reconstruction type, by normal ttJets lumi correction
ttJets_permcats = ['*right', '*matchable', '*unmatchable', '*other']
names = [dataset for dataset in list(set([key[0] for key in hdict[list(variables.keys())[0]].values().keys()]))] # get dataset names in hists
ttJets_cats = [name for name in names if any([fnmatch.fnmatch(name, cat) for cat in ttJets_permcats])] # gets ttJets(_PS)_other, ...
if len(ttJets_cats) > 0:
    ttJets_lumi_topo = '_'.join(ttJets_cats[0].split('_')[:-1]) # gets ttJets or ttJets_PS
    ttJets_eff_lumi = lumi_correction[ttJets_lumi_topo]
    lumi_correction.update({key: ttJets_eff_lumi for key in ttJets_cats})
for hname in hdict.keys():
    if hname == 'cutflow': continue
    hdict[hname].scale(lumi_correction, axis='dataset')


## make groups based on process
process = hist.Cat("process", "Process", sorting='placement')
process_cat = "dataset"
process_groups = plt_tools.make_dataset_groups(args.lepton, args.year, samples=names)
#set_trace()
for hname in hdict.keys():
    if hname == 'cutflow': continue
    hdict[hname] = hdict[hname].group(process_cat, process, process_groups)
    

    ## make plots
for hname in variables.keys():
    if hname not in hdict.keys():
        raise ValueError("%s not found in file" % hname)
    #set_trace()
    histo = hdict[hname][:, :, args.lepton, :, :, :].integrate('leptype') # process, jmult, leptype, lepcat, btag, object type

    if histo.dense_dim() == 1:
        xtitle, rebinning, x_lims = variables[hname]
        if rebinning != 1:
            xaxis_name = histo.dense_axes()[0].name
            histo = histo.rebin(xaxis_name, rebinning)

        #set_trace()
        ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
        for jmult in list(set([key[1] for key in histo.values().keys()])):
            for lepcat in list(set([key[3] for key in histo.values().keys()])):
                for btagregion in list(set([key[2] for key in histo.values().keys()])):
                    pltdir = '/'.join([outdir, args.lepton, jmult, lepcat, btagregion, hname.split('_')[0]])
                    if not os.path.isdir(pltdir):
                        os.makedirs(pltdir)

                    for obj in sorted(list(set([key[4] for key in histo.values().keys()]))):
                        objlabel, mass_range = objects[obj]
                        new_xtitle = xtitle.replace('obj', objlabel)
                        if ('mass' in hname) and (hname != 'Reso_mass'):
                            x_lims = mass_range

                        fig, ax = plt.subplots()
                        fig.subplots_adjust(hspace=.07)

                        #set_trace()
                        hslice = histo[:, jmult, btagregion, lepcat, obj].integrate('jmult').integrate('lepcat').integrate('btag').integrate('objtype')
                        ax = Plotter.plot_mc1d(ax, hslice, xlabel=new_xtitle, xlimits=x_lims)

                            # add lepton/jet multiplicity label
                        #set_trace()
                        ax.text(
                            0.02, 0.90, "%s, %s\n%s" % (lep_cats[lepcat], jet_mults[jmult], btag_cats[btagregion]),
                            fontsize=rcParams['font.size']*0.75, 
                            horizontalalignment='left', 
                            verticalalignment='bottom', 
                            transform=ax.transAxes
                        )
                        ax = hep.cms.cmslabel(ax=ax, data=False, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))

                        #set_trace()
                        figname = '%s/%s' % (pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, hname, obj]))
                        fig.savefig(figname)
                        print('%s written' % figname)
                        plt.close()


