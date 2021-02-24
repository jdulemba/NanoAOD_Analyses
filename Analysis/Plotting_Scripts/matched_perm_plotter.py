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
import Utilities.styles as styles
import Utilities.plot_tools as plt_tools
import Utilities.prettyjson as prettyjson
from coffea import hist
import numpy as np
import fnmatch
import Utilities.Plotter as Plotter

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
base_jobid = os.environ['base_jobid']
analyzer = 'matched_perm'

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('year', choices=['2016APV', '2016', '2017', '2018'] if base_jobid == 'ULnanoAOD' else ['2016', '2017', '2018'], help='Specify which year to run over')
parser.add_argument('lepton', choices=['Electron', 'Muon'], help='Choose which lepton to make plots for')
args = parser.parse_args()

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
    histo = hdict[hname][:, :, args.lepton, :, :, :].integrate('leptype') # process, jmult, leptype, lepcat, btag, object type

    if histo.dense_dim() == 1:
        xtitle, rebinning, x_lims = variables[hname]
        if rebinning != 1:
            xaxis_name = histo.dense_axes()[0].name
            histo = histo.rebin(xaxis_name, rebinning)

        ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
        for jmult in sorted(set([key[1] for key in histo.values().keys()])):
            for lepcat in sorted(set([key[3] for key in histo.values().keys()])):
                for btagregion in sorted(set([key[2] for key in histo.values().keys()])):
                    pltdir = os.path.join(outdir, args.lepton, jmult, lepcat, btagregion, hname.split('_')[0])
                    if not os.path.isdir(pltdir):
                        os.makedirs(pltdir)

                    for obj in sorted(sorted(set([key[4] for key in histo.values().keys()]))):
                        objlabel, mass_range = objects[obj]
                        new_xtitle = xtitle.replace('obj', objlabel)
                        if ('mass' in hname) and (hname != 'Reso_mass'):
                            x_lims = mass_range

                        fig, ax = plt.subplots()
                        fig.subplots_adjust(hspace=.07)

                        hslice = histo[:, jmult, btagregion, lepcat, obj].integrate('jmult').integrate('lepcat').integrate('btag').integrate('objtype')
                        Plotter.plot_mc1d(ax, hslice, xlabel=new_xtitle, xlimits=x_lims, ylabel='Events')

                            # add lepton/jet multiplicity label
                        ax.text(
                            0.02, 0.90, "%s, %s\n%s" % (lep_cats[lepcat], jet_mults[jmult], btag_cats[btagregion]),
                            fontsize=rcParams['font.size'], horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
                        )
                        hep.cms.label(ax=ax, data=False, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))

                        #set_trace()
                        figname = os.path.join(pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, hname, obj]))
                        fig.savefig(figname)
                        print('%s written' % figname)
                        plt.close()


