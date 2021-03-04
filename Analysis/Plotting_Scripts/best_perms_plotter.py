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
analyzer = 'best_perms'

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('year', choices=['2016APV', '2016', '2017', '2018'] if base_jobid == 'ULnanoAOD' else ['2016', '2017', '2018'], help='Specify which year to run over')
parser.add_argument('lepton', choices=['Electron', 'Muon'], help='Choose which lepton to make plots for')
parser.add_argument('--plot', default='all', choices=['bp', 'kin', 'all'], help='Make plots best perm objects, kinematic plots, or both.')
args = parser.parse_args()

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

bp_variables = {
    'BP_pt' : ('$p_{T}$(obj) [GeV]', 2, (0., 500.)),
    'BP_eta': ('$\\eta$(obj)', 2, (-2.6, 2.6)),
    'BP_phi': ('$\\phi$(obj)', 2, (-4, 4)),
    'BP_mass': ('m(obj) [GeV]', 2, (0., 300.)),
    'BP_energy': ('E(obj) [GeV]', 2, (0., 1000.)),
}

variables = {
    #'mtt' : ('m($t\\bar{t}$) [GeV]', 2, (200., 2000.), False),
    'mtt' : ('m($t\\bar{t}$) [GeV]', 2, (200., 2000.), True),
    'mthad' : ('m($t_{h}$) [GeV]', 2, (0., 300.), True),
    'pt_thad' : ('$p_{T}$($t_{h}$) [GeV]', 2, (0., 500.), True),
    'pt_tt' : ('$p_{T}$($t\\bar{t}$) [GeV]', 2, (0., 500.), True),
    'eta_thad' : ('$\\eta$($t_{h}$)', 2, (-4., 4.), True),
    'eta_tt' : ('$\\eta$($t\\bar{t}$)', 2, (-4., 4.), True),
    'tlep_ctstar' : ('cos($\\theta^{*}_{t_{l}}$)', 2, (-1., 1.), True),
    #'tlep_ctstar_abs' : ('|cos($\\theta^{*}_{t_{l}}$)|', 2, (0., 1.), False),
    'tlep_ctstar_abs' : ('|cos($\\theta^{*}_{t_{l}}$)|', 2, (0., 1.), True),
    'full_disc' : ('$\\lambda_{C}$', 2, (5, 25.), True),
    'mass_disc' : ('$\\lambda_{M}$', 2, (0, 20.), True),
    'ns_disc' : ('$\\lambda_{NS}$', 2, (3., 10.), True),
    'Jets_pt' : ('$p_{T}$(jets) [GeV]', 2, (0., 300.), True),
    'Jets_eta' : ('$\\eta$(jets)', 2, (-2.6, 2.6), True),
    'Jets_njets' : ('$n_{jets}$', 1, (0, 15), True),
    'Jets_LeadJet_pt' : ('$p_{T}$(leading jet) [GeV]', 2, (0., 300.), True),
    'Jets_LeadJet_eta' : ('$\\eta$(leading jet)', 2, (-2.6, 2.6), True),
    'MT' : ('$M_{T}$ [GeV]', 1, (0., 300.), True),
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
    if hname == 'cutflow': continue
    hdict[hname].scale(lumi_correction, axis='dataset') # scale hists
    hdict[hname] = hdict[hname].group(process_cat, process, process_groups) # group by process
    hdict[hname] = hdict[hname][:, :, args.lepton].integrate('leptype') # only pick out specified lepton
    

    ## make bp plots
if ((args.plot == 'all') or (args.plot == 'bp')):
    for hname in bp_variables.keys():
        if hname not in hdict.keys():
            raise ValueError("%s not found in file" % hname)
        histo = hdict[hname]
    
        xtitle, rebinning, x_lims = bp_variables[hname]
        if rebinning != 1:
            xaxis_name = histo.dense_axes()[0].name
            histo = histo.rebin(xaxis_name, rebinning)
    
        ## hists should have 3 category axes (process, jet multiplicity, lepton category) followed by variable
        for jmult in sorted(set([key[1] for key in histo.values().keys()])):
            for lepcat in sorted(set([key[3] for key in histo.values().keys()])):
                for btagregion in sorted(set([key[2] for key in histo.values().keys()])):
                    pltdir = os.path.join(outdir, args.lepton, jmult, lepcat, btagregion)
                    if not os.path.isdir(pltdir):
                        os.makedirs(pltdir)
    
                    for obj in sorted(set([key[4] for key in histo.values().keys()])):
                        objlabel, mass_range = objects[obj]
                        new_xtitle = xtitle.replace('obj', objlabel)
                        if ('mass' in hname) and (hname != 'Reso_mass'):
                            x_lims = mass_range
    
                        fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                        fig.subplots_adjust(hspace=.07)
    
                        hslice = histo[:, jmult, btagregion, lepcat, obj].integrate('jmult').integrate('lepcat').integrate('btag').integrate('objtype')
                        Plotter.plot_stack1d(ax, rax, hslice, xlabel=new_xtitle, xlimits=x_lims, **{'maskData': False})
    
                            # add lepton/jet multiplicity label
                        ax.text(
                            0.02, 0.85, "%s, %s\n%s" % (lep_cats[lepcat], jet_mults[jmult], btag_cats[btagregion]),
                            horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
                        )
                        hep.cms.label(ax=ax, data=True, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
    
                        #set_trace()
                        figname = os.path.join(pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, hname, obj]))
                        fig.savefig(figname)
                        print('%s written' % figname)
                        plt.close()

    ## make kinematic plots
if ((args.plot == 'all') or (args.plot == 'kin')):
    for hname in variables.keys():
        if hname not in hdict.keys():
            raise ValueError("%s not found in file" % hname)
        histo = hdict[hname]
    
        xtitle, rebinning, x_lims, withData = variables[hname]
        if rebinning != 1:
            xaxis_name = histo.dense_axes()[0].name
            histo = histo.rebin(xaxis_name, rebinning)
    
        ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
        for jmult in sorted(set([key[1] for key in histo.values().keys()])):
            for lepcat in sorted(set([key[3] for key in histo.values().keys()])):
                for btagregion in sorted(set([key[2] for key in histo.values().keys()])):
                    pltdir = os.path.join(outdir, args.lepton, jmult, lepcat, btagregion)
                    if not os.path.isdir(pltdir):
                        os.makedirs(pltdir)
    
                    fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                    fig.subplots_adjust(hspace=.07)
    
                    hslice = histo[:, jmult, btagregion, lepcat].integrate('jmult').integrate('lepcat').integrate('btag')
                    Plotter.plot_stack1d(ax, rax, hslice, xlabel=xtitle, xlimits=x_lims, **{'maskData': not withData})
    
                    if hname == 'Jets_njets':
                        print(jmult)
                        yields_txt, yields_json = Plotter.get_samples_yield_and_frac(hslice, data_lumi_year['%ss' % args.lepton]/1000.)
                        yields_fname = os.path.join(pltdir, '%s_%s_yields_and_fracs' % (jmult, args.lepton))
                        plt_tools.print_table(yields_txt, filename='%s.txt' % yields_fname, print_output=True)
                        with open('%s.json' % yields_fname, 'w') as out:
                            out.write(prettyjson.dumps(yields_json))

                        # add lepton/jet multiplicity label
                    ax.text(
                        0.02, 0.85, "%s, %s\n%s" % (lep_cats[lepcat], jet_mults[jmult], btag_cats[btagregion]),
                        horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
                    )
                    hep.cms.label(ax=ax, data=withData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
    
                    #set_trace()
                    figname = os.path.join(pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, hname]))
                    fig.savefig(figname)
                    print('%s written' % figname)
                    plt.close()
