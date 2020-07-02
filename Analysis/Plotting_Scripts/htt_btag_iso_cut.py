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
analyzer = 'htt_btag_iso_cut'

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
    'mtt' : ('m($t\\bar{t}$) [GeV]', 2, (200., 2000.), True),
    'mthad' : ('m($t_{h}$) [GeV]', 2, (0., 300.), True),
    'pt_thad' : ('$p_{T}$($t_{h}$) [GeV]', 2, (0., 500.), True),
    'pt_tlep' : ('$p_{T}$($t_{l}$) [GeV]', 2, (0., 500.), True),
    'pt_tt' : ('$p_{T}$($t\\bar{t}$) [GeV]', 2, (0., 500.), True),
    'eta_thad' : ('$\\eta$($t_{h}$)', 2, (-4., 4.), True),
    'eta_tlep' : ('$\\eta$($t_{l}$)', 2, (-4., 4.), True),
    'eta_tt' : ('$\\eta$($t\\bar{t}$)', 2, (-4., 4.), True),
    'tlep_ctstar' : ('cos($\\theta^{*}_{t_{l}}$)', 2, (-1., 1.), True),
    'tlep_ctstar_abs' : ('|cos($\\theta^{*}_{t_{l}}$)|', 2, (0., 1.), True),
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
}


    ## get plotting colors/settings
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
    for tt_cat in ttJets_cats:
        ttJets_lumi_topo = '_'.join(tt_cat.split('_')[:-1]) # gets ttJets[SL, Had, DiLep] or ttJets_PS
        ttJets_eff_lumi = lumi_correction[ttJets_lumi_topo]
        lumi_correction.update({tt_cat: ttJets_eff_lumi})
for hname in hdict.keys():
    if hname == 'cutflow': continue
    hdict[hname].scale(lumi_correction, axis='dataset')

## make groups based on process
process = hist.Cat("process", "Process", sorting='placement')
process_cat = "dataset"
process_groups = plt_tools.make_dataset_groups(args.lepton, args.year, samples=names)
for hname in hdict.keys():
    if hname == 'cutflow': continue
    hdict[hname] = hdict[hname].group(process_cat, process, process_groups)


    ## make plots
for hname in variables.keys():
    if hname not in hdict.keys():
        raise ValueError("%s not found in file" % hname)
    #set_trace()
    histo = hdict[hname][:, :, args.lepton, :, :].integrate('leptype') # process, jmult, leptype, btag, lepcat

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
                    pltdir = '/'.join([outdir, args.lepton, jmult, lepcat, btagregion])
                    if not os.path.isdir(pltdir):
                        os.makedirs(pltdir)

                    if withData:
                        fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                    else:
                        fig, ax = plt.subplots()
                    fig.subplots_adjust(hspace=.07)

                    hslice = histo[:, jmult, btagregion, lepcat].integrate('jmult').integrate('lepcat').integrate('btag')

                    if hname == 'Lep_iso':
                        if args.lepton == 'Muon':
                            x_lims = (0., 0.15) if lepcat == 'Tight' else (0.15, 1.)
                        if args.lepton == 'Electron':
                            x_lims = (0., 0.1) if lepcat == 'Tight' else (0., 0.5)

                    #if hname == 'mass_disc':
                    #    x_lims = (3., 13.) if jmult == '3Jets' else (5., 15.)
                    #if hname == 'full_disc':
                    #    x_lims = (6., 20.) if jmult == '3Jets' else (10., 22.)

                    mc_opts = {
                    #    'mcorder' : ['QCD', 'EWK', 'singlet', 'ttJets'] if not ttJets_cats else ['QCD', 'EWK', 'singlet', 'ttJets_other', 'ttJets_unmatchable', 'ttJets_matchable', 'ttJets_right']
                    }

                    if withData:
                        ax, rax = Plotter.plot_stack1d(ax, rax, hslice, xlabel=xtitle, xlimits=x_lims, **mc_opts)
                    else:
                        ax = Plotter.plot_mc1d(ax, hslice, xlabel=xtitle, xlimits=x_lims, **mc_opts)

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
                        0.02, 0.88 if withData else 0.90, "%s, %s\n%s" % (lep_cats[lepcat], jet_mults[jmult], btag_cats[btagregion]),
                        fontsize=rcParams['font.size']*0.75, 
                        horizontalalignment='left', 
                        verticalalignment='bottom', 
                        transform=ax.transAxes
                    )
                    ax = hep.cms.cmslabel(ax=ax, data=withData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))

                    #set_trace()
                    figname = '%s/%s' % (pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, hname]))
                    fig.savefig(figname)
                    print('%s written' % figname)
                    plt.close()


                        # plot with QCD estimation
                    if (btagregion  == 'btagPass') and (lepcat == 'Tight') and args.qcd_est:
                        shape_reg = 'BTAG'
                        #for norm in ['Sideband']:
                        for norm in ['ABCD', 'Sideband']:
                            qcd_name = '%s%s_Norm' % (shape_reg, norm) if norm == 'Sideband' else '%s_Norm' % norm
                            qcd_dir = '/'.join([outdir, args.lepton, jmult, 'QCD_Est', qcd_name])
                            if not os.path.isdir(qcd_dir):
                                os.makedirs(qcd_dir)

                            if withData:
                                fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                            else:
                                fig, ax = plt.subplots()
                            fig.subplots_adjust(hspace=.07)

                            iso_sb = histo[:, jmult, 'btagPass', 'Loose'].integrate('jmult').integrate('lepcat').integrate('btag')
                            btag_sb = histo[:, jmult, 'btagFail', 'Tight'].integrate('jmult').integrate('lepcat').integrate('btag')
                            double_sb = histo[:, jmult, 'btagFail', 'Loose'].integrate('jmult').integrate('lepcat').integrate('btag')
                            hslice = Plotter.QCD_Est(sig_reg=hslice, iso_sb=iso_sb, btag_sb=btag_sb, double_sb=double_sb, norm_type=norm, shape_region=shape_reg, norm_region=shape_reg if norm=='Sideband' else None)

                            if withData:
                                ax, rax = Plotter.plot_stack1d(ax, rax, hslice, xlabel=xtitle, xlimits=x_lims, **mc_opts)
                            else:
                                ax = Plotter.plot_mc1d(ax, hslice, xlabel=xtitle, xlimits=x_lims, **mc_opts)

                            if hname == 'Jets_njets':
                                print(jmult)
                                yields_txt, yields_json = Plotter.get_samples_yield_and_frac(hslice, data_lumi_year['%ss' % args.lepton]/1000., promptmc=True)
                                frac_name = '%s_yields_and_fracs_QCD_Est_%s' % ('_'.join([jmult, args.lepton, lepcat, btagregion]), qcd_name)
                                plt_tools.print_table(yields_txt, filename='%s/%s.txt' % (qcd_dir, frac_name), print_output=True)
                                print('%s/%s.txt written' % (qcd_dir, frac_name))
                                with open('%s/%s.json' % (qcd_dir, frac_name), 'w') as out:
                                    out.write(prettyjson.dumps(yields_json))

                                # add lepton/jet multiplicity label
                            ax.text(
                                0.02, 0.88 if withData else 0.90, "%s, %s\n%s" % (lep_cats[lepcat], jet_mults[jmult], btag_cats[btagregion]),
                                fontsize=rcParams['font.size']*0.75, 
                                horizontalalignment='left', 
                                verticalalignment='bottom', 
                                transform=ax.transAxes
                            )
                            ax = hep.cms.cmslabel(ax=ax, data=withData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))

                            #set_trace()
                            figname = '%s/%s_QCD_Est_%s' % (qcd_dir, '_'.join([jmult, args.lepton, lepcat, btagregion, hname]), qcd_name)
                            fig.savefig(figname)
                            print('%s written' % figname)
                            plt.close()

    

