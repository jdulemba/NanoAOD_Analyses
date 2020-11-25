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
#parser.add_argument('--nosys', action='store_true', help='Make plots without systematics and no qcd estimation')

args = parser.parse_args()

sys_to_name = systematics.sys_to_name[args.year]

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'htt_sigregion_ttdecay_fractions'

input_dir = os.path.join(proj_dir, 'results', '%s_%s' % (args.year, jobid), analyzer)
f_ext = 'TOT.coffea'
outdir = os.path.join(proj_dir, 'plots', '%s_%s' % (args.year, jobid), analyzer)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

fnames = sorted(['%s/%s' % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)])

#set_trace()
hdict = plt_tools.add_coffea_files(fnames) if len(fnames) > 1 else load(fnames[0])

ttdecay_frac_dir = os.path.join(proj_dir, 'results', '%s_%s' % (args.year, jobid), 'ttdecay_fractions')
ttdecay_fnames = sorted([os.path.join(ttdecay_frac_dir, fname) for fname in os.listdir(ttdecay_frac_dir) if fname.endswith(f_ext)])
tt_decay_dict = load(ttdecay_fnames[0])

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

variables = {
    #'mtt_vs_tlep_ctstar_abs' : ('m($t\\bar{t}$)', '|cos($\\theta^{*}_{t_{l}}$)|', linearize_binning[0], linearize_binning[1], (linearize_binning[0][0], linearize_binning[0][-1]), (linearize_binning[1][0], linearize_binning[1][-1]), True),
    'Jets_njets' : ('$n_{jets}$', 1, (0, 15), True),
    ####'mtt' : ('m($t\\bar{t}$) [GeV]', linearize_binning[0], (200., 2000.), True),
    ####'tlep_ctstar_abs' : ('|cos($\\theta^{*}_{t_{l}}$)|', linearize_binning[1], (0., 1.), True),
    #'mtt' : ('m($t\\bar{t}$) [GeV]', 4, (200., 2000.), True),
    #'tlep_ctstar_abs' : ('|cos($\\theta^{*}_{t_{l}}$)|', 1, (0., 1.), True),
    #'mthad' : ('m($t_{h}$) [GeV]', 2, (0., 300.), True),
    #'pt_thad' : ('$p_{T}$($t_{h}$) [GeV]', 2, (0., 500.), True),
    #'pt_tlep' : ('$p_{T}$($t_{l}$) [GeV]', 2, (0., 500.), True),
    #'pt_tt' : ('$p_{T}$($t\\bar{t}$) [GeV]', 2, (0., 500.), True),
    #'eta_thad' : ('$\\eta$($t_{h}$)', 2, (-4., 4.), True),
    #'eta_tlep' : ('$\\eta$($t_{l}$)', 2, (-4., 4.), True),
    #'eta_tt' : ('$\\eta$($t\\bar{t}$)', 2, (-4., 4.), True),
    #'tlep_ctstar' : ('cos($\\theta^{*}_{t_{l}}$)', 2, (-1., 1.), True),
    #'full_disc' : ('$\\lambda_{C}$', 2, (5, 25.), True),
    #'mass_disc' : ('$\\lambda_{M}$', 2, (0, 20.), True),
    #'ns_disc' : ('$\\lambda_{NS}$', 2, (3., 10.), True),
    #'Jets_pt' : ('$p_{T}$(jets) [GeV]', 2, (0., 300.), True),
    #'Jets_eta' : ('$\\eta$(jets)', 2, (-2.6, 2.6), True),
    #'Jets_LeadJet_pt' : ('$p_{T}$(leading jet) [GeV]', 2, (0., 300.), True),
    #'Jets_LeadJet_eta' : ('$\\eta$(leading jet)', 2, (-2.6, 2.6), True),
    #'Lep_pt' : ('$p_{T}$(%s) [GeV]' % objtypes['Lep'][args.lepton], 2, (0., 300.), True),
    #'Lep_eta' : ('$\\eta$(%s)' % objtypes['Lep'][args.lepton], 2, (-2.6, 2.6), True),
    #'Lep_iso' : ('pfRelIso, %s' % objtypes['Lep'][args.lepton], 1, (0., 1.), True),
    #'MT' : ('$M_{T}$ [GeV]', 1, (0., 300.), True),
    #'MET_pt' : ('$p_{T}$(MET) [GeV]', 1, (0., 300.), True),
    #'MET_phi' : ('$\phi$(MET)', 1, (-3.2, 3.2), True),
    'tt_decay' : ('$t\\bar{t}$ decay', 1, (1, 3), False),
}


#set_trace()

    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open('%s/inputs/lumis_data.json' % proj_dir).read())[args.year]
lumi_correction = load('%s/Corrections/%s/test_MC_LumiWeights.coffea' % (proj_dir, jobid))[args.year]['%ss' % args.lepton]
        # scale ttJets events, split by reconstruction type, by normal ttJets lumi correction
ttJets_permcats = ['*right', '*matchable', '*unmatchable', '*other', '*sl_tau']
names = [dataset for dataset in sorted(set([key[0] for key in hdict[sorted(variables.keys())[0]].values().keys()]))] # get dataset names in hists
ttJets_cats = [name for name in names if any([fnmatch.fnmatch(name, cat) for cat in ttJets_permcats])] # gets ttJets(_PS)_other, ...
if len(ttJets_cats) > 0:
    for tt_cat in ttJets_cats:
        ttJets_lumi_topo = '_'.join(tt_cat.split('_')[:-2]) if 'sl_tau' in tt_cat else '_'.join(tt_cat.split('_')[:-1]) # gets ttJets[SL, Had, DiLep] or ttJets_PS
        ttJets_eff_lumi = lumi_correction[ttJets_lumi_topo]
        lumi_correction.update({tt_cat: ttJets_eff_lumi})
for hname in hdict.keys():
    if 'cutflow' in hname: continue
    #set_trace()
    hdict[hname].scale(lumi_correction, axis='dataset')

## make groups based on process
process = hist.Cat("process", "Process", sorting='placement')
process_cat = "dataset"
process_groups = plt_tools.make_dataset_groups(args.lepton, args.year, samples=names)
for hname in hdict.keys():
    if 'cutflow' in hname: continue
    #set_trace()
    hdict[hname] = hdict[hname].group(process_cat, process, {'ttJets' : ['ttJetsSL', 'ttJetsDiLep', 'ttJetsHad']}) if (hname == 'tt_decay') else hdict[hname].group(process_cat, process, process_groups)


tt_decay_possibilities = {
    'Had' : [
        'Total'
    ],
    'SL' : [
        'e', 'mu', 'tau', 'tau->l', 'tau->h', 'Total'
    ],
    'DL' : [
        'e e',
        'e mu', 
        'e tau->l',
        'e tau->h',
        'mu mu', 
        'mu tau->l',
        'mu tau->h',
        'tau tau->ll',
        'tau tau->lh',
        'tau tau->hh',
        'Total',
    ],
}

    ## make plots
for hname in variables.keys():
    if hname not in hdict.keys():
        raise ValueError("%s not found in file" % hname)
    histo = hdict[hname][:, 'nosys', :, args.lepton, :, :].integrate('sys').integrate('leptype') # process, sys, jmult, leptype, btag, lepcat
    #set_trace()

    if hname == 'tt_decay':
        xtitle, rebinning, x_lims, withData = variables[hname]
        if rebinning != 1:
            xaxis_name = histo.dense_axes()[0].name
            histo = histo.rebin(xaxis_name, rebinning)
    
        for jmult in sorted(set([key[1] for key in histo.values().keys()])):
            for lepcat in sorted(set([key[3] for key in histo.values().keys()])):
                for btagregion in sorted(set([key[2] for key in histo.values().keys()])):
                    pltdir = os.path.join(outdir, args.lepton, jmult, lepcat, btagregion)
                    if not os.path.isdir(pltdir):
                        os.makedirs(pltdir)

                    hslice = histo[:, jmult, btagregion, lepcat].integrate('jmult').integrate('lepcat').integrate('btag').integrate('process')
                    plt_histo = hslice['SL e'].copy()+hslice['SL mu'].copy()+hslice['SL tau->l'].copy()+hslice['SL tau->h'].copy() + hslice['Had Total'].copy() + hslice['DL e e'].copy()+hslice['DL e mu'].copy()+hslice['DL mu mu'].copy()+hslice['DL e tau->l'].copy()+hslice['DL e tau->h'].copy()+hslice['DL mu tau->l'].copy()+hslice['DL mu tau->h'].copy()+hslice['DL tau tau->ll'].copy()+hslice['DL tau tau->lh'].copy()+hslice['DL tau tau->hh'].copy()
                    #set_trace() 
    
                    mc_opts = {
                        #'mcorder' : ['QCD', 'EWK', 'singlet', 'ttJets'] if not ttJets_cats else ['QCD', 'EWK', 'singlet', 'ttJets_other', 'ttJets_unmatchable', 'ttJets_matchable', 'ttJets_right']
                        #'maskData' : not withData
                    }
    
                    fig, ax = plt.subplots()
                    fig.subplots_adjust(hspace=.07)
                    ax = Plotter.plot_mc1d(ax, plt_histo, xlabel=xtitle, xlimits=x_lims, **mc_opts)
    
                        # add lepton/jet multiplicity label
                    #set_trace()
                    ax.text(
                        0.02, 0.88, "%s, %s\n%s" % (lep_cats[lepcat], jet_mults[jmult], btag_cats[btagregion]),
                        fontsize=rcParams['font.size']*0.75, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
                    )
                    hep.cms.cmslabel(ax=ax, data=withData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
    
                    #if hname == 'tt_decay': set_trace()
                    #set_trace()
                    figname = '%s/%s' % (pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, hname]))
                    fig.savefig(figname)
                    print('%s written' % figname)
                    plt.close()
    
                    rows = [("Lumi: %s fb^-1" % format(data_lumi_year['%ss' % args.lepton]/1000., '.1f'), args.lepton, jmult, lepcat, btagregion)]
                    rows += [("", "Decay", "Yield", "Error", "Frac of ttbar decays")]
                    
                    print(', '.join([jmult, lepcat, btagregion, hname])) 
                    sigreg_ntot = hslice.values()[('SL Total',)].sum()+hslice.values()[('DL Total',)].sum()+hslice.values()[('Had Total',)].sum()
                    for proc in tt_decay_possibilities.keys():
                        for decay in tt_decay_possibilities[proc]:
                            name, proc_yield, proc_err, proc_frac = '%s %s' % (proc, decay), hslice.values(overflow='all', sumw2=True)[('%s %s' % (proc, decay),)][0].sum(), np.sqrt(hslice.values(overflow='all', sumw2=True)[('%s %s' % (proc, decay),)][1].sum()), hslice.values(overflow='all', sumw2=True)[('%s %s' % (proc, decay),)][0].sum()/sigreg_ntot
                            rows += [("", name, format(proc_yield, '.2f'), format(proc_err, '.2f'), format(proc_frac, '.6f'))]
                        rows += [("", "", "", "", "")]
                                

                    #set_trace()
                        # nevents
                    ttdecay_name = '%s_ttdecay_yields' % '_'.join([jmult, args.lepton, lepcat, btagregion])
                    plt_tools.print_table(rows, filename='%s/%s.txt' % (pltdir, ttdecay_name), header_line=1, print_output=True)
                    print('%s/%s.txt written' % (pltdir, ttdecay_name))


        ## combine jet channels, so only depends on lepton
        #set_trace()
        lep_histo = histo[:, :, 'btagPass', 'Tight'].integrate('process').integrate('jmult').integrate('lepcat').integrate('btag')

        pltdir = os.path.join(outdir, args.lepton)
        if not os.path.isdir(pltdir):
            os.makedirs(pltdir)

        rows = [("Lumi: %s fb^-1" % format(data_lumi_year['%ss' % args.lepton]/1000., '.1f'), args.lepton, "3+ jets", "Tight", "btagPass")]
        rows += [("", "Decay", "Yield", "Error", "Frac of ttbar decays")]
                    
        print(args.lepton) 
        sigreg_ntot = lep_histo.values()[('SL Total',)].sum()+lep_histo.values()[('DL Total',)].sum()+lep_histo.values()[('Had Total',)].sum()
        for proc in tt_decay_possibilities.keys():
            ntot = lep_histo.values()[('%s Total' % proc,)].sum()
            for decay in tt_decay_possibilities[proc]:
                name, proc_yield, proc_err, proc_frac = '%s %s' % (proc, decay), lep_histo.values(overflow='all', sumw2=True)[('%s %s' % (proc, decay),)][0].sum(), np.sqrt(lep_histo.values(overflow='all', sumw2=True)[('%s %s' % (proc, decay),)][1].sum()), lep_histo.values(overflow='all', sumw2=True)[('%s %s' % (proc, decay),)][0].sum()/sigreg_ntot
                rows += [("", name, format(proc_yield, '.2f'), format(proc_err, '.2f'), format(proc_frac, '.6f'))]
            rows += [("", "", "", "", "")]
                    

        #set_trace()
            # nevents
        ttdecay_name = '%s_ttdecay_yields' % args.lepton
        plt_tools.print_table(rows, filename='%s/%s.txt' % (pltdir, ttdecay_name), header_line=1, print_output=True)
        print('%s/%s.txt written' % (pltdir, ttdecay_name))


    else:
        if histo.dense_dim() == 1:
            xtitle, rebinning, x_lims, withData = variables[hname]
            if rebinning != 1:
                xaxis_name = histo.dense_axes()[0].name
                histo = histo.rebin(xaxis_name, rebinning)
    
            #set_trace()
            ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
            for jmult in sorted(set([key[1] for key in histo.values().keys()])):
                #for lepcat in ['Tight']:
                #    for btagregion in ['btagPass']:
                for lepcat in sorted(set([key[3] for key in histo.values().keys()])):
                    for btagregion in sorted(set([key[2] for key in histo.values().keys()])):
                        pltdir = os.path.join(outdir, args.lepton, jmult, lepcat, btagregion)
                        if not os.path.isdir(pltdir):
                            os.makedirs(pltdir)
    
                        print(', '.join([jmult, lepcat, btagregion, hname])) 
    
                        hslice = histo[:, jmult, btagregion, lepcat].integrate('jmult').integrate('lepcat').integrate('btag')
                        #set_trace() 
    
                        if hname == 'Lep_iso':
                            if args.lepton == 'Muon':
                                x_lims = (0., 0.15) if lepcat == 'Tight' else (0.15, 1.)
                            if args.lepton == 'Electron':
                                x_lims = (0., 0.1) if lepcat == 'Tight' else (0., 0.5)
    
                        mc_opts = {
                            #'mcorder' : ['QCD', 'EWK', 'singlet', 'ttJets'] if not ttJets_cats else ['QCD', 'EWK', 'singlet', 'ttJets_other', 'ttJets_unmatchable', 'ttJets_matchable', 'ttJets_right']
                            #'maskData' : not withData
                        }
    
                        if withData:
                            fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                            fig.subplots_adjust(hspace=.07)
                            ax, rax = Plotter.plot_stack1d(ax, rax, hslice, xlabel=xtitle, xlimits=x_lims, **mc_opts)
                        else:
                            fig, ax = plt.subplots()
                            fig.subplots_adjust(hspace=.07)
                            ax = Plotter.plot_mc1d(ax, hslice, xlabel=xtitle, xlimits=x_lims, **mc_opts)
    
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
                        hep.cms.cmslabel(ax=ax, data=withData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
    
                        #set_trace()
                        figname = '%s/%s' % (pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, hname]))
                        fig.savefig(figname)
                        print('%s written' % figname)
                        plt.close()
    
        
        if histo.dense_dim() == 2:
            xtitle, ytitle, xrebinning, yrebinning, x_lims, y_lims, withData = variables[hname]
    
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
            for jmult in sorted(set([key[1] for key in histo.values().keys()])):
                #for lepcat in ['Tight']:
                #    for btagregion in ['btagPass']:
                for lepcat in sorted(set([key[3] for key in histo.values().keys()])):
                    for btagregion in sorted(set([key[2] for key in histo.values().keys()])):
                        pltdir = '/'.join([outdir, args.lepton, jmult, lepcat, btagregion])
                        if not os.path.isdir(pltdir):
                            os.makedirs(pltdir)
    
                        if (lepcat, btagregion) == ('Tight', 'btagPass'):
                            withData = False
                        else:
                            withData = variables[hname][-1]
                        print(', '.join([jmult, lepcat, btagregion, hname])) 
                        hslice = histo[:, jmult, btagregion, lepcat].integrate('jmult').integrate('lepcat').integrate('btag')
    
                        mc_opts = {
                            'maskData' : not withData
                        #    'mcorder' : ['QCD', 'EWK', 'singlet', 'ttJets'] if not ttJets_cats else ['QCD', 'EWK', 'singlet', 'ttJets_other', 'ttJets_unmatchable', 'ttJets_matchable', 'ttJets_right']
                        }
    
                            # make 1D projection along dense axes
                        for dax in range(2):
                            if dax == 0:
                                xlabel = ytitle
                                xlimits = y_lims
                                figname = '%s/%s' % (pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, yaxis_name, 'proj']))
                            else:
                                xlabel = '%s [GeV]' % xtitle
                                xlimits = x_lims
                                figname = '%s/%s' % (pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, xaxis_name, 'proj']))
    
                            fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                            fig.subplots_adjust(hspace=.07)
    
                            #set_trace()
                            hproj = hslice.integrate(hslice.dense_axes()[dax].name)
    
                            ax, rax = Plotter.plot_stack1d(ax, rax, hproj, xlabel=xlabel, xlimits=xlimits, **mc_opts)
    
                                # add lepton/jet multiplicity label
                            #set_trace()
                            ax.text(
                                0.02, 0.88, "%s, %s\n%s" % (lep_cats[lepcat], jet_mults[jmult], btag_cats[btagregion]),
                                fontsize=rcParams['font.size']*0.75, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
                            )
                            ax = hep.cms.cmslabel(ax=ax, data=withData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
    
                            #set_trace()
                            fig.savefig(figname)
                            print('%s written' % figname)
                            plt.close()
    
    
                            ## make 1D plots of mtt for each ctstar bin
                        for ybin in range(len(hslice.dense_axes()[1].edges())-1):
                            fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                            fig.subplots_adjust(hspace=.07)
    
                            hproj = hslice[:, :, hslice.dense_axes()[1].edges()[ybin]:hslice.dense_axes()[1].edges()[ybin+1]].integrate(hslice.dense_axes()[1].name)
    
                            ax, rax = Plotter.plot_stack1d(ax, rax, hproj, xlabel='%s [GeV]' % xtitle, xlimits=x_lims, **mc_opts)
    
                                # add lepton/jet multiplicity label, add ctstar range label
                            binlabel = '%s $\leq$ %s < %s' % (hslice.dense_axes()[1].edges()[ybin], ytitle, hslice.dense_axes()[1].edges()[ybin+1])
                            ax.text(
                                0.02, 0.88, "%s, %s\t%s\n%s" % (lep_cats[lepcat], jet_mults[jmult], binlabel, btag_cats[btagregion]),
                                fontsize=rcParams['font.size']*0.75, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
                            )
                            ax = hep.cms.cmslabel(ax=ax, data=withData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
    
                            #set_trace()
                            bintitle = '%sctstar%s' % (hslice.dense_axes()[1].edges()[ybin], hslice.dense_axes()[1].edges()[ybin+1])
                            bintitle = bintitle.replace('.', 'p')
                            figname = '%s/%s' % (pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, bintitle, 'mtt']))
                            fig.savefig(figname)
                            print('%s written' % figname)
                            plt.close()
    
    
                            # plot linearized view 
                        fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                        fig.subplots_adjust(hspace=.07)
    
                        #set_trace()
                        hline = Plotter.linearize_hist(hslice)
                        #set_trace()
                        
                        ax, rax = Plotter.plot_stack1d(ax, rax, hline, xlabel='%s $\otimes$ %s' % (xtitle, ytitle), xlimits=(0, len(hline.axis(hline.dense_axes()[0].name).edges())-1), **mc_opts)
    
                        #set_trace()
                            # draw vertical lines separating ctstar bins
                        bin_sep_lines = [hslice.values()[('EWK',)].shape[0]*ybin for ybin in range(1, hslice.values()[('EWK',)].shape[1])]
                        for binline in bin_sep_lines:
                            ax.axvline(binline, color='k', linestyle='--')
                            if rax is not None: rax.axvline(binline, color='k', linestyle='--')
                            # add lepton/jet multiplicity label
                        ax.text(
                            0.02, 0.88, "%s, %s\n%s" % (lep_cats[lepcat], jet_mults[jmult], btag_cats[btagregion]),
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



