# matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend('agg')
from matplotlib import rcParams
rcParams['font.size'] = 20
#rcParams["savefig.format"] = 'pdf'
rcParams["savefig.format"] = 'png'
rcParams["savefig.bbox"] = 'tight'

from coffea.util import load#, save
from pdb import set_trace
import numpy as np
import os, re, fnmatch
from coffea.lookup_tools.root_converters import convert_histo_root_file
import Utilities.Plotter as Plotter
from Utilities.styles import styles as styles
import Utilities.prettyjson as prettyjson

from argparse import ArgumentParser
parser = ArgumentParser()
#parser.add_argument('SM_prod', choices=['gg', 'pp'], help='Which production of SM ttbar to use.')
#parser.add_argument('genParts', choices=['status', 'lastcopy'], help='Which selection of top partons to use (LastCopyBeforeFSR or restricting range of status).')
parser.add_argument('--comp_signal', action='store_true', help='Make plots comparing signal distributions.')
#parser.add_argument('--comp_weights', action='store_true', help='Make plots comparing signal distributions resulting from weights using pp->tt and gg->tt.')
parser.add_argument('--comp_SM', action='store_true', help='Make plots comparing SM ttbar distributions.')
#parser.add_argument('--save_rfiles', action='store_true', help='Make plots comparing SM ttbar distributions.')
args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']


widthTOname = lambda width : str(width).replace('.', 'p')
nameTOwidth = lambda width : str(width).replace('p', '.').strip('W')

setup_opts = {
    'nanoAODv6_LastCopyBeforeFSR' : ('#e41a1c', 'v6, LastCopyBeforeFSR'), # red
    'nanoAODv6_DiffPartons' : ('#984ea3', 'v6, GenPartStatus'), # purple
    'nanoAODv7_LastCopyBeforeFSR' : ('#ff7f00', 'v7, LastCopyBeforeFSR'), # orange
    'nanoAODv7_DiffPartons' : ('#4daf4a', 'v7, GenPartStatus'), # green
}

genParts_dict = {
    'status' : 'GenPartStatus',
    'lastcopy' : 'LastCopyBeforeFSR',
} 

input_dir_desy = os.path.join(proj_dir, 'signal_scripts', 'inputs', 'Validation', 'desy_files')

hdicts = {
    'nanoAODv6_LastCopyBeforeFSR' : load(os.path.join(proj_dir, 'signal_scripts', 'inputs', 'Validation', 'my_files', 'LastCopyBeforeFSR', 'BATCH_SigValidation_SMpp_ttJetsSL_2017_26Oct2020_jpt30_ljpt50_MT40_cutBasedEl_TOT.coffea')),
    'nanoAODv6_DiffPartons' : load(os.path.join(proj_dir, 'signal_scripts', 'inputs', 'Validation', 'my_files', 'GenPartStatus', 'BATCH_SigValidation_DiffPartons_SMpp_ttJetsSL_2017_28Oct2020_jpt30_ljpt50_MT40_cutBasedEl_TOT.coffea')),
    'nanoAODv7_LastCopyBeforeFSR' : load(os.path.join(proj_dir, 'signal_scripts', 'inputs', 'Validation', 'my_files', 'NanoAODv7', 'BATCH_SigValidation_LastCopy_SMpp_ttJetsSL_2017_29Oct2020_NanoAODv7_TOT.coffea')),
    'nanoAODv7_DiffPartons' : load(os.path.join(proj_dir, 'signal_scripts', 'inputs', 'Validation', 'my_files', 'NanoAODv7', 'BATCH_SigValidation_DiffPartons_SMpp_ttJetsSL_2017_28Oct2020_NanoAODv7_TOT.coffea')),
}

outdir = os.path.join(proj_dir, 'signal_scripts', 'results', 'Validation', 'nanoAOD_Comp')
if not os.path.isdir(outdir): os.makedirs(outdir)

#set_trace()

xsecs_desy = prettyjson.loads(open(os.path.join(input_dir_desy, 'xsecs.json')).read())

    # values taken from https://pdg.lbl.gov/2019/reviews/rpp2019-rev-top-quark.pdf
tt_branching_ratios = {
    'Had' : 0.457,
    'SL' : 0.438,
    'DiLep' : 0.105,
}

    # scale by lumi (and SL branching ratio in 2016)
fake_data_lumi = 50000. # 50 fb-1
data_lumis = prettyjson.loads(open(os.path.join(proj_dir, 'inputs', 'lumis_data.json')).read())
lumi_correction_nanoAODv6 = load(os.path.join(proj_dir, 'Corrections', 'NanoAODv6', 'MC_LumiWeights_IgnoreSigEvts.coffea'))
lumi_correction_nanoAODv7 = load(os.path.join(proj_dir, 'Corrections', 'NanoAODv7', 'MC_LumiWeights_IgnoreSigEvts.coffea'))
for setup in hdicts.keys():
    lumi_correction = lumi_correction_nanoAODv6 if 'nanoAODv6' in setup else lumi_correction_nanoAODv7
    datasets = sorted(set([key[0] for key in hdicts[setup][sorted(hdicts[setup].keys())[0]].values().keys()]))
    signal_datasets = [name for name in datasets if any([fnmatch.fnmatch(name, sig) for sig in ['*AtoTT*', '*HtoTT*'] ])]
    ttJets_cats = [name for name in datasets if name not in signal_datasets]
    scale_dict = {name:(fake_data_lumi/data_lumis['2017']['Muons'])*lumi_correction['2017']['Muons']['%s_signal' % ttJets_cat] for name in datasets for ttJets_cat in ttJets_cats if ttJets_cat in name} 
    for histo in hdicts[setup].values():
        histo.scale(scale_dict, axis='dataset')

#set_trace()

signals = {
    'AtoTT_M600_W2p5_Int' : 'A_int-m600-relw2p5_lj_hist_gen_tt_no_cut.root',
    'AtoTT_M600_W2p5_Res' : 'A_res-m600-relw2p5_lj_hist_gen_tt_no_cut.root',
    'AtoTT_M600_W5_Int' : 'A_int-m600-relw5_lj_hist_gen_tt_no_cut.root',
    'AtoTT_M600_W5_Res' : 'A_res-m600-relw5_lj_hist_gen_tt_no_cut.root',
    'HtoTT_M600_W2p5_Int' : 'H_int-m600-relw2p5_lj_hist_gen_tt_no_cut.root',
    'HtoTT_M600_W2p5_Res' : 'H_res-m600-relw2p5_lj_hist_gen_tt_no_cut.root',
    'HtoTT_M600_W5_Int' : 'H_int-m600-relw5_lj_hist_gen_tt_no_cut.root',
    'HtoTT_M600_W5_Res' : 'H_res-m600-relw5_lj_hist_gen_tt_no_cut.root',
}

variables = {
    'top_pt' : ('$\mathsf{p_{T}(t)}$ [GeV]', 2, (0., 600.), 'top_pt', 'gen_top_pt_no_cut'),
    'top_rapidity' : ('$\mathsf{y}(t)}$', 2, (-4., 4.), 'top_rapidity', 'gen_top_rapidity_no_cut'), 
    'top_mass' : ('m$\mathsf{_{t}}$ [GeV]', 1, (160., 185.), 'top_mass', 'gen_top_m_no_cut'),
    'antitop_pt' : ('$\mathsf{p_{T}(\\bar t)}$ [GeV]', 2, (0., 600.), 'antitop_pt', 'gen_antitop_pt_no_cut'),
    'antitop_rapidity' : ('$\mathsf{y}(\\bar t)}$', 2, (-4., 4.), 'antitop_rapidity', 'gen_antitop_rapidity_no_cut'),
    'antitop_mass' : ('m$\mathsf{_{\\bar t}}$ [GeV]', 1, (160., 185.), 'antitop_mass', 'gen_antitop_m_no_cut'),
    'mtt' : ('$\mathsf{m_{t \\bar t}}$ [GeV]', 1, (300., 800.), 'mtt', 'gen_TT_m_no_cut'), # binning is np.linspace(250., 1500., 171)
    'top_ctstar' : ('$\mathsf{cos\\theta^{*}_{t}}$', 1, (-1., 1.), 'top_ctstar', 'gen_cpTTT_no_cut'), # binning is np.linspace(-1, 1, 25)
    'ttbar_pt' : ('$\mathsf{p_{T}(t \\bar t)}$ [GeV]', 1, (0., 600.), 'ttbar_pt', 'gen_TT_pt_no_cut'),
    'ttbar_rapidity' : ('$\mathsf{y}(t\\bar t)}$', 2, (-4., 4.), 'ttbar_rapidity', 'gen_TT_rapidity_no_cut'),
}
variables_2d = {
    'top_rap_vs_ctstar' : ('$\mathsf{cos\\theta^{*}_{t}}$', '$\mathsf{y}(t)}$', 1, (-1., 1.), 2, (-4., 4.)), 
}

   

if args.comp_signal:
    for signal, desy_rfile in signals.items():
        if not os.path.isfile(os.path.join(input_dir_desy, desy_rfile)): raise ValueError('%s is not found.' % os.path.join(input_dir_desy, desy_rfile))
        desy_file = convert_histo_root_file(os.path.join(input_dir_desy, desy_rfile))
        print(signal)
    
        boson, mass, width, shape = signal.split('_')
        figtitle = '%s, m$\mathsf{_{%s}}$=%s GeV, $\mathsf{\Gamma}$/m$\mathsf{_{%s}}$=%s%%, %s' % (styles[boson]['name'], boson[0], mass.strip('M'), boson[0], nameTOwidth(width), shape)

        dname = os.path.join(outdir, 'Signal')
        if not os.path.isdir(dname): os.makedirs(dname)
        for var in variables.keys():
            print('   ', var)
            xlabel, rebinning, xlims, myVar, desyVar = variables[var]
            figname = os.path.join(dname, '%s_%s' % (signal, var))
    
            #set_trace()
            desyHist = desy_file[(desyVar, 'dense_lookup')]
            sumWeights_desy = desy_file[('absolute_weight_sum', 'dense_lookup')][0].sum()
            norm_desy = abs(xsecs_desy[signal]/sumWeights_desy)*fake_data_lumi
    
            if shape == 'Res':
                fig, ax = plt.subplots()
                fig.subplots_adjust(hspace=.07)
    
                fig_norm, ax_norm = plt.subplots()
                fig_norm.subplots_adjust(hspace=.07)
    
                    # plot expected yield vals from DESY hist
                ax = hep.plot.histplot(desyHist[0]*norm_desy, desyHist[1], ax=ax, label='DESY', color='#377eb8', histtype='errorbar') # color is blue
                    # plot normalized vals from DESY hist
                ax_norm = hep.plot.histplot(desyHist[0]/desyHist[0].sum(), desyHist[1], ax=ax_norm, label='DESY', color='#377eb8', histtype='errorbar') # color is blue
    
                #set_trace()
                sig_pattern = re.compile(r'(ttJets.*_%s*)' % signal)
                for setup in sorted(hdicts.keys()):
                    Hist = hdicts[setup][var]
                    histo = Hist[sig_pattern].integrate('dataset')
                    histo = histo.rebin(histo.dense_axes()[0].name, rebinning)
                    
                        # plot expected yield vals from pp hist
                    ax = hep.plot.histplot(histo.values()[()], histo.dense_axes()[0].edges(), ax=ax, color=setup_opts[setup][0], label=setup_opts[setup][1], linestyle='-', histtype='step')
                        # plot normalized vals from pp hist
                    ax_norm = hep.plot.histplot(histo.values()[()]/histo.values()[()].sum(), histo.dense_axes()[0].edges(), ax=ax_norm, color=setup_opts[setup][0], label=setup_opts[setup][1], linestyle='-', histtype='step')
    
    
                ax.set_xlim(xlims)
                ax.set_xlabel(xlabel)
                ax.set_ylabel('Events')
                ax.legend(loc='upper right')
                ax = hep.cms.cmslabel(ax=ax, llabel='Prelim.',  rlabel='%s, %.1f fb$^{-1}$' % (figtitle, fake_data_lumi/1000.))
                #fig.savefig('test')
                #set_trace()
                fig.savefig(figname)
                print('   %s written' % figname)
                plt.close(fig)
    
                ax_norm.set_xlim(xlims)
                ax_norm.set_xlabel(xlabel)
                ax_norm.set_ylabel('A.U.')
                ax_norm.legend(loc='upper right')
                ax_norm = hep.cms.cmslabel(ax=ax_norm, llabel='Prelim.', rlabel=figtitle)
                #fig_norm.savefig('test')
                #set_trace()
                fig_norm.savefig('%s_norm' % figname)
                print('   %s_norm written' % figname)
                plt.close(fig_norm)
    
            else:
                desy_vals_pos = desyHist[0].copy()
                desy_vals_pos[desy_vals_pos < 0] = 0.
                desy_vals_neg = desyHist[0].copy()
                desy_vals_neg[desy_vals_neg > 0] = 0.

                    ## total interference
                fig, ax = plt.subplots()
                fig.subplots_adjust(hspace=.07)
    
                    ## positive interference
                fig_pos, ax_pos = plt.subplots()
                fig_pos.subplots_adjust(hspace=.07)
    
                    ## negative interference
                fig_neg, ax_neg = plt.subplots()
                fig_neg.subplots_adjust(hspace=.07)
    
                    # plot expected yield vals from DESY hist
                ax = hep.plot.histplot(desyHist[0]*norm_desy, desyHist[1], ax=ax, label='DESY', color='#377eb8', histtype='errorbar') # color is blue
                    # plot positive vals from DESY hist
                if desy_vals_pos[desy_vals_pos > 0].size > 0:
                    ax_pos = hep.plot.histplot(desy_vals_pos/abs(desy_vals_pos.sum()), desyHist[1], ax=ax_pos, label='DESY', color='#377eb8', histtype='errorbar') # color is blue
                    # plot negative vals from DESY hist
                if desy_vals_neg[desy_vals_neg < 0].size > 0:
                    ax_neg = hep.plot.histplot(desy_vals_neg/abs(desy_vals_neg.sum()), desyHist[1], ax=ax_neg, label='DESY', color='#377eb8', histtype='errorbar') # color is blue
    
                ymin = (desyHist[0]*norm_desy).min()
                sig_pattern_pos = re.compile(r'(ttJets.*_%s_pos)' % signal)
                sig_pattern_neg = re.compile(r'(ttJets.*_%s_neg)' % signal)
                #set_trace()
                for setup in sorted(hdicts.keys()):
                    Hist = hdicts[setup][var]
                    histo_pos = Hist[sig_pattern_pos].integrate('dataset')
                    histo_pos = histo_pos.rebin(histo_pos.dense_axes()[0].name, rebinning)
                    histo_neg = Hist[sig_pattern_neg].integrate('dataset')
                    histo_neg = histo_neg.rebin(histo_neg.dense_axes()[0].name, rebinning)
                    histo = histo_pos.copy()
                    histo.add(histo_neg)
                    
                    vals_pos = histo.values()[()].copy()
                    vals_pos[vals_pos < 0] = 0
                    vals_neg = histo.values()[()].copy()
                    vals_neg[vals_neg > 0] = 0

                        # plot expected yield vals from pp hist
                    ax = hep.plot.histplot(histo.values()[()], histo.dense_axes()[0].edges(), ax=ax, color=setup_opts[setup][0], label=setup_opts[setup][1], linestyle='-', histtype='step')
                        # plot normalized positive vals from pp hist
                    if vals_pos[vals_pos > 0].size > 0:
                        ax_pos = hep.plot.histplot(vals_pos/abs(vals_pos.sum()), histo.dense_axes()[0].edges(), ax=ax_pos, color=setup_opts[setup][0], label=setup_opts[setup][1], linestyle='-', histtype='step')
                        # plot normalized negative vals from pp hist
                    if vals_neg[vals_neg < 0].size > 0:
                        ax_neg = hep.plot.histplot(vals_neg/abs(vals_neg.sum()), histo.dense_axes()[0].edges(), ax=ax_neg, color=setup_opts[setup][0], label=setup_opts[setup][1], linestyle='-', histtype='step')
    
                    new_ymin = histo.values()[()].min()
                    ymin = new_ymin if new_ymin < ymin else ymin
    
                ax.set_xlim(xlims)
                ax.set_ylim(ymin*1.1, ax.get_ylim()[1])
                ax.set_xlabel(xlabel)
                ax.set_ylabel('Events')
                ax.legend(loc='upper right')
                ax = hep.cms.cmslabel(ax=ax, llabel='Prelim.',  rlabel='%s, %.1f fb$^{-1}$' % (figtitle, fake_data_lumi/1000.))
                #fig.savefig('test')
                #set_trace()
                fig.savefig('%s' % figname)
                print('   %s written' % figname)
                plt.close(fig)

                    # positive interference   
                ax_pos.set_xlim(xlims)
                ax_pos.set_xlabel(xlabel)
                ax_pos.set_ylabel('A.U.')
                ax_pos.legend(loc='upper right', title='w $>$ 0')
                ax_pos = hep.cms.cmslabel(ax=ax_pos, llabel='Prelim.',  rlabel=figtitle)
                #fig_pos.savefig('test')
                #set_trace()
                fig_pos.savefig('%s_pos' % figname)
                print('   %s_pos written' % figname)
                plt.close(fig_pos)
   
                    # negative interference   
                ax_neg.set_xlim(xlims)
                ax_neg.set_xlabel(xlabel)
                ax_neg.set_ylabel('A.U.')
                ax_neg.legend(loc='upper right', title='w $<$ 0')
                ax_neg = hep.cms.cmslabel(ax=ax_neg, llabel='Prelim.',  rlabel=figtitle)
                #fig_neg.savefig('test')
                #set_trace()
                fig_neg.savefig('%s_neg' % figname)
                print('   %s_neg written' % figname)
                plt.close(fig_neg)
   

        for var in variables_2d.keys():
            print('   ', var)
            xlabel, ylabel, xrebinning, xlims, yrebinning, ylims = variables_2d[var]
            figname = os.path.join(dname, '%s_%s' % (signal, var))
    
            #set_trace()
            if shape == 'Res':
                fig, ax = plt.subplots()
                fig.subplots_adjust(hspace=.07)
    
                fig_norm, ax_norm = plt.subplots()
                fig_norm.subplots_adjust(hspace=.07)
    
                #set_trace()
                sig_pattern = re.compile(r'(ttJets.*_%s*)' % signal)
                for setup in sorted(hdicts.keys()):
                    Hist = hdicts[setup][var]
                    histo = Hist[sig_pattern].integrate('dataset')
                    histo = histo.rebin(histo.dense_axes()[0].name, xrebinning)
                    histo = histo.rebin(histo.dense_axes()[1].name, yrebinning)
                    proj = histo[-0.5:0.5, :].integrate('ctstar')
                    
                        # plot expected yield vals from pp hist
                    ax = hep.plot.histplot(proj.values()[()], proj.dense_axes()[0].edges(), ax=ax, color=setup_opts[setup][0], label=setup_opts[setup][1], linestyle='-', histtype='step')
                        # plot normalized vals from pp hist
                    ax_norm = hep.plot.histplot(proj.values()[()]/proj.values()[()].sum(), proj.dense_axes()[0].edges(), ax=ax_norm, color=setup_opts[setup][0], label=setup_opts[setup][1], linestyle='-', histtype='step')
    
    
                ax.set_xlim(ylims)
                ax.set_xlabel(ylabel)
                ax.set_ylabel('Events')
                ax.legend(loc='upper right', title='-0.5 $<$ %s $<$ 0.5' % xlabel)
                hep.cms.cmslabel(ax=ax, llabel='Prelim.',  rlabel='%s, %.1f fb$^{-1}$' % (figtitle, fake_data_lumi/1000.))
                #hep.label.lumitext(text='%s, %.1f fb$^{-1}$' % (figtitle, fake_data_lumi/1000.), ax=ax)
                #hep.label._exptext(text=genParts_dict[args.genParts], ax=ax, loc=1)
                #fig.savefig('test')
                #set_trace()
                fig.savefig('%s0p5' % figname)
                print('   %s0p5 written' % figname)
                plt.close(fig)
    
                ax_norm.set_xlim(ylims)
                ax_norm.set_xlabel(ylabel)
                ax_norm.set_ylabel('A.U.')
                ax_norm.legend(loc='upper right', title='-0.5 $<$ %s $<$ 0.5' % xlabel)
                hep.cms.cmslabel(ax=ax_norm, llabel='Prelim.',  rlabel=figtitle)
                #hep.label.lumitext(text=figtitle, ax=ax_norm)
                #hep.label._exptext(text=genParts_dict[args.genParts], ax=ax_norm, loc=1)
                #fig_norm.savefig('test')
                #set_trace()
                fig_norm.savefig('%s0p5_norm' % figname)
                print('   %s0p5_norm written' % figname)
                plt.close(fig_norm)
    
            else:

                    ## total interference
                fig, ax = plt.subplots()
                fig.subplots_adjust(hspace=.07)
    
                    ## positive interference
                fig_pos, ax_pos = plt.subplots()
                fig_pos.subplots_adjust(hspace=.07)
    
                    ## negative interference
                fig_neg, ax_neg = plt.subplots()
                fig_neg.subplots_adjust(hspace=.07)
    
                ymin = 0
                sig_pattern_pos = re.compile(r'(ttJets.*_%s_pos)' % signal)
                sig_pattern_neg = re.compile(r'(ttJets.*_%s_neg)' % signal)
                #set_trace()
                for setup in sorted(hdicts.keys()):
                    Hist = hdicts[setup][var]
                    histo_pos = Hist[sig_pattern_pos].integrate('dataset')
                    histo_pos = histo_pos.rebin(histo_pos.dense_axes()[0].name, xrebinning)
                    histo_pos = histo_pos.rebin(histo_pos.dense_axes()[1].name, yrebinning)
                    histo_neg = Hist[sig_pattern_neg].integrate('dataset')
                    histo_neg = histo_neg.rebin(histo_neg.dense_axes()[0].name, xrebinning)
                    histo_neg = histo_neg.rebin(histo_neg.dense_axes()[1].name, yrebinning)
                    histo = histo_pos.copy()
                    histo.add(histo_neg)
                    proj = histo[-0.5:0.5, :].integrate('ctstar')
                    
                    vals_pos = proj.values()[()].copy()
                    vals_pos[vals_pos < 0] = 0
                    vals_neg = proj.values()[()].copy()
                    vals_neg[vals_neg > 0] = 0
                    
                        # plot expected yield vals from pp hist
                    ax = hep.plot.histplot(proj.values()[()], proj.dense_axes()[0].edges(), ax=ax, color=setup_opts[setup][0], label=setup_opts[setup][1], linestyle='-', histtype='step')
                        # plot normalized positive vals from pp hist
                    if vals_pos[vals_pos > 0].size > 0:
                        ax_pos = hep.plot.histplot(vals_pos/abs(vals_pos.sum()), proj.dense_axes()[0].edges(), ax=ax_pos, color=setup_opts[setup][0], label=setup_opts[setup][1], linestyle='-', histtype='step')
                        # plot normalized negative vals from pp hist
                    if vals_neg[vals_neg < 0].size > 0:
                        ax_neg = hep.plot.histplot(vals_neg/abs(vals_neg.sum()), proj.dense_axes()[0].edges(), ax=ax_neg, color=setup_opts[setup][0], label=setup_opts[setup][1], linestyle='-', histtype='step')
    
    
                    new_ymin = proj.values()[()].min()
                    ymin = new_ymin if new_ymin < ymin else ymin
    
                ax.set_xlim(ylims)
                ax.set_ylim(ymin*1.1, ax.get_ylim()[1])
                ax.set_xlabel(ylabel)
                ax.set_ylabel('Events')
                ax.legend(loc='upper right', title='-0.5 $<$ %s $<$ 0.5' % xlabel)
                hep.cms.cmslabel(ax=ax, llabel='Prelim.',  rlabel='%s, %.1f fb$^{-1}$' % (figtitle, fake_data_lumi/1000.))
                #hep.label.lumitext(text='%s, %.1f fb$^{-1}$' % (figtitle, fake_data_lumi/1000.), ax=ax)
                #hep.label._exptext(text=genParts_dict[args.genParts], ax=ax, loc=1)
                #fig.savefig('test')
                #set_trace()
                fig.savefig('%s0p5' % figname)
                print('   %s0p5 written' % figname)
                plt.close(fig)

                    # positive interference   
                ax_pos.set_xlim(ylims)
                ax_pos.set_xlabel(ylabel)
                ax_pos.set_ylabel('A.U.')
                ax_pos.legend(loc='upper right', title='-0.5 $<$ %s $<$ 0.5\nw $>$ 0' % xlabel)
                hep.cms.cmslabel(ax=ax_pos, llabel='Prelim.',  rlabel=figtitle)
                #hep.label.lumitext(text='%s, %.1f fb$^{-1}$' % (figtitle, fake_data_lumi/1000.), ax=ax_pos)
                #hep.label._exptext(text=genParts_dict[args.genParts], ax=ax_pos, loc=1)
                #fig_pos.savefig('test')
                #set_trace()
                fig_pos.savefig('%s0p5_pos' % figname)
                print('   %s0p5_pos written' % figname)
                plt.close(fig_pos)
   
                    # negative interference   
                ax_neg.set_xlim(ylims)
                ax_neg.set_xlabel(ylabel)
                ax_neg.set_ylabel('A.U.')
                ax_neg.legend(loc='upper right', title='-0.5 $<$ %s $<$ 0.5\nw $<$ 0' % xlabel)
                hep.cms.cmslabel(ax=ax_neg, llabel='Prelim.',  rlabel=figtitle)
                #hep.label.lumitext(text='%s, %.1f fb$^{-1}$' % (figtitle, fake_data_lumi/1000.), ax=ax_neg)
                #hep.label._exptext(text=genParts_dict[args.genParts], ax=ax_neg, loc=1)
                #fig_neg.savefig('test')
                #set_trace()
                fig_neg.savefig('%s0p5_neg' % figname)
                print('   %s0p5_neg written' % figname)
                plt.close(fig_neg)
   

 
if args.comp_SM:
    desy_fname = os.path.join(input_dir_desy, 'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_hist_gen_tt_no_cut.root')
    desy_file = convert_histo_root_file(desy_fname)
    figtitle = '2017 SM $\mathsf{t \\bar t}$, l+j'
    dname = os.path.join(outdir, 'SMttbar')
    if not os.path.isdir(dname): os.makedirs(dname)
    for var in variables.keys():
        print('   ', var)
        xlabel, rebinning, xlims, myVar, desyVar = variables[var]
        figname = os.path.join(dname, 'SMttbar_%s' % var)
    
        #set_trace()
        desyHist = desy_file[(desyVar, 'dense_lookup')]
        sumWeights_desy = desy_file[('absolute_weight_sum', 'dense_lookup')][0].sum()
        norm_desy = abs(364.31088/sumWeights_desy)*fake_data_lumi
    
        fig, ax = plt.subplots()
        fig.subplots_adjust(hspace=.07)
    
        fig_norm, ax_norm = plt.subplots()
        fig_norm.subplots_adjust(hspace=.07)
    
            # plot expected yield vals from DESY hist
        ax = hep.plot.histplot(desyHist[0]*norm_desy, desyHist[1], ax=ax, label='DESY', color='#377eb8', histtype='errorbar') # color is blue
            # plot normalized vals from DESY hist
        ax_norm = hep.plot.histplot(desyHist[0]/desyHist[0].sum(), desyHist[1], ax=ax_norm, label='DESY', color='#377eb8', histtype='errorbar') # color is blue
    
        #set_trace()
        for setup in sorted(hdicts.keys()):
            myHist = hdicts[setup][var]
            histo = myHist['ttJetsSL'].integrate('dataset')
            histo = histo.rebin(histo.dense_axes()[0].name, rebinning)
            
                # plot expected yield vals from my hist
            ax = hep.plot.histplot(histo.values()[()], histo.dense_axes()[0].edges(), ax=ax, color=setup_opts[setup][0], label=setup_opts[setup][1], histtype='errorbar')
                # plot normalized vals from my hist
            ax_norm = hep.plot.histplot(histo.values()[()]/histo.values()[()].sum(), histo.dense_axes()[0].edges(), ax=ax_norm, color=setup_opts[setup][0], label=setup_opts[setup][1], histtype='errorbar')
    
        ax.set_xlim(xlims)
        ax.set_xlabel(xlabel)
        ax.set_ylabel('Events')
        ax.legend(loc='upper right')
        ax = hep.cms.cmslabel(ax=ax, rlabel='%s, %.1f fb$^{-1}$' % (figtitle, fake_data_lumi/1000.))
        #fig.savefig('test')
        #set_trace()
        fig.savefig(figname)
        print('   %s written' % figname)
        plt.close(fig)
    
        ax_norm.set_xlim(xlims)
        ax_norm.set_xlabel(xlabel)
        ax_norm.set_ylabel('A.U.')
        ax_norm.legend(loc='upper right')
        ax_norm = hep.cms.cmslabel(ax=ax_norm, rlabel=figtitle)
        #fig_norm.savefig('test')
        #set_trace()
        fig_norm.savefig('%s_norm' % figname)
        print('   %s_norm written' % figname)
        plt.close(fig_norm)


