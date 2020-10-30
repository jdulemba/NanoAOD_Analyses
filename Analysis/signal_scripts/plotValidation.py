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
parser.add_argument('SM_prod', choices=['gg', 'pp'], help='Which production of SM ttbar to use.')
parser.add_argument('genParts', choices=['status', 'lastcopy'], help='Which selection of top partons to use (LastCopyBeforeFSR or restricting range of status).')
parser.add_argument('--comp_signal', action='store_true', help='Make plots comparing signal distributions.')
parser.add_argument('--comp_weights', action='store_true', help='Make plots comparing signal distributions resulting from weights using pp->tt and gg->tt.')
parser.add_argument('--comp_SM', action='store_true', help='Make plots comparing SM ttbar distributions.')
parser.add_argument('--save_rfiles', action='store_true', help='Make plots comparing SM ttbar distributions.')
args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']


widthTOname = lambda width : str(width).replace('.', 'p')
nameTOwidth = lambda width : str(width).replace('p', '.').strip('W')

year_opts = {
    '2016' : '#e41a1c', # red
    '2017' : '#984ea3', # purple
    #'2017' : '#ff7f00', # orange
    '2018' : '#4daf4a', # green
}

genParts_dict = {
    'status' : 'GenPartStatus',
    'lastcopy' : 'LastCopyBeforeFSR',
} 

input_dir_desy = os.path.join(proj_dir, 'signal_scripts', 'inputs', 'Validation', 'desy_files')
input_dir_ur = os.path.join(proj_dir, 'signal_scripts', 'inputs', 'Validation', 'my_files', genParts_dict[args.genParts])

outdir = os.path.join(proj_dir, 'signal_scripts', 'results', 'Validation', genParts_dict[args.genParts])
if not os.path.isdir(outdir): os.makedirs(outdir)

#set_trace()
if args.comp_weights:
    gg_coffea_files = [os.path.join(input_dir_ur, fname) for fname in os.listdir(input_dir_ur) if (fname.endswith('_TOT.coffea') and ('SMgg' in fname))]
    gg_hdicts = {year:load(fname) for fname in gg_coffea_files for year in ['2016', '2017', '2018'] if re.search(year, fname)} # make dict of hists for each year
    pp_coffea_files = [os.path.join(input_dir_ur, fname) for fname in os.listdir(input_dir_ur) if (fname.endswith('_TOT.coffea') and ('SMpp' in fname))]
    pp_hdicts = {year:load(fname) for fname in pp_coffea_files for year in ['2016', '2017', '2018'] if re.search(year, fname)} # make dict of hists for each year

coffea_files = [os.path.join(input_dir_ur, fname) for fname in os.listdir(input_dir_ur) if (fname.endswith('_TOT.coffea') and ('SM%s' % args.SM_prod in fname))]
hdicts = {year:load(fname) for fname in coffea_files for year in ['2016', '2017', '2018'] if re.search(year, fname)} # make dict of hists for each year

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
lumi_correction = load(os.path.join(proj_dir, 'Corrections', jobid, 'MC_LumiWeights_IgnoreSigEvts.coffea'))
for year in hdicts.keys():
    datasets = sorted(set([key[0] for key in hdicts[year][sorted(hdicts[year].keys())[0]].values().keys()]))
    signal_datasets = [name for name in datasets if any([fnmatch.fnmatch(name, sig) for sig in ['*AtoTT*', '*HtoTT*'] ])]
    ttJets_cats = [name for name in datasets if name not in signal_datasets]
    scale_dict = {name:(fake_data_lumi/data_lumis[year]['Muons'])*lumi_correction[year]['Muons']['%s_signal' % ttJets_cat]*tt_branching_ratios['SL'] for name in datasets for ttJets_cat in ttJets_cats if ttJets_cat in name} if year == '2016' else {name:(fake_data_lumi/data_lumis[year]['Muons'])*lumi_correction[year]['Muons']['%s_signal' % ttJets_cat] for name in datasets for ttJets_cat in ttJets_cats if ttJets_cat in name}
    for histo in hdicts[year].values():
        histo.scale(scale_dict, axis='dataset')

    if args.comp_weights:
        for histo in gg_hdicts[year].values():
            histo.scale(scale_dict, axis='dataset')
        for histo in pp_hdicts[year].values():
            histo.scale(scale_dict, axis='dataset')



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
        figtitle = '$\mathsf{w_{%s \\rightarrow t \\bar t}}$, %s, m$\mathsf{_{%s}}$=%s GeV, $\mathsf{\Gamma}$/m$\mathsf{_{%s}}$=%s%%, %s' % (args.SM_prod, styles[boson]['name'], boson[0], mass.strip('M'), boson[0], nameTOwidth(width), shape)
        for var in variables.keys():
            print('   ', var)
            xlabel, rebinning, xlims, myVar, desyVar = variables[var]
            dname = os.path.join(outdir, 'SM%sTOtt_Weights' % args.SM_prod)
            if not os.path.isdir(dname): os.makedirs(dname)
            figname = os.path.join(dname, 'SM%sTOtt_%s_%s' % (args.SM_prod, signal, var))
    
            #set_trace()
            if shape == 'Res':
                desyHist = desy_file[(desyVar, 'dense_lookup')]
                sumWeights_desy = desy_file[('absolute_weight_sum', 'dense_lookup')][0].sum()
                norm_desy = abs(xsecs_desy[signal]/sumWeights_desy)*fake_data_lumi
    
                fig, ax = plt.subplots()
                fig.subplots_adjust(hspace=.07)
    
                fig_norm, ax_norm = plt.subplots()
                fig_norm.subplots_adjust(hspace=.07)
    
                    # plot expected yield vals from DESY hist
                ax = hep.plot.histplot(desyHist[0]*norm_desy, desyHist[1], ax=ax, label='DESY, l+j', color='#377eb8', histtype='errorbar') # color is blue
                    # plot normalized vals from DESY hist
                ax_norm = hep.plot.histplot(desyHist[0]/desyHist[0].sum(), desyHist[1], ax=ax_norm, label='DESY, l+j', color='#377eb8', histtype='errorbar') # color is blue
    
                #set_trace()
                sig_pattern = re.compile(r'(ttJets.*_%s*)' % signal)
                for year in sorted(hdicts.keys()):
                    myHist = hdicts[year][var]
                    histo = myHist[sig_pattern].integrate('dataset')
                    histo = histo.rebin(histo.dense_axes()[0].name, rebinning)
                    
                        # plot expected yield vals from my hist
                    ax = hep.plot.histplot(histo.values()[()], histo.dense_axes()[0].edges(), ax=ax, label='%s, x %.3f' % (year, tt_branching_ratios['SL']) if year == '2016' else '%s, l+j' % year, color=year_opts[year], histtype='errorbar')
                        # plot normalized vals from my hist
                    ax_norm = hep.plot.histplot(histo.values()[()]/histo.values()[()].sum(), histo.dense_axes()[0].edges(), ax=ax_norm, label=year if year == '2016' else '%s, l+j' % year, color=year_opts[year], histtype='errorbar')
    
                ax.set_xlim(xlims)
                ax.set_xlabel(xlabel)
                ax.set_ylabel('Events')
                ax.legend(loc='upper right')
                hep.label.lumitext(text='%s, %.1f fb$^{-1}$' % (figtitle, fake_data_lumi/1000.), ax=ax)
                hep.label._exptext(text=genParts_dict[args.genParts], ax=ax, loc=1)
                #fig.savefig('test')
                #set_trace()
                fig.savefig(figname)
                print('   %s written' % figname)
                plt.close(fig)
    
                ax_norm.set_xlim(xlims)
                ax_norm.set_xlabel(xlabel)
                ax_norm.legend(loc='upper right')
                hep.label.lumitext(text=figtitle, ax=ax_norm)
                hep.label._exptext(text=genParts_dict[args.genParts], ax=ax_norm, loc=1)
                #fig_norm.savefig('test')
                #set_trace()
                fig_norm.savefig('%s_norm' % figname)
                print('   %s_norm written' % figname)
                plt.close(fig_norm)
    
            else:
                desyHist = desy_file[(desyVar, 'dense_lookup')]
                sumWeights_desy = desy_file[('absolute_weight_sum', 'dense_lookup')][0].sum()
                norm_desy = abs(xsecs_desy[signal]/sumWeights_desy)*fake_data_lumi
    
                fig, ax = plt.subplots()
                fig.subplots_adjust(hspace=.07)
    
                        # plot expected yield vals from DESY hist
                ax = hep.plot.histplot(desyHist[0]*norm_desy, desyHist[1], ax=ax, label='DESY, l+j', color='b', histtype='errorbar')
                #    # plot unscaled vals from DESY hist
                #ax = hep.plot.histplot(desyHist[0], desyHist[1], ax=ax, label='DESY, l+j', color='b', histtype='errorbar')
    
                ymin = (desyHist[0]*norm_desy).min()
                sig_pattern_pos = re.compile(r'(ttJets.*_%s_pos)' % signal)
                sig_pattern_neg = re.compile(r'(ttJets.*_%s_neg)' % signal)
                for year in sorted(hdicts.keys()):
                    myHist = hdicts[year][var]
                    histo_pos = myHist[sig_pattern_pos].integrate('dataset')
                    histo_pos = histo_pos.rebin(histo_pos.dense_axes()[0].name, rebinning)
                    histo_neg = myHist[sig_pattern_neg].integrate('dataset')
                    histo_neg = histo_neg.rebin(histo_neg.dense_axes()[0].name, rebinning)
    
                        # plot expected yield vals from my hist
                    ax = hep.plot.histplot((histo_pos.values()[()]+histo_neg.values()[()]), histo_pos.dense_axes()[0].edges(), ax=ax, label='%s, x %.3f' % (year, tt_branching_ratios['SL']) if year == '2016' else '%s, l+j' % year, color=year_opts[year], histtype='errorbar')
                    #ax = hep.plot.histplot((histo_pos.values()[()]+histo_neg.values()[()])*norm_ur, histo_pos.dense_axes()[0].edges(), ax=ax, label=year, color=year_opts[year], histtype='errorbar')
    
                    new_ymin = ((histo_pos.values()[()]+histo_neg.values()[()])).min()
                    ymin = new_ymin if new_ymin < ymin else ymin
                        # plot normalized vals from my hist
                    #ax = hep.plot.histplot((histo_pos.values()[()]+histo_neg.values()[()]), histo_pos.dense_axes()[0].edges(), ax=ax, label=year, color=year_opts[year], histtype='errorbar')
    
                ax.set_xlim(xlims)
                ax.set_ylim(ymin*1.1, ax.get_ylim()[1])
                ax.set_xlabel(xlabel)
                ax.set_ylabel('Events')
                ax.legend(loc='upper right')
                hep.label.lumitext(text='%s, %.1f fb$^{-1}$' % (figtitle, fake_data_lumi/1000.), ax=ax)
                hep.label._exptext(text=genParts_dict[args.genParts], ax=ax, loc=1)
                #fig.savefig('test')
                #set_trace()
                fig.savefig('%s' % figname)
                print('   %s written' % figname)
                plt.close(fig)
   

if args.comp_weights:
    for signal, desy_rfile in signals.items():
        if not os.path.isfile(os.path.join(input_dir_desy, desy_rfile)): raise ValueError('%s is not found.' % os.path.join(input_dir_desy, desy_rfile))
        desy_file = convert_histo_root_file(os.path.join(input_dir_desy, desy_rfile))
        print(signal)
    
        boson, mass, width, shape = signal.split('_')
        figtitle = '%s, m$\mathsf{_{%s}}$=%s GeV, $\mathsf{\Gamma}$/m$\mathsf{_{%s}}$=%s%%, %s' % (styles[boson]['name'], boson[0], mass.strip('M'), boson[0], nameTOwidth(width), shape)

        dname = os.path.join(outdir, 'SigWeights_Comp')
        if not os.path.isdir(dname): os.makedirs(dname)
        for var in variables.keys():
            print('   ', var)
            xlabel, rebinning, xlims, myVar, desyVar = variables[var]
            figname = os.path.join(dname, 'SigWeights_Comp_%s_%s' % (signal, var))
    
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
                ax = hep.plot.histplot(desyHist[0]*norm_desy, desyHist[1], ax=ax, label='DESY, l+j', color='#377eb8', histtype='errorbar') # color is blue
                    # plot normalized vals from DESY hist
                ax_norm = hep.plot.histplot(desyHist[0]/desyHist[0].sum(), desyHist[1], ax=ax_norm, label='DESY, l+j', color='#377eb8', histtype='errorbar') # color is blue
    
                #set_trace()
                sig_pattern = re.compile(r'(ttJets.*_%s*)' % signal)
                for year in sorted(hdicts.keys()):
                        # get hist from using gg -> tt weights
                    ggHist = gg_hdicts[year][var]
                    gg_histo = ggHist[sig_pattern].integrate('dataset')
                    gg_histo = gg_histo.rebin(gg_histo.dense_axes()[0].name, rebinning)
                    
                        # get hist from using pp -> tt weights
                    ppHist = pp_hdicts[year][var]
                    pp_histo = ppHist[sig_pattern].integrate('dataset')
                    pp_histo = pp_histo.rebin(pp_histo.dense_axes()[0].name, rebinning)
                    
                        # plot expected yield vals from pp hist
                    ax = hep.plot.histplot(pp_histo.values()[()], pp_histo.dense_axes()[0].edges(), ax=ax, label='%s, x %.3f, $\mathsf{pp \\rightarrow t \\bar t}$' % (year, tt_branching_ratios['SL']) if year == '2016' else '%s, l+j, $\mathsf{pp \\rightarrow t \\bar t}$' % year, color=year_opts[year], linestyle='-', histtype='step')
                        # plot normalized vals from pp hist
                    ax_norm = hep.plot.histplot(pp_histo.values()[()]/pp_histo.values()[()].sum(), pp_histo.dense_axes()[0].edges(), ax=ax_norm, label='%s, $\mathsf{pp \\rightarrow t \\bar t}$' % year if year == '2016' else '%s, l+j, $\mathsf{pp \\rightarrow t \\bar t}$' % year, color=year_opts[year], linestyle='-', histtype='step')
    
                        # plot expected yield vals from gg hist
                    ax = hep.plot.histplot(gg_histo.values()[()], gg_histo.dense_axes()[0].edges(), ax=ax, label='%s, x %.3f, $\mathsf{gg \\rightarrow t \\bar t}$' % (year, tt_branching_ratios['SL']) if year == '2016' else '%s, l+j, $\mathsf{gg \\rightarrow t \\bar t}$' % year, color=year_opts[year], linestyle='--', histtype='step')
                        # plot normalized vals from gg hist
                    ax_norm = hep.plot.histplot(gg_histo.values()[()]/gg_histo.values()[()].sum(), gg_histo.dense_axes()[0].edges(), ax=ax_norm, label='%s, $\mathsf{gg \\rightarrow t \\bar t}$' % year if year == '2016' else '%s, l+j, $\mathsf{gg \\rightarrow t \\bar t}$' % year, color=year_opts[year], linestyle='--', histtype='step')
    
    
                ax.set_xlim(xlims)
                ax.set_xlabel(xlabel)
                ax.set_ylabel('Events')
                ax.legend(loc='upper right')
                hep.label.lumitext(text='%s, %.1f fb$^{-1}$' % (figtitle, fake_data_lumi/1000.), ax=ax)
                hep.label._exptext(text=genParts_dict[args.genParts], ax=ax, loc=1)
                #fig.savefig('test')
                #set_trace()
                fig.savefig(figname)
                print('   %s written' % figname)
                plt.close(fig)
    
                ax_norm.set_xlim(xlims)
                ax_norm.set_xlabel(xlabel)
                ax_norm.legend(loc='upper right')
                hep.label.lumitext(text=figtitle, ax=ax_norm)
                hep.label._exptext(text=genParts_dict[args.genParts], ax=ax_norm, loc=1)
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
                ax = hep.plot.histplot(desyHist[0]*norm_desy, desyHist[1], ax=ax, label='DESY, l+j', color='#377eb8', histtype='errorbar') # color is blue
                    # plot positive vals from DESY hist
                if desy_vals_pos[desy_vals_pos > 0].size > 0:
                    ax_pos = hep.plot.histplot(desy_vals_pos/abs(desy_vals_pos.sum()), desyHist[1], ax=ax_pos, label='DESY, l+j', color='#377eb8', histtype='errorbar') # color is blue
                    # plot negative vals from DESY hist
                if desy_vals_neg[desy_vals_neg < 0].size > 0:
                    ax_neg = hep.plot.histplot(desy_vals_neg/abs(desy_vals_neg.sum()), desyHist[1], ax=ax_neg, label='DESY, l+j', color='#377eb8', histtype='errorbar') # color is blue
    
                ymin = (desyHist[0]*norm_desy).min()
                sig_pattern_pos = re.compile(r'(ttJets.*_%s_pos)' % signal)
                sig_pattern_neg = re.compile(r'(ttJets.*_%s_neg)' % signal)
                #set_trace()
                for year in sorted(hdicts.keys()):
                        # get hist from using gg -> tt weights
                    ggHist = gg_hdicts[year][var]
                    gg_histo_pos = ggHist[sig_pattern_pos].integrate('dataset')
                    gg_histo_pos = gg_histo_pos.rebin(gg_histo_pos.dense_axes()[0].name, rebinning)
                    gg_histo_neg = ggHist[sig_pattern_neg].integrate('dataset')
                    gg_histo_neg = gg_histo_neg.rebin(gg_histo_neg.dense_axes()[0].name, rebinning)
                    gg_histo = gg_histo_pos.copy()
                    gg_histo.add(gg_histo_neg)

                    gg_vals_pos = gg_histo.values()[()].copy()
                    gg_vals_pos[gg_vals_pos < 0] = 0
                    gg_vals_neg = gg_histo.values()[()].copy()
                    gg_vals_neg[gg_vals_neg > 0] = 0
                    
                        # get hist from using pp -> tt weights
                    ppHist = pp_hdicts[year][var]
                    pp_histo_pos = ppHist[sig_pattern_pos].integrate('dataset')
                    pp_histo_pos = pp_histo_pos.rebin(pp_histo_pos.dense_axes()[0].name, rebinning)
                    pp_histo_neg = ppHist[sig_pattern_neg].integrate('dataset')
                    pp_histo_neg = pp_histo_neg.rebin(pp_histo_neg.dense_axes()[0].name, rebinning)
                    pp_histo = pp_histo_pos.copy()
                    pp_histo.add(pp_histo_neg)
                    
                    pp_vals_pos = pp_histo.values()[()].copy()
                    pp_vals_pos[pp_vals_pos < 0] = 0
                    pp_vals_neg = pp_histo.values()[()].copy()
                    pp_vals_neg[pp_vals_neg > 0] = 0
                    
                        # plot expected yield vals from pp hist
                    ax = hep.plot.histplot(pp_histo.values()[()], pp_histo.dense_axes()[0].edges(), ax=ax, label='%s, x %.3f, $\mathsf{pp \\rightarrow t \\bar t}$' % (year, tt_branching_ratios['SL']) if year == '2016' else '%s, l+j, $\mathsf{pp \\rightarrow t \\bar t}$' % year, color=year_opts[year], linestyle='-', histtype='step')
                        # plot normalized positive vals from pp hist
                    if pp_vals_pos[pp_vals_pos > 0].size > 0:
                        ax_pos = hep.plot.histplot(pp_vals_pos/abs(pp_vals_pos.sum()), pp_histo.dense_axes()[0].edges(), ax=ax_pos, label='%s, $\mathsf{pp \\rightarrow t \\bar t}$' % year if year == '2016' else '%s, l+j, $\mathsf{pp \\rightarrow t \\bar t}$' % year, color=year_opts[year], linestyle='-', histtype='step')
                        # plot normalized negative vals from pp hist
                    if pp_vals_neg[pp_vals_neg < 0].size > 0:
                        ax_neg = hep.plot.histplot(pp_vals_neg/abs(pp_vals_neg.sum()), pp_histo.dense_axes()[0].edges(), ax=ax_neg, label='%s, $\mathsf{pp \\rightarrow t \\bar t}$' % year if year == '2016' else '%s, l+j, $\mathsf{pp \\rightarrow t \\bar t}$' % year, color=year_opts[year], linestyle='-', histtype='step')
    
                        # plot expected yield vals from gg hist
                    ax = hep.plot.histplot(gg_histo.values()[()], gg_histo.dense_axes()[0].edges(), ax=ax, label='%s, x %.3f, $\mathsf{gg \\rightarrow t \\bar t}$' % (year, tt_branching_ratios['SL']) if year == '2016' else '%s, l+j, $\mathsf{gg \\rightarrow t \\bar t}$' % year, color=year_opts[year], linestyle='--', histtype='step')
                        # plot normalized positive vals from gg hist
                    if gg_vals_pos[gg_vals_pos > 0].size > 0:
                        ax_pos = hep.plot.histplot(gg_vals_pos/abs(gg_vals_pos.sum()), gg_histo.dense_axes()[0].edges(), ax=ax_pos, label='%s, $\mathsf{gg \\rightarrow t \\bar t}$' % year if year == '2016' else '%s, l+j, $\mathsf{gg \\rightarrow t \\bar t}$' % year, color=year_opts[year], linestyle='--', histtype='step')
                        # plot normalized negative vals from gg hist
                    if gg_vals_neg[gg_vals_neg < 0].size > 0:
                        ax_neg = hep.plot.histplot(gg_vals_neg/abs(gg_vals_neg.sum()), gg_histo.dense_axes()[0].edges(), ax=ax_neg, label='%s, $\mathsf{gg \\rightarrow t \\bar t}$' % year if year == '2016' else '%s, l+j, $\mathsf{gg \\rightarrow t \\bar t}$' % year, color=year_opts[year], linestyle='--', histtype='step')
    
                    new_ymin = np.minimum( gg_histo.values()[()].min(), pp_histo.values()[()].min() )
                    ymin = new_ymin if new_ymin < ymin else ymin
    
                ax.set_xlim(xlims)
                ax.set_ylim(ymin*1.1, ax.get_ylim()[1])
                ax.set_xlabel(xlabel)
                ax.set_ylabel('Events')
                ax.legend(loc='upper right')
                hep.label.lumitext(text='%s, %.1f fb$^{-1}$' % (figtitle, fake_data_lumi/1000.), ax=ax)
                hep.label._exptext(text=genParts_dict[args.genParts], ax=ax, loc=1)
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
                hep.label.lumitext(text='%s, %.1f fb$^{-1}$' % (figtitle, fake_data_lumi/1000.), ax=ax_pos)
                hep.label._exptext(text=genParts_dict[args.genParts], ax=ax_pos, loc=1)
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
                hep.label.lumitext(text='%s, %.1f fb$^{-1}$' % (figtitle, fake_data_lumi/1000.), ax=ax_neg)
                hep.label._exptext(text=genParts_dict[args.genParts], ax=ax_neg, loc=1)
                #fig_neg.savefig('test')
                #set_trace()
                fig_neg.savefig('%s_neg' % figname)
                print('   %s_neg written' % figname)
                plt.close(fig_neg)
   

        for var in variables_2d.keys():
            print('   ', var)
            xlabel, ylabel, xrebinning, xlims, yrebinning, ylims = variables_2d[var]
            figname = os.path.join(dname, 'SigWeights_Comp_%s_%s' % (signal, var))
    
            #set_trace()
            if shape == 'Res':
                fig, ax = plt.subplots()
                fig.subplots_adjust(hspace=.07)
    
                fig_norm, ax_norm = plt.subplots()
                fig_norm.subplots_adjust(hspace=.07)
    
                #set_trace()
                sig_pattern = re.compile(r'(ttJets.*_%s*)' % signal)
                for year in sorted(hdicts.keys()):
                        # get hist from using gg -> tt weights
                    ggHist = gg_hdicts[year][var]
                    gg_histo = ggHist[sig_pattern].integrate('dataset')
                    gg_histo = gg_histo.rebin(gg_histo.dense_axes()[0].name, xrebinning)
                    gg_histo = gg_histo.rebin(gg_histo.dense_axes()[1].name, yrebinning)
                    gg_proj = gg_histo[-0.5:0.5, :].integrate('ctstar')
                    
                        # get hist from using pp -> tt weights
                    ppHist = pp_hdicts[year][var]
                    pp_histo = ppHist[sig_pattern].integrate('dataset')
                    pp_histo = pp_histo.rebin(pp_histo.dense_axes()[0].name, xrebinning)
                    pp_histo = pp_histo.rebin(pp_histo.dense_axes()[1].name, yrebinning)
                    pp_proj = pp_histo[-0.5:0.5, :].integrate('ctstar')
                    
                        # plot expected yield vals from pp hist
                    ax = hep.plot.histplot(pp_proj.values()[()], pp_proj.dense_axes()[0].edges(), ax=ax, label='%s, x %.3f, $\mathsf{pp \\rightarrow t \\bar t}$' % (year, tt_branching_ratios['SL']) if year == '2016' else '%s, l+j, $\mathsf{pp \\rightarrow t \\bar t}$' % year, color=year_opts[year], linestyle='-', histtype='step')
                        # plot normalized vals from pp hist
                    ax_norm = hep.plot.histplot(pp_proj.values()[()]/pp_proj.values()[()].sum(), pp_proj.dense_axes()[0].edges(), ax=ax_norm, label='%s, $\mathsf{pp \\rightarrow t \\bar t}$' % year if year == '2016' else '%s, l+j, $\mathsf{pp \\rightarrow t \\bar t}$' % year, color=year_opts[year], linestyle='-', histtype='step')
    
                        # plot expected yield vals from gg hist
                    ax = hep.plot.histplot(gg_proj.values()[()], gg_proj.dense_axes()[0].edges(), ax=ax, label='%s, x %.3f, $\mathsf{gg \\rightarrow t \\bar t}$' % (year, tt_branching_ratios['SL']) if year == '2016' else '%s, l+j, $\mathsf{gg \\rightarrow t \\bar t}$' % year, color=year_opts[year], linestyle='--', histtype='step')
                        # plot normalized vals from gg hist
                    ax_norm = hep.plot.histplot(gg_proj.values()[()]/gg_proj.values()[()].sum(), gg_proj.dense_axes()[0].edges(), ax=ax_norm, label='%s, $\mathsf{gg \\rightarrow t \\bar t}$' % year if year == '2016' else '%s, l+j, $\mathsf{gg \\rightarrow t \\bar t}$' % year, color=year_opts[year], linestyle='--', histtype='step')
    
    
                ax.set_xlim(ylims)
                ax.set_xlabel(ylabel)
                ax.set_ylabel('Events')
                ax.legend(loc='upper right', title='-0.5 $<$ %s $<$ 0.5' % xlabel)
                hep.label.lumitext(text='%s, %.1f fb$^{-1}$' % (figtitle, fake_data_lumi/1000.), ax=ax)
                hep.label._exptext(text=genParts_dict[args.genParts], ax=ax, loc=1)
                #fig.savefig('test')
                #set_trace()
                fig.savefig('%s0p5' % figname)
                print('   %s0p5 written' % figname)
                plt.close(fig)
    
                ax_norm.set_xlim(ylims)
                ax_norm.set_xlabel(ylabel)
                ax_norm.set_ylabel('A.U.')
                ax_norm.legend(loc='upper right', title='-0.5 $<$ %s $<$ 0.5' % xlabel)
                hep.label.lumitext(text=figtitle, ax=ax_norm)
                hep.label._exptext(text=genParts_dict[args.genParts], ax=ax_norm, loc=1)
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
                for year in sorted(hdicts.keys()):
                        # get hist from using gg -> tt weights
                    ggHist = gg_hdicts[year][var]
                    gg_histo_pos = ggHist[sig_pattern_pos].integrate('dataset')
                    gg_histo_pos = gg_histo_pos.rebin(gg_histo_pos.dense_axes()[0].name, xrebinning)
                    gg_histo_pos = gg_histo_pos.rebin(gg_histo_pos.dense_axes()[1].name, yrebinning)
                    gg_histo_neg = ggHist[sig_pattern_neg].integrate('dataset')
                    gg_histo_neg = gg_histo_neg.rebin(gg_histo_neg.dense_axes()[0].name, xrebinning)
                    gg_histo_neg = gg_histo_neg.rebin(gg_histo_neg.dense_axes()[1].name, yrebinning)
                    gg_histo = gg_histo_pos.copy()
                    gg_histo.add(gg_histo_neg)
                    gg_proj = gg_histo[-0.5:0.5, :].integrate('ctstar')

                    gg_vals_pos = gg_proj.values()[()].copy()
                    gg_vals_pos[gg_vals_pos < 0] = 0
                    gg_vals_neg = gg_proj.values()[()].copy()
                    gg_vals_neg[gg_vals_neg > 0] = 0
                    
                        # get hist from using pp -> tt weights
                    ppHist = pp_hdicts[year][var]
                    pp_histo_pos = ppHist[sig_pattern_pos].integrate('dataset')
                    pp_histo_pos = pp_histo_pos.rebin(pp_histo_pos.dense_axes()[0].name, xrebinning)
                    pp_histo_pos = pp_histo_pos.rebin(pp_histo_pos.dense_axes()[1].name, yrebinning)
                    pp_histo_neg = ppHist[sig_pattern_neg].integrate('dataset')
                    pp_histo_neg = pp_histo_neg.rebin(pp_histo_neg.dense_axes()[0].name, xrebinning)
                    pp_histo_neg = pp_histo_neg.rebin(pp_histo_neg.dense_axes()[1].name, yrebinning)
                    pp_histo = pp_histo_pos.copy()
                    pp_histo.add(pp_histo_neg)
                    pp_proj = gg_histo[-0.5:0.5, :].integrate('ctstar')
                    
                    pp_vals_pos = pp_proj.values()[()].copy()
                    pp_vals_pos[pp_vals_pos < 0] = 0
                    pp_vals_neg = pp_proj.values()[()].copy()
                    pp_vals_neg[pp_vals_neg > 0] = 0
                    
                        # plot expected yield vals from pp hist
                    ax = hep.plot.histplot(pp_proj.values()[()], pp_proj.dense_axes()[0].edges(), ax=ax, label='%s, x %.3f, $\mathsf{pp \\rightarrow t \\bar t}$' % (year, tt_branching_ratios['SL']) if year == '2016' else '%s, l+j, $\mathsf{pp \\rightarrow t \\bar t}$' % year, color=year_opts[year], linestyle='-', histtype='step')
                        # plot normalized positive vals from pp hist
                    if pp_vals_pos[pp_vals_pos > 0].size > 0:
                        ax_pos = hep.plot.histplot(pp_vals_pos/abs(pp_vals_pos.sum()), pp_proj.dense_axes()[0].edges(), ax=ax_pos, label='%s, $\mathsf{pp \\rightarrow t \\bar t}$' % year if year == '2016' else '%s, l+j, $\mathsf{pp \\rightarrow t \\bar t}$' % year, color=year_opts[year], linestyle='-', histtype='step')
                        # plot normalized negative vals from pp hist
                    if pp_vals_neg[pp_vals_neg < 0].size > 0:
                        ax_neg = hep.plot.histplot(pp_vals_neg/abs(pp_vals_neg.sum()), pp_proj.dense_axes()[0].edges(), ax=ax_neg, label='%s, $\mathsf{pp \\rightarrow t \\bar t}$' % year if year == '2016' else '%s, l+j, $\mathsf{pp \\rightarrow t \\bar t}$' % year, color=year_opts[year], linestyle='-', histtype='step')
    
                        # plot expected yield vals from gg hist
                    ax = hep.plot.histplot(gg_proj.values()[()], gg_proj.dense_axes()[0].edges(), ax=ax, label='%s, x %.3f, $\mathsf{gg \\rightarrow t \\bar t}$' % (year, tt_branching_ratios['SL']) if year == '2016' else '%s, l+j, $\mathsf{gg \\rightarrow t \\bar t}$' % year, color=year_opts[year], linestyle='--', histtype='step')
                        # plot normalized positive vals from gg hist
                    if gg_vals_pos[gg_vals_pos > 0].size > 0:
                        ax_pos = hep.plot.histplot(gg_vals_pos/abs(gg_vals_pos.sum()), gg_proj.dense_axes()[0].edges(), ax=ax_pos, label='%s, $\mathsf{gg \\rightarrow t \\bar t}$' % year if year == '2016' else '%s, l+j, $\mathsf{gg \\rightarrow t \\bar t}$' % year, color=year_opts[year], linestyle='--', histtype='step')
                        # plot normalized negative vals from gg hist
                    if gg_vals_neg[gg_vals_neg < 0].size > 0:
                        ax_neg = hep.plot.histplot(gg_vals_neg/abs(gg_vals_neg.sum()), gg_proj.dense_axes()[0].edges(), ax=ax_neg, label='%s, $\mathsf{gg \\rightarrow t \\bar t}$' % year if year == '2016' else '%s, l+j, $\mathsf{gg \\rightarrow t \\bar t}$' % year, color=year_opts[year], linestyle='--', histtype='step')
    
                    new_ymin = np.minimum( gg_proj.values()[()].min(), pp_proj.values()[()].min() )
                    ymin = new_ymin if new_ymin < ymin else ymin
    
                ax.set_xlim(ylims)
                ax.set_ylim(ymin*1.1, ax.get_ylim()[1])
                ax.set_xlabel(ylabel)
                ax.set_ylabel('Events')
                ax.legend(loc='upper right', title='-0.5 $<$ %s $<$ 0.5' % xlabel)
                hep.label.lumitext(text='%s, %.1f fb$^{-1}$' % (figtitle, fake_data_lumi/1000.), ax=ax)
                hep.label._exptext(text=genParts_dict[args.genParts], ax=ax, loc=1)
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
                hep.label.lumitext(text='%s, %.1f fb$^{-1}$' % (figtitle, fake_data_lumi/1000.), ax=ax_pos)
                hep.label._exptext(text=genParts_dict[args.genParts], ax=ax_pos, loc=1)
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
                hep.label.lumitext(text='%s, %.1f fb$^{-1}$' % (figtitle, fake_data_lumi/1000.), ax=ax_neg)
                hep.label._exptext(text=genParts_dict[args.genParts], ax=ax_neg, loc=1)
                #fig_neg.savefig('test')
                #set_trace()
                fig_neg.savefig('%s0p5_neg' % figname)
                print('   %s0p5_neg written' % figname)
                plt.close(fig_neg)
   

 
if args.comp_SM:
    desy_fname = os.path.join(input_dir_desy, 'TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8_hist_gen_tt_no_cut.root')
    desy_file = convert_histo_root_file(desy_fname)
    figtitle = 'SM $\mathsf{t \\bar t}$'
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
        ax = hep.plot.histplot(desyHist[0]*norm_desy, desyHist[1], ax=ax, label='DESY, l+j', color='#377eb8', histtype='errorbar') # color is blue
            # plot normalized vals from DESY hist
        ax_norm = hep.plot.histplot(desyHist[0]/desyHist[0].sum(), desyHist[1], ax=ax_norm, label='DESY, l+j', color='#377eb8', histtype='errorbar') # color is blue
    
        #set_trace()
        for year in sorted(hdicts.keys()):
            SM_name = 'ttJets_PS' if year == '2016' else 'ttJetsSL'
            myHist = hdicts[year][var]
            histo = myHist[SM_name].integrate('dataset')
            histo = histo.rebin(histo.dense_axes()[0].name, rebinning)
            
                # plot expected yield vals from my hist
            ax = hep.plot.histplot(histo.values()[()], histo.dense_axes()[0].edges(), ax=ax, label='%s, x %.3f' % (year, tt_branching_ratios['SL']) if year == '2016' else '%s, l+j' % year, color=year_opts[year], histtype='errorbar')
                # plot normalized vals from my hist
            ax_norm = hep.plot.histplot(histo.values()[()]/histo.values()[()].sum(), histo.dense_axes()[0].edges(), ax=ax_norm, label=year if year == '2016' else '%s, l+j' % year, color=year_opts[year], histtype='errorbar')
    
        ax.set_xlim(xlims)
        ax.set_xlabel(xlabel)
        ax.set_ylabel('Events')
        ax.legend(loc='upper right')
        hep.label.lumitext(text='%s, %.1f fb$^{-1}$' % (figtitle, fake_data_lumi/1000.), ax=ax)
        hep.label._exptext(text=genParts_dict[args.genParts], ax=ax, loc=1)
        #fig.savefig('test')
        #set_trace()
        fig.savefig(figname)
        print('   %s written' % figname)
        plt.close(fig)
    
        ax_norm.set_xlim(xlims)
        ax_norm.set_xlabel(xlabel)
        ax_norm.legend(loc='upper right')
        hep.label.lumitext(text=figtitle, ax=ax_norm)
        hep.label._exptext(text=genParts_dict[args.genParts], ax=ax_norm, loc=1)
        #fig_norm.savefig('test')
        #set_trace()
        fig_norm.savefig('%s_norm' % figname)
        print('   %s_norm written' % figname)
        plt.close(fig_norm)


if args.save_rfiles:
    import uproot
    from coffea import hist
    tmp_rname = os.path.join(outdir, 'test.root')
    upfout = uproot.recreate(tmp_rname, compression=uproot.ZLIB(4)) if os.path.isfile(tmp_rname) else uproot.create(tmp_rname)

    for var in variables.keys():
        print('   ', var)
        xlabel, rebinning, xlims, myVar, desyVar = variables[var]

            # save SM dists
        for year in sorted(hdicts.keys()):
            SM_name = 'ttJets_PS' if year == '2016' else 'ttJetsSL'
            myHist = hdicts[year][var]
            histo = myHist[SM_name].integrate('dataset')
            histo = histo.rebin(histo.dense_axes()[0].name, rebinning)
                # find normalization for scaling to fake_data_lumi and l+jets BR
            norm_ur = (fake_data_lumi/data_lumis[year]['Muons'])*lumi_correction[year]['Muons']['%s_inds12' % SM_name]*tt_branching_ratios['SL'] if year == '2016' else (fake_data_lumi/data_lumis[year]['Muons'])*lumi_correction[year]['Muons']['%s_inds12' % SM_name]

            histo.scale(norm_ur)
            out_hname = '_'.join([SM_name, year, var])

            upfout[out_hname] = hist.export1d(histo)
            #set_trace()    


            # save signal dists
        for signal in signals.keys():
            if 'Res' in signal:
                for year in sorted(hdicts.keys()):
                    myHist = hdicts[year][var]
                    histo = myHist['%s_pos' % signal].integrate('dataset')
                    histo = histo.rebin(histo.dense_axes()[0].name, rebinning)
                        # find normalization for scaling to fake_data_lumi and l+jets BR
                    norm_ur = (fake_data_lumi/data_lumis[year]['Muons'])*lumi_correction[year]['Muons']['ttJets_PS_signal']*tt_branching_ratios['SL'] if year == '2016' else (fake_data_lumi/data_lumis[year]['Muons'])*lumi_correction[year]['Muons']['ttJetsSL_signal']
                    
                    histo.scale(norm_ur)
                    #set_trace()
                    out_hname = '_'.join([signal, year, var])

                    upfout[out_hname] = hist.export1d(histo)

            else:
                for year in sorted(hdicts.keys()):
                    myHist = hdicts[year][var]
                    histo_pos = myHist['%s_pos' % signal].integrate('dataset')
                    histo_pos = histo_pos.rebin(histo_pos.dense_axes()[0].name, rebinning)
                    histo_neg = myHist['%s_neg' % signal].integrate('dataset')
                    histo_neg = histo_neg.rebin(histo_neg.dense_axes()[0].name, rebinning)
                        # find normalization for scaling to fake_data_lumi and l+jets BR
                    norm_ur = (fake_data_lumi/data_lumis[year]['Muons'])*lumi_correction[year]['Muons']['ttJets_PS_signal']*tt_branching_ratios['SL'] if year == '2016' else (fake_data_lumi/data_lumis[year]['Muons'])*lumi_correction[year]['Muons']['ttJetsSL_signal']

                    #set_trace()
                    int_histo = histo_pos.copy()
                    int_histo.add(histo_neg)
                    int_histo.scale(norm_ur)
                    out_hname = '_'.join([signal, year, var])

                    upfout[out_hname] = hist.export1d(int_histo)

    upfout.close()
    print('%s written' % tmp_rname)
