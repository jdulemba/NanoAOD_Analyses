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
from coffea import hist
import numpy as np
import Utilities.Plotter as Plotter
import Utilities.systematics as systematics
from coffea.hist import plot

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('year', choices=['2016', '2017', '2018'], help='What year is the ntuple from.')
parser.add_argument('lepton', choices=['Electron', 'Muon'], help='Choose which lepton to make plots for')
parser.add_argument('--indiv', action='store_true', help='Plot individual sys variations')
parser.add_argument('--comp', action='store_true', help='Plot up/down sys variations with nominal mc')
#parser.add_argument('--ratio', action='store_true', help='Make sys variation/nominal mc ratio plots')

args = parser.parse_args()

#if (not args.comp) and (not args.ratio) and (not args.indiv):
#    raise ValueError("No plots chosen to make")

#sys_to_name = systematics.sys_to_name[args.year]

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'htt_btag_iso_cut'

input_dir = '/'.join([proj_dir, 'plots', '%s_%s/' % (args.year, jobid), analyzer, 'Templates'])
outdir = '/'.join([input_dir, 'plots', args.lepton])
if not os.path.isdir(outdir):
    os.makedirs(outdir)

jet_mults = {
    '3Jets' : '3 jets',
    '4PJets' : '4+ jets'
}

leptypes = {
    'Muon' : '$\\mu$',
    'Electron' : '$e$',
}
lepdir = 'mujets' if args.lepton == 'Muon' else 'ejets'

baseSys = lambda sys : '_'.join(sys.split('_')[:-1])

templates_names = {
    '3Jets' : {
        'bkg' : (load('%s/templates_lj_3Jets_bkg_%s_QCD_Est_%s.coffea' % (input_dir, jobid, args.year)), load('%s/templates_lj_3Jets_bkg_smoothed_%s_QCD_Est_%s.coffea' % (input_dir, jobid, args.year))),
        #'sig' : (load('%s/templates_lj_3Jets_sig_%s_QCD_Est_%s.coffea' % (input_dir, jobid, args.year))),
    },
    '4PJets' : {
        'bkg' : (load('%s/templates_lj_4PJets_bkg_%s_QCD_Est_%s.coffea' % (input_dir, jobid, args.year)), load('%s/templates_lj_4PJets_bkg_smoothed_%s_QCD_Est_%s.coffea' % (input_dir, jobid, args.year))),
        #'sig' : (load('%s/templates_lj_4PJets_sig_%s_QCD_Est_%s.coffea' % (input_dir, jobid, args.year))),
    },
}

data_lumi_year = prettyjson.loads(open('%s/inputs/lumis_data.json' % proj_dir).read())[args.year]

nom_styles = {'color':'k', 'linestyle':'-', 'label':'Nominal'}
orig_styles = {'color':'b', 'linestyle':'-', 'label':'Original'}
smooth_styles = {'color':'r', 'linestyle':'-', 'label':'Smoothed'}


    ## make plots for background templates
for jmult in templates_names.keys():
    orig_dict, smoothed_dict = templates_names[jmult]['bkg'][0][args.lepton], templates_names[jmult]['bkg'][1][args.lepton]

        # get all keys from both files to make sure they're the same    
    orig_keys = sorted(orig_dict.keys())
    smoothed_keys = sorted(smoothed_dict.keys())
    diff = list(set(orig_keys) - set(smoothed_keys))
    if diff:
        raise ValueError("Input templates for smoothed and original distributions not the same for %s" % jmult)

    topologies = sorted(set([key.split('_')[0] for key in orig_keys if not key == 'data_obs']))
    systs = sorted(set(['_'.join(key.split('_')[1:]) for key in orig_keys if not (key == 'data_obs' or len(key.split('_')) == 1)]))

    #set_trace()
    if args.indiv:
        #for sys in ['hdampDOWN', 'hdampUP', 'mtop1695', 'mtop1755', 'mtopDOWN', 'mtopUP', 'ueDOWN', 'ueUP']:
        #for sys in ['mtop1695', 'mtop1755']:
        for sys in systs:
            #if (('ue' not in sys) or ('hdamp' not in sys) or ('mtop' not in sys)): continue
            for topo in topologies:
                if (sys in systematics.ttJets_sys.values()) and not (topo == 'TT'): continue
                print(jmult, sys, topo)
                pltdir = '/'.join([outdir, jmult, 'Individual', sys])
                if not os.path.isdir(pltdir):
                    os.makedirs(pltdir)
    
                nominal = orig_dict[topo]
                orig_sys = orig_dict['%s_%s' % (topo, sys)]
                smooth_sys = smoothed_dict['%s_%s' % (topo, sys)]
                #set_trace()
    
                x_lims = (0, nominal.dense_axes()[0].centers().size)
    
                fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)            
                fig.subplots_adjust(hspace=.07)
    
                    ## plot normal hists
                ax = hep.plot.histplot(nominal.values()[()], nominal.dense_axes()[0].edges(), ax=ax, histtype='step', **nom_styles) # nosys template
                ax = hep.plot.histplot(orig_sys.values()[()], nominal.dense_axes()[0].edges(), ax=ax, histtype='step', **orig_styles) # original template
                if not np.array_equal(smooth_sys.values()[()], orig_sys.values()[()]):
                    ax = hep.plot.histplot(smooth_sys.values()[()], nominal.dense_axes()[0].edges(), ax=ax, histtype='step', **smooth_styles) # smoothed template
                ax.legend(loc='upper right', title='%s, %s' % (sys, topo))
                ax.set_ylabel('Events')
    
    
                    ## plot relative deviation 
                        ## find first and last valid bins for original variation
                orig_first_valid_bin, orig_last_valid_bin = np.where(~np.isnan((orig_sys.values()[()]-nominal.values()[()])/nominal.values()[()]))[0][0], np.where(~np.isnan((orig_sys.values()[()]-nominal.values()[()])/nominal.values()[()]))[0][-1]+1
                orig_masked_vals = np.ma.masked_where(np.isnan(((orig_sys.values()[()]-nominal.values()[()])/nominal.values()[()])[orig_first_valid_bin:orig_last_valid_bin]), ((orig_sys.values()[()]-nominal.values()[()])/nominal.values()[()])[orig_first_valid_bin:orig_last_valid_bin])
                        ## plot original and smoothed variations
                            ## bins = edges     vals = np.r_[vals, vals[-1]]
                rax.step(nominal.dense_axes()[0].edges()[orig_first_valid_bin:orig_last_valid_bin+1], np.r_[orig_masked_vals, orig_masked_vals[-1]], where='post', **orig_styles)

                if not np.array_equal(smooth_sys.values()[()], orig_sys.values()[()]):
                        ## find first and last valid bins for smoothed variation
                    smooth_first_valid_bin, smooth_last_valid_bin = np.where(~np.isnan((smooth_sys.values()[()]-nominal.values()[()])/nominal.values()[()]))[0][0], np.where(~np.isnan((smooth_sys.values()[()]-nominal.values()[()])/nominal.values()[()]))[0][-1]+1
                    smooth_masked_vals = np.ma.masked_where(np.isnan(((smooth_sys.values()[()]-nominal.values()[()])/nominal.values()[()])[smooth_first_valid_bin:smooth_last_valid_bin]), ((smooth_sys.values()[()]-nominal.values()[()])/nominal.values()[()])[smooth_first_valid_bin:smooth_last_valid_bin])
                        ## plot original and smoothed variations
                            ## bins = edges     vals = np.r_[vals, vals[-1]]
                    rax.step(nominal.dense_axes()[0].edges()[smooth_first_valid_bin:smooth_last_valid_bin+1], np.r_[smooth_masked_vals, smooth_masked_vals[-1]], where='post', **smooth_styles)
                
                rax.axhline(0, **{'linestyle': '--', 'color': (0, 0, 0, 0.5), 'linewidth': 1})
                rax.autoscale(axis='x', tight=True)
                rax.set_xlim(x_lims)
                rax.set_xlabel('m($t\\bar{t}$) $\otimes$ |cos($\\theta^{*}_{t_{l}}$)|')
                rax.set_ylabel('Relative Deviaton from Nominal')
                #    # find y limits using buffer of 1.5 times max deviation from variations
                #max_rat_val = max(max(abs((up_mc.values()[()]/nom_mc.values()[()])[np.isfinite(up_mc.values()[()]/nom_mc.values()[()])]-1)), max(abs((dw_mc.values()[()]/nom_mc.values()[()])[np.isfinite(dw_mc.values()[()]/nom_mc.values()[()])]-1)))
                #max_buffer = round(max_rat_val*1.5, 2)
                #ax.set_ylim(1-max_buffer, 1+max_buffer)
                
                    # add lepton/jet multiplicity label
                ax.text(
                    0.02, 0.94, "%s, %s" % (leptypes[args.lepton], jet_mults[jmult]),
                    fontsize=rcParams['font.size']*0.9, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
                )
                    ## draw vertical lines for distinguishing different ctstar bins
                vlines = [x_lims[1]*ybin/5 for ybin in range(1, 5)]
                for vline in vlines:
                    ax.axvline(vline, color='k', linestyle='--')
                    rax.axvline(vline, color='k', linestyle='--')
                hep.cms.cmslabel(ax=ax, data=False, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
                
                #set_trace()
                figname = '%s/%s' % (pltdir, '_'.join([jmult, args.lepton, sys, topo, 'SmoothedSys_Comp_Scaled'])) if ( (sys == 'mtop1695') or (sys == 'mtop1755')) else '%s/%s' % (pltdir, '_'.join([jmult, args.lepton, sys, topo, 'SmoothedSys_Comp']))
                #figname = '%s/%s' % (pltdir, '_'.join([jmult, args.lepton, sys, topo, 'SmoothedSys_Comp_Corrected_Scaled'])) # only for 2016
                #figname = '%s/%s' % (pltdir, '_'.join([jmult, args.lepton, sys, topo, 'SmoothedSys_Comp_Corrected'])) # only for 2016
                #figname = '%s/%s' % (pltdir, '_'.join([jmult, args.lepton, sys, topo, 'SmoothedSys_Comp']))
                fig.savefig(figname)
                print('%s written' % figname)
                plt.close()
    
    
    if args.comp:
        systypes = sorted(set([baseSys(systematics.sys_to_name[args.year][sys]) for sys in systs])) 
        #set_trace()
        #for sys in ['HDAMP', 'MTOP', 'MTOP3GeV', 'UE']:
        #for sys in ['MTOP3GeV']:
        for sys in systypes:
            up_sys = '%s_UP' % sys
            dw_sys = '%s_DW' % sys
            for topo in topologies:
                if (('%s_UP' % sys in systematics.ttJets_sys.keys()) or ('%s_DW' % sys in systematics.ttJets_sys.keys())) and not (topo == 'TT'): continue
                #if (sys == 'MTOP3GeV') and not (topo == 'TT'): continue
                print(jmult, sys, topo)
                pltdir = '/'.join([outdir, jmult, 'Comp', sys])
                if not os.path.isdir(pltdir):
                    os.makedirs(pltdir)

                up_sysname = [key for key, val in systematics.sys_to_name[args.year].items() if val == up_sys][0]
                dw_sysname = [key for key, val in systematics.sys_to_name[args.year].items() if val == dw_sys][0]

                nominal = orig_dict[topo]
                orig_up = orig_dict['%s_%s' % (topo, up_sysname)]
                orig_dw = orig_dict['%s_%s' % (topo, dw_sysname)]
                smooth_up = smoothed_dict['%s_%s' % (topo, up_sysname)]
                smooth_dw = smoothed_dict['%s_%s' % (topo, dw_sysname)]

                up_histos = [(orig_up, {'color': 'r', 'linestyle': '-', 'label': 'Up'}, False)] if np.array_equal(smooth_up.values()[()], orig_up.values()[()]) else [(orig_up, {'color': 'r', 'linestyle': '--', 'label': 'Original Up'}, True), (smooth_up, {'color': 'r', 'linestyle': '-', 'label': 'Up'}, False)]
                dw_histos = [(orig_dw, {'color': 'b', 'linestyle': '-', 'label': 'Down'}, False)] if np.array_equal(smooth_dw.values()[()], orig_dw.values()[()]) else [(orig_dw, {'color': 'b', 'linestyle': '--', 'label': 'Original Down'}, True), (smooth_dw, {'color': 'b', 'linestyle': '-', 'label': 'Down'}, False)]
                #set_trace()
    
                x_lims = (0, nominal.dense_axes()[0].centers().size)
    
                fig, ax = plt.subplots()
                fig.subplots_adjust(hspace=.07)
    
                    ## plot relative deviations
                for up_histo, up_style, use_fill_between in up_histos:
                            ## find first and last valid bins for up variation
                    up_first_valid_bin, up_last_valid_bin = np.where(~np.isnan((up_histo.values()[()]-nominal.values()[()])/nominal.values()[()]))[0][0], np.where(~np.isnan((up_histo.values()[()]-nominal.values()[()])/nominal.values()[()]))[0][-1]+1
                    up_masked_vals = np.ma.masked_where(np.isnan(((up_histo.values()[()]-nominal.values()[()])/nominal.values()[()])[up_first_valid_bin:up_last_valid_bin]), ((up_histo.values()[()]-nominal.values()[()])/nominal.values()[()])[up_first_valid_bin:up_last_valid_bin])
                    ax.fill_between(nominal.dense_axes()[0].edges()[up_first_valid_bin:up_last_valid_bin+1], np.r_[up_masked_vals, up_masked_vals[-1]], facecolor=up_style['color'], step='post', alpha=0.5) if use_fill_between else ax.step(nominal.dense_axes()[0].edges()[up_first_valid_bin:up_last_valid_bin+1], np.r_[up_masked_vals, up_masked_vals[-1]], where='post', **up_style)

                for dw_histo, dw_style, use_fill_between in dw_histos:
                            ## find first and last valid bins for down variation
                    dw_first_valid_bin, dw_last_valid_bin = np.where(~np.isnan((dw_histo.values()[()]-nominal.values()[()])/nominal.values()[()]))[0][0], np.where(~np.isnan((dw_histo.values()[()]-nominal.values()[()])/nominal.values()[()]))[0][-1]+1
                    dw_masked_vals = np.ma.masked_where(np.isnan(((dw_histo.values()[()]-nominal.values()[()])/nominal.values()[()])[dw_first_valid_bin:dw_last_valid_bin]), ((dw_histo.values()[()]-nominal.values()[()])/nominal.values()[()])[dw_first_valid_bin:dw_last_valid_bin])
                    ax.fill_between(nominal.dense_axes()[0].edges()[dw_first_valid_bin:dw_last_valid_bin+1], np.r_[dw_masked_vals, dw_masked_vals[-1]], facecolor=dw_style['color'], step='post', alpha=0.5) if use_fill_between else ax.step(nominal.dense_axes()[0].edges()[dw_first_valid_bin:dw_last_valid_bin+1], np.r_[dw_masked_vals, dw_masked_vals[-1]], where='post', **dw_style)

                
                ax.legend(loc='upper right', title='%s, %s' % (sys, topo))
                ax.axhline(0, **{'linestyle': '--', 'color': (0, 0, 0, 0.5), 'linewidth': 1})
                ax.autoscale(axis='x', tight=True)
                ax.set_xlim(x_lims)
                ax.set_xlabel('m($t\\bar{t}$) $\otimes$ |cos($\\theta^{*}_{t_{l}}$)|')
                ax.set_ylabel('Relative Deviaton from Nominal')
                
                    # add lepton/jet multiplicity label
                ax.text(
                    0.02, 0.94, "%s, %s" % (leptypes[args.lepton], jet_mults[jmult]),
                    fontsize=rcParams['font.size']*0.9, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
                )
                    ## draw vertical lines for distinguishing different ctstar bins
                vlines = [x_lims[1]*ybin/5 for ybin in range(1, 5)]
                for vline in vlines:
                    ax.axvline(vline, color='k', linestyle='--')
                hep.cms.cmslabel(ax=ax, data=False, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
                
                #set_trace()
                figname = '%s/%s' % (pltdir, '_'.join([jmult, args.lepton, sys, topo, 'SysTemplates_Comp_Scaled'])) if sys == 'MTOP3GeV' else '%s/%s' % (pltdir, '_'.join([jmult, args.lepton, sys, topo, 'SysTemplates_Comp']))
                #figname = '%s/%s' % (pltdir, '_'.join([jmult, args.lepton, sys, topo, 'SysTemplates_Comp_Corrected_Scaled'])) # only 2016
                #figname = '%s/%s' % (pltdir, '_'.join([jmult, args.lepton, sys, topo, 'SysTemplates_Comp_Corrected'])) # only 2016
                #figname = '%s/%s' % (pltdir, '_'.join([jmult, args.lepton, sys, topo, 'SysTemplates_Comp']))
                fig.savefig(figname)
                print('%s written' % figname)
                plt.close()
    
    
