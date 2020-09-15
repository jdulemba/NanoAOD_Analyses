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
#parser.add_argument('--ind', action='store_true', help='Plot individual sys variations')
#parser.add_argument('--comp', action='store_true', help='Plot sys variations with nominal mc')
#parser.add_argument('--ratio', action='store_true', help='Make sys variation/nominal mc ratio plots')

args = parser.parse_args()

#if (not args.comp) and (not args.ratio) and (not args.ind):
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

    for sys in systs:
        for topo in topologies:
            if (sys in systematics.ttJets_sys.values()) and not (topo == 'TT'): continue
            print(jmult, sys, topo)
            pltdir = '/'.join([outdir, jmult, sys])
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
            ax = hep.plot.histplot(smooth_sys.values()[()], nominal.dense_axes()[0].edges(), ax=ax, histtype='step', **smooth_styles) # smoothed template
            ax.legend(loc='upper right', title='%s, %s' % (sys, topo))
            ax.set_ylabel('Events')


                ## plot relative deviation 
                    ## find first and last valid bins for up variation
            orig_first_valid_bin, orig_last_valid_bin = np.where(~np.isnan((orig_sys.values()[()]-nominal.values()[()])/nominal.values()[()]))[0][0], np.where(~np.isnan((orig_sys.values()[()]-nominal.values()[()])/nominal.values()[()]))[0][-1]+1
            orig_masked_vals = np.ma.masked_where(np.isnan(((orig_sys.values()[()]-nominal.values()[()])/nominal.values()[()])[orig_first_valid_bin:orig_last_valid_bin]), ((orig_sys.values()[()]-nominal.values()[()])/nominal.values()[()])[orig_first_valid_bin:orig_last_valid_bin])
                ## find first and last valid bins for down variation
            smooth_first_valid_bin, smooth_last_valid_bin = np.where(~np.isnan((smooth_sys.values()[()]-nominal.values()[()])/nominal.values()[()]))[0][0], np.where(~np.isnan((smooth_sys.values()[()]-nominal.values()[()])/nominal.values()[()]))[0][-1]+1
            smooth_masked_vals = np.ma.masked_where(np.isnan(((smooth_sys.values()[()]-nominal.values()[()])/nominal.values()[()])[smooth_first_valid_bin:smooth_last_valid_bin]), ((smooth_sys.values()[()]-nominal.values()[()])/nominal.values()[()])[smooth_first_valid_bin:smooth_last_valid_bin])
                    ## plot original and smoothed variations
                        ## bins = edges     vals = np.r_[vals, vals[-1]]
            rax.step(nominal.dense_axes()[0].edges()[orig_first_valid_bin:orig_last_valid_bin+1], np.r_[orig_masked_vals, orig_masked_vals[-1]], where='post', **orig_styles)
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
            hep.cms.cmslabel(ax=ax, data=False, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
            
            #set_trace()
            figname = '%s/%s' % (pltdir, '_'.join([jmult, args.lepton, sys, topo, 'SmoothedSys_Comp']))
            fig.savefig(figname)
            print('%s written' % figname)
            plt.close()


