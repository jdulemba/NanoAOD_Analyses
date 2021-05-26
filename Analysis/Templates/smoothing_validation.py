#! /bin/env python

from coffea.util import load#, save
from pdb import set_trace
import os
import Utilities.prettyjson as prettyjson
import numpy as np
import fnmatch
import coffea.processor as processor    
import Utilities.systematics as systematics
import statsmodels.nonparametric.smoothers_lowess as sm
LOWESS = sm.lowess
from scipy.stats import chisquare

# matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend('agg')
from matplotlib import rcParams
rcParams['font.size'] = 20
rcParams["savefig.format"] = 'png'
rcParams["savefig.bbox"] = 'tight'
    
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('year', choices=['2016', '2017', '2018'], help='What year is the ntuple from.')
parser.add_argument('--njets', default='all', nargs='?', choices=['3', '4+', 'all'], help='Specify which jet multiplicity to use.')

args = parser.parse_args()


def smoothing(nominal, template, nbinsx, nbinsy, debug=False):
    if debug: set_trace()
        # np array of original bin values
    nom_vals = nominal.values()[()]
    template_vals = template.values()[()]
    xin = np.arange(nbinsx)
    total_array = np.zeros(nbinsx*nbinsy)

        # loop over each bin of cos theta
    for ybin in range(nbinsy):
        yin = (template_vals[ybin*nbinsx:(ybin+1)*nbinsx] - nom_vals[ybin*nbinsx:(ybin+1)*nbinsx])/nom_vals[ybin*nbinsx:(ybin+1)*nbinsx] # relative deviation from nominal
        total_array[ybin*nbinsx:(ybin+1)*nbinsx] = LOWESS(yin, xin, frac=2./3, it=1, return_sorted=False)

    if debug: set_trace()
        # substitute smoothed array into copy of original hist
    smoothed_histo = template.copy()
    for idx in range(nbinsx*nbinsy):
        smoothed_histo.values()[()][idx] = (1+total_array[idx])*nom_vals[idx]

    if debug: set_trace()
    return smoothed_histo

def smooth_array(nom_array, sys_array, nbinsx, nbinsy, debug=False):
    if debug: set_trace()
    xin = np.arange(nbinsx)
    smoothed_array = np.zeros(nbinsx*nbinsy)

        # loop over each bin of cos theta
    for ybin in range(nbinsy):
        yin = (sys_array[ybin*nbinsx:(ybin+1)*nbinsx] - nom_array[ybin*nbinsx:(ybin+1)*nbinsx])/nom_array[ybin*nbinsx:(ybin+1)*nbinsx] # relative deviation from nominal
        smoothed_array[ybin*nbinsx:(ybin+1)*nbinsx] = LOWESS(yin, xin, frac=2./3, it=1, return_sorted=False)

    if debug: set_trace()
    return smoothed_array

def get_rel_deviation_array(nominal, sys):
        # np array of original bin values
    nom_vals = nominal.values()[()]
    sys_vals = sys.values()[()]
    rel_dev_array = (sys_vals - nom_vals)/nom_vals # relative deviation from nominal

    return rel_dev_array




def plot_bkg_templates():
    '''
    Function that writes linearized mtt vs costheta distributions to root file.
    '''
        # define variables to get histogram for background    
    bkg_fnmatch = '%s.coffea' % base_template_name.replace('NJETS', njets_regex).replace('SIG', 'bkg')
    bkg_fnames = fnmatch.filter(os.listdir(inputdir), bkg_fnmatch)
    
    nbinsx, nbinsy = len(linearize_binning[0])-1, len(linearize_binning[1])-1
    nbins = nbinsx*nbinsy
    #ntoys = 10
    ntoys = 1000
    plus_one_sigma_ind, minus_one_sigma_ind = int((ntoys/2)+(ntoys*0.34))-1, int((ntoys/2)-(ntoys*0.34))-1
    plus_two_sigma_ind, minus_two_sigma_ind = int((ntoys/2)+(ntoys*0.475))-1, int((ntoys/2)-(ntoys*0.475))-1

    sys_to_check = ['mtop', 'ue', 'isr', 'fsr', 'hdamp']

    for bkg_file in bkg_fnames:
        hdict = load(os.path.join(inputdir, bkg_file))
        jmult = '3Jets' if '3Jets' in bkg_file else '4PJets'
        for lep in hdict.keys():
            lepdir = 'mujets' if lep == 'Muon' else 'ejets'
            for tname in hdict[lep].keys():
                #set_trace()
                proc = tname.split('_')[0] if not 'data_obs' in tname else 'data_obs'
                sys = '_'.join(tname.split('_')[1:]) if not 'data_obs' in tname else 'nosys'
                    # skip non-systematic variations
                if sys == 'nosys': continue
                    # skip any systematicss that don't have dedicated samples
                if not any([name for name in sys_to_check if name in sys]): continue
                if not sys in sys_to_use.keys():
                    continue
                    # skip templates that aren't supposed to be smoothed
                if not templates_to_smooth[proc]: continue

                print(jmult, lep, sys, proc)
                pltdir = os.path.join(outdir, lep, jmult, sys)
                if not os.path.isdir(pltdir):
                    os.makedirs(pltdir)


                    # get nosys values
                nosys_histo = hdict[lep]['%s_nosys' % proc]
                nosys_vals = nosys_histo.values()[()]
                    # get vals and errors of systematic variation
                sys_histo = hdict[lep][tname]
                sys_histo_vals, sys_histo_sumw2 = sys_histo.values(sumw2=True)[()]
                sys_histo_errs = np.sqrt(sys_histo_sumw2)

                    # get original relative deviation from nominal
                orig_rel_dev = get_rel_deviation_array(nominal=nosys_histo, sys=sys_histo)
                orig_smooth_rel_dev_array = smooth_array(nom_array=nosys_vals, sys_array=sys_histo_vals, nbinsx=nbinsx, nbinsy=nbinsy)

                    # make toys based on Gaussian distribution of mu=bin_val, sigma=bin_error
                toy_arrays = np.zeros((nbins, ntoys))
                for idx in range(nbins):
                    toy_arrays[idx] = np.random.normal(sys_histo_vals[idx], sys_histo_errs[idx], size=ntoys)

                    # get smoothed relative deviation distributions from toys
                smoothed_rel_dev_arrays = np.zeros((ntoys, nbins))
                chi2_pvals = np.zeros((ntoys, 2))
                for idx in range(ntoys):
                    smoothed_rel_dev_arrays[idx] = smooth_array(nom_array=nosys_vals, sys_array=(toy_arrays.T)[idx], nbinsx=nbinsx, nbinsy=nbinsy)
                    chi2_pval = chisquare(f_obs=(smoothed_rel_dev_arrays[idx]+1)*nosys_vals, f_exp=(orig_smooth_rel_dev_array+1)*nosys_vals) # convert to expected yields so inputs are greater than 5
                    chi2_pvals[idx] = np.array([chi2_pval.statistic, chi2_pval.pvalue])


                    ## find 68% and 95% intervals
                plus_one_sigma_smooth_vals, minus_one_sigma_smooth_vals = np.zeros(nbins), np.zeros(nbins)
                plus_two_sigma_smooth_vals, minus_two_sigma_smooth_vals = np.zeros(nbins), np.zeros(nbins)
                for bin in range(nbins):
                    plus_one_sigma_smooth_vals[bin] = np.sort(smoothed_rel_dev_arrays[:,bin])[plus_one_sigma_ind]
                    minus_one_sigma_smooth_vals[bin] = np.sort(smoothed_rel_dev_arrays[:,bin])[minus_one_sigma_ind]
                    plus_two_sigma_smooth_vals[bin] = np.sort(smoothed_rel_dev_arrays[:,bin])[plus_two_sigma_ind]
                    minus_two_sigma_smooth_vals[bin] = np.sort(smoothed_rel_dev_arrays[:,bin])[minus_two_sigma_ind]
                    
                #### make plots
                    ## plot of yields and relative deviation
                fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                fig.subplots_adjust(hspace=.07)

                    ## plot of relative deviation
                fig2, ax2 = plt.subplots()
                fig2.subplots_adjust(hspace=.07)


                    # plot nosys yield                
                ax.step(nosys_histo.dense_axes()[0].edges(), np.r_[nosys_vals, nosys_vals[-1]], where='post', **{'color' : 'k', 'linestyle' : '-', 'label' : 'Nominal'})

                    # plot original sys dist and relative deviation
                ax.fill_between(nosys_histo.dense_axes()[0].edges(), np.r_[sys_histo_vals, nosys_vals[-1]], facecolor='r', step='post', alpha=0.5, label='Original Systematic')
                rax.fill_between(nosys_histo.dense_axes()[0].edges(), np.r_[orig_rel_dev, orig_rel_dev[-1]], facecolor='r', step='post', alpha=0.5)
                        # for second fig
                ax2.fill_between(nosys_histo.dense_axes()[0].edges(), np.r_[orig_rel_dev, orig_rel_dev[-1]], facecolor='r', step='post', alpha=0.5, label='Original Systematic')

                    # plot original smoothed
                ax.step(nosys_histo.dense_axes()[0].edges(), np.r_[(orig_smooth_rel_dev_array+1)*nosys_vals, ((orig_smooth_rel_dev_array+1)*nosys_vals)[-1]], where='post', **{'color': 'r', 'linestyle': '-', 'label': 'Original Smoothing'})
                rax.step(nosys_histo.dense_axes()[0].edges(), np.r_[orig_smooth_rel_dev_array, orig_smooth_rel_dev_array[-1]], where='post', **{'color': 'r', 'linestyle': '-'})
                        # for second fig
                ax2.step(nosys_histo.dense_axes()[0].edges(), np.r_[orig_smooth_rel_dev_array, orig_smooth_rel_dev_array[-1]], where='post', **{'color': 'r', 'linestyle': '-', 'label': 'Original Smoothing'})

                    # plot 68 and 95% intervals for yields 
                ax.fill_between(nosys_histo.dense_axes()[0].edges(), np.r_[(1+minus_one_sigma_smooth_vals)*nosys_vals, ((1+minus_one_sigma_smooth_vals)*nosys_vals)[-1]], np.r_[(1+plus_one_sigma_smooth_vals)*nosys_vals, ((1+plus_one_sigma_smooth_vals)*nosys_vals)[-1]], where=np.r_[plus_one_sigma_smooth_vals, plus_one_sigma_smooth_vals[-1]] > np.r_[minus_one_sigma_smooth_vals, minus_one_sigma_smooth_vals[-1]], step='post', **{'label':'68%', 'facecolor':'#00f847', 'alpha':0.5})
                ax.fill_between(nosys_histo.dense_axes()[0].edges(), np.r_[(1+minus_two_sigma_smooth_vals)*nosys_vals, ((1+minus_two_sigma_smooth_vals)*nosys_vals)[-1]], np.r_[(1+plus_two_sigma_smooth_vals)*nosys_vals, ((1+plus_two_sigma_smooth_vals)*nosys_vals)[-1]], where=np.r_[plus_two_sigma_smooth_vals, plus_two_sigma_smooth_vals[-1]] > np.r_[minus_two_sigma_smooth_vals, minus_two_sigma_smooth_vals[-1]], step='post', **{'label':'95%', 'facecolor':'#fffc4d', 'alpha':0.5})
                    # plot 68 and 95% intervals for relative deviation
                rax.fill_between(nosys_histo.dense_axes()[0].edges(), np.r_[minus_one_sigma_smooth_vals, minus_one_sigma_smooth_vals[-1]], np.r_[plus_one_sigma_smooth_vals, plus_one_sigma_smooth_vals[-1]], where=np.r_[plus_one_sigma_smooth_vals, plus_one_sigma_smooth_vals[-1]] > np.r_[minus_one_sigma_smooth_vals, minus_one_sigma_smooth_vals[-1]], step='post', **{'facecolor':'#00f847', 'alpha':0.5})
                rax.fill_between(nosys_histo.dense_axes()[0].edges(), np.r_[minus_two_sigma_smooth_vals, minus_two_sigma_smooth_vals[-1]], np.r_[plus_two_sigma_smooth_vals, plus_two_sigma_smooth_vals[-1]], where=np.r_[plus_two_sigma_smooth_vals, plus_two_sigma_smooth_vals[-1]] > np.r_[minus_two_sigma_smooth_vals, minus_two_sigma_smooth_vals[-1]], step='post', **{'facecolor':'#fffc4d', 'alpha':0.5})
                        # for second fig
                ax2.fill_between(nosys_histo.dense_axes()[0].edges(), np.r_[minus_one_sigma_smooth_vals, minus_one_sigma_smooth_vals[-1]], np.r_[plus_one_sigma_smooth_vals, plus_one_sigma_smooth_vals[-1]], where=np.r_[plus_one_sigma_smooth_vals, plus_one_sigma_smooth_vals[-1]] > np.r_[minus_one_sigma_smooth_vals, minus_one_sigma_smooth_vals[-1]], step='post', **{'label':'68%', 'facecolor':'#00f847', 'alpha':0.5})
                ax2.fill_between(nosys_histo.dense_axes()[0].edges(), np.r_[minus_two_sigma_smooth_vals, minus_two_sigma_smooth_vals[-1]], np.r_[plus_two_sigma_smooth_vals, plus_two_sigma_smooth_vals[-1]], where=np.r_[plus_two_sigma_smooth_vals, plus_two_sigma_smooth_vals[-1]] > np.r_[minus_two_sigma_smooth_vals, minus_two_sigma_smooth_vals[-1]], step='post', **{'label':'95%', 'facecolor':'#fffc4d', 'alpha':0.5})


                ax.legend(loc='upper right', title='%s, %s' % (sys, proc))
                ax.autoscale(axis='x', tight=True)
                ax.set_xlim(0, nbins)
                ax.set_ylabel('Events')
                ax.set_yscale('log')
                rax.set_xlabel('m($t\\bar{t}$) $\otimes$ |cos($\\theta^{*}_{t_{l}}$)|')
                rax.set_ylabel('Relative Deviaton from Nominal')
                rax.set_xlim(0, nbins)
                    # add lepton/jet multiplicity label
                ax.text(
                    0.02, 0.94, "%s, %s" % (leptypes[lep], jet_mults[jmult]),
                    fontsize=rcParams['font.size']*0.9, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
                )
                    ## draw vertical lines for distinguishing different ctstar bins
                vlines = [nbins*ybin/5 for ybin in range(1, 5)]
                for vline in vlines:
                    ax.axvline(vline, color='k', linestyle='--')
                    rax.axvline(vline, color='k', linestyle='--')
                hep.cms.cmslabel(ax=ax, data=False, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % lep]/1000., 1))

                #fig.savefig('test')
                #set_trace()
                figname = '%s_yields_and_relative_deviation' % '_'.join([jmult, lep, sys, proc])
                fig.savefig(os.path.join(pltdir, figname))
                print('%s written' % os.path.join(pltdir, figname))


                    ## plot of relative deviation
                ax2.legend(loc='upper right', title='%s, %s' % (sys, proc))
                ax2.autoscale(axis='x', tight=True)
                ax2.set_xlim(0, nbins)
                ax2.set_xlabel('m($t\\bar{t}$) $\otimes$ |cos($\\theta^{*}_{t_{l}}$)|')
                ax2.set_ylabel('Relative Deviaton from Nominal')
                    # add lepton/jet multiplicity label
                ax2.text(
                    0.02, 0.94, "%s, %s" % (leptypes[lep], jet_mults[jmult]),
                    fontsize=rcParams['font.size']*0.9, horizontalalignment='left', verticalalignment='bottom', transform=ax2.transAxes
                )
                    ## draw vertical lines for distinguishing different ctstar bins
                vlines = [nbins*ybin/5 for ybin in range(1, 5)]
                for vline in vlines:
                    ax2.axvline(vline, color='k', linestyle='--')
                hep.cms.cmslabel(ax=ax2, data=False, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % lep]/1000., 1))

                #fig2.savefig('test')
                #set_trace()
                figname2 = '%s_relative_deviation' % '_'.join([jmult, lep, sys, proc])
                fig2.savefig(os.path.join(pltdir, figname2))
                print('%s written' % os.path.join(pltdir, figname2))

                #set_trace()




if __name__ == '__main__':
    proj_dir = os.environ['PROJECT_DIR']
    jobid = os.environ['jobid']

    inputdir = os.path.join(proj_dir, 'Templates', 'results', jobid, args.year)
    base_template_name = 'raw_templates_lj_NJETS_SIG_%s_%s' % (jobid, args.year)

    jet_mults = {
        '3Jets' : '3 jets',
        '4PJets' : '4+ jets'
    }
    
    leptypes = {
        'Muon' : '$\\mu$',
        'Electron' : '$e$',
    }

        # get matching pattern based on args.njets
    if args.njets == '3':
        njets_regex = '3Jets'
    elif args.njets == '4+':
        njets_regex = '4PJets'
    else:
        njets_regex = '*'

    outdir = os.path.join(proj_dir, 'Templates', 'plots', 'Smoothing_Checks', jobid, args.year)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    data_lumi_year = prettyjson.loads(open('%s/inputs/lumis_data.json' % proj_dir).read())[args.year]
    
    # get systematics to run
    sys_to_use = systematics.template_sys_to_name[args.year]
    # get systematics to smooth
    templates_to_smooth = {
        'QCD' : False,
        'TT' : True,
        'VV' : False,
        'TTV' : False,
        'WJets' : False,
        'ZJets' : False,
        'sChannel' : True,
        'tChannel' : True,
        'tWChannel' : True,
        'data_obs' : False,
    }
    
    
    linearize_binning = (
        np.array([300.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 625.0, 650.0, 675.0, 700.0, 730.0, 760.0, 800.0, 850.0, 900.0, 1000.0, 1200.0]),
        #np.array([300.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 620.0, 650., 700.0, 750.0, 800.0, 850.0, 900.0, 2000.0]),
        np.array([0.0, 0.4, 0.6, 0.75, 0.9, 1.0])
    )

    print("Plotting background templates")
    plot_bkg_templates()
