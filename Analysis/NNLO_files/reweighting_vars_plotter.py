#!/usr/bin/env python

# matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend('agg')
from matplotlib import rcParams
rcParams['font.size'] = 20
rcParams["savefig.format"] = 'png'
rcParams["savefig.bbox"] = 'tight'

import Utilities.Plotter as Plotter
from coffea.hist import plot
from coffea import hist
from coffea.util import load, save
from pdb import set_trace
import os
import Utilities.prettyjson as prettyjson
import Utilities.plot_tools as plt_tools
import numpy as np
from coffea.lookup_tools.root_converters import convert_histo_root_file
from coffea.lookup_tools.dense_lookup import dense_lookup
from fnmatch import fnmatch

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('--save_ratios', action='store_true', help='Save NNLO/tune ratios for all years')
args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
base_jobid = os.environ['base_jobid']
analyzer = 'reweighting_vars'

f_ext = 'TOT.coffea'
outdir = os.path.join(proj_dir, 'plots', jobid, analyzer)
if not os.path.isdir(outdir):
    os.makedirs(outdir)
    
style_dict = {
    '2016' : ('2016 CP5' if base_jobid == 'ULnanoAOD' else '2016 T4', '#377eb8'), ## blue
    '2017' : ('2017 CP5', '#e42a2c'), ## red
    '2018' : ('2018 CP5', '#4daf4a'), ## green
}

save_dict = {
    '2016': {'thad_pt': {}, 'mtt_vs_thad_ctstar': {}},
    '2017': {'thad_pt': {}, 'mtt_vs_thad_ctstar': {}},
    '2018': {'thad_pt': {}, 'mtt_vs_thad_ctstar': {}},
}

mtt_binning = np.array([250., 420., 520., 620., 800., 1000., 3500.])
mtt_binlabels = ['%s $\leq$ m($t\\bar{t}$) $\leq$ %s' % (mtt_binning[bin], mtt_binning[bin+1]) for bin in range(len(mtt_binning)-1)]
mtt_bin_locs = np.linspace(3, 33, 6)
ctstar_binning = np.array([-1., -0.65, -0.3, 0., 0.3, 0.65, 1.])
ctstar_binlabels = ['%s $\leq$ cos($\\theta^{*}_{t_{h}}$) $\leq$ %s' % (ctstar_binning[bin], ctstar_binning[bin+1]) for bin in range(len(ctstar_binning)-1)]*len(mtt_binlabels)

    ## get values from NNLO root file
nnlo_fname = 'matrixhists_NNPDF.root' # 'cen' dist has only statistical uncs
#nnlo_fname = 'MATRIX_17_abs.root' # has scale and statistical uncs
png_ext = 'StatUncs' if nnlo_fname == 'matrixhists_NNPDF.root' else 'Scale_StatUncs'
nnlo_leg = '(Stat.)' if nnlo_fname == 'matrixhists_NNPDF.root' else '(Scale+Stat.)'
variables = [
    ('thad_pt', 'thadpt_cen' if nnlo_fname == 'matrixhists_NNPDF.root' else 'thadpth_xs', '$p_{T}$($t_{h}$) [GeV]', '$\dfrac{d \\sigma}{d p_{T}(t_{h})}$', False, None),
    ('mtt_vs_thad_ctstar', 'ttm+cts_cen' if nnlo_fname == 'matrixhists_NNPDF.root' else 'ttm+cts_xs', 'm($t\\bar{t}$) $\otimes$ cos($\\theta^{*}_{t_{h}}$)', '$\dfrac{d^{2} \\sigma}{d m(t\\bar{t}) d cos(\\theta^{*}_{t_{h}})}$', True, [6*ybin for ybin in range(1, 6)]), # last element is hardcoded from binning
]

nnlo_file = convert_histo_root_file(os.path.join(proj_dir, 'NNLO_files', nnlo_fname))
nnlo_dict = Plotter.root_converters_dict_to_hist(nnlo_file, vars=[val[1] for val in variables],
    sparse_axes_list=[{'name': 'dataset', 'label' : "Event Process", 'fill' : 'nnlo'}]
)
#set_trace()


for tune_var, nnlo_var, xtitle, ytitle, linearize, vlines in variables:
    ## make plots and save distributions for top pt reweighting
    logy_min, logy_max = 0., 0.
    normed_logy_min, normed_logy_max = 0., 0.
    
        # orig
    fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
    fig.subplots_adjust(hspace=.07)
        # orig logy
    fig_logy, (ax_logy, rax_logy) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
    fig_logy.subplots_adjust(hspace=.07)
        # norm
    fig_norm, (ax_norm, rax_norm) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
    fig_norm.subplots_adjust(hspace=.07)
        # norm logy
    fig_norm_logy, (ax_norm_logy, rax_norm_logy) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
    fig_norm_logy.subplots_adjust(hspace=.07)
    
    
        ## get values from NNLO root file
    #set_trace()
    nnlo_binning = nnlo_dict[nnlo_var].axis('xaxis').edges()
    nnlo_bins = hist.Bin('xaxis', 'xaxis', nnlo_binning) # axis needs same name as nnlo hist
    
        # raw histo
    nnlo_histo = nnlo_dict[nnlo_var].integrate('dataset')
            # orig
    ax = plot.plot1d(nnlo_histo,
        ax=ax, clear=False,
        line_opts={'linestyle' : '-', 'color' : 'k'},
    )
    ax = plot.plot1d(nnlo_histo, ## need to plot errorbar separately
        ax=ax, clear=False,
        error_opts={'color': 'k', 'marker' : None},
    )
            # orig logy
    ax_logy = plot.plot1d(nnlo_histo,
        ax=ax_logy, clear=False,
        line_opts={'linestyle' : '-', 'color' : 'k'},
    )
    ax_logy = plot.plot1d(nnlo_histo, ## need to plot errorbar separately
        ax=ax_logy, clear=False,
        error_opts={'color': 'k', 'marker' : None},
    )
       # normalized
    nnlo_normed_histo = nnlo_histo.copy()
    nnlo_normed_histo.scale(1./nnlo_normed_histo.values(overflow='all')[()].sum())
            # norm
    ax_norm = plot.plot1d(nnlo_normed_histo,
        ax=ax_norm, clear=False,
        line_opts={'linestyle' : '-', 'color' : 'k'},
    )
    ax_norm = plot.plot1d(nnlo_normed_histo, ## need to plot errorbar separately
        ax=ax_norm, clear=False,
        error_opts={'color': 'k', 'marker' : None},
    )
            # logy
    ax_norm_logy = plot.plot1d(nnlo_normed_histo,
        ax=ax_norm_logy, clear=False,
        line_opts={'linestyle' : '-', 'color' : 'k'},
    )
    ax_norm_logy = plot.plot1d(nnlo_normed_histo, ## need to plot errorbar separately
        ax=ax_norm_logy, clear=False,
        error_opts={'color': 'k', 'marker' : None},
    )

    # update max/min values
    logy_min, logy_max = np.min(nnlo_histo.values()[()]), np.max(nnlo_histo.values()[()]) 
    normed_logy_min, normed_logy_max = np.min(nnlo_normed_histo.values()[()]), np.max(nnlo_normed_histo.values()[()]) 
    
    years_to_run = ['2016APV', '2016', '2017', '2018'] if base_jobid == 'ULnanoAOD' else ['2016', '2017', '2018']
    for year in years_to_run:
        input_dir = os.path.join(proj_dir, 'results', '%s_%s' % (year, jobid), analyzer)
        fnames = sorted(['%s/%s' % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)])
        hdict = plt_tools.add_coffea_files(fnames) if len(fnames) > 1 else load(fnames[0])
    
            ## get NLO values from hdict
        ttSL = 'ttJets_PS' if ((year == '2016') and (base_jobid == 'NanoAODv6')) else 'ttJetsSL'
        xsec_file = prettyjson.loads(open(os.path.join(proj_dir, 'inputs', 'samples_%s.json' % year if base_jobid == 'NanoAODv6' else 'samples_%s_%s.json' % (year, base_jobid))).read())
        tt_dataset = list(filter(lambda x: fnmatch(x['name'], ttSL), xsec_file))[0]
        xsec = tt_dataset['xsection']
        meta_json = prettyjson.loads(open(os.path.join(proj_dir, 'inputs', '%s_%s' % (year, base_jobid), '%s.meta.json' % ttSL)).read())
        sumGenWeights = meta_json["sumGenWeights"]
        
            # orig
        tune_histo = hdict[tune_var].integrate('dataset')
        tune_histo.scale(xsec/sumGenWeights)
        if linearize:
            #set_trace()
            tune_histo = tune_histo.rebin(tune_histo.dense_axes()[0].name, hist.Bin(tune_histo.dense_axes()[0].name, tune_histo.dense_axes()[0].name, mtt_binning))
            tune_histo = tune_histo.rebin(tune_histo.dense_axes()[1].name, hist.Bin(tune_histo.dense_axes()[1].name, tune_histo.dense_axes()[1].name, ctstar_binning))

            # save integral to make normalized hist
        tune_integral = tune_histo.values(overflow='all')[()].sum()
        if linearize:
            tune_histo = Plotter.linearize_hist(tune_histo, no_transpose=True)

            # make normalized hist
        tune_normed_histo = tune_histo.copy()
        tune_normed_histo.scale(1./tune_integral)
        #tune_normed_histo.scale(1./tune_normed_histo.values(overflow='all')[()].sum())
        # rebin
        xaxis_name = tune_histo.dense_axes()[0].name
        tune_histo = tune_histo.rebin(xaxis_name, nnlo_bins)
        tune_normed_histo = tune_normed_histo.rebin(xaxis_name, nnlo_bins)
    
            # plot yields    
        ax = plot.plot1d(tune_histo,
            ax=ax, clear=False,
            line_opts={'linestyle' : '-', 'color' : style_dict[year][1]},
        )
        ax = plot.plot1d(tune_histo, ## need to plot errorbar separately
            ax=ax, clear=False,
            error_opts={'color': style_dict[year][1], 'marker' : None},
        )
                # logy
        ax_logy = plot.plot1d(tune_histo,
            ax=ax_logy, clear=False,
            line_opts={'linestyle' : '-', 'color' : style_dict[year][1]},
        )
        ax_logy = plot.plot1d(tune_histo, ## need to plot errorbar separately
            ax=ax_logy, clear=False,
            error_opts={'color': style_dict[year][1], 'marker' : None},
        )

                ## plot ratios
        first_valid_bin, last_valid_bin = np.where(~np.isnan(nnlo_histo.values()[()]/tune_histo.values()[()]))[0][0], np.where(~np.isnan(nnlo_histo.values()[()]/tune_histo.values()[()]))[0][-1]+1
        ratio_masked_vals = np.ma.masked_where(np.isnan((nnlo_histo.values()[()]/tune_histo.values()[()])[first_valid_bin:last_valid_bin]), (nnlo_histo.values()[()]/tune_histo.values()[()])[first_valid_bin:last_valid_bin])
                    # orig
        rax.step(nnlo_histo.dense_axes()[0].edges()[first_valid_bin:last_valid_bin+1], np.r_[ratio_masked_vals, ratio_masked_vals[-1]], where='post', **{'linestyle' : '-', 'color' : style_dict[year][1]})
                    # logy
        rax_logy.step(nnlo_histo.dense_axes()[0].edges()[first_valid_bin:last_valid_bin+1], np.r_[ratio_masked_vals, ratio_masked_vals[-1]], where='post', **{'linestyle' : '-', 'color' : style_dict[year][1]})
    
            # plot normalized
        tune_normed_histo = tune_histo.copy()
        tune_normed_histo.scale(1./tune_normed_histo.values()[()].sum())
                    # orig
        ax_norm = plot.plot1d(tune_normed_histo,
            ax=ax_norm, clear=False,
            line_opts={'linestyle' : '-', 'color' : style_dict[year][1]},
        )
        ax_norm = plot.plot1d(tune_normed_histo, ## need to plot errorbar separately
            ax=ax_norm, clear=False,
            error_opts={'color': style_dict[year][1], 'marker' : None},
        )
                    # logy
        ax_norm_logy = plot.plot1d(tune_normed_histo,
            ax=ax_norm_logy, clear=False,
            line_opts={'linestyle' : '-', 'color' : style_dict[year][1]},
        )
        ax_norm_logy = plot.plot1d(tune_normed_histo, ## need to plot errorbar separately
            ax=ax_norm_logy, clear=False,
            error_opts={'color': style_dict[year][1], 'marker' : None},
        )
                ## plot ratios
        first_valid_bin_normed, last_valid_bin_normed = np.where(~np.isnan(nnlo_normed_histo.values()[()]/tune_normed_histo.values()[()]))[0][0], np.where(~np.isnan(nnlo_normed_histo.values()[()]/tune_normed_histo.values()[()]))[0][-1]+1
        ratio_masked_vals_normed = np.ma.masked_where(np.isnan((nnlo_normed_histo.values()[()]/tune_normed_histo.values()[()])[first_valid_bin_normed:last_valid_bin_normed]), (nnlo_normed_histo.values()[()]/tune_normed_histo.values()[()])[first_valid_bin_normed:last_valid_bin_normed])
                    # orig
        rax_norm.step(nnlo_normed_histo.dense_axes()[0].edges()[first_valid_bin_normed:last_valid_bin_normed+1], np.r_[ratio_masked_vals_normed, ratio_masked_vals_normed[-1]], where='post', **{'linestyle' : '-', 'color' : style_dict[year][1]})
                    # logy
        rax_norm_logy.step(nnlo_normed_histo.dense_axes()[0].edges()[first_valid_bin_normed:last_valid_bin_normed+1], np.r_[ratio_masked_vals_normed, ratio_masked_vals_normed[-1]], where='post', **{'linestyle' : '-', 'color' : style_dict[year][1]})
    
            ## update min/max values
        if np.min(tune_histo.values()[()]) < logy_min: logy_min = np.min(tune_histo.values()[()])
        if np.max(tune_histo.values()[()]) > logy_max: logy_max = np.max(tune_histo.values()[()])
        if np.min(tune_normed_histo.values()[()]) < normed_logy_min: normed_logy_min = np.min(tune_normed_histo.values()[()])
        if np.max(tune_normed_histo.values()[()]) > normed_logy_max: normed_logy_max = np.max(tune_normed_histo.values()[()])

        #set_trace()
        if args.save_ratios:
            save_dict[year][tune_var] = dense_lookup(ratio_masked_vals_normed, (nnlo_binning)) if vlines is None else dense_lookup(ratio_masked_vals_normed.reshape(6,6), (mtt_binning, ctstar_binning))


    # plot yields
        # format axes
    ax.autoscale(axis='x', tight=True)
    ax.set_ylim(None)
    ax.set_xlabel(None)
    ax.set_ylabel('%s [pb/GeV]' % ytitle)
    
    rax.axhline(1, **{'linestyle': '--', 'color': (0, 0, 0, 0.5), 'linewidth': 1})
    rax.set_ylabel('NNLO/Tune')
    if vlines is None: rax.set_xlabel(xtitle)
    
    # set plotting styles
       ## set legend and corresponding colors
    handles, labels = ax.get_legend_handles_labels()
    new_handles, new_labels = [], []
    for idx, handle in enumerate(handles):
        if isinstance(handle, mpl.container.ErrorbarContainer): #set_trace()
            continue
        hstyle = [val for val in style_dict.values() if handle.get_c() in val[1]]
        if not hstyle:
            labels[idx] = 'NNLO %s' % nnlo_leg
        else:
            labels[idx] = hstyle[0][0]
        new_handles.append(handles[idx])
        new_labels.append(labels[idx])
    # call ax.legend() with the new values
    ax.legend(new_handles, new_labels, loc='upper right')

    ax.text(
        0.02, 0.86,
        "$t\\bart \\rightarrow e/\mu + jets$\nparton level",
        fontsize=rcParams['font.size'], horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
    )

    if vlines is not None:
        # plot ctstar labels
        rax.set_xticks(np.arange(len(ctstar_binlabels)))
        rax.set_xticklabels(ctstar_binlabels)
        ax.tick_params(which='minor', bottom=False, top=False)
        rax.tick_params(which='minor', bottom=False, top=False)

        plt.setp(rax.get_xticklabels(), rotation=90, ha="left", fontsize=rcParams['font.size']*0.5)

        # plot vertical lines for mtt vs ctstar
        for vline in vlines:
                # orig
            ax.axvline(vline, color='k', linestyle='--')
            rax.axvline(vline, color='k', linestyle='--')
                # orig logy
            ax_logy.axvline(vline, color='k', linestyle='--')
            rax_logy.axvline(vline, color='k', linestyle='--')
                # orig norm
            ax_norm.axvline(vline, color='k', linestyle='--')
            rax_norm.axvline(vline, color='k', linestyle='--')
                # orig norm logy
            ax_norm_logy.axvline(vline, color='k', linestyle='--')
            rax_norm_logy.axvline(vline, color='k', linestyle='--')
    
        # plot mtt labels
        for idx, label in enumerate(mtt_binlabels):
            rax.annotate(label, xy=(mtt_bin_locs[idx], 0), xycoords=('data', 'axes fraction'),
                xytext=(0, -120), textcoords='offset points', va='top', ha='center', fontsize=rcParams['font.size']*0.5, rotation=45)
    
    hep.cms.cmslabel(ax=ax, year='')

    #set_trace()
    figname = os.path.join(outdir, '%s_ttJets_xsec_%s_Comp_%s' % (jobid, tune_var, png_ext))
    fig.savefig(figname)
    #set_trace()
    print('%s written' % figname)
    
    
    # plot orig logy
        # format axes
    ax_logy.autoscale(axis='x', tight=True)
    ax_logy.set_xlabel(None)
    ax_logy.set_ylabel('%s [pb/GeV]' % ytitle)
    ax_logy.set_yscale('log')
    ax_logy.set_ylim(logy_min, 15.*logy_max)
    
    rax_logy.axhline(1, **{'linestyle': '--', 'color': (0, 0, 0, 0.5), 'linewidth': 1})
    rax_logy.set_ylabel('NNLO/Tune')
    if vlines is None: rax_logy.set_xlabel(xtitle)
    
    ## set plotting styles
    ax_logy.legend(new_handles, new_labels, loc='upper right')
    
    ax_logy.text(
        0.02, 0.86,
        "$t\\bart \\rightarrow e/\mu + jets$\nparton level",
        fontsize=rcParams['font.size'], horizontalalignment='left', verticalalignment='bottom', transform=ax_logy.transAxes
    )

    if vlines is not None:
        # plot ctstar labels
        rax_logy.set_xticks(np.arange(len(ctstar_binlabels)))
        rax_logy.set_xticklabels(ctstar_binlabels)
        ax_logy.tick_params(which='minor', bottom=False, top=False)
        rax_logy.tick_params(which='minor', bottom=False, top=False)

        plt.setp(rax_logy.get_xticklabels(), rotation=90, ha="left", fontsize=rcParams['font.size']*0.5)

        # plot mtt labels
        for idx, label in enumerate(mtt_binlabels):
            rax_logy.annotate(label, xy=(mtt_bin_locs[idx], 0), xycoords=('data', 'axes fraction'),
                xytext=(0, -120), textcoords='offset points', va='top', ha='center', fontsize=rcParams['font.size']*0.5, rotation=45)
    
    hep.cms.cmslabel(ax=ax_logy, year='')
    
    figname_logy = os.path.join(outdir, '%s_ttJets_xsec_%s_Comp_Logy_%s' % (jobid, tune_var, png_ext))
    fig_logy.savefig(figname_logy)
    print('%s written' % figname_logy)
    
    
    # plot normalized
        # format axes
    ax_norm.autoscale(axis='x', tight=True)
    ax_norm.set_ylim(None)
    ax_norm.set_xlabel(None)
    ax_norm.set_ylabel('$\dfrac{1}{\\sigma}$%s [$GeV^{-1}$]' % ytitle)
    
    rax_norm.axhline(1, **{'linestyle': '--', 'color': (0, 0, 0, 0.5), 'linewidth': 1})
    rax_norm.set_ylabel('NNLO/Tune')
    if vlines is None: rax_norm.set_xlabel(xtitle)
    
    ## set plotting styles
    ax_norm.legend(new_handles, new_labels, loc='upper right')
    
    ax_norm.text(
        0.02, 0.86,
        "$t\\bart \\rightarrow e/\mu + jets$\nparton level",
        fontsize=rcParams['font.size'], horizontalalignment='left', verticalalignment='bottom', transform=ax_norm.transAxes
    )

    if vlines is not None:
        # plot ctstar labels
        rax_norm.set_xticks(np.arange(len(ctstar_binlabels)))
        rax_norm.set_xticklabels(ctstar_binlabels)
        ax_norm.tick_params(which='minor', bottom=False, top=False)
        rax_norm.tick_params(which='minor', bottom=False, top=False)

        plt.setp(rax_norm.get_xticklabels(), rotation=90, ha="left", fontsize=rcParams['font.size']*0.5)

        # plot mtt labels
        for idx, label in enumerate(mtt_binlabels):
            rax_norm.annotate(label, xy=(mtt_bin_locs[idx], 0), xycoords=('data', 'axes fraction'),
                xytext=(0, -120), textcoords='offset points', va='top', ha='center', fontsize=rcParams['font.size']*0.5, rotation=45)
    
    hep.cms.cmslabel(ax=ax_norm, year='')
    
    figname_norm = os.path.join(outdir, '%s_ttJets_xsec_%s_Comp_Norm_%s' % (jobid, tune_var, png_ext))
    fig_norm.savefig(figname_norm)
    print('%s written' % figname_norm)
    
    
    # plot normed logy
        # format axes
    ax_norm_logy.autoscale(axis='x', tight=True)
    ax_norm_logy.set_xlabel(None)
    ax_norm_logy.set_ylabel('$\dfrac{1}{\\sigma}$%s [$GeV^{-1}$]' % ytitle)
    ax_norm_logy.set_yscale('log')
    #set_trace()
    ax_norm_logy.set_ylim(normed_logy_min, 15.*normed_logy_max)
    
    rax_norm_logy.axhline(1, **{'linestyle': '--', 'color': (0, 0, 0, 0.5), 'linewidth': 1})
    rax_norm_logy.set_ylabel('NNLO/Tune')
    if vlines is None: rax_norm_logy.set_xlabel(xtitle)
    
        ## set legend and corresponding colors
    ax_norm_logy.legend(new_handles, new_labels, loc='upper right')
    
    ax_norm_logy.text(
        0.02, 0.86,
        "$t\\bart \\rightarrow e/\mu + jets$\nparton level",
        fontsize=rcParams['font.size'], horizontalalignment='left', verticalalignment='bottom', transform=ax_norm_logy.transAxes
    )

    if vlines is not None:
        # plot ctstar labels
        rax_norm_logy.set_xticks(np.arange(len(ctstar_binlabels)))
        rax_norm_logy.set_xticklabels(ctstar_binlabels)
        ax_norm_logy.tick_params(which='minor', bottom=False, top=False)
        rax_norm_logy.tick_params(which='minor', bottom=False, top=False)

        plt.setp(rax_norm_logy.get_xticklabels(), rotation=90, ha="left", fontsize=rcParams['font.size']*0.5)

        # plot mtt labels
        for idx, label in enumerate(mtt_binlabels):
            rax_norm_logy.annotate(label, xy=(mtt_bin_locs[idx], 0), xycoords=('data', 'axes fraction'),
                xytext=(0, -120), textcoords='offset points', va='top', ha='center', fontsize=rcParams['font.size']*0.5, rotation=45)
    
    hep.cms.cmslabel(ax=ax_norm_logy, year='')
    
    figname_norm_logy = os.path.join(outdir, '%s_ttJets_xsec_%s_Comp_Norm_Logy_%s' % (jobid, tune_var, png_ext))
    fig_norm_logy.savefig(figname_norm_logy)
    print('%s written' % figname_norm_logy)


# save weights
if args.save_ratios:
    ratios_fname = os.path.join(proj_dir, 'NNLO_files', 'NNLO_to_Tune_Ratios_%s.coffea' % jobid)
    save(save_dict, ratios_fname)
    print('\n', ratios_fname, 'written')    
