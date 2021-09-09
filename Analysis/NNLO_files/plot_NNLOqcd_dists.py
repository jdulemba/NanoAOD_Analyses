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
#from coffea.util import load, save
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
parser.add_argument('--save_ratios', action='store_true', help='Save NNLO/tune ratios for all years (default is to not save them).')
args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
base_jobid = os.environ['base_jobid']
analyzer = 'NNLOqcd_dists'

f_ext = 'TOT.coffea'
outdir = os.path.join(proj_dir, 'plots', base_jobid, analyzer)
if not os.path.isdir(outdir):
    os.makedirs(outdir)
    
    ## get values from NNLO root file
nnlo_fname = 'MATRIX_ttmVStheta.root' # 'xsec_central' dist has only statistical uncs
#nnlo_fname = 'matrixhists_NNPDF.root' # 'cen' dist has only statistical uncs
#nnlo_fname = 'MATRIX_17_abs.root' # has scale and statistical uncs
nnlo_file = convert_histo_root_file(os.path.join(proj_dir, 'NNLO_files', nnlo_fname))
nnlo_var = 'xsec_central'
nnlo_dict = Plotter.root_converters_dict_to_hist(nnlo_file, vars=[nnlo_var],
    sparse_axes_list=[{'name': 'dataset', 'label' : "Event Process", 'fill' : 'nnlo'}],
    #dense_axes_list=[{'name': 'mtt', 'idx' : 1}, {'name' : 'ctstar', 'idx' : 0}],
    #transpose_da=True,
    dense_axes_list=[{'name' : 'ctstar', 'idx' : 0}, {'name': 'mtt', 'idx' : 1}],
)

mtt_binning = nnlo_dict[nnlo_var].axis('mtt').edges()
ctstar_binning = nnlo_dict[nnlo_var].axis('ctstar').edges()
vlines = [(len(mtt_binning)-1)*ybin for ybin in range(1, len(ctstar_binning)-1)]
#set_trace()

## linearize nnlo_hist
nnlo_lin_dict = {name: Plotter.linearize_hist(nnlo_dict[name].integrate('dataset')) for name in nnlo_dict.keys()}

png_ext = 'StatUncs'
nnlo_leg = '(Stat.)'

axis_labels_dict = {
    'mtt' : '$m_{t\\bar{t}}$',
    'ctstar' : 'cos($\\theta^{*}_{t_{h}}$)'
}


############## Plot hists ################

### raw histo ###
nnlo_histo = nnlo_dict[nnlo_var]

        # histo plotting params
xtitle, ytitle = axis_labels_dict[nnlo_histo.dense_axes()[0].name], axis_labels_dict[nnlo_histo.dense_axes()[1].name]
ztitle = '$\dfrac{d^{2} \\sigma}{d m_{t\\bar{t}} d cos(\\theta^{*}_{t_{h}})}$'
opts = {'cmap_label' : ztitle}
norm_opts = {'cmap_label' : '$\dfrac{1}{\\sigma}$%s [$GeV^{-1}$]' % ztitle}
x_lims = (np.min(nnlo_histo.dense_axes()[0].edges()), np.max(nnlo_histo.dense_axes()[0].edges()))
y_lims = (np.min(nnlo_histo.dense_axes()[1].edges()), np.max(nnlo_histo.dense_axes()[1].edges()))

    # plot 2d version of original hist
fig, ax = plt.subplots()
fig.subplots_adjust(hspace=.07)

nnlo_vals = nnlo_histo.integrate('dataset').values()[()]
Plotter.plot_2d_norm(nnlo_histo, xaxis_name=nnlo_histo.dense_axes()[0].name, yaxis_name=nnlo_histo.dense_axes()[1].name,
    values=np.ma.masked_where(nnlo_vals <= 0.0, nnlo_vals), # mask nonzero probabilities for plotting
    xlimits=x_lims, ylimits=y_lims, xlabel=xtitle, ylabel=ytitle,
    ax=ax, **opts)

ax.text(
    0.98, 0.90, "$t\\bart \\rightarrow e/\mu + jets$\nparton level",
    fontsize=rcParams['font.size'], horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes, color="w",
    #0.02, 0.86, "$t\\bart \\rightarrow e/\mu + jets$\nparton level",
    #fontsize=rcParams['font.size'], horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
)
hep.label._exp_label(ax=ax, exp="MATRIX", rlabel="NNLO QCD")

figname = os.path.join(outdir, '%s_ttJets_xsec_%s_vs_%s_%s' % (base_jobid, nnlo_histo.dense_axes()[0].name, nnlo_histo.dense_axes()[1].name, png_ext))
fig.savefig(figname)
print('%s written' % figname)
    

    # plot 2d version of transposed hist
fig, ax = plt.subplots()
fig.subplots_adjust(hspace=.07)

Plotter.plot_2d_norm(nnlo_histo, xaxis_name=nnlo_histo.dense_axes()[1].name, yaxis_name=nnlo_histo.dense_axes()[0].name,
    values=np.ma.masked_where(nnlo_vals.T <= 0.0, nnlo_vals.T), # mask nonzero probabilities for plotting
    xlimits=y_lims, ylimits=x_lims, xlabel=ytitle, ylabel=xtitle,
    ax=ax, **opts)

ax.text(
    0.98, 0.90, "$t\\bart \\rightarrow e/\mu + jets$\nparton level",
    fontsize=rcParams['font.size'], horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes, color="w",
)
hep.label._exp_label(ax=ax, exp="MATRIX", rlabel="NNLO QCD")

figname = os.path.join(outdir, '%s_ttJets_xsec_%s_vs_%s_%s' % (base_jobid, nnlo_histo.dense_axes()[1].name, nnlo_histo.dense_axes()[0].name, png_ext))
fig.savefig(figname)
print('%s written' % figname)
    

    # plot normalized 2d version of original hist
fig_norm, ax_norm = plt.subplots()
fig_norm.subplots_adjust(hspace=.07)

nnlo_norm_vals = nnlo_vals/np.sum(nnlo_vals)
Plotter.plot_2d_norm(nnlo_histo, xaxis_name=nnlo_histo.dense_axes()[0].name, yaxis_name=nnlo_histo.dense_axes()[1].name,
    values=np.ma.masked_where(nnlo_norm_vals <= 0.0, nnlo_norm_vals), # mask nonzero probabilities for plotting
    xlimits=x_lims, ylimits=y_lims, xlabel=xtitle, ylabel=ytitle,
    ax=ax_norm, **norm_opts)

ax_norm.text(
    0.98, 0.90, "$t\\bart \\rightarrow e/\mu + jets$\nparton level",
    fontsize=rcParams['font.size'], horizontalalignment='right', verticalalignment='bottom', transform=ax_norm.transAxes, color="w",
)
hep.label._exp_label(ax=ax_norm, exp="MATRIX", rlabel="NNLO QCD")

figname_norm = os.path.join(outdir, '%s_ttJets_xsec_%s_vs_%s_%s_Norm' % (base_jobid, nnlo_histo.dense_axes()[0].name, nnlo_histo.dense_axes()[1].name, png_ext))
fig_norm.savefig(figname_norm)
print('%s written' % figname_norm)
    

    # plot normalized 2d version of transposed hist
fig_norm, ax_norm = plt.subplots()
fig_norm.subplots_adjust(hspace=.07)

Plotter.plot_2d_norm(nnlo_histo, xaxis_name=nnlo_histo.dense_axes()[1].name, yaxis_name=nnlo_histo.dense_axes()[0].name,
    values=np.ma.masked_where(nnlo_norm_vals.T <= 0.0, nnlo_norm_vals.T), # mask nonzero probabilities for plotting
    xlimits=y_lims, ylimits=x_lims, xlabel=ytitle, ylabel=xtitle,
    ax=ax_norm, **norm_opts)

ax_norm.text(
    0.98, 0.90, "$t\\bart \\rightarrow e/\mu + jets$\nparton level",
    fontsize=rcParams['font.size'], horizontalalignment='right', verticalalignment='bottom', transform=ax_norm.transAxes, color="w",
)
hep.label._exp_label(ax=ax_norm, exp="MATRIX", rlabel="NNLO QCD")

figname_norm = os.path.join(outdir, '%s_ttJets_xsec_%s_vs_%s_%s_Norm' % (base_jobid, nnlo_histo.dense_axes()[1].name, nnlo_histo.dense_axes()[0].name, png_ext))
fig_norm.savefig(figname_norm)
print('%s written' % figname_norm)
    


### plot linearized version of hist ###
nnlo_lin_hist = nnlo_lin_dict[nnlo_var]

        # original hist
fig, ax = plt.subplots()
fig.subplots_adjust(hspace=.07)

plot.plot1d(nnlo_lin_hist,
    ax=ax, clear=False,
    line_opts={'linestyle' : '-', 'color' : 'k'},
)
#plot.plot1d(nnlo_lin_hist, ## need to plot errorbar separately
#    ax=ax, clear=False,
#    error_opts={'color': 'k', 'marker' : None},
#)

    # format axes
ax.autoscale()
ax.set_ylim(None, ax.get_ylim()[1]*1.15)
ax.set_xlabel("%s $\otimes$ %s" % (xtitle, ytitle))
ax.set_ylabel('%s [pb/GeV]' % ztitle)
ax.get_legend().remove()

ax.text(
    0.02, 0.90, "$t\\bart \\rightarrow e/\mu + jets$\nparton level",
    fontsize=rcParams['font.size'], horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes,
)
hep.label._exp_label(ax=ax, exp="MATRIX", rlabel="NNLO QCD")

figname_lin = os.path.join(outdir, '%s_ttJets_xsec_%s_otimes_%s_%s' % (base_jobid, nnlo_histo.dense_axes()[0].name, nnlo_histo.dense_axes()[1].name, png_ext))
fig.savefig(figname_lin)
print('%s written' % figname_lin)
    

        # normalized hist
nnlo_lin_normed_hist = nnlo_lin_hist.copy()
nnlo_lin_normed_hist.scale(1./nnlo_lin_normed_hist.values(overflow='all')[()].sum())

fig_norm, ax_norm = plt.subplots()
fig_norm.subplots_adjust(hspace=.07)

plot.plot1d(nnlo_lin_normed_hist,
    ax=ax_norm, clear=False,
    line_opts={'linestyle' : '-', 'color' : 'k'},
)
#plot.plot1d(nnlo_lin_normed_hist, ## need to plot errorbar separately
#    ax=ax_norm, clear=False,
#    error_opts={'color': 'k', 'marker' : None},
#)

    # format axes
ax_norm.autoscale()
ax_norm.set_ylim(None, ax_norm.get_ylim()[1]*1.15)
ax_norm.set_xlabel("%s $\otimes$ %s" % (xtitle, ytitle))
ax_norm.set_ylabel(norm_opts['cmap_label'])
ax_norm.get_legend().remove()

ax_norm.text(
    0.02, 0.90, "$t\\bart \\rightarrow e/\mu + jets$\nparton level",
    fontsize=rcParams['font.size'], horizontalalignment='left', verticalalignment='bottom', transform=ax_norm.transAxes,
)
hep.label._exp_label(ax=ax_norm, exp="MATRIX", rlabel="NNLO QCD")

figname_lin_norm = os.path.join(outdir, '%s_ttJets_xsec_%s_otimes_%s_%s_Norm' % (base_jobid, nnlo_histo.dense_axes()[0].name, nnlo_histo.dense_axes()[1].name, png_ext))
fig_norm.savefig(figname_lin_norm)
print('%s written' % figname_lin_norm)
