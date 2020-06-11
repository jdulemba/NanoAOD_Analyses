from coffea.hist import plot
# matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
#plt.style.use(hep.cms.style.ROOT)
#plt.switch_backend('agg')
import styles
import Utilities.plot_tools as plt_tools
import re
from pdb import set_trace
import numpy as np
#from coffea import hist

## make data and mc categories for data/MC plotting
mc_samples = re.compile('(?!data*)')
data_samples = re.compile('(data*)')
qcd_samples = re.compile('(QCD*)')
prompt_mc_mask = re.compile(r'^(?!.*(\b(?:%s)\b))' % '|'.join(['data*', 'QCD*']))

hstyles = styles.styles
stack_fill_opts = {'alpha': 0.8, 'edgecolor':(0,0,0,.5)}
stack_error_opts = {'edgecolor':(0,0,0,.5)}

def plot_stack1d(ax, rax, hdict, xlabel='', ylabel='', xlimits=None, ylimits=None, **mc_opts):

    mcorder = mc_opts.get('mcorder')

    #set_trace()
        ## plot MC and data
    plot.plot1d(hdict[mc_samples],
        overlay=hdict.axes()[0].name,
        ax=ax,
        clear=False,
        stack=True,
        line_opts=None,
        fill_opts=stack_fill_opts,
        error_opts=stack_error_opts,
        order=mcorder,
    )
    plot.plot1d(hdict[data_samples],
        overlay=hdict.axes()[0].name,
        ax=ax,
        clear=False,
        error_opts=hstyles['data_err_opts']
    )
    ax.autoscale(axis='x', tight=True)
    ax.set_ylim(0, None)
    ax.set_xlabel(None)
    ax.set_xlim(xlimits)
    
        ## set legend and corresponding colors
    handles, labels = ax.get_legend_handles_labels()
    for idx, sample in enumerate(labels):
        if sample == 'data' or sample == 'Observed': continue
        facecolor, legname = plt_tools.get_styles(sample, hstyles)
        handles[idx].set_facecolor(facecolor)
        labels[idx] = legname
    # call ax.legend() with the new values
    ax.legend(handles,labels, loc='upper right')
    #set_trace()
    
        ## plot data/MC ratio
    plot.plotratio(hdict[data_samples].sum(hdict.axes()[0].name), hdict[mc_samples].sum(hdict.axes()[0].name),
        ax=rax,
        error_opts=hstyles['data_err_opts'],
        denom_fill_opts={},
        guide_opts={},
        unc='num'
    )
    rax.set_ylabel('data/MC')
    rax.set_ylim(0.5, 1.5)
    rax.set_xlim(xlimits)
    
        ## set axes labels and titles
    rax.set_xlabel(xlabel)

    return ax, rax



def plot_mc1d(ax, hdict, xlabel='', ylabel='', xlimits=None, ylimits=None, **mc_opts):

    #set_trace()
    mcorder = mc_opts['mcorder'] if 'mcorder' in mc_opts.keys() else None

        ## plot MC and data
    plot.plot1d(hdict[mc_samples],
        overlay=hdict.axes()[0].name,
        ax=ax,
        clear=False,
        stack=True,
        line_opts=None,
        fill_opts=stack_fill_opts,
        error_opts=stack_error_opts,
        order=mcorder,
    )
    ax.autoscale(axis='x', tight=True)
    ax.set_ylim(0, ylimits)
    ax.set_xlabel(xlabel)
    ax.set_xlim(xlimits)
    
        ## set legend and corresponding colors
    handles, labels = ax.get_legend_handles_labels()
    for idx, sample in enumerate(labels):
        if sample == 'data' or sample == 'Observed': continue
        facecolor, legname = plt_tools.get_styles(sample, hstyles)
        handles[idx].set_facecolor(facecolor)
        labels[idx] = legname
    # call ax.legend() with the new values
    ax.legend(handles,labels, loc='upper right')
    
    return ax

def plot_2d_norm(hdict, xaxis_name, yaxis_name, values, xlimits, ylimits, xlabel, ylabel, ax=None, **opts):
    if ax is None:
        ax = plt.gca()

    xaxis = hdict.axis(xaxis_name)
    yaxis = hdict.axis(yaxis_name)
    xedges = xaxis.edges()
    yedges = yaxis.edges()

    cmap=opts['cmap'] if 'cmap' in opts.keys() else 'viridis'
    pc = ax.pcolormesh(xedges, yedges, values.T, cmap=cmap)
    ax.add_collection(pc)

    cmap_label = opts['cmap_label'] if 'cmap_label' in opts.keys() else ''
    plt.colorbar(pc, ax=ax, label=cmap_label, pad=0.)

    ax.autoscale(axis='x', tight=True)
    ax.set_xlim(xlimits)
    ax.set_ylim(ylimits)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    
    return ax


def plot_1D(values, bins, density=False, weights=None, ax=None, label='', histtype='errorbar', **kwargs):
    if ax is None:
        ax = plt.gca()

    hep.plot.histplot(values, bins, weights=weights, density=density, ax=ax, label=label, histtype=histtype, **kwargs)

    return ax


def QCD_Est(sig_reg, iso_sb, btag_sb, double_sb, norm_type=None, shape_region=None, norm_region=None):
    if not norm_type:
        raise ValueError("Normalization type has to be specified for qcd estimation")
    if not shape_region:
        raise ValueError("Region to get qcd shape has to be specified for qcd estimation")

    sig_dmp = data_minus_prompt(sig_reg)
    iso_dmp = data_minus_prompt(iso_sb)
    btag_dmp = data_minus_prompt(btag_sb)
    double_dmp = data_minus_prompt(double_sb)

    dmp_dict = {
        'ISO' : iso_dmp,
        'BTAG' : btag_dmp,
        'DOUBLE' : double_dmp,
    }

        # get normalized qcd shape (bins < 0. not allowed)
    qcd_norm_shape = get_qcd_shape(dmp_dict[shape_region])

    normalization = 0
    if norm_type == 'ABCD':
        normalization = (np.sum(btag_dmp.values(overflow='all')[()])*np.sum(iso_dmp.values(overflow='all')[()]))/(np.sum(double_dmp.values(overflow='all')[()]))

    qcd_est_array = qcd_norm_shape*normalization
        # substitute qcd_est array into original hist
    for idx, val in enumerate(qcd_est_array):
        sig_reg.values(overflow='all')[('QCD',)][idx] = val
    
    return sig_reg

def data_minus_prompt(histo):
    data_hist = histo[data_samples].integrate('process')
    pmc_hist  = histo[prompt_mc_mask].integrate('process')

    neg_pmc_hist = pmc_hist.copy()
    data_minus_prompt = data_hist.copy()
    neg_pmc_hist.scale(-1.)
    data_minus_prompt.add(neg_pmc_hist)

    return data_minus_prompt


def get_qcd_shape(dmp_hist):
    shape_array = dmp_hist.values(overflow='all')[()]
    shape_array[shape_array < 0] = 0
    norm_shape_array = shape_array/np.sum(shape_array)

    return norm_shape_array
