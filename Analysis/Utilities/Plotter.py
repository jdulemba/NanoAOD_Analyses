from coffea.hist import plot
# matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
#plt.style.use(hep.cms.style.ROOT)
#plt.switch_backend('agg')
from Utilities import styles
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

def plot_stack1d(ax, rax, hdict, xlabel='', ylabel='', sys='nosys', xlimits=None, ylimits=None, **mc_opts):

    mcorder = mc_opts.get('mcorder')

    #set_trace()
    if hdict.sparse_dim() > 1:
        mc_dict = hdict[:, sys].integrate('sys')
        mc_dict = mc_dict[mc_samples]
        data_dict = hdict[:, 'nosys'].integrate('sys')
        data_dict = data_dict[data_samples]
    
    else:
        mc_dict = hdict[mc_samples]
        data_dict = hdict[data_samples]

        ## plot MC and data
    plot.plot1d(mc_dict,
        overlay=mc_dict.axes()[0].name,
        ax=ax,
        clear=False,
        stack=True,
        line_opts=None,
        fill_opts=stack_fill_opts,
        error_opts=stack_error_opts,
        order=mcorder,
    )
    plot.plot1d(data_dict,
        overlay=data_dict.axes()[0].name,
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
    plot.plotratio(data_dict.sum(data_dict.axes()[0].name), mc_dict.sum(mc_dict.axes()[0].name),
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
    mcorder = mc_opts.get('mcorder')
    stack = mc_opts.get('stack', True)

        ## plot MC and data
    plot.plot1d(hdict[mc_samples],
        overlay=hdict.axes()[0].name,
        ax=ax,
        clear=False,
        stack=stack,
        #stack=True,
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

def plot_2d_norm(hdict, values, xlimits, ylimits, xlabel='', ylabel='', mask=None, ax=None, xaxis_name=None, yaxis_name=None, xbins=None, ybins=None, **opts):
    if ax is None:
        ax = plt.gca()

    if xbins is not None:
        xedges = xbins
    else:
        xaxis = hdict.axis(xaxis_name)
        xedges = xaxis.edges()
    if ybins is not None:
        yedges = ybins
    else:
        yaxis = hdict.axis(yaxis_name)
        yedges = yaxis.edges()

    if mask is not None:
        values = np.ma.masked_where(mask, values)

    cmap = opts.get('cmap', 'viridis')
    pc = ax.pcolormesh(xedges, yedges, values.T, cmap=cmap)
    ax.add_collection(pc)

    cmap_label = opts.get('cmap_label', '')
    plt.colorbar(pc, ax=ax, label=cmap_label, pad=0.)

    ax.autoscale(axis='x', tight=True)
    ax.set_xlim(xlimits)
    ax.set_ylim(ylimits)
    xtitle = opts.get('xtitle', xlabel)
    ytitle = opts.get('ytitle', ylabel)
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)
    
    return ax


def plot_1D(values, bins, xlimits, xlabel='', ylabel='Events', linestyle='-', color='k', density=False, weights=None, ax=None, label='', histtype='errorbar', **kwargs):
    if ax is None:
        ax = plt.gca()

    ax = hep.plot.histplot(values, bins, weights=weights, density=density, ax=ax, label=label, linestyle=linestyle, color=color, histtype=histtype)
    ax.set_xlim(xlimits)
    xtitle = kwargs.get('xtitle', xlabel)
    ytitle = kwargs.get('ytitle', ylabel)
    ax.set_xlabel(xtitle)
    ax.set_ylabel(ytitle)

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

    qcd_dict = {
        'ISO' : iso_sb[qcd_samples].integrate('process'),
        'BTAG' : btag_sb[qcd_samples].integrate('process'),
        'DOUBLE' : double_sb[qcd_samples].integrate('process'),
        'SIG' : sig_reg[qcd_samples].integrate('process'),
    }

        # get normalized qcd shape (bins < 0. not allowed)
    qcd_norm_shape = get_qcd_shape(dmp_dict[shape_region])

    # find normalization
    '''
    sideband norm = (QCD MC in sig region)/(QCD MC in norm_region)*(data-prompt in norm_region)

    ABCD norm = (BTAG*ISO)/(DOUBLE) data-prompt yields
    '''
    normalization = 0
    if norm_type == 'ABCD':
        normalization = (np.sum(btag_dmp.values(overflow='all')[()])*np.sum(iso_dmp.values(overflow='all')[()]))/(np.sum(double_dmp.values(overflow='all')[()]))
    elif norm_type == 'Sideband':
        normalization = np.sum(qcd_dict['SIG'].values(overflow='all')[()])/np.sum(qcd_dict[norm_region].values(overflow='all')[()])*np.sum(dmp_dict[norm_region].values(overflow='all')[()])

    qcd_est_array = qcd_norm_shape*normalization
    output_qcd = sig_reg.copy()
        # substitute qcd_est array into original hist
    for idx, val in enumerate(qcd_est_array):
        output_qcd.values(overflow='all')[('QCD',)][idx] = val
    
    return output_qcd

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


def get_samples_yield_and_frac(histo, lumi, promptmc=False, overflow='all', sys='nosys'):
    '''
    Get the yield and relative fraction for each sample of MC, get data yield and compare data/MC
    Returns: list of tuples containing sample name, yield, and fraction
    '''
    if histo.sparse_dim() > 1:
        mc_dict = histo[:, sys].integrate('sys')
        mc_dict = mc_dict[mc_samples]
        data_dict = histo[:, 'nosys'].integrate('sys')
        data_dict = data_dict[data_samples]

    else:
        mc_dict = histo[mc_samples]
        data_dict = histo[data_samples]

    yields = mc_dict.integrate(mc_dict.dense_axes()[0].name).sum().values(overflow=overflow)
    yields.update(data_dict.integrate(data_dict.dense_axes()[0].name).sum().values(overflow=overflow))
    proc_yields_list = [(''.join(process), proc_yields) for process, proc_yields in yields.items()]
    mc_yield = sum([process[1] for process in proc_yields_list if not (process[0] == 'data')])
    data_yield = sum([process[1] for process in proc_yields_list if (process[0] == 'data')])

    rows = [("Lumi: %s fb^-1" % format(lumi, '.1f'), "Sample", "Yield", "Frac")]
    rows += [("", process, format(proc_yield, '.1f'), format((proc_yield/mc_yield)*100, '.1f')) for process, proc_yield in proc_yields_list if not process == 'data']
    rows += [("", "SIM", format(mc_yield, '.1f'), '100.0')]
    rows += [("", "", "", "")]
    rows += [("", "data", format(data_yield, '.1f'), "")]
    rows += [("", "", "", "")]
    rows += [("", "data/SIM", "", format(data_yield/mc_yield, '.3f'))]
    if promptmc:
        prompt_mc_yield = sum([process[1] for process in proc_yields_list if not ((process[0] == 'data') or (process[0] == 'QCD'))])
        rows += [("", "", "", "")]
        rows += [("", "Prompt MC", format(prompt_mc_yield, '.1f'), "")]
        rows += [("", "data-Prompt MC", format(data_yield-prompt_mc_yield, '.1f'), "")]

    yields_dict = {process:round(proc_yield, 2) for process, proc_yield in proc_yields_list}
    yields_dict.update({'SIM': round(mc_yield, 2)})
    yields_dict.update({'data/SIM': round(data_yield/mc_yield, 3)})

    return rows, yields_dict

