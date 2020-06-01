from coffea.hist import plot
# matplotlib
#import matplotlib.pyplot as plt
#import mplhep as hep
#plt.style.use(hep.cms.style.ROOT)
#plt.switch_backend('agg')
import styles
import Utilities.plot_tools as plt_tools
import re
from pdb import set_trace

## make data and mc categories for data/MC plotting
mc_samples = re.compile('(?!data*)')
data_samples = re.compile('(data*)')

hstyles = styles.styles
stack_fill_opts = {'alpha': 0.8, 'edgecolor':(0,0,0,.5)}
stack_error_opts = {'edgecolor':(0,0,0,.5)}

def plot_stack1d(ax, rax, hdict, xlabel='', ylabel='', xlimits=None, ylimits=None, **mc_opts):

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
