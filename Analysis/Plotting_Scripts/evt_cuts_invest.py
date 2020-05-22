from coffea.hist import plot
# matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend('agg')
from coffea.util import load
from pdb import set_trace
import os
import styles
import Utilities.plot_tools as plt_tools
import Utilities.prettyjson as prettyjson
import re
from coffea import hist
import numpy as np

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('year', choices=['2016', '2017', '2018'], help='What year is the ntuple from.')
parser.add_argument('lepton', choices=['Electron', 'Muon'], help='Choose which lepton to make plots for')

args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'evt_cuts_invest'

input_dir = '/'.join([proj_dir, 'results', '%s_%s' % (args.year, jobid), analyzer])
f_ext = 'TOT.coffea'
outdir = '/'.join([proj_dir, 'plots', '%s_%s' % (args.year, jobid), analyzer])
if not os.path.isdir(outdir):
    os.makedirs(outdir)

fnames = sorted(['%s/%s' % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)])

#set_trace()
hdict = plt_tools.add_coffea_files(fnames) if len(fnames) > 1 else load(fnames[0])

jet_mults = {
    '3Jets' : '3 jets',
    '4PJets' : '4+ jets'
}

objtypes = {
    'Jets' : 'jets',
    'Lep' :  {
        'Muon' : '$\\mu$',
        'Electron' : '$e$',
    }
}

variables = {
    'Jets_LeadJet_pt' : ('$p_{T}$(leading jet) [GeV]', 2, (0., 300.), True),
    'Jets_njets' : ('$n_{jets}$', 1, (0, 15), True),
    'Lep_iso' : ('pfRelIso, %s' % objtypes['Lep'][args.lepton], 1, (0., 1.), True),
    'MT' : ('$M_{T}$', 1, (0., 300.), True),
    #'MT_Iso' : ('$M_{T}$', 'pfRelIso, %s' % objtypes['Lep'][args.lepton], 1, (0., 300.), 1, (0., 1.), True),
    #'MT_LJpt' : ('$M_{T}$', '$p_{T}$(leading jet) [GeV]', 1, (0., 300.), 2, (0., 300.), True),
    #'Iso_LJpt' : ('pfRelIso, %s' % objtypes['Lep'][args.lepton], '$p_{T}$(leading jet) [GeV]', 1, (0., 1.), 2, (0., 300.), True),
}


    ## get plotting colors/settings
hstyles = styles.styles
stack_fill_opts = {'alpha': 0.8, 'edgecolor':(0,0,0,.5)}
stack_error_opts = {'edgecolor':(0,0,0,.5)}

    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open('%s/inputs/lumis_data.json' % proj_dir).read())[args.year]
lumi_correction = load('%s/Corrections/%s/MC_LumiWeights.coffea' % (proj_dir, jobid))
for hname in hdict.keys():
    if hname == 'cutflow': continue
    hdict[hname].scale(lumi_correction[args.year]['%ss' % args.lepton], axis='dataset')


## make data and mc categories for data/MC plotting
mc_samples = re.compile('(?!data*)')
data_samples = re.compile('(data*)')

## make groups based on process
process = hist.Cat("process", "Process", sorting='placement')
process_cat = "dataset"
process_groups = plt_tools.make_dataset_groups(args.lepton, args.year)
#set_trace()
for hname in hdict.keys():
    if hname == 'cutflow': continue
    hdict[hname] = hdict[hname].group(process_cat, process, process_groups)
    

def get_samples_yield_and_frac(histo, lep):
    '''
    Get the yield and relative fraction for each sample of MC, get data yield and compare data/MC
    Returns: list of tuples containing sample name, yield, and fraction
    '''
    yields = histo.integrate('njets').sum().values()
    proc_yields_list = [(''.join(process), proc_yields) for process, proc_yields in yields.items()]
    mc_yield = sum([process[1] for process in proc_yields_list if not (process[0] == 'data')])
    data_yield = sum([process[1] for process in proc_yields_list if (process[0] == 'data')])

    rows = [("Lumi: %s fb^-1" % format(data_lumi_year['%ss' % lep]/1000., '.1f'), "Sample", "Yield", "Frac")]
    rows += [("", process, format(proc_yield, '.1f'), format((proc_yield/mc_yield)*100, '.1f')) for process, proc_yield in proc_yields_list if not process == 'data']
    rows += [("", "SIM", format(mc_yield, '.1f'), '100.0')]
    rows += [("", "", "", "")]
    rows += [("", "data", format(data_yield, '.1f'), "")]
    rows += [("", "", "", "")]
    rows += [("", "data/SIM", "", format(data_yield/mc_yield, '.3f'))]
        
    yields_dict = {process:round(proc_yield, 2) for process, proc_yield in proc_yields_list}
    yields_dict.update({'SIM': round(mc_yield, 2)})
    yields_dict.update({'data/SIM': round(data_yield/mc_yield, 3)})

    return rows, yields_dict


def cumulative_fraction(data, mc, lep, topos_to_use=[]):
    mc_mask = re.compile(r'\b(?:%s)\b' % '|'.join(topos_to_use))
    data_vals = np.sum(np.array(list(data.values().values())), axis=0)
    all_mc_vals = np.sum(np.array(list(mc.values().values())), axis=0)
    tt_vals = np.sum(np.array(list(mc['ttJets'].values().values())), axis=0)
    mc_to_use_vals = np.sum(np.array(list(mc[mc_mask].values().values())), axis=0)

    edges = mc.dense_axes()[0].edges()
    effs = [0.95, 0.9, 0.85, 0.8, 0.75, 0.7]#, 0.65, 0.6, 0.55, 0.5]
    #set_trace()
    eff_bins = [np.min(np.where(mc_to_use_vals/data_vals < eff)[0]) for eff in effs]
    iso_vals = edges[eff_bins]

    rows = [tuple(['Prompt MC/data: ']+['%.2f' % eff for eff in effs])]
    rows += [tuple(['Iso: ']+['%.2f' % iso for iso in iso_vals])]

    return rows


    ## make plots
for hname in variables.keys():
    if hname not in hdict.keys():
        raise ValueError("%s not found in file" % hname)
    histo = hdict[hname][:, :, args.lepton].integrate('leptype')
    #set_trace()

    if histo.dense_dim() == 1:
        xtitle, rebinning, x_lims, withData = variables[hname]

        ## hists should have 3 category axes (dataset, jet multiplicity, lepton type) followed by variable
        for jmult in list(set([key[1] for key in histo.values().keys()])):
            pltdir = '/'.join([outdir, args.lepton, jmult])
            if not os.path.isdir(pltdir):
                os.makedirs(pltdir)

            if withData:
                fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
            else:
                fig, ax = plt.subplots()
            fig.subplots_adjust(hspace=.07)
            hslice = histo[:, jmult].integrate('jmult')

            #if 'Lep_eta' in hname: set_trace()
            if hname == 'Jets_njets':
                print(jmult)
                yields_txt, yields_json = get_samples_yield_and_frac(hslice, args.lepton)
                plt_tools.print_table(yields_txt, filename='%s/%s_%s_yields_and_fracs.txt' % (pltdir, jmult, args.lepton), print_output=True)
                with open('%s/%s_%s_yields_and_fracs.json' % (pltdir, jmult, args.lepton), 'w') as out:
                    out.write(prettyjson.dumps(yields_json))

            if rebinning != 1:
                xaxis_name = hslice.dense_axes()[0].name
                hslice = hslice.rebin(xaxis_name, rebinning)

            mc_dict = {key:np.sum(list(hslice[key].values().values()), axis=1) for key in [key[0] for key in list(hslice[mc_samples].values().keys())]} # make {topo: sum events} dict for sorting
            mc_order = sorted(mc_dict, key=mc_dict.get, reverse=False)
            #set_trace()
                ## plot MC and data
            plot.plot1d(hslice[mc_samples],
                overlay=hslice.axes()[0].name,
                ax=ax,
                clear=False,
                stack=True,
                line_opts=None,
                fill_opts=stack_fill_opts,
                error_opts=stack_error_opts,
                order=mc_order,
            )
            if withData:
                plot.plot1d(hslice[data_samples],
                    overlay=hslice.axes()[0].name,
                    ax=ax,
                    clear=False,
                    error_opts=hstyles['data_err_opts']
                )
            ax.autoscale(axis='x', tight=True)
            ax.set_ylim(0, None)
            ax.set_xlabel(None)
            ax.set_xlim(x_lims)

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

            if withData:
                    ## plot data/MC ratio
                plot.plotratio(hslice[data_samples].sum(hslice.axes()[0].name), hslice[mc_samples].sum(hslice.axes()[0].name), 
                    ax=rax,
                    error_opts=hstyles['data_err_opts'], 
                    #error_opts=data_err_opts, 
                    denom_fill_opts={},
                    guide_opts={},
                    unc='num'
                )
                rax.set_ylabel('data/MC')
                rax.set_ylim(0.5, 1.5)
                rax.set_xlim(x_lims)


                ## set axes labels and titles
            plt.xlabel(xtitle)
                # add lepton/jet multiplicity label
            ax.text(
                0.02, 0.92, "%s/%s" % (objtypes['Lep'][args.lepton], jet_mults[jmult]),
                fontsize=hep.styles_cms.CMS['font.size']*0.90, 
                horizontalalignment='left', 
                verticalalignment='bottom', 
                transform=ax.transAxes
            )
            ax = hep.cms.cmslabel(ax=ax, data=True if withData else False, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
            if hname == 'Lep_iso':
                print(jmult)
                prompt_mc_frac = cumulative_fraction(data=hslice[data_samples], mc=hslice[mc_samples], lep=args.lepton, topos_to_use=['ttJets', 'singlet','EWK'])
                #set_trace()
                plt_tools.print_table(prompt_mc_frac, filename='%s/%s_%s_Iso_Fracs.txt' % (pltdir, jmult, args.lepton), print_output=True)

                mc_vals = np.sum(np.array(list(hslice[mc_samples].values().values())), axis=0)
                data_vals = np.sum(np.array(list(hslice[data_samples].values().values())), axis=0)
                max_binval = np.maximum(np.max(mc_vals), np.max(data_vals))
                min_binval = np.minimum(np.min(mc_vals), np.min(data_vals))
                ax.set_yscale('log')
                ax.set_ylim(np.maximum(1e-1, min_binval), 10**np.ceil(np.log10(max_binval)))

            figname = '%s/%s.png' % (pltdir, '_'.join([jmult, args.lepton, hname]))
            fig.savefig(figname, bbox_inches='tight')
            print('%s written' % figname)
            plt.close()
            #set_trace()

    elif histo.dense_dim() == 2:
        xtitle, ytitle, x_rebinning, x_lims, y_rebinning, y_lims, withData = variables[hname]

        ## hists should have 3 category axes (dataset, jet multiplicity, lepton type) followed by variable
        for lep in [args.lepton]:
            for jmult in histo.axis('jmult')._sorted:
                pltdir = '/'.join([outdir, lep, jmult])
                if not os.path.isdir(pltdir):
                    os.makedirs(pltdir)

                hslice = histo[:, jmult, lep].integrate('jmult').integrate('leptype')

                if x_rebinning != 1:
                    xaxis_name = hslice.dense_axes()[0].name
                    hslice = hslice.rebin(xaxis_name, x_rebinning)
                if y_rebinning != 1:
                    yaxis_name = hslice.dense_axes()[1].name
                    hslice = hslice.rebin(yaxis_name, y_rebinning)

                    # make 1D projections along y axis
                bin_slices = [(x, x+10) for x in range(0, 100, 10)] # [(0, 10), (10, 20)..]
                for bin_xmin, bin_xmax in bin_slices:
                    if withData:
                        fig, (ax, rax) = plt.subplots(2, 1, figsize=(7,7), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                        fig.subplots_adjust(hspace=.07)
                    else:
                        fig, ax = plt.subplots(1, 1, figsize=(7,7))#, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                        fig.subplots_adjust(hspace=.07)

                    #set_trace()
                    y_slice = hslice[:, bin_xmin:bin_xmax, :].integrate(hslice.dense_axes()[0].name)
                        ## plot MC and data
                    plot.plot1d(y_slice[mc_samples],
                        overlay=y_slice.axes()[0].name,
                        ax=ax,
                        clear=False,
                        stack=True,
                        line_opts=None,
                        fill_opts=stack_fill_opts,
                        error_opts=stack_error_opts
                    )
                    if withData:
                        plot.plot1d(y_slice[data_samples],
                            overlay=y_slice.axes()[0].name,
                            ax=ax,
                            clear=False,
                            error_opts=hstyles['data_err_opts']
                        )
                    ax.autoscale(axis='x', tight=True)
                    ax.set_ylim(0, None)
                    ax.set_xlabel(None)
                    ax.set_xlim(y_lims)

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

                    if withData:
                            ## plot data/MC ratio
                        plot.plotratio(y_slice[data_samples].sum(y_slice.axes()[0].name), y_slice[mc_samples].sum(y_slice.axes()[0].name), 
                            ax=rax,
                            error_opts=hstyles['data_err_opts'], 
                            #error_opts=data_err_opts, 
                            denom_fill_opts={},
                            guide_opts={},
                            unc='num'
                        )
                        rax.set_ylabel('data/MC')
                        rax.set_ylim(0.5, 1.5)
                        rax.set_xlim(x_lims)

                    #    # add x axis range for projection
                    #ax.text(
                    #    0.02, 0.95, "%.1f $\leq$ %s $\leq$ %.1f" % (bin_xmin, xtitle, bin_xmax),
                    #    fontsize=8, 
                    #    horizontalalignment='left', 
                    #    verticalalignment='bottom', 
                    #    transform=ax.transAxes
                    #)

                        ## set axes labels and titles
                    plt.xlabel(ytitle)
                    ax.text(
                        0.02, 0.92, "%s/%s\n%.1f $\leq$ %s $\leq$ %.1f" % (objtypes['Lep'][lep], jet_mults[jmult], bin_xmin, xtitle, bin_xmax),
                        fontsize=hep.styles_cms.CMS['font.size']*0.90, 
                        horizontalalignment='left', 
                        verticalalignment='bottom', 
                        transform=ax.transAxes
                    )
                    ax = hep.cms.cmslabel(ax=ax, data=True if withData else False, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % lep]/1000., 1))

                    #figname = 'test.png'
                    figname = '%s/%s.png' % (pltdir, '_'.join([jmult, lep, hname, 'xrange%sto%s' % (bin_xmin, bin_xmax)]))
                    fig.savefig(figname, bbox_inches='tight')
                    print('%s written' % figname)
                    plt.close()
                    #set_trace()


                    # make 2D plots for each topology
                for top in hslice.sparse_axes()[0]._sorted:
                    if not withData and top == 'data': continue

                    if withData:
                        fig, (ax, rax) = plt.subplots(2, 1, figsize=(7,7), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                        fig.subplots_adjust(hspace=.07)
                    else:
                        fig, ax = plt.subplots(1, 1, figsize=(7,7))#, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                        fig.subplots_adjust(hspace=.07)

                    htop = hslice[top].integrate('process')
                    top_label = plt_tools.get_label(top, hstyles)
                    htop.label = '%s %s' % (top_label, htop.label)
                        ## plot MC and data
                    #set_trace()
                    plot.plot2d(htop,
                        xaxis=htop.axes()[0].name,
                        ax=ax,
                        patch_opts={'cmap' : 'OrRd'},
                        clear=True,
                    )
                    ax.autoscale(axis='x', tight=True)
                    ax.set_ylim(y_lims)
                    ax.set_xlim(x_lims)


                        ## set axes labels and titles
                    plt.xlabel(xtitle)
                    plt.ylabel(ytitle)
                    #kwargs = {
                    #    'Lumi_blurb' : "(13 TeV %.2f fb$^{-1}$, %s/%s)" % (data_lumi_year['%ss' % lep]/1000., objtypes['Lep'][lep], jet_mults[jmult]),
                    #}
                    #ax = plt_tools.make_cms_lumi_blurb(ax, **kwargs)
                    ax.text(
                        0.02, 0.92, "%s/%s" % (objtypes['Lep'][lep], jet_mults[jmult]),
                        fontsize=hep.styles_cms.CMS['font.size']*0.90, 
                        horizontalalignment='left', 
                        verticalalignment='bottom', 
                        transform=ax.transAxes
                    )
                    ax = hep.cms.cmslabel(ax=ax, data=True if withData else False, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % lep]/1000., 1))

                    figname = '%s/%s.png' % (pltdir, '_'.join([jmult, lep, hname, top]))
                    fig.savefig(figname, bbox_inches='tight')
                    print('%s written' % figname)
                    plt.close()
                    #set_trace()

