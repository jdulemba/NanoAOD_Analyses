from coffea.hist import plot
# matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
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
import styles
import Utilities.plot_tools as plt_tools
import Utilities.prettyjson as prettyjson
import re
from coffea import hist
import numpy as np
from equal_split import partition_list
import Plotter

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('year', choices=['2016', '2017', '2018'], help='What year is the ntuple from.')
parser.add_argument('lepton', choices=['Electron', 'Muon'], help='Choose which lepton to make plots for')

args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'htt_test_of_homogeneity'

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

lep_cats = {
    'Tight' : 'tight %s' % objtypes['Lep'][args.lepton],
    'Loose' : 'loose %s' % objtypes['Lep'][args.lepton],
}

mtt_binning = np.array([200., 350., 400., 450., 500., 550.,  600., 700., 800., 900., 1200., 1500., 2000.])
#mtt_binning = np.array([200., 360., 380, 400., 420., 440., 460., 480., 500., 520., 540., 560., 580., 610., 640., 680., 730., 800., 900., 2000.])
mtt_bins = hist.Bin('mtt', 'mtt', mtt_binning)
variables = {
    #'mtt_vs_iso' : ('m($t\\bar{t}$) [GeV]', 'pfRelIso, %s' % objtypes['Lep'][args.lepton], mtt_bins, (200., 2000.), 1, (0., 2.) if args.lepton == 'Muon' else (0., 1.5), True),
    'mtt_vs_iso' : ('m($t\\bar{t}$) [GeV]', 'pfRelIso, %s' % objtypes['Lep'][args.lepton], 2, (200., 2000.), 1, (0., 2.) if args.lepton == 'Muon' else (0., 1.5), True),
    'tlep_ctstar_vs_iso' : ('cos($\\theta^{*}_{t_{l}}$)', 'pfRelIso, %s' % objtypes['Lep'][args.lepton], 2, (-1., 1.), 1, (0., 2.) if args.lepton == 'Muon' else (0., 1.5), True),
}



    ## get plotting colors/settings
hstyles = styles.styles
stack_fill_opts = {'alpha': 0.8, 'edgecolor':(0,0,0,.5)}
stack_error_opts = {'edgecolor':(0,0,0,.5)}

    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open('%s/inputs/lumis_data.json' % proj_dir).read())[args.year]
lumi_correction = load('%s/Corrections/%s/MC_LumiWeights_IgnoreSigEvts.coffea' % (proj_dir, jobid))
for hname in hdict.keys():
    if hname == 'cutflow': continue
    hdict[hname].scale(lumi_correction[args.year]['%ss' % args.lepton], axis='dataset')


## make data and mc categories for data/MC plotting
mc_samples = re.compile('(?!data*)')
data_samples = re.compile('(data*)')
prompt_mc_mask = re.compile(r'\b(?:%s)\b' % '|'.join(['ttJets', 'singlet','EWK']))

## make groups based on process
process = hist.Cat("process", "Process", sorting='placement')
process_cat = "dataset"
process_groups = plt_tools.make_dataset_groups(args.lepton, args.year)
#set_trace()
for hname in hdict.keys():
    if hname == 'cutflow': continue
    hdict[hname] = hdict[hname].group(process_cat, process, process_groups)
    

make_plots = True
#make_plots = False
    ## make plots
if make_plots:
    for hname in variables.keys():
        if hname not in hdict.keys():
            raise ValueError("%s not found in file" % hname)
        #set_trace()
        histo = hdict[hname][:, :, args.lepton, :, :].integrate('leptype') # process, jmult, lepton, lepcat, btagregion
    
        if histo.dense_dim() == 2:
            xtitle, ytitle, x_rebinning, x_lims, y_rebinning, y_lims, withData = variables[hname]
            xvar, yvar = hname.split('_vs_')
    
            #set_trace()
            ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
            for jmult in list(set([key[1] for key in histo.values().keys()])):
                for btagregion in list(set([key[3] for key in histo.values().keys()])):
                    for lepcat in list(set([key[2] for key in histo.values().keys()])):
                        pltdir = '/'.join([outdir, args.lepton, jmult, lepcat, btagregion])
                        if not os.path.isdir(pltdir):
                            os.makedirs(pltdir)
    
                        hslice = histo[:, jmult, lepcat, btagregion].integrate('jmult').integrate('lepcat').integrate('btag')
    
                        xaxis_name = hslice.dense_axes()[0].name
                        hslice = hslice.rebin(xaxis_name, x_rebinning)
                        #if x_rebinning != 1:
                        #    xaxis_name = hslice.dense_axes()[0].name
                        #    hslice = hslice.rebin(xaxis_name, x_rebinning)
                        if y_rebinning != 1:
                            yaxis_name = hslice.dense_axes()[1].name
                            hslice = hslice.rebin(yaxis_name, y_rebinning)
    
    
                            # make 1D projection along dense axes
                        for dax in range(2):
                            fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                            fig.subplots_adjust(hspace=.07)
    
                            #set_trace()
                            hproj = hslice.integrate(hslice.dense_axes()[dax].name)
                                ## plot comparison of perm categories
                            plot.plot1d(hproj[mc_samples],
                                overlay=hproj.axes()[0].name,
                                ax=ax,
                                clear=False,
                                stack=True,
                                line_opts=None,
                                fill_opts=stack_fill_opts,
                                error_opts=stack_error_opts,
                                order=['QCD', 'singlet', 'EWK', 'ttJets'],
                            )
                            if withData:
                                plot.plot1d(hproj[data_samples],
                                    overlay=hproj.axes()[0].name,
                                    ax=ax,
                                    clear=False,
                                    error_opts=hstyles['data_err_opts']
                                )
                            ax.autoscale(axis='x', tight=True)
                            ax.set_ylim(0, None)
                            ax.set_ylabel(None)
                            ax.set_xlabel(None)
                            if hproj.dense_axes()[0].name == 'iso':
                                mc_vals = np.sum(np.array(list(hproj[mc_samples].values().values())), axis=0)
                                data_vals = np.sum(np.array(list(hproj[data_samples].values().values())), axis=0)
                                max_binval = np.maximum(np.max(mc_vals), np.max(data_vals))
                                min_binval = np.minimum(np.min(mc_vals), np.min(data_vals))
                                ax.set_yscale('log')
                                ax.set_ylim(np.maximum(1e-1, min_binval), 10**(np.ceil(np.log10(max_binval))+1))
                                if lepcat == 'Tight':
                                    ax.set_xlim((0., 0.2))
                            else:
                                ax.set_xlim(y_lims) if dax == 0 else ax.set_xlim(x_lims)
    
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
                                plot.plotratio(hproj[data_samples].sum(hproj.axes()[0].name), hproj[mc_samples].sum(hproj.axes()[0].name),
                                    ax=rax,
                                    error_opts=hstyles['data_err_opts'],
                                    denom_fill_opts={},
                                    guide_opts={},
                                    unc='num'
                                )
                                rax.set_ylabel('data/MC')
                                rax.set_ylim(0.5, 1.5)
                                if lepcat == 'Tight' and hproj.dense_axes()[0].name == 'iso':
                                    rax.set_xlim((0., 0.2))
                                else:
                                    rax.set_xlim(y_lims) if dax == 0 else rax.set_xlim(x_lims)
    
                            #set_trace()
                                # add lep category
                            ax.text(
                                0.02, 0.91, "%s, %s" % (lep_cats[lepcat], jet_mults[jmult]),
                                fontsize=18,
                                horizontalalignment='left',
                                verticalalignment='bottom',
                                transform=ax.transAxes
                            )
    
                                ## set axes labels and titles
                            plt.xlabel(ytitle) if dax == 0 else plt.xlabel(xtitle)
                            ax = hep.cms.cmslabel(ax=ax, data=withData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
    
                            #set_trace()
                            #figname = 'test.png'
                            figname = '%s/%s.png' % (pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, yvar if dax == 0 else xvar]))
                            fig.savefig(figname, bbox_inches='tight')
                            print('%s written' % figname)
                            plt.close()
                            #set_trace()


#lepIso_invest = False
lepIso_invest = True
if lepIso_invest:
    pltdir = '/'.join([outdir, args.lepton])
    if not os.path.isdir(pltdir):
        os.makedirs(pltdir)

    #set_trace()
    histo = hdict['tlep_ctstar_vs_iso'][:, :, args.lepton, 'Loose', 'btagPass'][data_samples].integrate('process').integrate('leptype').integrate('lepcat').integrate('btag').integrate('ctstar')
    xtitle, ytitle, x_rebinning, x_lims, y_rebinning, y_lims, withData = variables['tlep_ctstar_vs_iso']

    if y_rebinning != 1:
        xaxis_name = histo.dense_axes()[0].name
        histo = histo.rebin(xaxis_name, y_rebinning)

    iso_vals = histo.axis('iso').edges(overflow='all')
    max_iso_ind = 151 if args.lepton == 'Muon' else 101 # splitting doesn't work after a certain cutoff

    vals_to_save = {
        args.lepton : {'3Jets' : {}, '4PJets' : {}, '3PJets' : {}}
    }

    fig, ax = plt.subplots()
    fig.subplots_adjust(hspace=.07)
    for idx, jmult in enumerate(['3Jets', '4PJets', '3PJets']):
        hslice = histo.integrate('jmult') if jmult == '3PJets' else histo[jmult].integrate('jmult')

        ## construct iso sideband regions to investigate based on data in 'Loose' category
        data_array = hslice.values(sumw2=False, overflow='all')[()]
        first_nonzero_bin = np.where(data_array > 0)[0][0]
        #set_trace()
        array_to_partition = np.concatenate((data_array[0:max_iso_ind], np.array([sum(data_array[max_iso_ind:])]))) # add values over max ind as 'overflow'
        data_partitions, inds_for_splitting = partition_list(array_to_partition, 3) # get splitting of data and inds where that splitting occurs: returns list for partitions, ndarray for inds

        iso_region_edge_inds = [0]+inds_for_splitting.tolist()+[len(data_array)-1]
        iso_region_edges = [iso_vals[first_nonzero_bin]]+iso_vals[inds_for_splitting].tolist()+[iso_vals[-2]] # [min iso val, 1st upper, 2nd upper, last non-overflow val]
        iso_region_vals = [(iso_region_edges[idx], iso_region_edges[idx+1]) for idx in range(len(iso_region_edges)-1)]
        vals_to_save[args.lepton][jmult] = {
            'iso_binning' : iso_vals.tolist(),
            'orig_data_array' : data_array.tolist(),
            'partition_loc_inds' : iso_region_edge_inds,
            'iso_region_edges' : iso_region_vals,
            'data_partitions' : data_partitions,
            'data_partition_sums' : list(map(sum, data_partitions)),
        }

        plot.plot1d(hslice,
            ax=ax,
            clear=False,
            error_opts=hstyles[jmult],
        )
    ax.autoscale(axis='x', tight=True)
    #ax.set_ylim(0, None)
    ax.set_xlabel(ytitle)
    ax.set_xlim(y_lims)

        ## set legend labels
    handles, labels = ax.get_legend_handles_labels()
    labels[0] = '3 jets'
    labels[1] = '4+ jets'
    labels[2] = '3+ jets'
    # call ax.legend() with the new values
    ax.legend(handles,labels, loc='upper right')

        ## set axes labels and titles
        # add lepton/jet multiplicity label
    ax.text(
        0.02, 0.91, '%s data\n%s' % (lep_cats['Loose'], '$n_{btags} \geq$ 2'),
        fontsize=18,
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax.transAxes
    )
    ax = hep.cms.cmslabel(ax=ax, data=True, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
    data_vals = np.array(list(histo.values().values()))
    max_binval = np.max(data_vals)
    min_binval = np.min(data_vals[data_vals > 0])
    ax.set_yscale('log')
    ax.set_ylim(np.maximum(1e-1, min_binval), 10**np.ceil(np.log10(max_binval)))

    figname = '%s/%s' % (pltdir, '_'.join(['Loose', args.lepton, 'btagPass', 'data', 'Lep_iso']))
    fig.savefig(figname)
    print('%s written' % figname)
    #set_trace()
    plt.close()
    #set_trace()

    outjson_name = '%s/%s_iso_sideband_regions_construction.json' % (pltdir, '_'.join(['Loose', args.lepton]))
    with open(outjson_name, 'w') as out:
        out.write(prettyjson.dumps(vals_to_save))
    print('%s written' % outjson_name)



    ## make data-prompt plots
    for hname in variables.keys():
        if hname not in hdict.keys():
            raise ValueError("%s not found in file" % hname)
        #set_trace()
        xtitle, ytitle, x_rebinning, x_lims, y_rebinning, y_lims, withData = variables[hname]
        x_rebinning = 5
        xvar, yvar = hname.split('_vs_')

        histo = hdict[hname][:, :, args.lepton, 'Loose', 'btagPass'].integrate('leptype').integrate('lepcat').integrate('btag') # process, jmult, lepton, lepcat, btagregion
        if xvar == 'mtt':
            histo = histo.rebin('mtt', mtt_bins)
        else:
            xaxis_name = histo.dense_axes()[0].name
            histo = histo.rebin(xaxis_name, x_rebinning)
        #if x_rebinning != 1:
            #xaxis_name = histo.dense_axes()[0].name
            #histo = histo.rebin(xaxis_name, x_rebinning)
        #if y_rebinning != 1:
        #    yaxis_name = histo.dense_axes()[1].name
        #    histo = histo.rebin(yaxis_name, y_rebinning)

        iso_val_edges = vals_to_save[args.lepton]['3PJets']['iso_region_edges']
        #iso_val_edges = [(0.15, 0.24), (0.24, 0.43), (0.43, 20.)]
        for jmult in list(set([key[1] for key in histo.values().keys()])):
            pltdir = '/'.join([outdir, args.lepton, jmult, 'Loose', 'btagPass'])
            if not os.path.isdir(pltdir):
                os.makedirs(pltdir)

            fig, ax = plt.subplots()
            fig.subplots_adjust(hspace=.07)
            norm_fig, norm_ax = plt.subplots()
            norm_fig.subplots_adjust(hspace=.07)

            for idx in range(len(iso_val_edges)):
                iso_bin_min, iso_bin_max = iso_val_edges[idx] 
                    ## project mtt and costheta for different values of iso
                hproj = histo[:, jmult, :, iso_bin_min:iso_bin_max].integrate('jmult').integrate('iso') if idx != 2 else histo[:, jmult, :, iso_bin_min:].integrate('jmult').integrate('iso')
                    # get data and prompt mc hists
                data_hist = hproj[data_samples].integrate('process')
                pmc_hist = hproj[prompt_mc_mask].integrate('process')
                    # copy data and prompt mc hists to new ones for creating data-prompt hist without losing originals
                neg_pmc_hist = pmc_hist.copy()
                data_minus_prompt = data_hist.copy()
                neg_pmc_hist.scale(-1.)
                data_minus_prompt.add(neg_pmc_hist)
                #data_minus_prompt.scale(1./np.sum(list(data_minus_prompt.values().values()))) # normalize
    
                set_trace()
                dmp_vals = data_minus_prompt.values()[()]
                dmp_vals[np.where(dmp_vals < 0)] = 0
                #kwargs = {'color' : hstyles['iso_reg_%i' % idx]['color'], 'label' : '%s $\leq$ iso' % (iso_val_edges[idx][0]) if idx == 2 else '%s $\leq$ iso $<$ %s' % (iso_val_edges[idx][0], iso_val_edges[idx][1])}
                norm_vals = dmp_vals/np.sum(dmp_vals)
                bins = data_minus_prompt.axes()[0].edges()
                label = '%s $\leq$ iso' % (iso_val_edges[idx][0]) if idx == 2 else '%s $\leq$ iso $<$ %s' % (iso_val_edges[idx][0], iso_val_edges[idx][1])
                #ax = Plotter.plot_1D(dmp_vals, bins, density=True, ax=ax, label=label, **hstyles['iso_reg_%i' % idx])
                #ax = Plotter.plot_1D(dmp_vals, bins, weights=norm_vals, ax=ax, label=label, **hstyles['iso_reg_%i' % idx])
                norm_ax = Plotter.plot_1D(norm_vals, bins, ax=norm_ax, label=label, **hstyles['iso_reg_%i' % idx])

                plot.plot1d(data_minus_prompt,
                    ax=ax,
                    clear=False,
                    error_opts=hstyles['iso_reg_%i' % idx],
                )
            ax.autoscale(axis='x', tight=True)
            #ax.set_ylabel('A.U.')
            ax.set_xlabel(xtitle)
            ax.set_xlim(x_lims)

                ## set legend labels
            handles, labels = ax.get_legend_handles_labels()
            labels[0] = '%s $\leq$ iso $<$ %s' % (iso_val_edges[0][0], iso_val_edges[0][1])
            labels[1] = '%s $\leq$ iso $<$ %s' % (iso_val_edges[1][0], iso_val_edges[1][1])
            labels[2] = '%s $\leq$ iso' % (iso_val_edges[2][0])
            # call ax.legend() with the new values
            ax.legend(handles,labels, loc='upper right')
        
                ## set axes labels and titles
                # add lepton/jet multiplicity label
            ax.text(
                0.02, 0.86, '%s, %s\n%s\ndata-$prompt_{MC}$' % (lep_cats['Loose'], jet_mults[jmult], '$n_{btags} \geq$ 2'),
                fontsize=18,
                horizontalalignment='left',
                verticalalignment='bottom',
                transform=ax.transAxes
            )
            ax = hep.cms.cmslabel(ax=ax, data=True, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
       
                # norm axes 
            norm_ax.autoscale(axis='x', tight=True)
            norm_ax.set_ylabel('A.U.')
            norm_ax.set_xlabel(xtitle)
            norm_ax.set_xlim(x_lims)
            norm_ax.legend(loc='upper right')
        
                # add lepton/jet multiplicity label
            norm_ax.text(
                0.0, 0.80, '%s, %s\n%s\ndata-$prompt_{MC}$' % (lep_cats['Loose'], jet_mults[jmult], '$n_{btags} \geq$ 2'),
                fontsize=18,
                horizontalalignment='left',
                verticalalignment='bottom',
                transform=ax.transAxes
            )
            norm_ax = hep.cms.cmslabel(ax=norm_ax, data=True, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
        

            #set_trace()
            figname = '%s/%s' % (pltdir, '_'.join([jmult, args.lepton, 'Loose', 'btagPass', 'DataMinusPrompt', xvar]))
            fig.savefig(figname)
            print('%s written' % figname)
            norm_figname = '%s/%s_norm' % (pltdir, '_'.join([jmult, args.lepton, 'Loose', 'btagPass', 'DataMinusPrompt', xvar]))
            norm_fig.savefig(norm_figname)
            print('%s written' % norm_figname)
            #set_trace()
            plt.close()
    
