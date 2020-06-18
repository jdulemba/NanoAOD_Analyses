# matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend('agg')
from matplotlib import rcParams
rcParams['font.size'] = 20
rcParams["savefig.format"] = 'png'
rcParams["savefig.bbox"] = 'tight'
from coffea.util import load, save
from pdb import set_trace
import os
import styles
import Utilities.plot_tools as plt_tools
import Utilities.prettyjson as prettyjson
import re
from coffea import hist
from coffea.hist import plot
import numpy as np
import fnmatch
import Plotter as Plotter
from equal_split import partition_list
from scipy import stats
from scipy import interpolate
from coffea.lookup_tools.dense_lookup import dense_lookup

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('--year', type=str, help='Choose year(s) to run')

args = parser.parse_args()

years_to_run = [args.year] if args.year else ['2016', '2017', '2018']

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'ttbar_alpha_reco'

blurb = 'tight $e/\mu$+3 jets\n$n_{btags} \geq$ 2'


alpha_corrections = {
    '2016' : {
        'E' : {},
        'P' : {},
    },
    '2017' : {
        'E' : {},
        'P' : {},
    },
    '2018' : {
        'E' : {},
        'P' : {},
    },
}


variables = {
    'Alpha_THad_P' : ('172.5/m($t_{h}$)', 10, (0., 5.), '$\\alpha_{P}$=Gen P($t_{h}$)/Reco P($t_{h}$)', 5, (0., 10.), 'Reco m($t\\bar{t}$) [GeV]', 1, (200., 2000.)),
    'Alpha_THad_E' : ('172.5/m($t_{h}$)', 10, (0., 5.), '$\\alpha_{E}$=Gen E($t_{h}$)/Reco E($t_{h}$)', 5, (0., 10.), 'Reco m($t\\bar{t}$) [GeV]', 1, (200., 2000.)),
}


    ## get plotting colors/settings
hstyles = styles.styles
stack_fill_opts = {'alpha': 0.8, 'edgecolor':(0,0,0,.5)}
stack_error_opts = {'edgecolor':(0,0,0,.5)}

    ## get data lumi and scale MC by lumi
data_lumi_dict = prettyjson.loads(open('%s/inputs/lumis_data.json' % proj_dir).read())
lumi_correction = load('%s/Corrections/%s/MC_LumiWeights_IgnoreSigEvts.coffea' % (proj_dir, jobid))

        # scale ttJets events, split by reconstruction type, by normal ttJets lumi correction
ttJets_permcats = ['*right', '*matchable', '*unmatchable', '*other']


## make groups based on process
process = hist.Cat("process", "Process", sorting='placement')
process_cat = "dataset"



def get_mtt_bins(histo):
    mtt_vals = histo.axis('bp_mtt').edges(overflow='all')
    bin_cutoff = np.where(mtt_vals > 1000)[0][0]
    val_array = histo.values(sumw2=False, overflow='all')[()]
    array_to_partition = np.concatenate((val_array[0:bin_cutoff], np.array([np.sum(val_array[bin_cutoff:])])))

    partitions, inds_for_splitting = partition_list(array_to_partition, 4)
    mtt_region_edges = [mtt_vals[1]] + mtt_vals[inds_for_splitting].tolist() + [mtt_vals[-2]]
    mtt_region_yields = list(map(sum, partitions))

    mtt_region_splittings = [(mtt_region_edges[idx], mtt_region_edges[idx+1]) for idx in range(len(mtt_region_edges)-1)]
    #set_trace()

    return mtt_region_splittings

def FindMedianAndMedianError(histo):
    if histo.dense_dim() != 1:
        raise ValueError("Hist must be 1D in order to find median")

    vals = histo.values()[()]
    vals[vals < 0] = 0 # set bins with negative contents = 0
    bin_centers = histo.dense_axes()[0].centers()
    if np.sum(vals) <= 0.:
        median, median_error = -10., -10.
    else:
        #print('max:', max(vals), ', min:', min(vals[vals > 0]), ', integral:', np.sum(vals))
        if np.sum(vals) <= 100.:
            multiplier = 1e5
        else:
            multiplier = 1e4 if max(vals) <= 100. else 100
        pseudo_bincounts = np.repeat(bin_centers, (vals*multiplier).astype(int)) ## "recreate" number of counts populating each bin to caluculate median and error on median
        median = round(np.median(pseudo_bincounts)*100, 2)/100
        median_error = 1.2533*stats.sem(pseudo_bincounts) if stats.sem(pseudo_bincounts) != 0. else 1e-10

        #set_trace()
    return median, median_error


def get_median_from_2d(histo, xaxis_name, xmin=None, xmax=None):

    x_edges = histo.dense_axes()[0].edges()
    first_xbin = 0 if not xmin else np.where(histo.axis(xaxis_name).edges() >= xmin)[0][0]
    last_xbin = len(x_edges)-2 if not xmax else np.where(histo.axis(xaxis_name).edges() <= xmax)[0][-1]

    medians, median_errors = [], []
    #set_trace()
    for idx in range(first_xbin, last_xbin+1):
        bin_min, bin_max = x_edges[idx], x_edges[idx+1]

        #print(bin_min)
        hslice = histo[bin_min:bin_max, :].integrate(xaxis_name)
        median, median_error = FindMedianAndMedianError(hslice)

        medians.append(median)
        median_errors.append(median_error)

    #print('  medians:', medians, '  errors:', median_errors)
    return medians, median_errors    

def find_alpha_correction(medians, errors, xbins, output_xbins, ybins=None, output_ybins=None):

    #set_trace()
    if np.ndim(medians) == 1:
        if (medians < 0).sum() > 0:
            raise ValueError("Not all median input values are valid!")
        np_fit = np.polyfit(xbins, medians, 1, w=np.reciprocal(errors))
        fitvals = np.poly1d(np_fit)(output_xbins)
        lookup = dense_lookup(*(fitvals, output_xbins))

    else:
            # get x bincenter range by finding first and last bin that has valid medians for all mtt bins
        first_xbin, last_xbin = np.where((medians > 0).all(axis=0))[0][0], np.where((medians > 0).all(axis=0))[0][-1]
        fit_ybins = np.array( [(ybins[i+1]+ybins[i])/2 for i in range(len(ybins)-1)] ) # ybins to be used in interpolation
        valid_medians = medians[:, first_xbin:last_xbin+1]

        fit = interpolate.interp2d(xbins, fit_ybins, valid_medians, kind='linear')
        fitvals = fit(output_xbins, output_ybins)
        lookup = dense_lookup(*(fitvals, (output_xbins, output_ybins)))

    return lookup, fitvals

    
for year in years_to_run:
    input_dir = '/'.join([proj_dir, 'results', '%s_%s' % (year, jobid), analyzer])
    f_ext = 'TOT.coffea'
    outdir = '/'.join([proj_dir, 'plots', '%s_%s' % (year, jobid), analyzer])
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    fnames = sorted(['%s/%s' % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)])

    hdict = plt_tools.add_coffea_files(fnames) if len(fnames) > 1 else load(fnames[0])

    lumi_to_use = (data_lumi_dict[year]['Muons']+data_lumi_dict[year]['Electrons'])/2000.

        # scale ttJets events, split by reconstruction type, by normal ttJets lumi correction
    names = [dataset for dataset in list(set([key[0] for key in hdict[list(variables.keys())[0]].values().keys()]))] # get dataset names in hists
    ttJets_cats = [name for name in names if any([fnmatch.fnmatch(name, cat) for cat in ttJets_permcats])] # gets ttJets(_PS)_other, ...
    if len(ttJets_cats) > 0:
        ttJets_lumi_topo = '_'.join(ttJets_cats[0].split('_')[:-1]) # gets ttJets or ttJets_PS
        mu_lumi = lumi_correction[year]['Muons'][ttJets_lumi_topo]
        el_lumi = lumi_correction[year]['Electrons'][ttJets_lumi_topo]
        lumi_correction[year]['Muons'].update({key: mu_lumi for key in ttJets_cats})
        lumi_correction[year]['Electrons'].update({key: el_lumi for key in ttJets_cats})

        ## make groups based on process
    process_groups = plt_tools.make_dataset_groups('Muon', year, samples=names) # works when only MC present
    for hname in hdict.keys():
        if hname == 'cutflow': continue
        hdict[hname] = hdict[hname].group(process_cat, process, process_groups)


    mthad_fit_range = (0.9, 3.5)
    for hname in variables.keys():
        if not hname in hdict.keys():
            raise ValueError("Hist %s not found" % hname)

        #set_trace()
        histo = hdict[hname]
             ## rescale hist by lumi for muons and electrons separately and then combine
        h_mu = histo[:, 'Muon'].integrate('leptype')
        h_mu.scale(lumi_correction[year]['Muons'], axis='process')
        h_el = histo[:, 'Electron'].integrate('leptype')
        h_el.scale(lumi_correction[year]['Electrons'], axis='process')
        h_tot = h_mu+h_el

        alpha_axis_name = [h_tot.dense_axes()[idx].name for idx in range(len(h_tot.dense_axes())) if 'alpha' in h_tot.dense_axes()[idx].name][0]
        mthad_title, mthad_rebinning, mthad_lims, alpha_title, alpha_rebinning, alpha_lims, mtt_title, mtt_rebinning, mtt_lims = variables[hname]
        if mthad_rebinning != 1:
            h_tot = h_tot.rebin('norm_mthad', mthad_rebinning)
        if alpha_rebinning != 1:
            h_tot = h_tot.rebin(alpha_axis_name, alpha_rebinning)
        if mtt_rebinning != 1:
            h_tot = h_tot.rebin('bp_mtt', mtt_rebinning)

            # get mtt splitting from correct perms
        #mtt_bin_ranges = get_mtt_bins(h_tot['ttJets_right'].integrate('process').integrate('norm_mthad').integrate(alpha_axis_name))
        mtt_bin_ranges = [(200., 350.), (350., 400.), (400., 500.), (500., 700.), (700., 1000.), (1000., 2000.)]

            ## make plots for each perm category
        for cat in ['ttJets_right']:
        #for cat in list(set([key[0] for key in h_tot.values().keys()])):
            pltdir = '/'.join([outdir, cat])
            if not os.path.isdir(pltdir):
                os.makedirs(pltdir)
            
            histo = h_tot[cat].integrate('process')
            opts = {'cmap_label' : 'Events'}
            #set_trace()

            mthad_edges = histo.axis('norm_mthad').edges()
            mthad_bins = mthad_edges[np.where(mthad_edges == mthad_fit_range[0])[0][0]:np.where(mthad_edges == mthad_fit_range[1])[0][0]+1]
            output_mthad_bins = np.linspace(mthad_bins[0], mthad_bins[-1], (mthad_bins.size-1)*10+1)
            mtt_bins = np.unique(np.array(list(set(mtt_bin_ranges))))
            output_mtt_bins = np.linspace(mtt_bins[0], mtt_bins[-1], int((mtt_bins[-1]-mtt_bins[0])/10+1) )

                # plots for mtt ranges
            if cat == 'ttJets_right':
                binned_mtt_medians = np.zeros( (len(mtt_bin_ranges), mthad_bins.size) )
                binned_mtt_errors  = np.zeros( (len(mtt_bin_ranges), mthad_bins.size) )
            for idx in range(len(mtt_bin_ranges)):
                fig, ax = plt.subplots()
                fig.subplots_adjust(hspace=.07)
                bin_min, bin_max = mtt_bin_ranges[idx]

                hslice = histo[bin_min:bin_max, :, :].integrate('bp_mtt')
                if cat == 'ttJets_right':
                    alpha_median, alpha_median_errs = get_median_from_2d(hslice, 'norm_mthad', xmin=mthad_fit_range[0], xmax=mthad_fit_range[1])
                    alpha_medians, alpha_errors = np.array(alpha_median), np.array(alpha_median_errs)
                    binned_mtt_medians[idx] = alpha_medians
                    binned_mtt_errors[idx] = alpha_errors

                ax = Plotter.plot_2d_norm(hslice, xaxis_name='norm_mthad', yaxis_name=alpha_axis_name,
                    values=np.ma.masked_where(hslice.values()[()] <= 0., hslice.values()[()]),
                    xlimits=mthad_lims, ylimits=alpha_lims, xlabel=mthad_title, ylabel=alpha_title, ax=ax, **opts)

                   # add lepton/jet multiplicity label
                ax.text(
                    0.02, 0.90, blurb,
                    fontsize=rcParams['font.size']*0.75, 
                    horizontalalignment='left', 
                    verticalalignment='bottom', 
                    transform=ax.transAxes
                )
                    # add perm category and mtt region
                mtt_label = 'm($t\\bar{t}$) $\geq$ %s' % bin_min if idx == len(mtt_bin_ranges)-1 else '%s $\leq$ m($t\\bar{t}$) $<$ %s' % (bin_min, bin_max)
                ax.text(
                    0.98, 0.90, '%s\n%s' % (hstyles[cat]['name'].split(' ')[-1].capitalize(), mtt_label),
                    fontsize=rcParams['font.size']*0.75,
                    horizontalalignment='right',
                    verticalalignment='bottom',
                    transform=ax.transAxes
                )
                    ## add lumi/cms label
                ax = hep.cms.cmslabel(ax=ax, data=False, paper=False, year=year, lumi=round(lumi_to_use, 1), fontsize=18)

                #set_trace()
                figname = '%s/%s' % (pltdir, '_'.join([hname, 'Mtt%sto%s' % (int(bin_min), int(bin_max))]))
                fig.savefig(figname)
                print('%s written' % figname)
                plt.close()

                ## make 2d interpolated fit
            if cat == 'ttJets_right':
                lookup_mtt, fitvals_mtt = find_alpha_correction(medians=binned_mtt_medians, errors=binned_mtt_errors, xbins=mthad_bins, output_xbins=output_mthad_bins, ybins=mtt_bins, output_ybins=output_mtt_bins)
                alpha_corrections[year][hname.split('_')[-1]].update({'Mtt' : lookup_mtt})
                    # plot alpha correction
                fig, ax = plt.subplots()
                fig.subplots_adjust(hspace=.07)

                ax = Plotter.plot_2d_norm(histo, xlimits=(min(lookup_mtt._axes[0]), max(lookup_mtt._axes[0])), ylimits=mtt_lims, xlabel=mthad_title, ylabel='m($t\\bar{t}$) [GeV]', ax=ax,
                    values=fitvals_mtt.T, xbins=lookup_mtt._axes[0], ybins=lookup_mtt._axes[1], **{'cmap_label' : '%s Fit Values' % alpha_title.split('=')[0]})
                ax = hep.cms.cmslabel(ax=ax, data=False, paper=False, year=year, lumi=round(lumi_to_use, 1), fontsize=18)
                figname = '%s/%s_Mtt_FitVals' % (pltdir, hname)
                fig.savefig(figname)
                print('%s written' % figname)
                plt.close()


                # plots over entire mttbar range
            all_mtt = histo.integrate('bp_mtt')
            if cat == 'ttJets_right':
                alpha_median_all, alpha_median_errs_all = get_median_from_2d(all_mtt, 'norm_mthad', xmin=mthad_fit_range[0], xmax=mthad_fit_range[1])
                alpha_medians_all, alpha_errors_all = np.array(alpha_median_all), np.array(alpha_median_errs_all)
                lookup_all, fitvals_all = find_alpha_correction(medians=alpha_medians_all, errors=alpha_errors_all, xbins=mthad_bins, output_xbins=output_mthad_bins)
                alpha_corrections[year][hname.split('_')[-1]].update({'All' : lookup_all})
                    # plot alpha correction
                fig, ax = plt.subplots()
                fig.subplots_adjust(hspace=.07)

                #set_trace()
                plt.plot(lookup_all._axes, fitvals_all, color='black')
                ax.set_xlabel(mthad_title)
                ax.set_ylabel('%s Fit Values' % alpha_title.split('=')[0])
                ax = hep.cms.cmslabel(ax=ax, data=False, paper=False, year=year, lumi=round(lumi_to_use, 1), fontsize=18)
                figname = '%s/%s_All_FitVals' % (pltdir, hname)
                fig.savefig(figname)
                print('%s written' % figname)
                plt.close()


            fig, ax = plt.subplots()
            fig.subplots_adjust(hspace=.07)

            mtt_range = histo.axis('bp_mtt').edges()[0]
            ax = Plotter.plot_2d_norm(all_mtt, xaxis_name='norm_mthad', yaxis_name=alpha_axis_name,
                values=np.ma.masked_where(all_mtt.values()[()] <= 0., all_mtt.values()[()]),
                xlimits=mthad_lims, ylimits=alpha_lims, xlabel=mthad_title, ylabel=alpha_title, ax=ax, **opts)

               # add lepton/jet multiplicity label
            ax.text(
                0.02, 0.90, blurb,
                fontsize=rcParams['font.size']*0.75, 
                horizontalalignment='left', 
                verticalalignment='bottom', 
                transform=ax.transAxes
            )
                # add perm category and mtt region
            ax.text(
                0.98, 0.90, '%s\nm($t\\bar{t}$) $\geq$ %s' % (hstyles[cat]['name'].split(' ')[-1].capitalize(), mtt_range),
                fontsize=rcParams['font.size']*0.75,
                horizontalalignment='right',
                verticalalignment='bottom',
                transform=ax.transAxes
            )
                ## add lumi/cms label
            ax = hep.cms.cmslabel(ax=ax, data=False, paper=False, year=year, lumi=round(lumi_to_use, 1), fontsize=18)

            #set_trace()
            figname = '%s/%s' % (pltdir, '_'.join([hname, 'All']))
            fig.savefig(figname)
            print('%s written' % figname)
            plt.close()


    ## save corrections
if len(years_to_run) == 3:
    corrdir = '%s/Corrections/%s' % (proj_dir, jobid)
    if not os.path.isdir(corrdir):
        os.makedirs(corrdir)

    alpha_corr_name = '%s/alpha_correction_%s.coffea' % (corrdir, jobid)
    save(alpha_corrections, alpha_corr_name)
    print('\n', alpha_corr_name, 'written')
