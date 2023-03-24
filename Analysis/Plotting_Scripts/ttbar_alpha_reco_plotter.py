# matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend("agg")
from matplotlib import rcParams
rcParams["font.size"] = 20
rcParams["savefig.format"] = "pdf"
rcParams["savefig.bbox"] = "tight"

from coffea.util import load, save
from pdb import set_trace
import os
import Utilities.styles as styles
import Utilities.plot_tools as plt_tools
import Utilities.prettyjson as prettyjson
import re
from coffea import hist
from coffea.hist import plot
import numpy as np
import fnmatch
import Utilities.Plotter as Plotter
from equal_split import partition_list
from scipy import stats
from scipy import interpolate
from coffea.lookup_tools.dense_lookup import dense_lookup
from scipy.optimize import curve_fit
import Utilities.common_features as cfeatures

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
base_jobid = os.environ["base_jobid"]
analyzer = "ttbar_alpha_reco"
plot_outdir = os.environ["plots_dir"]
eos_dir = os.environ["eos_dir"]

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--year", type=str, help="Choose year(s) to run")
parser.add_argument("--force_save", action="store_true", help="Force the efficiencies to be saved.")
args = parser.parse_args()

years_to_run = [args.year] if args.year else ["2016APV", "2016", "2017", "2018"]
max_years = 4
#years_to_run = [args.year] if args.year else ["2017", "2018"]
#max_years = 2

blurb = cfeatures.channel_labels["Lepton_3Jets"]

alpha_corrections = {year : {"E" : {}, "P" : {}} for year in years_to_run}


variables = {
    #"Alpha_THad_P" : ("172.5/Reco m($t_{h}$)", 1, (0., 5.), "Gen P($t_{h}$)/Reco P($t_{h}$)", 1, (0., 5.), "Reco m($t\\bar{t}$) [GeV]", 1, (200., 2000.)),
    #"Alpha_THad_E" : ("172.5/Reco m($t_{h}$)", 1, (0., 5.), "Gen E($t_{h}$)/Reco E($t_{h}$)", 1, (0., 5.), "Reco m($t\\bar{t}$) [GeV]", 1, (200., 2000.)),
    "Alpha_THad_P" : ("172.5/$m_{t_{h}^{Reco}}$", 1, (0., 5.), "$P_{t_{h}^{Gen}}$/$P_{t_{h}^{Reco}}$", 1, (0., 5.), "Reco $m_{t\\bar{t}}$ [GeV]", 1, (200., 2000.)),
    "Alpha_THad_E" : ("172.5/$m_{t_{h}^{Reco}}$", 1, (0., 5.), "$E_{t_{h}^{Gen}}$/$E_{t_{h}^{Reco}}$", 1, (0., 5.), "Reco $m_{t\\bar{t}}$ [GeV]", 1, (200., 2000.)),
}


    ## get plotting colors/settings
hstyles = styles.styles

    ## get data lumi and scale MC by lumi
data_lumi_dict = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())
lumi_correction = load(os.path.join(proj_dir, "Corrections", base_jobid, "MC_LumiWeights.coffea"))

        # scale ttJets events, split by reconstruction type, by normal ttJets lumi correction
ttJets_permcats = ["*right", "*matchable", "*unmatchable", "*sl_tau", "*other"]

## make groups based on process
process = hist.Cat("process", "Process", sorting="placement")
process_cat = "dataset"


# define sigmoid function for fitting
def sigmoid(x, L ,x0, k, b):
    y = L / (1 + np.exp(-k*(x-x0)))+b
    return y

def sqrt_x(x, a, b, c):
    y = a*np.sqrt(x-b)+c
    return y

def get_mtt_bins(histo):
    mtt_vals = histo.axis("bp_mtt").edges(overflow="all")
    bin_cutoff = np.where(mtt_vals > 1000)[0][0]
    val_array = histo.values(sumw2=False, overflow="all")[()]
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

    #set_trace()
    vals = histo.values()[()]
    vals[vals < 0] = 0 # set bins with negative contents = 0
    if np.sum(vals) <= 100.:
        median = -10.
    else:
        norm_vals = vals/np.sum(histo.values(overflow="all")[()])
        cdf = np.cumsum(norm_vals)
        bin_centers = histo.dense_axes()[0].centers()
        median = bin_centers[np.where(cdf > 0.5)[0][0]]

        #set_trace()
    return median


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
        median = FindMedianAndMedianError(hslice)
        #median, median_error = FindMedianAndMedianError(hslice)

        medians.append(median)
        #median_errors.append(median_error)

    #print("  medians:", medians, "  errors:", median_errors)
    return medians    
    #return medians, median_errors    

def find_alpha_correction(medians, xbins, output_xbins, errors=None, ybins=None, output_ybins=None, degree=None, Fit=None):

    #set_trace()
    if np.ndim(medians) == 1:
        if (medians < 0).sum() > 0:
            print("Not all median input values are valid!")
        #    raise ValueError("Not all median input values are valid!")
        #set_trace()
        valid_medians = medians[medians != -10.]
        valid_xbins = xbins[medians != -10.]
        if Fit is None:
            deg = 1 if degree is None else degree
            #set_trace()
            np_fit = np.polyfit(valid_xbins, valid_medians, deg, w=np.reciprocal(errors[medians != -10.]), full=True) if errors is not None else np.polyfit(valid_xbins, valid_medians, deg, full=True)
            chisq_ndof = np_fit[1]/(len(valid_xbins)-(deg+1))
            fitvals = np.poly1d(np_fit[0])(output_xbins)
            lookup = dense_lookup(*(fitvals, output_xbins))

            return lookup, fitvals, chisq_ndof

        elif Fit == "Sqrt":
            #set_trace()
            p0 = [0.5, 0.5, 0.5]
            popt, pcov = curve_fit(sqrt_x, valid_xbins, valid_medians, p0, method="dogbox")
            #popt, pcov = curve_fit(sqrt_x, xbins, medians, p0, method="dogbox")
            #popt, pcov = curve_fit(sqrt_x, xbins, medians, p0, bounds=(xbins[0], xbins[-1]), method="dogbox")
            fitvals = sqrt_x(output_xbins, *popt)
            lookup = dense_lookup(*(fitvals, output_xbins))

            return lookup, fitvals

        elif Fit == "Sigmoid":
            #set_trace()
            set_trace()
            p0 = [1., 1., 1., 0.5] # this is an mandatory initial guess
            #p0 = [max(medians), min(xbins), 1, min(medians)] # this is an mandatory initial guess
            #p0 = [max(medians), np.median(xbins), 1, min(medians)] # this is an mandatory initial guess
            
            popt, pcov = curve_fit(sigmoid, xbins, medians, method="dogbox")
            #popt, pcov = curve_fit(sigmoid, xbins, medians, p0, bounds=(xbins[0], xbins[-1]), method="dogbox")
            fitvals = sigmoid(output_xbins, *popt)
            lookup = dense_lookup(*(fitvals, output_xbins))

            return lookup, fitvals

    else:
            # get x bincenter range by finding first and last bin that has valid medians for all mtt bins
        first_xbin, last_xbin = np.where((medians > 0).all(axis=0))[0][0], np.where((medians > 0).all(axis=0))[0][-1]
        fit_ybins = np.array( [(ybins[i+1]+ybins[i])/2 for i in range(len(ybins)-1)] ) # ybins to be used in interpolation
        valid_medians = medians[:, first_xbin:last_xbin+1]

        deg = "linear" if degree is None else degree
        fit = interpolate.interp2d(xbins, fit_ybins, valid_medians, kind=deg)
        fitvals = fit(output_xbins, output_ybins)
        lookup = dense_lookup(*(fitvals, (output_xbins, output_ybins)))

        return lookup, fitvals

    
for year in years_to_run:
    input_dir = os.path.join(eos_dir, "results", f"{year}_{jobid}", analyzer)
    f_ext = "TOT.coffea"
    fnames = sorted(["%s/%s" % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)])
    hdict = plt_tools.add_coffea_files(fnames) if len(fnames) > 1 else load(fnames[0])

    outdir = os.path.join(plot_outdir, f"{year}_{jobid}", analyzer)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    lumi_to_use = (data_lumi_dict[year]["Muons"]+data_lumi_dict[year]["Electrons"])/2000.

        # scale ttJets events, split by reconstruction type, by normal ttJets lumi correction
    names = [dataset for dataset in sorted(set([key[0] for key in hdict[sorted(variables.keys())[0]].values().keys()]))] # get dataset names in hists
    ttJets_cats = [name for name in names if any([fnmatch.fnmatch(name, cat) for cat in ttJets_permcats])] # gets ttJets(_PS)_other, ...
    if len(ttJets_cats) > 0:
        for tt_cat in ttJets_cats:
            ttJets_lumi_topo = "_".join(tt_cat.split("_")[:-2]) if "sl_tau" in tt_cat else "_".join(tt_cat.split("_")[:-1]) # gets ttJets[SL, Had, DiLep] or ttJets_PS
            mu_lumi = lumi_correction[year]["Muons"][ttJets_lumi_topo]
            el_lumi = lumi_correction[year]["Electrons"][ttJets_lumi_topo]
            lumi_correction[year]["Muons"].update({tt_cat: mu_lumi})
            lumi_correction[year]["Electrons"].update({tt_cat: el_lumi})

        ## make groups based on process
    process_groups = plt_tools.make_dataset_groups("Muon", year, samples=names) # works when only MC present


    mthad_fit_range = (0.9, 2.2) ## chosen because ~95% of events fall in this range
    #mthad_fit_range = (0.9, 3.5)
    for hname in variables.keys():
        if not hname in hdict.keys():
            raise ValueError(f"Hist {hname} not found")

        #set_trace()
        histo = hdict[hname]
             ## rescale hist by lumi for muons and electrons separately and then combine
        h_mu = histo[:, "Muon"].integrate("leptype")
        h_mu.scale(lumi_correction[year]["Muons"], axis="dataset")
        h_el = histo[:, "Electron"].integrate("leptype")
        h_el.scale(lumi_correction[year]["Electrons"], axis="dataset")
        h_tot = h_mu+h_el
        h_tot = h_tot.group(process_cat, process, process_groups)

        alpha_axis_name = [h_tot.dense_axes()[idx].name for idx in range(len(h_tot.dense_axes())) if "alpha" in h_tot.dense_axes()[idx].name][0]
        mthad_title, mthad_rebinning, mthad_lims, alpha_title, alpha_rebinning, alpha_lims, mtt_title, mtt_rebinning, mtt_lims = variables[hname]
        if mthad_rebinning != 1:
            h_tot = h_tot.rebin("norm_mthad", mthad_rebinning)
        if alpha_rebinning != 1:
            h_tot = h_tot.rebin(alpha_axis_name, alpha_rebinning)
        if mtt_rebinning != 1:
            h_tot = h_tot.rebin("bp_mtt", mtt_rebinning)

            # get mtt splitting from correct perms
        #mtt_bin_ranges = get_mtt_bins(h_tot["ttJets_right"].integrate("process").integrate("norm_mthad").integrate(alpha_axis_name))
        mtt_bin_ranges = [(200., 350.), (350., 400.), (400., 500.), (500., 700.), (700., 1000.), (1000., 2000.)]

            ## make plots for each perm category
        for cat in ["ttJets_right"]:
        #for cat in list(set([key[0] for key in h_tot.values().keys()])):
            pltdir = os.path.join(outdir, cat)
            if not os.path.isdir(pltdir):
                os.makedirs(pltdir)
            
            histo = h_tot[cat].integrate("process")
            opts = {"cmap_label" : "%s $t\\bar{t}$ Events" % hstyles[cat]["name"].split(" ")[-1].capitalize()}
            #opts = {"cmap_label" : "Events"}
            #set_trace()

            mthad_edges = histo.axis("norm_mthad").edges()
            mthad_bins = mthad_edges[np.where(np.around(mthad_edges, 5) == mthad_fit_range[0])[0][0]:np.where(np.around(mthad_edges, 5) == mthad_fit_range[1])[0][0]+1]
            output_mthad_bins = np.linspace(mthad_bins[0], mthad_bins[-1], (mthad_bins.size-1)*10+1)
            mtt_bins = np.unique(np.array(list(set(mtt_bin_ranges))))
            output_mtt_bins = np.linspace(mtt_bins[0], mtt_bins[-1], int((mtt_bins[-1]-mtt_bins[0])/10+1) )

                # plots for mtt ranges
            if cat == "ttJets_right":
                binned_mtt_medians = np.zeros( (len(mtt_bin_ranges), mthad_bins.size) )
                binned_mtt_errors  = np.zeros( (len(mtt_bin_ranges), mthad_bins.size) )
            for idx in range(len(mtt_bin_ranges)):
                fig, ax = plt.subplots()
                fig.subplots_adjust(hspace=.07)
                bin_min, bin_max = mtt_bin_ranges[idx]

                hslice = histo[bin_min:bin_max, :, :].integrate("bp_mtt")
                if cat == "ttJets_right":
                    #alpha_median, alpha_median_errs = get_median_from_2d(hslice, "norm_mthad", xmin=mthad_fit_range[0], xmax=mthad_fit_range[1])
                    #alpha_medians, alpha_errors = np.array(alpha_median), np.array(alpha_median_errs)
                    #binned_mtt_medians[idx] = alpha_medians
                    #binned_mtt_errors[idx] = alpha_errors
                    alpha_median = get_median_from_2d(hslice, "norm_mthad", xmin=mthad_fit_range[0], xmax=mthad_fit_range[1])
                    alpha_medians = np.array(alpha_median)
                    binned_mtt_medians[idx] = alpha_medians

                Plotter.plot_2d_norm(hdict=hslice, xaxis_name="norm_mthad", yaxis_name=alpha_axis_name,
                    values=np.ma.masked_where(hslice.values()[()] <= 0., hslice.values()[()]),
                    xlimits=mthad_lims, ylimits=alpha_lims, xlabel=mthad_title, ylabel=alpha_title, ax=ax, **opts)

                mtt_label = "$m_{t\\bar{t}}^{Reco}$ $\geq$ %s" % bin_min if idx == len(mtt_bin_ranges)-1 else "%s $\leq$ $m_{t\\bar{t}}^{Reco}$ $<$ %s" % (bin_min, bin_max)
                   # add lepton/jet multiplicity label
                ax.text(
                    0.02, 0.88, "%s\n%s" % (blurb, mtt_label),
                    #0.02, 0.84, "%s\n%s" % (blurb, mtt_label),
                    #0.02, 0.88, blurb,
                    horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                )
                #    # add perm category and mtt region
                #ax.text(
                #    0.98, 0.90, mtt_label,
                #    horizontalalignment="right", verticalalignment="bottom", transform=ax.transAxes
                #)
                    ## add lumi/cms label
                hep.cms.label(ax=ax, data=False, label="Preliminary", year=cfeatures.year_labels[year], lumi=round(lumi_to_use, 1))
                #hep.cms.label(ax=ax, data=False, year=year, lumi=round(lumi_to_use, 1))

                #set_trace()
                figname = os.path.join(pltdir, "_".join([year, jobid, hname, "Mtt%sto%s" % (int(bin_min), int(bin_max))]))
                fig.savefig(figname)
                print(f"{figname} written")
                plt.close(fig)

                ## make 2d interpolated fit
            if cat == "ttJets_right":
                #set_trace()
                valid_binned_mtt_medians = binned_mtt_medians[np.where((binned_mtt_medians == -10.).sum(axis=1) == 0)]
                valid_mtt_bins = mtt_bins[np.where((binned_mtt_medians == -10.).sum(axis=1) == 0)]
                valid_mtt_bins = np.array(valid_mtt_bins.tolist()+[mtt_bins[np.where((binned_mtt_medians == -10.).sum(axis=1) == 0)[0][-1]+1]]) # add next bin for bin edges
                lookup_mtt, fitvals_mtt = find_alpha_correction(medians=valid_binned_mtt_medians, xbins=mthad_bins, output_xbins=output_mthad_bins, ybins=valid_mtt_bins, output_ybins=output_mtt_bins)
                #lookup_mtt, fitvals_mtt = find_alpha_correction(medians=binned_mtt_medians, errors=binned_mtt_errors, xbins=mthad_bins, output_xbins=output_mthad_bins, ybins=mtt_bins, output_ybins=output_mtt_bins)
                alpha_corrections[year][hname.split("_")[-1]].update({"Mtt" : lookup_mtt})
                    # plot alpha correction
                fig, ax = plt.subplots()
                fig.subplots_adjust(hspace=.07)

                Plotter.plot_2d_norm(hdict=histo, xlimits=(min(lookup_mtt._axes[0]), max(lookup_mtt._axes[0])), ylimits=mtt_lims, xlabel=mthad_title, ylabel="$m_{t\\bar{t}}^{Reco}$) [GeV]", ax=ax,
                    values=fitvals_mtt.T, xbins=lookup_mtt._axes[0], ybins=lookup_mtt._axes[1], **{"cmap_label" : "%s Fit Values" % alpha_title.split("=")[0]})
                #hep.cms.label(ax=ax, data=False, year=year, lumi=round(lumi_to_use, 1))
                hep.cms.label(ax=ax, data=False, label="Preliminary", year=cfeatures.year_labels[year], lumi=round(lumi_to_use, 1))
                figname = os.path.join(pltdir, "_".join([year, jobid, hname, "Mtt_FitVals"]))
                fig.savefig(figname)
                print(f"{figname} written")
                plt.close(fig)


                # plots over entire mttbar range
            all_mtt = histo.integrate("bp_mtt")
            if cat == "ttJets_right":
                #alpha_median_all = get_median_from_2d(all_mtt, "norm_mthad")
                alpha_median_all = get_median_from_2d(all_mtt, "norm_mthad", xmin=mthad_fit_range[0], xmax=mthad_fit_range[1])
                #alpha_median_all, alpha_median_errs_all = get_median_from_2d(all_mtt, "norm_mthad", xmin=mthad_fit_range[0], xmax=mthad_fit_range[1])
                #alpha_medians_all, alpha_errors_all = np.array(alpha_median_all), np.array(alpha_median_errs_all)
                alpha_medians_all = np.array(alpha_median_all)
                valid_medians_all = alpha_medians_all[alpha_medians_all != -10.]
                valid_mthad_bins = mthad_bins[alpha_medians_all != -10.]
                    # linear fit
                lookup_all_1d, fitvals_all_1d, chisq_1d = find_alpha_correction(medians=valid_medians_all, xbins=valid_mthad_bins, output_xbins=output_mthad_bins)
                alpha_corrections[year][hname.split("_")[-1]].update({"All_1D" : lookup_all_1d})
                    # quadratic fit
                lookup_all_2d, fitvals_all_2d, chisq_2d = find_alpha_correction(medians=valid_medians_all, xbins=valid_mthad_bins, output_xbins=output_mthad_bins, degree=2)
                alpha_corrections[year][hname.split("_")[-1]].update({"All_2D" : lookup_all_2d})

                    # plot alpha correction
                fig, ax = plt.subplots()
                fig.subplots_adjust(hspace=.07)

                #set_trace()
                #plt.plot(lookup_all_1d._axes, fitvals_all_1d, color="r", label="Linear Fit, $\\chi^2$/ndof=%.4f" % chisq_1d)
                #plt.plot(lookup_all_2d._axes, fitvals_all_2d, color="b", label="Quadratic Fit, $\\chi^2$/ndof=%.4f" % chisq_2d)
                plt.plot(lookup_all_2d._axes, fitvals_all_2d, color="b", label="$\\alpha_E$ Quadratic Fit" if hname == "Alpha_THad_E" else "$\\alpha_P$ Quadratic Fit")
                #plt.plot(lookup_all_2d._axes, fitvals_all_2d, color="b", label="$\\alpha_E$ Quadratic Fit, $\\chi^2$/ndof=%.4f" % chisq_2d if hname == "Alpha_THad_E" else "$\\alpha_P$ Quadratic Fit, $\\chi^2$/ndof=%.4f" % chisq_2d)
                plt.errorbar(valid_mthad_bins, valid_medians_all, marker="o", color="black", fmt=".", label="Medians")
                ax.autoscale()
                ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.05)
                ax.legend(loc="upper right")
                ax.set_xlabel(mthad_title)
                ax.set_ylabel(alpha_title.split("=")[0])
                hep.cms.label(ax=ax, data=False, label="Preliminary", year=cfeatures.year_labels[year], lumi=round(lumi_to_use, 1))
                #hep.cms.label(ax=ax, data=False, year=year, lumi=round(lumi_to_use, 1))
                figname = os.path.join(pltdir, "_".join([year, jobid, hname, "All_FitVals"]))
                fig.savefig(figname)
                print(f"{figname} written")
                plt.close(fig)


            fig, ax = plt.subplots()
            fig.subplots_adjust(hspace=.07)

            mtt_range = histo.axis("bp_mtt").edges()[0]
            Plotter.plot_2d_norm(hdict=all_mtt, xaxis_name="norm_mthad", yaxis_name=alpha_axis_name,
                values=np.ma.masked_where(all_mtt.values()[()] <= 0., all_mtt.values()[()]),
                xlimits=mthad_lims, ylimits=alpha_lims, xlabel=mthad_title, ylabel=alpha_title, ax=ax, **opts)

               # add lepton/jet multiplicity label
            ax.text(
                0.02, 0.88, "%s\n$m_{t\\bar{t}}^{Reco}$ $\geq$ %s" % (blurb, mtt_range),
                #0.02, 0.84, "%s\n$m_{t\\bar{t}}^{Reco}$ $\geq$ %s" % (blurb, mtt_range),
                #0.02, 0.88, blurb,
                horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
            )
            #    # add perm category and mtt region
            #ax.text(
            #    0.98, 0.90, "$m_{t\\bar{t}}^{Reco}$ $\geq$ %s" % mtt_range,
            #    horizontalalignment="right", verticalalignment="bottom", transform=ax.transAxes
            #)
                ## add lumi/cms label
            hep.cms.label(ax=ax, data=False, label="Preliminary", year=cfeatures.year_labels[year], lumi=round(lumi_to_use, 1))
            #hep.cms.label(ax=ax, data=False, year=year, lumi=round(lumi_to_use, 1))

            #set_trace()
            figname = os.path.join(pltdir, "_".join([year, jobid, hname, "All"]))
            fig.savefig(figname)
            print(f"{figname} written")
            plt.close(fig)


    ## save corrections
if (len(years_to_run) == max_years) or (args.force_save):
    #set_trace()
    corrdir = os.path.join(proj_dir, "Corrections", jobid)
    if not os.path.isdir(corrdir):
        os.makedirs(corrdir)

    alpha_corr_name = os.path.join(corrdir, f"alpha_correction_{jobid}.coffea")
    save(alpha_corrections, alpha_corr_name)
    print("\n", alpha_corr_name, "written")
