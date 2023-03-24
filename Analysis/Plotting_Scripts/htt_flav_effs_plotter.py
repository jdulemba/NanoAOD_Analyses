import numpy as np
from coffea.lookup_tools.dense_lookup import dense_lookup
from coffea.util import load, save
from pdb import set_trace
import os
from coffea import hist
from coffea.hist import plot
# matplotlib
import matplotlib.pyplot as plt
import matplotlib as mpl
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend("agg")
from matplotlib import rcParams
rcParams["font.size"] = 20
rcParams["savefig.format"] = "pdf"
rcParams["savefig.bbox"] = "tight"
import Utilities.prettyjson as prettyjson
import Utilities.common_features as cfeatures

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--year", type=str, help="What year is the ntuple from.")
parser.add_argument("--construct_btag", action="store_false", help="Makes btag SF constructor (default is True)")
parser.add_argument("--no_plots", action="store_true", help="Don't plot efficienies (default is False)")
parser.add_argument("--force_save", action="store_true", help="Force the efficiencies to be saved.")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
base_jobid = os.environ["base_jobid"]
analyzer = "htt_flav_effs"
plot_outdir = os.environ["plots_dir"]
eos_dir = os.environ["eos_dir"]

outdir = os.path.join(proj_dir, "Corrections", jobid)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

years_to_run = [args.year] if args.year else ["2016APV", "2016", "2017", "2018"]
max_years = 4

flav_effs = {year : {"DeepJet" : {"3Jets" : {}, "4PJets" : {}}, "DeepCSV" : {"3Jets" : {},"4PJets" : {},} } for year in years_to_run}

if args.construct_btag:
    from copy import deepcopy
    btag_constructs_dict = deepcopy(flav_effs)

flav_to_name = {"bjet" : "bottom", "cjet" : "charm", "ljet" : "light"}
hname = "Jets_pt_eta"
lumi_correction = load(os.path.join(proj_dir, "Corrections", base_jobid, "MC_LumiWeights.coffea"))

#pt_binning = np.array([30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 125.0, 150.0,170.0, 200.0, 250.0, 1000.0])
#eta_binning = np.array([-2.5, -1.5, -0.5, 0.0, 0.5, 1.5, 2.5])
pt_binning = np.array([30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 125.0, 150.0,170.0, 200.0, 1000.0])
eta_binning = np.array([-2.5, -2., -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5])
pt_bins = hist.Bin("pt", "pt", pt_binning)
eta_bins = hist.Bin("eta", "eta", eta_binning)

working_points = []

def plot_effs(heff, edges, lumi_to_use, year, jmult, btagger, wp, flav, plotdir, clear=True):
    fig, ax = plt.subplots()
    fig.subplots_adjust(hspace=.07)

    opts = {}
    #opts = {"cmap" : "OrRd"}
    xedges, yedges = edges[0], edges[1]
    sumw = heff._values
    pc = ax.pcolormesh(xedges, yedges, sumw.T, **opts)
    ax.add_collection(pc)
    if clear:
        fig.colorbar(pc, ax=ax, label=f"{flav_to_name[flav]} jet Efficiency, {btagger} {wp}", pad=0.)
    ax.autoscale()
    ax.set_xlim(xedges[0], xedges[-2]+50.)
    #ax.set_xlim(xedges[0], xedges[-1])
    ax.set_ylim(yedges[0], yedges[-1])

    ax.set_xscale("log")
    ax.xaxis.set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.xaxis.set_minor_formatter(mpl.ticker.ScalarFormatter())
    ax.ticklabel_format(style="plain", axis="x",useOffset=False)
    
        ## set axes labels and titles
    ax.set_xlabel("$p_{T}$ [GeV]")
    ax.set_ylabel("$\\eta$")
        # add lepton/jet multiplicity label
    ax.text(
        0.02, 0.95, cfeatures.channel_labels[f"Lepton_{jmult}"],
        fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
    )
    hep.cms.label(ax=ax, fontsize=rcParams["font.size"], data=False, year=cfeatures.year_labels[year], lumi=round(lumi_to_use, 1))

    figname = os.path.join(plotdir, f"{btagger}_{wp}_{jmult}_{flav}_Efficiency")
    fig.savefig(figname)
    print(f"{figname} written")
    plt.close(fig)


data_lumi_dict = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())
combine_2016 = ("2016" in years_to_run) and ("2016APV" in years_to_run) and (base_jobid == "ULnanoAOD")
computed_combined_2016 = False

    # ZJets Summer20UL samples have too many negative contributions
import re
non_ZJets_samples = re.compile("(?!ZJets*)")

for year in years_to_run:
    print(year)

    f_ext = "TOT.coffea"
    if combine_2016 and ("2016" in year):
        if computed_combined_2016:
            computed_combined_2016_year_to_copy = year
            continue
        input_dir_2016 = os.path.join(proj_dir, "results", f"2016_{jobid}", analyzer)
        input_dir_2016APV = os.path.join(proj_dir, "results", f"2016APV_{jobid}", analyzer)
        pltdir = os.path.join(proj_dir, "plots", f"2016Combined_{jobid}", analyzer)

        fnames_2016 = sorted([os.path.join(input_dir_2016, fname) for fname in os.listdir(input_dir_2016) if fname.endswith(f_ext)])
        fnames_2016APV = sorted([os.path.join(input_dir_2016APV, fname) for fname in os.listdir(input_dir_2016APV) if fname.endswith(f_ext)])

        hdict_2016 = plt_tools.add_coffea_files(fnames_2016) if len(fnames_2016) > 1 else load(fnames_2016[0])
        hdict_2016APV = plt_tools.add_coffea_files(fnames_2016APV) if len(fnames_2016APV) > 1 else load(fnames_2016APV[0])

        lumi_to_use_2016 = (data_lumi_dict["2016"]["Muons"]+data_lumi_dict["2016"]["Electrons"])/2000.
        lumi_to_use_2016APV = (data_lumi_dict["2016APV"]["Muons"]+data_lumi_dict["2016APV"]["Electrons"])/2000.
        lumi_to_use = lumi_to_use_2016+lumi_to_use_2016APV

        computed_combined_2016 = True
        computed_combined_2016_year_key = year

        hpass_16 = hdict_2016[f"{hname}_pass"]
        hall_16 = hdict_2016[f"{hname}_all"]
        hpass_16APV = hdict_2016APV[f"{hname}_pass"]
        hall_16APV = hdict_2016APV[f"{hname}_all"]

        #set_trace()
            ## rescale hist by lumi for muons and electrons separately and then combine
        hpass_mu_16 = hpass_16[:, :, :, "Muon", :].integrate("leptype")
        hpass_mu_16.scale(lumi_correction["2016"]["Muons"], axis="dataset")
        hpass_mu_16APV = hpass_16APV[:, :, :, "Muon", :].integrate("leptype")
        hpass_mu_16APV.scale(lumi_correction["2016APV"]["Muons"], axis="dataset")
        hpass_el_16 = hpass_16[:, :, :, "Electron", :].integrate("leptype")
        hpass_el_16.scale(lumi_correction["2016"]["Electrons"], axis="dataset")
        hpass_el_16APV = hpass_16APV[:, :, :, "Electron", :].integrate("leptype")
        hpass_el_16APV.scale(lumi_correction["2016APV"]["Electrons"], axis="dataset")
        hpass_tot = hpass_mu_16+hpass_el_16 + hpass_mu_16APV+hpass_el_16APV

        hall_mu_16 = hall_16[:, :, :, "Muon", :].integrate("leptype")
        hall_mu_16.scale(lumi_correction["2016"]["Muons"], axis="dataset")
        hall_mu_16APV = hall_16APV[:, :, :, "Muon", :].integrate("leptype")
        hall_mu_16APV.scale(lumi_correction["2016APV"]["Muons"], axis="dataset")
        hall_el_16 = hall_16[:, :, :, "Electron", :].integrate("leptype")
        hall_el_16.scale(lumi_correction["2016"]["Electrons"], axis="dataset")
        hall_el_16APV = hall_16APV[:, :, :, "Electron", :].integrate("leptype")
        hall_el_16APV.scale(lumi_correction["2016APV"]["Electrons"], axis="dataset")
        hall_tot = hall_mu_16+hall_el_16 + hall_mu_16APV+hall_el_16APV

        if not hpass_tot.compatible(hall_tot):
            raise ValueError("Passing and All hists don't have the same binning!")

    else:
        input_dir = os.path.join(eos_dir, "results", f"{year}_{jobid}", analyzer)
        fnames = [f"{input_dir}/{fname}" for fname in os.listdir(input_dir) if fname.endswith("TOT.coffea")]
        if len(fnames) > 1:
            raise ValueError("Only one TOT file should be used")
        hdict = load(fnames[0])

            ## get data lumi and scale MC by lumi
        lumi_to_use = (data_lumi_dict[year]["Muons"]+data_lumi_dict[year]["Electrons"])/2000.
    
        pltdir = os.path.join(plot_outdir, f"{year}_{jobid}", analyzer)

        hpass = hdict[f"{hname}_pass"]
        hall = hdict[f"{hname}_all"]

            ## rescale hist by lumi for muons and electrons separately and then combine
        hpass_mu = hpass[:, :, :, "Muon", :].integrate("leptype")
        hpass_mu.scale(lumi_correction[year]["Muons"], axis="dataset")
        hpass_el = hpass[:, :, :, "Electron", :].integrate("leptype")
        hpass_el.scale(lumi_correction[year]["Electrons"], axis="dataset")
        hpass_tot = hpass_mu+hpass_el

        hall_mu = hall[:, :, :, "Muon", :].integrate("leptype")
        hall_mu.scale(lumi_correction[year]["Muons"], axis="dataset")
        hall_el = hall[:, :, :, "Electron", :].integrate("leptype")
        hall_el.scale(lumi_correction[year]["Electrons"], axis="dataset")
        hall_tot = hall_mu+hall_el

        # rebin
    hpass_tot = hpass_tot.rebin("pt", pt_bins).rebin("eta", eta_bins)
    hall_tot = hall_tot.rebin("pt", pt_bins).rebin("eta", eta_bins)

    #    # remove ZJets
    #hpass_tot = hpass_tot[:, non_ZJets_samples]
    #hall_tot = hall_tot[:, non_ZJets_samples]

    if not hpass_tot.compatible(hall_tot):
        raise ValueError("Passing and All hists don't have the same binning!")

    if not os.path.isdir(pltdir):
        os.makedirs(pltdir)

    for wp in hpass_tot.axis("btagger")._sorted: # [DEEPCSVMEDIUM, DEEPJETMEDIUM]
        for jmult in hpass_tot.axis("jmult")._sorted: #["3Jets", "4PJets"]
            for flav in hpass_tot.axis("hFlav")._sorted: # [bjet, cjet, ljet]
                print(wp, jmult, flav)
                #set_trace()
                    # get passing and all hists for 3 and 4+ jets separately, only as a function of pT and eta
                h_pass = hpass_tot[wp, :, jmult, flav].integrate("btagger").integrate("dataset").integrate("jmult").integrate("hFlav")
                h_all = hall_tot[wp, :, jmult, flav].integrate("btagger").integrate("dataset").integrate("jmult").integrate("hFlav")

                edges = (h_pass.axis("pt").edges(), h_pass.axis("eta").edges())
                pass_lookup = dense_lookup(*(h_pass.values().values()), edges)
                    ## check that number of entries in passing bins > 20.
                if min([min(val) for val in pass_lookup._values]) < 20.:
                    print(f"bin for hist {wp}, {flav}, {jmult} has bin entries < 20: %s" % [min(val) for val in pass_lookup._values])
                    #raise ValueError("bin for hist %s, %s, %s has bin entries < 20: %s" % (wp, flav, jmult, [min(val) for val in pass_lookup._values]))

                all_lookup = dense_lookup(*(h_all.values().values()), edges)
                eff_lookup = dense_lookup(pass_lookup._values/all_lookup._values, edges)

                tagger = "DeepCSV" if wp.upper().startswith("DEEPCSV") else "DeepJet"
                working_points.append(wp.upper().split(tagger.upper())[-1])
                flav_effs[year][tagger][jmult].update({flav_to_name[flav] : eff_lookup})

                if not args.no_plots: plot_effs(eff_lookup, edges, lumi_to_use, year, jmult, tagger, wp.upper().split(tagger.upper())[-1][0], flav, pltdir)

wp_name = list(set(working_points))[0]

    # save files
if (len(years_to_run) == max_years) or (args.force_save):
    flav_effs_name = os.path.join(outdir, f"htt_3PJets_{wp_name}_flavour_efficiencies_{jobid}.coffea")
    
        # 2016 and 2016APV are the same if they"re computed together
    if combine_2016:
        flav_effs[computed_combined_2016_year_to_copy] = flav_effs[computed_combined_2016_year_key]

    save(flav_effs, flav_effs_name)
    print("\n", flav_effs_name, "written")


if args.construct_btag:
    import python.BTagScaleFactors as btagSF

    cfg_file = prettyjson.loads(open(os.path.join(proj_dir, "cfg_files", f"cfg_pars_{jobid}.json")).read())

    #set_trace()
    for year in flav_effs.keys():
        for btagger in flav_effs[year].keys():
            csv_path = os.path.join(proj_dir, "inputs", "data", base_jobid, "btagSFs", btagSF.btag_csvFiles[year][btagger])
            if not os.path.isfile(csv_path):
                raise IOError(f"BTagging csv file {csv_path} not found.")
            print(f"csv file: {csv_path}")

            for njets_cat in flav_effs[year][btagger].keys():
                eff_dict = flav_effs[year][btagger][njets_cat]
                sf_computer = btagSF.BTagSF(
                    csv = csv_path,
                    wp_key = (btagger, "used", wp_name.lower().capitalize()),
                    effs = eff_dict
                )

                print(f"BTag SF constructed for {year}, {btagger} {wp_name.lower().capitalize()} wp, {njets_cat}")
                btag_constructs_dict[year][btagger][njets_cat].update({wp_name.lower().capitalize() : sf_computer})

    #set_trace()
        # save files
    if (len(years_to_run) == max_years) or (args.force_save):
        btagSFs_name = os.path.join(outdir, f"htt_3PJets_{wp_name}_btag_scalefactors_Subsources_{jobid}.coffea")
        #btagSFs_name = os.path.join(outdir, f"htt_3PJets_{wp_name}_btag_scalefactors_{jobid}.coffea")
        save(btag_constructs_dict, btagSFs_name)
        print("\n", btagSFs_name, "written")
