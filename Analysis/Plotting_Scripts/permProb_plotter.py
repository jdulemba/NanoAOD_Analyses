from coffea.hist import plot
# matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend("agg")
from matplotlib import rcParams
rcParams["font.size"] = 20
rcParams["savefig.format"] = "png"
rcParams["savefig.bbox"] = "tight"

from coffea.util import load, save
from pdb import set_trace
import os
from Utilities import styles
import Utilities.plot_tools as plt_tools
import Utilities.prettyjson as prettyjson
from coffea import hist
from coffea.lookup_tools.dense_lookup import dense_lookup
import numpy as np
from Utilities import Plotter as Plotter
from scipy import interpolate

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("--njets", type=str, help="Choose jet multiplicity.")
parser.add_argument("--combine_2016", action="store_true", help="Combine 2016APV and 2016 distributions into one")
parser.add_argument("--force_save", action="store_true", help="Force the dict of probabilities to be saved.")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
base_jobid = os.environ["base_jobid"]
analyzer = "permProbComputer"

years_to_run = [args.year] if args.year else ["2016APV", "2016", "2017", "2018"]
max_years = 4

njets_to_run = ["3", "4+"]
if args.njets == "3":
    njets_to_run = ["3"]
elif args.njets == "4+":
    njets_to_run = ["4+"]
else:
    njets_to_run = ["3", "4+"]

combine_2016 = ("2016" in years_to_run) and ("2016APV" in years_to_run) and (base_jobid != "NanoAODv6") and (args.combine_2016)
permProbs = {year : {"3Jets" : {},"4PJets" : {}}  for year in years_to_run}

jet_mults = {
    "3Jets" : "3 jets",
    "4PJets" : "4+ jets",
    #"4Jets" : "4 jets",
    #"5Jets" : "5 jets",
    #"6PJets" : "5+ jets"
}

lep_cats = {
    "LoT" : "loose or tight $e/\mu$",
    "Tight" : "$e/\mu$",
    #"Tight" : "tight $e/\mu$",
}

perm_cats = {
    "3Jets" : {
        "Correct" : ["Correct"],
        "Wrong" : ["Wrong"],
    },
    "4PJets" : {
        "1D" : {
            "Correct_BLep" : ["Correct", "Right_TLep"],
            "Wrong_BLep" : ["Right_THad", "Wrong"]
        },
        "2D" : {
            "Correct_THad" : ["Correct", "Right_THad"],
            "Wrong_THad" : ["Right_TLep", "Wrong"],
        },
    },
}

variables = {}
if "3" in njets_to_run:
    variables.update({
        "Lost_nusolver_chi2" : ("$\\chi_{\\nu}^{2}$", 5, (0., 1000.), True, True),
        "Lost_nusolver_dist" : ("$D_{\\nu, min}$ [GeV]", 1, (0., 150.), True, True),
        "Lost_mTHadProxy" : ("$m_{t_{h}^{proxy}}$ [GeV]", 1, (0., 500.), True, True),
        "Merged_nusolver_chi2" : ("$\\chi_{\\nu}^{2}$", 5, (0., 1000.), True, True),
        "Merged_nusolver_dist" : ("$D_{\\nu, min}$ [GeV]", 1, (0., 150.), True, True),
        "Merged_mTHadProxy_vs_maxmjet" : ("max m(jet) [GeV]", "$m_{t_{h}^{proxy}}$ [GeV]", 1, (0., 150.), 1, (0., 500.), True, True),
    })
if "4+" in njets_to_run:
    variables.update({
        "nusolver_chi2" : ("$\\chi_{\\nu}^{2}$", 5, (0., 1000.), True, False),
        "nusolver_dist" : ("$D_{\\nu, min}$ [GeV]", 1, (0., 150.), True, False),
        "mWHad_vs_mTHad" : ("$m_{t_{h}}$ [GeV]", "$m_{W_{h}}$ [GeV]", 10, (0., 500.), 10, (0., 500.), True, False),
    })


    ## get plotting colors/settings
hstyles = styles.styles

    ## get data lumi and scale MC by lumi
data_lumi_dict = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())
lumi_correction = load(os.path.join(proj_dir, "Corrections", jobid, "MC_LumiWeights.coffea"))

    ## make groups based on perm category
pcat = hist.Cat("permcat", "Perm Category", sorting="placement")
pcat_cat = "permcat"


computed_combined_2016 = False
for year in years_to_run:
    f_ext = "TOT.coffea"
    if combine_2016 and ("2016" in year):
        if computed_combined_2016:
            computed_combined_2016_year_to_copy = year
            continue
        input_dir_2016 = os.path.join(proj_dir, "results", f"2016_{jobid}", analyzer)
        input_dir_2016APV = os.path.join(proj_dir, "results", "2016APV_{jobid}", analyzer)
        outdir = os.path.join(proj_dir, "plots", "2016Combined_{jobid}", analyzer)

        fnames_2016 = sorted([os.path.join(input_dir_2016, fname) for fname in os.listdir(input_dir_2016) if fname.endswith(f_ext)])
        fnames_2016APV = sorted([os.path.join(input_dir_2016APV, fname) for fname in os.listdir(input_dir_2016APV) if fname.endswith(f_ext)])
        
        hdict_2016 = plt_tools.add_coffea_files(fnames_2016) if len(fnames_2016) > 1 else load(fnames_2016[0])
        hdict_2016APV = plt_tools.add_coffea_files(fnames_2016APV) if len(fnames_2016APV) > 1 else load(fnames_2016APV[0])

        lumi_to_use_2016 = (data_lumi_dict["2016"]["Muons"]+data_lumi_dict["2016"]["Electrons"])/2000.
        lumi_to_use_2016APV = (data_lumi_dict["2016APV"]["Muons"]+data_lumi_dict["2016APV"]["Electrons"])/2000.
        lumi_to_use = lumi_to_use_2016+lumi_to_use_2016APV

        computed_combined_2016 = True
        computed_combined_2016_year_key = year

    else:
        input_dir = os.path.join(proj_dir, "results", f"{year}_{jobid}", analyzer)
        outdir = os.path.join(proj_dir, "plots", f"{year}_{jobid}", analyzer)
        
        fnames = [os.path.join(input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)]
        fnames = sorted(fnames)
        
        hdict = plt_tools.add_coffea_files(fnames) if len(fnames) > 1 else load(fnames[0])

        lumi_to_use = (data_lumi_dict[year]["Muons"]+data_lumi_dict[year]["Electrons"])/2000.

    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    for hname in variables.keys():
        #set_trace()
        if combine_2016 and ("2016" in year):
            #set_trace()
            #if hname not in hdict_2016.keys():
            #    raise ValueError("%s not found in file" % hname)
            histo_16 = hdict_2016[hname]
            histo_16APV = hdict_2016APV[hname]

                ## rescale hist by lumi for muons and electrons separately and then combine
                    # 2016
            h_mu_16 = histo_16[:, :, "Muon", :, :, :].integrate("leptype").integrate("jmult")
            h_mu_16.scale(lumi_correction["2016"]["Muons"], axis="dataset")
            h_el_16 = histo_16[:, :, "Electron", :, :, :].integrate("leptype").integrate("jmult")
            h_el_16.scale(lumi_correction["2016"]["Electrons"], axis="dataset")
                    # 2016APV
            h_mu_16APV = histo_16APV[:, :, "Muon", :, :, :].integrate("leptype").integrate("jmult")
            h_mu_16APV.scale(lumi_correction["2016APV"]["Muons"], axis="dataset")
            h_el_16APV = histo_16APV[:, :, "Electron", :, :, :].integrate("leptype").integrate("jmult")
            h_el_16APV.scale(lumi_correction["2016APV"]["Electrons"], axis="dataset")

            h_tot = h_mu_16+h_el_16 + h_mu_16APV+h_el_16APV
            h_tot = h_tot[:, :, "MTHigh", :].integrate("dataset").integrate("mtregion")

        else:
            if hname not in hdict.keys():
                raise ValueError(f"{hname} not found in file")

            histo = hdict[hname]
                ## rescale hist by lumi for muons and electrons separately and then combine
            h_mu = histo[:, :, "Muon", :, :, :].integrate("leptype").integrate("jmult")
            h_mu.scale(lumi_correction[year]["Muons"], axis="dataset")
            h_el = histo[:, :, "Electron", :, :, :].integrate("leptype").integrate("jmult")
            h_el.scale(lumi_correction[year]["Electrons"], axis="dataset")
            h_tot = h_mu+h_el
            h_tot = h_tot[:, :, "MTHigh", :].integrate("dataset").integrate("mtregion")

        print(year, hname)
        var_pars = variables[hname]
        if var_pars[-1]: # distribution is for 3 jets

            topo = hname.split("_")[0] # Merged/Lost
    
                ## make groups based on perm category
            pcat_groups = perm_cats["3Jets"]

            if h_tot.dense_dim() == 1:
                xtitle, rebinning, x_lims, save_dist, is3jet_dist = var_pars
        
                ## hists should have 3 category axes (jet multiplicity, lepton type, perm category) followed by variable
                for lepcat in sorted(set([key[0] for key in h_tot.values().keys()])):
                    pltdir = os.path.join(outdir, "3Jets", lepcat, topo)
                    if not os.path.isdir(pltdir):
                        os.makedirs(pltdir)
        
                    fig, ax = plt.subplots()
                    fig.subplots_adjust(hspace=.07)
                    hslice = h_tot[lepcat].integrate("lepcat")
                    hslice = hslice.group(pcat_cat, pcat, pcat_groups) ## make groups based on perm category
        
                    if rebinning != 1:
                        xaxis_name = hslice.dense_axes()[0].name
                        hslice = hslice.rebin(xaxis_name, rebinning)

                        # save distribution to dict
                    if save_dist and (lepcat == "Tight"):
                        hcor = hslice["Correct"].integrate("permcat")
                        edges = hcor.dense_axes()[0].edges()
                        lookup = dense_lookup(*(hcor.values().values()), edges) # not normalized
                        permProbs[year]["3Jets"].update({hname : lookup})
    
                        ## plot comparison of perm categories
                    plot.plot1d(hslice,
                        overlay=hslice.axes()[0].name,
                        ax=ax,
                        clear=False,
                        line_opts={"linestyle" : "-"},
                        density=True, # normalized to 1,
                    )
                    ax.autoscale()#, tight=True)
                    ax.set_ylim(0, ax.get_ylim()[1]*1.15)
                    ax.set_ylabel("Probability Density")
                    ax.set_xlabel(xtitle)
                    ax.set_xlim(x_lims)

                    #set_trace()        
                       ## set legend and corresponding colors
                    if hname == "Lost_mTHadProxy":
                        h_opts = {key: hstyles[f"{key}_THad"] for key in pcat_groups.keys()}
                    else:
                        h_opts = {key: hstyles[f"{key}_BLep" if "nusolver" in hname else key] for key in pcat_groups.keys()}
                    handles, labels = ax.get_legend_handles_labels()
                    for idx, cat in enumerate(labels):
                        labels[idx] = h_opts[cat]["name"]
                        handles[idx].set_color(h_opts[cat]["color"])
                        handles[idx].set_linestyle(h_opts[cat]["linestyle"])
                    # call ax.legend() with the new values
                    ax.legend(handles,labels, loc="upper right", title=f"{topo} Jet" if topo == "Lost" else f"{topo} Jets")
        
                        # add lep category
                    ax.text(
                        0.02, 0.91, "%s\n%s" % (lep_cats[lepcat], jet_mults["3Jets"]),
                        fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                    )
                        ## set axes labels and titles
                    hep.cms.label(ax=ax, data=False, year="2016APV+2016" if (combine_2016 and ("2016" in year)) else year, lumi=round(lumi_to_use, 1))
        
                    #set_trace()
                    figname = os.path.join(pltdir, "_".join([year, jobid, "3Jets", lepcat, hname]))
                    fig.savefig(figname)
                    print(f"{figname} written")
                    plt.close()
    
            else:
                xtitle, ytitle, x_rebinning, x_lims, y_rebinning, y_lims, save_dist, is3jet_dist = variables[hname]
                topo, yvar, useless, xvar = hname.split("_") # Merged/Lost
        
                #set_trace()
                ## hists should have 3 category axes (jet multiplicity, lepton category, permutation category) followed by variable(s)
                for lepcat in sorted(set([key[0] for key in h_tot.values().keys()])):
                    pltdir = os.path.join(outdir, "3Jets", lepcat, topo)
                    if not os.path.isdir(pltdir):
                        os.makedirs(pltdir)
        
                    hslice = h_tot[lepcat].integrate("lepcat")
                    hslice = hslice.group(pcat_cat, pcat, pcat_groups) ## make groups based on perm category
        
                    if x_rebinning != 1:
                        xaxis_name = hslice.dense_axes()[0].name
                        hslice = hslice.rebin(xaxis_name, x_rebinning)
                    if y_rebinning != 1:
                        yaxis_name = hslice.dense_axes()[1].name
                        hslice = hslice.rebin(yaxis_name, y_rebinning)
        
                        # save distribution to dict
                    if save_dist and (lepcat == "Tight"):
                        hcor = hslice["Correct"].integrate("permcat")
                        edges = (hcor.dense_axes()[0].edges(), hcor.dense_axes()[1].edges())
                        lookup = dense_lookup(*(hcor.values().values()), edges) # not normalized
                        permProbs[year]["3Jets"].update({hname : lookup})
        
        
                        # make 1D projection along dense axes
                    for dax in range(2):
                        fig, ax = plt.subplots()
                        fig.subplots_adjust(hspace=.07)
        
                        hproj = hslice.integrate(hslice.dense_axes()[dax].name)
                            ## plot comparison of perm categories
                        plot.plot1d(hproj,
                            overlay=hproj.axes()[0].name,
                            ax=ax,
                            clear=False,
                            line_opts={"linestyle" : "-"},
                            density=True, # normalized to 1,
                        )
                        ax.autoscale()#axis="x", tight=True)
                        ax.set_ylim(0, ax.get_ylim()[1]*1.15)
                        ax.set_ylabel("Probability Density")
                        ax.set_xlabel(ytitle) if dax == 0 else ax.set_xlabel(xtitle)
                        ax.set_xlim(y_lims) if dax == 0 else ax.set_xlim(x_lims)
        
                           ## set legend and corresponding colors
                        h_opts ={key: hstyles[key] for key in pcat_groups.keys()}
                        handles, labels = ax.get_legend_handles_labels()
                        for idx, cat in enumerate(labels):
                            labels[idx] = h_opts[cat]["name"]
                            handles[idx].set_color(h_opts[cat]["color"])
                            handles[idx].set_linestyle(h_opts[cat]["linestyle"])
                        # call ax.legend() with the new values
                        ax.legend(handles,labels, loc="upper right", title=f"{topo} Jet" if topo == "Lost" else f"{topo} Jets")
        
                            # add lep category
                        ax.text(
                            0.02, 0.91, "%s\n%s" % (lep_cats[lepcat], jet_mults["3Jets"]),
                            fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                        )
        
                            ## set axes labels and titles
                        hep.cms.label(ax=ax, fontsize=rcParams["font.size"], data=False, year="2016APV+2016" if (combine_2016 and ("2016" in year)) else year, lumi=round(lumi_to_use, 1))
        
                        figname = os.path.join(pltdir, "_".join([year, jobid, "3Jets", lepcat, topo, yvar if dax == 0 else xvar]))
                        fig.savefig(figname)
                        print(f"{figname} written")
                        plt.close()
    
                        # make 2D plots for each permutation category
                    for cat in sorted(set([key[0] for key in hslice.values().keys()])):            
                        hcat = hslice[cat].integrate("permcat")
        
                            ## normalized plots
                        fig, ax = plt.subplots()
                        fig.subplots_adjust(hspace=.07)
        
                        norm_values = hcat.values()[()]/np.sum(hcat.values()[()]) ## normalized array of hist values
                        opts = {"cmap_label" : "$P_{M}$"}
                        Plotter.plot_2d_norm(xaxis_name=hcat.axes()[0].name, yaxis_name=hcat.axes()[1].name,
                            values=np.ma.masked_where(norm_values <= 0.0, norm_values), # mask nonzero probabilities for plotting
                            xlimits=x_lims, ylimits=y_lims, xlabel=xtitle, ylabel=ytitle,
                            hdict=hcat, ax=ax, **opts)
    
                            # add lep category
                        ax.text(
                            0.02, 0.91, "%s\n%s" % (lep_cats[lepcat], jet_mults["3Jets"]),
                            fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                        )
                            # add perm category
                        ax.text(
                            1, 0.95, cat,
                            fontsize=rcParams["font.size"], 
                            horizontalalignment="right", verticalalignment="bottom", transform=ax.transAxes
                        )
                        hep.cms.label(ax=ax, fontsize=rcParams["font.size"], data=False, year="2016APV+2016" if (combine_2016 and ("2016" in year)) else year, lumi=round(lumi_to_use, 1))
        
                        figname = os.path.join(pltdir, "_".join([year, jobid, "3Jets", lepcat, hname, cat, "norm"]))
                        fig.savefig(figname)
                        print(f"{figname} written")
                        plt.close()


        else: # distribution is for 4+ jets
            if h_tot.dense_dim() == 1:
                xtitle, rebinning, x_lims, save_dist, is3jet_dist = var_pars

                ## hists should have 3 category axes (jet multiplicity, lepton type, perm category) followed by variable
                for lepcat in sorted(set([key[0] for key in h_tot.values().keys()])):
                    pltdir = os.path.join(outdir, "4PJets", lepcat)
                    if not os.path.isdir(pltdir):
                        os.makedirs(pltdir)
        
                    fig, ax = plt.subplots()
                    fig.subplots_adjust(hspace=.07)
                    hslice = h_tot[lepcat].integrate("lepcat")
                        ## make groups based on perm category
                    pcat_groups = perm_cats["4PJets"]["1D"]
                    hslice = hslice.group(pcat_cat, pcat, pcat_groups)
        
                    if rebinning != 1:
                        xaxis_name = hslice.dense_axes()[0].name
                        hslice = hslice.rebin(xaxis_name, rebinning)
        
                        # save distribution to dict
                    if save_dist and (lepcat == "Tight"):
                        hcor = hslice["Correct_BLep"].integrate("permcat")
                        edges = hcor.dense_axes()[0].edges()
                        lookup = dense_lookup(*(hcor.values().values()), edges) # not normalized
                        permProbs[year]["4PJets"].update({hname : lookup})
    
                        ## plot comparison of perm categories
                    plot.plot1d(hslice,
                        overlay=hslice.axes()[0].name,
                        ax=ax,
                        clear=False,
                        line_opts={"linestyle" : "-"},
                        density=True, # normalized to 1,
                    )
                    ax.autoscale()#axis="x", tight=True)
                    ax.set_ylim(0, ax.get_ylim()[1]*1.15)
                    ax.set_ylabel("Probability Density")
                    ax.set_xlabel(xtitle)
                    ax.set_xlim(x_lims)
        
                       ## set legend and corresponding colors
                    h_opts ={key: hstyles[key] for key in pcat_groups.keys()}
                    handles, labels = ax.get_legend_handles_labels()
                    for idx, cat in enumerate(labels):
                        labels[idx] = h_opts[cat]["name"]
                        handles[idx].set_color(h_opts[cat]["color"])
                        handles[idx].set_linestyle(h_opts[cat]["linestyle"])
                    # call ax.legend() with the new values
                    ax.legend(handles,labels, loc="upper right")
        
                        # add lep category
                    ax.text(
                        0.02, 0.91, "%s\n%s" % (lep_cats[lepcat], jet_mults["4PJets"]),
                        fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                    )
                        ## set axes labels and titles
                    hep.cms.label(ax=ax, fontsize=rcParams["font.size"], data=False, year="2016APV+2016" if (combine_2016 and ("2016" in year)) else year, lumi=round(lumi_to_use, 1))
        
                    figname = os.path.join(pltdir, "_".join([year, jobid, "4PJets", lepcat, hname]))
                    fig.savefig(figname)
                    print(f"{figname} written")
                    plt.close()
        
            elif h_tot.dense_dim() == 2:
                xtitle, ytitle, x_rebinning, x_lims, y_rebinning, y_lims, save_dist, is3jet_dist = var_pars
    
                ## hists should have 3 category axes (jet multiplicity, lepton category, permutation category) followed by variable(s)
                for lepcat in sorted(set([key[0] for key in h_tot.values().keys()])):
                    pltdir = os.path.join(outdir, "4PJets", lepcat)
                    if not os.path.isdir(pltdir):
                        os.makedirs(pltdir)
        
                    hslice = h_tot[lepcat].integrate("lepcat")
                        ## make groups based on perm category
                    pcat_groups = perm_cats["4PJets"]["2D"]
                    hslice = hslice.group(pcat_cat, pcat, pcat_groups)
        
                    if x_rebinning != 1:
                        xaxis_name = hslice.dense_axes()[0].name
                        hslice = hslice.rebin(xaxis_name, x_rebinning)
                    if y_rebinning != 1:
                        yaxis_name = hslice.dense_axes()[1].name
                        hslice = hslice.rebin(yaxis_name, y_rebinning)
        
                        # make 1D projection along dense axes
                    for dax in range(2):
                        fig, ax = plt.subplots()
                        fig.subplots_adjust(hspace=.07)
        
                        #set_trace()
                        hproj = hslice.integrate(hslice.dense_axes()[dax].name)
                            ## plot comparison of perm categories
                        plot.plot1d(hproj,
                            overlay=hproj.axes()[0].name,
                            ax=ax,
                            clear=False,
                            line_opts={"linestyle" : "-"},
                            density=True, # normalized to 1,
                        )
                        ax.autoscale()#axis="x", tight=True)
                        ax.set_ylim(0, ax.get_ylim()[1]*1.15)
                        ax.set_ylabel("Probability Density")
                        ax.set_xlabel(ytitle) if dax == 0 else ax.set_xlabel(xtitle)
                        ax.set_xlim(y_lims) if dax == 0 else ax.set_xlim(x_lims)
        
                        #set_trace()
                           ## set legend and corresponding colors
                        h_opts ={key: hstyles[key] for key in pcat_groups.keys()}
                        handles, labels = ax.get_legend_handles_labels()
                        for idx, cat in enumerate(labels):
                            labels[idx] = h_opts[cat]["name"]
                            handles[idx].set_color(h_opts[cat]["color"])
                            handles[idx].set_linestyle(h_opts[cat]["linestyle"])
                        # call ax.legend() with the new values
                        ax.legend(handles,labels, loc="upper right")
        
                        #set_trace()
                            # add lep category
                        ax.text(
                            0.02, 0.91, "%s\n%s" % (lep_cats[lepcat], jet_mults["4PJets"]),
                            fontsize=rcParams["font.size"], 
                            horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                        )
        
                            ## set axes labels and titles
                        hep.cms.label(ax=ax, fontsize=rcParams["font.size"], data=False, year="2016APV+2016" if (combine_2016 and ("2016" in year)) else year, lumi=round(lumi_to_use, 1))
    
                        figname = os.path.join(pltdir, "_".join([year, jobid, "4PJets", lepcat, hname.split("_vs_")[0] if dax == 0 else hname.split("_vs_")[1]]))
                        fig.savefig(figname)
                        print(f"{figname} written")
                        plt.close()
        
                        # make 2D plots for each permutation category
                    for cat in sorted(set([key[0] for key in hslice.values().keys()])):            
                        hcat = hslice[cat].integrate("permcat")
    
                            ## normalized plots before interpolation
                        fig, ax = plt.subplots()
                        fig.subplots_adjust(hspace=.07)
        
                        norm_values = hcat.values()[()]/np.sum(hcat.values()[()]) ## normalized array of hist values
                        opts = {"cmap_label" : "$P_{M}$"}
                        Plotter.plot_2d_norm(xaxis_name=hcat.axes()[0].name, yaxis_name=hcat.axes()[1].name,
                            values=np.ma.masked_where(norm_values <= 0.0, norm_values), # mask nonzero probabilities for plotting
                            xlimits=x_lims, ylimits=y_lims, xlabel=xtitle, ylabel=ytitle,
                            hdict=hcat, ax=ax, **opts)
    
                            # add lep category
                        ax.text(
                            0.02, 0.91, "%s\n%s" % (lep_cats[lepcat], jet_mults["4PJets"]),
                            fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                        )
                            # add perm category
                        ax.text(
                            0.95, 0.95, hstyles[cat]["name"],
                            fontsize=rcParams["font.size"], horizontalalignment="right", verticalalignment="bottom", transform=ax.transAxes
                        )
                        hep.cms.label(ax=ax, fontsize=rcParams["font.size"], data=False, year="2016APV+2016" if (combine_2016 and ("2016" in year)) else year, lumi=round(lumi_to_use, 1))
        
                        figname = os.path.join(pltdir, "_".join([year, jobid, "4PJets", lepcat, hname, cat, "norm", "orig"]))
                        fig.savefig(figname)
                        print(f"{figname} written")
                        plt.close()
    
    
                            ## normalized plots after interpolation
                        fig, ax = plt.subplots()
                        fig.subplots_adjust(hspace=.07)
       
                        #set_trace() 
                            # get interpolated values (using bilinear) 
                        orig_xbins, orig_ybins, orig_vals = hcat.axes()[0].edges(), hcat.axes()[1].edges(), hcat.values()[()]
                        orig_xcenters = [(orig_xbins[idx]+orig_xbins[idx+1])/2 for idx in range(len(orig_xbins)-1)] 
                        orig_ycenters = [(orig_ybins[idx]+orig_ybins[idx+1])/2 for idx in range(len(orig_ybins)-1)] 
                        output_xbins = np.arange(min(orig_xbins), max(orig_xbins)+1, 1)
                        output_ybins = np.arange(min(orig_ybins), max(orig_ybins)+1, 1)
                        fit = interpolate.interp2d(orig_xcenters, orig_ycenters, orig_vals, kind="linear")
                        interped_array = fit(output_xbins, output_ybins)

                            # save distribution to dict
                        if save_dist and (lepcat == "Tight") and (cat == "Correct_THad"):
                            lookup = dense_lookup(*(interped_array, (output_xbins, output_ybins))) # not normalized
                            permProbs[year]["4PJets"].update({hname : lookup})
        
                        norm_intval = interped_array/np.sum(interped_array) ## normalized array of values after interpolating
                        values = np.ma.masked_where(norm_intval <= 0., norm_intval)
                        Plotter.plot_2d_norm(xbins=output_xbins, ybins=output_ybins,
                            values=np.ma.masked_where(values <= 0.0, values), # mask nonzero probabilities for plotting
                            xlimits=x_lims, ylimits=y_lims, xlabel=xtitle, ylabel=ytitle,
                            hdict=hcat, ax=ax, **opts)
    
                            # add lep category
                        ax.text(
                            0.02, 0.91, "%s\n%s" % (lep_cats[lepcat], jet_mults["4PJets"]),
                            fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                        )
                            # add perm category
                        ax.text(
                            0.95, 0.95, hstyles[cat]["name"],
                            fontsize=rcParams["font.size"], horizontalalignment="right", verticalalignment="bottom", transform=ax.transAxes
                        )
                        hep.cms.label(ax=ax, fontsize=rcParams["font.size"], data=False, year="2016APV+2016" if (combine_2016 and ("2016" in year)) else year, lumi=round(lumi_to_use, 1))
        
                        figname = os.path.join(pltdir, "_".join([year, jobid, "4PJets", lepcat, hname, cat, "norm", "interp"]))
                        fig.savefig(figname)
                        print(f"{figname} written")
                        plt.close()



    # write corrections to coffea file
if ((len(njets_to_run) == 2) and (len(years_to_run) == max_years)) or (args.force_save):
    corrdir = os.path.join(proj_dir, "Corrections", jobid)
    if not os.path.isdir(corrdir):
        os.makedirs(corrdir)

        # 2016 and 2016APV are the same if they"re computed together
    if combine_2016:
        permProbs[computed_combined_2016_year_to_copy] = permProbs[computed_combined_2016_year_key]

    permprobs_name = os.path.join(corrdir, f"prob_{jobid}.coffea")
    save(permProbs, permprobs_name)
    print("\n", permprobs_name, "written")
