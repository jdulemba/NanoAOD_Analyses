# matplotlib
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend("agg")
from matplotlib import rcParams
rcParams["font.size"] = 20
rcParams["savefig.format"] = "png"
rcParams["savefig.bbox"] = "tight"

from coffea.hist import plot
from coffea.util import load
from pdb import set_trace
import os
import Utilities.plot_tools as plt_tools
import Utilities.prettyjson as prettyjson
import re
from coffea import hist
import numpy as np
import fnmatch
import Utilities.Plotter as Plotter
import Utilities.final_analysis_binning as final_binning
import uproot3
import Utilities.HistExport as HistExport
import Utilities.systematics as systematics

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("topology", choices=["TT", "Signal"], help="Choose the topology to make plots for")
parser.add_argument("plots", choices=["RECO", "GEN", "All"], help="Choose which hists to make plots for")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
base_jobid = os.environ["base_jobid"]
jobid = os.environ["jobid"]
analyzer = "ResponseMatrix"

input_dir = os.path.join(proj_dir, "results", f"{args.year}_{jobid}", analyzer)
#outdir = os.path.join(proj_dir, "plots", f"{args.year}_{jobid}", analyzer, "MTOPcut")
#outdir = os.path.join(proj_dir, "plots", f"{args.year}_{jobid}", analyzer)
#outdir = os.path.join(proj_dir, "plots", f"{args.year}_{jobid}", analyzer, "GenKinCuts")
#outdir = os.path.join(proj_dir, "plots", f"{args.year}_{jobid}", analyzer, "NoCuts")
outdir = os.path.join(proj_dir, "plots", f"{args.year}_{jobid}", analyzer, "NoCuts_HigherST")
if not os.path.isdir(outdir):
    os.makedirs(outdir)

jet_mults = {
    "3Jets" : "3 jets",
    "4PJets" : "4+ jets"
}

objtypes = {
    "Lep" :  {
        "Muon" : "$\\mu$",
        "Electron" : "$e$",
    }
}

reco_variables = {
    "Reco_mtt_x_tlep_ctstar_abs_x_st_inds" : "$S_T$(tops) $\otimes$ $m_{t\\bar{t}}$ $\otimes$ |cos($\\theta^{*}_{t_{l}}$)|",
}
if args.topology == "TT":
    reco_variables["Reco_VS_Gen_mtt_x_tlep_ctstar_abs_x_st_inds"] = ("Gen $S_T$(tops) $\otimes$ $m_{t\\bar{t}}$ $\otimes$ |cos($\\theta^{*}_{t_{l}}$)|", "Reco")

gen_variables = {
    "Gen_mtt_x_tlep_ctstar_abs_x_st_inds" : "$S_T$(tops) $\otimes$ $m_{t\\bar{t}}$ $\otimes$ |cos($\\theta^{*}_{t_{l}}$)|",
}

signal_LHEscale_wts_name_dict = {
    "AH_FACTORDown" : "uF_down",
    "AH_FACTORUp"   : "uF_up",
    "AH_RENORMDown" : "uR_down",
    "AH_RENORMUp"   : "uR_up",
}

tt_LHEscale_wts_name_dict = {
    "FACTORDown" : "uF_down",
    "FACTORUp"   : "uF_up",
    "RENORMDown" : "uR_down",
    "RENORMUp"   : "uR_up",
}
widthTOname = lambda width : str(width).replace(".", "p")

    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())[args.year]
lumi_correction = load(os.path.join(proj_dir, "Corrections", base_jobid, "MC_LumiWeights_kfactors.coffea"))[args.year]

if args.plots == "All":
    #root_fname = os.path.join(outdir, f"ResponseMatrix_{args.topology}_MTOPcut_nosys_{args.year}.root")
    #root_fname = os.path.join(outdir, f"ResponseMatrix_{args.topology}_nosys_{args.year}.root")
    #root_fname = os.path.join(outdir, f"ResponseMatrix_{args.topology}_NoCuts_{args.year}.root")
    root_fname = os.path.join(outdir, f"ResponseMatrix_{args.topology}_NoCuts_HigherST_{args.year}.root")
    #root_fname = os.path.join(outdir, f"ResponseMatrix_{args.topology}_GenKinCuts_{args.year}.root")
    #root_fname = os.path.join(outdir, f"ResponseMatrix_{args.topology}_{args.year}.root")
    fout = uproot3.recreate(root_fname, compression=uproot3.ZLIB(4)) if os.path.isfile(root_fname) else uproot3.create(root_fname)


    ## make gen plots
if (args.plots == "GEN") or (args.plots == "All"):
    #set_trace()
    gen_fname_fnmatch = f"*GenLevel*{args.topology}*NoCuts_HigherST*TOT.coffea"
    #gen_fname_fnmatch = f"*GenLevel*{args.topology}*NoCuts*TOT.coffea"
    #gen_fname_fnmatch = f"*GenLevel*{args.topology}*KinCuts*TOT.coffea"
    #gen_fname_fnmatch = f"*GenLevel*{args.topology}*TOT.coffea"
    gen_fnames = fnmatch.filter(os.listdir(input_dir), gen_fname_fnmatch)
    gen_fnames = [os.path.join(input_dir, fname) for fname in gen_fnames]
    if len(gen_fnames) > 1: raise ValueError("More than 1 file found for GenLevel")
    gen_hdict = load(gen_fnames[0])

        ## get data lumi and scale MC by lumi
    lumi_to_use = (data_lumi_year["Muons"]+data_lumi_year["Electrons"])/2000.

    for hname in gen_variables.keys():
        if hname not in gen_hdict.keys():
            raise ValueError(f"{hname} not found in file")
        histo = gen_hdict[hname].copy()
        procs = sorted(set([key[0] for key in histo.values().keys()]))
        for proc in procs:
            pltdir = os.path.join(outdir, "Gen", proc)
            if not os.path.isdir(pltdir):
                os.makedirs(pltdir)

            hproc = (histo[proc, :].integrate("dataset")).copy()
            hproc.scale((lumi_correction["Muons"][proc]+lumi_correction["Electrons"][proc])/2)
    
            if hproc.dense_dim() == 1:
                xtitle = gen_variables[hname]
                bin_edges = hproc.dense_axes()[0].edges()   

                for sys in sorted(set([key[0] for key in hproc.values().keys()])):
                    #if sys != "nosys": continue
                    if "RENORM_FACTOR" in sys: continue
                    sys_histo = (hproc[sys].integrate("sys")).copy()

                    # rescale LHEscale systematics correctly
                    if args.topology == "TT":
                        if sys in tt_LHEscale_wts_name_dict.keys():
                            lhe_scale = (lumi_correction["Muons"][f"{proc}_%s" % tt_LHEscale_wts_name_dict[sys]]+lumi_correction["Muons"][f"{proc}_%s" % tt_LHEscale_wts_name_dict[sys]])/(lumi_correction["Muons"][proc]+lumi_correction["Electrons"][proc])
                            sys_histo.scale(lhe_scale)
                    else:
                        if sys in signal_LHEscale_wts_name_dict.keys():
                            lhe_scale = (lumi_correction["Muons"][f"{proc}_{signal_LHEscale_wts_name_dict[sys]}"]+lumi_correction["Electrons"][f"{proc}_{signal_LHEscale_wts_name_dict[sys]}"])/(lumi_correction["Muons"][proc]+lumi_correction["Electrons"][proc])
                            sys_histo.scale(lhe_scale)

                    if proc == "ttJetsSL":
                        sysname = "TT" if sys == "nosys" else f"TT_{systematics.template_sys_to_name[args.year][sys]}"
                    else:
                        if "Int" in proc:
                            boson, mass, width, pI, wt = tuple(proc.split("_"))
                        else:
                            boson, mass, width, pI = tuple(proc.split("_"))
                        sub_name = "_".join(["%s%s" % (boson[0], mass[1:]), "relw%s" % widthTOname(width).split("W")[-1], pI.lower(), wt]) if pI == "Int" else "_".join(["%s%s" % (boson[0], mass[1:]), "relw%s" % widthTOname(width).split("W")[-1], pI.lower()])

                        sysname = sub_name if sys == "nosys" else f"{sub_name}_{systematics.template_sys_to_name[args.year][sys]}"

                    if "LEP" in sysname: sysname = sysname.replace("LEP", lep.lower())
                    if args.plots == "All":
                        fout[f"Gen_{sysname}"] = HistExport.export1d(sys_histo)

                    # plots
                    fig, ax = plt.subplots()
                    fig.subplots_adjust(hspace=.07)

                    Plotter.plot_1D(sys_histo.values()[()], bin_edges, xlabel=xtitle, ax=ax, label=sysname)
                            
                    ax.legend(loc="upper right")
                    ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.15)
    
                        # add lepton/jet multiplicity label
                    ax.text(
                        0.02, 0.92, "parton level",
                        fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                    )

                    ## draw vertical lines for distinguishing different ctstar bins
                    ctstar_lines = [idx*(final_binning.mtt_binning.size-1) for idx in range(1, (final_binning.st_binning.size-1)*(final_binning.ctstar_abs_binning.size-1))]
                    for ct_line in ctstar_lines:
                        ax.axvline(ct_line, color="k", linestyle="--")
                    ## draw vertical lines for distinguishing different S_T bins
                    st_lines = [idx*(final_binning.mtt_binning.size-1)*(final_binning.ctstar_abs_binning.size-1) for idx in range(1, final_binning.st_binning.size-1)]
                    for st_line in st_lines:
                        ax.axvline(st_line, color="k", linestyle="-", linewidth=2)
                    hep.cms.label(ax=ax, data=False, year=args.year, lumi=round(lumi_to_use, 1))
    
                    #set_trace()
                    figname = os.path.join(pltdir, f"{hname}_{sysname}")
                    fig.savefig(figname)
                    print(f"{figname} written")
                    plt.close()


    ## make reco plots
if (args.plots == "RECO") or (args.plots == "All"):
    #set_trace()
    reco_fname_fnmatch = f"*RecoLevel*{args.topology}*NoCuts_HigherST*TOT.coffea"
    #reco_fname_fnmatch = f"*RecoLevel*{args.topology}*NoCuts*TOT.coffea"
    #reco_fname_fnmatch = f"*RecoLevel*{args.topology}*KinCuts*TOT.coffea"
    #reco_fname_fnmatch = f"*RecoLevel*{args.topology}*TOT.coffea"
    #reco_fname_fnmatch = f"*RecoLevel*{args.topology}*MTOPcut*TOT.coffea"
    reco_fnames = fnmatch.filter(os.listdir(input_dir), reco_fname_fnmatch)
    reco_fnames = [os.path.join(input_dir, fname) for fname in reco_fnames]
    if len(reco_fnames) > 1: raise ValueError("More than 1 file found for RecoLevel")
    reco_hdict = load(reco_fnames[0])

    for hname in reco_variables.keys():
        if hname not in reco_hdict.keys():
            raise ValueError(f"{hname} not found in file")
        histo = reco_hdict[hname].copy()
        procs = sorted(set([key[0] for key in histo.values().keys()]))
        for proc in procs:
            hproc = (histo[proc, :].integrate("dataset")).copy()

            for lep in sorted(set([key[2] for key in hproc.values().keys()])):
                lep_histo = hproc.copy()
                lep_histo.scale(lumi_correction[f"{lep}s"][proc])

                for jmult in sorted(set([key[1] for key in lep_histo.values().keys()])):
                    hslice = (lep_histo[:, jmult, lep].integrate("jmult").integrate("leptype")).copy()
    
                    if hslice.dense_dim() == 1:
                        pltdir = os.path.join(outdir, lep, jmult, "RECO", proc)
                        if not os.path.isdir(pltdir):
                            os.makedirs(pltdir)

                        xtitle = reco_variables[hname]
                        bin_edges = hslice.dense_axes()[0].edges()   

                        #set_trace()
                        for sys in sorted(set([key[0] for key in hslice.values().keys()])):
                            #if sys != "nosys": continue
                            if "RENORM_FACTOR" in sys: continue
                            print(lep, jmult, sys)
                            sys_histo = (hslice[sys].integrate("sys")).copy()

                            # rescale LHEscale systematics correctly
                            if args.topology == "TT":
                                if sys in tt_LHEscale_wts_name_dict.keys():
                                    lhe_scale = lumi_correction[f"{lep}s"][f"{proc}_{tt_LHEscale_wts_name_dict[sys]}"]/lumi_correction[f"{lep}s"][proc]
                                    sys_histo.scale(lhe_scale)
                            else:
                                if sys in signal_LHEscale_wts_name_dict.keys():
                                    lhe_scale = lumi_correction[f"{lep}s"][f"{proc}_{signal_LHEscale_wts_name_dict[sys]}"]/lumi_correction[f"{lep}s"][proc]
                                    sys_histo.scale(lhe_scale)

                            if proc == "ttJetsSL":
                                sysname = "TT" if sys == "nosys" else f"TT_{systematics.template_sys_to_name[args.year][sys]}"
                            else:
                                if "Int" in proc:
                                    boson, mass, width, pI, wt = tuple(proc.split("_"))
                                else:
                                    boson, mass, width, pI = tuple(proc.split("_"))
                                sub_name = "_".join(["%s%s" % (boson[0], mass[1:]), "relw%s" % widthTOname(width).split("W")[-1], pI.lower(), wt]) if pI == "Int" else "_".join(["%s%s" % (boson[0], mass[1:]), "relw%s" % widthTOname(width).split("W")[-1], pI.lower()])
                                sysname = sub_name if sys == "nosys" else f"{sub_name}_{systematics.template_sys_to_name[args.year][sys]}"

                            if "LEP" in sysname: sysname = sysname.replace("LEP", lep.lower())
                            if args.plots == "All":
                                fout[f"Reco_{lep}_{jmult}_{sysname}"] = HistExport.export1d(sys_histo)

                            # plots
                            fig, ax = plt.subplots()
                            fig.subplots_adjust(hspace=.07)
   
                            Plotter.plot_1D(sys_histo.values()[()], bin_edges, xlabel=xtitle, ax=ax, label=sysname)
                                    
                            ax.legend(loc="upper right")
                            #ax.set_ylim(0, ax.get_ylim()[1]*1.15)
                            ax.set_ylim(ax.get_ylim()[0], ax.get_ylim()[1]*1.15)
    
                                # add lepton/jet multiplicity label
                            ax.text(
                                0.02, 0.88, "%s, %s\nparticle level" % (objtypes["Lep"][lep], jet_mults[jmult]),
                                fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                            )

                            ## draw vertical lines for distinguishing different ctstar bins
                            ctstar_lines = [idx*(final_binning.mtt_binning.size-1) for idx in range(1, (final_binning.st_binning.size-1)*(final_binning.ctstar_abs_binning.size-1))]
                            for ct_line in ctstar_lines:
                                ax.axvline(ct_line, color="k", linestyle="--")
                            ## draw vertical lines for distinguishing different S_T bins
                            st_lines = [idx*(final_binning.mtt_binning.size-1)*(final_binning.ctstar_abs_binning.size-1) for idx in range(1, final_binning.st_binning.size-1)]
                            for st_line in st_lines:
                                ax.axvline(st_line, color="k", linestyle="-", linewidth=2)
                            hep.cms.label(ax=ax, data=False, year=args.year, lumi=round(data_lumi_year[f"{lep}s"]/1000., 1))
    
                            figname = os.path.join(pltdir, f"{hname}_{lep}_{jmult}_{sysname}")
                            fig.savefig(figname)
                            print(f"{figname} written")
                            plt.close()

                    elif hslice.dense_dim() == 2:
                        pltdir = os.path.join(outdir, lep, jmult, "GEN_vs_RECO", proc)
                        if not os.path.isdir(pltdir):
                            os.makedirs(pltdir)
                        xtitle, ytitle = reco_variables[hname]
                        bin_edges = hslice.dense_axes()[0].edges()   

                        for sys in sorted(set([key[0] for key in hslice.values().keys()])):
                            if "RENORM_FACTOR" in sys: continue
                            print(lep, jmult, sys)
                            sys_histo = (hslice[sys].integrate("sys")).copy()

                            # rescale LHEscale systematics correctly
                            if args.topology == "TT":
                                if sys in tt_LHEscale_wts_name_dict.keys():
                                    lhe_scale = lumi_correction[f"{lep}s"][f"{proc}_{tt_LHEscale_wts_name_dict[sys]}"]/lumi_correction[f"{lep}s"][proc]
                                    sys_histo.scale(lhe_scale)
                            else:
                                raise ValueError(f"{hname} should not exist for signal")

                            if proc == "ttJetsSL":
                                sysname = "TT" if sys == "nosys" else f"TT_{systematics.template_sys_to_name[args.year][sys]}"
                            else:
                                raise ValueError(f"{hname} should not exist for signal")
                            if "LEP" in sysname: sysname = sysname.replace("LEP", lep.lower())

                            if args.plots == "All":
                                fout[f"Gen_vs_Reco_{lep}_{jmult}_{sysname}"] = HistExport.export2d(sys_histo)

                            # plots
                            fig, ax = plt.subplots()
                            fig.subplots_adjust(hspace=.07)
   
                            Plotter.plot_2d_norm(xaxis_name=sys_histo.axes()[0].name, yaxis_name=sys_histo.axes()[1].name,
                                values=np.ma.masked_where(sys_histo.values()[()] <= 0.0, sys_histo.values()[()]), # mask nonzero probabilities for plotting
                                xlimits=(sys_histo.axes()[0].edges()[0], sys_histo.axes()[0].edges()[-1]),
                                ylimits=(sys_histo.axes()[1].edges()[0], sys_histo.axes()[1].edges()[-1]), xlabel=xtitle, ylabel=ytitle,
                                hdict=sys_histo, ax=ax)
                                    
                            ax.legend(loc="upper right", title=sysname)
                                # add lepton/jet multiplicity label
                            ax.text(
                                0.02, 0.88, "%s, %s, $\mathrm{t\\bart}_{\ell j}$" % (objtypes["Lep"][lep], jet_mults[jmult]),
                                fontsize=rcParams["font.size"], horizontalalignment="left", verticalalignment="bottom", transform=ax.transAxes
                            )
                            hep.cms.label(ax=ax, data=False, year=args.year, lumi=round(data_lumi_year[f"{lep}s"]/1000., 1))
    
                            figname = os.path.join(pltdir, f"{hname}_{lep}_{jmult}_{sys}")
                            fig.savefig(figname)
                            print(f"{figname} written")
                            plt.close()


if args.plots == "All":
    fout.close()
    print(f"{root_fname} written")
