#! /bin/env python

import time
tic = time.time()

from coffea.util import load, save
from pdb import set_trace
import os
import Utilities.plot_tools as plt_tools
from coffea import hist
import numpy as np
import fnmatch
import Utilities.Plotter as Plotter
import uproot3
import coffea.processor as processor    
import Utilities.systematics as systematics
import Utilities.final_analysis_binning as final_binning

base_jobid = os.environ["base_jobid"]

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("--only_bkg", action="store_true", help="Make background templates only.")
parser.add_argument("--only_sig", action="store_true", help="Make signal templates only.")
parser.add_argument("--maskData", action="store_false", help="Mask templates for data, default is True.")
parser.add_argument("--kfactors", action="store_true", help="Apply signal k-factors to signal")
parser.add_argument("--scale_mtop3gev", action="store_true", help="Scale 3GeV mtop variations by 1/6")
args = parser.parse_args()

if args.year == "2016APV":
    btag_reg_names_dict = {
        "Signal" : {"reg" : "btagPass"},
        "Down"   : {"reg" : "p00p20", "label" : "Down (0.0-0.2)", "color" : "b"},
        "Central": {"reg" : "p20p40", "label" : "Cen (0.2-0.4)", "color" : "k"},
        "Up"     : {"reg" : "p40p60", "label" : "Up (0.4-0.6)", "color" : "r"},
    }
if args.year == "2016":
    btag_reg_names_dict = {
        "Signal" : {"reg" : "btagPass"},
        "Down"   : {"reg" : "p00p1949", "label" : "Down (0.0-0.1949)", "color" : "b"},
        "Central": {"reg" : "p1949p3898", "label" : "Cen (0.1949-0.3898)", "color" : "k"},
        "Up"     : {"reg" : "p3898p5847", "label" : "Up (0.3898-0.5847)", "color" : "r"},
    }
if args.year == "2017":
    btag_reg_names_dict = {
        "Signal" : {"reg" : "btagPass"},
        "Down"   : {"reg" : "p00p1502", "label" : "Down (0.0-0.1502)", "color" : "b"},
        "Central": {"reg" : "p1502p3004", "label" : "Cen (0.1502-0.3004)", "color" : "k"},
        "Up"     : {"reg" : "p3004p4506", "label" : "Up (0.3004-0.4506)", "color" : "r"},
    }
if args.year == "2018":
    btag_reg_names_dict = {
        "Signal" : {"reg" : "btagPass"},
        "Down"   : {"reg" : "p00p1389", "label" : "Down (0.0-0.1389)", "color" : "b"},
        "Central": {"reg" : "p1389p2779", "label" : "Cen (0.1389-0.2779)", "color" : "k"},
        "Up"     : {"reg" : "p2779p4168", "label" : "Up (0.2779-0.4168)", "color" : "r"},
    }



def get_bkg_templates():
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """
    ## variables that only need to be defined/evaluated once
    templates_to_check = systematics.template_sys_to_name if args.scale_mtop3gev else systematics.combine_template_sys_to_name
    hdict = plt_tools.add_coffea_files(bkg_fnames) if len(bkg_fnames) > 1 else load(bkg_fnames[0])

        # get correct hist and rebin
    hname_to_use = "mtt_vs_tlep_ctstar_abs"
    if hname_to_use not in hdict.keys():
        raise ValueError(f"{hname_to_use} not found in file")
    xrebinning, yrebinning = linearize_binning
    histo = hdict[hname_to_use][Plotter.nonsignal_samples] # process, sys, jmult, leptype, btag, lepcat
    
    xaxis_name = histo.dense_axes()[0].name
    yaxis_name = histo.dense_axes()[1].name
        ## rebin x axis
    if isinstance(xrebinning, np.ndarray):
        new_xbins = hist.Bin(xaxis_name, xaxis_name, xrebinning)
    elif isinstance(xrebinning, float) or isinstance(xrebinning, int):
        new_xbins = xrebinning
    histo = histo.rebin(xaxis_name, new_xbins)
        ## rebin y axis
    if isinstance(yrebinning, np.ndarray):
        new_ybins = hist.Bin(yaxis_name, yaxis_name, yrebinning)
    elif isinstance(yrebinning, float) or isinstance(yrebinning, int):
        new_ybins = yrebinning
    rebin_histo = histo.rebin(yaxis_name, new_ybins)
    
        ## scale ttJets events, split by reconstruction type, by normal ttJets lumi correction
    ttJets_permcats = ["*right", "*matchable", "*unmatchable", "*sl_tau", "*other"]
    names = [dataset for dataset in sorted(set([key[0] for key in histo.values().keys()]))] # get dataset names in hists
    ttJets_cats = [name for name in names if any([fnmatch.fnmatch(name, cat) for cat in ttJets_permcats])] # gets ttJets(_PS)_other, ...

        ## make groups based on process
    process = hist.Cat("process", "Process", sorting="placement")
    process_cat = "dataset"

    histo_dict = processor.dict_accumulator({njets : {"Muon" : {}, "Electron" :{}} for njets in njets_to_run})

    for lep in ["Muon", "Electron"]:
        #set_trace()    
        ## make groups based on process
        process_groups = plt_tools.make_dataset_groups(lep, args.year, samples=names, gdict="templates")
        
        lumi_correction = lumi_corr_dict[args.year][f"{lep}s"]
                # scale ttJets events, split by reconstruction type, by normal ttJets lumi correction
        if len(ttJets_cats) > 0:
            for tt_cat in ttJets_cats:
                ttJets_lumi_topo = "_".join(tt_cat.split("_")[:-2]) if "sl_tau" in tt_cat else "_".join(tt_cat.split("_")[:-1]) # gets ttJets[SL, Had, DiLep] or ttJets_PS
                ttJets_eff_lumi = lumi_correction[ttJets_lumi_topo]
                lumi_correction.update({tt_cat: ttJets_eff_lumi})
    
        histo = rebin_histo.copy()
        histo.scale(lumi_correction, axis="dataset")
        for nom_tt in ["ttJetsSL", "ttJetsDiLep", "ttJetsHad"]:
            for cat in ttJets_permcats:
                # rescale LHEscale systematics correctly
                for hname, dname in tt_LHEscale_wts_name_dict.items():
                    if f"{nom_tt}_{cat[1:]}" not in ttJets_cats: continue
                    #print(f"{nom_tt}_{cat[1:]}_{dname}")
                    lhe_scale = lumi_correction[f"{nom_tt}_{dname}"]/lumi_correction[f"{nom_tt}_{cat[1:]}"]
                    histo.scale({(f"{nom_tt}_{cat[1:]}", hname) : lhe_scale}, axis=("dataset", "sys"))

        histo = histo.group(process_cat, process, process_groups)[:, :, :, lep, :, :].integrate("leptype")

        #set_trace()
        systs = sorted(set([key[1] for key in histo.values().keys()]))
        systs.insert(0, systs.pop(systs.index("nosys"))) # move "nosys" to the front

            # loop over each jet multiplicity
        for jmult in njets_to_run:
                # get sideband and signal region hists
            cen_sb_histo = Plotter.linearize_hist(histo[:, "nosys", jmult, btag_reg_names_dict["Central"]["reg"]].integrate("jmult").integrate("btag").integrate("sys"))
            up_sb_histo = Plotter.linearize_hist(histo[:, "nosys", jmult, btag_reg_names_dict["Up"]["reg"]].integrate("jmult").integrate("btag").integrate("sys"))
            dw_sb_histo = Plotter.linearize_hist(histo[:, "nosys", jmult, btag_reg_names_dict["Down"]["reg"]].integrate("jmult").integrate("btag").integrate("sys"))
            sig_histo = Plotter.linearize_hist(histo[:, :, jmult, btag_reg_names_dict["Signal"]["reg"]].integrate("jmult").integrate("btag"))

                # loop over each systematic
            for sys in systs:
                #set_trace()
                if sys not in templates_to_check[args.year].keys(): continue

                    # EWQCD background estimation only needed for 'nosys'
                if "data_obs" not in sorted(set([key[0] for key in sig_histo.values().keys()])):
                    sys_histo = sig_histo[:, sys].integrate("sys")
                else:
                    sys_histo = Plotter.BKG_Est(sig_reg=sig_histo[:, sys].integrate("sys"), sb_reg=cen_sb_histo, norm_type="SigMC", sys=sys, ignore_uncs=True, isForTemplates=True) if sys == "nosys" else sig_histo[:, sys].integrate("sys")

                    ## write nominal and systematic variations for each topology to file
                for proc in sorted(set([key[0] for key in sys_histo.values().keys()])):
                    if ("TT" not in proc) and (sys in systematics.ttJets_sys.values()): continue
                    if (proc == "data_obs") and not (sys == "nosys"): continue
                    if not sys_histo[proc].values().keys():
                        print(f"Systematic {sys} for {lep} {jmult} {proc} not found, skipping")
                        continue
                    if (proc == "EWQCD"):
                        if sys != "nosys": continue
                        print(args.year, lep, jmult, sys, proc)

                            # get original nominal distribution
                        template_histo = sys_histo[proc].integrate("process")

                        # get shape variations from btag sb regions
                            # Up region
                        bkg_shapeUp = Plotter.data_minus_top(up_sb_histo, isForTemplates=True) # find data - (tt+st) distributions
                        bkg_shapeUp.scale(np.sum(template_histo.values()[()])/np.sum(bkg_shapeUp.values()[()])) # normalize to mc yield in signal region
                        print(args.year, lep, jmult, "shapeUp", proc)
                            # Down region
                        bkg_shapeDown = Plotter.data_minus_top(dw_sb_histo, isForTemplates=True) # find data - (tt+st) distributions
                        bkg_shapeDown.scale(np.sum(template_histo.values()[()])/np.sum(bkg_shapeDown.values()[()])) # normalize to mc yield in signal region
                        print(args.year, lep, jmult, "shapeDown", proc)

                            ## save template histos to coffea dict
                        histo_dict[jmult][lep][f"{proc}_{sys}"] = template_histo.copy()
                        histo_dict[jmult][lep][f"{proc}_shapeUp"] = bkg_shapeUp.copy()
                        histo_dict[jmult][lep][f"{proc}_shapeDown"] = bkg_shapeDown.copy()
                    else:
                        print(args.year, lep, jmult, sys, proc)
                        #if "EWcorr" in sys: set_trace()
                        template_histo = sys_histo[proc].integrate("process")

                            # scale relative deviation for mtop3GeV by 1/6
                        if ("MTOP3GEV" in templates_to_check[args.year][sys].upper()) and (args.scale_mtop3gev):
                            #set_trace()
                            nominal_histo = histo_dict[jmult][lep][f"{proc}_nosys"].copy()
                            MTOP3GEV_scaled_var_yields = nominal_histo.values()[()] + (template_histo.values()[()] - nominal_histo.values()[()]) * 1./6.
                            template_histo = Plotter.np_array_TO_hist(sumw=MTOP3GEV_scaled_var_yields, sumw2=np.zeros(MTOP3GEV_scaled_var_yields.size), hist_template=nominal_histo)

                        #set_trace()
                            ## save template histos to coffea dict
                        histo_dict[jmult][lep][f"{proc}_{sys}"] = template_histo.copy()

    #set_trace()
    coffea_out = os.path.join(outdir, f"raw_templates_lj_bkg_mtopscaled_{args.year}_{jobid}.coffea" if args.scale_mtop3gev else f"raw_templates_lj_bkg_{args.year}_{jobid}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")


def get_sig_templates():
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """
    ## variables that only need to be defined/evaluated once
    templates_to_check = systematics.template_sys_to_name if args.kfactors else systematics.combine_template_sys_to_name

    hdict = plt_tools.add_coffea_files(sig_fnames) if len(sig_fnames) > 1 else load(sig_fnames[0])

        # get correct hist and rebin
    hname_to_use = "mtt_vs_tlep_ctstar_abs"
    if hname_to_use not in hdict.keys():
        raise ValueError(f"{hname_to_use} not found in file")
    xrebinning, yrebinning = linearize_binning
    histo = hdict[hname_to_use] # process, sys, jmult, leptype, btag, lepcat

    #set_trace()    
    xaxis_name = histo.dense_axes()[0].name
    yaxis_name = histo.dense_axes()[1].name
        ## rebin x axis
    if isinstance(xrebinning, np.ndarray):
        new_xbins = hist.Bin(xaxis_name, xaxis_name, xrebinning)
    elif isinstance(xrebinning, float) or isinstance(xrebinning, int):
        new_xbins = xrebinning
    histo = histo.rebin(xaxis_name, new_xbins)

        ## rebin y axis
    if isinstance(yrebinning, np.ndarray):
        new_ybins = hist.Bin(yaxis_name, yaxis_name, yrebinning)
    elif isinstance(yrebinning, float) or isinstance(yrebinning, int):
        new_ybins = yrebinning
    histo = histo.rebin(yaxis_name, new_ybins)
    rebin_histo = histo[Plotter.signal_samples, :, :, :, "btagPass"].integrate("btag")

    names = [dataset for dataset in sorted(set([key[0] for key in rebin_histo.values().keys()]))] # get dataset names in hists

    signals = sorted(set([key[0] for key in rebin_histo.values().keys()]))    
    signals = [sig for sig in signals if "TTJetsSL" in sig] # only use SL decays

    systs = sorted(set([key[1] for key in rebin_histo.values().keys()]))
    systs.insert(0, systs.pop(systs.index("nosys"))) # move "nosys" to the front

    histo_dict = processor.dict_accumulator({njets : {"Muon" : {}, "Electron" :{}} for njets in njets_to_run})

    #set_trace()
        # write signal dists to temp file        
    for lep in ["Muon", "Electron"]:
            # scale by lumi
        lumi_correction = lumi_corr_dict[args.year][f"{lep}s"]
        histo = rebin_histo.copy()
        histo.scale(lumi_correction, axis="dataset")
        process_groups = plt_tools.make_dataset_groups(lep, args.year, samples=names, gdict="templates")
        histo = histo.group("dataset", hist.Cat("process", "Process", sorting="placement"), process_groups)
    
        for jmult in njets_to_run:
            #set_trace()
            for signal in signals:
                #set_trace()
                for sys in systs:
                    if sys not in templates_to_check[args.year].keys(): continue

                    print(args.year, lep, jmult, signal, sys)
                    template_histo = Plotter.linearize_hist(histo[signal, sys, jmult, lep].integrate("jmult").integrate("leptype").integrate("process").integrate("sys"))
                    if ("RENORM" in sys.upper()) or ("FACTOR" in sys.upper()):
                        #set_trace()
                        lhe_scale = lumi_correction[f"{signal}_{signal_LHEscale_wts_name_dict[sys]}"]/lumi_correction[signal]
                        template_histo.scale(lhe_scale)

                        ## save template histos to coffea dict
                    histo_dict[jmult][lep][f"{signal}_{sys}"] = template_histo.copy()

    coffea_out = os.path.join(outdir, f"raw_templates_lj_sig_kfactors_{args.year}_{jobid}.coffea" if args.kfactors else  f"raw_templates_lj_sig_{args.year}_{jobid}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")


if __name__ == "__main__":
    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]

    njets_to_run = ["3Jets", "4PJets"]

        ## initialize lumi scaling files 
    lumi_corr_dict = load(os.path.join(proj_dir, "Corrections", base_jobid, "MC_LumiWeights_kfactors.coffea" if args.kfactors else "MC_LumiWeights.coffea"))

    #set_trace()
        # define variables to get histogram for signal
    sig_analyzer = "htt_btag_sb_regions"
    sig_input_dir = os.path.join(proj_dir, "results", f"{args.year}_{jobid}", sig_analyzer)
    if os.path.isdir(sig_input_dir):
        sig_fnames = fnmatch.filter(os.listdir(sig_input_dir), "*SIG*TOT.coffea")
        sig_fnames = [os.path.join(sig_input_dir, fname) for fname in sig_fnames]
    else: print("No signal file found.")
    
        # define variables to get histogram for background    
    bkg_analyzer = "htt_btag_sb_regions"
    bkg_input_dir = os.path.join(proj_dir, "results", f"{args.year}_{jobid}", bkg_analyzer)
    if os.path.isdir(bkg_input_dir):
        bkg_fnames = fnmatch.filter(os.listdir(bkg_input_dir), "*BKG*TOT.coffea")
        bkg_fnames = [os.path.join(bkg_input_dir, fname) for fname in bkg_fnames]
    else: print("No background file found.")
    
    outdir = os.path.join(proj_dir, "results", f"{args.year}_{jobid}", f"Templates_{bkg_analyzer}")
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
  
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

    linearize_binning = (
        final_binning.mtt_binning,
        final_binning.ctstar_abs_binning
    ) 

    #set_trace()
    if not args.only_sig:
        print("Creating background templates")
        get_bkg_templates()

    if not args.only_bkg:
        print("Creating signal templates")
        #set_trace()
        get_sig_templates()

    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
