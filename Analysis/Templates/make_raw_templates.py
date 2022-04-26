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
import coffea.processor as processor    
import Utilities.systematics as systematics
import Utilities.final_analysis_binning as final_binning
import Utilities.btag_sideband_regions as btag_sidebands
import Utilities.prettyjson as prettyjson

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("templates_to_run", type=str, help="Choose which type of templates to run, multiple options can be input as ':' separated strings.")
parser.add_argument("--maskData", action="store_false", help="Mask templates for data, default is True.")
parser.add_argument("--kfactors", action="store_true", help="Apply signal k-factors to signal")
parser.add_argument("--scale_mtop3gev", action="store_true", help="Scale 3GeV mtop variations by 1/6")
args = parser.parse_args()


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

    #set_trace()    
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
                        #set_trace()
                        template_histo = sys_histo[proc].integrate("process")
                        sumw, sumw2 = template_histo.values(sumw2=True)[()]
                        rel_err = np.sqrt(sumw2)/np.abs(sumw)
                        rel_err_mask = rel_err > 10
                        if np.any(rel_err_mask):
                            #set_trace()
                            if (sys == "nosys"):
                                print(f"\tRelative error > 10 for this process! Setting bin {np.where(rel_err_mask)[0]} to 0")
                                sumw[rel_err_mask], sumw2[rel_err_mask] = 0., 0.
                            else:
                            # check if nosys hist has same behavior (if sys isn't nosys)
                                #set_trace()
                                nosys_template = histo_dict[jmult][lep][f"{proc}_nosys"].copy()
                                nosys_sumw, nosys_sumw2 = nosys_template.values(sumw2=True)[()]
                                nosys_rel_err_mask = (nosys_sumw == 0.) & (nosys_sumw2 == 0.)
                                if np.any(rel_err_mask & nosys_rel_err_mask):
                                    print(f"\tRelative error > 10 for this process! Setting bin {np.where(rel_err_mask & nosys_rel_err_mask)[0]} to 0")
                                    sumw[rel_err_mask & nosys_rel_err_mask], sumw2[rel_err_mask & nosys_rel_err_mask] = 0., 0.

                            if proc != "TB": set_trace()

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

    #set_trace()
    signals = sorted(plt_tools.create_sig_groups("MC_combined")[args.year].keys()) # get combined signal to loop over

    systs = sorted(set([key[1] for key in rebin_histo.values().keys()]))
    systs.insert(0, systs.pop(systs.index("nosys"))) # move "nosys" to the front

    histo_dict = processor.dict_accumulator({njets : {"Muon" : {}, "Electron" :{}} for njets in njets_to_run})

        # write signal dists to temp file        
    for lep in ["Muon", "Electron"]:
            # scale by lumi
        lumi_correction = lumi_corr_dict[args.year][f"{lep}s"]
        histo = rebin_histo.copy()
        histo = histo[:, :, :, lep].integrate("leptype")
        histo.scale(lumi_correction, axis="dataset")
        process_groups = plt_tools.make_dataset_groups(lep, args.year, samples=names, bkgdict="templates", sigdict="MC_indiv")
        histo = histo.group("dataset", hist.Cat("process", "Process", sorting="placement"), process_groups)

        for jmult in njets_to_run:
            for signal in signals:
                for sys in systs:
                    if sys not in templates_to_check[args.year].keys(): continue
                    #set_trace()

                    sl_template_histo = Plotter.linearize_hist(histo[signal.replace("TT", "TTJetsSL"), sys, jmult].integrate("jmult").integrate("process").integrate("sys"))
                    if not sl_template_histo.values().keys(): continue
                    dl_template_histo = Plotter.linearize_hist(histo[signal.replace("TT", "TTJetsDiLep"), sys, jmult].integrate("jmult").integrate("process").integrate("sys"))
                    if not dl_template_histo.values().keys(): continue
                    print(args.year, lep, jmult, signal, sys)
                    if ("RENORM" in sys.upper()) or ("FACTOR" in sys.upper()):
                        #set_trace()
                        sl_lhe_scale = lumi_correction[f"{signal.replace('TT', 'TTJetsSL')}_{signal_LHEscale_wts_name_dict[sys]}"]/lumi_correction[signal.replace("TT", "TTJetsSL")]
                        sl_template_histo.scale(sl_lhe_scale)
                        dl_lhe_scale = lumi_correction[f"{signal.replace('TT', 'TTJetsDiLep')}_{signal_LHEscale_wts_name_dict[sys]}"]/lumi_correction[signal.replace("TT", "TTJetsDiLep")]
                        dl_template_histo.scale(dl_lhe_scale)

                    template_histo = sl_template_histo.copy()
                    template_histo.add(dl_template_histo.copy())

                        ## save template histos to coffea dict
                    histo_dict[jmult][lep][f"{signal}_{sys}"] = template_histo.copy()

    coffea_out = os.path.join(outdir, f"raw_templates_lj_sig_kfactors_{args.year}_{jobid}.coffea" if args.kfactors else  f"raw_templates_lj_sig_{args.year}_{jobid}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")


def get_MEreweight_sig_templates():
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """
    ## variables that only need to be defined/evaluated once
    templates_to_check = systematics.template_sys_to_name if args.kfactors else systematics.combine_template_sys_to_name
    histo_dict = processor.dict_accumulator({njets : {"Muon" : {}, "Electron" :{}} for njets in njets_to_run})
    signal_kfactors = prettyjson.loads(open(os.path.join(proj_dir, "inputs", "signal_kfactors_ulkfactor_final_220129.json")).read())
    hname_to_use = "mtt_vs_tlep_ctstar_abs"

    allowed_masses = ["365", "380", "400", "425", "450", "475", "500", "525", "550", "575", "600", "625", "650", "675", "700", "725", "750", "775", "800", "825", "850", "875", "900", "925", "950", "975", "1000"]
    #allowed_widths = ["2p5", "5p0", "10p0", "25p0"]
    allowed_widths = ["2p5", "3p0", "4p0", "5p0", "8p0", "10p0", "25p0"]

    fname_dir = f"/eos/user/{os.environ['USER'][0]}/{os.environ['USER']}/NanoAOD_Analyses/results/{args.year}_{jobid}/signal_ME_evtReweighting/RecoLevel"
    fnames_to_run = [fname for fname in os.listdir(fname_dir) if fname.startswith("m")]
    if not fnames_to_run: raise ValueError(f"No file found in {fname_dir}")
    #set_trace()

    for fname in fnames_to_run:
        hdict = load(os.path.join(fname_dir, fname))

            # get correct hist and rebin
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
        rebin_histo = histo.rebin(yaxis_name, new_ybins)

        names = [dataset for dataset in sorted(set([key[0] for key in rebin_histo.values().keys()]))] # get dataset names in hists
        process_groups = plt_tools.make_dataset_groups("Muon", args.year, samples=names, bkgdict="templates", sigdict="MEreweight_combined")
        rebin_histo = rebin_histo.group("dataset", hist.Cat("process", "Process", sorting="placement"), process_groups)

        #set_trace()
        signals = sorted(set([key[0] for key in rebin_histo.values().keys()]))

        systs = sorted(set([key[1] for key in rebin_histo.values().keys()]))
        systs.insert(0, systs.pop(systs.index("nosys"))) # move "nosys" to the front

            # write signal dists to temp file        
        for lep in ["Muon", "Electron"]:
            histo = rebin_histo[:, :, :, lep].integrate("leptype")

            for jmult in njets_to_run:
                for signal in signals:
                        # only run over allowed masses and width points for now
                    mass, width = signal.split("_")[1], signal.split("_")[2]
                    if not ((mass.strip("M") in allowed_masses) and (width.strip("W") in allowed_widths)): continue
                    for sys in systs:
                        if sys not in templates_to_check[args.year].keys(): continue
                        #set_trace()

                        template_histo = Plotter.linearize_hist(histo[signal, sys, jmult].integrate("jmult").integrate("process").integrate("sys"))
                        if not template_histo.values().keys(): continue
                        print(args.year, lep, jmult, signal, sys)

                        if args.kfactors:
                            if ("RENORM" in sys.upper()) or ("FACTOR" in sys.upper()):
                                #set_trace()
                                if "Int_neg" in signal:
                                    kfactor = signal_kfactors[f"{signal.replace('TT', 'TTJetsSL').strip('_neg')}_{signal_LHEscale_wts_name_dict[sys]}"] # same as DiLep
                                if "Int_pos" in signal:
                                    kfactor = signal_kfactors[f"{signal.replace('TT', 'TTJetsSL').strip('_pos')}_{signal_LHEscale_wts_name_dict[sys]}"]
                                if "Res" in signal:
                                    kfactor = signal_kfactors[f"{signal.replace('TT', 'TTJetsSL')}_{signal_LHEscale_wts_name_dict[sys]}"]
                            else:
                                if "Int_neg" in signal:
                                    kfactor = signal_kfactors[signal.replace("TT", "TTJetsSL").strip("_neg")] # same as DiLep
                                if "Int_pos" in signal:
                                    kfactor = signal_kfactors[signal.replace("TT", "TTJetsSL").strip("_pos")] # same as DiLep
                                if "Res" in signal:
                                    kfactor = signal_kfactors[signal.replace("TT", "TTJetsSL")]
                            template_histo.scale(kfactor)

                            ## save template histos to coffea dict
                        histo_dict[jmult][lep][f"{signal}_{sys}"] = template_histo.copy()

    #set_trace()
    coffea_out = os.path.join(outdir, f"raw_templates_lj_MEsig_kfactors_{args.year}_{jobid}.coffea" if args.kfactors else  f"raw_templates_lj_MEsig_{args.year}_{jobid}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")


if __name__ == "__main__":
    allowed_template_options = ["bkg", "sig", "MEreweight_sig"]
    templates_to_run = (args.templates_to_run).split(":")
    templates_to_run = [template for template in templates_to_run if template in allowed_template_options]

    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]
    base_jobid = os.environ["base_jobid"]

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

    btag_reg_names_dict = btag_sidebands.btag_reg_names_dict[args.year]

    linearize_binning = (
        final_binning.mtt_binning,
        final_binning.ctstar_abs_binning
    ) 

    #set_trace()
    if "bkg" in templates_to_run:
        print("Creating background templates")
        get_bkg_templates()

    if "sig" in templates_to_run:
        print("Creating signal templates")
        get_sig_templates()

    if "MEreweight_sig" in templates_to_run:
        print("Creating ME reweighting signal templates")
        get_MEreweight_sig_templates()

    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
