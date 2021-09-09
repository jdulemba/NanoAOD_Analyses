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
from rootpy.io import root_open
import coffea.processor as processor    
import Utilities.systematics as systematics

base_jobid = os.environ["base_jobid"]

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016", "2017", "2018"] if base_jobid == "NanoAODv6" else ["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("--njets", default="all", nargs="?", choices=["3", "4+", "all"], help="Specify which jet multiplicity to use.")
parser.add_argument("--only_bkg", action="store_true", help="Make background templates only.")
parser.add_argument("--only_sig", action="store_true", help="Make signal templates only.")
parser.add_argument("--maskData", action="store_false", help="Mask templates for data, default is True.")
parser.add_argument("--naming", choices=["dir", "name"], default="name", help="Naming in output root files has separate names within one directory (name) or same naming in separate directories (dir). Default is name.")
args = parser.parse_args()

njets_to_run = []
if (args.njets == "3") or (args.njets == "all"):
    njets_to_run += ["3Jets"]
if (args.njets == "4+") or (args.njets == "all"):
    njets_to_run += ["4PJets"]


btag_reg_names_dict = {
    "Signal" : {"reg" : "btagPass"},
    "Central": {"reg" : "p15p30", "label" : "Cen (0.15-0.3)", "color" : "k"},
    "Up"     : {"reg" : "p30p45", "label" : "Up (0.3-0.45)", "color" : "r"},
    "Down"   : {"reg" : "p00p15", "label" : "Down (0.0-0.15)", "color" : "b"}
}


def get_bkg_templates(tmp_rname):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """
    ## variables that only need to be defined/evaluated once
    hdict = plt_tools.add_coffea_files(bkg_fnames) if len(bkg_fnames) > 1 else load(bkg_fnames[0])

        # get correct hist and rebin
    hname_to_use = "mtt_vs_tlep_ctstar_abs"
    if hname_to_use not in hdict.keys():
        raise ValueError("%s not found in file" % hname_to_use)
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

        # need to save coffea hist objects to file so they can be opened by uproot in the proper format
    upfout = uproot3.recreate(tmp_rname, compression=uproot3.ZLIB(4)) if os.path.isfile(tmp_rname) else uproot3.create(tmp_rname)

    if "3Jets" in njets_to_run:
        histo_dict_3j = processor.dict_accumulator({"Muon" : {}, "Electron" :{}})
    if "4PJets" in njets_to_run:
        histo_dict_4pj = processor.dict_accumulator({"Muon" : {}, "Electron" :{}})

    for lep in ["Muon", "Electron"]:
        orig_lepdir = "muNJETS" if lep == "Muon" else "eNJETS"

        #set_trace()    
        ## make groups based on process
        process_groups = plt_tools.make_dataset_groups(lep, args.year, samples=names, gdict="templates")
        #process_groups = plt_tools.make_dataset_groups(lep, args.year, samples=names, gdict="dataset")
        
        lumi_correction = lumi_corr_dict[args.year]["%ss" % lep]
                # scale ttJets events, split by reconstruction type, by normal ttJets lumi correction
        if len(ttJets_cats) > 0:
            for tt_cat in ttJets_cats:
                ttJets_lumi_topo = "_".join(tt_cat.split("_")[:-2]) if "sl_tau" in tt_cat else "_".join(tt_cat.split("_")[:-1]) # gets ttJets[SL, Had, DiLep] or ttJets_PS
                ttJets_eff_lumi = lumi_correction[ttJets_lumi_topo]
                lumi_correction.update({tt_cat: ttJets_eff_lumi})
    
        histo = rebin_histo.copy()
        histo.scale(lumi_correction, axis="dataset")
        histo = histo.group(process_cat, process, process_groups)[:, :, :, lep, :, :].integrate("leptype")

        #set_trace()
        systs = sorted(set([key[1] for key in histo.values().keys()]))
        systs.insert(0, systs.pop(systs.index("nosys"))) # move "nosys" to the front

            # loop over each jet multiplicity
        for jmult in njets_to_run:
            lepdir = orig_lepdir.replace("NJETS", jmult.lower())

                # get sideband and signal region hists
            cen_sb_histo = Plotter.linearize_hist(histo[:, "nosys", jmult, btag_reg_names_dict["Central"]["reg"]].integrate("jmult").integrate("btag").integrate("sys"))
            #up_sb_histo = histo[:, "nosys", jmult, btag_reg_names_dict["Up"]["reg"]].integrate("jmult").integrate("btag")
            #dw_sb_histo = histo[:, "nosys", jmult, btag_reg_names_dict["Down"]["reg"]].integrate("jmult").integrate("btag")
            sig_histo = Plotter.linearize_hist(histo[:, :, jmult, btag_reg_names_dict["Signal"]["reg"]].integrate("jmult").integrate("btag"))

                # loop over each systematic
            for sys in systs:
                if sys not in systematics.template_sys_to_name[args.year].keys(): continue

                sys_histo = sig_histo[:, sys].integrate("sys") if sys in systematics.ttJets_sys.values() else Plotter.BKG_Est(sig_reg=sig_histo[:, sys].integrate("sys"), sb_reg=cen_sb_histo, norm_type="SigMC", sys=sys, ignore_uncs=True)

                    ## write nominal and systematic variations for each topology to file
                #for proc in sorted(set([key[0] for key in sig_histo.values().keys()])):
                for proc in sorted(set([key[0] for key in sys_histo.values().keys()])):
                    if ("tt" not in proc) and (sys in systematics.ttJets_sys.values()): continue
                    #if (proc != "tt") and (sys in systematics.ttJets_sys.values()): continue
                    if (proc == "data_obs") and not (sys == "nosys"): continue
                    if not sys_histo[proc].values().keys():
                    #if not sig_histo[proc, sys].values().keys():
                        print(f"Systematic {sys} for {lep} {jmult} {proc} not found, skipping")
                        continue

                    print(args.year, lep, jmult, sys, proc)
                    #set_trace()
                    outhname = "_".join(list(filter(None, [proc, systematics.template_sys_to_name[args.year][sys][0], lepdir, (args.year)[-2:]])))
                    if "LEP" in outhname: outhname = outhname.replace("LEP", "muon") if lep == "Muon" else outhname.replace("LEP", "electron")

                    template_histo = sys_histo[proc].integrate("process")
                    #template_histo = sig_histo[proc, sys].integrate("process").integrate("sys")

                    #set_trace()
                        ## save template histos to coffea dict
                    if jmult == "3Jets":
                        histo_dict_3j[lep][f"{proc}_{sys}"] = template_histo.copy()
                    if jmult == "4PJets":
                        histo_dict_4pj[lep][f"{proc}_{sys}"] = template_histo.copy()

                        ## save template histo to root file
                    upfout[outhname] = hist.export1d(template_histo)

    if "3Jets" in njets_to_run:
        coffea_out_3j = os.path.join(outdir, f"test_raw_templates_lj_3Jets_bkg_{args.year}_{jobid}.coffea")
        save(histo_dict_3j, coffea_out_3j)
        print(f"{coffea_out_3j} written")
    if "4PJets" in njets_to_run:
        coffea_out_4pj = os.path.join(outdir, f"test_raw_templates_lj_4PJets_bkg_{args.year}_{jobid}.coffea")
        save(histo_dict_4pj, coffea_out_4pj)
        print(f"{coffea_out_4pj} written")

    upfout.close()
    print(f"{tmp_rname} written")


def get_sig_templates(tmp_rname):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """
    widthTOname = lambda width : str(width).replace(".", "p")
    nameTOwidth = lambda width : str(width).replace("p", ".")

    ## variables that only need to be defined/evaluated once
    hdict = plt_tools.add_coffea_files(sig_fnames) if len(sig_fnames) > 1 else load(sig_fnames[0])

        # get correct hist and rebin
    hname_to_use = "mtt_vs_tlep_ctstar_abs"
    if hname_to_use not in hdict.keys():
        raise ValueError(f"{hname_to_use} not found in file")
    xrebinning, yrebinning = linearize_binning
    #xrebinning, yrebinning = mtt_ctstar_2d_binning
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

        # need to save coffea hist objects to file so they can be opened by uproot in the proper format
    upfout = uproot3.recreate(tmp_rname, compression=uproot3.ZLIB(4)) if os.path.isfile(tmp_rname) else uproot3.create(tmp_rname)

    if "3Jets" in njets_to_run:
        histo_dict_3j = processor.dict_accumulator({"Muon" : {}, "Electron" :{}})
    if "4PJets" in njets_to_run:
        histo_dict_4pj = processor.dict_accumulator({"Muon" : {}, "Electron" :{}})

        # write signal dists to temp file        
    for lep in ["Muon", "Electron"]:
        orig_lepdir = "muNJETS" if lep == "Muon" else "eNJETS"

            # scale by lumi
        lumi_correction = lumi_corr_dict[args.year]["%ss" % lep]
        histo = rebin_histo.copy()
        histo.scale(lumi_correction, axis="dataset")
        process_groups = plt_tools.make_dataset_groups(lep, args.year, samples=names, gdict="templates")
        histo = histo.group("dataset", hist.Cat("process", "Process", sorting="placement"), process_groups)
    
        for jmult in njets_to_run:
            lepdir = orig_lepdir.replace("NJETS", jmult.lower())

            #set_trace()
            lin_histo = Plotter.linearize_hist(histo[:, :, jmult, lep].integrate("jmult").integrate("leptype"))
            for signal in signals:
                if "Int" in signal:
                    boson, mass, width, pI, wt = tuple(signal.split("_"))
                else:
                    boson, mass, width, pI = tuple(signal.split("_"))
                sub_name = "_".join(["%s%s" % (boson[0], mass[1:]), "relw%s" % widthTOname(width).split("W")[-1], pI.lower(), wt]) if pI == "Int" else "_".join(["%s%s" % (boson[0], mass[1:]), "relw%s" % widthTOname(width).split("W")[-1], pI.lower()])
    
                #set_trace()
                for sys in systs:
                    if sys not in systematics.template_sys_to_name[args.year].keys(): continue
                    if not lin_histo[signal, sys].values().keys():
                        print(f"Systematic {sys} for {lep} {jmult} {signal} not found, skipping")
                        continue
    
                    print(args.year, lep, jmult, sub_name, sys)
                    outhname = "_".join(list(filter(None, [sub_name, systematics.template_sys_to_name[args.year][sys][0], lepdir, (args.year)[-2:]])))
                    if "LEP" in outhname: outhname = outhname.replace("LEP", "muon") if lep == "Muon" else outhname.replace("LEP", "electron")

                    template_histo = lin_histo[signal, sys].integrate("process").integrate("sys")

                        ## save template histos to coffea dict
                    if jmult == "3Jets":
                        histo_dict_3j[lep][f"{signal}_{sys}"] = template_histo.copy()
                    if jmult == "4PJets":
                        histo_dict_4pj[lep][f"{signal}_{sys}"] = template_histo.copy()

                        ## save template histo to root file
                    upfout[outhname] = hist.export1d(template_histo)

    if "3Jets" in njets_to_run:
        coffea_out_3j = os.path.join(outdir, f"test_raw_templates_lj_3Jets_sig_{args.year}_{jobid}.coffea")
        save(histo_dict_3j, coffea_out_3j)
        print(f"{coffea_out_3j} written")
    if "4PJets" in njets_to_run:
        coffea_out_4pj = os.path.join(outdir, f"test_raw_templates_lj_4PJets_sig_{args.year}_{jobid}.coffea")
        save(histo_dict_4pj, coffea_out_4pj)
        print(f"{coffea_out_4pj} written")

    upfout.close()
    print(f"{tmp_rname} written")


def write_correct_template_format(in_fname, isSignal=None):
    """
    Opens temporary root file where template distributions are and then saves them with the correct structure/naming
    """
    if isSignal is None:
        raise ValueError("isSignal needs to be set to True to write signal templates, False for background")
    signame = "sig" if isSignal else "bkg"

    rfile = root_open(in_fname) if in_fname.endswith(".root") else root_open("%s.root" % in_fname)

    if args.naming == "dir":
        out_rfile = os.path.join(outdir, "%s.root" % "_".join(["test", "raw", "templates", "lj", signame, args.year, jobid, "SeparateDirs"]))

        mu_3j_dirs = sorted(set(["_".join((key.name).split("_")[-2:]) for key in rfile.keys() if "mu3jets" in key.name])) if "3Jets" in njets_to_run else None
        el_3j_dirs = sorted(set(["_".join((key.name).split("_")[-2:]) for key in rfile.keys() if "e3jets" in key.name])) if "3Jets" in njets_to_run else None
        mu_4pj_dirs = sorted(set(["_".join((key.name).split("_")[-2:]) for key in rfile.keys() if "mu4pjets" in key.name])) if "4PJets" in njets_to_run else None
        el_4pj_dirs = sorted(set(["_".join((key.name).split("_")[-2:]) for key in rfile.keys() if "e4pjets" in key.name])) if "4PJets" in njets_to_run else None

        #set_trace()
        with root_open(out_rfile, "w") as rout:
                # mu, 3jets
            for dirname in mu_3j_dirs:
                mu_3j_dir = rout.mkdir(dirname)
                mu_3j_dir.cd()
                keys = [key.name for key in rfile.keys() if dirname in key.name]
                for key in keys:
                    hname = key.split("_%s" % dirname)[0]
                    histo = rfile.Get(key)
                    histo.name = hname
                    if (hname == "data_obs") and (args.maskData):
                        histo.Reset()
                    mu_3j_dir.WriteTObject(histo, hname)
        
                # el, 3jets
            for dirname in el_3j_dirs:
                el_3j_dir = rout.mkdir(dirname)
                el_3j_dir.cd()
                keys = [key.name for key in rfile.keys() if dirname in key.name]
                for key in keys:
                    hname = key.split("_%s" % dirname)[0]
                    histo = rfile.Get(key)
                    histo.name = hname
                    if (hname == "data_obs") and (args.maskData):
                        histo.Reset()
                    el_3j_dir.WriteTObject(histo, hname)
        
                # mu, 4pjets
            for dirname in mu_4pj_dirs:
                mu_4pj_dir = rout.mkdir(dirname)
                mu_4pj_dir.cd()
                keys = [key.name for key in rfile.keys() if dirname in key.name]
                for key in keys:
                    hname = key.split("_%s" % dirname)[0]
                    histo = rfile.Get(key)
                    histo.name = hname
                    if (hname == "data_obs") and (args.maskData):
                        histo.Reset()
                    mu_4pj_dir.WriteTObject(histo, hname)
        
                # el, 4pjets
            for dirname in el_4pj_dirs:
                el_4pj_dir = rout.mkdir(dirname)
                el_4pj_dir.cd()
                keys = [key.name for key in rfile.keys() if dirname in key.name]
                for key in keys:
                    hname = key.split("_%s" % dirname)[0]
                    histo = rfile.Get(key)
                    histo.name = hname
                    if (hname == "data_obs") and (args.maskData):
                        histo.Reset()
                    el_4pj_dir.WriteTObject(histo, hname)

    if args.naming == "name":
        out_rfile = os.path.join(outdir, "%s.root" % "_".join(["test", "raw", "templates", "lj", signame, args.year, jobid, "SeparateNames"]))

        mu_keys = [key.name for key in rfile.keys() if (("mu3jets" in key.name) or ("mu4pjets" in key.name))]
        el_keys = [key.name for key in rfile.keys() if (("e3jets" in key.name) or ("e4pjets" in key.name))]

        #set_trace()
        with root_open(out_rfile, "w") as rout:
                # mu
            mu_dir = rout.mkdir("mujets")
            mu_dir.cd()
            for key in mu_keys:
                histo = rfile.Get(key)
                if ("data_obs" in key) and (args.maskData):
                    histo.Reset()
                mu_dir.WriteTObject(histo, key)
        
                # el
            el_dir = rout.mkdir("ejets")
            el_dir.cd()
            for key in el_keys:
                histo = rfile.Get(key)
                if ("data_obs" in key) and (args.maskData):
                    histo.Reset()
                el_dir.WriteTObject(histo, key)
    
    print(f"{out_rfile} written")



if __name__ == "__main__":
    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]
    
    f_ext = "TOT.coffea"
        ## initialize lumi scaling files 
    lumi_corr_dict = load(os.path.join(proj_dir, "Corrections", base_jobid, "MC_LumiWeights.coffea"))

    #set_trace()
        # define variables to get histogram for signal
    #sig_analyzer = "htt_signal_reweight"
    sig_analyzer = "htt_btag_sb_regions"
    sig_input_dir = os.path.join(proj_dir, "results", "%s_%s" % (args.year, jobid), sig_analyzer)
    if os.path.isdir(sig_input_dir):
        sig_fnames = sorted(["%s/%s" % (sig_input_dir, fname) for fname in os.listdir(sig_input_dir) if fname.endswith(f_ext)])
    else: print("No signal file found.")
    
        # define variables to get histogram for background    
    #bkg_analyzer = "htt_btag_iso_cut"
    bkg_analyzer = "htt_btag_sb_regions"
    bkg_input_dir = os.path.join(proj_dir, "results", "%s_%s" % (args.year, jobid), bkg_analyzer)
    if os.path.isdir(bkg_input_dir):
        bkg_fnames = sorted(["%s/%s" % (bkg_input_dir, fname) for fname in os.listdir(bkg_input_dir) if fname.endswith(f_ext)])
    else: print("No background file found.")
    
    outdir = os.path.join(proj_dir, "results", "%s_%s" % (args.year, jobid), "Templates_%s" % bkg_analyzer)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    
    linearize_binning = (
        #np.array([300.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 625.0, 650.0, 675.0, 700.0, 730.0, 760.0, 800.0, 850.0, 900.0, 1000.0, 1200.0]), # orig
        np.array([300.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 625.0, 650.0, 675.0,
            700.0, 730.0, 760.0, 800.0, 850.0, 900.0, 950., 1000.0, 1050., 1100.0, 1150., 1200., 1300., 1500., 2000.]),
        #np.array([300.0, 340.0, 360.0, 380.0, 400.0, 420.0, 440.0, 460.0, 480.0, 500.0, 520.0, 540.0, 560.0, 580.0, 600.0, 620.0, 650., 700.0, 750.0, 800.0, 850.0, 900.0, 2000.0]),
        np.array([0.0, 0.4, 0.6, 0.75, 0.9, 1.0])
    )
    #    # binning for signal 2d dists
    #mtt_ctstar_2d_binning = (
    #    np.arange(300., 2005., 5.),
    #    #np.arange(300., 1205., 5.),
    #    #1,
    #    np.array([0.0, 0.4, 0.6, 0.75, 0.9, 1.0])
    #)

    if not args.only_sig:
        try:
            temp_bkg_rname = "tmp_bkg.root"
            print("Creating background templates")
            get_bkg_templates(temp_bkg_rname)
            write_correct_template_format(temp_bkg_rname, isSignal=False)
            os.system("rm %s" % temp_bkg_rname)
            print(f"{temp_bkg_rname} deleted")
            
        except:
            print("Could not write background templates to file")

    if not args.only_bkg:
        try:
            temp_sig_rname = "tmp_sig.root"
            print("Creating signal templates")
            get_sig_templates(temp_sig_rname)
            write_correct_template_format(temp_sig_rname, isSignal=True)
            os.system("rm %s" % temp_sig_rname)
            print(f"{temp_sig_rname} deleted")
            
        except:
            print("Could not write signal templates to file")

    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
