#! /bin/env python

import time
tic = time.time()

from coffea.util import load, save
from pdb import set_trace
import os
import Utilities.plot_tools as plt_tools
import Utilities.prettyjson as prettyjson
from coffea import hist
import numpy as np
import fnmatch
import Utilities.Plotter as Plotter
import uproot3
from rootpy.io import root_open
import coffea.processor as processor    
import Utilities.systematics as systematics
import Utilities.Plotter as Plotter

base_jobid = os.environ["base_jobid"]

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016", "2017", "2018"] if base_jobid == "NanoAODv6" else ["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("--njets", default="all", nargs="?", choices=["3", "4+", "all"], help="Specify which jet multiplicity to use.")
parser.add_argument("--only_bkg", action="store_true", help="Make background templates only.")
parser.add_argument("--only_sig", action="store_true", help="Make signal templates only.")
parser.add_argument("--maskData", action="store_false", help="Mask templates for data, default is True.")
args = parser.parse_args()

njets_to_run = []
if (args.njets == "3") or (args.njets == "all"):
    njets_to_run += ["3Jets"]
if (args.njets == "4+") or (args.njets == "all"):
    njets_to_run += ["4PJets"]



def scale_mtop3gev(nominal, template):
    ratio_vals = (template.values()[()] - nominal.values()[()])/nominal.values()[()]
    scaled_vals = ratio_vals*(1./6.)
    interped_vals = (scaled_vals+1)*nominal.values()[()]
    interped_histo = nominal.copy()
    for idx in range(len(interped_vals)):
        interped_histo.values()[()][idx]  = interped_vals[idx]

    return interped_histo



def substitute_ttJets(sys_histo, ttJets_histo, ttJets_PS_histo):
    sys_vals = sys_histo.values()[()]
    ttJets_vals = ttJets_histo.values()[()]
    ttJetsPS_vals = ttJets_PS_histo.values()[()]

        # find (sys-ttJets)/ttJets relative deviation and then scale by ttJets_PS
    scaled_vals = ((sys_vals-ttJets_vals)/ttJets_vals + 1)*ttJetsPS_vals
    scaled_histo = sys_histo.copy()
    for idx in range(len(scaled_vals)):
        scaled_histo.values()[()][idx]  = scaled_vals[idx]

    #set_trace()
    return scaled_histo



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
    histo = hdict[hname_to_use] # process, sys, jmult, leptype, btag, lepcat
    
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
    
    nbins = (len(xrebinning)-1)*(len(yrebinning)-1)
    
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
        lepdir = "mujets" if lep == "Muon" else "ejets"
    
        ## make groups based on process
        process_groups = plt_tools.make_dataset_groups(lep, args.year, samples=names, gdict="templates")
        
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
        histo = histo[Plotter.nonsignal_samples]
        systs = sorted(set([key[1] for key in histo.values().keys()]))
        systs.insert(0, systs.pop(systs.index("nosys"))) # move "nosys" to the front


            # loop over each jet multiplicity
        for jmult in njets_to_run:
            #iso_sb    = Plotter.linearize_hist(histo[:, "nosys", jmult, "btagPass", "Loose"].integrate("sys").integrate("jmult").integrate("lepcat").integrate("btag"))
            #btag_sb   = Plotter.linearize_hist(histo[:, "nosys", jmult, "btagFail", "Tight"].integrate("sys").integrate("jmult").integrate("lepcat").integrate("btag"))
            #double_sb = Plotter.linearize_hist(histo[:, "nosys", jmult, "btagFail", "Loose"].integrate("sys").integrate("jmult").integrate("lepcat").integrate("btag"))
            #sig_histo = Plotter.linearize_hist(histo[:, :, jmult, "btagPass", "Tight"].integrate("jmult").integrate("lepcat").integrate("btag"))
            sig_histo = Plotter.linearize_hist(histo[:, :, jmult, "btagPass"].integrate("jmult").integrate("btag"))

            #set_trace()
                # loop over each systematic
            for sys in systs:
                #qcd_est_histo = Plotter.QCD_Est(sig_reg=sig_histo, iso_sb=iso_sb, btag_sb=btag_sb, double_sb=double_sb, norm_type="Sideband", shape_region="BTAG", norm_region="BTAG", sys=sys)

                #set_trace()
                if sys not in systematics.template_sys_to_name[args.year].keys(): continue

                    ## write nominal and systematic variations for each topology to file
                #for proc in sorted(set([key[0] for key in qcd_est_histo.values().keys()])):
                for proc in sorted(set([key[0] for key in sig_histo.values().keys()])):
                    if (proc != "TT") and (sys in systematics.ttJets_sys.values()): continue
                    if (proc == "data_obs") and not (sys == "nosys"): continue
                    name = proc+lepdir if proc == "QCD" else proc
                    print(lep, jmult, sys, name)
                    #set_trace()
                    outhname = "_".join([jmult, lepdir, name, sys]) if sys == "nosys" else "_".join([jmult, lepdir, name, systematics.template_sys_to_name[args.year][sys][0]])
                    if "LEP" in outhname: outhname = outhname.replace("LEP", "m") if lep == "Muon" else outhname.replace("LEP", "e")

                    #template_histo = qcd_est_histo[proc].integrate("process")
                    if not sig_histo[proc, sys].values().keys():
                        print(f"Systematic {sys} for {lep} {jmult} {proc} not found, skipping")
                        continue
                    template_histo = sig_histo[proc, sys].integrate("process").integrate("sys")

                    if ((sys == "mtop1695") or (sys == "mtop1755")):
                        template_histo = scale_mtop3gev(nominal=histo_dict_3j[lep]["%s_nosys" % proc] if jmult == "3Jets" else histo_dict_4pj[lep]["%s_nosys" % proc], template=template_histo)
                        #set_trace()

                        ## save template histos to coffea dict
                    if jmult == "3Jets":
                        histo_dict_3j[lep]["%s_%s" % (proc, sys)] = template_histo
                        #histo_dict_3j[lep][proc if sys == "nosys" else "%s_%s" % (proc, sys)] = template_histo
                    if jmult == "4PJets":
                        histo_dict_4pj[lep]["%s_%s" % (proc, sys)] = template_histo
                        #histo_dict_4pj[lep][proc if sys == "nosys" else "%s_%s" % (proc, sys)] = template_histo

                        ## save template histo to root file
                    upfout[outhname] = hist.export1d(template_histo)

    if "3Jets" in njets_to_run:
        coffea_out_3j = os.path.join(outdir, "raw_templates_lj_3Jets_bkg_%s_%s.coffea" % (args.year, jobid))
        save(histo_dict_3j, coffea_out_3j)
        print(f"{coffea_out_3j} written")
    if "4PJets" in njets_to_run:
        coffea_out_4pj = os.path.join(outdir, "raw_templates_lj_4PJets_bkg_%s_%s.coffea" % (args.year, jobid))
        save(histo_dict_4pj, coffea_out_4pj)
        print(f"{coffea_out_4pj} written")

    upfout.close()
    print(f"{tmp_rname} written")


def get_sig_templates(tmp_rname):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """
    from rootpy.plotting import Hist2D

    #set_trace()
    widthTOname = lambda width : str(width).replace(".", "p")
    nameTOwidth = lambda width : str(width).replace("p", ".")

    ## variables that only need to be defined/evaluated once
    hdict = plt_tools.add_coffea_files(sig_fnames) if len(sig_fnames) > 1 else load(sig_fnames[0])

        # get correct hist and rebin
    hname_to_use = "mtt_vs_tlep_ctstar_abs"
    if hname_to_use not in hdict.keys():
        raise ValueError(f"{hname_to_use} not found in file")
        #raise ValueError("%s not found in file" % hname_to_use)
    xrebinning, yrebinning = mtt_ctstar_2d_binning
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
    rebin_histo = histo[:, :, :, :, "btagPass"].integrate("btag")

    rebin_histo = rebin_histo[Plotter.signal_samples]
    signals = sorted(set([key[0] for key in rebin_histo.values().keys()]))    
        # only use SL decays
    signals = [sig for sig in signals if "TTJetsSL" in sig]

    systs = sorted(set([key[1] for key in rebin_histo.values().keys()]))
    systs.insert(0, systs.pop(systs.index("nosys"))) # move "nosys" to the front

        # create 2D signal hists and write to temp file        
    with root_open(tmp_rname, "w") as out:
        for lep in ["Muon", "Electron"]:
            lepdir = "mujets" if lep == "Muon" else "ejets"

                # scale by lumi
            lumi_correction = lumi_corr_dict[args.year]["%ss" % lep]
            scaled_histo = rebin_histo.copy()
            scaled_histo.scale(lumi_correction, axis="dataset")
    
            for jmult in njets_to_run:
                histo = scaled_histo[:, :, jmult, lep].integrate("jmult").integrate("leptype")
    
                for signal in signals:
                    if "Int" in signal:
                        boson, mass, width, pI, wt = tuple(signal.split("_"))
                    else:
                        boson, mass, width, pI = tuple(signal.split("_"))
                        wt = "pos"
                    samtype = "int" if pI == "Int" else "sgn"
                    bostype = "ggA" if "AtoTT" in boson else "ggH"
    
                    sub_name = "%s_%s-%s-%s-%s" % (bostype, wt, samtype, widthTOname(width).split("W")[-1]+"pc", mass) if pI == "Int" else "%s_pos-%s-%s-%s" % (bostype, samtype, widthTOname(width).split("W")[-1]+"pc", mass)
    
                    #set_trace()
                    for sys in systs:
                        if sys not in systematics.template_sys_to_name[args.year].keys(): continue
                        template_histo = histo[signal, sys].integrate("dataset").integrate("sys")
                        if wt == "neg":
                            template_histo.scale(-1.)
                        print(lep, jmult, sub_name, sys)
                        sumw, sumw2 = template_histo.values(sumw2=True, overflow="all")[()] # get vals and errors for all bins (including under/overflow)

                            ## create rootpy hist and rename
                        rtpy_h2d = Hist2D(template_histo.dense_axes()[0].edges(), template_histo.dense_axes()[1].edges())
                        outhname = "_".join([jmult, lepdir, sub_name]) if sys == "nosys" else "_".join([jmult, lepdir, sub_name, systematics.template_sys_to_name[args.year][sys][0]])
                        if "LEP" in outhname: outhname = outhname.replace("LEP", "m") if lep == "Muon" else outhname.replace("LEP", "e")
                        rtpy_h2d.name = outhname
                            # set bin content for rootpy hist
                        for binx in range(0, rtpy_h2d.GetNbinsX()+2):
                            for biny in range(0, rtpy_h2d.GetNbinsY()+2):
                                rtpy_h2d[binx, biny] = sumw[binx, biny], sumw2[binx, biny]
                        #set_trace()
                        rtpy_h2d.Write()
        
    print(f"{tmp_rname} written")


def write_correct_template_format(in_fname, isSignal=None):
    """
    Opens temporary root file where template distributions are and then saves them with the correct structure/naming
    """
    if isSignal is None:
        raise ValueError("isSignal needs to be set to True to write signal templates, False for background")
    signame = "sig" if isSignal else "bkg"

    rfile = root_open(in_fname) if in_fname.endswith(".root") else root_open("%s.root" % in_fname)

    #set_trace()
    if "3Jets" in njets_to_run:    
        mu_3j_keys = [key.name for key in rfile.keys() if "3Jets_mujets" in key.name]
        el_3j_keys = [key.name for key in rfile.keys() if "3Jets_ejets" in key.name]

        fname_3j = os.path.join(outdir, "%s.root" % "_".join(["raw", "templates", "lj", "3Jets", signame, args.year, jobid]))
        with root_open(fname_3j, "w") as rout:
            mu_dir = rout.mkdir("mujets")
            mu_dir.cd()

            for key in mu_3j_keys:
                hname = key.split("3Jets_mujets_")[-1]
                if "nosys" in hname: hname = "_".join(hname.split("_")[:-1])
                histo = rfile.Get(key)
                histo.name = hname
                if (hname == "data_obs") and (args.maskData):
                    histo.Reset()
                mu_dir.WriteTObject(histo, hname)
        
            el_dir = rout.mkdir("ejets")
            el_dir.cd()
        
            for key in el_3j_keys:
                hname = key.split("3Jets_ejets_")[-1]
                if "nosys" in hname: hname = "_".join(hname.split("_")[:-1])
                histo = rfile.Get(key)
                histo.name = hname
                if (hname == "data_obs") and (args.maskData):
                    histo.Reset()
                el_dir.WriteTObject(histo, hname)

        print(f"{fname_3j} written")

    if "4PJets" in njets_to_run:
        mu_4pj_keys = [key.name for key in rfile.keys() if "4PJets_mujets" in key.name]
        el_4pj_keys = [key.name for key in rfile.keys() if "4PJets_ejets" in key.name]
        
        fname_4pj = os.path.join(outdir, "%s.root" % "_".join(["raw", "templates", "lj", "4PJets", signame, args.year, jobid]))
        with root_open(fname_4pj, "w") as rout:
            mu_dir = rout.mkdir("mujets")
            mu_dir.cd()
            for key in mu_4pj_keys:
                hname = key.split("4PJets_mujets_")[-1]
                if "nosys" in hname: hname = "_".join(hname.split("_")[:-1])
                histo = rfile.Get(key)
                histo.name = hname
                if (hname == "data_obs") and (args.maskData):
                    histo.Reset()
                mu_dir.WriteTObject(histo, hname)
        
            el_dir = rout.mkdir("ejets")
            el_dir.cd()
            for key in el_4pj_keys:
                hname = key.split("4PJets_ejets_")[-1]
                if "nosys" in hname: hname = "_".join(hname.split("_")[:-1])
                histo = rfile.Get(key)
                histo.name = hname
                if (hname == "data_obs") and (args.maskData):
                    histo.Reset()
                el_dir.WriteTObject(histo, hname)

        print(f"{fname_4pj} written")



if __name__ == "__main__":
    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]
    
    f_ext = "TOT.coffea"
        ## initialize lumi scaling files 
    lumi_corr_dict = load(os.path.join(proj_dir, "Corrections", base_jobid, "MC_LumiWeights.coffea"))

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
        # binning for signal 2d dists
    mtt_ctstar_2d_binning = (
        np.arange(300., 2005., 5.),
        #np.arange(300., 1205., 5.),
        #1,
        np.array([0.0, 0.4, 0.6, 0.75, 0.9, 1.0])
    )

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
