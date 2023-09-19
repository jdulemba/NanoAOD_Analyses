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

import argparse
class ParseKwargs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, dict())
        for value in values:
            key, value = value.split("=")
            getattr(namespace, self.dest)[key] = value


from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("templates_to_run", type=str, help="Choose which type of templates to run, multiple options can be input as ':' separated strings.")
parser.add_argument("--nomSMTTxsec", action="store_true", help="Apply nominal SM cross sections to top mass and LHE scale weights")
args = parser.parse_args()


def get_bkg_templates():
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """
    ## variables that only need to be defined/evaluated once
    templates_to_check = systematics.combine_template_sys_to_name
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
        process_groups = plt_tools.make_dataset_groups(lep, args.year, samples=names, bkgdict="templates")
        
        lumi_correction = lumi_corr_dict[args.year][f"{lep}s"]
                # scale ttJets events, split by reconstruction type, by normal ttJets lumi correction
        if len(ttJets_cats) > 0:
            for tt_cat in ttJets_cats:
                ttJets_lumi_topo = "_".join(tt_cat.split("_")[:-2]) if "sl_tau" in tt_cat else "_".join(tt_cat.split("_")[:-1]) # gets ttJets[SL, Had, DiLep] or ttJets_PS
                ttJets_eff_lumi = lumi_correction[ttJets_lumi_topo]
                lumi_correction.update({tt_cat: ttJets_eff_lumi})
    
        histo = rebin_histo.copy()
        #set_trace()
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

        systs = sorted(set([key[1] for key in histo.values().keys()]))
        systs.insert(0, systs.pop(systs.index("nosys"))) # move "nosys" to the front
        #set_trace()

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
                #if ("CR1" in sys) or ("CR2" in sys) or ("erdON" in sys): set_trace()
                #if sys not in templates_to_check[args.year].keys(): continue

                    # EWQCD background estimation only needed for 'nosys'
                if "data_obs" not in sorted(set([key[0] for key in sig_histo.values().keys()])):
                    sys_histo = sig_histo[:, sys].integrate("sys")
                else:
                    sys_histo = Plotter.BKG_Est(sig_reg=sig_histo[:, sys].integrate("sys"), sb_reg=cen_sb_histo, norm_type="Sideband", sys=sys, uncs_percentage=None, isForTemplates=True) if sys == "nosys" else sig_histo[:, sys].integrate("sys")

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
                        shapeUP_name, shapeDW_name = systematics.sys_groups[args.year]["SHAPE"]
                            # Up region
                        bkg_shapeUp = Plotter.data_minus_top(up_sb_histo, isForTemplates=True) # find data - (tt+st) distributions
                        bkg_shapeUp.scale(np.sum(template_histo.values()[()])/np.sum(bkg_shapeUp.values()[()])) # normalize to yield in signal region
                        print(args.year, lep, jmult, shapeUP_name, proc)
                            # Down region
                        bkg_shapeDown = Plotter.data_minus_top(dw_sb_histo, isForTemplates=True) # find data - (tt+st) distributions
                        bkg_shapeDown.scale(np.sum(template_histo.values()[()])/np.sum(bkg_shapeDown.values()[()])) # normalize to yield in signal region
                        print(args.year, lep, jmult, shapeDW_name, proc)

                            ## save template histos to coffea dict
                        histo_dict[jmult][lep][f"{proc}_{sys}"] = template_histo.copy()
                        histo_dict[jmult][lep][f"{proc}_{shapeUP_name}"] = bkg_shapeUp.copy()
                        histo_dict[jmult][lep][f"{proc}_{shapeDW_name}"] = bkg_shapeDown.copy()

                        #set_trace()
                        ## estimate unc related to scaling top+tt contribution in central control region
                        TTsubUP_name, TTsubDW_name = systematics.sys_groups[args.year]["EWQCD_TTsub"]

                        cen_ratio = cen_sb_histo[Plotter.data_samples].integrate("process").values()[()]/cen_sb_histo[Plotter.mc_samples].integrate("process").values()[()]
                        cen_ratio_err = np.sqrt(cen_sb_histo[Plotter.data_samples].integrate("process").values(sumw2=True)[()][1]/np.square(cen_sb_histo[Plotter.data_samples].integrate("process").values()[()]) + \
                            cen_sb_histo[Plotter.mc_samples].integrate("process").values(sumw2=True)[()][1]/np.square(cen_sb_histo[Plotter.mc_samples].integrate("process").values()[()]))
                        #set_trace()
                        cen_sf = np.polyfit(np.arange(cen_ratio[np.isfinite(cen_ratio_err) & np.isfinite(cen_ratio)].size), cen_ratio[np.isfinite(cen_ratio_err) & np.isfinite(cen_ratio)], deg=0, w=np.reciprocal(cen_ratio_err[np.isfinite(cen_ratio_err) & np.isfinite(cen_ratio)]))[0]

                            # find original dmt shape
                        orig_top_mc_hist  = cen_sb_histo[Plotter.template_top_samples].integrate("process").copy()
                        orig_top_mc_hist.scale(-1)
                        orig_data_minus_top = cen_sb_histo[Plotter.data_samples].integrate("process").copy()
                        orig_data_minus_top.add(orig_top_mc_hist)
                        orig_cen_dmt_sumw, orig_cen_dmt_sumw2 = Plotter.get_qcd_shape(orig_data_minus_top)

                            # find shape of scaled dmt dist
                        scaled_top_mc_hist  = cen_sb_histo[Plotter.template_top_samples].integrate("process").copy()
                        scaled_top_mc_hist.scale(-1*cen_sf)
                        scaled_data_minus_top = cen_sb_histo[Plotter.data_samples].integrate("process").copy()
                        scaled_data_minus_top.add(scaled_top_mc_hist)
                        scaled_cen_dmt_sumw, scaled_cen_dmt_sumw2 = Plotter.get_qcd_shape(scaled_data_minus_top)

                            # symmetrize relative differences
                        TTsubUP_rel_sumw = scaled_cen_dmt_sumw/orig_cen_dmt_sumw
                        TTsubUP_rel_sumw2 = scaled_cen_dmt_sumw2/np.square(scaled_cen_dmt_sumw) + orig_cen_dmt_sumw2/np.square(orig_cen_dmt_sumw)
                        TTsubDW_rel_sumw = 2 - TTsubUP_rel_sumw

                        # create histos
                            # up
                        TTsubUP_histo = template_histo.copy() 
                        TTsubUP_histo.values(overflow="all", sumw2=True)[()][0][:] = TTsubUP_rel_sumw * template_histo.values(overflow="all")[()]
                        TTsubUP_histo.values(overflow="all", sumw2=True)[()][1][:] = TTsubUP_rel_sumw2 * np.square(TTsubUP_rel_sumw * template_histo.values(overflow="all")[()])

                            # down
                        TTsubDW_histo = template_histo.copy() 
                        TTsubDW_histo.values(overflow="all", sumw2=True)[()][0][:] = TTsubDW_rel_sumw * template_histo.values(overflow="all")[()]
                        TTsubDW_histo.values(overflow="all", sumw2=True)[()][1][:] = TTsubUP_rel_sumw2 * np.square(TTsubDW_rel_sumw * template_histo.values(overflow="all")[()])
            
                        histo_dict[jmult][lep][f"{proc}_{TTsubUP_name}"] = TTsubUP_histo.copy()
                        histo_dict[jmult][lep][f"{proc}_{TTsubDW_name}"] = TTsubDW_histo.copy()

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

                            ## save template histos to coffea dict
                        histo_dict[jmult][lep][f"{proc}_{sys}"] = template_histo.copy()

    #set_trace()
    coffea_out = os.path.join(outdir, f"raw_templates_lj_bkg_nomSMTTxsec_{args.year}_{jobid}.coffea" if args.nomSMTTxsec else f"raw_templates_lj_bkg_{args.year}_{jobid}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")



def get_pdf_templates():
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """
    ## variables that only need to be defined/evaluated once
    hdict = plt_tools.add_coffea_files(pdf_fnames) if len(pdf_fnames) > 1 else load(pdf_fnames[0])
    histo_dict = processor.dict_accumulator({njets : {"Muon" : {}, "Electron" :{}} for njets in njets_to_run})

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

    for lep in ["Muon", "Electron"]:
        #set_trace()
        ## make groups based on process
        process_groups = plt_tools.make_dataset_groups(lep, args.year, samples=names, bkgdict="templates")

        lumi_correction = lumi_corr_dict[args.year][f"{lep}s"]
                # scale ttJets events, split by reconstruction type, by normal ttJets lumi correction
        if len(ttJets_cats) > 0:
            for tt_cat in ttJets_cats:
                ttJets_lumi_topo = "_".join(tt_cat.split("_")[:-2]) if "sl_tau" in tt_cat else "_".join(tt_cat.split("_")[:-1]) # gets ttJets[SL, Had, DiLep] or ttJets_PS
                ttJets_eff_lumi = lumi_correction[ttJets_lumi_topo]
                lumi_correction.update({tt_cat: ttJets_eff_lumi})

        histo = rebin_histo.copy()
        histo.scale(lumi_correction, axis="dataset")
        histo = histo.group(process_cat, process, process_groups)[:, :, :, lep].integrate("leptype")

            # remove 0th and last 2 pdf replicas 0th is base set compatible with 1, last two sets are variations in alpha_S
        histo = histo.remove(["pdf_0"], "sys")
        systs = sorted(set([key[1] for key in histo.values().keys()]))
        systs.insert(0, systs.pop(systs.index("nosys"))) # move "nosys" to the front

            # loop over each jet multiplicity
        for jmult in njets_to_run:
            sig_histo = Plotter.linearize_hist(histo[:, :, jmult].integrate("jmult"))

                # loop over each systematic
            for sys in systs:
                sys_histo = sig_histo[:, sys].integrate("sys")

                    ## write nominal and systematic variations for each topology to file
                for proc in sorted(set([key[0] for key in sys_histo.values().keys()])):
                    if not sys_histo[proc].values().keys():
                        print(f"Systematic {sys} for {lep} {jmult} {proc} not found, skipping")
                        continue

                    print(args.year, lep, jmult, sys, proc)
                    template_histo = sys_histo[proc].integrate("process")

                        ## save template histos to coffea dict
                    histo_dict[jmult][lep][f"{proc}_{sys}"] = template_histo.copy()

    coffea_out = os.path.join(pdf_outdir, f"raw_pdf_templates_lj_bkg_{args.year}_{jobid}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")


def make_output_dir(analyzer):
    outdir = os.path.join(eos_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}")
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    return outdir
  

if __name__ == "__main__":
    allowed_template_options = ["bkg", "PDF"]
    templates_to_run = [template for template in (args.templates_to_run).split(":") if template in allowed_template_options]
    templates_to_not_run = [template for template in (args.templates_to_run).split(":") if template not in allowed_template_options]
    if templates_to_not_run:
        print(f"{templates_to_not_run} are not valid options for making templates, will be skipped")

    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]
    base_jobid = os.environ["base_jobid"]
    eos_dir = os.environ["eos_dir"]

    njets_to_run = sorted(["3Jets", "4PJets"])

        ## initialize lumi scaling files
    lumi_name = "MC_LumiWeights"
    if args.nomSMTTxsec:
        lumi_name += "_nomSMTTxsec"
    lumi_corr_dict = load(os.path.join(proj_dir, "Corrections", base_jobid, f"{lumi_name}.coffea"))
  
    tt_LHEscale_wts_name_dict = {
        "FACTORDown" : "uF_down",
        "FACTORUp"   : "uF_up",
        "RENORMDown" : "uR_down",
        "RENORMUp"   : "uR_up",
    }

    jet_pars = prettyjson.loads(open(os.path.join(proj_dir, "cfg_files", f"cfg_pars_{jobid}.json")).read())["Jets"]
    btagger = jet_pars["btagger"]
    wps_to_use = list(set([jet_pars["permutations"]["tightb"],jet_pars["permutations"]["looseb"]]))
    if not( len(wps_to_use) == 1):
        raise IOError("Only 1 unique btag working point supported now")
    btag_wp = btagger+wps_to_use[0]
    btag_reg_names_dict = btag_sidebands.btag_reg_names_dict[args.year][btag_wp]

    linearize_binning = (
        final_binning.mtt_binning,
        final_binning.ctstar_abs_binning
    ) 

    #set_trace()
    if "bkg" in templates_to_run:
        analyzer = "htt_btag_sb_regions"
            # define variables to get histogram for background    
        bkg_input_dir = os.path.join(eos_dir, "results", f"{args.year}_{jobid}", analyzer)
        if os.path.isdir(bkg_input_dir):
            bkg_fnames = fnmatch.filter(os.listdir(bkg_input_dir), "*BKG*TOT.coffea")
            bkg_fnames = [os.path.join(bkg_input_dir, fname) for fname in bkg_fnames]
        else: raise ValueError("No background file found.")
    
        outdir = make_output_dir(analyzer)

        print("Creating background templates")
        get_bkg_templates()

    if "PDF" in templates_to_run:
            # define variables to get histogram for background
        analyzer = "htt_pdfUncs"
        pdf_outdir = make_output_dir(analyzer)

        pdf_input_dir = os.path.join(eos_dir, "results", f"{args.year}_{jobid}", analyzer)
        if os.path.isdir(pdf_input_dir):
            pdf_fnames = fnmatch.filter(os.listdir(pdf_input_dir), "*TOT.coffea")
            pdf_fnames = [os.path.join(pdf_input_dir, fname) for fname in pdf_fnames]
        else: print("No PDF file found.")

        print("Creating PDF templates")
        get_pdf_templates()


    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
