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
import uproot

import argparse
class ParseKwargs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, dict())
        for value in values:
            key, value = value.split("=")
            getattr(namespace, self.dest)[key] = value


from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("years_to_run", type=str, help="Choose which year to run, multiple options can be input as ':' separated strings.")
parser.add_argument("--no_est", action="store_true", help="Do not perform background estimation, default is True.")
args = parser.parse_args()

version="V29"

def get_bkg_templates():
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """
    #set_trace()
    if len(years_to_run) == 4:
        outdir = make_output_dir(analyzer, year="Run2")
        rname = os.path.join(outdir, f"control_plot_templates_run2_{jobid}_noEst_{version}.root" if args.no_est else f"control_plot_templates_run2_{jobid}_{version}.root")
        upfout = uproot.recreate(rname, compression=uproot.ZLIB(4)) if os.path.isfile(rname) else uproot.create(rname)


    ## variables that only need to be defined/evaluated once
    templates_to_check = systematics.combine_template_sys_to_name

    for year in years_to_run:
        if len(years_to_run) < 4:
            outdir = make_output_dir(analyzer, year=year)
            #rname = os.path.join(outdir, f"test.root")
            rname = os.path.join(outdir, f"control_plot_templates_{year}_{jobid}_noEst_{version}.root" if args.no_est else f"control_plot_templates_{year}_{jobid}_{version}.root")
            upfout = uproot.recreate(rname, compression=uproot.ZLIB(4)) if os.path.isfile(rname) else uproot.create(rname)

        bkg_input_dir = os.path.join(eos_dir, "results", f"{year}_{jobid}", analyzer)
        if os.path.isdir(bkg_input_dir):
            bkg_fnames = fnmatch.filter(os.listdir(bkg_input_dir), "*BKG*TOT.coffea")
            bkg_fnames = [os.path.join(bkg_input_dir, fname) for fname in bkg_fnames]
        else: raise ValueError("No background file found.")
        
        hdict = plt_tools.add_coffea_files(bkg_fnames) if len(bkg_fnames) > 1 else load(bkg_fnames[0])

        btag_reg_names_dict = btag_sidebands.btag_reg_names_dict[year][btag_wp]

        if year == "2016APV": year_to_use = "2016pre"
        elif year == "2016": year_to_use = "2016post"
        else: year_to_use = year

        #set_trace()
            # get correct hist and rebin
        for hname_to_use, rebinning in variables.items():
            if hname_to_use not in hdict.keys():
                print(f"{hname_to_use} not found in file")
                continue

            print(hname_to_use)
            orig_histo = hdict[hname_to_use][Plotter.nonsignal_samples].copy() # process, sys, jmult, leptype, btag, lepcat
            #set_trace()

            if orig_histo.dense_dim() == 1:
                xaxis_name = orig_histo.dense_axes()[0].name
                rebin_histo = orig_histo.rebin(xaxis_name, rebinning)
            else:        
                xrebinning, yrebinning = rebinning
                xaxis_name = orig_histo.dense_axes()[0].name
                yaxis_name = orig_histo.dense_axes()[1].name
                    ## rebin x axis
                if isinstance(xrebinning, np.ndarray):
                    new_xbins = hist.Bin(xaxis_name, xaxis_name, xrebinning)
                elif isinstance(xrebinning, float) or isinstance(xrebinning, int):
                    new_xbins = xrebinning
                    ## rebin y axis
                if isinstance(yrebinning, np.ndarray):
                    new_ybins = hist.Bin(yaxis_name, yaxis_name, yrebinning)
                elif isinstance(yrebinning, float) or isinstance(yrebinning, int):
                    new_ybins = yrebinning
                rebin_histo = orig_histo.rebin(xaxis_name, new_xbins).rebin(yaxis_name, new_ybins)

            #set_trace()    
                ## scale ttJets events, split by reconstruction type, by normal ttJets lumi correction
            ttJets_permcats = ["*right", "*matchable", "*unmatchable", "*sl_tau", "*other"]
            names = [dataset for dataset in sorted(set([key[0] for key in rebin_histo.values().keys()]))] # get dataset names in hists
            ttJets_cats = [name for name in names if any([fnmatch.fnmatch(name, cat) for cat in ttJets_permcats])] # gets ttJets(_PS)_other, ...

                ## make groups based on process
            process = hist.Cat("process", "Process", sorting="placement")
            process_cat = "dataset"

            for lep in ["Muon", "Electron"]:
                #set_trace()    
                ## make groups based on process
                process_groups = plt_tools.make_dataset_groups(lep, year, samples=names, bkgdict="control_plots" if args.no_est else "templates")
                
                lumi_correction = lumi_corr_dict[year][f"{lep}s"]
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

                systs = sorted(set([key[1] for key in histo.values().keys()]))
                systs.insert(0, systs.pop(systs.index("nosys"))) # move "nosys" to the front
                #set_trace()

                    # loop over each jet multiplicity
                for jmult in njets_to_run:
                    orig_lepdir = "muNJETS" if lep == "Muon" else "eNJETS"
                    lepdir = orig_lepdir.replace("NJETS", jmult.lower())

                    dirname = f"{lepdir}_{year_to_use}"

                    upfout.mkdir(dirname+"/"+hname_to_use)
                    tmp_histo_dict = {}

                        # get sideband and signal region hists
                    if histo.dense_dim() == 1:
                        cen_sb_histo = histo[:, "nosys", jmult, btag_reg_names_dict["Central"]["reg"]].integrate("jmult").integrate("btag").integrate("sys")
                        up_sb_histo = histo[:, "nosys", jmult, btag_reg_names_dict["Up"]["reg"]].integrate("jmult").integrate("btag").integrate("sys")
                        dw_sb_histo = histo[:, "nosys", jmult, btag_reg_names_dict["Down"]["reg"]].integrate("jmult").integrate("btag").integrate("sys")
                        sig_histo = histo[:, :, jmult, btag_reg_names_dict["Signal"]["reg"]].integrate("jmult").integrate("btag")
                    else:
                        cen_sb_histo = Plotter.linearize_hist(histo[:, "nosys", jmult, btag_reg_names_dict["Central"]["reg"]].integrate("jmult").integrate("btag").integrate("sys"))
                        up_sb_histo = Plotter.linearize_hist(histo[:, "nosys", jmult, btag_reg_names_dict["Up"]["reg"]].integrate("jmult").integrate("btag").integrate("sys"))
                        dw_sb_histo = Plotter.linearize_hist(histo[:, "nosys", jmult, btag_reg_names_dict["Down"]["reg"]].integrate("jmult").integrate("btag").integrate("sys"))
                        sig_histo = Plotter.linearize_hist(histo[:, :, jmult, btag_reg_names_dict["Signal"]["reg"]].integrate("jmult").integrate("btag"))

                        # loop over each systematic
                    for sys in systs:
                        if (jobid == "Summer20UL_POG_lepSFs") or (jobid == "Summer20UL_DeepJet"):
                            #if ("RECO" in sys) and (lep == "Muon"):
                            #    #set_trace()
                            #    continue
                            if (lep == "Muon") and (("IDtot" in sys) or ("ISOtot" in sys) or ("TRIGtot" in sys)):
                                continue
                        #set_trace()
                        if sys not in templates_to_check[year].keys():
                            if not ((sys == "CR1") or (sys == "CR2") or (sys == "erdON")): continue
                            

                            # EWQCD background estimation only needed for 'nosys'
                        if "data_obs" not in sorted(set([key[0] for key in sig_histo.values().keys()])):
                            sys_histo = sig_histo[:, sys].integrate("sys")
                        else:
                            #set_trace()
                            if args.no_est:
                                sys_histo = sig_histo[:, sys].integrate("sys")
                            else:
                                sys_histo = Plotter.BKG_Est(sig_reg=sig_histo[:, sys].integrate("sys"), sb_reg=cen_sb_histo, norm_type="Sideband", sys=sys, isForTemplates=True) if sys == "nosys" else sig_histo[:, sys].integrate("sys")

                            ## write nominal and systematic variations for each topology to file
                        for proc in sorted(set([key[0] for key in sys_histo.values().keys()])):
                            if (("TT" and "ttJets") not in proc) and (sys in systematics.ttJets_sys.values()): continue
                            if (proc == "data_obs") and not (sys == "nosys"): continue
                            if not sys_histo[proc].values().keys():
                                print(f"Systematic {sys} for {lep} {jmult} {proc} not found, skipping")
                                continue

                            if (sys == "CR1") or (sys == "CR2") or (sys == "erdON"):
                                #set_trace()
                                outhname = "_".join(list(filter(None, [proc, systematics.combine_template_sys_to_name[year][sys+"Up"]])))
                            #set_trace()
                            else:
                                outhname = "_".join(list(filter(None, [proc, systematics.combine_template_sys_to_name[year][sys]])))
                            if "LEP" in outhname: outhname = outhname.replace("LEP", lep[0].lower())
                            if "_tot" in outhname: outhname = outhname.replace("_tot", "")
                            if sys == "SHAPE":
                                #set_trace()
                                outhname = outhname.replace("CHAN", lepdir)
                            if (year == "2016APV") and ("2016APV" in outhname):
                                #set_trace()
                                outhname = outhname.replace("2016APV", "2016pre")
                            if (year == "2016") and ("2016" in outhname):
                                #set_trace()
                                outhname = outhname.replace("2016", "2016post")

                            if (proc == "EWQCD"):
                                if sys != "nosys": continue
                                print(year, hname_to_use, lep, jmult, sys, proc)

                                    # get original nominal distribution
                                template_histo = sys_histo[proc].integrate("process")

                                # get shape variations from btag sb regions
                                shapeUP_name, shapeDW_name = systematics.sys_groups[year]["SHAPE"]
                                    # Up region
                                bkg_shapeUp = Plotter.data_minus_top(up_sb_histo, isForTemplates=True) # find data - (tt+st) distributions
                                bkg_shapeUp.scale(np.sum(template_histo.values()[()])/np.sum(bkg_shapeUp.values()[()])) # normalize to mc yield in signal region
                                print(year, hname_to_use, lep, jmult, shapeUP_name, proc)
                                    # Down region
                                bkg_shapeDown = Plotter.data_minus_top(dw_sb_histo, isForTemplates=True) # find data - (tt+st) distributions
                                bkg_shapeDown.scale(np.sum(template_histo.values()[()])/np.sum(bkg_shapeDown.values()[()])) # normalize to mc yield in signal region
                                print(year, hname_to_use, lep, jmult, shapeDW_name, proc)

                                    ## save template histos to coffea dict
                                tmp_histo_dict[f"{proc}_{sys}"] = template_histo.copy()
                                upfout[dirname][hname_to_use][outhname] = (template_histo.copy()).to_hist()
                                upfout[dirname][hname_to_use][f"{proc}_{shapeUP_name}"] = (bkg_shapeUp.copy()).to_hist()
                                upfout[dirname][hname_to_use][f"{proc}_{shapeDW_name}"] = (bkg_shapeDown.copy()).to_hist()
                            else:
                                print(year, hname_to_use, lep, jmult, sys, proc)
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
                                        nosys_template = tmp_histo_dict[f"{proc}_nosys"].copy()
                                        nosys_sumw, nosys_sumw2 = nosys_template.values(sumw2=True)[()]
                                        nosys_rel_err_mask = (nosys_sumw == 0.) & (nosys_sumw2 == 0.)
                                        if np.any(rel_err_mask & nosys_rel_err_mask):
                                            print(f"\tRelative error > 10 for this process! Setting bin {np.where(rel_err_mask & nosys_rel_err_mask)[0]} to 0")
                                            sumw[rel_err_mask & nosys_rel_err_mask], sumw2[rel_err_mask & nosys_rel_err_mask] = 0., 0.

                                    #if proc != "TB": set_trace()
                                if (sys == "CR1") or (sys == "CR2") or (sys == "erdON"):
                                    #set_trace()
                                    dw_histo = template_histo.copy()
                                        # get nominal vals
                                    orig_sumw = sig_histo.values()[(proc, "nosys")]
                                        # get ratios
                                    ratio_vals = sumw/orig_sumw
                                        # set 0 values to nan so the yields aren't doubled in the end
                                    ratio_vals[ratio_vals == 0.] = np.nan
                                    #if np.any(np.isnan(ratio_vals)): set_trace()
                                    #ratio_vals = np.nan_to_num(sumw/orig_sumw)
                                        # find relative deviation from nominal vals
                                    up_rel_vals = ratio_vals - 1.
                                    dw_rel_vals = -1. * up_rel_vals

                                    symmetrized_ratios_dw = np.nan_to_num(dw_rel_vals + 1.)
                                    dw_histo.values()[()][:] = symmetrized_ratios_dw*orig_sumw
                                    upfout[dirname][hname_to_use][outhname] = (template_histo.copy()).to_hist()
                                    upfout[dirname][hname_to_use][outhname.replace("Up", "Down")] = (dw_histo.copy()).to_hist()

                                else:
                                        ## save template histos to coffea dict
                                    tmp_histo_dict[f"{proc}_{sys}"] = template_histo.copy()
                                    upfout[dirname][hname_to_use][outhname] = (template_histo.copy()).to_hist()

    upfout.close()
    print(f"{rname} written")


def make_output_dir(analyzer, year):
    outdir = os.path.join(eos_dir, "results", jobid, f"Control_Plot_Templates_{analyzer}", year)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    return outdir
  

if __name__ == "__main__":
    allowed_year_options_ = ["2016APV", "2016", "2017", "2018"]
    years_to_run = (args.years_to_run).split(":")
    years_to_run = [year for year in years_to_run if year in allowed_year_options_]

    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]
    base_jobid = os.environ["base_jobid"]
    eos_dir = os.environ["eos_dir"]

    njets_to_run = sorted(["3Jets", "4PJets"])

        ## initialize lumi scaling files 
    lumi_corr_dict = load(os.path.join(proj_dir, "Corrections", base_jobid, "MC_LumiWeights.coffea"))
  
    tt_LHEscale_wts_name_dict = {
        "FACTORDown" : "uF_down",
        "FACTORUp"   : "uF_up",
        "RENORMDown" : "uR_down",
        "RENORMUp"   : "uR_up",
    }

    #set_trace()
    jet_pars = prettyjson.loads(open(os.path.join(proj_dir, "cfg_files", f"cfg_pars_{jobid}.json")).read())["Jets"]
    btagger = jet_pars["btagger"]
    wps_to_use = list(set([jet_pars["permutations"]["tightb"],jet_pars["permutations"]["looseb"]]))
    if not( len(wps_to_use) == 1):
        raise IOError("Only 1 unique btag working point supported now")
    btag_wp = btagger+wps_to_use[0]

    linearize_binning = (
        final_binning.mtt_binning,
        final_binning.ctstar_abs_binning
    ) 

    phi_eta_binning = (
        np.array([-3.2, -2.4, -1.57, -0.87, 0., 0.87, 1.57, 2.4, 3.2]), # phi
        np.array([-2.5, -1.3, -0.7, 0., 0.7, 1.3, 2.5]) # eta binning
    )

    variables = {
        "mtt_vs_tlep_ctstar_abs" : (linearize_binning[0], linearize_binning[1]),
        "Jets_phi_vs_eta" : (phi_eta_binning[0], phi_eta_binning[1]),
        "Lep_phi_vs_eta" : (phi_eta_binning[0], phi_eta_binning[1]),
        "Jets_njets" : 1,
        "mtt" : 4,
        "tlep_ctstar_abs" : 1,
        "mthad" : 2,
        "mWHad" : 2,
        "mWLep" : 2,
        "pt_thad" : 2,
        "pt_tlep" : 2,
        "pt_tt" : 2,
        "eta_thad" : 2,
        "eta_tlep" : 2,
        "eta_tt" : 2,
        "tlep_ctstar" : 2,
        "full_disc" : 2,
        "mass_disc" : 2,
        "ns_disc" : 2,
        "ns_dist" : 1,
        "Jets_pt" : 1,
        "Jets_eta" : 2,
        "Jets_phi" : 2,
        "Jets_LeadJet_pt" : 1,
        "Jets_LeadJet_eta" : 2,
        #"Jets_DeepCSV_bDisc" : 1,
        #"Jets_DeepJet_bDisc" : 1,
        "Lep_pt" : 1,
        "Lep_eta" : 2,
        "Lep_phi" : 2,
        "Lep_iso" : 1,
        "MT" : 1,
        "MET_pt" : 1,
        "MET_phi" : 1,
    }

    analyzer = "htt_btag_sb_regions"
        # define variables to get histogram for background    
    print("Creating background templates")
    get_bkg_templates()


    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
