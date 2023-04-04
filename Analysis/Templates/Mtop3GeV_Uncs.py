#! /bin/env python

import time
tic = time.time()

from coffea.util import load, save
from pdb import set_trace
import os
import numpy as np
import fnmatch
import Utilities.Plotter as Plotter
import coffea.processor as processor    
import Utilities.final_analysis_binning as final_binning

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("temp_type", choices=["raw", "comb_era_lepton", "comb_lepton"], help="What type of combination do you want to be performed.")
args = parser.parse_args()

def make_raw_templates(fnames_dict):
    """
    Function that writes original linearized mtt vs costheta distributions to coffea file.
    """
    ## variables that only need to be defined/evaluated once
    hname_to_use = "mtt_vs_tlep_ctstar_abs"
    ttJets_permcats = ["*right", "*matchable", "*unmatchable", "*sl_tau", "*other"]

        ## make groups based on process
    process = hist.Cat("process", "Process", sorting="placement")
    process_cat = "dataset"
            ## hardcoded!!!!
    process_groups = {
        "TT_mtop1695" : ["ttJetsDiLep_mtop1695_other", "ttJetsHad_mtop1695_other", "ttJetsSL_mtop1695_matchable", "ttJetsSL_mtop1695_right", "ttJetsSL_mtop1695_sl_tau", "ttJetsSL_mtop1695_unmatchable"],
        "TT_mtop1755" : ["ttJetsDiLep_mtop1755_other", "ttJetsHad_mtop1755_other", "ttJetsSL_mtop1755_matchable", "ttJetsSL_mtop1755_right", "ttJetsSL_mtop1755_sl_tau", "ttJetsSL_mtop1755_unmatchable"],
    }

    histo_dict = processor.dict_accumulator({njets : {"Muon" : {}, "Electron" :{}} for njets in njets_to_run})

    for year, fname in fnames_dict.items():
        hdict = load(fname)

            # get correct hist and rebin
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
        names = [dataset for dataset in sorted(set([key[0] for key in histo.values().keys()]))] # get dataset names in hists
        ttJets_cats = [name for name in names if any([fnmatch.fnmatch(name, cat) for cat in ttJets_permcats])] # gets ttJets(_PS)_other, ...

        for lep in ["Muon", "Electron"]:
            lumi_correction = lumi_corr_dict[year][f"{lep}s"]
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

            #set_trace()
            if sorted(set([key[4] for key in histo.values().keys()])) != ["btagPass"]:
                raise ValueError("'btagPass' is not the only region for the 'btag' axis")
            histo = histo.group(process_cat, process, process_groups)[:, :, :, lep, :].integrate("leptype").integrate("btag")

            if sorted(set([key[0] for key in histo.values().keys()])) != ["TT_mtop1695", "TT_mtop1755"]:
                raise ValueError("'TT' is not the only process")

            systs = sorted(set([key[1] for key in histo.values().keys()]))
            systs.insert(0, systs.pop(systs.index("nosys"))) # move "nosys" to the front
            #set_trace()

                # loop over each jet multiplicity
            for jmult in njets_to_run:
                sig_histo = Plotter.linearize_hist(histo[:, :, jmult].integrate("jmult"))

                for (proc, sys) in sig_histo.values().keys():
                    print(year, lep, jmult, sys, proc)
                    #set_trace()
                    template_histo = sig_histo[(proc, sys)].integrate("process").integrate("sys").copy()
                    sumw, sumw2 = template_histo.values(sumw2=True)[()]
                    rel_err = np.sqrt(sumw2)/np.abs(sumw)
                    rel_err_mask = rel_err > 10
                    if np.any(rel_err_mask):
                        set_trace()
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

                        ## save template histos to coffea dict
                    histo_dict[jmult][lep][f"{proc}_{sys}"] = template_histo.copy()

        #set_trace()
        coffea_out = os.path.join(output_dir, f"raw_templates_{year}.coffea")
        save(histo_dict, coffea_out)
        print(f"{coffea_out} written")


def make_combined_year_and_lepton_templates(fnames_dict):
    """
    Function that writes linearized mtt vs costheta distributions that have been combined across eras and lepton channels to coffea file.
    """
    #set_trace()
    histo_dict = processor.dict_accumulator({njets : {} for njets in njets_to_run})
    for jmult in njets_to_run:
            # find common names across lepton channels
        mu_tnames_dict = {year : sorted(fnames_dict[year][jmult]["Muon"].keys()) for year in fnames_dict.keys() }
        el_tnames_dict = {year : sorted(fnames_dict[year][jmult]["Electron"].keys()) for year in fnames_dict.keys() }
        tnames = sorted( set.intersection(*map(set, sorted(mu_tnames_dict.values()))) & set.intersection(*map(set, sorted(el_tnames_dict.values()))) )
        for tname in tnames:
            print(jmult, tname)

            combined_year_template = None
            idx = 0
            for year in fnames_dict.keys():
                for lep in sorted(["Muon", "Electron"]):
                    combined_year_template = fnames_dict[year][jmult][lep][tname].copy() if idx == 0 else combined_year_template.add(fnames_dict[year][jmult][lep][tname].copy())
                    idx += 1

                ## save template histos to coffea dict
            histo_dict[jmult][tname] = combined_year_template.copy()

    #set_trace()
    coffea_out = os.path.join(input_dir, "raw_combined_year_and_lepton_templates.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")

def make_combined_lep_templates(fnames_dict):
    """
    Function that writes linearized mtt vs costheta distributions that have been combined across lepton channels for each era to coffea file.
    """
    #set_trace()
    for year, hdict in fnames_dict.items():
        histo_dict = processor.dict_accumulator({njets : {} for njets in njets_to_run})
        for jmult in hdict.keys():
                # find common names across lepton channels
            mu_tnames = sorted(hdict[jmult]["Muon"].keys())
            el_tnames = sorted(hdict[jmult]["Electron"].keys())
            tnames = sorted(set(mu_tnames) & set(el_tnames))
            for tname in tnames:
                print(year, jmult, tname)

                combined_lep_template = hdict[jmult]["Muon"][tname].copy()
                combined_lep_template = combined_lep_template.add(hdict[jmult]["Electron"][tname].copy())

                    ## save template histos to coffea dict
                histo_dict[jmult][tname] = combined_lep_template.copy()

        #set_trace()
        coffea_out = os.path.join(input_dir, f"raw_combined_lep_templates_{year}.coffea")
        save(histo_dict, coffea_out)
        print(f"{coffea_out} written")




if __name__ == "__main__":
    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]
    base_jobid = os.environ["base_jobid"]
    eos_dir = os.environ["eos_dir"]
    analyzer = "Mtop3GeV_Uncs"

    njets_to_run = sorted(["3Jets", "4PJets"])

        ## initialize lumi scaling files
    lumi_corr_dict = load(os.path.join(proj_dir, "Corrections", base_jobid, "MC_LumiWeights.coffea"))
  
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

    if args.temp_type == "raw":
        from coffea import hist
        input_dir = os.path.join(eos_dir, "results", jobid, analyzer)
        output_dir = os.path.join(eos_dir, "results", jobid, analyzer, "Templates")
        fnames_dict = {year : os.path.join(input_dir, year, fnmatch.filter(os.listdir(os.path.join(input_dir, year)), "*TOT.coffea")[0]) for year in ["2016APV", "2016", "2017", "2018"]}
        #set_trace()
        print("Creating raw templates")
        make_raw_templates(fnames_dict)

    if args.temp_type == "comb_era_lepton":
        input_dir = os.path.join(eos_dir, "results", jobid, analyzer, "Templates")
        fnames_dict = {year : load(os.path.join(input_dir, f"raw_templates_{year}.coffea")) for year in ["2016APV", "2016", "2017", "2018"]}
        print("Combining templates across era+lepton channels")
        make_combined_year_and_lepton_templates(fnames_dict)

    if args.temp_type == "comb_lepton":
        input_dir = os.path.join(eos_dir, "results", jobid, analyzer, "Templates")
        fnames_dict = {year : load(os.path.join(input_dir, f"raw_templates_{year}.coffea")) for year in ["2016APV", "2016", "2017", "2018"]}
        print("Combining templates across lepton channels")
        make_combined_lep_templates(fnames_dict)


    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
