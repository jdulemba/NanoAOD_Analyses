#! /bin/env python

import time
tic = time.time()

from coffea.util import load, save
from pdb import set_trace
import os
import coffea.processor as processor    
    
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("templates_to_run", type=str, help="Choose which type of templates to run, multiple options can be input as ':' separated strings.")
parser.add_argument("channels_to_combine", type=str, help="Choose how to combine systematic templates (by lepton or eras), multiple options can be input as ':' separated strings.")
parser.add_argument("--nomSMTTxsec", action="store_true", help="Apply nominal SM cross sections to top mass and LHE scale weights")
args = parser.parse_args()


def combine_lep_templates(fname, process):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """
    #set_trace()
    histo_dict = processor.dict_accumulator({njets : {} for njets in njets_to_run})
    hdict = load(fname)
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

    coffea_out = os.path.join(input_dir, f"raw_combined_lep_templates_lj_{process}_{year}_{jobid}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")


def combine_year_and_lepton_templates(fnames_dict, process):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """
    #set_trace()
    histo_dict = processor.dict_accumulator({njets : {} for njets in njets_to_run})
    years_hdict = {key : load(fnames_dict[key]) for key in fnames_dict.keys()}
    for jmult in njets_to_run:
            # find common names across lepton channels
        mu_tnames_dict = {year : sorted(years_hdict[year][jmult]["Muon"].keys()) for year in years_hdict.keys() }
        el_tnames_dict = {year : sorted(years_hdict[year][jmult]["Electron"].keys()) for year in years_hdict.keys() }
        tnames = sorted( set.intersection(*map(set, sorted(mu_tnames_dict.values()))) & set.intersection(*map(set, sorted(el_tnames_dict.values()))) )
        for tname in tnames:
            print(jmult, tname)

            combined_year_template = None
            idx = 0
            for year in years_hdict.keys():
                for lep in sorted(["Muon", "Electron"]):
                    combined_year_template = years_hdict[year][jmult][lep][tname].copy() if idx == 0 else combined_year_template.add(years_hdict[year][jmult][lep][tname].copy())
                    idx += 1

                ## save template histos to coffea dict
            histo_dict[jmult][tname] = combined_year_template.copy()

    coffea_out = os.path.join(output_dir, f"raw_combined_year_and_lepton_templates_lj_{process}_nomSMTTxsec_{jobid}.coffea" if args.nomSMTTxsec else f"raw_combined_year_and_lepton_templates_lj_{process}_{jobid}.coffea")
    save(histo_dict, coffea_out)
    print(f"{coffea_out} written")



if __name__ == "__main__":
    allowed_template_options = ["bkg", "PDF"]
    templates_to_run = (args.templates_to_run).split(":")
    templates_to_run = [template for template in templates_to_run if template in allowed_template_options]

    allowed_combination_options = ["lepton", "era_lepton"]
    combinations_to_run = (args.channels_to_combine).split(":")
    combinations_to_run = [combination for combination in combinations_to_run if combination in allowed_combination_options]

    jobid = os.environ["jobid"]
    eos_dir = os.environ["eos_dir"]

    njets_to_run = ["3Jets", "4PJets"]
    analyzer = "htt_btag_sb_regions"

    for combination in combinations_to_run:
        if combination == "lepton":
            for year in ["2016APV", "2016", "2017", "2018"]:
                input_dir = os.path.join(eos_dir, "results", f"{year}_{jobid}", f"Templates_{analyzer}")
                if not os.path.isdir(input_dir): raise ValueError("No background file found.")

                if "bkg" in templates_to_run:
                    base_bkg_template_name = f"raw_templates_lj_bkg_{year}_{jobid}.coffea"
                    bkg_fname = os.path.join(input_dir, base_bkg_template_name)
                    if not os.path.isfile(bkg_fname): raise ValueError(f"{bkg_fname} not found.")
                    print(f"Combining e+mu channels in {year} for background templates")
                    combine_lep_templates(bkg_fname, "bkg")

        if combination == "era_lepton":
            years_to_combine = ["2016APV", "2016", "2017", "2018"]

            if "PDF" in templates_to_run:
                analyzer = "htt_pdfUncs"
                output_dir = os.path.join(eos_dir, "results", jobid, f"Templates_{analyzer}")
                if not os.path.isdir(output_dir): os.makedirs(output_dir)

                fnames_dict = {year : os.path.join(eos_dir, "results", f"{year}_{jobid}", f"Templates_{analyzer}", f"raw_pdf_templates_lj_bkg_{year}_{jobid}.coffea") for year in years_to_combine}
                for year in years_to_combine:
                    if not os.path.isfile(fnames_dict[year]): raise ValueError(f"{fnames_dict[year]} not found.")
                print(f"Combining channels across years and leptons for PDF templates")
                combine_year_and_lepton_templates(fnames_dict, "PDF")

            else:
                output_dir = os.path.join(eos_dir, "results", jobid, f"Templates_{analyzer}")
                if not os.path.isdir(output_dir): os.makedirs(output_dir)

                if "bkg" in templates_to_run:
                    bkg_fnames_dict = {year : os.path.join(eos_dir, "results", f"{year}_{jobid}", f"Templates_{analyzer}", f"raw_templates_lj_bkg_nomSMTTxsec_{year}_{jobid}.coffea" if args.nomSMTTxsec \
                        else f"raw_templates_lj_bkg_{year}_{jobid}.coffea") for year in years_to_combine}
                    for year in years_to_combine:
                        if not os.path.isfile(bkg_fnames_dict[year]): raise ValueError(f"{bkg_fnames_dict[year]} not found.")
                    print(f"Combining channels across years and leptons for background templates")
                    combine_year_and_lepton_templates(bkg_fnames_dict, "bkg")


    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
