#! /bin/env python

import time
tic = time.time()

from coffea.util import load
from pdb import set_trace
import os
import numpy as np
import Utilities.systematics as systematics
from coffea import hist
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
parser.add_argument("templates_to_run", type=str, help="Choose which type of templates to run, multiple options can be input as ':' separated strings.")
parser.add_argument("--MEreweight_opts", nargs="*", action=ParseKwargs, help="Options to specify which mass and width points to produce for ME reweighted signal.")
parser.add_argument("--maskData", action="store_false", help="Mask templates for data, default is True.")
args = parser.parse_args()


def final_bkg_templates(years_to_run):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """

    #set_trace()
    if len(years_to_run) == 4:
        outdir = global_outdir
        rname = os.path.join(outdir, f"final_templates_lj_bkg_{jobid}.root")
        upfout = uproot.recreate(rname, compression=uproot.ZLIB(4)) if os.path.isfile(rname) else uproot.create(rname)

    for year in years_to_run:
        if len(years_to_run) < 4:
            outdir = os.path.join(proj_dir, "results", f"{year}_{jobid}", f"Templates_{analyzer}", "FINAL")
            rname = os.path.join(outdir, f"final_templates_lj_bkg_{year}_{jobid}.root")
            upfout = uproot.recreate(rname, compression=uproot.ZLIB(4)) if os.path.isfile(rname) else uproot.create(rname)

        base_bkg_template_name = f"final_templates_lj_bkg_{year}_{jobid}.coffea"
        input_dir = os.path.join(proj_dir, "results", f"{year}_{jobid}", f"Templates_{analyzer}", "FINAL")
        bkg_fname = os.path.join(input_dir, base_bkg_template_name)
        if not os.path.isfile(bkg_fname): raise ValueError("No background file found.")
        hdict = load(bkg_fname)

        #set_trace()
        for jmult in hdict.keys():
            for lep, histo in hdict[jmult].items():
                orig_lepdir = "muNJETS" if lep == "Muon" else "eNJETS"
                lepdir = orig_lepdir.replace("NJETS", jmult.lower())

                if year == "2016APV": year_to_use = "2016pre"
                elif year == "2016": year_to_use = "2016post"
                else: year_to_use = year
                dirname = f"{lepdir}_{year_to_use}"

                upfout.mkdir(dirname)

                #set_trace()
                systypes = ["nosys"] + sorted(systematics.sys_groups[year].keys())
                for sys in systypes:
                    if (jobid == "Summer20UL_POG_lepSFs"):
                        #if ("RECO" in sys) and (lep == "Muon"):
                        #    #set_trace()
                        #    continue
                        if ("syst" in sys) or ("stat" in sys):
                            #set_trace()
                            continue

                    if jobid == "Summer20UL_regroupedJECs":
                        if (sys == "nosys") or (sys == "JES_FlavorQCD") or (sys == "JES_RelativeBal"): continue
                        #if (sys == "nosys") or (sys == "JES_FlavorQCD") or (sys == "JES_RelativeBal") or (sys == f"JES_RelativeSample_{args.year}"): continue

                        # find histograms of associated systematics and their processes
                    #if sys == "EWQCD_SHAPE": set_trace()

                    if sys == "nosys":
                        procs = sorted(set([key.split(f"_{sys}")[0] for key in histo.keys() if sys in key]))
                    else:
                        up_sysname = systematics.sys_groups[year][sys][0]
                        dw_sysname = systematics.sys_groups[year][sys][1]
                        procs = sorted(set([key.split(f"_{up_sysname}")[0] for key in sorted(histo.keys()) if up_sysname in key])) if not dw_sysname \
                            else sorted(set([key.split(f"_{up_sysname}")[0] for key in sorted(histo.keys()) if up_sysname in key] + [key.split(f"_{dw_sysname}")[0] for key in sorted(histo.keys()) if dw_sysname in key]))

                    if not procs: continue

                    for proc in procs:
                        print(year, lep, jmult, sys, proc)
                        #if (proc == "TB") and (lep == "Muon") and (jmult == "4PJets"): set_trace()
                        if sys == "nosys":
                            template, treatment = histo[f"{proc}_{sys}"]
                            if np.any(np.isnan(template.values()[()])):
                                print(f"\tSome bins contain nan values!")
                                set_trace()
                            if proc == "data_obs": template.clear()
                            upfout[dirname][proc] = template.to_hist()

                        else:
                            up_template, treatment = histo[f"{proc}_{up_sysname}"]
                            dw_template, treatment = histo[f"{proc}_{dw_sysname}"]

                            up_outhname = "_".join(list(filter(None, [proc, systematics.combine_template_sys_to_name[year][up_sysname]])))
                            if "LEP" in up_outhname: up_outhname = up_outhname.replace("LEP", lep[0].lower())
                            dw_outhname = "_".join(list(filter(None, [proc, systematics.combine_template_sys_to_name[year][dw_sysname]])))
                            if "LEP" in dw_outhname: dw_outhname = dw_outhname.replace("LEP", lep[0].lower())

                            #set_trace()
                                # replace '2016APV' with '2016pre' and '2016' with '2016post'
                            if (year == "2016APV") and ("2016APV" in up_outhname):
                                #set_trace()
                                up_outhname = up_outhname.replace("2016APV", "2016pre")
                                dw_outhname = dw_outhname.replace("2016APV", "2016pre")
                            if (year == "2016") and ("2016" in up_outhname):
                                up_outhname = up_outhname.replace("2016", "2016post")
                                dw_outhname = dw_outhname.replace("2016", "2016post")

                            if np.any(np.isnan(up_template.values()[()])):
                                print(f"\tSome bins in up_template contain nan values!")
                                set_trace()
                            if np.any(np.isnan(dw_template.values()[()])):
                                print(f"\tSome bins in dw_template contain nan values!")
                                set_trace()

                            upfout[dirname][up_outhname] = up_template.to_hist()
                            upfout[dirname][dw_outhname] = dw_template.to_hist()
                    
    upfout.close()
    print(f"{rname} written")


def final_sig_templates(years_to_run):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """

    if len(years_to_run) == 4:
        outdir = global_outdir
        rname = os.path.join(outdir, f"final_templates_lj_sig_{jobid}.root")
        upfout = uproot.recreate(rname, compression=uproot.ZLIB(4)) if os.path.isfile(rname) else uproot.create(rname)

    #set_trace()
    for year in years_to_run:
        if len(years_to_run) < 4:
            outdir = os.path.join(proj_dir, "results", f"{year}_{jobid}", f"Templates_{analyzer}", "FINAL")
            rname = os.path.join(outdir, f"final_templates_lj_sig_{year}_{jobid}.root")
            upfout = uproot.recreate(rname, compression=uproot.ZLIB(4)) if os.path.isfile(rname) else uproot.create(rname)

        base_sig_template_name = f"final_templates_lj_sig_{year}_{jobid}.coffea"
        input_dir = os.path.join(proj_dir, "results", f"{year}_{jobid}", f"Templates_{analyzer}", "FINAL")
        sig_fname = os.path.join(input_dir, base_sig_template_name)
        if not os.path.isfile(sig_fname): raise ValueError("No signal file found.")
        hdict = load(sig_fname)
        
        for jmult in hdict.keys():
            for lep, histo in hdict[jmult].items():
                orig_lepdir = "muNJETS" if lep == "Muon" else "eNJETS"
                lepdir = orig_lepdir.replace("NJETS", jmult.lower())

                if year == "2016APV": year_to_use = "2016pre"
                elif year == "2016": year_to_use = "2016post"
                else: year_to_use = year
                dirname = f"{lepdir}_{year_to_use}"

                upfout.mkdir(dirname)

                systs = []
                tnames = sorted(histo.keys())
                for tname in tnames:
                    if "Res" in tname:
                        signal = "_".join([tname.split("_Res")[0], "Res"])
                    else:
                        signal = "_".join([tname.split("_neg")[0], "neg"]) if "neg" in tname else "_".join([tname.split("_pos")[0], "pos"])
                    sys = sorted(filter(None, tname.split(f"{signal}_")))[0]
                    systs.append(sys)

                #set_trace()
                systs = sorted(set(systs))
                systs.remove("nosys")
                systypes = ["nosys"]+sorted(filter(None, sorted(set([baseSys(systematics.sys_to_name[year][sys]) for sys in systs]))))
                for sys in systypes:
                    if (jobid == "Summer20UL_POG_lepSFs"):
                        #if ("RECO" in sys) and (lep == "Muon"):
                        #    #set_trace()
                        #    continue
                        if ("syst" in sys) or ("stat" in sys):
                            #set_trace()
                            continue

                    if jobid == "Summer20UL_regroupedJECs":
                        if (sys == "nosys") or (sys == "JES_FlavorQCD") or (sys == "JES_RelativeBal"): continue
                        #if (sys == "nosys") or (sys == "JES_FlavorQCD") or (sys == "JES_RelativeBal") or (sys == f"JES_RelativeSample_{args.year}"): continue

                        # find histograms of associated systematics and their processes
                    if sys == "nosys":
                        signals = sorted(set([key.split(f"_{sys}")[0] for key in histo.keys() if sys in key]))
                    else:
                        up_sysname = [key for key, val in systematics.sys_to_name[year].items() if val == f"{sys}_UP"][0]
                        dw_sysname = [key for key, val in systematics.sys_to_name[year].items() if val == f"{sys}_DW"][0]
                        signals = sorted(set([key.split(f"_{up_sysname}")[0] for key in histo.keys() if up_sysname in key] + [key.split(f"_{dw_sysname}")[0] for key in histo.keys() if dw_sysname in key]))

                    #set_trace()
                    for signal in signals:
                        print(year, lep, jmult, sys, signal)
                        if "Int" in signal:
                            boson, mass, width, pI, wt = tuple(signal.split("_"))
                        else:
                            boson, mass, width, pI = tuple(signal.split("_"))

                        sub_name = "_".join([boson[0], mass.lower(), width.lower(), wt]) if pI == "Int" else "_".join([boson[0], mass.lower(), width.lower(), pI.lower()])

                        if sys == "nosys":
                            template, treatment = histo[f"{signal}_{sys}"]

                            if np.any(np.isnan(template.values()[()])):
                                print(f"\tSome bins contain nan values!")
                                set_trace()
                            upfout[dirname][sub_name] = template.to_hist()

                        else:
                            #set_trace()
                            up_template, treatment = histo[f"{signal}_{up_sysname}"]
                            dw_template, treatment = histo[f"{signal}_{dw_sysname}"]

                            up_outhname = "_".join(list(filter(None, [sub_name, systematics.combine_template_sys_to_name[year][up_sysname]])))
                            if "LEP" in up_outhname: up_outhname = up_outhname.replace("LEP", lep[0].lower())
                            dw_outhname = "_".join(list(filter(None, [sub_name, systematics.combine_template_sys_to_name[year][dw_sysname]])))
                            if "LEP" in dw_outhname: dw_outhname = dw_outhname.replace("LEP", lep[0].lower())

                                # replace '2016APV' with '2016pre' and '2016' with '2016post' 
                            if (year == "2016APV") and ("2016APV" in up_outhname):
                                up_outhname = up_outhname.replace("2016APV", "2016pre")
                                dw_outhname = dw_outhname.replace("2016APV", "2016pre")
                            if (year == "2016") and ("2016" in up_outhname):
                                up_outhname = up_outhname.replace("2016", "2016post")
                                dw_outhname = dw_outhname.replace("2016", "2016post")

                            if np.any(np.isnan(up_template.values()[()])):
                                print(f"\tSome bins in up_template contain nan values!")
                                set_trace()
                            if np.any(np.isnan(dw_template.values()[()])):
                                print(f"\tSome bins in dw_template contain nan values!")
                                set_trace()

                            upfout[dirname][up_outhname] = up_template.to_hist()
                            upfout[dirname][dw_outhname] = dw_template.to_hist()

    upfout.close()
    print(f"{rname} written")


def final_MEsig_templates(years_to_run):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """

    possible_masses = [365, 380, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000]
    possible_widths = [2.5, 3.0, 4.0, 5.0, 8.0, 10.0, 13.0, 15.0, 18.0, 21.0, 25.0]
    #possible_widths = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 8.0, 10.0, 13.0, 15.0, 18.0, 21.0, 25.0]
    possible_masses = [str(mass) for mass in possible_masses]
    possible_widths = [str(width) for width in possible_widths]

    #set_trace()
    masses_to_run = args.MEreweight_opts.get("allowed_masses", "All")
    widths_to_run = args.MEreweight_opts.get("allowed_widths", "All")

    if masses_to_run == "All":
        allowed_masses = possible_masses
    else:
        allowed_masses = masses_to_run.split(":")
        allowed_masses = [mass for mass in allowed_masses if mass in possible_masses]

    if widths_to_run == "All":
        allowed_widths = possible_widths
    else:
        allowed_widths = widths_to_run.split(":")
        allowed_widths = [width for width in allowed_widths if width in possible_widths]

    #set_trace()
    allowed_masses = ["M"+mass for mass in allowed_masses]
    allowed_widths = ["W"+width.replace(".", "p") for width in allowed_widths]
    print(f"Making target signal points for {allowed_masses} and {allowed_widths}")



    if len(years_to_run) == 4:
        outdir = global_outdir
        if masses_to_run == "All":
            rname = os.path.join(outdir, f"final_templates_lj_MEsig_{jobid}.root" if (widths_to_run == "All") else f"final_templates_lj_MEsig_{''.join(allowed_widths).lower()}_{jobid}.root")
        elif widths_to_run == "All":
            rname = os.path.join(outdir, f"final_templates_lj_MEsig_{jobid}.root" if (masses_to_run == "All") else f"final_templates_lj_MEsig_{''.join(allowed_masses).lower()}_{jobid}.root")
        else:
            rname = os.path.join(outdir, f"final_templates_lj_MEsig_{''.join(allowed_masses).lower()}_{''.join(allowed_widths).lower()}_{jobid}.root")

        upfout = uproot.recreate(rname, compression=uproot.ZLIB(4)) if os.path.isfile(rname) else uproot.create(rname)

    #set_trace()
    for year in years_to_run:
        if len(years_to_run) < 4:
            outdir = os.path.join(proj_dir, "results", f"{year}_{jobid}", f"Templates_{analyzer}", "FINAL")
            if masses_to_run == "All":
                rname = os.path.join(outdir, f"final_templates_lj_MEsig_{year}_{jobid}.root" if (widths_to_run == "All") else f"final_templates_lj_MEsig_{''.join(allowed_widths).lower()}_{year}_{jobid}.root")
            elif widths_to_run == "All":
                rname = os.path.join(outdir, f"final_templates_lj_MEsig_{year}_{jobid}.root" if (masses_to_run == "All") else f"final_templates_lj_MEsig_{''.join(allowed_masses).lower()}_{year}_{jobid}.root")
            else:
                rname = os.path.join(outdir, f"final_templates_lj_MEsig_{year}_{''.join(allowed_masses).lower()}_{''.join(allowed_widths).lower()}_{jobid}.root")

            upfout = uproot.recreate(rname, compression=uproot.ZLIB(4)) if os.path.isfile(rname) else uproot.create(rname)

        base_sig_template_name = f"final_templates_lj_MEsig_{year}_{jobid}.coffea"
        input_dir = os.path.join(proj_dir, "results", f"{year}_{jobid}", f"Templates_{analyzer}", "FINAL")
        sig_fname = os.path.join(input_dir, base_sig_template_name)
        if not os.path.isfile(sig_fname): raise ValueError("No signal file found.")
        hdict = load(sig_fname)
        
        #set_trace()
        for jmult in hdict.keys():
            for lep, histo in hdict[jmult].items():
                orig_lepdir = "muNJETS" if lep == "Muon" else "eNJETS"
                lepdir = orig_lepdir.replace("NJETS", jmult.lower())

                if year == "2016APV": year_to_use = "2016pre"
                elif year == "2016": year_to_use = "2016post"
                else: year_to_use = year
                dirname = f"{lepdir}_{year_to_use}"

                upfout.mkdir(dirname)
                systs = []
                tnames = sorted(histo.keys())
                for tname in tnames:
                    if "Res" in tname:
                        signal = "_".join([tname.split("_Res")[0], "Res"])
                    else:
                        signal = "_".join([tname.split("_neg")[0], "neg"]) if "neg" in tname else "_".join([tname.split("_pos")[0], "pos"])
                    sys = sorted(filter(None, tname.split(f"{signal}_")))[0]
                    systs.append(sys)

                #set_trace()
                systs = sorted(set(systs))
                systs.remove("nosys")
                systypes = ["nosys"]+sorted(filter(None, sorted(set([baseSys(systematics.sys_to_name[year][sys]) for sys in systs]))))
                for sys in systypes:
                    if jobid == "Summer20UL_regroupedJECs":
                        if (sys == "nosys") or (sys == "JES_FlavorQCD") or (sys == "JES_RelativeBal"): continue
                        #if (sys == "nosys") or (sys == "JES_FlavorQCD") or (sys == "JES_RelativeBal") or (sys == f"JES_RelativeSample_{year}"): continue

                    if (jobid == "Summer20UL_POG_lepSFs"):
                        #if ("RECO" in sys) and (lep == "Muon"):
                        #    #set_trace()
                        #    continue
                        if ("syst" in sys) or ("stat" in sys):
                            #set_trace()
                            continue

                        # find histograms of associated systematics and their processes
                    if sys == "nosys":
                        signals = sorted(set([key.split(f"_{sys}")[0] for key in histo.keys() if sys in key]))
                    else:
                        up_sysname = [key for key, val in systematics.sys_to_name[year].items() if val == f"{sys}_UP"][0]
                        dw_sysname = [key for key, val in systematics.sys_to_name[year].items() if val == f"{sys}_DW"][0]
                        signals = sorted(set([key.split(f"_{up_sysname}")[0] for key in histo.keys() if up_sysname in key] + [key.split(f"_{dw_sysname}")[0] for key in histo.keys() if dw_sysname in key]))

                    #set_trace()
                    for signal in signals:
                        ### make sure mass and width point for this signal point is allowed
                        mass, width = signal.split("_")[1], signal.split("_")[2]
                        if not ((mass in allowed_masses) and (width in allowed_widths)): continue
                        ###

                        print(year, lep, jmult, sys, signal)
                        if "Int" in signal:
                            boson, mass, width, pI, wt = tuple(signal.split("_"))
                        else:
                            boson, mass, width, pI = tuple(signal.split("_"))

                        sub_name = "_".join([boson[0], mass.lower(), width.lower(), wt]) if pI == "Int" else "_".join([boson[0], mass.lower(), width.lower(), pI.lower()])

                        if sys == "nosys":
                            template, treatment = histo[f"{signal}_{sys}"]

                            if np.any(np.isnan(template.values()[()])):
                                print(f"\tSome bins contain nan values!")
                                set_trace()
                            upfout[dirname][sub_name] = template.to_hist()

                        else:
                            #set_trace()
                            up_template, treatment = histo[f"{signal}_{up_sysname}"]
                            dw_template, treatment = histo[f"{signal}_{dw_sysname}"]

                            up_outhname = "_".join(list(filter(None, [sub_name, systematics.combine_template_sys_to_name[year][up_sysname]])))
                            if "LEP" in up_outhname: up_outhname = up_outhname.replace("LEP", lep[0].lower())
                            dw_outhname = "_".join(list(filter(None, [sub_name, systematics.combine_template_sys_to_name[year][dw_sysname]])))
                            if "LEP" in dw_outhname: dw_outhname = dw_outhname.replace("LEP", lep[0].lower())

                                # replace '2016APV' with '2016pre' and '2016' with '2016post' for files sent to Afiq
                            if (year == "2016APV") and ("2016APV" in up_outhname):
                                #set_trace()
                                up_outhname = up_outhname.replace("2016APV", "2016pre")
                                dw_outhname = dw_outhname.replace("2016APV", "2016pre")
                            if (year == "2016") and ("2016" in up_outhname):
                                up_outhname = up_outhname.replace("2016", "2016post")
                                dw_outhname = dw_outhname.replace("2016", "2016post")

                            if np.any(np.isnan(up_template.values()[()])):
                                print(f"\tSome bins in up_template contain nan values!")
                                set_trace()
                            if np.any(np.isnan(dw_template.values()[()])):
                                print(f"\tSome bins in dw_template contain nan values!")
                                set_trace()

                            upfout[dirname][up_outhname] = up_template.to_hist()
                            upfout[dirname][dw_outhname] = dw_template.to_hist()


    upfout.close()
    print(f"{rname} written")


def final_pdf_templates(years_to_run):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """

    #set_trace()
    if len(years_to_run) == 4:
        outdir = global_outdir
        rname = os.path.join(outdir, f"final_pdf_templates_lj_bkg_{jobid}.root")
        upfout = uproot.recreate(rname, compression=uproot.ZLIB(4)) if os.path.isfile(rname) else uproot.create(rname)

    for year in years_to_run:
        if len(years_to_run) < 4:
            outdir = os.path.join(proj_dir, "results", f"{year}_{jobid}", f"Templates_{analyzer}", "FINAL")
            rname = os.path.join(outdir, f"final_pdf_templates_lj_bkg_{year}_{jobid}.root")
            upfout = uproot.recreate(rname, compression=uproot.ZLIB(4)) if os.path.isfile(rname) else uproot.create(rname)

        base_bkg_template_name = f"final_pdf_templates_lj_bkg_{year}_{jobid}.coffea"
        input_dir = os.path.join(proj_dir, "results", f"{year}_{jobid}", f"Templates_{analyzer}", "FINAL")
        bkg_fname = os.path.join(input_dir, base_bkg_template_name)
        if not os.path.isfile(bkg_fname): raise ValueError("No background file found.")
        hdict = load(bkg_fname)

        for jmult in hdict.keys():
            for lep, histo in hdict[jmult].items():
                orig_lepdir = "muNJETS" if lep == "Muon" else "eNJETS"
                lepdir = orig_lepdir.replace("NJETS", jmult.lower())

                if year == "2016APV": year_to_use = "2016pre"
                elif year == "2016": year_to_use = "2016post"
                else: year_to_use = year
                dirname = f"{lepdir}_{year_to_use}"
                upfout.mkdir(dirname)

                #set_trace()
                systs = sorted(set(["_".join(key.split("_")[1:]) for key in histo.keys() if not ("data_obs" in key or len(key.split("_")) == 1 or "shape" in key)]))
                for sys in systs:
                    if sys == "nosys": continue
                    #set_trace()
                        # find histograms of associated systematics and their processes
                    procs = sorted(set([key.split(f"_{sys}")[0] for key in histo.keys() if sys in key]))
                    for proc in procs:
                        print(year, lep, jmult, sys, proc)
                        if sys == "nosys":
                            template, treatment = histo[f"{proc}_{sys}"]
                            upfout[dirname][proc] = template.to_hist()

                        else:
                            template, treatment = histo[f"{proc}_{sys}"]
                            outhname = "_".join([proc, "CMS", "PDF", sys.split("_")[-1]])
                            upfout[dirname][outhname] = template.to_hist()

    upfout.close()
    print(f"{rname} written")


if __name__ == "__main__":
    allowed_template_options = ["bkg", "sig", "MEreweight_sig", "PDF"]
    templates_to_run = (args.templates_to_run).split(":")
    templates_to_run = [template for template in templates_to_run if template in allowed_template_options]

    allowed_year_options = ["2016APV", "2016", "2017", "2018"]
    years_to_run = (args.years_to_run).split(":")
    years_to_run = [year for year in years_to_run if year in allowed_year_options]

    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]

    global_outdir = os.path.join(proj_dir, "results", jobid, f"Templates_FINAL")
    if not os.path.isdir(global_outdir):
        os.makedirs(global_outdir)

    baseSys = lambda sys : "_".join(sys.split("_")[:-1])

    if "bkg" in templates_to_run:
        analyzer = "htt_btag_sb_regions"

        print("Creating final background templates")
        final_bkg_templates(years_to_run)

    if "sig" in templates_to_run:
        analyzer = "htt_btag_sb_regions"

        print("Creating final signal templates")
        final_sig_templates(years_to_run)

    if "MEreweight_sig" in templates_to_run:
        analyzer = "htt_btag_sb_regions"

        print("Creating final ME reweighted signal templates")
        final_MEsig_templates(years_to_run)

    if "PDF" in templates_to_run:
        analyzer = "htt_pdfUncs"

        print("Creating final background templates")
        final_pdf_templates(years_to_run)


    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
