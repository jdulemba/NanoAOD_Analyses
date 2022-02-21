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
    
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("year", choices=["2016APV", "2016", "2017", "2018"], help="What year is the ntuple from.")
parser.add_argument("--maskData", action="store_false", help="Mask templates for data, default is True.")
parser.add_argument("--only_bkg", action="store_true", help="Make background templates only.")
parser.add_argument("--only_sig", action="store_true", help="Make signal templates only.")
parser.add_argument("--kfactors", action="store_true", help="Apply signal k-factors to signal")
parser.add_argument("--scale_mtop3gev", action="store_true", help="Scale 3GeV mtop variations by 1/6")
args = parser.parse_args()


def final_bkg_templates(hdict):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """

    rname = os.path.join(outdir, f"final_templates_lj_bkg_mtopscaled_{args.year}_{jobid}.root" if args.scale_mtop3gev else f"final_templates_lj_bkg_{args.year}_{jobid}.root")
    upfout = uproot.recreate(rname, compression=uproot.ZLIB(4)) if os.path.isfile(rname) else uproot.create(rname)

    #set_trace()
    for jmult in hdict.keys():
        for lep, histo in hdict[jmult].items():
            orig_lepdir = "muNJETS" if lep == "Muon" else "eNJETS"
            lepdir = orig_lepdir.replace("NJETS", jmult.lower())

            upfout.mkdir(lepdir)
            systypes = ["nosys"] + sorted(systematics.sys_groups[args.year].keys())
            for sys in systypes:
                if jobid == "Summer20UL_regroupedJECs":
                    if (sys == "nosys") or (sys == "JES_FlavorQCD") or (sys == "JES_RelativeBal"): continue
                    #if (sys == "nosys") or (sys == "JES_FlavorQCD") or (sys == "JES_RelativeBal") or (sys == f"JES_RelativeSample_{args.year}"): continue

                    # find histograms of associated systematics and their processes
                if sys == "nosys":
                    procs = sorted(set([key.split(f"_{sys}")[0] for key in histo.keys() if sys in key]))
                else:
                    up_sysname = systematics.sys_groups[args.year][sys][0]
                    dw_sysname = systematics.sys_groups[args.year][sys][1]
                    procs = sorted(set([key.split(f"_{up_sysname}")[0] for key in sorted(histo.keys()) if up_sysname in key])) if not dw_sysname \
                        else sorted(set([key.split(f"_{up_sysname}")[0] for key in sorted(histo.keys()) if up_sysname in key] + [key.split(f"_{dw_sysname}")[0] for key in sorted(histo.keys()) if dw_sysname in key]))

                if not procs: continue

                for proc in procs:
                    print(lep, jmult, sys, proc)
                    if sys == "nosys":
                        template, treatment = histo[f"{proc}_{sys}"]
                        if proc == "data_obs": template.clear()
                        upfout[lepdir][proc] = template.to_hist()

                    elif (sys == "deltaQCDdeltaEW") and (args.scale_mtop3gev):
                        outhname = "_".join(list(filter(None, [proc, systematics.template_sys_to_name[args.year][up_sysname]])))
                        template, treatment = histo[f"{proc}_{up_sysname}"]
                        upfout[lepdir][outhname] = template.to_hist()
                        if treatment != "raw":
                            set_trace()

                    else:
                        #set_trace()
                        up_template, treatment = histo[f"{proc}_{up_sysname}"]
                        dw_template, treatment = histo[f"{proc}_{dw_sysname}"]

                        up_outhname = "_".join(list(filter(None, [proc, systematics.template_sys_to_name[args.year][up_sysname]]))) if args.scale_mtop3gev \
                            else "_".join(list(filter(None, [proc, systematics.combine_template_sys_to_name[args.year][up_sysname]])))
                        if "LEP" in up_outhname: up_outhname = up_outhname.replace("LEP", lep.lower()) if args.scale_mtop3gev else up_outhname.replace("LEP", lep[0].lower())
                        dw_outhname = "_".join(list(filter(None, [proc, systematics.template_sys_to_name[args.year][dw_sysname]]))) if args.scale_mtop3gev \
                            else "_".join(list(filter(None, [proc, systematics.combine_template_sys_to_name[args.year][dw_sysname]])))
                        if "LEP" in dw_outhname: dw_outhname = dw_outhname.replace("LEP", lep.lower()) if args.scale_mtop3gev else dw_outhname.replace("LEP", lep[0].lower())

                        #set_trace()
                            # replace '2016APV' with '2016pre' and '2016' with '2016post' for files sent to Afiq
                        if (args.year == "2016APV") and ("2016APV" in up_outhname) and (not args.scale_mtop3gev):
                            #set_trace()
                            up_outhname = up_outhname.replace("2016APV", "2016pre")
                            dw_outhname = dw_outhname.replace("2016APV", "2016pre")
                        if (args.year == "2016") and ("2016" in up_outhname) and (not args.scale_mtop3gev):
                            up_outhname = up_outhname.replace("2016", "2016post")
                            dw_outhname = dw_outhname.replace("2016", "2016post")

                        upfout[lepdir][up_outhname] = up_template.to_hist()
                        upfout[lepdir][dw_outhname] = dw_template.to_hist()
                
                        # add extra 'chi2' histogram for variations that were smoothed/flattened
                        if treatment != "raw":
                            #print("\t", sys, treatment)
                            #if treatment == "flat": set_trace()
                            treated_up_val = up_template.values(overflow="over")[()][-1]
                            treated_dw_val = dw_template.values(overflow="over")[()][-1]
            
                            chi2_outhname = "_".join(list(filter(None, [proc, systematics.template_sys_to_name[args.year][dw_sysname].split("Down")[0], "chi2"]))) if args.scale_mtop3gev \
                                else "_".join(list(filter(None, [proc, systematics.combine_template_sys_to_name[args.year][dw_sysname].split("Down")[0], "chi2"])))
                            if "LEP" in chi2_outhname: chi2_outhname = chi2_outhname.replace("LEP", lep.lower()) if args.scale_mtop3gev else chi2_outhname.replace("LEP", lep[0].lower())

                                # replace '2016APV' with '2016pre' and '2016' with '2016post' for files sent to Afiq
                            if (args.year == "2016APV") and ("2016APV" in chi2_outhname) and (not args.scale_mtop3gev):
                                #set_trace()
                                chi2_outhname = up_outhname.replace("2016APV", "2016pre")
                            if (args.year == "2016") and ("2016" in chi2_outhname) and (not args.scale_mtop3gev):
                                chi2_outhname = up_outhname.replace("2016", "2016post")

                            sumw = np.array([0., 0., 0., 0., treated_up_val, treated_dw_val])
                            tmp_chi2_histo = orig_chi2_histo.copy()
                                ## fill bins
                            for xbin in range(sumw.size):
                                tmp_chi2_histo.values()[()][xbin] = sumw[xbin]
                            
                            #set_trace()
                            upfout[lepdir][chi2_outhname] = tmp_chi2_histo.to_hist()

    upfout.close()
    print(f"{rname} written")


def final_sig_templates(hdict):
    """
    Function that writes linearized mtt vs costheta distributions to root file.
    """

    rname = os.path.join(outdir, f"final_templates_lj_sig_kfactors_{args.year}_{jobid}.root" if args.kfactors else f"final_templates_lj_sig_{args.year}_{jobid}.root")
    upfout = uproot.recreate(rname, compression=uproot.ZLIB(4)) if os.path.isfile(rname) else uproot.create(rname)

    #set_trace()
    for jmult in hdict.keys():
        for lep, histo in hdict[jmult].items():
            orig_lepdir = "muNJETS" if lep == "Muon" else "eNJETS"
            lepdir = orig_lepdir.replace("NJETS", jmult.lower())

            upfout.mkdir(lepdir)
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
            systypes = ["nosys"]+sorted(filter(None, sorted(set([baseSys(systematics.sys_to_name[args.year][sys]) for sys in systs]))))
            for sys in systypes:
                #set_trace()
                if jobid == "Summer20UL_regroupedJECs":
                    if (sys == "nosys") or (sys == "JES_FlavorQCD") or (sys == "JES_RelativeBal"): continue
                    #if (sys == "nosys") or (sys == "JES_FlavorQCD") or (sys == "JES_RelativeBal") or (sys == f"JES_RelativeSample_{args.year}"): continue

                    # find histograms of associated systematics and their processes
                if sys == "nosys":
                    signals = sorted(set([key.split(f"_{sys}")[0] for key in histo.keys() if sys in key]))
                else:
                    up_sysname = [key for key, val in systematics.sys_to_name[args.year].items() if val == f"{sys}_UP"][0]
                    dw_sysname = [key for key, val in systematics.sys_to_name[args.year].items() if val == f"{sys}_DW"][0]
                    signals = sorted(set([key.split(f"_{up_sysname}")[0] for key in histo.keys() if up_sysname in key] + [key.split(f"_{dw_sysname}")[0] for key in histo.keys() if dw_sysname in key]))

                #set_trace()
                for signal in signals:
                    print(lep, jmult, sys, signal)
                    if "Int" in signal:
                        boson, mass, width, pI, wt = tuple(signal.split("_"))
                    else:
                        boson, mass, width, pI = tuple(signal.split("_"))

                    if args.kfactors:
                        sub_name = "_".join(["%s%s" % (boson[0], mass[1:]), "relw%s" % widthTOname(width).split("W")[-1], pI.lower(), wt]) if pI == "Int" else "_".join(["%s%s" % (boson[0], mass[1:]), "relw%s" % widthTOname(width).split("W")[-1], pI.lower()])
                    else:
                        sub_name = "_".join([boson[0], mass.lower(), width.lower(), wt]) if pI == "Int" else "_".join([boson[0], mass.lower(), width.lower(), pI.lower()])

                    if sys == "nosys":
                        template, treatment = histo[f"{signal}_{sys}"]

                        upfout[lepdir][sub_name] = template.to_hist()

                    else:
                        #set_trace()
                        up_template, treatment = histo[f"{signal}_{up_sysname}"]
                        dw_template, treatment = histo[f"{signal}_{dw_sysname}"]

                        up_outhname = "_".join(list(filter(None, [sub_name, systematics.template_sys_to_name[args.year][up_sysname]]))) if args.kfactors \
                            else "_".join(list(filter(None, [sub_name, systematics.combine_template_sys_to_name[args.year][up_sysname]])))
                        if "LEP" in up_outhname: up_outhname = up_outhname.replace("LEP", lep.lower()) if args.kfactors else up_outhname.replace("LEP", lep[0].lower())
                        dw_outhname = "_".join(list(filter(None, [sub_name, systematics.template_sys_to_name[args.year][dw_sysname]]))) if args.kfactors \
                            else "_".join(list(filter(None, [sub_name, systematics.combine_template_sys_to_name[args.year][dw_sysname]])))
                        if "LEP" in dw_outhname: dw_outhname = dw_outhname.replace("LEP", lep.lower()) if args.kfactors else dw_outhname.replace("LEP", lep[0].lower())

                            # replace '2016APV' with '2016pre' and '2016' with '2016post' for files sent to Afiq
                        if (args.year == "2016APV") and ("2016APV" in up_outhname) and (not args.kfactors):
                            #set_trace()
                            up_outhname = up_outhname.replace("2016APV", "2016pre")
                            dw_outhname = dw_outhname.replace("2016APV", "2016pre")
                        if (args.year == "2016") and ("2016" in up_outhname) and (not args.kfactors):
                            up_outhname = up_outhname.replace("2016", "2016post")
                            dw_outhname = dw_outhname.replace("2016", "2016post")

                        upfout[lepdir][up_outhname] = up_template.to_hist()
                        upfout[lepdir][dw_outhname] = dw_template.to_hist()
                
                        # add extra 'chi2' histogram for variations that were smoothed/flattened
                        if treatment != "raw":
                            #print("\t", sys, treatment)
                            #if treatment == "flat": set_trace()
                            treated_up_val = up_template.values(overflow="over")[()][-1]
                            treated_dw_val = dw_template.values(overflow="over")[()][-1]
            
                            chi2_outhname = "_".join(list(filter(None, [sub_name, systematics.template_sys_to_name[args.year][dw_sysname].split("Down")[0], "chi2"]))) if args.kfactors \
                                else "_".join(list(filter(None, [sub_name, systematics.combine_template_sys_to_name[args.year][dw_sysname].split("Down")[0], "chi2"])))
                            if "LEP" in chi2_outhname: chi2_outhname = chi2_outhname.replace("LEP", lep.lower()) if args.kfactors else chi2_outhname.replace("LEP", lep[0].lower())

                                # replace '2016APV' with '2016pre' and '2016' with '2016post' for files sent to Afiq
                            if (args.year == "2016APV") and ("2016APV" in chi2_outhname) and (not args.kfactors):
                                #set_trace()
                                chi2_outhname = up_outhname.replace("2016APV", "2016pre")
                            if (args.year == "2016") and ("2016" in chi2_outhname) and (not args.kfactors):
                                chi2_outhname = up_outhname.replace("2016", "2016post")


                            sumw = np.array([0., 0., 0., 0., treated_up_val, treated_dw_val])
                            tmp_chi2_histo = orig_chi2_histo.copy()
                                ## fill bins
                            for xbin in range(sumw.size):
                                tmp_chi2_histo.values()[()][xbin] = sumw[xbin]
                            
                            #set_trace()
                            upfout[lepdir][chi2_outhname] = tmp_chi2_histo.to_hist()

    upfout.close()
    print(f"{rname} written")


if __name__ == "__main__":
    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]

    analyzer = "htt_btag_sb_regions"
    base_bkg_template_name = f"final_templates_lj_bkg_mtopscaled_{args.year}_{jobid}.coffea" if args.scale_mtop3gev else f"final_templates_lj_bkg_{args.year}_{jobid}.coffea"
    base_sig_template_name = f"final_templates_lj_sig_kfactors_{args.year}_{jobid}.coffea" if args.kfactors else f"final_templates_lj_sig_{args.year}_{jobid}.coffea"

    input_dir = os.path.join(proj_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}", "FINAL")
    outdir = os.path.join(proj_dir, "results", f"{args.year}_{jobid}", f"Templates_{analyzer}", "FINAL")
    if not os.path.isdir(outdir):
        os.makedirs(outdir)

    if not args.only_sig:
        bkg_fname = os.path.join(input_dir, base_bkg_template_name)
        if not os.path.isfile(bkg_fname): raise ValueError("No background file found.")
        bkg_hdict = load(bkg_fname)

        orig_chi2_histo = hist.Hist("Events", hist.Bin("x_y", "x_y", np.arange(7)))
        orig_chi2_histo.fill(x_y=np.zeros(0))
        
        print("Creating final background templates")
        final_bkg_templates(bkg_hdict)

    if not args.only_bkg:
        widthTOname = lambda width : str(width).replace(".", "p")

        sig_fname = os.path.join(input_dir, base_sig_template_name)
        if not os.path.isfile(sig_fname): raise ValueError("No signal file found.")
        sig_hdict = load(sig_fname)
    
        orig_chi2_histo = hist.Hist("Events", hist.Bin("x_y", "x_y", np.arange(7)))
        orig_chi2_histo.fill(x_y=np.zeros(0))
        
        baseSys = lambda sys : "_".join(sys.split("_")[:-1])
    
        print("Creating final signal templates")
        final_sig_templates(sig_hdict)

    toc = time.time()
    print("Total time: %.1f" % (toc - tic))
