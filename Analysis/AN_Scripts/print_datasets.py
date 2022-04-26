import os
from pdb import set_trace
import Utilities.prettyjson as prettyjson

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--type", default="data", nargs="?", choices=["SM", "signal", "data", "all"], help="specify which sample type to print")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
base_jobid = os.environ["base_jobid"]
jobid = os.environ["jobid"]

    # only print signal MC info
if (args.type == "signal") or (args.type == "all"):
    import itertools
    from string import Template
    name2val = lambda x: float(x.replace("pc","").replace("p", "."))

    samples_list = prettyjson.loads(open("%s/inputs/samples_2016.json" % proj_dir).read())
    
    samples_dict = {}
    for sample in samples_list:
        name = sample.pop("name")
        samples_dict[name] = sample

        ## make table for A
    A_output = "\multirow{2}{*}{Parity} & \multirow{2}{*}{$\mathsf{m_{A}}$ [GeV]} & \multirow{2}{*}{$\Gamma_{\mathsf{A}}$ [\% $\mathsf{m_{A}}$]} & \multicolumn{2}{c |}{LO $\sigma$ [pb]} & \multirow{2}{*}{$\mathsf{k_{R}}$} \\\ \n"
    A_output += " & & & Resonance & Interference & \\\ \n\hline \n"
    H_output = "\multirow{2}{*}{Parity} & \multirow{2}{*}{$\mathsf{m_{H}}$ [GeV]} & \multirow{2}{*}{$\Gamma_{\mathsf{H}}$ [\% $\mathsf{m_{H}}$]} & \multicolumn{2}{c |}{LO $\sigma$ [pb]} & \multirow{2}{*}{$\mathsf{k_{R}}$} \\\ \n"
    H_output += " & & & Resonance & Interference & \\\ \n\hline \n"
    for sig_point in itertools.product(["M400","M500", "M600", "M750"], ["W2p5", "W5", "W10", "W25"]):
        signal = "_".join(sig_point)

        mass = sig_point[0][1:]
        width = name2val(sig_point[1][1:])

        A_res_xsec = samples_dict["AtoTT_%s_Res" % signal]["xsection"]
        A_int_xsec = samples_dict["AtoTT_%s_Int" % signal]["xsection"]
        H_res_xsec = samples_dict["HtoTT_%s_Res" % signal]["xsection"]
        H_int_xsec = samples_dict["HtoTT_%s_Int" % signal]["xsection"]

        A_template = Template("\multirow{16}{*}{A} & \multirow{4}{*}{$MASS} & $WIDTH & $XSEC_A_RES & $XSEC_A_INT & $KFACTOR \\\\ \n")
        A_output += A_template.substitute(MASS=mass, WIDTH=width, XSEC_A_RES="{:.3f}".format(A_res_xsec), XSEC_A_INT="{:.3f}".format(A_int_xsec), KFACTOR="---")
        H_template = Template("\multirow{16}{*}{H} & \multirow{4}{*}{$MASS} & $WIDTH & $XSEC_H_RES & $XSEC_H_INT & $KFACTOR \\\\ \n")
        H_output += H_template.substitute(MASS=mass, WIDTH=width, XSEC_H_RES="{:.3f}".format(H_res_xsec), XSEC_H_INT="{:.3f}".format(H_int_xsec), KFACTOR="---")

    print(A_output)
    print(H_output)
    sig_datasets_info = open("%s/AN_Scripts/sig_datasets.txt" % proj_dir, "w")
    sig_datasets_info.write(A_output)
    sig_datasets_info.write(H_output)
    sig_datasets_info.close()
    print("%s/AN_Scripts/sig_datasets.txt written" % proj_dir)


lumi_list = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"{base_jobid}_lumis_data.json")).read())

years = ["2016APV", "2016", "2017", "2018"]
for year in years:
    samples_list = prettyjson.loads(open(os.path.join(proj_dir, "inputs", f"samples_{year}_{base_jobid}.json")).read())
    
    samples_dict = {}
    for sample in samples_list:
        name = sample.pop("name")
        samples_dict[name] = sample

    if (args.type == "data") or (args.type == "all"):
        lumiMasks = {}
        set_trace()
        data_output = "Data set & Integrated luminosity ($\\fbinv$)\\\\ \n\hline \n"
        data_output += f"{year} & {lumi_list[year]['Muons']/1000:.2f} ({lumi_list[year]['Electrons']/1000:.2f})\\\\ \n"
    
    if (args.type == "SM") or (args.type == "all"):
        sm_output = "{\small Data set} & {\small $\sigma$ [pb]} & {\small $\Sigma$genweights} \\\\ \n\hline\hline \n"
        #sm_output = "{\small Data set} & {\small $\sigma$ [pb]} & {\small Weights} \\\\ \n\hline \n"


    for sample in samples_dict.keys():
            # only print SM MC info
        if (args.type == "SM") or (args.type == "all"):
            if (sample.startswith("AtoTT")) or (sample.startswith("HtoTT")): continue
            if "data" in sample: continue
            #set_trace()
            dbs_name = samples_dict[sample]["DBSName"] # get full DBS name
            name_to_use = dbs_name.split("/")[1] if dbs_name.startswith("/") else dbs_name.split("/")[0]
            name_to_use = name_to_use.replace("_", "\_")
            xsec = samples_dict[sample]["xsection"]
            #nwts = "-" if not os.path.isfile(os.path.join(proj_dir, "inputs", f"{year}_{base_jobid}", f"{sample}.meta.json")) else prettyjson.loads(open(os.path.join(proj_dir, "inputs",f"{year}_{base_jobid}", f"{sample}.meta.json")).read())["nEvents"]
            nwts = "-" if not os.path.isfile(os.path.join(proj_dir, "inputs", f"{year}_{base_jobid}", f"{sample}.meta.json")) else prettyjson.loads(open(os.path.join(proj_dir, "inputs",f"{year}_{base_jobid}", f"{sample}.meta.json")).read())["sumGenWeights"]
        
            sm_output += f"{{\small {name_to_use}}} & {{\small {xsec}}} & {{\small {nwts}}} \\\\ \n"


            # only print data info
        if (args.type == "data") or (args.type == "all"):
            if not "data" in sample: continue
            #set_trace()
            lepton = sample.split("data_Single")[-1].split("_")[0]+"s"
            period = sample.split("2016")[-1] if year == "2016APV" else sample.split(year)[-1]
            if period not in lumi_list["Ind"][year][lepton].keys():
                print(f"{sample} not found in list of lumis used. Make sure this is correct.")
                continue
            lumi = lumi_list["Ind"][year][lepton][period]/1000
            dbs_name = samples_dict[sample]["DBSName"] # get full DBS name
            dbs_name = ("/".join(dbs_name.split("/")[:-1])).replace("_", "\_") # remove "NANOAOD" part of dbs name
            set_trace()
            data_output += f"{dbs_name} & {lumi:.2f} \\\\ \n"

                # add lumimask used for datset
            lumiMasks[sample] = samples_dict[sample]["lumimask"]
    
    if (args.type == "data") or (args.type == "all"):
            # add lumimask to output
        if len(sorted(set(lumiMasks.values()))) == 1:
            lumimask_val = os.path.basename(sorted(set(lumiMasks.values()))[0]).replace("_", "\_")
            data_output += f"\n\n{lumimask_val}\n\n"

        print(data_output)
        data_fname = os.path.join(proj_dir, "AN_Scripts", f"{year}_data_datasets.txt")
        data_datasets_lumi = open(data_fname, "w")
        data_datasets_lumi.write(data_output)
        data_datasets_lumi.close()
        print(f"{data_fname} written")

    if (args.type == "SM") or (args.type == "all"):
        sm_fname = os.path.join(proj_dir, "AN_Scripts", f"{year}_SM_datasets.txt")
        SM_datasets_info = open(sm_fname, "w")
        SM_datasets_info.write(sm_output)
        SM_datasets_info.close()
        print(f"{sm_fname} written")
