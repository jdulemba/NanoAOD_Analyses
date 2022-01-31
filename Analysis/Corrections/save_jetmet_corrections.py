from coffea.lookup_tools import extractor
from coffea.jetmet_tools import JECStack
from coffea.jetmet_tools import CorrectedJetsFactory, CorrectedMETFactory
import os
from pdb import set_trace
from coffea.util import save

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("split_uncs", choices=["total", "split", "regrouped"], help="Choose which jec uncertainty sources file to use")
#parser.add_argument("--split_uncs", action="store_true", help="Use individual jec uncertainty sources file")
parser.add_argument("--test", action="store_true", help="Output txt file named test_jetmet.txt")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
base_jobid = os.environ["base_jobid"]

def make_name_map(jstack, isMC):
    # values are hardcoded based on names in IDJet
    name_map = jstack.blank_name_map
        # set name_map keys needed for CorrectedJetsFactory
    name_map["JetPt"] = "pt"
    name_map["JetMass"] = "mass"
    name_map["JetEta"] = "eta"
    name_map["JetA"] = "area"
    if isMC: name_map["ptGenJet"] = "pt_gen"
    name_map["ptRaw"] = "pt_raw"
    name_map["massRaw"] = "mass_raw"
    name_map["Rho"] = "rho"
        # set name_map keys needed for CorrectedMETFactory
    name_map["METpt"] = "pt"
    name_map["METphi"] = "phi"
    name_map["JetPhi"] = "phi"
    name_map["UnClusteredEnergyDeltaX"] = "MetUnclustEnUpDeltaX"
    name_map["UnClusteredEnergyDeltaY"] = "MetUnclustEnUpDeltaY"

    return name_map

#jec_levels_MC = ["L1FastJet", "L2Relative", "L3Absolute"]
#jec_levels_Data = ["L1FastJet", "L2L3Residual", "L2Relative", "L3Absolute"]


years_to_run = ["2016APV", "2016", "2017", "2018"]
if args.test: years_to_run = ["2018"]

jet_corrections = {year:{"MC" : {}, "DATA" : {}} for year in years_to_run}
cfg_dict = {
    "2016" : {
        "MC" : {
            "ERA" : [""],
        },
        "DATA" : {
            "ERA" : ["FGH"],
        }
    },
    "2016APV" : {
        "MC" : {"ERA" : [""]},
        "DATA" : {
            "ERA" : ["BCD", "EF"],
        }
    },
    "2017" : {
        "MC" : {
            "ERA" : [""],
        },
        "DATA" : {
            "ERA" : ["B", "C", "D", "E", "F"],
        }
    },
    "2018" : {
        "MC" : {
            "ERA" : [""],
        },
        "DATA" : {
            "ERA" : ["A", "B", "C", "D"],
        }
    },
}

#set_trace()
jec_uncs_dict = {"total" : "Uncertainty", "split" : "UncertaintySources", "regrouped" : "RegroupedV2"}
jec_unc_type = jec_uncs_dict[args.split_uncs]
#jec_unc_type = "UncertaintySources" if args.split_uncs else "Uncertainty"

dtypes_to_run = ["MC"] if args.split_uncs == "regrouped" else ["MC", "DATA"]

for year in years_to_run:
    for dtype in dtypes_to_run:
    #for dtype in ["DATA", "MC"]:
    #    if (args.split_uncs == "regrouped") and (dtype == "DATA"): continue
        for era in cfg_dict[year][dtype]["ERA"]:
            print(f"\tMaking corrections for {year} {dtype} {era}")
            directory = os.path.join(proj_dir, "inputs", "data", base_jobid, "Jet_Corrections", year, dtype, era)
            jec_stack_names = [os.path.basename(fname).split(".")[0] for fname in os.listdir(directory)]

                # find jec uncertainty files and remove whichever one isn't wanted based on args.split_uncs
            unc_files = [fname for fname in jec_stack_names if "Uncertainty" in fname]
            junc_to_remove = [fname for fname in unc_files if not any(jec_unc_type == substring for substring in fname.split("_"))]
            for jtr in junc_to_remove:
                jec_stack_names.remove(jtr)
            #if len(junc_to_remove) > 1: raise ValueError(f"Multiple jec uncertainty files found to be removed for {year} {dtype} {era}")
            #jec_stack_names.remove(junc_to_remove[0])

            fnames = [os.path.join(directory, fname) for fname in os.listdir(directory) if os.path.basename(fname).split(".")[0] in jec_stack_names] # add files that are only in jec_stack_names

            ext = extractor()
            ext.add_weight_sets([f"* * {fname}" for fname in fnames])
            ext.finalize()
            evaluator = ext.make_evaluator()

            jec_inputs = {name: evaluator[name] for name in dir(evaluator)} 
            jec_stack = JECStack(jec_inputs)

                # make jet and met factory
            name_map = make_name_map(jec_stack, isMC=(dtype == "MC"))
            jet_factory = CorrectedJetsFactory(name_map, jec_stack)
            met_factory = CorrectedMETFactory(name_map)
            if dtype == "DATA":
                jet_corrections[year][dtype][era] = {
                    "JetsFactory" : jet_factory,
                    "METFactory" : met_factory,
                }
            else:
                jet_corrections[year][dtype] = {
                    "JetsFactory" : jet_factory,
                    "METFactory" : met_factory,
                }
            print(f"Jet corrections for {year} {dtype} {era} saved")

#set_trace()
if args.split_uncs == "total":
    fname = os.path.join(proj_dir, "Corrections", base_jobid, "JetMETCorrections_Total.coffea")
elif args.split_uncs == "split":
    fname = os.path.join(proj_dir, "Corrections", base_jobid, "JetMETCorrections_UncSources.coffea")
elif args.split_uncs == "regrouped":
    fname = os.path.join(proj_dir, "Corrections", base_jobid, "JetMETCorrections_RegroupedV2.coffea")
#fname = os.path.join(proj_dir, "Corrections", base_jobid, "JetMETCorrections_UncSources.coffea") if args.split_uncs else os.path.join(proj_dir, "Corrections", base_jobid, "JetMETCorrections.coffea")
if args.test: fname = os.path.join(proj_dir, "test_jetmet.coffea")

save(jet_corrections, fname)
print(f"\n{fname} written")
