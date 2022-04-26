from pdb import set_trace
import fnmatch
import os
import matplotlib.pyplot as plt

base_jobid = os.environ["base_jobid"]


def create_bkg_groups(mode):
    if mode == "dataset":
        #dataset_groups = {
        groupings = {
            year : {
                "EWK" : ["[WZ][WZ]", "ZJets*", "WJets_HT*", "tt[WZ]*"] if base_jobid == "Summer20UL" else ["[WZ][WZ]", "[WZ]Jets*", "tt[WZ]*"],
                "ttJets_right" : ["ttJets*_right"],
                "ttJets_matchable" : ["ttJets*_matchable"],
                "ttJets_unmatchable" : ["ttJets*_unmatchable"],
                "ttJets_sl_tau" : ["ttJets*_sl_tau"],
                "ttJets_other" : ["ttJets*_other"],
                "ttJets" : ["ttJets", "ttJets_PS"] if (base_jobid == "NanoAODv6" and year == "2016") else ["ttJetsSL", "ttJetsDiLep", "ttJetsHad"],
                "ttJets_Sys" : ["ttJets_*DOWN", "ttJets_*UP"] if (base_jobid == "NanoAODv6" and year == "2016") else ["ttJetsSL_*DOWN", "ttJetsDiLep_*DOWN", "ttJetsHad_*DOWN", "ttJetsSL_*UP", "ttJetsDiLep_*UP", "ttJetsHad_*UP"],
                "singlet" : ["single*"],
                "QCD" : ["QCD*"],
                "data" : ["data_Single*"],
            }
            for year in ["2016APV", "2016", "2017", "2018"]
        }

    elif mode == "templates":
            ## dataset groupings for making templates to be used in fit
        groupings = {
        #template_groups = {
            year: {
                "EWQCD" : ["QCD*", "[WZ][WZ]", "ZJets*", "WJets_HT*", "tt[WZ]*"] if base_jobid == "Summer20UL" else ["QCD*", "[WZ][WZ]", "[WZ]Jets*", "tt[WZ]*"],
                "TT" : ["ttJets*_right", "ttJets*_matchable", "ttJets*_unmatchable", "ttJets*_sl_tau", "ttJets*_other"],
                #"QCD" : ["QCD*"],
                #"VV" : ["[WZ][WZ]"],
                #"TTV" : ["tt[WZ]*"],
                #"W" : ["WJets_HT*"],
                #"DY" : ["ZJets*"],
                "TQ" : ["single*_tchannel*"],
                "TW" : ["single*_tW*"],
                "TB" : ["single*_schannel*"],
                "data_obs" : ["data_Single*"],
            }
            for year in ["2016APV", "2016", "2017", "2018"]
        }
    else:
        raise ValueError("Can only choose between 'dataset' and 'templates' for background groupings options")

    return groupings

def create_sig_groups(mode):
    from itertools import product
    import numpy as np

    groupings = {year: {} for year in ["2016APV", "2016", "2017", "2018"]}
    if mode == "MC_indiv":
        #MC_masses = np.array([800.])
        MC_masses = np.array([365., 400., 500., 600., 800., 1000.])
        MC_widths = np.array([2.5, 10.0, 25.0])
        #MC_widths = np.array([10.0])
        for year in ["2016APV", "2016", "2017", "2018"]:
            for bundle in product(MC_masses, [str(width).replace(".", "p") for width in sorted(MC_widths)]):
                groupings[year]["AtoTTJetsSL_M%d_W%s_Int_neg" % bundle] = ["AtoTTJetsSL_M%d_W%s_Int_neg" % bundle]
                groupings[year]["AtoTTJetsSL_M%d_W%s_Int_pos" % bundle] = ["AtoTTJetsSL_M%d_W%s_Int_pos" % bundle]
                groupings[year]["AtoTTJetsSL_M%d_W%s_Res" % bundle]     = ["AtoTTJetsSL_M%d_W%s_Res" % bundle]
                groupings[year]["HtoTTJetsSL_M%d_W%s_Int_neg" % bundle] = ["HtoTTJetsSL_M%d_W%s_Int_neg" % bundle]
                groupings[year]["HtoTTJetsSL_M%d_W%s_Int_pos" % bundle] = ["HtoTTJetsSL_M%d_W%s_Int_pos" % bundle]
                groupings[year]["HtoTTJetsSL_M%d_W%s_Res" % bundle]     = ["HtoTTJetsSL_M%d_W%s_Res" % bundle]
                groupings[year]["AtoTTJetsDiLep_M%d_W%s_Int_neg" % bundle] = ["AtoTTJetsDiLep_M%d_W%s_Int_neg" % bundle]
                groupings[year]["AtoTTJetsDiLep_M%d_W%s_Int_pos" % bundle] = ["AtoTTJetsDiLep_M%d_W%s_Int_pos" % bundle]
                groupings[year]["AtoTTJetsDiLep_M%d_W%s_Res" % bundle]     = ["AtoTTJetsDiLep_M%d_W%s_Res" % bundle]
                groupings[year]["HtoTTJetsDiLep_M%d_W%s_Int_neg" % bundle] = ["HtoTTJetsDiLep_M%d_W%s_Int_neg" % bundle]
                groupings[year]["HtoTTJetsDiLep_M%d_W%s_Int_pos" % bundle] = ["HtoTTJetsDiLep_M%d_W%s_Int_pos" % bundle]
                groupings[year]["HtoTTJetsDiLep_M%d_W%s_Res" % bundle]     = ["HtoTTJetsDiLep_M%d_W%s_Res" % bundle]
            for bundle in product(np.array([400.]), ["5p0"]):
                groupings[year]["AtoTTJetsSL_M%d_W%s_Int_neg" % bundle] = ["AtoTTJetsSL_M%d_W%s_Int_neg" % bundle]
                groupings[year]["AtoTTJetsSL_M%d_W%s_Int_pos" % bundle] = ["AtoTTJetsSL_M%d_W%s_Int_pos" % bundle]
                groupings[year]["AtoTTJetsSL_M%d_W%s_Res" % bundle]     = ["AtoTTJetsSL_M%d_W%s_Res" % bundle]
                groupings[year]["HtoTTJetsSL_M%d_W%s_Int_neg" % bundle] = ["HtoTTJetsSL_M%d_W%s_Int_neg" % bundle]
                groupings[year]["HtoTTJetsSL_M%d_W%s_Int_pos" % bundle] = ["HtoTTJetsSL_M%d_W%s_Int_pos" % bundle]
                groupings[year]["HtoTTJetsSL_M%d_W%s_Res" % bundle]     = ["HtoTTJetsSL_M%d_W%s_Res" % bundle]
                groupings[year]["AtoTTJetsDiLep_M%d_W%s_Int_neg" % bundle] = ["AtoTTJetsDiLep_M%d_W%s_Int_neg" % bundle]
                groupings[year]["AtoTTJetsDiLep_M%d_W%s_Int_pos" % bundle] = ["AtoTTJetsDiLep_M%d_W%s_Int_pos" % bundle]
                groupings[year]["AtoTTJetsDiLep_M%d_W%s_Res" % bundle]     = ["AtoTTJetsDiLep_M%d_W%s_Res" % bundle]
                groupings[year]["HtoTTJetsDiLep_M%d_W%s_Int_neg" % bundle] = ["HtoTTJetsDiLep_M%d_W%s_Int_neg" % bundle]
                groupings[year]["HtoTTJetsDiLep_M%d_W%s_Int_pos" % bundle] = ["HtoTTJetsDiLep_M%d_W%s_Int_pos" % bundle]
                groupings[year]["HtoTTJetsDiLep_M%d_W%s_Res" % bundle]     = ["HtoTTJetsDiLep_M%d_W%s_Res" % bundle]
        
    elif mode == "MC_combined":
        MC_masses = np.array([365., 400., 500., 600., 800., 1000.])
        #MC_widths = np.array([2.5, 5.0, 10.0, 25.0])
        MC_widths = np.array([2.5, 10.0, 25.0])
        for year in ["2016APV", "2016", "2017", "2018"]:
            for bundle in product(MC_masses, [str(width).replace(".", "p") for width in sorted(MC_widths)]):
                groupings[year]["AtoTT_M%d_W%s_Int_neg" % bundle] = ["AtoTT*_M%d_W%s_Int_neg" % bundle]
                groupings[year]["AtoTT_M%d_W%s_Int_pos" % bundle] = ["AtoTT*_M%d_W%s_Int_pos" % bundle]
                groupings[year]["AtoTT_M%d_W%s_Res" % bundle]     = ["AtoTT*_M%d_W%s_Res" % bundle]
                groupings[year]["HtoTT_M%d_W%s_Int_neg" % bundle] = ["HtoTT*_M%d_W%s_Int_neg" % bundle]
                groupings[year]["HtoTT_M%d_W%s_Int_pos" % bundle] = ["HtoTT*_M%d_W%s_Int_pos" % bundle]
                groupings[year]["HtoTT_M%d_W%s_Res" % bundle]     = ["HtoTT*_M%d_W%s_Res" % bundle]
            for bundle in product(np.array([400.]), ["5p0"]):
                groupings[year]["AtoTT_M%d_W%s_Int_neg" % bundle] = ["AtoTT*_M%d_W%s_Int_neg" % bundle]
                groupings[year]["AtoTT_M%d_W%s_Int_pos" % bundle] = ["AtoTT*_M%d_W%s_Int_pos" % bundle]
                groupings[year]["AtoTT_M%d_W%s_Res" % bundle]     = ["AtoTT*_M%d_W%s_Res" % bundle]
                groupings[year]["HtoTT_M%d_W%s_Int_neg" % bundle] = ["HtoTT*_M%d_W%s_Int_neg" % bundle]
                groupings[year]["HtoTT_M%d_W%s_Int_pos" % bundle] = ["HtoTT*_M%d_W%s_Int_pos" % bundle]
                groupings[year]["HtoTT_M%d_W%s_Res" % bundle]     = ["HtoTT*_M%d_W%s_Res" % bundle]

    elif mode == "MEreweight_indiv":
        MEreweight_masses = np.array([365, 380, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000])
        MEreweight_widths = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 8.0, 10.0, 13.0, 15.0, 18.0, 21.0, 25.0])
        #MEreweight_masses = np.arange(365., 1005., 5.)
        #MEreweight_widths = np.arange(0.5, 25.5, 0.5)
        for year in ["2016APV", "2016", "2017", "2018"]:
            for bundle in product(MEreweight_masses, [str(width).replace(".", "p") for width in sorted(MEreweight_widths)]):
                groupings[year]["AtoTTJetsSL_M%d_W%s_Int_neg" % bundle] = ["AtoTTJetsSL_M%d_W%s_Int_neg" % bundle]
                groupings[year]["AtoTTJetsSL_M%d_W%s_Int_pos" % bundle] = ["AtoTTJetsSL_M%d_W%s_Int_pos" % bundle]
                groupings[year]["AtoTTJetsSL_M%d_W%s_Res" % bundle]     = ["AtoTTJetsSL_M%d_W%s_Res" % bundle]
                groupings[year]["HtoTTJetsSL_M%d_W%s_Int_neg" % bundle] = ["HtoTTJetsSL_M%d_W%s_Int_neg" % bundle]
                groupings[year]["HtoTTJetsSL_M%d_W%s_Int_pos" % bundle] = ["HtoTTJetsSL_M%d_W%s_Int_pos" % bundle]
                groupings[year]["HtoTTJetsSL_M%d_W%s_Res" % bundle]     = ["HtoTTJetsSL_M%d_W%s_Res" % bundle]
                groupings[year]["AtoTTJetsDiLep_M%d_W%s_Int_neg" % bundle] = ["AtoTTJetsDiLep_M%d_W%s_Int_neg" % bundle]
                groupings[year]["AtoTTJetsDiLep_M%d_W%s_Int_pos" % bundle] = ["AtoTTJetsDiLep_M%d_W%s_Int_pos" % bundle]
                groupings[year]["AtoTTJetsDiLep_M%d_W%s_Res" % bundle]     = ["AtoTTJetsDiLep_M%d_W%s_Res" % bundle]
                groupings[year]["HtoTTJetsDiLep_M%d_W%s_Int_neg" % bundle] = ["HtoTTJetsDiLep_M%d_W%s_Int_neg" % bundle]
                groupings[year]["HtoTTJetsDiLep_M%d_W%s_Int_pos" % bundle] = ["HtoTTJetsDiLep_M%d_W%s_Int_pos" % bundle]
                groupings[year]["HtoTTJetsDiLep_M%d_W%s_Res" % bundle]     = ["HtoTTJetsDiLep_M%d_W%s_Res" % bundle]
    
    elif mode == "MEreweight_combined":
        MEreweight_masses = np.array([365, 380, 400, 425, 450, 475, 500, 525, 550, 575, 600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000])
        MEreweight_widths = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 8.0, 10.0, 13.0, 15.0, 18.0, 21.0, 25.0])
        #MEreweight_widths = np.arange(0.5, 25.5, 0.5)
        #MEreweight_masses = np.arange(365., 605., 5.)
        #MEreweight_widths = np.arange(0.5, 8.5, 0.5)
        #MEreweight_masses = np.array([425., 575.])
        #MEreweight_widths = np.array([0.5, 8.0])
        #MEreweight_masses = np.arange(600., 1005., 25.)
        #MEreweight_masses = np.array([600, 625, 650, 675, 700, 725, 750, 775, 800, 825, 850, 875, 900, 925, 950, 975, 1000])
        #MEreweight_widths = np.array([2.5])
        for year in ["2016APV", "2016", "2017", "2018"]:
            for bundle in product(MEreweight_masses, [str(width).replace(".", "p") for width in sorted(MEreweight_widths)]):
                groupings[year]["AtoTT_M%d_W%s_Int_neg" % bundle] = ["AtoTT*_M%d_W%s_Int_neg" % bundle]
                groupings[year]["AtoTT_M%d_W%s_Int_pos" % bundle] = ["AtoTT*_M%d_W%s_Int_pos" % bundle]
                groupings[year]["AtoTT_M%d_W%s_Res" % bundle]     = ["AtoTT*_M%d_W%s_Res" % bundle]
                groupings[year]["HtoTT_M%d_W%s_Int_neg" % bundle] = ["HtoTT*_M%d_W%s_Int_neg" % bundle]
                groupings[year]["HtoTT_M%d_W%s_Int_pos" % bundle] = ["HtoTT*_M%d_W%s_Int_pos" % bundle]
                groupings[year]["HtoTT_M%d_W%s_Res" % bundle]     = ["HtoTT*_M%d_W%s_Res" % bundle]

    else:
        raise ValueError("Can only choose between 'MC_indiv', 'MC_combined', 'MEreweight_indiv', and 'MEreweight_combined' for signal groupings options")

    return groupings


def get_styles(sample, styles):
    best_pattern = ""
    for pattern, style_dict in styles.items():
        if fnmatch.fnmatch(sample, pattern):
            if len(pattern) > len(best_pattern):
                best_pattern = pattern
    if best_pattern:
        return styles[best_pattern]["facecolor"], styles[best_pattern]["name"]
    else:
        return "r", sample

def get_label(sample, styles):
    best_pattern = ""
    for pattern, style_dict in styles.items():
        if fnmatch.fnmatch(sample, pattern):
            if len(pattern) > len(best_pattern):
                best_pattern = pattern
    if best_pattern:
        return styles[best_pattern]["name"]
    else:
        return sample

#def get_group(sample, styles=dataset_groups):
def get_group(sample, styles=create_bkg_groups("dataset")):
    best_pattern = ""
    for group_name, patterns in styles.items():
        if any([fnmatch.fnmatch(sample, pattern) for pattern in patterns]):
            best_pattern = group_name

    if best_pattern:
        return best_pattern
    else:
        print(f"Pattern not found for {sample}")
        return sample

def make_dataset_groups(lepton, year, samples=[], bkgdict="dataset", sigdict="MC_indiv"):
    proj_dir = os.environ["PROJECT_DIR"]
    jobid = os.environ["jobid"]

    if not samples:
            ## get file that has names for all datasets to use
        fpath = os.path.join(proj_dir, "inputs", "%s_%s" % (year, jobid), "analyzer_inputs.txt")
        if not os.path.isfile(fpath):
            raise IOError(f"File {fpath} not found")

        txt_file = open(fpath, "r")
        samples = [sample.strip("\n") for sample in txt_file if not sample.startswith("#")]
        samples = [sample for sample in samples if not (sample.startswith("data") and lepton not in sample)] # get rid of data samples that don"t correspond to lepton chosen

    groups_dict = {}
    bkg_groupings = create_bkg_groups(mode=bkgdict)
    sig_groupings = create_sig_groups(mode=sigdict)
    groups_dict.update(bkg_groupings)
    [groups_dict[year].update(sig_groupings[year]) for year in groups_dict.keys()]

    groupings = {}
    for group_name, patterns in groups_dict[year].items():
        flist = []
        for sample in samples:
            if any([fnmatch.fnmatch(sample, pattern) for pattern in patterns]):
                flist.append(sample)
        if flist:
            groupings[group_name] = flist

    return groupings


def add_coffea_files(input_files):
    from coffea.util import load
    from tqdm import tqdm

    output_acc = load(input_files[0])
    for idx in tqdm(range(1, len(input_files))):
        try:
            output_acc.add(load(input_files[idx]))
        except:
            raise ValueError("File %s (number %i) could not be added" % (input_files[idx], idx))

    return output_acc

def save_accumulator(accumulator, output_fname):
    from coffea.util import save
    save(accumulator, output_fname)
    print(f"{output_fname} written")


def print_table(lines, filename="", separate_header=True, header_line=0, print_output=False):
    """ Prints a formatted table given a 2 dimensional array """
    #Count the column width
    widths = []
    for line in lines:
        for i, size in enumerate( [ len(x) for x in line ] ):
            while i >= len(widths):
                widths.append(0)
            if size > widths[i]:
                widths[i] = size
    
    #Generate the format string to pad the columns
    print_string = ""
    for i, width in enumerate(widths):
        print_string += "{" + str(i) + ":" + str(width) + "} | "
    if (len(print_string) == 0):
        return
    print_string = print_string[:-3]
    
    if print_output: ## print data
        for i,line in enumerate(lines):
            print(print_string.format(*line))
            if (i == header_line and separate_header):
                print("-"*(sum(widths)+3*(len(widths)-1)))

    if filename != "": # save to filename    
        with open(filename, "w") as f:
            #Write the data into filename
            for i,line in enumerate(lines):
                print(print_string.format(*line), file=f)
                if (i == header_line and separate_header):
                    print("-"*(sum(widths)+3*(len(widths)-1)), file=f)

