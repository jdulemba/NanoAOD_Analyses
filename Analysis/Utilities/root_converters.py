"""
This file was created by Joseph Dulemba on Thursday 3 February 2022.

It's based on the coffea file https://github.com/CoffeaTeam/coffea/blob/master/coffea/lookup_tools/root_converters.py
but converts TGraph objects into dictionaries.
"""
import uproot
import re
from pdb import set_trace

cycle = re.compile(r";\d+")
def killcycle(s, cycle):
    return cycle.sub("", s)


histTypes = ["TH1D", "TH1F", "TH2D", "TH2F", "TH3D", "TH3F"]
graphTypes = ["TGraphAsymmErrors", "TGraph2D"]


def convert_TGraph_root_file(file):
    #set_trace()
    converted_file = {}
    fin = uproot.open(file.strip())
    for path, item in fin.iteritems(recursive=True):
        if isinstance(item, uproot.ReadOnlyDirectory):
            continue
        nicepath = killcycle(path, cycle)
        rootclass = item.classname
        if rootclass in graphTypes:
            json = item.tojson()
            points = {(json["fX"][idx], json["fY"][idx]): json["fZ"][idx] for idx in range(len(json["fX"]))}
            converted_file[json["fName"]] = points
        elif rootclass in histTypes:
            print("Use coffea's convert_histo_root_file for TH1/2D histos.")
            continue
        else:
            print(f"Class {rootclass} not supported.")
    return converted_file
