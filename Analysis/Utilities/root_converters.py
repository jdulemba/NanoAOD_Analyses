"""
This file was created by Joseph Dulemba on Thursday 3 February 2022.

It's based on the coffea file https://github.com/CoffeaTeam/coffea/blob/master/coffea/lookup_tools/root_converters.py
but converts TGraph objects into dictionaries.
"""
import uproot
import re
from pdb import set_trace
import numpy as np

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
            xcenters, ycenters = (np.unique(json["fX"]), np.unique(json["fY"]))
                # format values into (xdim, ydim) array so vals_2d[0] corresponds to all of the ybin values for a single xbin value
            vals_2d = np.reshape(np.array(json["fZ"]),  (xcenters.size, ycenters.size))
            xdiffs = np.array([(xcenters[idx+1]-xcenters[idx])/2 for idx in range(xcenters.size - 1)])
            ydiffs = np.array([(ycenters[idx+1]-ycenters[idx])/2 for idx in range(ycenters.size - 1)])
                # make correct bin edges by making original fX and fY values the centers
            xedges = np.array( [xcenters[0]-xdiffs[0]] + [xcenters[idx+1]-xdiffs[idx] for idx in range(xcenters.size -1)] + [xcenters[-1]+xdiffs[-1]])
            yedges = np.array( [ycenters[0]-ydiffs[0]] + [ycenters[idx+1]-ydiffs[idx] for idx in range(ycenters.size -1)] + [ycenters[-1]+ydiffs[-1]])
            converted_file[(nicepath, "dense_lookup")] = vals_2d, (xedges, yedges)
        elif rootclass in histTypes:
            print("Use coffea's convert_histo_root_file for TH1/2D histos.")
            continue
        else:
            print(f"Class {rootclass} not supported.")
    return converted_file
