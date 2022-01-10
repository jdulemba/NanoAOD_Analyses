from coffea.util import load
import uproot3
from pdb import set_trace
import os
from Utilities import HistExport
from Utilities import Plotter
from coffea import hist
import numpy as np

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("input_file", help="Input filename")
args = parser.parse_args()

proj_dir = os.environ["PROJECT_DIR"]
input_dir = input_dir = os.path.join(proj_dir, "NNLO_files")

input_file = os.path.join(input_dir, args.input_file)
if not os.path.isfile(input_file):
    raise ValueError(f"File {input_file} not found")

cfile = load(input_file)
# hardcoded values known beforehand
var = "mtt_vs_thad_ctstar"
hist_template = hist.Hist(
    "Events",
    hist.Bin("mtt", "mtt", np.around(np.concatenate((np.arange(300., 1010., 10), np.arange(1050., 3550., 50))), decimals=0)),
    hist.Bin("thad_cstar", "ctstar", np.around(np.linspace(-1., 1., 21), decimals=2))
)
    
set_trace()
outrname = os.path.join(input_dir, args.input_file.replace("coffea", "root"))
upfout = uproot3.recreate(outrname, compression=uproot3.ZLIB(4)) if os.path.isfile(outrname) else uproot3.create(outrname)
for year in cfile.keys():
    corr = cfile[year][var]
    tmp_hist = Plotter.np_array_TO_hist(sumw=corr._values, sumw2=np.zeros(corr._values.shape), hist_template=hist_template)
    upfout[f"NNLO_to_{year}_Ratios_Summer20UL"] = HistExport.export2d(tmp_hist)
upfout.close()
print(f"{outrname} written")
