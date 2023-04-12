import numpy as np
from coffea.util import save, load
from pdb import set_trace
import os
import python.LeptonSF as lepSF

proj_dir = os.environ["PROJECT_DIR"]
jobid = os.environ["jobid"]
base_jobid = os.environ["base_jobid"]

outdir = os.path.join(proj_dir, "Corrections", jobid)
if not os.path.isdir(outdir):
    os.makedirs(outdir)


import Utilities.prettyjson as prettyjson
cfg_pars = prettyjson.loads(open(os.path.join(proj_dir, "cfg_files", f"cfg_pars_{jobid}.json")).read())
inputLepSFs = load(os.path.join(proj_dir, "Corrections", base_jobid, cfg_pars["corrections"]["input_lepton"]))

#set_trace()
outputSFs = {}
for year in inputLepSFs.keys():
    print(year)
    outputSFs[year] = lepSF.LeptonSF(inputLepSFs[year], debug=True)

#set_trace()
lepSF_name = os.path.join(outdir, f"test_lepton_{jobid}.coffea")
save(outputSFs, lepSF_name)    
print(f"{lepSF_name} written")

