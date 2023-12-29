import time
tic = time.time()

# matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend("agg")
from matplotlib import rcParams
rcParams["font.size"] = 20
rcParams["savefig.format"] = "pdf"
rcParams["savefig.bbox"] = "tight"

import numpy as np
from pdb import set_trace
import uproot
import os
from coffea import hist

jobid = os.environ["jobid"]
proj_dir = os.environ["PROJECT_DIR"]
plot_outdir = os.environ["plots_dir"]

basedir = "root://cmseos.fnal.gov//store/user/jdulemba/HeavyHiggsFitter/Thesis_constrainedYt_EtaT_v1"

outdir = os.path.join(plot_outdir, jobid, "Limits_htt_btag_sb_regions", os.path.basename(basedir).split(f"{jobid}_")[-1], "GOF")
if not os.path.isdir(outdir):
    os.makedirs(outdir)


def make_plot(rname, label):
    rfile = uproot.open(os.path.join(basedir, rname))
    gof = rfile["Hgof"].to_hist()

    fig, ax = plt.subplots(figsize=(10,10))
    fig.subplots_adjust(hspace=0.0, wspace=0.0)
    
    histo = hist.Hist("Events", hist.Bin("x", "x", gof.axes.edges[0]))
    histo.fill(**{"x": np.zeros(0), "weight" : np.zeros(0)})
    histo.values()[()][:] = gof.values()
    
    pvalue = np.sum(histo.values()[()][np.where(histo.dense_axes()[0].edges() < 0)[0]])/np.sum(histo.values()[()])
    rebin_histo = histo.rebin("x", 50)
    hep.plot.histplot(rebin_histo.values()[()], rebin_histo.dense_axes()[0].edges(), ax=ax, color="k", histtype="step")
        # fill values up to zero
    hep.plot.histplot(rebin_histo.values()[()][np.where(rebin_histo.dense_axes()[0].edges() <= 0)[0][:-1]], rebin_histo.dense_axes()[0].edges()[np.where(rebin_histo.dense_axes()[0].edges() <= 0)[0]], ax=ax, color="b", alpha=0.5, histtype="fill")
    
    ax.axvline(x=0, color="r")
    ax.autoscale()
    ax.legend(title="Saturated model\n%s, $\ell$+jets\np-value: %.3f" % (label, pvalue), title_fontsize=30, loc="upper right")
    ax.set_xlabel("-2log($L_{toy}/L_{data}$)", loc="center")
    ax.set_ylabel("Number of pseudo-experiments")
    figname = os.path.join(outdir, f"GOF_{rname.replace('.root', '')}")
    fig.savefig(figname)
    print(f"{figname} written")
    plt.close(fig)

## background-only
rname = f"toys_bkg.root"
make_plot(rname, "background-only")

## best-fit signal plots
print(f"\nMaking plots for: A 750 5.0\n")
rname = f"toys_A_750_5p0.root"
make_plot(rname, "$A_{750\ GeV}^{5.0\%}$")

print(f"\nMaking plots for: H 725 3.0\n")
rname = f"toys_H_725_3p0.root"
make_plot(rname, "$H_{725\ GeV}^{3.0\%}$")

## standard signal points
for parity in ["A", "H"]:
    for mass in ["400", "800"]:
        for width in ["5.0"]:
            print(f"\nMaking plots for: {parity} {mass} {width}\n")
            
            rname = f"toys_{parity}_{mass}_{width.replace('.', 'p')}.root"
            make_plot(rname, "$%s_{%s\ GeV}^{%s\%%}$" % (parity, mass, width))


toc = time.time()
print("Total time: %.1f" % (toc - tic))
