# matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend('agg')
from matplotlib import rcParams
rcParams['font.size'] = 20
rcParams["savefig.format"] = 'png'
rcParams["savefig.bbox"] = 'tight'
from coffea.util import load
from pdb import set_trace
import os
import styles
import Utilities.plot_tools as plt_tools
import Utilities.prettyjson as prettyjson
#import re
from coffea import hist
import numpy as np
import fnmatch
import Plotter as Plotter
from coffea.hist import plot

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('year', choices=['2016', '2017', '2018'], help='What year is the ntuple from.')

args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'ttbar_post_alpha_reco'

input_dir = '/'.join([proj_dir, 'results', '%s_%s' % (args.year, jobid), analyzer])
f_ext = 'TOT.coffea'
outdir = '/'.join([proj_dir, 'plots', '%s_%s' % (args.year, jobid), analyzer])
if not os.path.isdir(outdir):
    os.makedirs(outdir)

fnames = sorted(['%s/%s' % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)])

#set_trace()
hdict = plt_tools.add_coffea_files(fnames) if len(fnames) > 1 else load(fnames[0])

blurb = 'tight $e/\mu$+3 jets\n$n_{btags} \geq$ 2'

alpha_corrections = {
    'E_All' : ('$\\alpha_{E}$ All', '#e41a1c'), ## red,
    'E_Mtt' : ('$\\alpha_{E}$ Mtt', '#377eb8'), ## blue,
    'P_All' : ('$\\alpha_{P}$ All', '#4daf4a'), ## green,
    'P_Mtt' : ('$\\alpha_{P}$ Mtt', '#ff7f00'), ## orange
    'Uncorrected' : ('Uncorrected', 'black'),
}


variables = {
    'Reco_mtt': ('m($t\\bar{t}$) [GeV]', 2, (200., 2000.)),
    'Reco_mthad': ('m($t_{h}$) [GeV]', 2, (0., 300.)),
    'Reco_thad_ctstar': ('cos($\\theta^{*}_{t_{h}}$)', 2, (-1., 1.)),
    'Reso_mtt': ('m($t\\bar{t}$) Resolution [GeV]', 1, (-300., 300.)),
    'Reso_mthad': ('m($t_{h}$) Resolution [GeV]', 2, (-200., 200.)),
    'Reso_thad_ctstar': ('cos($\\theta^{*}_{t_{h}}$) Resolution', 2, (-1., 1.)),
}


    ## get plotting colors/settings
hstyles = styles.styles
stack_fill_opts = {'alpha': 0.8, 'edgecolor':(0,0,0,.5)}
stack_error_opts = {'edgecolor':(0,0,0,.5)}

    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open('%s/inputs/lumis_data.json' % proj_dir).read())[args.year]
lumi_correction = load('%s/Corrections/%s/MC_LumiWeights_IgnoreSigEvts.coffea' % (proj_dir, jobid))[args.year]
lumi_to_use = (data_lumi_year['Muons']+data_lumi_year['Electrons'])/2000.

        # scale ttJets events, split by reconstruction type, by normal ttJets lumi correction
ttJets_permcats = ['*right', '*matchable', '*unmatchable', '*other']
names = [dataset for dataset in list(set([key[0] for key in hdict[list(variables.keys())[0]].values().keys()]))] # get dataset names in hists
ttJets_cats = [name for name in names if any([fnmatch.fnmatch(name, cat) for cat in ttJets_permcats])] # gets ttJets(_PS)_other, ...
if len(ttJets_cats) > 0:
    ttJets_lumi_topo = '_'.join(ttJets_cats[0].split('_')[:-1]) # gets ttJets or ttJets_PS
    mu_lumi = lumi_correction['Muons'][ttJets_lumi_topo]
    el_lumi = lumi_correction['Electrons'][ttJets_lumi_topo]
    lumi_correction['Muons'].update({key: mu_lumi for key in ttJets_cats})
    lumi_correction['Electrons'].update({key: el_lumi for key in ttJets_cats})
for hname in hdict.keys():
    if hname == 'cutflow': continue
    hdict[hname].scale(lumi_correction, axis='dataset')


## make groups based on process
process = hist.Cat("process", "Process", sorting='placement')
process_cat = "dataset"
process_groups = plt_tools.make_dataset_groups('Muon', args.year, samples=names)
for hname in hdict.keys():
    if hname == 'cutflow': continue
    hdict[hname] = hdict[hname].group(process_cat, process, process_groups)
    

    ## make plots
for hname in variables.keys():
    if hname not in hdict.keys():
        raise ValueError("%s not found in file" % hname)

    histo = hdict[hname]
    h_mu = histo[:, 'Muon', :].integrate('leptype')
    h_mu.scale(lumi_correction['Muons'], axis='process')
    h_el = histo[:, 'Electron', :].integrate('leptype')
    h_el.scale(lumi_correction['Electrons'], axis='process')
    h_tot = h_mu+h_el

    if h_tot.dense_dim() == 1:
        xtitle, rebinning, x_lims = variables[hname]
        if rebinning != 1:
            xaxis_name = h_tot.dense_axes()[0].name
            h_tot = h_tot.rebin(xaxis_name, rebinning)

        for cat in list(set([key[0] for key in h_tot.values().keys()])):
            pltdir = '/'.join([outdir, cat])
            if not os.path.isdir(pltdir):
                os.makedirs(pltdir)

            histo = h_tot[cat].integrate('process')

            fig, ax = plt.subplots()
            fig.subplots_adjust(hspace=.07)

            plot.plot1d(histo,
                overlay=histo.axes()[0].name,
                ax=ax,
                clear=False,
                line_opts={'linestyle' : '-'},
            )
            ax.autoscale(axis='x', tight=True)
            ax.set_ylim(0, None)
            ax.set_xlabel(xtitle)
            ax.set_xlim(x_lims)
        
            #set_trace()
                ## set legend and corresponding colors
            handles, labels = ax.get_legend_handles_labels()
            for idx, label in enumerate(labels):
                labels[idx] = alpha_corrections[label][0]
                handles[idx].set_color(alpha_corrections[label][1])
            # call ax.legend() with the new values
            ax.legend(handles,labels, loc='upper right')

                # add perm category 
            #set_trace()
            ax.text(
                0.02, 0.95, hstyles[cat]['name'].split(' ')[-1].capitalize(),#blurb,
                fontsize=rcParams['font.size']*0.75, 
                horizontalalignment='left', 
                verticalalignment='bottom', 
                transform=ax.transAxes
            )
            ax = hep.cms.cmslabel(ax=ax, data=False, paper=False, year=args.year, lumi=round(lumi_to_use, 1), fontsize=18)

            figname = '%s/%s' % (pltdir, '_'.join([cat, hname]))
            fig.savefig(figname)
            print('%s written' % figname)
            plt.close()


