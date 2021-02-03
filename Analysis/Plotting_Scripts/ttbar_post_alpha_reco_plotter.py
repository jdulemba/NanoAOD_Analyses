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
import Utilities.styles as styles
import Utilities.plot_tools as plt_tools
import Utilities.prettyjson as prettyjson
import re
from coffea import hist
import numpy as np
import fnmatch
import Utilities.Plotter as Plotter
from coffea.hist import plot

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('year', choices=['2016', '2017', '2018'], help='What year is the ntuple from.')
parser.add_argument('--plot', default='all', choices=['nosys', 'uncs', 'all'], help='Make plots for no systematics, variations of JES/JER systematics, or both.')
args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'ttbar_post_alpha_reco'

input_dir = os.path.join(proj_dir, 'results', '%s_%s' % (args.year, jobid), analyzer)
f_ext = 'TOT.coffea'
outdir = os.path.join(proj_dir, 'plots', '%s_%s' % (args.year, jobid), analyzer)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

fnames = sorted(['%s/%s' % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)])

#set_trace()
hdict = plt_tools.add_coffea_files(fnames) if len(fnames) > 1 else load(fnames[0])

jet_mults = {
    '3Jets' : '3 jets',
    '4Jets' : '4 jets',
    '5PJets' : '5+ jets',
}

blurb = 'tight $e/\mu$+3 jets\n$n_{btags} \geq$ 2'

alpha_corrections = {
    'E_All_1D' : ('$\\alpha_{E}$ All Lin.', '#377eb8', 2), ## blue
    'E_All_2D' : ('$\\alpha_{E}$ All Quad.', '#e42a2c', 2), ## red
    'E_Mtt' : ('$\\alpha_{E}$ Mtt', '#4daf4a', 2), ## green
    'P_All_1D' : ('$\\alpha_{P}$ All Lin.', '#984ea3', 2), ## purple
    'P_All_2D' : ('$\\alpha_{P}$ All Quad.', '#ff7f00', 2), ## orange
    'P_Mtt' : ('$\\alpha_{P}$ Mtt', '#a65628', 2), ## brown
    'Uncorrected' : ('Uncorrected', 'k', 2),
}

corr_to_use = 'E_All_2D'
comp_3j_mask = re.compile(r'((?:%s))' % '|'.join(['E_All_2D*', 'Uncorrected*']))

systematics = {
    'nosys' : ('Nominal', 'k', 2),
    'JES_UP' : ('JES Up', '#e42a2c', 2), ## red
    'JES_DW' : ('JES Down', '#377eb8', 2), ## blue
    'JER_UP' : ('JER Up', '#4daf4a', 2), ## green
    'JER_DW' : ('JER Down', '#984ea3', 2), ## purple
}

variables = {
    #'Reco_mtt': ('m($t\\bar{t}$) [GeV]', 2, (200., 2000.)),
    'Reco_mtt': ('m($t\\bar{t}$) [GeV]', 2, (200., 1000.)),
    'Reco_mthad': ('m($t_{h}$) [GeV]', 2, (0., 300.)),
    'Reco_thad_ctstar': ('cos($\\theta^{*}_{t_{h}}$)', 2, (-1., 1.)),
    'Reco_thad_ctstar_abs': ('|cos($\\theta^{*}_{t_{h}}$)|', 2, (0., 1.)),
    'Reso_mtt': ('m($t\\bar{t}$) Resolution [GeV]', 1, (-300., 300.)),
    'Reso_mthad': ('m($t_{h}$) Resolution [GeV]', 2, (-200., 200.)),
    'Reso_thad_ctstar': ('cos($\\theta^{*}_{t_{h}}$) Resolution', 2, (-1., 1.)),
    'Reso_thad_ctstar_abs': ('|cos($\\theta^{*}_{t_{h}}$)| Resolution', 2, (-1., 1.)),
}


    ## get plotting colors/settings
hstyles = styles.styles
stack_fill_opts = {'alpha': 0.8, 'edgecolor':(0,0,0,.5)}
stack_error_opts = {'edgecolor':(0,0,0,.5)}

    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, 'inputs', 'lumis_data.json')).read())[args.year]
lumi_correction = load(os.path.join(proj_dir, 'Corrections', jobid, 'MC_LumiWeights_allTTJets.coffea'))[args.year]
lumi_to_use = (data_lumi_year['Muons']+data_lumi_year['Electrons'])/2000.

        # scale ttJets events, split by reconstruction type, by normal ttJets lumi correction
ttJets_permcats = ['*right', '*matchable', '*unmatchable', '*other']
names = [dataset for dataset in list(set([key[0] for key in hdict[list(variables.keys())[0]].values().keys()]))] # get dataset names in hists

ttJets_cats = [name for name in names if any([fnmatch.fnmatch(name, cat) for cat in ttJets_permcats])] # gets ttJets(_PS)_other, ...
if len(ttJets_cats) > 0:
    for tt_cat in ttJets_cats:
        ttJets_lumi_topo = '_'.join(tt_cat.split('_')[:-1]) # gets ttJets[SL, Had, DiLep] or ttJets_PS
        mu_lumi = lumi_correction['Muons'][ttJets_lumi_topo]
        el_lumi = lumi_correction['Electrons'][ttJets_lumi_topo]
        lumi_correction['Muons'].update({tt_cat: mu_lumi})
        lumi_correction['Electrons'].update({tt_cat: el_lumi})


## make groups based on process
process = hist.Cat("process", "Process", sorting='placement')
process_cat = "dataset"
process_groups = plt_tools.make_dataset_groups('Muon', args.year, samples=names)
    

    ## make plots
for hname in variables.keys():
    if hname not in hdict.keys():
        raise ValueError("%s not found in file" % hname)

    histo = hdict[hname]
    h_mu = histo[:, :, :, 'Muon', :].integrate('leptype')
    h_mu.scale(lumi_correction['Muons'], axis='dataset')
    h_el = histo[:, :, :, 'Electron', :].integrate('leptype')
    h_el.scale(lumi_correction['Electrons'], axis='dataset')
    h_tot = h_mu+h_el
    h_tot = h_tot.group(process_cat, process, process_groups)
    #set_trace()    

    if (args.plot == 'nosys') or (args.plot == 'all'):
        xtitle, rebinning, x_lims = variables[hname]
        if rebinning != 1:
            xaxis_name = h_tot.dense_axes()[0].name
            h_tot = h_tot.rebin(xaxis_name, rebinning)
    
        nosys_histo = h_tot[:, 'nosys', :, :].integrate('sys')
            # make plot for each jet multiplicity
        #for jmult in ['3Jets']:
        for jmult in sorted(set([key[1] for key in nosys_histo.values().keys()])):
            for cat in sorted(set([key[0] for key in nosys_histo.values().keys()])):
                pltdir = os.path.join(outdir, jmult, cat)
                if not os.path.isdir(pltdir):
                    os.makedirs(pltdir)

                hslice = nosys_histo[cat, jmult].integrate('process').integrate('jmult')
    
                fig, ax = plt.subplots()
                fig.subplots_adjust(hspace=.07)
    
                plot.plot1d(hslice,
                    overlay=hslice.axes()[0].name,
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
                    handles[idx].set_linewidth(alpha_corrections[label][2])
                # call ax.legend() with the new values
                ax.legend(handles,labels, loc='upper right')
    
                    # add perm category 
                #set_trace()
                ax.text(
                    0.02, 0.85, 'tight $e/\mu$, %s\n%s' % (jet_mults[jmult], hstyles[cat]['name']),
                    fontsize=rcParams['font.size'], horizontalalignment='left', verticalalignment='bottom', 
                    transform=ax.transAxes
                )
                hep.cms.cmslabel(ax=ax, data=False, paper=False, year=args.year, lumi=round(lumi_to_use, 1), fontsize=rcParams['font.size'])
    
                figname = os.path.join(pltdir, '_'.join([jmult, cat, hname]))
                fig.savefig(figname)
                print('%s written' % figname)
                plt.close()
    
            # compare plots for each jet multiplicity 
        for cat in sorted(set([key[0] for key in nosys_histo.values().keys()])):
            pltdir = os.path.join(outdir, 'Comp_JMults', cat)
            if not os.path.isdir(pltdir):
                os.makedirs(pltdir)
   
                # yields 
            fig, ax = plt.subplots()
            fig.subplots_adjust(hspace=.07)
    
                # norm 
            fig_norm, ax_norm = plt.subplots()
            fig_norm.subplots_adjust(hspace=.07)
    
            for jmult in sorted(set([key[1] for key in nosys_histo.values().keys()])):
                hslice = nosys_histo[cat, jmult, :].integrate('process').integrate('jmult') if jmult == '3Jets' else nosys_histo[cat, jmult, :].integrate('process').integrate('corrtype')
                if jmult == '3Jets':
                    hslice = hslice[comp_3j_mask]
                #set_trace()

                    # yields 
                plot.plot1d(hslice,
                    overlay=hslice.axes()[0].name,
                    ax=ax,
                    clear=False,
                    line_opts={'linestyle' : '-'},
                )

                    # norm
                #set_trace()
                for corr in hslice.values().keys():
                    ax_norm = Plotter.plot_1D(values=hslice.values()[corr]/np.sum(hslice.values()[corr]), bins=hslice.dense_axes()[0].edges(),
                        ax=ax_norm, xlimits=x_lims, xlabel=xtitle, ylabel='A.U.', label=corr[0], histtype='step')

                ## set legend and corresponding colors
            handles, labels = ax.get_legend_handles_labels()
            for idx, label in enumerate(labels):
                if label == '4Jets':
                    labels[idx] = jet_mults[label]
                    handles[idx].set_color('b')
                    handles[idx].set_linewidth(2)
                elif label == '5PJets':
                    labels[idx] = jet_mults[label]
                    handles[idx].set_color('g')
                    handles[idx].set_linewidth(2)
                else:
                    labels[idx] = '3 jets, %s' % alpha_corrections[label][0]
                    handles[idx].set_color(alpha_corrections[label][1])
                    handles[idx].set_linewidth(alpha_corrections[label][2])

            # set axes and call ax.legend() with the new values
            ax.autoscale(axis='x', tight=True)
            ax.set_ylim(0, None)
            ax.set_xlabel(xtitle)
            ax.set_xlim(x_lims)
                
            ax.legend(handles,labels, loc='upper right')
    
                # add perm category 
            ax.text(
                0.02, 0.85, 'tight $e/\mu$\n%s' % hstyles[cat]['name'],
                fontsize=rcParams['font.size'], horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
            )
            hep.cms.cmslabel(ax=ax, data=False, paper=False, year=args.year, lumi=round(lumi_to_use, 1), fontsize=rcParams['font.size'])

            figname = os.path.join(pltdir, '_'.join(['Comp_JMults', cat, hname]))
            fig.savefig(figname)
            print('%s written' % figname)
            plt.close(fig)
    
                # for norm hists
            handles, labels = ax_norm.get_legend_handles_labels()
            for idx, label in enumerate(labels):
                if label == '4Jets':
                    labels[idx] = jet_mults[label]
                    handles[idx].set_color('b')
                    handles[idx].set_linewidth(2)
                elif label == '5PJets':
                    labels[idx] = jet_mults[label]
                    handles[idx].set_color('g')
                    handles[idx].set_linewidth(2)
                else:
                    labels[idx] = '3 jets, %s' % alpha_corrections[label][0]
                    handles[idx].set_color(alpha_corrections[label][1])
                    handles[idx].set_linewidth(alpha_corrections[label][2])

            # set axes and call ax.legend() with the new values
            ax_norm.autoscale(axis='x', tight=True)
            ax_norm.set_ylim(0, None)
            ax_norm.set_xlabel(xtitle)
            ax_norm.set_ylabel('A.U.')
            ax_norm.set_xlim(x_lims)
                
            ax_norm.legend(handles,labels, loc='upper right')
    
                # add perm category 
            ax_norm.text(
                0.02, 0.85, 'tight $e/\mu$\n%s\n' % hstyles[cat]['name'],
                fontsize=rcParams['font.size'], horizontalalignment='left', verticalalignment='bottom', transform=ax_norm.transAxes
            )
            hep.cms.cmslabel(ax=ax_norm, data=False, paper=False, year=args.year, lumi=round(lumi_to_use, 1), fontsize=rcParams['font.size'])
            #set_trace()
    
            figname_norm = os.path.join(pltdir, '_'.join(['Comp_JMults', cat, hname, 'Norm']))
            fig_norm.savefig(figname_norm)
            print('%s written' % figname_norm)
            plt.close(fig_norm)
    
    
    if (args.plot == 'uncs') or (args.plot == 'all'):
        xtitle, rebinning, x_lims = variables[hname]
        if rebinning != 1:
            xaxis_name = h_tot.dense_axes()[0].name
            h_tot = h_tot.rebin(xaxis_name, rebinning)
    
        uncs_histo = h_tot[:, :, '3Jets', :].integrate('jmult')
            # make plot for each category
        for cat in sorted(set([key[0] for key in uncs_histo.values().keys()])):
            pltdir = os.path.join(outdir, '3Jets', 'Sys_Vars', cat)
            if not os.path.isdir(pltdir):
                os.makedirs(pltdir)
            for corr in sorted(set([key[2] for key in uncs_histo.values().keys()])):

                hslice = uncs_histo[cat, :, corr].integrate('process').integrate('corrtype')

                fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                #fig, ax = plt.subplots()
                fig.subplots_adjust(hspace=.07)
    
                    # plot yields
                plot.plot1d(hslice,
                    overlay=hslice.axes()[0].name,
                    ax=ax,
                    clear=False,
                    line_opts={'linestyle' : '-'},
                )
                ax.autoscale(axis='x', tight=True)
                ax.set_ylim(0, None)
                ax.set_xlabel(None)
                ax.set_xlim(x_lims)

                    # plot ratios
                for sys in sorted([key[0] for key in hslice.values().keys()]):
                    if sys == 'nosys': continue
                    nom_histo = hslice['nosys'].integrate('sys')
                    sys_histo = hslice[sys].integrate('sys')

                    sys_first_valid_bin, sys_last_valid_bin = np.where(~np.isnan(sys_histo.values()[()]/nom_histo.values()[()]))[0][0], np.where(~np.isnan(sys_histo.values()[()]/nom_histo.values()[()]))[0][-1]+1
                    ratio_masked_vals = np.ma.masked_where(np.isnan((sys_histo.values()[()]/nom_histo.values()[()])[sys_first_valid_bin:sys_last_valid_bin]), (sys_histo.values()[()]/nom_histo.values()[()])[sys_first_valid_bin:sys_last_valid_bin])
                    rax.step(nom_histo.dense_axes()[0].edges()[sys_first_valid_bin:sys_last_valid_bin+1], np.r_[ratio_masked_vals, ratio_masked_vals[-1]], where='post', **{'linestyle' : '-', 'color' : systematics[sys][1], 'linewidth' :  systematics[sys][2]})

                rax.set_xlabel(xtitle)
                rax.set_ylabel('Sys/Nominal')
                rax.set_xlim(x_lims)
                rax.set_ylim(0.5, 1.5)
                rax.axhline(1, **{'linestyle': '--', 'color': (0, 0, 0, 0.5), 'linewidth': 1})
            
                    ## set legend and corresponding colors
                handles, labels = ax.get_legend_handles_labels()
                for idx, label in enumerate(labels):
                    labels[idx] = systematics[label][0]
                    handles[idx].set_color(systematics[label][1])
                    handles[idx].set_linewidth(systematics[label][2])
                # call ax.legend() with the new values
                ax.legend(handles,labels, loc='upper right', title=alpha_corrections[corr][0])
    
                    # add lepton/jet mult, and tt perm category 
                ax.text(
                    0.02, 0.87, 'tight $e/\mu$, %s\n%s' % (jet_mults['3Jets'], hstyles[cat]['name']),
                    fontsize=rcParams['font.size']*0.9, horizontalalignment='left', verticalalignment='bottom', 
                    transform=ax.transAxes
                )
                hep.cms.cmslabel(ax=ax, data=False, paper=False, year=args.year, lumi=round(lumi_to_use, 1), fontsize=rcParams['font.size'])
    
                #set_trace()
                figname = os.path.join(pltdir, '_'.join(['3Jets', cat, corr, hname, 'Sys_Comp']))
                fig.savefig(figname)
                print('%s written' % figname)
                plt.close()
    
