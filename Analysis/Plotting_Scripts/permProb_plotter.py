from coffea.hist import plot
# matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend('agg')
from matplotlib import rcParams
rcParams['font.size'] = 18
rcParams["savefig.format"] = 'png'
rcParams["savefig.bbox"] = 'tight'
from coffea.util import load, save
from pdb import set_trace
import os
import styles
import Utilities.plot_tools as plt_tools
import Utilities.prettyjson as prettyjson
from coffea import hist
from coffea.lookup_tools.dense_lookup import dense_lookup
import numpy as np
import Plotter as Plotter
from scipy import interpolate

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('--year', type=str, help='What year is the ntuple from.')
parser.add_argument('--njets', type=str, help='Choose jet multiplicity.')
args = parser.parse_args()

years_to_run = [args.year] if args.year else ['2016', '2017', '2018']
njets_to_run = ['3', '4+']
if args.njets == '3':
    njets_to_run = ['3']
elif args.njets == '4+':
    njets_to_run = ['4+']
else:
    njets_to_run = ['3', '4+']

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'permProbComputer'

permProbs = {
    '2016' : {
        '3Jets' : {},
        '4PJets' : {},
    },
    '2017' : {
        '3Jets' : {},
        '4PJets' : {},
    },
    '2018' : {
        '3Jets' : {},
        '4PJets' : {},
    },
}

jet_mults = {
    '3Jets' : '3 jets',
    '4PJets' : '4+ jets',
    #'4Jets' : '4 jets',
    #'5Jets' : '5 jets',
    #'6PJets' : '5+ jets'
}

lep_cats = {
    'LoT' : 'loose or tight $e/\mu$',
    'Tight' : 'tight $e/\mu$',
}

perm_cats = {
    '3Jets' : {
        'Correct' : ['Correct'],
        'Wrong' : ['Wrong'],
    },
    '4PJets' : {
        '1D' : {
            'Correct_BLep' : ['Correct', 'Right_TLep'],
            'Wrong_BLep' : ['Right_THad', 'Wrong']
        },
        '2D' : {
            'Correct_THad' : ['Correct', 'Right_THad'],
            'Wrong_THad' : ['Right_TLep', 'Wrong'],
        },
    },
}

variables_3jets = {
    'Lost_nusolver_chi2' : ('$\\chi_{\\nu}^{2}$', 5, (0., 1000.), True),
    'Lost_nusolver_dist' : ('$D_{\\nu, min}$ [GeV]', 1, (0., 150.), True),
    'Lost_mbpjet' : ('m(b+j) [GeV]', 1, (0., 500.), True),
    'Merged_nusolver_chi2' : ('$\\chi_{\\nu}^{2}$', 5, (0., 1000.), True),
    'Merged_nusolver_dist' : ('$D_{\\nu, min}$ [GeV]', 2, (0., 150.), True),
    'Merged_mbpjet_vs_maxmjet' : ('max m(jet) [GeV]', 'm(b+j) [GeV]', 2, (0., 150.), 2, (0., 500.), True),
}

variables_4pjets = {
    'nusolver_chi2' : ('$\\chi_{\\nu}^{2}$', 5, (0., 1000.), True),
    'nusolver_dist' : ('$D_{\\nu, min}$ [GeV]', 1, (0., 150.), True),
    'mWHad_vs_mTHad' : ('m($t_{h}$) [GeV]', 'm($W_{h}$) [GeV]', 8, (0., 500.), 8, (0., 500.), True),
}

    ## get plotting colors/settings
hstyles = styles.styles

    ## get data lumi and scale MC by lumi
data_lumi_dict = prettyjson.loads(open('%s/inputs/lumis_data.json' % proj_dir).read())
lumi_correction = load('%s/Corrections/%s/MC_LumiWeights_IgnoreSigEvts.coffea' % (proj_dir, jobid))

    ## make groups based on perm category
pcat = hist.Cat("permcat", "Perm Category", sorting='placement')
pcat_cat = "permcat"

    ## make groups based on process
process = hist.Cat("process", "Process", sorting='placement')
process_cat = "dataset"

for year in years_to_run:
#for year in ['2016']:#, '2017', '2018']:
    input_dir = '/'.join([proj_dir, 'results', '%s_%s' % (year, jobid), analyzer])
    f_ext = 'TOT.coffea'
    outdir = '/'.join([proj_dir, 'plots', '%s_%s' % (year, jobid), analyzer])
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    
    fnames = ['%s/%s' % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)]
    fnames = sorted(fnames)
    
    hdict = plt_tools.add_coffea_files(fnames) if len(fnames) > 1 else load(fnames[0])
    
    lumi_to_use = (data_lumi_dict[year]['Muons']+data_lumi_dict[year]['Electrons'])/2000.

        ## make groups based on process
    process_groups = plt_tools.make_dataset_groups('Muon', year) # works when only MC present
    for hname in hdict.keys():
        if hname == 'cutflow': continue
        hdict[hname] = hdict[hname].group(process_cat, process, process_groups)


    if '3' in njets_to_run:
            ## make 3 jets plots
        for hname in variables_3jets.keys():
            if not hname in hdict.keys():
                raise ValueError("Hist %s not found" % hname)
        
            topo = hname.split('_')[0] # Merged/Lost
    
            histo = hdict[hname]
                ## rescale hist by lumi for muons and electrons separately and then combine
            h_mu = histo[:, :, 'Muon', :, :, :].integrate('leptype')
            h_mu.scale(lumi_correction[year]['Muons'], axis='process')
            h_el = histo[:, :, 'Electron', :, :, :].integrate('leptype')
            h_el.scale(lumi_correction[year]['Electrons'], axis='process')
            h_tot = h_mu+h_el
            h_tot = h_tot[:, :, :, 'MTHigh', :].integrate('process').integrate('mtregion')
        
            if h_tot.dense_dim() == 1:
                xtitle, rebinning, x_lims, save_dist = variables_3jets[hname]
        
                ## hists should have 3 category axes (jet multiplicity, lepton type, perm category) followed by variable
                for jmult in list(set([key[0] for key in h_tot.values().keys()])):
                    for lepcat in list(set([key[1] for key in h_tot.values().keys()])):
                        pltdir = '/'.join([outdir, jmult, lepcat, topo])
                        if not os.path.isdir(pltdir):
                            os.makedirs(pltdir)
        
                        fig, ax = plt.subplots()
                        fig.subplots_adjust(hspace=.07)
                        hslice = h_tot[jmult, lepcat].integrate('jmult').integrate('lepcat')
        
                            ## make groups based on perm category
                        pcat_groups = perm_cats[jmult]
                        hslice = hslice.group(pcat_cat, pcat, pcat_groups)
        
                        if rebinning != 1:
                            xaxis_name = hslice.dense_axes()[0].name
                            hslice = hslice.rebin(xaxis_name, rebinning)
    
                            # save distribution to dict
                        if save_dist and (lepcat == 'Tight'):
                            hcor = hslice['Correct'].integrate('permcat')
                            edges = hcor.dense_axes()[0].edges()
                            lookup = dense_lookup(*(hcor.values().values()), edges) # not normalized
                            permProbs[year][jmult].update({hname : lookup})
    
                            ## plot comparison of perm categories
                        plot.plot1d(hslice,
                            overlay=hslice.axes()[0].name,
                            ax=ax,
                            clear=False,
                            line_opts={'linestyle' : '-'},
                            density=True, # normalized to 1,
                        )
                        ax.autoscale(axis='x', tight=True)
                        ax.set_ylim(0, None)
                        ax.set_ylabel(None)
                        ax.set_xlabel(None)
                        ax.set_xlim(x_lims)
        
                           ## set legend and corresponding colors
                        h_opts = {key: hstyles['%s_BLep' % key if 'nusolver' in hname else key] for key in pcat_groups.keys()}
                        handles, labels = ax.get_legend_handles_labels()
                        for idx, cat in enumerate(labels):
                            labels[idx] = h_opts[cat]['name']
                            handles[idx].set_color(h_opts[cat]['color'])
                            handles[idx].set_linestyle(h_opts[cat]['linestyle'])
                        # call ax.legend() with the new values
                        ax.legend(handles,labels, loc='upper right', title='%s Jet' % topo if topo == 'Lost' else '%s Jets' % topo)
        
                            # add lep category
                        ax.text(
                            0.02, 0.91, "%s\n%s" % (lep_cats[lepcat], jet_mults[jmult]),
                            fontsize=18, 
                            horizontalalignment='left', 
                            verticalalignment='bottom', 
                            transform=ax.transAxes
                        )
                            ## set axes labels and titles
                        plt.xlabel(xtitle)
                        ax = hep.cms.cmslabel(ax=ax, data=False, paper=False, year=year, lumi=round(lumi_to_use, 1))
        
                        #set_trace()
                        figname = '%s/%s' % (pltdir, '_'.join([jmult, lepcat, hname]))
                        fig.savefig(figname)
                        print('%s written' % figname)
                        plt.close()
    
            elif h_tot.dense_dim() == 2:
                xtitle, ytitle, x_rebinning, x_lims, y_rebinning, y_lims, save_dist = variables_3jets[hname]
                topo, yvar, useless, xvar = hname.split('_') # Merged/Lost
        
                #set_trace()
                ## hists should have 3 category axes (jet multiplicity, lepton category, permutation category) followed by variable(s)
                for jmult in list(set([key[0] for key in h_tot.values().keys()])):
                    for lepcat in list(set([key[1] for key in h_tot.values().keys()])):
                        pltdir = '/'.join([outdir, jmult, lepcat, topo])
                        if not os.path.isdir(pltdir):
                            os.makedirs(pltdir)
        
                        hslice = h_tot[jmult, lepcat].integrate('jmult').integrate('lepcat')
        
                            ## make groups based on perm category
                        pcat_groups = perm_cats[jmult]
                        hslice = hslice.group(pcat_cat, pcat, pcat_groups)
        
                        if x_rebinning != 1:
                            xaxis_name = hslice.dense_axes()[0].name
                            hslice = hslice.rebin(xaxis_name, x_rebinning)
                        if y_rebinning != 1:
                            yaxis_name = hslice.dense_axes()[1].name
                            hslice = hslice.rebin(yaxis_name, y_rebinning)
        
                            # save distribution to dict
                        if save_dist and (lepcat == 'Tight'):
                            hcor = hslice['Correct'].integrate('permcat')
                            edges = (hcor.dense_axes()[0].edges(), hcor.dense_axes()[1].edges())
                            lookup = dense_lookup(*(hcor.values().values()), edges) # not normalized
                            permProbs[year][jmult].update({hname : lookup})
    
    
                        #set_trace()
                            # make 1D projection along dense axes
                        for dax in range(2):
                            fig, ax = plt.subplots()
                            fig.subplots_adjust(hspace=.07)
        
                            #set_trace()
                            hproj = hslice.integrate(hslice.dense_axes()[dax].name)
                                ## plot comparison of perm categories
                            plot.plot1d(hproj,
                                overlay=hproj.axes()[0].name,
                                ax=ax,
                                clear=False,
                                line_opts={'linestyle' : '-'},
                                density=True, # normalized to 1,
                            )
                            ax.autoscale(axis='x', tight=True)
                            ax.set_ylim(0, None)
                            ax.set_ylabel(None)
                            ax.set_xlabel(None)
                            ax.set_xlim(y_lims) if dax == 0 else ax.set_xlim(x_lims)
        
                            #set_trace()
                               ## set legend and corresponding colors
                            h_opts ={key: hstyles[key] for key in pcat_groups.keys()}
                            handles, labels = ax.get_legend_handles_labels()
                            for idx, cat in enumerate(labels):
                                labels[idx] = h_opts[cat]['name']
                                handles[idx].set_color(h_opts[cat]['color'])
                                handles[idx].set_linestyle(h_opts[cat]['linestyle'])
                            # call ax.legend() with the new values
                            ax.legend(handles,labels, loc='upper right', title='%s Jet' % topo if topo == 'Lost' else '%s Jets' % topo)
        
                            #set_trace()
                                # add lep category
                            ax.text(
                                0.02, 0.91, "%s\n%s" % (lep_cats[lepcat], jet_mults[jmult]),
                                fontsize=18, 
                                horizontalalignment='left', 
                                verticalalignment='bottom', 
                                transform=ax.transAxes
                            )
        
                                ## set axes labels and titles
                            plt.xlabel(ytitle) if dax == 0 else plt.xlabel(xtitle)
                            ax = hep.cms.cmslabel(ax=ax, fontsize=18, data=False, paper=False, year=year, lumi=round(lumi_to_use, 1))
       
                            #set_trace() 
                            #figname = 'test.png'
                            figname = '%s/%s' % (pltdir, '_'.join([jmult, lepcat, topo, yvar if dax == 0 else xvar]))
                            fig.savefig(figname)
                            print('%s written' % figname)
                            plt.close()
                            #set_trace()
        
                            # make 2D plots for each permutation category
                        for cat in list(set([key[0] for key in hslice.values().keys()])):            
                                # raw hist
                            fig, ax = plt.subplots()
                            fig.subplots_adjust(hspace=.07)
        
                            hcat = hslice[cat].integrate('permcat')
                                ## plot 
                            plot.plot2d(hcat,
                                xaxis=hcat.axes()[0].name,
                                ax=ax,
                                clear=True,
                            )
                            ax.autoscale(axis='x', tight=True)
                            ax.set_ylim(y_lims)
                            ax.set_xlim(x_lims)
        
                            #set_trace()
                                # add lep category
                            ax.text(
                                0.02, 0.91, "%s\n%s" % (lep_cats[lepcat], jet_mults[jmult]),
                                fontsize=18, 
                                horizontalalignment='left', 
                                verticalalignment='bottom', 
                                transform=ax.transAxes
                            )
                                # add perm category
                            ax.text(
                                1., 0.95, cat,#hstyles[cat]['name'],
                                fontsize=18, 
                                horizontalalignment='right', 
                                verticalalignment='bottom', 
                                transform=ax.transAxes
                            )
                                ## set axes labels and titles
                            plt.xlabel(xtitle)
                            plt.ylabel(ytitle)
                            ax = hep.cms.cmslabel(ax=ax, data=False, paper=False, year=year, lumi=round(lumi_to_use, 1), fontsize=18)
       
                            #set_trace() 
                            figname = '%s/%s' % (pltdir, '_'.join([jmult, lepcat, hname, cat]))
                            fig.savefig(figname)
                            print('%s written' % figname)
                            plt.close()
    
                                ## normalized plots
                            fig, ax = plt.subplots()
                            fig.subplots_adjust(hspace=.07)
        
                            norm_values = hcat.values()[()]/np.sum(hcat.values()[()]) ## normalized array of hist values
                            opts = {'cmap_label' : '$P_{M}$'}
                            ax = Plotter.plot_2d_norm(hcat, xaxis_name=hcat.axes()[0].name, yaxis_name=hcat.axes()[1].name,
                                values=np.ma.masked_where(norm_values <= 0.0, norm_values), # mask nonzero probabilities for plotting
                                xlimits=x_lims, ylimits=y_lims, xlabel=xtitle, ylabel=ytitle,
                                ax=ax, **opts)

                                # add lep category
                            ax.text(
                                0.98, 0.02, "%s\n%s" % (lep_cats[lepcat], jet_mults[jmult]),
                                #0.02, 0.91, "%s\n%s" % (lep_cats[lepcat], jet_mults[jmult]),
                                fontsize=18, 
                                #horizontalalignment='left', 
                                horizontalalignment='right', 
                                verticalalignment='bottom', 
                                transform=ax.transAxes
                            )
                                # add perm category
                            ax.text(
                                1, 0.95, cat,
                                fontsize=18, 
                                horizontalalignment='right', 
                                verticalalignment='bottom', 
                                transform=ax.transAxes
                            )
                            ax = hep.cms.cmslabel(ax=ax, data=False, paper=False, year=year, lumi=round(lumi_to_use, 1), fontsize=18)
        
                            figname = '%s/%s_norm' % (pltdir, '_'.join([jmult, lepcat, hname, cat]))
                            fig.savefig(figname)
                            print('%s written' % figname)
                            plt.close()

                            ## plot distribution of mass disc values
                        fig, ax = plt.subplots()
                        fig.subplots_adjust(hspace=.07)
                        correct_nl_norm_values = -1*np.log(hslice['Correct'].integrate('permcat').values()[()]/np.sum(hslice['Correct'].integrate('permcat').values()[()]))
                        wrong_nl_norm_values = -1*np.log(hslice['Wrong'].integrate('permcat').values()[()]/np.sum(hslice['Wrong'].integrate('permcat').values()[()]))
                        correct_mass_disc_freq = correct_nl_norm_values.ravel()[np.isfinite(correct_nl_norm_values.ravel())] # get non-infinite mass disc values
                        correct_weight = correct_mass_disc_freq/np.sum(correct_mass_disc_freq)
                        wrong_mass_disc_freq = wrong_nl_norm_values.ravel()[np.isfinite(wrong_nl_norm_values.ravel())] # get non-infinite mass disc values
                        wrong_weight = wrong_mass_disc_freq/np.sum(wrong_mass_disc_freq)
                        xlim = (0, 20)
                        bins = np.linspace(xlim[0], xlim[1], 101)
                        ax.hist(correct_mass_disc_freq, bins, weights=correct_weight, linestyle='-', color='r', label='Correct', histtype='step')
                        ax.hist(wrong_mass_disc_freq, bins, weights=wrong_weight, linestyle='-', color='b', label='Wrong', histtype='step')

                        # call ax.legend() with the new values
                        handles, labels = ax.get_legend_handles_labels()
                        new_handles = [Line2D([], [], c=h.get_edgecolor()) for h in handles] ## needed to make legend lines instead of empty boxes
                        ax.legend(new_handles,labels, loc='upper right')

                        ax.autoscale(axis='x', tight=True)
                        ax.set_xlim(xlim)
                        ax.set_xlabel('$\\lambda_{M}$ = -log($P_{M}$)')

                            # add lep category
                        ax.text(
                            0.02, 0.91, "%s\n%s" % (lep_cats[lepcat], jet_mults[jmult]),
                            fontsize=18, 
                            horizontalalignment='left', 
                            verticalalignment='bottom', 
                            transform=ax.transAxes
                        )
                        ax = hep.cms.cmslabel(ax=ax, data=False, paper=False, year=year, lumi=round(lumi_to_use, 1))
        
                        figname = '%s/%s_massdisc' % (pltdir, '_'.join([jmult, lepcat, hname]))
                        fig.savefig(figname)
                        print('%s written' % figname)
                        #set_trace()
                        plt.close()


    
        

    if '4+' in njets_to_run:    
            ## make 4+ jets plots
        for hname in variables_4pjets.keys():
            if not hname in hdict.keys():
                raise ValueError("Hist %s not found" % hname)
        
            histo = hdict[hname]
                ## rescale hist by lumi for muons and electrons separately and then combine
            h_mu = histo[:, :, 'Muon', :, :, :].integrate('leptype')
            h_mu.scale(lumi_correction[year]['Muons'], axis='process')
            h_el = histo[:, :, 'Electron', :, :, :].integrate('leptype')
            h_el.scale(lumi_correction[year]['Electrons'], axis='process')
            h_tot = h_mu+h_el
            h_tot = h_tot[:, :, :, 'MTHigh', :].integrate('process').integrate('mtregion')
        
            if h_tot.dense_dim() == 1:
                xtitle, rebinning, x_lims, save_dist = variables_4pjets[hname]
        
                ## hists should have 3 category axes (jet multiplicity, lepton type, perm category) followed by variable
                for jmult in list(set([key[0] for key in h_tot.values().keys()])):
                    for lepcat in list(set([key[1] for key in h_tot.values().keys()])):
                        pltdir = '/'.join([outdir, jmult, lepcat])
                        if not os.path.isdir(pltdir):
                            os.makedirs(pltdir)
        
                        fig, ax = plt.subplots()
                        fig.subplots_adjust(hspace=.07)
                        hslice = h_tot[jmult, lepcat, :].integrate('jmult').integrate('lepcat')
        
                            ## make groups based on perm category
                        pcat_groups = perm_cats[jmult]['1D']
                        hslice = hslice.group(pcat_cat, pcat, pcat_groups)
        
                        if rebinning != 1:
                            xaxis_name = hslice.dense_axes()[0].name
                            hslice = hslice.rebin(xaxis_name, rebinning)
        
                            # save distribution to dict
                        if save_dist and (lepcat == 'Tight'):
                            hcor = hslice['Correct_BLep'].integrate('permcat')
                            edges = hcor.dense_axes()[0].edges()
                            lookup = dense_lookup(*(hcor.values().values()), edges) # not normalized
                            permProbs[year][jmult].update({hname : lookup})
    
                            ## plot comparison of perm categories
                        plot.plot1d(hslice,
                            overlay=hslice.axes()[0].name,
                            ax=ax,
                            clear=False,
                            line_opts={'linestyle' : '-'},
                            density=True, # normalized to 1,
                        )
                        ax.autoscale(axis='x', tight=True)
                        ax.set_ylim(0, None)
                        ax.set_ylabel(None)
                        ax.set_xlabel(None)
                        ax.set_xlim(x_lims)
        
                        #set_trace()
                           ## set legend and corresponding colors
                        h_opts ={key: hstyles[key] for key in pcat_groups.keys()}
                        handles, labels = ax.get_legend_handles_labels()
                        for idx, cat in enumerate(labels):
                            labels[idx] = h_opts[cat]['name']
                            handles[idx].set_color(h_opts[cat]['color'])
                            handles[idx].set_linestyle(h_opts[cat]['linestyle'])
                        # call ax.legend() with the new values
                        ax.legend(handles,labels, loc='upper right')
                        #set_trace()
        
                            # add lep category
                        ax.text(
                            0.02, 0.91, "%s\n%s" % (lep_cats[lepcat], jet_mults[jmult]),
                            fontsize=18, 
                            horizontalalignment='left', 
                            verticalalignment='bottom', 
                            transform=ax.transAxes
                        )
                            ## set axes labels and titles
                        plt.xlabel(xtitle)
                        ax = hep.cms.cmslabel(ax=ax, fontsize=18, data=False, paper=False, year=year, lumi=round(lumi_to_use, 1))
        
                        #set_trace()
                        figname = '%s/%s' % (pltdir, '_'.join([jmult, lepcat, hname]))
                        fig.savefig(figname)
                        print('%s written' % figname)
                        plt.close()
        
            elif h_tot.dense_dim() == 2:
                xtitle, ytitle, x_rebinning, x_lims, y_rebinning, y_lims, save_dist = variables_4pjets[hname]

                #set_trace()
                ## hists should have 3 category axes (jet multiplicity, lepton category, permutation category) followed by variable(s)
                for jmult in list(set([key[0] for key in h_tot.values().keys()])):
                    for lepcat in list(set([key[1] for key in h_tot.values().keys()])):
                        pltdir = '/'.join([outdir, jmult, lepcat])
                        if not os.path.isdir(pltdir):
                            os.makedirs(pltdir)
        
                        hslice = h_tot[jmult, lepcat, :].integrate('jmult').integrate('lepcat')
        
                            ## make groups based on perm category
                        pcat_groups = perm_cats[jmult]['2D']
                        hslice = hslice.group(pcat_cat, pcat, pcat_groups)
        
                        if x_rebinning != 1:
                            xaxis_name = hslice.dense_axes()[0].name
                            hslice = hslice.rebin(xaxis_name, x_rebinning)
                        if y_rebinning != 1:
                            yaxis_name = hslice.dense_axes()[1].name
                            hslice = hslice.rebin(yaxis_name, y_rebinning)
        
                            # make 1D projection along dense axes
                        for dax in range(2):
                            fig, ax = plt.subplots()
                            fig.subplots_adjust(hspace=.07)
        
                            #set_trace()
                            hproj = hslice.integrate(hslice.dense_axes()[dax].name)
                                ## plot comparison of perm categories
                            plot.plot1d(hproj,
                                overlay=hproj.axes()[0].name,
                                ax=ax,
                                clear=False,
                                line_opts={'linestyle' : '-'},
                                density=True, # normalized to 1,
                            )
                            ax.autoscale(axis='x', tight=True)
                            ax.set_ylim(0, None)
                            ax.set_ylabel(None)
                            ax.set_xlabel(None)
                            ax.set_xlim(y_lims) if dax == 0 else ax.set_xlim(x_lims)
        
                            #set_trace()
                               ## set legend and corresponding colors
                            h_opts ={key: hstyles[key] for key in pcat_groups.keys()}
                            handles, labels = ax.get_legend_handles_labels()
                            for idx, cat in enumerate(labels):
                                labels[idx] = h_opts[cat]['name']
                                handles[idx].set_color(h_opts[cat]['color'])
                                handles[idx].set_linestyle(h_opts[cat]['linestyle'])
                            # call ax.legend() with the new values
                            ax.legend(handles,labels, loc='upper right')
        
                            #set_trace()
                                # add lep category
                            ax.text(
                                0.02, 0.91, "%s\n%s" % (lep_cats[lepcat], jet_mults[jmult]),
                                fontsize=18, 
                                horizontalalignment='left', 
                                verticalalignment='bottom', 
                                transform=ax.transAxes
                            )
        
                                ## set axes labels and titles
                            plt.xlabel(ytitle) if dax == 0 else plt.xlabel(xtitle)
                            ax = hep.cms.cmslabel(ax=ax, fontsize=18, data=False, paper=False, year=year, lumi=round(lumi_to_use, 1))
        
                            #figname = 'test.png'
                            figname = '%s/%s' % (pltdir, '_'.join([jmult, lepcat, hname.split('_vs_')[1] if dax == 0 else hname.split('_vs_')[0]]))
                            fig.savefig(figname)
                            print('%s written' % figname)
                            plt.close()
                            #set_trace()
        
                            # make 2D plots for each permutation category
                        for cat in list(set([key[0] for key in hslice.values().keys()])):            
                            hcat = hslice[cat].integrate('permcat')

                                ## plot raw values
                            fig, ax = plt.subplots()
                            fig.subplots_adjust(hspace=.07)
        
                            plot.plot2d(hcat,
                                xaxis=hcat.axes()[0].name,
                                ax=ax,
                                clear=True,
                            )
                            ax.autoscale(axis='x', tight=True)
                            ax.set_ylim(y_lims)
                            ax.set_xlim(x_lims)
        
                            #set_trace()
                                # add lep category
                            ax.text(
                                0.02, 0.91, "%s\n%s" % (lep_cats[lepcat], jet_mults[jmult]),
                                #fontsize=rcParams['font.size']*0.75, 
                                fontsize=18, 
                                horizontalalignment='left', 
                                verticalalignment='bottom', 
                                transform=ax.transAxes
                            )
                                # add perm category
                            ax.text(
                                0.98, 0.95, hstyles[cat]['name'],
                                #fontsize=rcParams['font.size']*0.75, 
                                fontsize=18, 
                                horizontalalignment='right', 
                                verticalalignment='bottom', 
                                transform=ax.transAxes
                            )
                                ## set axes labels and titles
                            plt.xlabel(xtitle)
                            plt.ylabel(ytitle)
                            ax = hep.cms.cmslabel(ax=ax, data=False, paper=False, year=year, lumi=round(lumi_to_use, 1), fontsize=18)
        
                            figname = '%s/%s' % (pltdir, '_'.join([jmult, lepcat, hname, cat]))
                            fig.savefig(figname)
                            print('%s written' % figname)
                            plt.close()
                            #set_trace()

                                ## normalized plots before interpolation
                            fig, ax = plt.subplots()
                            fig.subplots_adjust(hspace=.07)
        
                            norm_values = hcat.values()[()]/np.sum(hcat.values()[()]) ## normalized array of hist values
                            opts = {'cmap_label' : '$P_{M}$'}
                            ax = Plotter.plot_2d_norm(hcat, xaxis_name=hcat.axes()[0].name, yaxis_name=hcat.axes()[1].name,
                                values=np.ma.masked_where(norm_values <= 0.0, norm_values), # mask nonzero probabilities for plotting
                                xlimits=x_lims, ylimits=y_lims, xlabel=xtitle, ylabel=ytitle,
                                ax=ax, **opts)

                                # add lep category
                            ax.text(
                                0.02, 0.91, "%s\n%s" % (lep_cats[lepcat], jet_mults[jmult]),
                                fontsize=18, 
                                horizontalalignment='left', 
                                verticalalignment='bottom', 
                                transform=ax.transAxes
                            )
                                # add perm category
                            ax.text(
                                0.95, 0.95, hstyles[cat]['name'],
                                fontsize=18, 
                                horizontalalignment='right', 
                                verticalalignment='bottom', 
                                transform=ax.transAxes
                            )
                            ax = hep.cms.cmslabel(ax=ax, data=False, paper=False, year=year, lumi=round(lumi_to_use, 1), fontsize=18)
        
                            figname = '%s/%s_norm_orig' % (pltdir, '_'.join([jmult, lepcat, hname, cat]))
                            fig.savefig(figname)
                            print('%s written' % figname)
                            plt.close()


                                ## normalized plots after interpolation
                            fig, ax = plt.subplots()
                            fig.subplots_adjust(hspace=.07)
        
                            #set_trace()
                            xcenter = hcat.axes()[0].centers()
                            ycenter = hcat.axes()[1].centers()
                            f = interpolate.interp2d(xcenter, ycenter, hcat.values()[()])
                            new_xcenters = np.arange(0, 505, 5)
                            new_ycenters = np.arange(0, 505, 5)
                            interp_vals = f(new_xcenters, new_ycenters)
                            norm_intval = interp_vals/np.sum(interp_vals) ## normalized array of values after interpolating
                                # save distribution to dict
                            if save_dist and (lepcat == 'Tight') and (cat == 'Correct_THad'):
                                lookup = dense_lookup(*(interp_vals, (new_xcenters, new_ycenters))) # not normalized
                                permProbs[year][jmult].update({hname : lookup})
    
                            values = np.ma.masked_where(norm_intval <= 0., norm_intval)
                            ax = Plotter.plot_2d_norm(hcat, xbins=new_xcenters, ybins=new_ycenters,
                                values=np.ma.masked_where(values <= 0.0, values), # mask nonzero probabilities for plotting
                                xlimits=x_lims, ylimits=y_lims, xlabel=xtitle, ylabel=ytitle,
                                ax=ax, **opts)

                                # add lep category
                            ax.text(
                                0.02, 0.91, "%s\n%s" % (lep_cats[lepcat], jet_mults[jmult]),
                                fontsize=18, 
                                horizontalalignment='left', 
                                verticalalignment='bottom', 
                                transform=ax.transAxes
                            )
                                # add perm category
                            ax.text(
                                0.95, 0.95, hstyles[cat]['name'],
                                fontsize=18, 
                                horizontalalignment='right', 
                                verticalalignment='bottom', 
                                transform=ax.transAxes
                            )
                            ax = hep.cms.cmslabel(ax=ax, data=False, paper=False, year=year, lumi=round(lumi_to_use, 1), fontsize=18)
        
                            figname = '%s/%s_norm_interp' % (pltdir, '_'.join([jmult, lepcat, hname, cat]))
                            fig.savefig(figname)
                            print('%s written' % figname)
                            plt.close()



                            ## plot distribution of mass disc values
                        fig, ax = plt.subplots()
                        fig.subplots_adjust(hspace=.07)
                        correct_nl_norm_values = -1*np.log(hslice['Correct_THad'].integrate('permcat').values()[()]/np.sum(hslice['Correct_THad'].integrate('permcat').values()[()]))
                        wrong_nl_norm_values = -1*np.log(hslice['Wrong_THad'].integrate('permcat').values()[()]/np.sum(hslice['Wrong_THad'].integrate('permcat').values()[()]))
                        correct_mass_disc_freq = correct_nl_norm_values.ravel()[np.isfinite(correct_nl_norm_values.ravel())] # get non-infinite mass disc values
                        correct_weight = correct_mass_disc_freq/np.sum(correct_mass_disc_freq)
                        wrong_mass_disc_freq = wrong_nl_norm_values.ravel()[np.isfinite(wrong_nl_norm_values.ravel())] # get non-infinite mass disc values
                        wrong_weight = wrong_mass_disc_freq/np.sum(wrong_mass_disc_freq)
                        bins = np.linspace(0, 20, 101)
                        ax.hist(correct_mass_disc_freq, bins, weights=correct_weight, linestyle='-', color='r', label='Correct $t_{h}$', histtype='step')
                        ax.hist(wrong_mass_disc_freq, bins, weights=wrong_weight, linestyle='-', color='b', label='Wrong $t_{h}$', histtype='step')

                        # call ax.legend() with the new values
                        handles, labels = ax.get_legend_handles_labels()
                        new_handles = [Line2D([], [], c=h.get_edgecolor()) for h in handles] ## needed to make legend lines instead of empty boxes
                        ax.legend(new_handles,labels, loc='upper right')

                        ax.autoscale(axis='x', tight=True)
                        ax.set_xlim(0, 20)
                        ax.set_xlabel('$\\lambda_{M}$ = -log($P_{M}$)')

                            # add lep category
                        ax.text(
                            0.02, 0.91, "%s\n%s" % (lep_cats[lepcat], jet_mults[jmult]),
                            fontsize=18, 
                            horizontalalignment='left', 
                            verticalalignment='bottom', 
                            transform=ax.transAxes
                        )
                        ax = hep.cms.cmslabel(ax=ax, data=False, paper=False, year=year, lumi=round(lumi_to_use, 1))
        
                        figname = '%s/%s_massdisc' % (pltdir, '_'.join([jmult, lepcat, hname]))
                        fig.savefig(figname)
                        print('%s written' % figname)
                        #set_trace()
                        plt.close()


    # write corrections to coffea file
if (len(njets_to_run) == 2) and (len(years_to_run) == 3):
    corrdir = '%s/Corrections/%s' % (proj_dir, jobid)
    if not os.path.isdir(corrdir):
        os.makedirs(corrdir)
    
    permprobs_name = '%s/prob_htt_%s.coffea' % (corrdir, jobid)
    save(permProbs, permprobs_name)
    print('\n', permprobs_name, 'written')
