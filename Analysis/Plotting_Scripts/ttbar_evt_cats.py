from coffea.hist import plot
# matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend('agg')
from matplotlib import rcParams
rcParams['font.size'] = 20
rcParams["savefig.format"] = 'png'
rcParams["savefig.bbox"] = 'tight'
from coffea.util import load, save
from pdb import set_trace
import os
import Utilities.plot_tools as plt_tools
import Utilities.prettyjson as prettyjson
from coffea import hist
import numpy as np
import Utilities.Plotter as Plotter

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
base_jobid = os.environ['base_jobid']
analyzer = 'ttbar_evt_cats'

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('year', choices=['2016', '2017', '2018'] if base_jobid == 'NanoAODv6' else ['2016APV', '2016', '2017', '2018'], help='What year is the ntuple from.')
args = parser.parse_args()

input_dir = os.path.join(proj_dir, 'results', '%s_%s' % (args.year, jobid), analyzer)
f_ext = 'TOT.coffea'
outdir = os.path.join(proj_dir, 'plots', '%s_%s' % (args.year, jobid), analyzer)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

fnames = sorted([os.path.join(input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)])
hdict = plt_tools.add_coffea_files(fnames) if len(fnames) > 1 else load(fnames[0])

jet_mults = {
    '3Jets' : '3 jets',
    '4PJets' : '4+ jets',
}

lep_cats = {
    'Tight' : '$e/\mu$',
}

indiv_cat_styles = {
    'Lost_BHad' : {'facecolor' : '#fb9a99', 'name' : 'Lost $b_{h}$'},
    'Lost_BLep' : {'facecolor' : '#e31a1c', 'name' : 'Lost $b_{l}$'},
    'Lost_WJa' : {'facecolor' : '#fdbf6f', 'name' : 'Lost $W_{ja}$'},
    'Lost_WJb' : {'facecolor' : '#ff7f00', 'name' : 'Lost $W_{jb}$'},
    'Merged_BHadBLep' : {'facecolor' : '#a6cee3', 'name' : 'Merged $b_{h}$, $b_{l}$'},
    'Merged_BHadWJa' : {'facecolor' : '#1f78b4', 'name' : 'Merged $b_{h}$, $W_{ja}$'},
    'Merged_BHadWJb' : {'facecolor' : '#b2df8a', 'name' : 'Merged $b_{h}$, $W_{jb}$'},
    'Merged_BLepWJa' : {'facecolor' : '#33a02c', 'name' : 'Merged $b_{l}$, $W_{ja}$'},
    'Merged_BLepWJb' : {'facecolor' : '#cab2d6', 'name' : 'Merged $b_{l}$, $W_{jb}$'},
    'Merged_WJets' : {'facecolor' : '#6a3d9a', 'name' : 'Merged $W_{ja}$, $W_{jb}$'},
    'Other' : {'facecolor' : '#808080', 'name' : 'Other'},
    'Sum unc.' : {
        'facecolor' : 'none',
        'linewidth' : 0,
        'name' : 'Stat. Unc.',
        'hatch' : '///',
    },
}
general_cat_styles = {
    'Lost' : {'facecolor' : '#e41a1c', 'name' : 'Lost'},
    'Merged' : {'facecolor' : '#377eb8', 'name' : 'Partially Merged'},
    'Other' : {'facecolor' : '#808080', 'name' : 'Other'},
    'Sum unc.' : {
        'facecolor' : 'none',
        'linewidth' : 0,
        'name' : 'Stat. Unc.',
        'hatch' : '///',
    },
}

indiv_evt_groups = {
    'Lost_BHad' : ['Lost_BHad'],
    'Lost_BLep' : ['Lost_BLep'],
    'Lost_WJa' : ['Lost_WJa'],
    'Lost_WJb' : ['Lost_WJb'],
    'Merged_BHadBLep' : ['Merged_BHadBLep'],
    'Merged_BHadWJa' : ['Merged_BHadWJa'],
    'Merged_BHadWJb' : ['Merged_BHadWJb'],
    'Merged_BLepWJa' : ['Merged_BLepWJa'],
    'Merged_BLepWJb' : ['Merged_BLepWJb'],
    'Merged_WJets' : ['Merged_WJets'],
    'Other' : ['Other'],
}

evt_groups = {
    'Lost' : ['Lost_BHad', 'Lost_BLep', 'Lost_WJa', 'Lost_WJb'],
    'Merged':['Merged_BHadBLep', 'Merged_BHadWJa', 'Merged_BHadWJb', 'Merged_BLepWJa', 'Merged_BLepWJb', 'Merged_WJets'],
    'Other' : ['Other'],
}

variables = {
    'mtt' : ('$m_{t\\bar{t}}$ [GeV]', 4, (300., 2000.)),
    #'pass_mtt' : ('m($t\\bar{t}$) [GeV]', 4, (200., 2000.)),
    #'fail_mtt' : ('m($t\\bar{t}$) [GeV]', 4, (200., 2000.)),
}

    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open(os.path.join(proj_dir, 'inputs', '%s_lumis_data.json' % base_jobid)).read())[args.year]
lumi_to_use = (data_lumi_year['Muons']+data_lumi_year['Electrons'])/2000.
lumi_correction = load(os.path.join(proj_dir, 'Corrections', base_jobid, 'MC_LumiWeights.coffea'))[args.year]

    ## make groups based on perm category
orig_axis = "cat"
        #all 
all_evt_cat = hist.Cat("all_evt_cat", "All Event Category", sorting='placement')
        # general
g_evt_cat = hist.Cat("gen_evt_cat", "General Event Category", sorting='placement')


    ## make plots
for hname in variables.keys():
    if not hname in hdict.keys():
        raise ValueError("Hist %s not found" % hname)

    #set_trace()
    histo = hdict[hname]
        ## rescale hist by lumi for muons and electrons separately and then combine
    h_mu = histo[:, :, 'Muon', :].integrate('leptype')
    h_mu.scale(lumi_correction['Muons'], axis='dataset')
    h_el = histo[:, :, 'Electron', :].integrate('leptype')
    h_el.scale(lumi_correction['Electrons'], axis='dataset')
    h_tot = h_mu+h_el
    h_tot = h_tot.integrate('dataset')

    if h_tot.dense_dim() == 1:
        xtitle, rebinning, x_lims = variables[hname]

        ## hists should have 3 category axes (jet multiplicity, MT region, perm category) followed by variable
        for jmult in sorted(set([key[0] for key in h_tot.values().keys()])):
            pltdir = os.path.join(outdir, jmult)
            if not os.path.isdir(pltdir):
                os.makedirs(pltdir)
            print(hname, jmult)

            hslice = h_tot[jmult, :].integrate('jmult')
            if rebinning != 1:
                xaxis_name = hslice.dense_axes()[0].name
                hslice = hslice.rebin(xaxis_name, rebinning)

                # make histo with event all individual categories for 'Lost' and 'Merged'
            indiv_cat_histo = hslice.group(orig_axis, all_evt_cat, indiv_evt_groups)
                # hists with yields
            fig, ax = plt.subplots()
            fig.subplots_adjust(hspace=.07)

            Plotter.plot_mc1d(ax=ax, hdict=indiv_cat_histo, xlabel=xtitle, ylabel='Events', xlimits=x_lims, hist_styles=indiv_cat_styles, **{'error_opts':None, 'leg_ncols':2})
            ax.set_ylim(0, ax.get_ylim()[1]*1.1)
                # add lep category
            ax.text(
                0.02, 0.90, "$e/\mu$, %s\nParton Level" % jet_mults[jmult],
                fontsize=rcParams['font.size'], horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
            )
                ## set axes labels and titles
            hep.cms.label(ax=ax, data=False, paper=False, year=args.year, lumi=round(lumi_to_use, 1))

            #set_trace()
            figname = os.path.join(pltdir, '_'.join([args.year, jobid, jmult, 'Indiv_Cats', hname]))
            fig.savefig(figname)
            print('%s written' % figname)

            #set_trace()
                # hists normalized
            fig_norm, ax_norm = plt.subplots()
            fig_norm.subplots_adjust(hspace=.07)

                    # get vals normalized by bin total
            tot_all_vals = indiv_cat_histo.integrate('all_evt_cat').values()[()]
            indiv_norm_histo = indiv_cat_histo.copy()
            for cat in indiv_evt_groups.keys():
                indiv_norm_vals = indiv_cat_histo.values()[(cat,)]/tot_all_vals
                indiv_norm_vals[np.isnan(indiv_norm_vals)] = 0.
                # substitute norm values into hist
                for idx in range(len(indiv_norm_vals)):
                    indiv_norm_histo.values()[(cat,)][idx] = indiv_norm_vals[idx]

            Plotter.plot_mc1d(ax=ax_norm, hdict=indiv_norm_histo, xlabel=xtitle, ylabel='Fraction', xlimits=x_lims, ylimits=(0., 1.), hist_styles=indiv_cat_styles, **{'error_opts':None, 'leg_ncols':2})
            #set_trace()
            ax_norm.set_ylim(0, ax_norm.get_ylim()[1]*1.5)
                # add lep category
            ax_norm.text(
                0.02, 0.90, "$e/\mu$, %s\nParton Level" % jet_mults[jmult],
                fontsize=rcParams['font.size'], horizontalalignment='left', verticalalignment='bottom', transform=ax_norm.transAxes
            )
                ## set axes labels and titles
            hep.cms.label(ax=ax_norm, data=False, paper=False, year=args.year, lumi=round(lumi_to_use, 1))

            #set_trace()
            figname_norm = os.path.join(pltdir, '_'.join([args.year, jobid, jmult, 'Indiv_Cats', hname, 'Norm']))
            fig_norm.savefig(figname_norm)
            print('%s written' % figname_norm)
            plt.close()


            # make histo with event categories 'Lost', 'Merged', 'Other'
            general_cat_histo = hslice.group(orig_axis, g_evt_cat, evt_groups)   

                # hists with yields
            fig, ax = plt.subplots()
            fig.subplots_adjust(hspace=.07)

            Plotter.plot_mc1d(ax=ax, hdict=general_cat_histo, xlabel=xtitle, ylabel='Events', xlimits=x_lims, hist_styles=general_cat_styles, **{'error_opts':None})
            ax.set_ylim(0, ax.get_ylim()[1]*1.1)
                # add lep category
            ax.text(
                0.02, 0.90, "$e/\mu$, %s\nParton Level" % jet_mults[jmult],
                fontsize=rcParams['font.size'], horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes
            )
                ## set axes labels and titles
            hep.cms.label(ax=ax, data=False, paper=False, year=args.year, lumi=round(lumi_to_use, 1))

            #set_trace()
            figname = os.path.join(pltdir, '_'.join([args.year, jobid, jmult, 'Group_Cats', hname]))
            fig.savefig(figname)
            print('%s written' % figname)


                # hists normalized
            fig_norm, ax_norm = plt.subplots()
            fig_norm.subplots_adjust(hspace=.07)

                    # get vals normalized by bin total
            tot_gen_vals = general_cat_histo.integrate('gen_evt_cat').values()[()]
                # substitute norm values into hist
            gen_norm_histo = general_cat_histo.copy()
            for cat in evt_groups.keys():
                gen_norm_vals = general_cat_histo.values()[(cat,)]/tot_gen_vals
                gen_norm_vals[np.isnan(gen_norm_vals)] = 0.
                # substitute norm values into hist
                for idx in range(len(gen_norm_vals)):
                    gen_norm_histo.values()[(cat,)][idx] = gen_norm_vals[idx]

            Plotter.plot_mc1d(ax=ax_norm, hdict=gen_norm_histo, xlabel=xtitle, ylabel='Fraction', xlimits=x_lims, ylimits=(0., 1.), hist_styles=general_cat_styles, **{'error_opts':None})
            ax_norm.set_ylim(0, ax_norm.get_ylim()[1]*1.2)
                # add lep category
            ax_norm.text(
                0.02, 0.90, "$e/\mu$, %s\nParton Level" % jet_mults[jmult],
                fontsize=rcParams['font.size'], horizontalalignment='left', verticalalignment='bottom', transform=ax_norm.transAxes
            )
                ## set axes labels and titles
            hep.cms.label(ax=ax_norm, data=False, paper=False, year=args.year, lumi=round(lumi_to_use, 1))

            #set_trace()
            figname_norm = os.path.join(pltdir, '_'.join([args.year, jobid, jmult, 'Group_Cats', hname, 'Norm']))
            fig_norm.savefig(figname_norm)
            print('%s written' % figname_norm)
            plt.close()


                    ## make table
            rows = [("Lumi: %s fb^-1" % format(round(lumi_to_use, 1), '.1f'), 'e/mu', jmult, "")]
            rows += [("Event Type", "Yield", "Error", "Frac of l+jets ttbar events")]
            for cat in evt_groups.keys():
                rows += [(cat, format(general_cat_histo.values(overflow='all', sumw2=True)[(cat,)][0].sum(), '.3f'), format(np.sqrt(general_cat_histo.values(overflow='all', sumw2=True)[(cat,)][1].sum()), '.3f'), format(general_cat_histo.values(overflow='all', sumw2=True)[(cat,)][0].sum()/general_cat_histo.integrate('gen_evt_cat').values(overflow='all')[()].sum(), '.5f'))]

            rows += [("", "", "", "")]
            for cat in indiv_evt_groups.keys():
                if cat == 'Other': continue
                rows += [(cat, format(indiv_cat_histo.values(overflow='all', sumw2=True)[(cat,)][0].sum(), '.3f'), format(np.sqrt(indiv_cat_histo.values(overflow='all', sumw2=True)[(cat,)][1].sum()), '.3f'), format(indiv_cat_histo.values(overflow='all', sumw2=True)[(cat,)][0].sum()/indiv_cat_histo.integrate('all_evt_cat').values(overflow='all')[()].sum(), '.5f'))]

            #set_trace()
            cat_fractions_name = os.path.join(pltdir, '%s_%s_Event_Categories_%s_yields_and_fractions.txt' % (args.year, jobid, jmult))
            plt_tools.print_table(rows, filename=cat_fractions_name, header_line=1, print_output=True)
            print('%s written' % cat_fractions_name)

