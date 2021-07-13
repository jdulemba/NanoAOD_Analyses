#!/usr/bin/env python

# matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend('agg')
from matplotlib import rcParams
rcParams['font.size'] = 20
rcParams["savefig.format"] = 'png'
rcParams["savefig.bbox"] = 'tight'

from Utilities.Plotter import plot_1D
from coffea.hist import plot
import coffea
from coffea.util import load#, save
from pdb import set_trace
import os
import Utilities.prettyjson as prettyjson
import Utilities.plot_tools as plt_tools
import numpy as np

proj_dir = os.environ['PROJECT_DIR']
base_jobid = os.environ['base_jobid']
analyzer = 'meta_info'

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('year', choices=['2016', '2017', '2018'] if base_jobid == 'NanoAODv6' else ['2016APV', '2016', '2017', '2018'], help='Specify which year to run over')
parser.add_argument('data_mc', choices=['Data', 'MC', 'All'], help='Make plots for MC (MC), make meta.json files for Data (Data), or both (All).')

args = parser.parse_args()

input_dir = os.path.join(proj_dir, 'results', '%s_%s' % (args.year, base_jobid), analyzer)
f_ext = 'TOT.coffea'
outdir = os.path.join(proj_dir, 'plots', '%s_%s' % (args.year, base_jobid), analyzer)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

fnames = sorted(['%s/%s' % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)])
hdict = plt_tools.add_coffea_files(fnames) if len(fnames) > 1 else load(fnames[0])
#set_trace()

isSample = lambda x: ('runs_to_lumis' not in x) and ('PU_nTrueInt' not in x) and ('PU_nPU' not in x)
samples = sorted([key for key in hdict.keys() if isSample(key)])


if (args.data_mc == 'MC') or (args.data_mc == 'All'):
    # get hists
    PU_nTrueInt_histo = hdict['PU_nTrueInt']
    PU_nPU_histo = hdict['PU_nPU']

if (args.data_mc == 'Data') or (args.data_mc == 'All'):
    el_lumimask_check_comb = {}
    el_lumimask_check_ind = {}
    mu_lumimask_check_comb = {}
    mu_lumimask_check_ind = {}
    def as_range(iterable): # not sure how to do this part elegantly
        l = list(iterable)
        if len(l) > 1:
            return [l[0], l[-1]]
        else:
            return [l[0], l[0]]

for sample in samples:
    #print("    %s is being analyzed" % sample)
    if 'data_Single' in sample:
        if not ((args.data_mc == 'Data') or (args.data_mc == 'All')): continue
        print("    %s is being analyzed" % sample)
        #set_trace()
        run_lumi_list = hdict['%s_runs_to_lumis' % sample].value
        lumi_map = {}
        for run, lumis in run_lumi_list:
            lumi_map.setdefault(run, []).append(lumis)

        from itertools import groupby, count
            ## format lumi_map in same way as lumimask golden json
        for run in lumi_map.keys():
            lumis = sorted(list(set([int(item) for sublist in lumi_map[run] for item in sublist])))
            lumi_ranges = [as_range(g) for _, g in groupby(lumis, key=lambda n, c=count(): n-next(c))]
            lumi_map[run] = lumi_ranges

        el_lumimask_check_comb.update(lumi_map) if 'data_SingleElectron' in sample else mu_lumimask_check_comb.update(lumi_map)
        el_lumimask_check_ind.update({sample : lumi_map}) if 'data_SingleElectron' in sample else mu_lumimask_check_ind.update({sample : lumi_map})
        #set_trace()
        continue

    else: 
        if not ((args.data_mc == 'MC') or (args.data_mc == 'All')): continue
        print("    %s is being analyzed" % sample)

            ## write meta info to json
        meta_dict = {}
        for key, val in hdict[sample].items():
            if isinstance(val, (int, float, list)):
                meta_dict[key] = val
            else:
                meta_dict[key] = val.tolist()
        meta_fname = os.path.join(proj_dir, 'inputs', '%s_%s' % (args.year, base_jobid), '%s.meta.json' % sample)
        with open(meta_fname, 'w') as out:
            out.write(prettyjson.dumps(meta_dict))
        print('%s written' % meta_fname)
    
            ## plot histograms
                ## pileup nTrueInt distribution
        pu_nTrueInt_histo = PU_nTrueInt_histo[sample].integrate('dataset')
        pu_nTrueInt_bins = pu_nTrueInt_histo.axis('pu_nTrueInt').edges()
        fig_nTrueInt, ax_nTrueInt = plt.subplots()
        plot_1D(pu_nTrueInt_histo.values()[()], pu_nTrueInt_bins, xlimits=(0., 100.), xlabel=('$\mathsf{%s}$' % pu_nTrueInt_histo.axes()[-1].label), ax=ax_nTrueInt, label=sample)
        ax_nTrueInt.legend(loc='upper right')
        hep.cms.label(ax=ax_nTrueInt, rlabel=args.year)
        figname_nTrueInt = os.path.join(outdir, '%s_PU_nTrueInt.png' % sample)
        fig_nTrueInt.savefig(figname_nTrueInt)
        print('%s written' % figname_nTrueInt)
        plt.close(fig_nTrueInt)

                ## pileup nPU distribution
        pu_nPU_histo = PU_nPU_histo[sample].integrate('dataset')
        pu_nPU_bins = pu_nPU_histo.axis('pu_nPU').edges()
        fig_nPU, ax_nPU = plt.subplots()
        plot_1D(pu_nPU_histo.values()[()], pu_nPU_bins, xlimits=(0., 100.), xlabel=('$\mathsf{%s}$' % pu_nPU_histo.axes()[-1].label), ax=ax_nPU, label=sample)
        ax_nPU.legend(loc='upper right')
        hep.cms.label(ax=ax_nPU, rlabel=args.year)
        figname_nPU = os.path.join(outdir, '%s_PU_nPU.png' % sample)
        fig_nPU.savefig(figname_nPU)
        print('%s written' % figname_nPU)
        plt.close(fig_nPU)


if (args.data_mc == 'Data') or (args.data_mc == 'All'):
    #set_trace()
        ## save individual run.json files for each data period
    for period in sorted(el_lumimask_check_ind.keys()):
        el_lumi_map_ind_dict = {}
        for key, val in sorted(el_lumimask_check_ind[period].items()):
            if isinstance(val, (int, float, list)):
                el_lumi_map_ind_dict[key] = val
            else:
                el_lumi_map_ind_dict[key] = val.tolist()
        with open(os.path.join(proj_dir,'inputs', '%s_%s' % (args.year, base_jobid), '%s.run.json' % period), 'w') as out:
            out.write(prettyjson.dumps(el_lumi_map_ind_dict))
        print('%s.run.json written' % period)

    for period in sorted(mu_lumimask_check_ind.keys()):
        mu_lumi_map_ind_dict = {}
        for key, val in sorted(mu_lumimask_check_ind[period].items()):
            if isinstance(val, (int, float, list)):
                mu_lumi_map_ind_dict[key] = val
            else:
                mu_lumi_map_ind_dict[key] = val.tolist()
        with open(os.path.join(proj_dir,'inputs', '%s_%s' % (args.year, base_jobid), '%s.run.json' % period), 'w') as out:
            out.write(prettyjson.dumps(mu_lumi_map_ind_dict))
        print('%s.run.json written' % period)

        ## save combined run.json files
    el_lumi_map_comb_dict = {}
    for key, val in sorted(el_lumimask_check_comb.items()):
        if isinstance(val, (int, float, list)):
            el_lumi_map_comb_dict[key] = val
        else:
            el_lumi_map_comb_dict[key] = val.tolist()
    with open(os.path.join(proj_dir,'inputs', '%s_%s' % (args.year, base_jobid), 'data_SingleElectron_%s.run.json' % args.year), 'w') as out:
        out.write(prettyjson.dumps(el_lumi_map_comb_dict))
    print('data_SingleElectron_%s.run.json written' % args.year)

    mu_lumi_map_comb_dict = {}
    for key, val in sorted(mu_lumimask_check_comb.items()):
        if isinstance(val, (int, float, list)):
            mu_lumi_map_comb_dict[key] = val
        else:
            mu_lumi_map_comb_dict[key] = val.tolist()
    with open(os.path.join(proj_dir,'inputs', '%s_%s' % (args.year, base_jobid), 'data_SingleMuon_%s.run.json' % args.year), 'w') as out:
        out.write(prettyjson.dumps(mu_lumi_map_comb_dict))
    print('data_SingleMuon_%s.run.json written' % args.year)
