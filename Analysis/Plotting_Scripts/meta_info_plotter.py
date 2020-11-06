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

import Utilities.Plotter as Plotter
from coffea.hist import plot
import coffea
from coffea.util import load, save
from pdb import set_trace
import os
import Utilities.prettyjson as prettyjson
import Utilities.plot_tools as plt_tools
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('year', choices=['2016', '2017', '2018'], help='Specify which year to run over')
parser.add_argument('--dump_lumi', action='store_true', help='Crosscheck resulting lumi map from data to golden json.')

args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'meta_info'

input_dir = os.path.join(proj_dir, 'results', '%s_%s' % (args.year, jobid), analyzer)
f_ext = 'TOT.coffea'
outdir = os.path.join(proj_dir, 'plots', '%s_%s' % (args.year, jobid), analyzer)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

fnames = sorted(['%s/%s' % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)])
hdict = plt_tools.add_coffea_files(fnames) if len(fnames) > 1 else load(fnames[0])
#set_trace()

isSample = lambda x: ('runs_to_lumis' not in x) and ('PUDistribution' not in x)
samples = sorted([key for key in hdict.keys() if isSample(key)])

    # get hists
PU_histo = hdict['PUDistribution']

if args.dump_lumi:
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
    print("    %s is being analyzed" % sample)
    if 'data_Single' in sample:
        if not args.dump_lumi: continue
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
        #if args.dump_lumi: continue

            ## write meta info to json
        meta_dict = {}
        for key, val in hdict[sample].items():
            if isinstance(val, (int, float, list)):
                meta_dict[key] = val
            else:
                meta_dict[key] = val.tolist()
        meta_fname = os.path.join(proj_dir, 'inputs', '%s_%s' % (args.year, jobid), '%s.meta.json' % sample)
        with open(meta_fname, 'w') as out:
            out.write(prettyjson.dumps(meta_dict))
        print('%s written' % meta_fname)
    
            ## plot histograms
                ## pileup distribution
        pu_histo = PU_histo[sample].integrate('dataset')
        pu_bins = pu_histo.axis('pu').edges()
        fig_pu, ax_pu = plt.subplots()
        Plotter.plot_1D(pu_histo.values()[()], pu_bins, (pu_bins.min(), pu_bins.max()), xlabel=('$\mathrm{%s}$' % pu_histo.axes()[-1].label), ax=ax_pu, label=sample, histtype='step')
        ax_pu.legend(loc='upper right')
        figname_pu = os.path.join(outdir, '%s_PUDistribution.png' % sample)
        fig_pu.savefig(figname_pu)
        print('%s written' % figname_pu)
        plt.close(fig_pu)

if args.dump_lumi:
    #set_trace()
        ## save individual run.json files for each data period
    for period in sorted(el_lumimask_check_ind.keys()):
        el_lumi_map_ind_dict = {}
        for key, val in sorted(el_lumimask_check_ind[period].items()):
            if isinstance(val, (int, float, list)):
                el_lumi_map_ind_dict[key] = val
            else:
                el_lumi_map_ind_dict[key] = val.tolist()
        with open(os.path.join(proj_dir,'inputs', '%s_%s' % (args.year, jobid), '%s.run.json' % period), 'w') as out:
            out.write(prettyjson.dumps(el_lumi_map_ind_dict))
        print('%s.run.json written' % period)

    for period in sorted(mu_lumimask_check_ind.keys()):
        mu_lumi_map_ind_dict = {}
        for key, val in sorted(mu_lumimask_check_ind[period].items()):
            if isinstance(val, (int, float, list)):
                mu_lumi_map_ind_dict[key] = val
            else:
                mu_lumi_map_ind_dict[key] = val.tolist()
        with open(os.path.join(proj_dir,'inputs', '%s_%s' % (args.year, jobid), '%s.run.json' % period), 'w') as out:
            out.write(prettyjson.dumps(mu_lumi_map_ind_dict))
        print('%s.run.json written' % period)

        ## save combined run.json files
    el_lumi_map_comb_dict = {}
    for key, val in sorted(el_lumimask_check_comb.items()):
        if isinstance(val, (int, float, list)):
            el_lumi_map_comb_dict[key] = val
        else:
            el_lumi_map_comb_dict[key] = val.tolist()
    with open(os.path.join(proj_dir,'inputs', '%s_%s' % (args.year, jobid), 'data_SingleElectron_%s.run.json' % args.year), 'w') as out:
        out.write(prettyjson.dumps(el_lumi_map_comb_dict))
    print('data_SingleElectron_%s.run.json written' % args.year)

    mu_lumi_map_comb_dict = {}
    for key, val in sorted(mu_lumimask_check_comb.items()):
        if isinstance(val, (int, float, list)):
            mu_lumi_map_comb_dict[key] = val
        else:
            mu_lumi_map_comb_dict[key] = val.tolist()
    with open(os.path.join(proj_dir,'inputs', '%s_%s' % (args.year, jobid), 'data_SingleMuon_%s.run.json' % args.year), 'w') as out:
        out.write(prettyjson.dumps(mu_lumi_map_comb_dict))
    print('data_SingleMuon_%s.run.json written' % args.year)
