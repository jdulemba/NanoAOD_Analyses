from coffea.hist import plot
import coffea
# matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from coffea.util import load, save
from pdb import set_trace
import os
import Utilities.prettyjson as prettyjson
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('year', choices=['2016', '2017', '2018'], help='Specify which year to run over')
parser.add_argument('--sample', type=str, help='Input sample to use.')
parser.add_argument('--testing', action='store_true', help='Determines where input file is.')
parser.add_argument('--dump_lumi', action='store_true', help='Crosscheck resulting lumi map from data to golden json.')

args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'meta_info'

variables = {
    'mtt' : '$m_{t\\bart}$ [GeV]',
    'ctstar' : 'cos($\\theta_{t}^{*}$)'
}

input_dir = proj_dir if args.testing else '/'.join([proj_dir, 'results', '%s_%s' % (args.year, jobid), analyzer])
f_ext = '%s.test.coffea' % analyzer if args.testing else '.coffea'
outdir = '/'.join([proj_dir, 'plots', '%s_%s' % (args.year, jobid), analyzer, 'Test']) if args.testing else '/'.join([proj_dir, 'plots', '%s_%s' % (args.year, jobid), analyzer])

if args.testing:
    fnames = ['%s/%s_%s' % (input_dir, args.sample, f_ext)] if args.sample else ['%s/%s' % (input_dir, f_ext)]
else:
    fnames = ['%s/%s%s' % (input_dir, args.sample, f_ext)] if args.sample else ['%s/%s' % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)]
fnames = sorted(fnames)
#set_trace()

if not os.path.isdir(outdir):
    os.makedirs(outdir)

if args.dump_lumi:
    el_lumimask_check = {}
    mu_lumimask_check = {}
    def as_range(iterable): # not sure how to do this part elegantly
        l = list(iterable)
        if len(l) > 1:
            return [l[0], l[-1]]
        else:
            return [l[0], l[0]]

for fname in fnames:
    if not os.path.isfile(fname):
        raise IOError("%s not found" % fname)
    hists = load(fname)

    #set_trace()
    if 'data_Single' in fname:
        print(fname.split('/')[-1].split('.')[0])
        if not args.dump_lumi: continue
        run_lumi_list = hists['%s_runs_to_lumis' % fname.split('/')[-1].split('.')[0]].value
        lumi_map = {}
        for run, lumis in run_lumi_list:
            lumi_map.setdefault(run, []).append(lumis)

        from itertools import groupby, count

        #set_trace()
            ## format lumi_map in same way as lumimask golden json
        for run in lumi_map.keys():
            lumis = sorted(list(set([int(item) for sublist in lumi_map[run] for item in sublist])))
            lumi_ranges = [as_range(g) for _, g in groupby(lumis, key=lambda n, c=count(): n-next(c))]
            lumi_map[run] = lumi_ranges

        el_lumimask_check.update(lumi_map) if 'data_SingleElectron' in fname else mu_lumimask_check.update(lumi_map)
        #set_trace()
        continue
    #set_trace()

    else: 
        if args.dump_lumi: continue
        for hname in hists.keys():
            if 'runs_to_lumis' in hname: continue
                ## plot histograms
            if isinstance(hists[hname], coffea.hist.hist_tools.Hist):
                histo = hists[hname]
    
                if histo.dense_dim() == 1:
                    ## make plot for separate samples
                    for sample in histo.axes()[0]._sorted:
                        sample_histo = histo[sample]
        
                        fig = plt.figure()
                        ax = plot.plot1d(sample_histo)
                        if 'mtt' in hname:
                            plt.xlabel(variables['mtt'])
                        elif 'ctstar' in hname:
                            plt.xlabel(variables['ctstar'])
                        else:
                            plt.xlabel('$%s$' % sample_histo.axes()[-1].label)
                        plt.ylabel('Events')
                        figname = '%s/%s_%s.png' % (outdir, sample, hname)
                        fig.savefig(figname)
                        print('%s written' % figname)
                        plt.close()
        
        
                elif histo.dense_dim() == 2:
                    xvar, yvar = histo.axes()[-2].name, histo.axes()[-1].name
        
                    ## make plots for different ttJets samples
                    for sample in histo.axes()[0]._sorted:
                        if not (sample == 'ttJets' or sample == 'ttJets_PS'): continue
                        sample_histo = histo[sample].project(histo.axes()[1].name, xvar, yvar)
    
                        #set_trace()    
                            ## plot x projection
                        fig = plt.figure()
                        x_proj_histo = sample_histo.sum(yvar)
                        x_ax = plot.plot1d(x_proj_histo)
                        x_ax.set_xlabel(variables[xvar])
                        x_ax.set_ylabel('Events')
                        xfigname = '%s/%s_%s.png' % (outdir, sample, xvar)
                        fig.savefig(xfigname)
                        print('%s written' % xfigname)
                        plt.close()
        
                            ## plot y projection
                        fig = plt.figure()
                        y_proj_histo = sample_histo.sum(xvar)
                        y_ax = plot.plot1d(y_proj_histo)
                        y_ax.set_xlabel(variables[yvar])
                        y_ax.set_ylabel('Events')
                        yfigname = '%s/%s_%s.png' % (outdir, sample, yvar)
                        fig.savefig(yfigname)
                        print('%s written' % yfigname)
                        plt.close()
        
                ## write other meta info to meta.json files
            else:
                if args.testing: continue
                #set_trace()
                #if args.testing: set_trace()
                meta_dict = {}
                for key, val in hists[hname].items():
                    if isinstance(val, (int, float, list)):
                        meta_dict[key] = val
                    else:
                        meta_dict[key] = val.tolist()
                with open('%s/inputs/%s/%s.meta.json' % (proj_dir, '%s_%s' % (args.year, jobid), hname), 'w') as out:
                    out.write(prettyjson.dumps(meta_dict))
                    print('%s/inputs/%s/%s.meta.json written' % (proj_dir, '%s_%s' % (args.year, jobid), hname))
    

if args.dump_lumi:
    el_lumi_map_dict = {}
    for key, val in el_lumimask_check.items():
        if isinstance(val, (int, float, list)):
            el_lumi_map_dict[key] = val
        else:
            el_lumi_map_dict[key] = val.tolist()
    with open('%s/inputs/%s_%s/data_SingleElectron_%s.run.json' % (proj_dir, args.year, jobid, args.year), 'w') as out:
        out.write(prettyjson.dumps(el_lumi_map_dict))
    print('data_SingleElectron_%s.run.json written' % args.year)

    mu_lumi_map_dict = {}
    for key, val in mu_lumimask_check.items():
        if isinstance(val, (int, float, list)):
            mu_lumi_map_dict[key] = val
        else:
            mu_lumi_map_dict[key] = val.tolist()
    with open('%s/inputs/%s_%s/data_SingleMuon_%s.run.json' % (proj_dir, args.year, jobid, args.year), 'w') as out:
        out.write(prettyjson.dumps(mu_lumi_map_dict))
    print('data_SingleMuon_%s.run.json written' % args.year)


