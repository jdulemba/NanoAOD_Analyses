#!/usr/bin/env python

import os, subprocess, time
from pdb import set_trace
import Run_Jobs.tools as tools
import Utilities.prettyjson as prettyjson
import fnmatch

import argparse
class ParseKwargs(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, dict())
        for value in values:
            key, value = value.split('=')
            getattr(namespace, self.dest)[key] = value

parser = argparse.ArgumentParser('submit analyzer to the batch queues')
parser.add_argument('analyzer', help='Analyzer to use.')
parser.add_argument('frange', type=str, help='Specify start:stop indices for files, inclusive. 0:1 means valid files in indices 0 and 1 will be used')
parser.add_argument('year', choices=['2016APV', '2016', '2017', '2018'], help='Specify which year to run over')
parser.add_argument('--opts', nargs='*', action=ParseKwargs, help='Options to pass to analyzers.')
#parser.add_argument('--sample', type=str, help='String used to specify which samples to run, can be name or regex.')
#parser.add_argument('--outfname', type=str, help='Specify output filename, including directory and file extension')
#parser.add_argument('--debug', action='store_true', help='Uses iterative_executor for debugging purposes, otherwise futures_excutor will be used (faster)')
#parser.add_argument('--signal', type=str, default='.*', help='Signal sample to use, regex formatting.')
#parser.add_argument('--evt_sys', type=str, default='NONE', help='Specify event systematics to run, will be capitalized. Default is NONE')
#parser.add_argument('--rewt_sys', type=str, default='NONE', help='Specify reweighting systematics to run, will be capitalized. Default is NONE')
#parser.add_argument('--only_sys', type=int, default=0, help='Only run specified systematics and not nominal weights (nosys)')
args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
base_jobid = os.environ['base_jobid']
analyzer=args.analyzer
#set_trace()

# define dictionary of options to pass
opts_dict = args.opts

    ## get samples to use
indir = os.path.join(proj_dir, 'inputs', '%s_%s' % (args.year, base_jobid))
samples_to_use = tools.get_sample_list(indir=indir, sample=opts_dict['sample']) if 'sample' in opts_dict.keys() else tools.get_sample_list(indir=indir, text_file='analyzer_inputs.txt')

fileset = {}
for sample in samples_to_use:
    if not os.path.isfile(sample):
        raise IOError("Sample file %s not found" % sample)

    sample_name = sample.split('/')[-1].split('.')[0]
    sfiles = open(sample, 'r')
    valid_inds_fnames = [(idx, fname.strip('\n')) for idx, fname in enumerate(sfiles) if not fname.startswith('#')]
    valid_inds = [x[0] for x in valid_inds_fnames]
    valid_files = [x[1] for x in valid_inds_fnames]

    if ':' in args.frange:
        file_start, file_stop = int((args.frange).split(':')[0]), int((args.frange).split(':')[1])
    else:
        file_start = 0
        file_stop = len(valid_files)-1 if (args.frange).lower() == 'all' else int(args.frange)-1

    if file_start >= 0 and file_stop < len(valid_files):
        valid_files = valid_files[file_start:file_stop+1]
        inds_to_use = valid_inds[file_start:file_stop+1]
    else:
        raise IOError("The number of root files available for the %s sample is %i. args.frange must be less than this." % (sample, len(valid_files) ) )

    #print(inds_to_use)
    fileset[sample_name] = valid_files

print(fileset)

    ## save output to coffea pkl file
#set_trace()
if (args.frange).lower() == 'all':
    outdir = os.path.join(proj_dir, 'results', '%s_%s' % (args.year, jobid), analyzer)
    if 'outfname' in opts_dict.keys():
        cfname = opts_dict['outfname']
    elif 'sample' in opts_dict.keys():
        cfname = os.path.join(outdir, '%s.coffea' % opts_dict['sample'])
    else:
        cfname = os.path.join(outdir,'test_%s.coffea' % analyzer)
else:
    if ':' in args.frange:
        outdir = os.path.join(proj_dir, 'results', '%s_%s' % (args.year, jobid), analyzer)
        if 'outfname' in opts_dict.keys():
            cfname = opts_dict['outfname']
        elif ('signal' in opts_dict.keys()) and ('sample' in opts_dict.keys()) and ((analyzer == 'htt_signal_reweight') or (analyzer == 'signal_validation')):
            cfname = os.path.join(outdir,'%s_%s_%sto%s.coffea' % (opts_dict['signal'], opts_dict['sample'], file_start, file_stop)) if not opts_dict['signal'] == ".*" else os.path.join(outdir, 'AHtoTT_%s_%sto%s.coffea' % (opts_dict['sample'], file_start, file_stop))
        #elif args.signal and args.sample and ((analyzer == 'htt_signal_reweight') or (analyzer == 'signal_validation')):
        ##elif args.signal and args.sample and analyzer == 'htt_signal_reweight':
        #    cfname = '%s/%s_%s_%sto%s.coffea' % (outdir, args.signal, args.sample, file_start, file_stop) if not args.signal == parser.get_default('signal') else '%s/AHtoTT_%s_%sto%s.coffea' % (outdir, args.sample, file_start, file_stop)
        elif 'sample' in opts_dict.keys():
            cfname = os.path.join(outdir, '%s_%sto%s.coffea' % (opts_dict['sample'], file_start, file_stop))
        else:
            cfname = os.path.join(outdir, 'test_%sto%s.coffea' % (file_start, file_stop))
    else:
        outdir = proj_dir
        if 'outfname' in opts_dict.keys():
            cfname = opts_dict['outfname']
        elif ('signal' in opts_dict.keys()) and ('sample' in opts_dict.keys()) and ((analyzer == 'htt_signal_reweight') or (analyzer == 'signal_validation')):
            cfname = os.path.join(outdir, '%s_%s_%s_%s.test.coffea' % (opts_dict['signal'], opts_dict['sample'], args.year, analyzer)) if not opts_dict['signal'] == ".*" else os.path.join(outdir, 'AHtoTT_%s_%s_%s.test.coffea' % (opts_dict['sample'], args.year, analyzer))
        #elif args.signal and args.sample and ((analyzer == 'htt_signal_reweight') or (analyzer == 'signal_validation')):
        ##elif args.signal and args.sample and analyzer == 'htt_signal_reweight':
        #    cfname = '%s/%s_%s_%s_%s.test.coffea' % (outdir, args.signal, args.sample, args.year, analyzer) if not args.signal == parser.get_default('signal') else '%s/AHtoTT_%s_%s_%s.test.coffea' % (outdir, args.sample, args.year, analyzer)
        elif 'sample' in opts_dict.keys():
            cfname = os.path.join(outdir, '%s_%s_%s.test.coffea' % (opts_dict['sample'], args.year, analyzer))
        else:
            cfname = os.path.join(outdir, '%s_%s.test.coffea' % (args.year, analyzer))
if not os.path.isdir(outdir):
    os.makedirs(outdir)

##set_trace()
#    ## get string for passed options
#to_debug = bool(opts_dict.get('debug'))
#evt_sys = opts_dict.get('evt_sys', 'NONE').upper()
#rewt_sys = opts_dict.get('rewt_sys', 'NONE').upper()
#only_sys = opts_dict.get('only_sys', 'False')

#opts = ""
##if args.sample: opts += " --sample=%s" % args.sample
#if args.debug: opts += " --debug"
#if to_debug: opts += " --debug"

if (analyzer == 'htt_signal_reweight') or (analyzer == 'signal_validation'):
#if analyzer == 'htt_signal_reweight':
    signal = opts_dict.get('signal', None)
    if signal is None:
    #if args.signal is None:
        raise ValueError("Signal sample must be specified when running %s" % analyzer)
    opts += " --evt_sys={EVTSYS} --rewt_sys={REWTSYS}".format(
            #EVTSYS=args.evt_sys.upper(),
            #REWTSYS=args.rewt_sys.upper(),
            EVTSYS=evt_sys.upper(),
            REWTSYS=rewt_sys.upper(),
        )
    opts += " --only_sys=%i" % only_sys
    #opts += " --only_sys=%i" % args.only_sys
    run_cmd = """python {PROJDIR}/bin/{ANALYZER}.py "{FSET}" {YEAR} {OUTFNAME} --signal={SIGNAL} {OPTS}""".format(
            PROJDIR=proj_dir,
            ANALYZER=analyzer,
            FSET=fileset,
            YEAR=args.year,
            SIGNAL=signal,
            #SIGNAL=args.signal,
            OUTFNAME=cfname,
            OPTS=opts
    )
#elif analyzer == 'htt_btag_iso_cut':
#elif (analyzer == 'htt_btag_iso_cut') or (analyzer == 'evtWeights'):
elif (analyzer == 'htt_btag_iso_cut') or (analyzer == 'evtWeights') or (analyzer == 'htt_btag_sb_regions'):
    run_cmd = """python {PROJDIR}/bin/{ANALYZER}.py "{FSET}" {YEAR} {OUTFNAME} "{OPTS}" """.format(
            PROJDIR=proj_dir,
            ANALYZER=analyzer,
            FSET=fileset,
            YEAR=args.year,
            OUTFNAME=cfname,
            OPTS=opts_dict
    )
elif analyzer == 'ttbar_post_alpha_reco':
    run_cmd = """python {PROJDIR}/bin/{ANALYZER}.py "{FSET}" {YEAR} {OUTFNAME} "{OPTS}" """.format(
            PROJDIR=proj_dir,
            ANALYZER=analyzer,
            FSET=fileset,
            YEAR=args.year,
            OUTFNAME=cfname,
            OPTS=opts_dict
    )
else:
    run_cmd = """python {PROJDIR}/bin/{ANALYZER}.py "{FSET}" {YEAR} {OUTFNAME} "{OPTS}" """.format(
            PROJDIR=proj_dir,
            ANALYZER=analyzer,
            FSET=fileset,
            YEAR=args.year,
            OUTFNAME=cfname,
            OPTS=opts_dict
    )

print('Running command: %s' % run_cmd)
os.system(run_cmd)
