from coffea.hist import plot
# matplotlib
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from coffea.util import load
from pdb import set_trace
import os
import styles
import Utilities.plot_tools as plt_tools
import Utilities.prettyjson as prettyjson
import re
from coffea import hist

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('year', choices=['2016', '2017', '2018'], help='What year is the ntuple from.')
parser.add_argument('lepton', choices=['Electron', 'Muon'], help='Choose which lepton to make plots for')
parser.add_argument('--sample', type=str, help='Input sample to use.')
parser.add_argument('--testing', action='store_true', help='Determines where input file is.')
parser.add_argument('--use_combined', action='store_true', help='Used file that has multiple datasets.')

args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'htt_simple'

input_dir = proj_dir if args.testing else '/'.join([proj_dir, 'results', '%s_%s' % (args.year, jobid), analyzer])
f_ext = '%s.test.coffea' % analyzer if args.testing else 'TOT.coffea'
outdir = '/'.join([proj_dir, 'plots', '%s_%s' % (args.year, jobid), analyzer, 'Test']) if args.testing else '/'.join([proj_dir, 'plots', '%s_%s' % (args.year, jobid), analyzer])
if not os.path.isdir(outdir):
    os.makedirs(outdir)

if args.testing:
    if args.use_combined:
        fnames = ['%s/%s' % (input_dir, f_ext)]
    else:
        fnames = ['%s/%s' % (input_dir, fname) for fname in os.listdir(input_dir) if fname == '%s_%s' % (args.sample, f_ext)] if args.sample else ['%s/%s' % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)]
        fnames.remove('%s/%s' % (input_dir, f_ext))
else:
    fnames = ['%s/%s%s' % (input_dir, args.sample, f_ext)] if args.sample else ['%s/%s' % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)]
fnames = sorted(fnames)

#set_trace()
hdict = plt_tools.add_coffea_files(fnames) if len(fnames) > 1 else load(fnames[0])

jet_mults = {
    '3Jets' : '3 jets',
    '4PJets' : '4+ jets'
}

objtypes = {
    'Jets' : 'jets',
    'Lep' :  {
        'Muon' : '$\\mu$',
        'Electron' : '$e$',
    }
}

variables = {
    'Jets_LeadJet_pt' : ('$p_{T}$(leading jet) [GeV]', 2, (0., 300.), True),
    'Jets_LeadJet_eta' : ('$\\eta$(leading jet)', 1, (-2.6, 2.6), True),
    'Jets_LeadJet_phi' : ('$\\phi$(leading jet)', 1, (-4., 4.), True),
    'Jets_LeadJet_energy' : ('E(leading jet) [GeV]', 2, (0., 500.), True),
    'Jets_pt' : ('$p_{T}$(jets) [GeV]', 2, (0., 300.), True),
    'Jets_eta' : ('$\\eta$(jets)', 1, (-2.6, 2.6), True),
    'Jets_phi' : ('$\\phi$(jets)', 1, (-4., 4.), True),
    'Jets_energy' : ('E(jets) [GeV]', 2, (0., 500.), True),
    'Jets_njets' : ('$n_{jets}$', 1, (0, 15), True),
    'Lep_pt' : ('$p_{T}$(%s) [GeV]' % objtypes['Lep'][args.lepton], 2, (0., 300.), True),
    'Lep_eta' : ('$\\eta$(%s)' % objtypes['Lep'][args.lepton], 1, (-2.6, 2.6), True),
    'Lep_phi' : ('$\\phi$(%s)' % objtypes['Lep'][args.lepton], 1, (-4., 4.), True),
    'Lep_energy' : ('E(%s) [GeV]' % objtypes['Lep'][args.lepton], 2, (0., 500.), True),
    'Lep_iso' : ('pfRelIso, %s' % objtypes['Lep'][args.lepton], 1, (0., 1.), True),
    'BTag_SF' : ('$SF_{btag}$', 1, (0.5, 1.5), True),
    'Lep_SF' : ('$SF_{l}$', 1, (0.8, 1.1), True),
}
if args.lepton == 'Electron':
    variables.update({'Lep_etaSC' : ('$\\eta_{SC}$(%s)' % objtypes['Lep'][args.lepton], 1, (-2.6, 2.6), True)})


    ## get plotting colors/settings
hstyles = styles.styles
stack_fill_opts = {'alpha': 0.8, 'edgecolor':(0,0,0,.5)}
stack_error_opts = {'edgecolor':(0,0,0,.5)}

    ## get data lumi and scale MC by lumi
data_lumi_year = prettyjson.loads(open('%s/inputs/lumis_data.json' % proj_dir).read())[args.year]
lumi_correction = load('%s/Corrections/%s/MC_LumiWeights_IgnoreSigEvts.coffea' % (proj_dir, jobid))
for hname in hdict.keys():
    if hname == 'cutflow': continue
    hdict[hname].scale(lumi_correction[args.year]['%ss' % args.lepton], axis='dataset')


## make data and mc categories for data/MC plotting
mc_samples = re.compile('(?!data*)')
data_samples = re.compile('(data*)')

## make groups based on process
process = hist.Cat("process", "Process", sorting='placement')
process_cat = "dataset"
process_groups = plt_tools.make_dataset_groups(args.lepton, args.year)
#set_trace()
for hname in hdict.keys():
    if hname == 'cutflow': continue
    hdict[hname] = hdict[hname].group(process_cat, process, process_groups)
    

def get_samples_yield_and_frac(histo, lep):
    '''
    Get the yield and relative fraction for each sample of MC, get data yield and compare data/MC
    Returns: list of tuples containing sample name, yield, and fraction
    '''
    yields = histo.integrate('njets').sum().values()
    proc_yields_list = [(''.join(process), proc_yields) for process, proc_yields in yields.items()]
    mc_yield = sum([process[1] for process in proc_yields_list if not (process[0] == 'data')])
    data_yield = sum([process[1] for process in proc_yields_list if (process[0] == 'data')])

    rows = [("Lumi: %s fb^-1" % format(data_lumi_year['%ss' % lep]/1000., '.1f'), "Sample", "Yield", "Frac")]
    rows += [("", process, format(proc_yield, '.1f'), format((proc_yield/mc_yield)*100, '.1f')) for process, proc_yield in proc_yields_list if not process == 'data']
    rows += [("", "SIM", format(mc_yield, '.1f'), '100.0')]
    rows += [("", "", "", "")]
    rows += [("", "data", format(data_yield, '.1f'), "")]
    rows += [("", "", "", "")]
    rows += [("", "data/SIM", "", format(data_yield/mc_yield, '.3f'))]
        
    yields_dict = {process:round(proc_yield, 2) for process, proc_yield in proc_yields_list}
    yields_dict.update({'SIM': round(mc_yield, 2)})
    yields_dict.update({'data/SIM': round(data_yield/mc_yield, 3)})

    return rows, yields_dict



    ## make plots
for hname in hdict.keys():
    if (hname not in variables.keys()): continue
    #if hname == 'cutflow': continue
    histo = hdict[hname]
    #set_trace()

    if histo.dense_dim() == 1:
        xtitle, rebinning, x_lims, withData = variables[hname]

        ## hists should have 4 category axes (dataset, jet multiplicity, btagger, lepton type) followed by variable 
        for lep in [args.lepton]:
            for btagger in histo.axis('bdisc')._sorted:
                for jmult in histo.axis('jmult')._sorted:
                    pltdir = outdir if args.testing else '/'.join([outdir, lep, jmult, btagger])
                    if not os.path.isdir(pltdir):
                        os.makedirs(pltdir)

                    if withData:
                        fig, (ax, rax) = plt.subplots(2, 1, figsize=(7,7), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                        fig.subplots_adjust(hspace=.07)
                    else:
                        fig, ax = plt.subplots(1, 1, figsize=(7,7))#, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                        fig.subplots_adjust(hspace=.07)
                    hslice = histo[:, jmult, lep, btagger].integrate('jmult').integrate('leptype').integrate('bdisc')

                    #set_trace()
                    #if 'Lep_eta' in hname: set_trace()
                    if hname == 'Jets_njets':
                        print(jmult)
                        yields_txt, yields_json = get_samples_yield_and_frac(hslice, lep)
                        txt_name = '%s/%s_%s_%s_yields_and_fracs.txt' % (pltdir, btagger, jmult, lep)
                        plt_tools.print_table(yields_txt, filename=txt_name, print_output=True)
                        with open(txt_name, 'w') as out:
                            out.write(prettyjson.dumps(yields_json))

                    if rebinning != 1:
                        xaxis_name = hslice.dense_axes()[0].name
                        hslice = hslice.rebin(xaxis_name, rebinning)

                        ## plot MC and data
                    plot.plot1d(hslice[mc_samples],
                        overlay=hslice.axes()[0].name,
                        ax=ax,
                        clear=False,
                        stack=True,
                        line_opts=None,
                        fill_opts=stack_fill_opts,
                        error_opts=stack_error_opts
                    )
                    #set_trace()
                    if withData:
                        plot.plot1d(hslice[data_samples],
                            overlay=hslice.axes()[0].name,
                            ax=ax,
                            clear=False,
                            error_opts=hstyles['data_err_opts']
                        )
                    ax.autoscale(axis='x', tight=True)
                    ax.set_ylim(0, None)
                    ax.set_xlabel(None)
                    if hname == 'BTag_SF':
                        x_lims = (0.5, 2.5) if btagger == 'DeepJet' else (0.7, 1.2)
                    ax.set_xlim(x_lims)

                        ## set legend and corresponding colors
                    handles, labels = ax.get_legend_handles_labels()
                    for idx, sample in enumerate(labels):
                        if sample == 'data' or sample == 'Observed': continue
                        facecolor, legname = plt_tools.get_styles(sample, hstyles)
                        handles[idx].set_facecolor(facecolor)
                        labels[idx] = legname
                    # call ax.legend() with the new values
                    ax.legend(handles,labels, loc='upper right')
                    #set_trace()

                    if withData:
                            ## plot data/MC ratio
                        plot.plotratio(hslice[data_samples].sum(hslice.axes()[0].name), hslice[mc_samples].sum(hslice.axes()[0].name), 
                            ax=rax,
                            error_opts=hstyles['data_err_opts'], 
                            #error_opts=data_err_opts, 
                            denom_fill_opts={},
                            guide_opts={},
                            unc='num'
                        )
                        rax.set_ylabel('data/MC')
                        rax.set_ylim(0.5, 1.5)
                        rax.set_xlim(x_lims)

                        # add btagger name
                    ax.text(
                        0.02, 0.95, btagger,
                        fontsize=10,
                        horizontalalignment='left',
                        verticalalignment='bottom',
                        transform=ax.transAxes
                    )

                        ## set axes labels and titles
                    plt.xlabel(xtitle)
                    cms_blurb = plt.text(
                        0., 1., r"CMS Preliminary",
                        fontsize=12, 
                        horizontalalignment='left', 
                        verticalalignment='bottom', 
                        transform=ax.transAxes,
                        style='italic'
                    )
                    lumi_blurb = plt.text(
                        1., 1., r"(13 TeV %.2f fb$^{-1}$, %s/%s)" % (data_lumi_year['%ss' % lep]/1000., objtypes['Lep'][lep], jet_mults[jmult]),
                        fontsize=12, 
                        horizontalalignment='right', 
                        verticalalignment='bottom', 
                        transform=ax.transAxes
                    )

                    figname = '%s/%s.png' % (pltdir, '_'.join([btagger, jmult, lep, hname]))
                    fig.savefig(figname)
                    print('%s written' % figname)
                    plt.close()
                    #set_trace()


