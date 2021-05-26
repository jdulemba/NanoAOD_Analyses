from coffea.hist import plot
# matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend('agg')
from matplotlib import rcParams
rcParams["savefig.format"] = 'png'
rcParams["savefig.bbox"] = 'tight'
from coffea.util import load
from pdb import set_trace
import os
import Utilities.plot_tools as plt_tools
import Utilities.prettyjson as prettyjson
import re
from coffea import hist
import Utilities.Plotter as Plotter
from Utilities import styles

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
analyzer = 'presel_analyzer'

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
    'MT' : ('$M_{T}$', 1, (0, 300), True),
}
if args.lepton == 'Electron':
    variables.update({'Lep_etaSC' : ('$\\eta_{SC}$(%s)' % objtypes['Lep'][args.lepton], 1, (-2.6, 2.6), True)})



    ## get plotting colors/settings
hstyles = styles.styles

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
    #if (hname not in variables.keys()) and (hname not in variables_2d.keys()): continue
    #if hname == 'cutflow': continue
    histo = hdict[hname]
    #set_trace()

    if histo.dense_dim() == 1:
        xtitle, rebinning, x_lims, withData = variables[hname]

        ## hists should have 3 category axes (dataset, jet multiplicity, lepton type) followed by variable
        for lep in [args.lepton]:
            for jmult in histo.axis('jmult')._sorted:
                pltdir = outdir if args.testing else '/'.join([outdir, lep, jmult])
                if not os.path.isdir(pltdir):
                    os.makedirs(pltdir)

                if withData:
                    fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                else:
                    fig, ax = plt.subplots()
                fig.subplots_adjust(hspace=.07)
                hslice = histo[:, jmult, lep].integrate('jmult').integrate('leptype')

                #if 'Lep_eta' in hname: set_trace()
                if hname == 'Jets_njets':
                    print(jmult)
                    yields_txt, yields_json = get_samples_yield_and_frac(hslice, lep)
                    plt_tools.print_table(yields_txt, filename='%s/%s_%s_yields_and_fracs.txt' % (pltdir, jmult, lep), print_output=True)
                    with open('%s/%s_%s_yields_and_fracs.json' % (pltdir, jmult, lep), 'w') as out:
                        out.write(prettyjson.dumps(yields_json))

                if rebinning != 1:
                    xaxis_name = hslice.dense_axes()[0].name
                    hslice = hslice.rebin(xaxis_name, rebinning)

                    ## plot mc+data
                if withData:
                    ax, rax = Plotter.plot_stack1d(ax, rax, hslice, xlabel=xtitle, xlimits=x_lims)
                else:
                    ax = Plotter.plot_mc1d(ax, hslice, xlabel=xtitle, xlimits=x_lims)

                    ## set axes labels and titles
                    # add lepton/jet multiplicity label
                ax.text(
                    0.02, 0.92, "%s, %s" % (objtypes['Lep'][lep], jet_mults[jmult]),
                    fontsize=hep.styles_cms.CMS['font.size']*0.75, 
                    horizontalalignment='left', 
                    verticalalignment='bottom', 
                    transform=ax.transAxes
                )
                ax = hep.cms.cmslabel(ax=ax, data=True if withData else False, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % lep]/1000., 1))

                figname = '%s/%s' % (pltdir, '_'.join([jmult, lep, hname]))
                fig.savefig(figname)
                print('%s written' % figname)
                plt.close()
                #set_trace()

    elif histo.dense_dim() == 2:
        xtitle, ytitle, x_rebinning, x_lims, y_rebinning, y_lims, withData = variables_2d[hname]

        ## hists should have 3 category axes (dataset, jet multiplicity, lepton type) followed by variable
        for lep in [args.lepton]:
            for jmult in histo.axis('jmult')._sorted:
                pltdir = outdir if args.testing else '/'.join([outdir, lep, jmult])
                if not os.path.isdir(pltdir):
                    os.makedirs(pltdir)

                hslice = histo[:, jmult, lep].integrate('jmult').integrate('leptype')

                if x_rebinning != 1:
                    xaxis_name = hslice.dense_axes()[0].name
                    hslice = hslice.rebin(xaxis_name, x_rebinning)
                if y_rebinning != 1:
                    yaxis_name = hslice.dense_axes()[1].name
                    hslice = hslice.rebin(yaxis_name, y_rebinning)

                    # make 1D projections along y axis
                bin_slices = [(x, x+10) for x in range(0, 100, 10)] # [(0, 10), (10, 20)..]
                for bin_xmin, bin_xmax in bin_slices:
                    if withData:
                        fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                    else:
                        fig, ax = plt.subplots()
                    fig.subplots_adjust(hspace=.07)

                    #set_trace()
                    y_slice = hslice[:, bin_xmin:bin_xmax, :].integrate(hslice.dense_axes()[0].name)

                        ## plot mc+data
                    if withData:
                        ax, rax = Plotter.plot_stack1d(ax, rax, yslice, xlabel=xtitle, xlimits=x_lims)
                    else:
                        ax = Plotter.plot_mc1d(ax, yslice, xlabel=xtitle, xlimits=x_lims)


                        ## set axes labels and titles
                    ax.text(
                        0.02, 0.92, "%s, %s\n%.1f $\leq$ %s $\leq$ %.1f" % (objtypes['Lep'][lep], jet_mults[jmult], bin_xmin, xtitle, bin_xmax),
                        fontsize=hep.styles_cms.CMS['font.size']*0.75, 
                        horizontalalignment='left', 
                        verticalalignment='bottom', 
                        transform=ax.transAxes
                    )
                    ax = hep.cms.cmslabel(ax=ax, data=True if withData else False, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % lep]/1000., 1))

                    #figname = 'test'
                    figname = '%s/%s' % (pltdir, '_'.join([jmult, lep, hname, 'xrange%sto%s' % (bin_xmin, bin_xmax)]))
                    fig.savefig(figname)
                    print('%s written' % figname)
                    plt.close()
                    #set_trace()


                    # make 2D plots for each topology
                for top in hslice.sparse_axes()[0]._sorted:
                    if not withData and top == 'data': continue

                    if withData:
                        fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                    else:
                        fig, ax = plt.subplots()
                    fig.subplots_adjust(hspace=.07)

                    htop = hslice[top].integrate('process')
                    top_label = plt_tools.get_label(top, hstyles)
                    htop.label = '%s %s' % (top_label, htop.label)
                        ## plot MC and data
                    #set_trace()
                    plot.plot2d(htop,
                        xaxis=htop.axes()[0].name,
                        ax=ax,
                        #patch_opts={'cmap' : 'OrRd'},
                        clear=True,
                    )
                    ax.autoscale(axis='x', tight=True)
                    ax.set_ylim(y_lims)
                    ax.set_xlim(x_lims)


                        ## set axes labels and titles
                    plt.xlabel(xtitle)
                    plt.ylabel(ytitle)
                    ax.text(
                        0.02, 0.92, "%s, %s" % (objtypes['Lep'][lep], jet_mults[jmult]),
                        fontsize=hep.styles_cms.CMS['font.size']*0.75, 
                        horizontalalignment='left', 
                        verticalalignment='bottom', 
                        transform=ax.transAxes
                    )
                    ax = hep.cms.cmslabel(ax=ax, data=True if withData else False, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % lep]/1000., 1))

                    figname = '%s/%s' % (pltdir, '_'.join([jmult, lep, hname, top]))
                    fig.savefig(figname)
                    print('%s written' % figname)
                    plt.close()
                    #set_trace()

