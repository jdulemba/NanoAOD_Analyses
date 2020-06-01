from coffea.hist import plot
# matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend('agg')
from matplotlib import rcParams
rcParams['font.size'] = 20
from coffea.util import load
from pdb import set_trace
import os
import styles
import Utilities.plot_tools as plt_tools
import Utilities.prettyjson as prettyjson
import re
from coffea import hist
import numpy as np
from equal_split import partition_list

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('year', choices=['2016', '2017', '2018'], help='What year is the ntuple from.')
parser.add_argument('lepton', choices=['Electron', 'Muon'], help='Choose which lepton to make plots for')

args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'evt_cuts_invest'

input_dir = '/'.join([proj_dir, 'results', '%s_%s' % (args.year, jobid), analyzer])
f_ext = 'TOT.coffea'
outdir = '/'.join([proj_dir, 'plots', '%s_%s' % (args.year, jobid), analyzer])
if not os.path.isdir(outdir):
    os.makedirs(outdir)

fnames = sorted(['%s/%s' % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith(f_ext)])

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

lep_cats = {
    #'Loose_or_Tight' : 'loose or tight %s' % objtypes['Lep'][args.lepton],
    'Tight' : 'tight %s' % objtypes['Lep'][args.lepton],
    'Loose' : 'loose %s' % objtypes['Lep'][args.lepton],
}

variables = {
    'Jets_LeadJet_pt' : ('$p_{T}$(leading jet) [GeV]', 2, (0., 300.), True),
    'Jets_njets' : ('$n_{jets}$', 1, (0, 15), True),
    'Lep_iso' : ('pfRelIso, %s' % objtypes['Lep'][args.lepton], 1, (0., 5.) if args.lepton == 'Muon' else (0., 1.5), True),
    'MT' : ('$M_{T}$', 1, (0., 300.), True),
    'BTagSF' : ('$SF_{btag}$', 1, (0.7, 1.5), False),
    'LepSF' : ('$SF_{lep}$', 1, (0.8, 1.1), False),
    'EvtWeight' : ('Event Weight', 1, (0., 2.), False),
    #'MT_Iso' : ('$M_{T}$', 'pfRelIso, %s' % objtypes['Lep'][args.lepton], 1, (0., 300.), 1, (0., 1.), True),
    #'MT_LJpt' : ('$M_{T}$', '$p_{T}$(leading jet) [GeV]', 1, (0., 300.), 2, (0., 300.), True),
    #'Iso_LJpt' : ('pfRelIso, %s' % objtypes['Lep'][args.lepton], '$p_{T}$(leading jet) [GeV]', 1, (0., 1.), 2, (0., 300.), True),
}
if args.lepton == 'Electron':
    variables.update({
        'El_iso_barrel' : ('pass/fail tight Iso, |$\eta_{SC}$| $\leq$ 1.479', 1, (0., 2.), True),
        'El_iso_endcap' : ('pass/fail tight Iso, |$\eta_{SC}$| $>$ 1.479', 1, (0., 2.), True),
    })


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


def cumulative_fraction(data, mc, lep, topos_to_use=[]):
    mc_mask = re.compile(r'\b(?:%s)\b' % '|'.join(topos_to_use))
    data_vals = np.sum(np.array(list(data.values().values())), axis=0)
    all_mc_vals = np.sum(np.array(list(mc.values().values())), axis=0)
    tt_vals = np.sum(np.array(list(mc['ttJets'].values().values())), axis=0)
    mc_to_use_vals = np.sum(np.array(list(mc[mc_mask].values().values())), axis=0)

    edges = mc.dense_axes()[0].edges()
    effs = [0.95, 0.9, 0.85, 0.8, 0.75, 0.7]#, 0.65, 0.6, 0.55, 0.5]
    #set_trace()
    eff_bins = [np.min(np.where(mc_to_use_vals/data_vals < eff)[0]) for eff in effs]
    iso_vals = edges[eff_bins]

    rows = [tuple(['Prompt MC/data: ']+['%.2f' % eff for eff in effs])]
    rows += [tuple(['Iso: ']+['%.2f' % iso for iso in iso_vals])]

    return rows


    ## make plots
for hname in variables.keys():
    if hname not in hdict.keys():
        raise ValueError("%s not found in file" % hname)
    #set_trace()
    histo = hdict[hname][:, :, args.lepton, :, :].integrate('leptype') # process, jmult, lepton, lepcat, btagregion

    if histo.dense_dim() == 1:
        xtitle, rebinning, x_lims, withData = variables[hname]

        #set_trace()
        ## hists should have 3 category axes (dataset, jet multiplicity, lepton category) followed by variable
        for jmult in list(set([key[1] for key in histo.values().keys()])):
            for btagregion in list(set([key[3] for key in histo.values().keys()])):
                for lepcat in list(set([key[2] for key in histo.values().keys()])):
                    pltdir = '/'.join([outdir, args.lepton, jmult, lepcat, btagregion])
                    if not os.path.isdir(pltdir):
                        os.makedirs(pltdir)

                    if withData:
                        fig, (ax, rax) = plt.subplots(2, 1, gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
                    else:
                        fig, ax = plt.subplots()
                    fig.subplots_adjust(hspace=.07)
                    hslice = histo[:, jmult, lepcat, btagregion].integrate('jmult').integrate('lepcat').integrate('btag')

                    #if 'Lep_eta' in hname: set_trace()
                    if hname == 'Jets_njets':
                        print(jmult)
                        yields_txt, yields_json = get_samples_yield_and_frac(hslice, args.lepton)
                        plt_tools.print_table(yields_txt, filename='%s/%s_%s_yields_and_fracs.txt' % (pltdir, jmult, args.lepton), print_output=True)
                        with open('%s/%s_%s_yields_and_fracs.json' % (pltdir, jmult, args.lepton), 'w') as out:
                            out.write(prettyjson.dumps(yields_json))

                    if rebinning != 1:
                        xaxis_name = hslice.dense_axes()[0].name
                        hslice = hslice.rebin(xaxis_name, rebinning)

                    #set_trace()
                    #mc_dict = {key:np.sum(list(hslice[key].values().values()), axis=1) for key in [key[0] for key in list(hslice[mc_samples].values().keys())]} # make {topo: sum events} dict for sorting
                    #mc_order = sorted(mc_dict, key=mc_dict.get, reverse=False)
                        ## plot MC and data
                    plot.plot1d(hslice[mc_samples],
                        overlay=hslice.axes()[0].name,
                        ax=ax,
                        clear=False,
                        stack=True,
                        line_opts=None,
                        fill_opts=stack_fill_opts,
                        error_opts=stack_error_opts,
                        order=['QCD', 'EWK', 'singlet', 'ttJets'],
                        #order=mc_order,
                    )
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
                            denom_fill_opts={},
                            guide_opts={},
                            unc='num'
                        )
                        rax.set_ylabel('data/MC')
                        rax.set_ylim(0.5, 1.5)
                        rax.set_xlim(x_lims)


                        ## set axes labels and titles
                    plt.xlabel(xtitle)
                        # add lepton/jet multiplicity label
                    #set_trace()
                    ax.text(
                        0.02, 0.92 if withData else 0.95, "%s, %s" % (lep_cats[lepcat], jet_mults[jmult]),
                        fontsize=rcParams['font.size']*0.9, 
                        horizontalalignment='left', 
                        verticalalignment='bottom', 
                        transform=ax.transAxes
                    )
                    ax = hep.cms.cmslabel(ax=ax, data=withData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
                    if hname == 'Lep_iso':
                        #print(jmult)
                        #prompt_mc_frac = cumulative_fraction(data=hslice[data_samples], mc=hslice[mc_samples], lep=args.lepton, topos_to_use=['ttJets', 'singlet','EWK'])
                        #plt_tools.print_table(prompt_mc_frac, filename='%s/%s_%s_Iso_Fracs.txt' % (pltdir, jmult, args.lepton), print_output=True)

                        mc_vals = np.sum(np.array(list(hslice[mc_samples].values().values())), axis=0)
                        data_vals = np.sum(np.array(list(hslice[data_samples].values().values())), axis=0)
                        max_binval = np.maximum(np.max(mc_vals), np.max(data_vals))
                        min_binval = np.minimum(np.min(mc_vals), np.min(data_vals))
                        ax.set_yscale('log')
                        ax.set_ylim(np.maximum(1e-1, min_binval), 10**(np.ceil(np.log10(max_binval))+1))
                        if lepcat == 'Tight':
                            ax.set_xlim((0., 0.2))

                    figname = '%s/%s.png' % (pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, hname]))
                    fig.savefig(figname, bbox_inches='tight')
                    print('%s written' % figname)
                    #set_trace()
                    plt.close()

    


## save values used to determine isolation sidebands (only for Loose leptons in data)
lepIso_invest = False
#lepIso_invest = True
if lepIso_invest:
    pltdir = '/'.join([outdir, args.lepton])
    if not os.path.isdir(pltdir):
        os.makedirs(pltdir)

    histo = hdict['Lep_iso'][:, :, args.lepton, 'Loose', :].integrate('leptype').integrate('lepcat')[data_samples].integrate('process')
    xtitle, rebinning, x_lims, withData = variables['Lep_iso']

    if rebinning != 1:
        xaxis_name = histo.dense_axes()[0].name
        hslice = histo.rebin(xaxis_name, rebinning)
        
    iso_vals = histo.axis('iso').edges(overflow='all')

    vals_to_save = {
        args.lepton : {'3Jets' : {}, '4PJets' : {}, '3PJets' : {}}
    }    
    for btag in ['btagPass', 'btagFail', 'btagComb']:
        histo = histo.integrate('btag') if btag == 'btagComb' else histo[btag].integrate('btag')
        fig, ax = plt.subplots()
        fig.subplots_adjust(hspace=.07)
        for idx, jmult in enumerate(['3Jets', '4PJets', '3PJets']):
            hslice = histo.integrate('jmult') if jmult == '3PJets' else histo[jmult].integrate('jmult')

            #set_trace()        
            ## construct iso sideband regions to investigate based on data in 'Loose' category
            data_array = hslice.values(sumw2=False, overflow='all')[()]
            data_partitions, inds_for_splitting = partition_list(data_array, 3) # get splitting of data and inds where that splitting occurs: returns list for partitions, ndarray for inds
            vals_to_save[args.lepton][jmult] = {
                'iso_binning' : iso_vals.tolist(),
                'orig_data_array' : data_array.tolist(),
                'partition_loc_inds' : inds_for_splitting.tolist(),
                'iso_region_edges' : iso_vals[inds_for_splitting].tolist(),
                'data_partitions' : data_partitions,
                'data_partition_sums' : list(map(sum, data_partitions)),
            }

            plot.plot1d(hslice,
                ax=ax,
                clear=False,
                line_opts={'linestyle' : hstyles[jmult]['linestyle'], 'color' : hstyles[jmult]['color'], 'linewidth' : hstyles[jmult]['elinewidth']},
            )
        ax.autoscale(axis='x', tight=True)
        ax.set_ylim(0, None)
        ax.set_xlabel(None)
        ax.set_xlim(x_lims)
        
            ## set legend labels
        handles, labels = ax.get_legend_handles_labels()
        labels[0] = hstyles['3Jets']['name']
        labels[1] = hstyles['4PJets']['name']
        labels[2] = hstyles['3PJets']['name']
        # call ax.legend() with the new values
        ax.legend(handles,labels, loc='upper right')
        
        #set_trace()
        
            ## set axes labels and titles
        plt.xlabel(xtitle)
            # add lepton/jet multiplicity label
        ax.text(
            0.02, 0.95, lep_cats['Loose'],
            fontsize=hep.styles_cms.CMS['font.size']*0.90, 
            horizontalalignment='left', 
            verticalalignment='bottom', 
            transform=ax.transAxes
        )
        ax = hep.cms.cmslabel(ax=ax, data=True if withData else False, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))
        data_vals = np.array(list(histo.values().values()))
        max_binval = np.max(data_vals)
        min_binval = np.min(data_vals[data_vals > 0])
        ax.set_yscale('log')
        ax.set_ylim(np.maximum(1e-1, min_binval), 10**np.ceil(np.log10(max_binval)))
        
        figname = '%s/%s.png' % (pltdir, '_'.join(['Loose', args.lepton, 'data', 'Lep_iso']))
        fig.savefig(figname, bbox_inches='tight')
        print('%s written' % figname)
        #set_trace()
        plt.close()
        #set_trace()


    outjson_name = '%s/%s_iso_sideband_regions_construction.json' % (pltdir, '_'.join(['Loose', args.lepton]))
    with open(outjson_name, 'w') as out:
        out.write(prettyjson.dumps(vals_to_save))
    print('%s written' % outjson_name)
    
