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

import Plotter as Plotter

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('year', choices=['2016', '2017', '2018'], help='What year is the ntuple from.')
parser.add_argument('lepton', choices=['Electron', 'Muon'], help='Choose which lepton to make plots for')

args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'btagging_invest'

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

btag_cats = {
    'zero_btags' : '0 btags',
    'one_btag' : '1 btag',
    'two_btags' : '2 btags',
    'threePlus_btags' : '3+ btags',
}

lep_cats = {
    'Tight' : 'tight %s' % objtypes['Lep'][args.lepton],
    'Loose' : 'loose %s' % objtypes['Lep'][args.lepton],
}

variables = {
    'Jets_LeadJet_pt' : ('$p_{T}$(leading jet) [GeV]', 2, (0., 300.), True),
    'Jets_pt' : ('$p_{T}$(jets) [GeV]', 2, (0., 300.), True),
    'Jets_njets' : ('$n_{jets}$', 1, (0, 15), True),
    'Jets_DeepCSV' : ('DeepCSV disc (all jets)', 1, (0, 1), True),
    'Lep_pt' : ('$p_{T}$(%s) [GeV]' % objtypes['Lep'][args.lepton], 2, (0., 300.), True),
    'MT' : ('$M_{T}$', 1, (0., 300.), True),
    'nTrueInt_puweight' : ('Reweighted Pileup', 1, (0., 100.), False),
    'nTrueInt_noweight' : ('Unweighted Pileup', 1, (0., 100.), False),
    'rho_puweight' : ('Reweighted $\\rho$', 1, (0., 100.), True),
    'rho_noweight' : ('Unweighted $\\rho$', 1, (0., 100.), True),
    'nvtx_puweight' : ('Reweighted n vertices', 1, (0., 100.), True),
    'nvtx_noweight' : ('Unweighted n vertices', 1, (0., 100.), True),
    'BTagSF' : ('$SF_{btag}$', 1, (0.7, 1.5), False),
    'LepSF' : ('$SF_{lep}$', 1, (0.8, 1.1), False),
    'PileupWeight' : ('Pileup Weight', 1, (0., 2.), False),
    'EvtWeight' : ('Event Weight', 1, (0., 2.), False),
}


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

                    if hname == 'Jets_njets':
                        print(jmult)
                        yields_txt, yields_json = get_samples_yield_and_frac(hslice, args.lepton)
                        plt_tools.print_table(yields_txt, filename='%s/%s_%s_yields_and_fracs.txt' % (pltdir, jmult, args.lepton), print_output=True)
                        with open('%s/%s_%s_yields_and_fracs.json' % (pltdir, jmult, args.lepton), 'w') as out:
                            out.write(prettyjson.dumps(yields_json))

                    if rebinning != 1:
                        xaxis_name = hslice.dense_axes()[0].name
                        hslice = hslice.rebin(xaxis_name, rebinning)

                    mc_opts = {
                        'mcorder' : ['QCD', 'EWK', 'singlet', 'ttJets']
                    }
                    if withData:
                        ax, rax = Plotter.plot_stack1d(ax, rax, hslice, xlabel=xtitle, xlimits=x_lims, **mc_opts)
                    else:
                        ax = Plotter.plot_mc1d(ax, hslice, xlabel=xtitle, xlimits=x_lims, **mc_opts)

                        # add lepton/jet multiplicity label
                    #set_trace()
                    ax.text(
                        0.02, 0.88 if withData else 0.90, "%s, %s\n%s" % (lep_cats[lepcat], jet_mults[jmult], btag_cats[btagregion]),
                        fontsize=rcParams['font.size']*0.75, 
                        horizontalalignment='left', 
                        verticalalignment='bottom', 
                        transform=ax.transAxes
                    )
                    ax = hep.cms.cmslabel(ax=ax, data=withData, paper=False, year=args.year, lumi=round(data_lumi_year['%ss' % args.lepton]/1000., 1))

                    #set_trace()
                    figname = '%s/%s.png' % (pltdir, '_'.join([jmult, args.lepton, lepcat, btagregion, hname]))
                    fig.savefig(figname, bbox_inches='tight')
                    print('%s written' % figname)
                    plt.close()

    


