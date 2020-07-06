import numpy as np
from coffea.lookup_tools.dense_lookup import dense_lookup
from coffea.util import load, save
from pdb import set_trace
import os
from coffea import hist
from coffea.hist import plot
# matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
plt.style.use(hep.cms.style.ROOT)
plt.switch_backend('agg')
from matplotlib import rcParams
rcParams['font.size'] = 18
rcParams["savefig.format"] = 'png'
rcParams["savefig.bbox"] = 'tight'
import Utilities.styles as styles
import Utilities.plot_tools as plt_tools
import Utilities.prettyjson as prettyjson

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('--construct_btag', action='store_false', help='Makes btag SF constructor (default is True)')
args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'htt_flav_effs_analyzer'

jobid = os.environ['jobid']
outdir = '%s/Corrections/%s' % (proj_dir, jobid)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

flav_effs = {
    '2016' : {
        'DeepJet' : {
            '3Jets' : {},
            '4PJets' : {},
        },
        'DeepCSV' : {
            '3Jets' : {},
            '4PJets' : {},
        },
    },
    '2017' : {
        'DeepJet' : {
            '3Jets' : {},
            '4PJets' : {},
        },
        'DeepCSV' : {
            '3Jets' : {},
            '4PJets' : {},
        },
    },
    '2018' : {
        'DeepJet' : {
            '3Jets' : {},
            '4PJets' : {},
        },
        'DeepCSV' : {
            '3Jets' : {},
            '4PJets' : {},
        },
    },
}

if args.construct_btag:
    from copy import deepcopy
    btag_contructs_dict = deepcopy(flav_effs)

jet_mults = {
    '3Jets' : '3 jets',
    '4PJets' : '4+ jets',
}

flav_to_name = {'bjet' : 'bottom', 'cjet' : 'charm', 'ljet' : 'light'}
hname = 'Jets_pt_eta'
lumi_correction = load('%s/Corrections/%s/MC_LumiWeights_IgnoreSigEvts.coffea' % (proj_dir, jobid))

#pt_binning = np.array([30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 125.0, 150.0,170.0, 200.0, 250.0, 1000.0])
#eta_binning = np.array([-2.5, -1.5, -0.5, 0.0, 0.5, 1.5, 2.5])
pt_binning = np.array([30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 125.0, 150.0,170.0, 200.0, 1000.0])
eta_binning = np.array([-2.5, -2., -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5])
pt_bins = hist.Bin('pt', 'pt', pt_binning)
eta_bins = hist.Bin('eta', 'eta', eta_binning)

working_points = []

def plot_effs(heff, edges, lumi_to_use, year, jmult, btagger, wp, flav, plotdir, clear=True):
    fig, ax = plt.subplots()
    fig.subplots_adjust(hspace=.07)

    opts = {'cmap' : 'OrRd'}
    xedges, yedges = edges[0], edges[1]
    sumw = heff._values
    pc = ax.pcolormesh(xedges, yedges, sumw.T, **opts)
    ax.add_collection(pc)
    if clear:
        fig.colorbar(pc, ax=ax, label='%s Efficiency, %s %s' % (flav_to_name[flav], btagger, wp))
    ax.autoscale(axis='x', tight=True)
    ax.set_xlim(xedges[0], xedges[-1])
    ax.set_ylim(yedges[0], yedges[-1])
    
        ## set axes labels and titles
    plt.xlabel('$p_{T}$ [GeV]')
    plt.ylabel('$\\eta$')
        # add lepton/jet multiplicity label
    ax.text(
        0.02, 0.95, "%s" % (jet_mults[jmult]),
        fontsize=18,
        horizontalalignment='left',
        verticalalignment='bottom',
        transform=ax.transAxes
    )
    ax = hep.cms.cmslabel(ax=ax, fontsize=18, data=False, paper=False, year=year, lumi=round(lumi_to_use, 1))

    #set_trace()
    #figname = 'test.png'    
    figname = '%s/%s_Efficiency' % (plotdir, '_'.join([btagger, wp, jmult, flav]))
    fig.savefig(figname)
    print('%s written' % figname)
    plt.close()
    #set_trace()


data_lumi_dict = prettyjson.loads(open('%s/inputs/lumis_data.json' % proj_dir).read())

for year in ['2016', '2017', '2018']:
    input_dir = '/'.join([proj_dir, 'results', '%s_%s' % (year, jobid), analyzer])
    fnames = ['%s/%s' % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith('TOT.coffea')]
    if len(fnames) > 1:
        raise ValueError("Only one TOT file should be used")
    hdict = load(fnames[0])
    hpass = hdict['%s_pass' % hname]
    hall = hdict['%s_all' % hname]

    if not hpass.compatible(hall):
        raise ValueError("Passing and All hists don't have the same binning!")

        ## rescale hist by lumi for muons and electrons separately and then combine
    hpass_mu = hpass[:, :, :, 'Muon', :].integrate('leptype')
    hpass_mu.scale(lumi_correction[year]['Muons'], axis='dataset')
    hpass_el = hpass[:, :, :, 'Electron', :].integrate('leptype')
    hpass_el.scale(lumi_correction[year]['Electrons'], axis='dataset')
    hpass_lep = hpass_mu+hpass_el

    hall_mu = hall[:, :, :, 'Muon', :].integrate('leptype')
    hall_mu.scale(lumi_correction[year]['Muons'], axis='dataset')
    hall_el = hall[:, :, :, 'Electron', :].integrate('leptype')
    hall_el.scale(lumi_correction[year]['Electrons'], axis='dataset')
    hall_lep = hall_mu+hall_el

        ## get data lumi and scale MC by lumi
    lumi_to_use = (data_lumi_dict[year]['Muons']+data_lumi_dict[year]['Electrons'])/2000.

    pltdir = '/'.join([proj_dir, 'plots', '%s_%s' % (year, jobid), analyzer])
    if not os.path.isdir(pltdir):
        os.makedirs(pltdir)

    for wp in hpass_lep.axis('btagger')._sorted: # [DEEPCSVMEDIUM, DEEPJETMEDIUM]
        for jmult in hpass_lep.axis('jmult')._sorted: #['3Jets', '4PJets']
            for flav in hpass_lep.axis('hFlav')._sorted: # [bjet, cjet, ljet]
                #set_trace()
                    # get passing and all hists for 3 and 4+ jets separately, only as a function of pT and eta
                h_pass = hpass_lep[wp, :, jmult, flav].integrate('btagger').integrate('dataset').integrate('jmult').integrate('hFlav')
                h_pass = h_pass.rebin('pt', pt_bins).rebin('eta', eta_bins)
                h_all = hall_lep[wp, :, jmult, flav].integrate('btagger').integrate('dataset').integrate('jmult').integrate('hFlav')
                h_all = h_all.rebin('pt', pt_bins).rebin('eta', eta_bins)

                edges = (h_pass.axis('pt').edges(), h_pass.axis('eta').edges())
                pass_lookup = dense_lookup(*(h_pass.values().values()), edges)
                    ## check that number of entries in passing bins > 20.
                if min([min(val) for val in pass_lookup._values]) < 20.:
                    print("bin for hist %s, %s, %s has bin entries < 20: %s" % (wp, flav, jmult, [min(val) for val in pass_lookup._values]))
                    #raise ValueError("bin for hist %s, %s, %s has bin entries < 20: %s" % (wp, flav, jmult, [min(val) for val in pass_lookup._values]))

                all_lookup = dense_lookup(*(h_all.values().values()), edges)
                eff_lookup = dense_lookup(pass_lookup._values/all_lookup._values, edges)

                tagger = 'DeepCSV' if wp.upper().startswith('DEEPCSV') else 'DeepJet'
                working_points.append(wp.upper().split(tagger.upper())[-1])
                flav_effs[year][tagger][jmult].update({flav_to_name[flav] : eff_lookup})

                plot_effs(eff_lookup, edges, lumi_to_use, year, jmult, tagger, wp.upper().split(tagger.upper())[-1][0], flav, pltdir)

#set_trace()
wp_name = list(set(working_points))[0]
    # save files
flav_effs_name = '%s/htt_3PJets_%s_flavour_efficiencies_%s.coffea' % (outdir, wp_name, jobid)
save(flav_effs, flav_effs_name)
print('\n', flav_effs_name, 'written')


if args.construct_btag:
    import Utilities.prettyjson as prettyjson
    import python.BTagScaleFactors as btagSF

    cfg_file = prettyjson.loads(open('%s/cfg_files/cfg_pars_%s.json' % (proj_dir, jobid)).read())

    for year in flav_effs.keys():
        for btagger in flav_effs[year].keys():
            csv_path = '/'.join([proj_dir, 'inputs', 'data', btagSF.btag_csvFiles[year][btagger]])
            if not os.path.isfile(csv_path):
                raise IOError('BTagging csv file %s not found.' % csv_path)

            for njets_cat in flav_effs[year][btagger].keys():
                eff_dict = flav_effs[year][btagger][njets_cat]
                sf_computer = btagSF.BTagSF(
                    csv = csv_path,
                    wp_key = (btagger, 'used', wp_name.lower().capitalize()),
                    effs = eff_dict
                )

                print('BTag SF constructed for %s, %s %s wp, %s' % (year, btagger, wp_name.lower().capitalize(), njets_cat))
                btag_contructs_dict[year][btagger][njets_cat].update({wp_name.lower().capitalize() : sf_computer})

        # save files
    btagSFs_name = '%s/htt_3PJets_%s_btag_scalefactors_%s.coffea' % (outdir, wp_name, jobid)
    save(btag_contructs_dict, btagSFs_name)
    print('\n', btagSFs_name, 'written')
