import numpy as np
from coffea.lookup_tools.dense_lookup import dense_lookup
from coffea.util import load, save
from pdb import set_trace
import os
from coffea import hist

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'htt_flav_effs_analyzer'

jobid = 'Testing' ## only temporary!!
outdir = '%s/Corrections/%s' % (proj_dir, jobid)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

flav_effs = {
    '2016' : {
        'DEEPJET' : {
            'bottom' : {},
            'charm' : {},
            'light' : {},
        },
        'DEEPCSV' : {
            'bottom' : {},
            'charm' : {},
            'light' : {},
        },
    },
    '2017' : {
        'DEEPJET' : {
            'bottom' : {},
            'charm' : {},
            'light' : {},
        },
        'DEEPCSV' : {
            'bottom' : {},
            'charm' : {},
            'light' : {},
        },
    },
    '2018' : {
        'DEEPJET' : {
            'bottom' : {},
            'charm' : {},
            'light' : {},
        },
        'DEEPCSV' : {
            'bottom' : {},
            'charm' : {},
            'light' : {},
        },
    },
}

flav_to_name = {'bjet' : 'bottom', 'cjet' : 'charm', 'ljet' : 'light'}
hname = 'Jets_pt_eta'
lumi_correction = load('%s/Corrections/%s/MC_LumiWeights.coffea' % (proj_dir, jobid))

pt_bins_3j = np.array([20., 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 125.0, 150.0,170.0, 200.0, 250.0, 1000.0])
eta_bins_3j = np.array([-2.5, -1.8, -1.2, -0.6, 0.0, 0.6, 1.2, 1.8, 2.5])
pt_3j = hist.Bin('pt', 'pt', pt_bins_3j)
eta_3j = hist.Bin('eta', 'eta', eta_bins_3j)

for year in ['2016', '2017', '2018']:
    input_dir = '/'.join([proj_dir, 'results', '%s_%s' % (year, jobid), analyzer])
    fnames = ['%s/%s' % (input_dir, fname) for fname in os.listdir(input_dir) if fname.endswith('TOT.coffea')]
    if len(fnames) > 1:
        raise ValueError("Only one TOT file should be used")
    hdict = load(fnames[0])
    h_pass = hdict['%s_pass' % hname]
    h_all = hdict['%s_all' % hname]

    if not h_pass.compatible(h_all):
        raise ValueError("Passing and All hists don't have the same binning!")

        ## rescale hist by lumi for muons and electrons separately and then combine
    h_pass_mu = h_pass[:, :, :, 'Muon', :].integrate('leptype')
    h_pass_mu.scale(lumi_correction[year]['Muons'], axis='dataset')
    h_pass_el = h_pass[:, :, :, 'Electron', :].integrate('leptype')
    h_pass_el.scale(lumi_correction[year]['Electrons'], axis='dataset')
    h_pass_lep = h_pass_mu+h_pass_el

    h_all_mu = h_all[:, :, :, 'Muon', :].integrate('leptype')
    h_all_mu.scale(lumi_correction[year]['Muons'], axis='dataset')
    h_all_el = h_all[:, :, :, 'Electron', :].integrate('leptype')
    h_all_el.scale(lumi_correction[year]['Electrons'], axis='dataset')
    h_all_lep = h_all_mu+h_all_el

    for wp in h_pass_lep.axis('btagger')._sorted: # [DEEPCSVMEDIUM, DEEPJETMEDIUM]
        for flav in h_pass_lep.axis('hFlav')._sorted: # [bjet, cjet, ljet]
            for jmult in h_pass_lep.axis('jmult')._sorted: #['3Jets', '4PJets']
                set_trace()
                    # get passing and all hists for 3 and 4+ jets separately, only as a function of pT and eta
                h_pass = h_pass_lep[wp, :, jmult, flav].integrate('btagger').integrate('dataset').integrate('jmult').integrate('hFlav')
                h_pass = h_pass.rebin('pt', pt_3j).rebin('eta', 10)
                #h_pass = h_pass.rebin('pt', pt_3j).rebin('eta', eta_3j)
                h_all = h_all_lep[wp, :, jmult, flav].integrate('btagger').integrate('dataset').integrate('jmult').integrate('hFlav')
                h_all = h_all.rebin('pt', pt_3j).rebin('eta', 10)
                #h_all = h_all.rebin('pt', pt_3j).rebin('eta', eta_3j)

                edges = (h_pass.axis('pt').edges(), h_pass.axis('eta').edges())
                pass_lookup = dense_lookup(*(h_pass.values().values()), edges)
                    ## check that number of entries in passing bins > 20.
                if min([min(val) for val in pass_lookup._values]) < 20.:
                    raise ValueError("bin for hist %s, %s, %s has bin entries < 20: %s" % (wp, flav, jmult, [min(val) for val in pass_lookup._values]))

                all_lookup = dense_lookup(*(h_all.values().values()), edges)
                eff_lookup = dense_lookup(pass_lookup._values/all_lookup._values, edges)

                tagger = 'DEEPCSV' if wp.startswith('DEEPCSV') else 'DEEPJET'
                flav_effs[year][tagger][flav_to_name[flav]].update({jmult : eff_lookup})
   
        
    # save files
mcweights_name = '%s/MC_PU_Weights.coffea' % outdir
save(mc_pu_weights, mcweights_name)
print('\n', mcweights_name, 'written')

data_pu_name = '%s/data_PU_dists.coffea' % outdir
save(data_pu_dists, data_pu_name)
print('\n', data_pu_name, 'written')

