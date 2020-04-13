import numpy as np
from coffea.lookup_tools.dense_lookup import dense_lookup
from coffea.util import load, save
from pdb import set_trace
import os
from coffea import hist

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
analyzer = 'htt_flav_effs_analyzer'

jobid = os.environ['jobid']
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

pt_binning = np.array([30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0, 100.0, 105.0, 110.0, 125.0, 150.0,170.0, 200.0, 250.0, 1000.0])
eta_binning = np.array([-2.5, -1.5, -0.5, 0.0, 0.5, 1.5, 2.5])
#eta_binning = np.array([-2.5, -2., -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5])
pt_bins = hist.Bin('pt', 'pt', pt_binning)
eta_bins = hist.Bin('eta', 'eta', eta_binning)

working_points = []

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

    for wp in hpass_lep.axis('btagger')._sorted: # [DEEPCSVMEDIUM, DEEPJETMEDIUM]
        for flav in hpass_lep.axis('hFlav')._sorted: # [bjet, cjet, ljet]
            for jmult in hpass_lep.axis('jmult')._sorted: #['3Jets', '4PJets']
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

                tagger = 'DEEPCSV' if wp.startswith('DEEPCSV') else 'DEEPJET'
                working_points.append(wp.split(tagger)[-1])
                flav_effs[year][tagger][flav_to_name[flav]].update({jmult : eff_lookup})
   
wp_name = list(set(working_points))[0]
    # save files
flav_effs_name = '%s/htt_3PJets_%s_flavour_efficiencies_%s.coffea' % (outdir, wp_name, jobid)
save(flav_effs, flav_effs_name)
print('\n', flav_effs_name, 'written')

