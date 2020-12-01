import numpy as np
from coffea.util import save
from coffea.lookup_tools.root_converters import convert_histo_root_file
from coffea.lookup_tools.dense_lookup import dense_lookup
from pdb import set_trace
import os

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']
base_jobid = jobid.split('_')[0]

outdir = os.path.join(proj_dir, 'Corrections', jobid)
if not os.path.isdir(outdir):
    os.makedirs(outdir)

## pt and eta lists are hardcoded for now
leptons = {
    'Electrons' : {
        'pt' : {
            'reco_id' : ['SF_ElReco_ElTOT_0', 'SF_ElReco_ElTOT_1', 'SF_ElReco_ElTOT_2', 'SF_ElReco_ElTOT_3', 'SF_ElReco_ElTOT_4', 'SF_ElReco_ElTOT_5', 'SF_ElReco_ElTOT_6', 'SF_ElReco_ElTOT_7', 'SF_ElReco_ElTOT_8', 'SF_ElReco_ElTOT_9'],
            'trig' : ['SF_ElISOTRG_0', 'SF_ElISOTRG_1', 'SF_ElISOTRG_2', 'SF_ElISOTRG_3', 'SF_ElISOTRG_4', 'SF_ElISOTRG_5', 'SF_ElISOTRG_6', 'SF_ElISOTRG_7', 'SF_ElISOTRG_8', 'SF_ElISOTRG_9'],
        },
        'eta' : 'eta_ElTOT',
        'fnames' : {
            #'2016' : 'SF_El_V16010b.root',
            '2017' : 'SF_El_V170400b.root',
            '2018' : 'SF_El_V180400b.root',
        }
    },
    'Muons' : {
        'pt' : {
            'reco_id' : ['SF_MuTRK_MuTOT_0', 'SF_MuTRK_MuTOT_1', 'SF_MuTRK_MuTOT_2', 'SF_MuTRK_MuTOT_3', 'SF_MuTRK_MuTOT_4', 'SF_MuTRK_MuTOT_5', 'SF_MuTRK_MuTOT_6', 'SF_MuTRK_MuTOT_7'],
            'trig' : ['SF_MuISOTRG_0', 'SF_MuISOTRG_1', 'SF_MuISOTRG_2', 'SF_MuISOTRG_3', 'SF_MuISOTRG_4', 'SF_MuISOTRG_5', 'SF_MuISOTRG_6', 'SF_MuISOTRG_7'],
        },
        'eta' : 'eta_MuTOT',
        'fnames' : {
            #'2016' : 'SF_Mu_V16010.root', 
            '2017' : 'SF_Mu_V170400.root',
            '2018' : 'SF_Mu_V180400.root',
        }
    }
}

sf_output = {
    #'2016' : {
    #    'Electrons' : {
    #        'Reco_ID' : {
    #            'Central' : {},
    #            'Error' : {},
    #        },
    #        'Trig' : {
    #            'Central' : {},
    #            'Error' : {},
    #        },
    #    },
    #    'Muons' : {
    #        'Reco_ID' : {
    #            'Central' : {},
    #            'Error' : {},
    #        },
    #        'Trig' : {
    #            'Central' : {},
    #            'Error' : {},
    #        },
    #    },
    #},
    '2017' : {
        'Electrons' : {
            'Reco_ID' : {
                'Central' : {},
                'Error' : {},
            },
            'Trig' : {
                'Central' : {},
                'Error' : {},
            },
        },
        'Muons' : {
            'Reco_ID' : {
                'Central' : {},
                'Error' : {},
            },
            'Trig' : {
                'Central' : {},
                'Error' : {},
            },
        },
    },
    '2018' : {
        'Electrons' : {
            'Reco_ID' : {
                'Central' : {},
                'Error' : {},
            },
            'Trig' : {
                'Central' : {},
                'Error' : {},
            },
        },
        'Muons' : {
            'Reco_ID' : {
                'Central' : {},
                'Error' : {},
            },
            'Trig' : {
                'Central' : {},
                'Error' : {},
            },
        },
    },
}

for lep in leptons.keys():
    for year, fname in leptons[lep]['fnames'].items():
        sf_file = convert_histo_root_file(os.path.join(proj_dir, 'inputs', 'data', base_jobid, 'lepSFs', fname))
        eta_binning = sf_file[(leptons[lep]['eta'], 'dense_lookup')][1]

        #if lep == 'Muons': set_trace()
        sf_output[year][lep]['eta_ranges'] = [(eta_binning[idx], eta_binning[idx+1]) for idx in range(len(eta_binning)-1)]
        for idx in range(len(eta_binning)-1):
                # reco/ID SFs
            sf_output[year][lep]['Reco_ID']['Central']['eta_bin%i' % idx] = dense_lookup(*sf_file[(leptons[lep]['pt']['reco_id'][idx], 'dense_lookup')])
            sf_output[year][lep]['Reco_ID']['Error']['eta_bin%i' % idx] = dense_lookup(*sf_file[('%s_error' % leptons[lep]['pt']['reco_id'][idx], 'dense_lookup')])
                # trigger SFs
            sf_output[year][lep]['Trig']['Central']['eta_bin%i' % idx] = dense_lookup(*sf_file[(leptons[lep]['pt']['trig'][idx], 'dense_lookup')])
            sf_output[year][lep]['Trig']['Error']['eta_bin%i' % idx] = dense_lookup(*sf_file[('%s_error' % leptons[lep]['pt']['trig'][idx], 'dense_lookup')])

lepSF_name = os.path.join(outdir, 'leptonSFs.coffea')
save(sf_output, lepSF_name)    
print('%s written' % lepSF_name)
