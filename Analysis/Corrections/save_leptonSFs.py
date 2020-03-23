import numpy as np
from coffea.util import save
from coffea.lookup_tools.root_converters import convert_histo_root_file
from coffea.lookup_tools.dense_lookup import dense_lookup
from pdb import set_trace
import os

proj_dir = os.environ['PROJECT_DIR']

outdir = '%s/Corrections' % proj_dir
if not os.path.isdir(outdir):
    os.makedirs(outdir)

## pt and eta lists are hardcoded for now
leptons = {
    'Electrons' : {
        'pt' : ['SF_ElTOT_0', 'SF_ElTOT_1', 'SF_ElTOT_2', 'SF_ElTOT_3', 'SF_ElTOT_4', 'SF_ElTOT_5', 'SF_ElTOT_6', 'SF_ElTOT_7', 'SF_ElTOT_8', 'SF_ElTOT_9'],
        'eta' : 'eta_ElTOT',
        'fnames' : {
            '2016' : 'SF_El_V16010b.root',
            '2017' : 'SF_El_V170301b.root',
            '2018' : 'SF_El_V180301b.root',
        }
    },
    'Muons' : {
        'pt' : ['SF_MuTOT_0', 'SF_MuTOT_1', 'SF_MuTOT_2', 'SF_MuTOT_3', 'SF_MuTOT_4', 'SF_MuTOT_5', 'SF_MuTOT_6', 'SF_MuTOT_7'],
        'eta' : 'eta_MuTOT',
        'fnames' : {
            '2016' : 'SF_Mu_V16010.root', 
            '2017' : 'SF_Mu_V170301.root',
            '2018' : 'SF_Mu_V180301.root',
        }
    }
}

sf_output = {
    '2016' : {
        'Electrons' : {
            'Central' : {},
            'Error' : {},
        },
        'Muons' : {
            'Central' : {},
            'Error' : {},
        },
    },
    '2017' : {
        'Electrons' : {
            'Central' : {},
            'Error' : {},
        },
        'Muons' : {
            'Central' : {},
            'Error' : {},
        },
    },
    '2018' : {
        'Electrons' : {
            'Central' : {},
            'Error' : {},
        },
        'Muons' : {
            'Central' : {},
            'Error' : {},
        },
    },
}

for lep in leptons.keys():
    for year, fname in leptons[lep]['fnames'].items():
        sf_file = convert_histo_root_file('%s/inputs/data/%s' % (proj_dir, fname))
        eta_binning = sf_file[(leptons[lep]['eta'], 'dense_lookup')][1]

        #set_trace()
        sf_output[year][lep]['eta_ranges'] = [(eta_binning[idx], eta_binning[idx+1]) for idx in range(len(eta_binning)-1)]
        for idx, pt_hist in enumerate(leptons[lep]['pt']):
            sf_output[year][lep]['Central']['eta_bin%i' % idx] = dense_lookup(*sf_file[(pt_hist, 'dense_lookup')])
            sf_output[year][lep]['Error']['eta_bin%i' % idx] = dense_lookup(*sf_file[('%s_error' % pt_hist, 'dense_lookup')])

lepSF_name = '%s/leptonSFs.coffea' % outdir
save(sf_output, lepSF_name)    
print('%s written' % lepSF_name)
