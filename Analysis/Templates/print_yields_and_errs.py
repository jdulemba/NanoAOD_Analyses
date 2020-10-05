#! /bin/env python

from coffea.util import load#, save
from pdb import set_trace
import os
import Utilities.prettyjson as prettyjson
import numpy as np
#import fnmatch
#import Utilities.systematics as systematics

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('year', choices=['2016', '2017', '2018'], help='What year is the ntuple from.')
parser.add_argument('--njets', default='all', nargs='?', choices=['3', '4+', 'all'], help='Specify which jet multiplicity to use.')
parser.add_argument('--only_bkg', action='store_true', help='Make background templates only.')
parser.add_argument('--only_sig', action='store_true', help='Make signal templates only.')

args = parser.parse_args()

proj_dir = os.environ['PROJECT_DIR']
jobid = os.environ['jobid']

njets_to_run = []
if (args.njets == '3') or (args.njets == 'all'):
    njets_to_run += ['3Jets']
if (args.njets == '4+') or (args.njets == 'all'):
    njets_to_run += ['4PJets']


input_dir = os.path.join(proj_dir, 'Templates', 'results', jobid, args.year)

fname_3j = os.path.join(input_dir, 'raw_templates_lj_3Jets_bkg_%s_%s.coffea' % (jobid, args.year))
if not os.path.isfile(fname_3j):
    raise ValueError("File %s with bkg templates for %s not found" % (fname_3j, args.year))
fname_4pj = os.path.join(input_dir, 'raw_templates_lj_4PJets_bkg_%s_%s.coffea' % (jobid, args.year))
if not os.path.isfile(fname_4pj):
    raise ValueError("File %s with bkg templates for %s not found" % (fname_4pj, args.year))

hdict_3j = load(fname_3j)
mu_dict_3j = hdict_3j['Muon']
el_dict_3j = hdict_3j['Electron']

hdict_4pj = load(fname_4pj)
mu_dict_4pj = hdict_4pj['Muon']
el_dict_4pj = hdict_4pj['Electron']

nosys_names = [key for key in mu_dict_3j.keys() if 'nosys' in key]

proc_to_names = {
    'TT_nosys' : "\\ttbar",
    'QCD_nosys' : 'QCD',
    'WJets_nosys' : 'W+jets',
    'ZJets_nosys' : "Z/$\gamma^{*}$+jets",
    'TTV_nosys' : "{\\ttbar}V",
    'VV_nosys' : 'VV',
    'sChannel_nosys' : 'Single top',
    'tChannel_nosys' : 'Single top',
    'tWChannel_nosys' : 'Single top',
    'data_obs_nosys' : 'Data'
}

yields_out = "& \multicolumn{2}{c}{Muon Channel} & \multicolumn{2}{c}{Electron Channel} \\\ \n\hline\hline \n"
yields_out += "Process & 3 jets & 4+ jets & 3 jets & 4+ jets \\\ \n\hline \n"

yields = {}
sumw2s = {}

for name, proc in proc_to_names.items():
    mu_3j_histo = mu_dict_3j[name]
    el_3j_histo = el_dict_3j[name]

    mu_3j_sumw, mu_3j_sumw2 = mu_3j_histo.values(sumw2=True)[()]
    mu_3j_yield = mu_3j_sumw.sum()
    mu_3j_err = np.sqrt(mu_3j_sumw2.sum())
    el_3j_sumw, el_3j_sumw2 = el_3j_histo.values(sumw2=True)[()]
    el_3j_yield = el_3j_sumw.sum()
    el_3j_err = np.sqrt(el_3j_sumw2.sum())

    mu_4pj_histo = mu_dict_4pj[name]
    el_4pj_histo = el_dict_4pj[name]

    mu_4pj_sumw, mu_4pj_sumw2 = mu_4pj_histo.values(sumw2=True)[()]
    mu_4pj_yield = mu_4pj_sumw.sum()
    mu_4pj_err = np.sqrt(mu_4pj_sumw2.sum())
    el_4pj_sumw, el_4pj_sumw2 = el_4pj_histo.values(sumw2=True)[()]
    el_4pj_yield = el_4pj_sumw.sum()
    el_4pj_err = np.sqrt(el_4pj_sumw2.sum())

    yields[name] = {'m3' : mu_3j_yield, 'm4' : mu_4pj_yield, 'e3' : el_3j_yield, 'e4' : el_4pj_yield}
    sumw2s[name] = {'m3' : mu_3j_sumw2.sum(), 'm4' : mu_4pj_sumw2.sum(), 'e3' : el_3j_sumw2.sum(), 'e4' : el_4pj_sumw2.sum()}
    if proc != 'Single top':
        yields_out += "{PROC} & {MU_3J_YIELD:.1f} $\pm$ {MU_3J_ERR:.1f} & {MU_4PJ_YIELD:.1f} $\pm$ {MU_4PJ_ERR:.1f} & {EL_3J_YIELD:.1f} $\pm$ {EL_3J_ERR:.1f} & {EL_4PJ_YIELD:.1f} $\pm$ {EL_4PJ_ERR:.1f} \\\\ \n".format(PROC=proc, MU_3J_YIELD=mu_3j_yield, MU_3J_ERR=mu_3j_err, MU_4PJ_YIELD=mu_4pj_yield, MU_4PJ_ERR=mu_4pj_err, EL_3J_YIELD=el_3j_yield, EL_3J_ERR=el_3j_err, EL_4PJ_YIELD=el_4pj_yield, EL_4PJ_ERR=el_4pj_err)



    # get single top yields from all processes
singlet_yields = {key:(yields['sChannel_nosys'][key]+yields['tChannel_nosys'][key]+yields['tWChannel_nosys'][key]) for key in yields['sChannel_nosys'].keys()}
singlet_errs = {key:np.sqrt(sumw2s['sChannel_nosys'][key]+sumw2s['tChannel_nosys'][key]+sumw2s['tWChannel_nosys'][key]) for key in yields['sChannel_nosys'].keys()}

yields_out += "{PROC} & {MU_3J_YIELD:.1f} $\pm$ {MU_3J_ERR:.1f} & {MU_4PJ_YIELD:.1f} $\pm$ {MU_4PJ_ERR:.1f} & {EL_3J_YIELD:.1f} $\pm$ {EL_3J_ERR:.1f} & {EL_4PJ_YIELD:.1f} $\pm$ {EL_4PJ_ERR:.1f} \\\\ \n".format(PROC="Single top", MU_3J_YIELD=singlet_yields['m3'], MU_3J_ERR=singlet_errs['m3'], MU_4PJ_YIELD=singlet_yields['m4'], MU_4PJ_ERR=singlet_errs['m4'], EL_3J_YIELD=singlet_yields['e3'], EL_3J_ERR=singlet_errs['e3'], EL_4PJ_YIELD=singlet_yields['e4'], EL_4PJ_ERR=singlet_errs['e4'])

    # get total background yields
bkg_yields = {chan:sum([yields[key][chan] for key in yields.keys() if not 'data_obs' in key]) for chan in yields['TT_nosys'].keys()}
bkg_errs = {chan:np.sqrt(sum([sumw2s[key][chan] for key in sumw2s.keys() if not 'data_obs' in key])) for chan in sumw2s['TT_nosys'].keys()}

yields_out += "Total Background & {MU_3J_YIELD:.1f} $\pm$ {MU_3J_ERR:.1f} & {MU_4PJ_YIELD:.1f} $\pm$ {MU_4PJ_ERR:.1f} & {EL_3J_YIELD:.1f} $\pm$ {EL_3J_ERR:.1f} & {EL_4PJ_YIELD:.1f} $\pm$ {EL_4PJ_ERR:.1f} \\\\ \n".format(MU_3J_YIELD=bkg_yields['m3'], MU_3J_ERR=bkg_errs['m3'], MU_4PJ_YIELD=bkg_yields['m4'], MU_4PJ_ERR=bkg_errs['m4'], EL_3J_YIELD=bkg_yields['e3'], EL_3J_ERR=bkg_errs['e3'], EL_4PJ_YIELD=bkg_yields['e4'], EL_4PJ_ERR=bkg_errs['e4'])


print(yields_out)
obs_yields = open('%s/yields_and_errs_%s.txt' % (input_dir, args.year), 'w')
obs_yields.write(yields_out)
obs_yields.close()
print('%s/yields_and_errs_%s.txt written' % (input_dir, args.year))
