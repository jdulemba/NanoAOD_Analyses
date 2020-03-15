from pdb import set_trace
import os
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('year', choices=['2016', '2017', '2018'], help='Specify which year to run over')
parser.add_argument('--group', type=str, help='Run samples from specific group')
parser.add_argument('--sample', type=str, help='Use specific sample')

args = parser.parse_args()

if args.group and args.sample:
    raise IOError("Either group or sample can be specified, not both.")

#proj_dir = os.environ['PROJECT_DIR']
#jobid = os.environ['jobid']


## define sample groups
QCD_EM = [
    "QCD_EM_20to30",
    "QCD_EM_30to50",
    "QCD_EM_50to80",
    "QCD_EM_80to120",
    "QCD_EM_120to170",
    "QCD_EM_170to300",
    "QCD_EM_300toInf",
]
if args.year == '2017': QCD_EM.append('QCD_EM_15to20')

QCD_Mu = [
    "QCD_Mu_15to20",
    "QCD_Mu_20to30",
    "QCD_Mu_30to50",
    "QCD_Mu_50to80",
    "QCD_Mu_80to120",
    "QCD_Mu_120to170",
    "QCD_Mu_170to300",
    "QCD_Mu_300to470",
    "QCD_Mu_470to600",
    "QCD_Mu_600to800",
    "QCD_Mu_800to1000",
    "QCD_Mu_1000toInf",
]

ttV = [
    "ttWlnu",
    "ttWqq",
    "ttZll",
    "ttZqq",
]

diboson = [
    "WW",
    "WZ",
    "ZZ",
]

WJets = [
    "W1Jets",
    "W2Jets",
    "W3Jets",
    "W4Jets",
    "WJets",
]

ZJets = [
    "ZJets",
]

if args.year == '2016':
    singlet = [
        "singlet_schannel_PS",
        "singlet_tW_PS",
        "singletbar_tW_PS",
        "singlet_tchannel_PS",
        "singletbar_tchannel_PS"
    ]
else:
    singlet = [
        "singlet_schannel",
        "singlet_tchannel",
        "singlet_tW",
        "singletbar_tW",
        "singletbar_tchannel",
    ]

ttSL_sys = [
    "ttJetsSL_hdampDOWN",
    "ttJetsSL_hdampUP",
    "ttJetsSL_mtopDOWN",
    "ttJetsSL_mtopUP",
    "ttJetsSL_mtop1695",
    "ttJetsSL_mtop1755",
    "ttJetsSL_ueDOWN",
    "ttJetsSL_ueUP",
]

ttHad_sys = [
    "ttJetsHad_hdampDOWN",
    "ttJetsHad_hdampUP",
    "ttJetsHad_mtopDOWN",
    "ttJetsHad_mtopUP",
    "ttJetsHad_mtop1695",
    "ttJetsHad_mtop1755",
    "ttJetsHad_ueDOWN",
    "ttJetsHad_ueUP",
]

ttDiLep_sys = [
    "ttJetsDiLep_hdampDOWN",
    "ttJetsDiLep_hdampUP",
    "ttJetsDiLep_mtopDOWN",
    "ttJetsDiLep_mtopUP",
    "ttJetsDiLep_mtop1695",
    "ttJetsDiLep_mtop1755",
    "ttJetsDiLep_ueDOWN",
    "ttJetsDiLep_ueUP",
]

tt_sys = [
    "ttJets_fsrDOWN",
    "ttJets_fsrUP",
    "ttJets_hdampDOWN",
    "ttJets_hdampUP",
    "ttJets_isrDOWN",
    "ttJets_isrUP",
    "ttJets_mtopDOWN",
    "ttJets_mtopUP",
    "ttJets_mtop1695",
    "ttJets_mtop1755",
    "ttJets_ueDOWN",
    "ttJets_ueUP",
]

if args.year == '2016':
    data_el = [
        'data_SingleElectron_2016Bv2',
        'data_SingleElectron_2016C',
        'data_SingleElectron_2016D',
        'data_SingleElectron_2016E',
        'data_SingleElectron_2016F',
        'data_SingleElectron_2016G',
        'data_SingleElectron_2016H'
    ]
    data_mu = [
        'data_SingleMuon_2016Bv2',
        'data_SingleMuon_2016C',
        'data_SingleMuon_2016D',
        'data_SingleMuon_2016E',
        'data_SingleMuon_2016F',
        'data_SingleMuon_2016G',
        'data_SingleMuon_2016H'
    ]

elif args.year == '2017':
    data_el = [
        'data_SingleElectron_2017B',
        'data_SingleElectron_2017C',
        'data_SingleElectron_2017D',
        'data_SingleElectron_2017E',
        'data_SingleElectron_2017F',
    ]
    data_mu = [
        'data_SingleMuon_2017B',
        'data_SingleMuon_2017C',
        'data_SingleMuon_2017D',
        'data_SingleMuon_2017E',
        'data_SingleMuon_2017F',
    ]
else:
    data_el = [
        'data_SingleElectron_2018A',
        'data_SingleElectron_2018B',
        'data_SingleElectron_2018C',
        'data_SingleElectron_2018D',
    ]
    data_mu = [
        'data_SingleMuon_2018A',
        'data_SingleMuon_2018B',
        'data_SingleMuon_2018C',
        'data_SingleMuon_2018D',
    ]


groups_dict = {
    'QCD' : QCD_EM+QCD_Mu, 'QCD_EM' : QCD_EM, 'QCD_Mu' : QCD_Mu,
    'ttV' : ttV, 'diboson' : diboson,
    'singlet' : singlet,
    'WJets' : WJets, 'ZJets' : ZJets, 'EWK' : WJets+ZJets,
    'tt_sys' : tt_sys, 'ttJets' : ["ttJets"], 'ttJets_PS' : ["ttJets_PS"],
    'ttSL_sys' : ttSL_sys, 'ttHad_sys' : ttHad_sys, 'ttDiLep_sys' : ttDiLep_sys,
    'ttJetsSL' : ["ttJetsSL"], 'ttJetsHad' : ["ttJetsHad"], 'ttJetsDiLep' : ["ttJetsDiLep"],
    'data_el' : data_el, 'data_mu' : data_mu, 'data' : data_el+data_mu,
}

if args.year == '2016':
    if args.group in ['ttJetsSL', 'ttJetsHad', 'ttJetsDiLep', 'ttSL_sys', 'ttHad_sys', 'ttDiLep_sys']:
        raise IOError("%s is not a valid group for running 2016" % args.group)

#set_trace()
if args.group:
    for sample in groups_dict[args.group]:
        run_cmd = "python bin/test_analyzer.py all {YEAR} --sample={SAMPLE}".format(YEAR=args.year, SAMPLE=sample)
        print("Executing: %s" % run_cmd)
        os.system(run_cmd)

if args.sample:
    run_cmd = "python bin/test_analyzer.py all {YEAR} --sample={SAMPLE}".format(YEAR=args.year, SAMPLE=args.sample)
    print("Executing: %s" % run_cmd)
    os.system(run_cmd)

