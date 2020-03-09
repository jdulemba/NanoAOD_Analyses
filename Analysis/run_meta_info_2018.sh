#!/bin/bash

process=$1

declare -a QCD=(
"QCD_EM_20to30"
"QCD_EM_30to50"
"QCD_EM_50to80"
"QCD_EM_80to120"
"QCD_EM_120to170"
"QCD_EM_170to300"
"QCD_EM_300toInf"
"QCD_Mu_15to20"
"QCD_Mu_20to30"
"QCD_Mu_30to50"
"QCD_Mu_50to80"
"QCD_Mu_80to120"
"QCD_Mu_120to170"
"QCD_Mu_170to300"
"QCD_Mu_300to470"
"QCD_Mu_470to600"
"QCD_Mu_600to800"
"QCD_Mu_800to1000"
"QCD_Mu_1000toInf"
)

declare -a ttV=(
"ttWlnu"
"ttWqq"
"ttZll"
"ttZqq"
)

declare -a diboson=(
"WW"
"WZ"
"ZZ"
)

declare -a singlet=(
"singlet_schannel"
"singlet_tW"
"singlet_tchannel"
"singletbar_tW"
"singletbar_tchannel"
)

declare -a WJets=(
"W1Jets"
"W2Jets"
"W3Jets"
"W4Jets"
"WJets"
)

declare -a ttSL_sys=(
"ttJetsSL_hdampDOWN"
"ttJetsSL_hdampUP"
"ttJetsSL_mtopDOWN"
"ttJetsSL_mtopUP"
"ttJetsSL_mtop1695"
"ttJetsSL_mtop1755"
"ttJetsSL_ueDOWN"
"ttJetsSL_ueUP"
)

declare -a ttDiLep_sys=(
"ttJetsDiLep_hdampDOWN"
"ttJetsDiLep_hdampUP"
"ttJetsDiLep_mtopDOWN"
"ttJetsDiLep_mtopUP"
"ttJetsDiLep_mtop1695"
"ttJetsDiLep_mtop1755"
"ttJetsDiLep_ueDOWN"
"ttJetsDiLep_ueUP"
)

declare -a ttHad_sys=(
"ttJetsHad_hdampDOWN"
"ttJetsHad_hdampUP"
"ttJetsHad_mtopDOWN"
"ttJetsHad_mtopUP"
"ttJetsHad_mtop1695"
"ttJetsHad_mtop1755"
"ttJetsHad_ueDOWN"
"ttJetsHad_ueUP"
)

declare -a data_el=(
"data_SingleElectron_2018A"
"data_SingleElectron_2018B"
"data_SingleElectron_2018C"
"data_SingleElectron_2018D"
)

declare -a data_mu=(
"data_SingleMuon_2018A"
"data_SingleMuon_2018B"
"data_SingleMuon_2018C"
"data_SingleMuon_2018D"
)

if [ $process == "QCD" ]; then
    samples=${QCD[@]}
elif [ $process == "ttV" ]; then
    samples=${ttV[@]}
elif [ $process == "diboson" ]; then
    samples=${diboson[@]}
elif [ $process == "singlet" ]; then
    samples=${singlet[@]}
elif [ $process == "WJets" ]; then
    samples=${WJets[@]}
elif [ $process == "ttSL_sys" ]; then
    samples=${ttSL_sys[@]}
elif [ $process == "ttDiLep_sys" ]; then
    samples=${ttDiLep_sys[@]}
elif [ $process == "ttHad_sys" ]; then
    samples=${ttHad_sys[@]}
elif [ $process == "data_el" ]; then
    samples=${data_el[@]}
elif [ $process == "data_mu" ]; then
    samples=${data_mu[@]}
else
    echo "Current choices are QCD, ttV, diboson, singlet, or WJets"
    exit
fi


for sample in ${samples[@]}; do
    echo "Running meta info for $sample"
    python bin/meta_info.py all 2018 --sample=$sample
done
