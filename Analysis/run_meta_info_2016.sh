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
"singlet_schannel_PS"
"singlet_tW"
"singlet_tW_PS"
"singletbar_tW"
"singletbar_tW_PS"
"singlet_tchannel"
"singlet_tchannel_PS"
"singletbar_tchannel"
"singletbar_tchannel_PS"
)

declare -a WJets=(
"W1Jets"
"W2Jets"
"W4Jets"
"WJets"
)

declare -a tt_sys=(
"ttJets_fsrDOWN"
"ttJets_fsrUP"
"ttJets_hdampDOWN"
"ttJets_hdampUP"
"ttJets_isrDOWN"
"ttJets_isrUP"
"ttJets_mtopDOWN"
"ttJets_mtopUP"
"ttJets_mtop1695"
"ttJets_mtop1755"
"ttJets_ueDOWN"
"ttJets_ueUP"
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
elif [ $process == "tt_sys" ]; then
    samples=${tt_sys[@]}
else
    echo "Current choices are QCD, ttV, diboson, singlet, or WJets"
    exit
fi


for sample in ${samples[@]}; do
    echo "Running meta info for $sample"
    python bin/meta_info.py all 2016 --sample=$sample
done
