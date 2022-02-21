#export jobid=NanoAODv6
#export jobid=NanoAODv6_jpt30_ljpt50_MT40_cutBasedEl
#export jobid=NanoAODv6_pt_thad_reweighting
#export jobid=NanoAODv6_mtt_thadctstar_reweighting
#export jobid=NanoAODv6_mtt_thadctstar_Interp_reweighting
#export jobid=NanoAODv6_signal_testing

#export jobid=ULnanoAOD
#export jobid=ULnanoAOD_jpt30_ljpt50_MT40_cutBasedEl
#export jobid=ULnanoAOD_noMT_noLJpt
#export jobid=ULnanoAOD_noMT_noLJpt_NNLO
#export jobid=ULnanoAOD_mujets_btagSFs

#export jobid=Summer20UL
#export jobid=Summer20UL_mujets_btagSFs
#export jobid=Summer20UL_OttoEWCorr
#export jobid=Summer20UL_OttoEWCorr_newBTagSFs
#export jobid=Summer20UL_regroupedJECs
#export jobid=Summer20UL_new_lepSFs
export jobid=Summer20UL_POG_lepSFs

if [[ $jobid == *"NanoAODv6"* ]]; then
    export base_jobid=NanoAODv6
elif [[ $jobid == *"Summer20UL"* ]]; then
    export base_jobid=Summer20UL
elif [[ $jobid == *"ULnanoAOD"* ]]; then
    export base_jobid=ULnanoAOD
else
    echo "base_jobid NOT SET"
fi
