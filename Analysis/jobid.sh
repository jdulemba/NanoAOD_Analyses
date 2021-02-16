export jobid=NanoAODv6
#export jobid=NanoAODv6_jpt30_ljpt50_MT40_cutBasedEl
#export jobid=NanoAODv6_newJCs
#export jobid=NanoAODv6_pt_thad_reweighting
#export jobid=NanoAODv6_mtt_thadctstar_reweighting
#export jobid=NanoAODv6_mtt_thadctstar_Interp_reweighting

#export jobid=ULnanoAOD
#export jobid=ULnanoAOD_jpt30_ljpt50_MT40_cutBasedEl

if [[ $jobid == *"NanoAODv6"* ]]; then
    export base_jobid=NanoAODv6
else
    export base_jobid=ULnanoAOD
fi
