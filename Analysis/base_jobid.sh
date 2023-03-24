if [[ $jobid == *"NanoAODv6"* ]]; then
    export base_jobid=NanoAODv6
elif [[ $jobid == *"Summer20UL"* ]]; then
    export base_jobid=Summer20UL
elif [[ $jobid == *"ULnanoAODv9"* ]]; then
    export base_jobid=ULnanoAODv9
#elif [[ $jobid == *"ULnanoAOD"* ]]; then
#    export base_jobid=ULnanoAOD
else
    echo "base_jobid NOT SET"
fi
