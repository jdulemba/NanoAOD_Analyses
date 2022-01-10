#!/bin/bash

# Project-related environment variables
export PROJECT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P)

string_to_check="jdulemba"
if [[ $PROJECT_DIR == *"$string_to_check"* ]]
then
    source /cvmfs/sft.cern.ch/lcg/views/LCG_97python3/x86_64-centos7-gcc8-opt/setup.sh
fi 

#This part should be changed by the user(s)
export jobid=NOTSET
if [ -e $PROJECT_DIR/jobid.sh ]
then
        source $PROJECT_DIR/jobid.sh
        echo "set jobid: $jobid"
        echo "set nanoAOD version: $base_jobid"
else
        echo "I did not find jobid.sh, are you sure you do not want to set the jobid and leave it to $jobid?"
fi

#HERE ARE LIONS!
export LIBRARY_PATH=$LD_LIBRARY_PATH

export PYTHONPATH="$PROJECT_DIR:$PYTHONPATH"
export LD_LIBRARY_PATH=$PROJECT_DIR/compiled:$LD_LIBRARY_PATH

python --version
