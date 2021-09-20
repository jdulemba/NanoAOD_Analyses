# Project-related environment variables
whereIam=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P)

#This part should be changed by the user(s)
export jobid=NOTSET
if [ -e $whereIam/jobid.sh ]
then
        source $whereIam/jobid.sh
        echo "set jobid: $jobid"
        echo "set nanoAOD version: $base_jobid"
else
        echo "I did not find jobid.sh, are you sure you do not want to set the jobid and leave it to $jobid?"
fi

#HERE ARE LIONS!
export PROJECT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P)

export LIBRARY_PATH=$LD_LIBRARY_PATH

export PYTHONPATH="$PROJECT_DIR:$PYTHONPATH"
export LD_LIBRARY_PATH=$PROJECT_DIR/compiled:$LD_LIBRARY_PATH

python --version
