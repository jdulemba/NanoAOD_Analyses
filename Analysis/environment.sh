# Project-related environment variables
whereIam=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P)

#This part should be changed by the user(s)
export jobid=NOTSET
if [ -e $whereIam/jobid.sh ]
then
        source $whereIam/jobid.sh
        echo "set jobid: $jobid"
        source $whereIam/nanoTag.sh
        echo "set nanoAOD version: $base_jobid"
else
        echo "I did not find jobid.sh, are you sure you do not want to set the jobid and leave it to $jobid?"
fi
export URA_NTHREADS=1
#export URA_PROJECT_LIBS='-lwhatever_you_need_to_make_it_compile' #<-- Add here needed libraries for the project compilation

#HERE ARE LIONS!
#This part should be handled automatically by the scripts,
#touch it carefully
export PROJECT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P)

#export PYTHONPATH=$PROJECT_DIR
