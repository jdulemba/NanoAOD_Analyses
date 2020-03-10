# Project-related environment variables
NANODIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P)

source /cvmfs/sft.cern.ch/lcg/views/dev3python3/Tue/x86_64-centos7-gcc8-opt/setup.sh
export LIBRARY_PATH=$LD_LIBRARY_PATH

##unset PYTHONPATH
#if [[ $NANODIR ==  "/uscms_data/d3/jdulemba/NanoAOD_Analyses" ]]; then
#    echo "Adding coffea to PYTHONPATH to run on the LPC normal nodes"
#fi
#
#export PYTHONPATH="$NANODIR/coffea:$PYTHONPATH"
export PYTHONPATH="$NANODIR/Analysis:$PYTHONPATH"
export LD_LIBRARY_PATH=$NANODIR/Analysis/compiled:$LD_LIBRARY_PATH
