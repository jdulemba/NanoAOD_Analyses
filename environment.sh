# Project-related environment variables
export NANODIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P)

source /cvmfs/sft.cern.ch/lcg/views/dev3python3/Tue/x86_64-centos7-gcc8-opt/setup.sh
#source /cvmfs/sft.cern.ch/lcg/views/dev3python3/Thu/x86_64-centos7-gcc8-opt/setup.sh
export LIBRARY_PATH=$LD_LIBRARY_PATH

export PYTHONPATH="$NANODIR/coffea:$PYTHONPATH"
export PYTHONPATH="$NANODIR/Analysis:$PYTHONPATH"
export LD_LIBRARY_PATH=$NANODIR/Analysis/compiled:$LD_LIBRARY_PATH
