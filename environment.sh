# Project-related environment variables
export NANODIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P)

source /cvmfs/sft.cern.ch/lcg/views/LCG_97python3/x86_64-centos7-gcc8-opt/setup.sh
export LIBRARY_PATH=$LD_LIBRARY_PATH

export PYTHONPATH="$NANODIR/Analysis:$PYTHONPATH"
export PYTHONPATH=/afs/cern.ch/user/j/jdulemba/.local/lib/python3.7/site-packages:$PYTHONPATH
export PYTHONPATH="$NANODIR/coffea:$PYTHONPATH"
export LD_LIBRARY_PATH=$NANODIR/Analysis/compiled:$LD_LIBRARY_PATH
export PATH=/afs/cern.ch/user/j/jdulemba/.local/bin:$PATH

python --version
