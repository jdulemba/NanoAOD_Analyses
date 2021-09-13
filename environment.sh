# Project-related environment variables
export NANODIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P)

export LIBRARY_PATH=$LD_LIBRARY_PATH

export PYTHONPATH="$NANODIR/Analysis:$PYTHONPATH"
export PYTHONPATH="$NANODIR/coffea:$PYTHONPATH"
export LD_LIBRARY_PATH=$NANODIR/Analysis/compiled:$LD_LIBRARY_PATH

python --version
