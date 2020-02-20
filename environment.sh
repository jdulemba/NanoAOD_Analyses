# Project-related environment variables
NANODIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P)

unset PYTHONPATH
export PYTHONPATH="$NANODIR/coffea:$PYTHONPATH"
export PYTHONPATH="$NANODIR/Analysis:$PYTHONPATH"
