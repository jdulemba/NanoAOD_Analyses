# Project-related environment variables
NANODIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd -P)

#unset PYTHONPATH
if [[ $NANODIR ==  "/uscms_data/d3/jdulemba/NanoAOD_Analyses" ]]; then
    echo "Adding coffea to PYTHONPATH to run on the LPC normal nodes"
    export PYTHONPATH="$NANODIR/coffea:$PYTHONPATH"
fi

export PYTHONPATH="$NANODIR/Analysis:$PYTHONPATH"
export LD_LIBRARY_PATH=$NANODIR/Analysis/compiled
