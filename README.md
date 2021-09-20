Creating to commit

To run singularity from this directory (on lxplus)

Interactively through shell
```
singularity shell  -B /tmp/x509up_u81826 -B ${PWD}:/srv -B /afs/cern.ch/work/j/jdulemba/Test_Coffea:/scratch --pwd /srv   /cvmfs/unpacked.cern.ch/registry.hub.docker.com/coffeateam/coffea-base:latest   /bin/bash
```

Execute commands
```
singularity exec  -B /tmp/x509up_u81826 -B ${PWD}:/srv -B /afs/cern.ch/work/j/jdulemba/Test_Coffea:/scratch --pwd /srv   /cvmfs/unpacked.cern.ch/registry.hub.docker.com/coffeateam/coffea-base:latest   /bin/bash -c "command(s) to run"
```

In order to set up analysis environment:

```
cd Analysis & source environment.sh

```
