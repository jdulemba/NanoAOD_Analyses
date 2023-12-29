Creating to commit

To run singularity from this directory (on lxplus)

Interactively through shell
```
singularity shell  -B /tmp/x509up_u81826 -B ${PWD}:/srv -B /afs/cern.ch/work/j/jdulemba/Test_Coffea:/scratch -B /eos/user/j/jdulemba/NanoAOD_Analyses --pwd /srv   /cvmfs/unpacked.cern.ch/registry.hub.docker.com/coffeateam/coffea-base:0.7.20-fastjet-3.4.0.1   /bin/bash
```

Execute commands
```
singularity exec  -B /tmp/x509up_u81826 -B ${PWD}:/srv -B /afs/cern.ch/work/j/jdulemba/Test_Coffea:/scratch --pwd /srv -B /eos/user/j/jdulemba/NanoAOD_Analyses  /cvmfs/unpacked.cern.ch/registry.hub.docker.com/coffeateam/coffea-base:0.7.20-fastjet-3.4.0.1   /bin/bash -c "command(s) to run"
```

In order to set up analysis environment:

```
cd Analysis & source environment.sh

```
