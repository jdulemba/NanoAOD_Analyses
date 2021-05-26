```
LUMIMASK Golden Jsons
The lumimasks for each year were taken from the following links.

2016: https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2016Analysis under Re-reco datasets 07Aug17

2017 & 2018: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TWikiLUM, 'Golden JSON' links
```


```
PILEUP

find latest pileup file at https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions[YEAR(16/17/18)]/13TeV/PileUp/pileup_latest.txt

Create pileup files for each year by running make_nanoAOD_data_pileup.sh from CMSSW_9_4_10/src/Analyses/URTTbar/data_scripts (as of March 6, 2020)

```

```
Lepton Scale Factors:

Files used are from Otto. Which files to use are found in his cfg files (/uscms_data/d3/ohindric/RunAnalyzer/config_1[678]UL.cfg)
and the files are then copied from /uscms_data/d1/ohindric/RunAnalyzer/INUL1[678]/INPUT/
```
