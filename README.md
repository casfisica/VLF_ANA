# VLF_ANA

## Installation instructions:
```bash
export SCRAM_ARCH=slc6_amd64_gcc530

cmsrel CMSSW_8_0_25

cd CMSSW_8_0_25/src

cmsenv

git cms-init  
git cms-addpkg RecoMET/METProducers  
scram b -j 10  
git cms-merge-topic -u cms-met:fromCMSSW_8_0_20_postICHEPfilter  
scram b -j 10  
#git clone git@github.com:HENGEX/VLF_ANA.git
git clone git@github.com:casfisica/VLF_ANA.git
cd VLF_ANA/
git checkout TagAndSave
cd ..
scram b -j 10

```
## Running with crab:

<par>Create a new crab task and submit the jobs:  </par>
```bash
crab submit -c python/YOUR_CRAB_CONFIG_FILE.py  
#Checking status:  
crab status -d crab_projects/YOUR_CRAB_TASK_DIRECTORY/  
#Resubmitting failed jobs (check in the status output if the jobs can be resubmitted):  
crab resubmit -d crab_projects/YOUR_CRAB_TASK_DIRECTORY/  
```
<par>A good tutorial on crab can be found here: https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookCRAB3Tutorial  </par>


## Runing the submit scripts
<par>To send many samples at the same time you can use the <i>Submit_*.py</i> scripts in the <i>test</i> folder </par>

```bash
./Submit_DYPlusJets.py

```
