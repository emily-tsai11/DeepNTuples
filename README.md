# DeepNTuples
NTuple framework for DeepFlavour

Installation CMSSW 11.1.2 patch3 with new track timing variables, PV3D+timing
============
```
cmsrel CMSSW_11_1_2_patch3
cd CMSSW_11_1_2_patch3/src/
cmsenv
git cms-init
git cms-merge-topic mneukum:11_1_2_patch3_dnt_pv3d
git clone https://github.com/mneukum/DeepNTuples
cd DeepNTuples
git checkout timing_variables_pv3d
# Add JetToolBox
git submodule init
git submodule update
# copy files explicitly. Read additional_files/info for more information
git cms-addpkg DataFormats/BTauReco
. additional_files/copy_files.sh
# The following line is necessary for a consistent setup, since TaggingVariables.h/cc require recompilation, esp. for track related variables! Fetch a coffee while compiling
git cms-checkdeps -a
cd $CMSSW_BASE/src
scram b -j 10
```

Installation CMSSW 11.3.0 pre3, PV3D+timing, for RelVal
============
```
cmsrel CMSSW_11_3_0_pre3
cd CMSSW_11_3_0_pre3/src/
cmsenv
git cms-init
git cms-merge-topic mneukum:CMSSW_11_3_0_pre3_dnt_pv3d
git clone https://github.com/mneukum/DeepNTuples
cd DeepNTuples
git checkout timing_variables_pv3d
# Add JetToolBox
git submodule init
git submodule update
# copy files explicitly. Read additional_files/info for more information
git cms-addpkg DataFormats/BTauReco
. additional_files/copy_files.sh
# The following line is necessary for a consistent setup, since TaggingVariables.h/cc require recompilation, esp. for track related variables! Fetch a coffee while compiling
git cms-checkdeps -a
cd $CMSSW_BASE/src
scram b -j 10
```

Caveats
============

- Compilation of a lot of CMSSW modules is necessary, because of changed TaggingVariable.h/cc for time related variables.
- In 11_3_0_pre3: differences in muon_isHighPt and muon_energy. These variables are not used in the training.
- For TagVarCSV_vertexCategory==1: set TagVarCSV_flightDistance* to 0. Script will be uploaded soon and should be included in DeepNTuplizer code later. This is an old bug to recover the same values as in BTagAnalyzer (only few jets).


Further settings
============

It is important to create your grid proxy in a location that is accessible by other nodes (there is no security issue, your full credentials are still needed for access). For this purpose, redirect the grid proxy location by adding the following to your login script:

```
export X509_USER_PROXY=${HOME}/.gridproxy.pem
```

Production
==========

Before doing a batch submission you can test the ntuplizer locally in the production directory with:
```
cmsRun DeepNtuplizer.py inputFiles=/path/to/file.root
```
The jobs can be submitted using the following syntax
```
jobSub.py --file <sample file> DeepNtuplizer.py <batch directory> --outpath /path/to/output/directory/
```
For an example of sample files, please refer to the .cfg files already in the production directory. You first specify the number of jobs to be submitted, then the input dataset name, which should then be followed by the name of the output. Other arguments such as gluonReduction can then be specified if needed. Each argument need to be separted by at least two whitespaces.
 
The large job output (root files) will NOT be stored in the batch directory. The storage directory is specified by the --outpath argument. The batch directory will contain a symlink to this directory. If the outpath is not specified the ntuples are stored in the deepjet directory, where you need write permission.

The status of the jobs can be checked with
```
cd <batch directory>
check.py <sample subdirectories to be checked>
```

The check.py script provides additional options to resubmit failed jobs or to create sample lists in case a satisfying fraction of jobs ended successfully. 
In this case do:
```
check.py <sample subdirectories to be checked> --action filelist
```
This will create file lists that can be further processed by the DeepJet framework
For resubmitting failed jobs, do:
```
check.py <sample subdirectories to be checked> --action resubmit
```

When the file lists are created, the part used for training of the ttbar and QCD samples (or in principle any other process) can be merged using the executable:
```
mergeSamples.py <no of jets per file> <output dir> <file lists 1> <file lists 2> <file lists 3> ...
```
For example:
```
mergeSamples.py 400000 /path/to/dir/merged ntuple_*/train_val_samples.txt
```
This will take a significant amount of time - likely more than the ntuple production itself. It is therefore recommended to run the command within 'screen'. In the 94X branch you can also submit via batch by doing --batch. This will create a batch directory in the folder the command is called from.

```
mergeSamples.py 400000 /path/to/dir/merged ntuple_*/train_val_samples.txt --batch
```
