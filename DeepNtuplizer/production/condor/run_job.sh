#!/bin/bash

inputFile=$1
sample=$2
index=$3

echo "Starting job on `date`" # date/time of start of job
echo "Running on: `uname -a`" # condor job is running on this node
# echo "System software: `cat /etc/redhat-release`" # operating system on that node

source /cvmfs/cms.cern.ch/cmsset_default.sh

export X509_USER_PROXY=/afs/cern.ch/user/e/etsai/.globus/x509up_u125665
export EOS_MGM_URL=root://eosuser.cern.ch
export SCRAM_ARCH=el9_amd64_gcc12

cp /afs/cern.ch/user/e/etsai/workspace/SecondaryVertexing_CMSSW_14_1_0_pre0.tgz ./SecondaryVertexing_CMSSW_14_1_0_pre0.tgz
tar -xf SecondaryVertexing_CMSSW_14_1_0_pre0.tgz
rm SecondaryVertexing_CMSSW_14_1_0_pre0.tgz
cd SecondaryVertexing_CMSSW_14_1_0_pre0/src/
scramv1 b ProjectRename # this handles linking the already compiled code - do NOT recompile
cmsenv
echo $CMSSW_BASE "is the CMSSW on the local worker node"

cd DeepNTuples/DeepNtuplizer/production/
echo "`pwd` is the current directory"
ls -alrth

export HOME=/afs/cern.ch/user/e/etsai
echo Home is $HOME
echo Output directory: /eos/cms/store/group/phys_btag/etsai/output_$sample
echo Input file: $inputFile
cmsRun DeepNtuplizer.py inputFiles=root://cms-xrd-global.cern.ch/$inputFile outputFile=/eos/cms/store/group/phys_btag/etsai/output_$sample/output_$index
