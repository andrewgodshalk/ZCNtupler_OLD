#!/bin/bash
# 2014-04-28 - Written as wrapper script for CONDOR jobs to use
# VHbbNtupler on VHbb PATTuple files to make ZcNtuples.
#
# Last modified: 2015-01-16
#
# Variables:
# $1: job number
#
# Things to change in each job:
# Original run directory (CMSSW Release directory where you compiled Ntupler)
# fileName <- list of file locations for Ntupler to run over.
#    - many located in formattedFileLists folder in tar included in condor job.

# Set up env and unpack all files
pwd
#source /uscmst1/prod/sw/cms/bashrc prod
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd ${_CONDOR_SCRATCH_DIR}
pwd 
# original run directory
cd /uscms_data/d2/godshalk/2015-01_ZHbb_cp/CMSSW_5_3_23_patch1/src
pwd
eval `scramv1 runtime -sh`
cd ${_CONDOR_SCRATCH_DIR}
pwd
tar -zxf zcntuplerFiles.tgz


# Set up variables to select files to run cmsRun, ZCAnalyzer over
let jobNum=$1
let filesPerJob=10
let firstFile=$filesPerJob*$jobNum+1
let lastFile=$firstFile+$filesPerJob-1
fileName="formattedFileLists/dy.txt"

echo ""
echo "files from $fileName"
echo "job number $jobNum"
echo "files/job  $filesPerJob"
echo "from $firstFile to $lastFile"
echo ""


# Set up configuration file for running the set of files it needs.
sed -n "1,9 p" <ntuple_for_condor.py > new_config.py
#sed -n "$firstFile,$lastFile p" <muonA.txt >> new_config.py
sed -n "$firstFile,$lastFile p" <$fileName >> new_config.py
sed -n "10,47 p" <ntuple_for_condor.py >> new_config.py
echo "fname = 'Test' + channel + '_$1.root'" >> new_config.py
sed -n "49,400 p" <ntuple_for_condor.py >> new_config.py
mv new_config.py ntuple_for_condor.py

# Print a mess of information to CMSSW_#.stdout
echo ""
echo ""
echo "== PROCESSING FILES :"
#echo `sed -n "$firstFile,$lastFile p" <muonA.txt`
echo `sed -n "$firstFile,$lastFile p" <$fileName`
echo ""
echo ""
 
# Run the main process
./Ntupler ntuple_for_condor.py

# Clean up

rm Ntupler
rm Cert_190456-203002_8TeV_PromptReco_Collisions12_JSON.txt
rm btag_generic.txt
mv ntuple_for_condor.py ntuple_job_$1.py
