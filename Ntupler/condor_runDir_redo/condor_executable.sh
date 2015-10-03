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
tar -zxf configFiles.tgz

# Set up variables to select files to run cmsRun, ZCAnalyzer over
let jobNum=$1+1
fileName=`sed -n "$jobNum p" < configFileList.txt`

# Run the main process
./Ntupler $fileName

# Clean up
rm Ntupler
rm Cert_*
rm btag_generic.txt
rm ntuple_job_*.py
