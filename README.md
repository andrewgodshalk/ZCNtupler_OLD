```
VBHF
====
Code for CMS SMP study of Z+c vs Z+b and Z+j XS's.
Developed from code written by CMS VHbb Group.
Git repository (andrewgodshalk/VBHF) created 2014-12-10

Andrew Godshalk
godshalk@buffalo.edu
andrewgodshalk@gmail.com
Created: Some time in 2014 or whatever.
Last updated: 2015-04-22


==Step 1 PATtuples== - created by VHbb Analysis group.
TODO - add lists containing file locations, directions to access
TODO - add description of code (or code itself) used to create PATTuples
TODO - add reference to information, TWIKIs, anything else that might be successful.

==Step 2 Ntuples==
Created with Ntupler executable.
Source: VHbbAnalysis/VHbbAnalyzer/bin/Ntupler.cc

Best done in three-day scratch due to output file size.
GIT repo out of date. For latest version, copy from:
    ~godshalk/nobackup/2015-01_ZHbb_cp/CMSSW_5_3_23_patch1/src/
	- Only need "Ntupler", "VHbbAnalysis", and "ZSV" folders.

Directions for use of Ntupler in condor on LPC:
  Compile - in src, type "scram b"
  Move to "Ntupler" directory
  Copy compiled executable (in release's bin directory. From Ntupler directory: ../../bin/<arch>/Ntupler)
  Archive files for transfer:
    "tar czv --file=zcntuplerFiles.tgz --files-from=fileList_for_condorTar.txt"
  Choose list of files from "formattedFileLists" to run on.
    (Ex: To run over the WW set, chose formatedFileLists/WW.txt)
  Create new directory to run in by copying "condor_runDir_template", cd into.
    (Ex: cp -r condor_runDir_template condor_runDir_WW)
  Configure ntuple_for_condor.py, condor_executable.sh, condor_config.script
    ntuple_for_condor.py - Remember change to run on MC, etc.
    condor_executable.sh - Change "fileName" variable to your formattedFileList choice.
    condor_config.script - Queue = (number of files in file list)/10 (rounded up)
  Run "voms-proxy-init -voms cms" to set up grid credentials.
  Run with "condor_submit condor_config.script".
  Use condor_q <username> to check status.
  When complete, rerun any missing jobs.
  Combine output files using hadd.

TODO - modify Ntupler to extract extra btag variables
TODO - create README for Ntupler
TODO - write description of files needed, how they are used, etc.
TODO - see if you can compile Ntupler.cc from Ntupler folder


==Step 3 Selection Plots==
TODO - write readme
TODO - transfer code to the server
TODO - commit latest plotting macros to git

==Step 4 Combined Plots & Measurement==
TODO - write readme
TODO - transfer code to the server


Ntupler.cc Dependencies contained in repo
VHbbAnalysis/VHbbDataFormats/
     BuildFile.xml
     bin/
          BuildFile.xml
          Ntupler.cc
     interface/
          HbbCandidateFinderAlgo.h
          TopMassReco.h
          VHbbCandidate.h
          VHbbEvent.h
          JECFWLite.h
          TriggerReader.h
          VHbbCandidateTools.h
          VHbbEventAuxInfo.h
          BTagWeight.h
          BTagReshaping.h
          TriggerWeight.h
          MultiThresholdEfficiency.h
     src/HbbCandidateFinderAlgo.cc
ZSV/BAnalysis/interface/
     SimBHadron.h
     SimSecondaryVertex.h

Dependencies: VHbbAnalysis/VHbbDataFormats/interface/BTagReshaping.h
VHbbAnalysis/VHbbDataFormats/interface/
     btag_payload_b.h
     btag_payload_light.h



-2015-04-20-
-Changed Ntupler to run on newer, correct JSON file (files_to_tar, ntupler_for_condor.py
    in template folder)

-2015-10-
Added a *MESS* of scripts to help run the Ntupler, as well as a few extra data points to the Ntuples themselves.
Full documentation of scripts and the process to come.

```

