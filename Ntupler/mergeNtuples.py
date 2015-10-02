from os import system
import datetime
import time

# Create an ntuple directory with a timestamp
ts = time.time()
ntupleDir = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d') + "_step2ntuples"
system("mkdir "+ntupleDir)

# Get list of datasets
countFile = open('formattedFileLists/file_counts_per_list.txt', 'r')
datasets = []
for line in countFile :
    datasets.append(line.split()[0])

# For each dataset in file counts list...
for ds in datasets:
  # Get the appropriate condor directory
    condDir = "condor_runDir_"+ds
  # Combine files using root's "hadd" function
    system("hadd "+ntuplerDir+"/"+ds+".root "+condDir+"/TestMuon*.root")

# Modify ttbar titles accordingly
system("mv "+ntupleDir+"/ttbar_lep.root "+ntupleDir+"/ttlep.root ")
system("mv "+ntupleDir+"/ttbar_semi.root "+ntupleDir+"/ttsemi.root ")
system("mv "+ntupleDir+"/ttbar_had.root "+ntupleDir+"/tthad.root ")
# Combine muon and electron ntuples into one place.
system("hadd "+ntuplerDir+"/muon.root "+ntupleDir+"/muon*.root")
system("hadd "+ntuplerDir+"/elec.root "+ntupleDir+"/elec*.root")

