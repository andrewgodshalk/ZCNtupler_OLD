from os import system
import datetime
import time

# Create an ntuple directory with a timestamp
ts = time.time()
ntuplerDir = datetime.datetime.fromtimestamp(ts).strftime('%Y-%m-%d') + "_step2ntuples"
system("mkdir "+ntuplerDir)

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
#system("mv "+ntuplerDir+"/ttbar_lep.root " +ntuplerDir+"/ttlep.root" )
#system("mv "+ntuplerDir+"/ttbar_semi.root "+ntuplerDir+"/ttsemi.root")
#system("mv "+ntuplerDir+"/ttbar_had.root " +ntuplerDir+"/tthad.root" )
# Combine muon and electron ntuples into one place.
system("hadd "+ntuplerDir+"/muon.root "+ntuplerDir+"/muon*.root")
system("hadd "+ntuplerDir+"/elec.root "+ntuplerDir+"/elec*.root")

