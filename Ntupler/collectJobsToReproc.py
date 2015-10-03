
from os import system

# String to find
THE_STRING = "Aborted                 ./Ntupler ntuple_for_condor.py"

# Get list of datasets
countFile = open('formattedFileLists/file_counts_per_list.txt', 'r')
datasets = []
for line in countFile :
    datasets.append(line.split()[0])

system("mkdir condor_runDir_do-over")

# In each direcotry...
for ds in datasets:
  # Check each stderr for THE STRING
    dirName = "condor_runDir_"+ds
  # Create a list of stderr files
    fileList = listdir(dirName)
    fileList = [ item for item in fileList if item.endswith(".stderr")]
  # If the file contains THE STRING...
    for errFile in fileList :
        if THE_STRING in open(errFile).read() :
          # Save the corresponding config file to a new folder.
            jobNumber = errFile.split('_')[2].split('.')[0]
            system("cp "+dirName+"/ntuple_job_"+jobNumber+".py condor_runDir_do-over/")
