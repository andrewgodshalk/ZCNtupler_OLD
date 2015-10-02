
from os import listdir
from math import ceil

# Get list of datasets and the number of jobs from file_counts_per_list.txt
countFile = open('formattedFileLists/file_counts_per_list.txt', 'r')
datasets = []
for line in countFile :
    datasets.append(line.split())

# For each dataset...
for ds in datasets :
  # Get the number of jobs and the directory to look at.
    numJobs = int(ceil(float(ds[1])/10.0))
    dirName = "condor_runDir_"+ds[0]
  # Create a list of job
    fileList = listdir(dirName)
    fileList = [ item for item in fileList if not item.startswith("TestMuon_")]
    listOfIncompleteJobs = []
    strListOfIncompletes = ""
  # For each job number, check if there is a finished product, and add the mess to a list if no product is present.
    for i in range(numJobs) :
        if fileList.index("TestMuon_"+str(i)+".root") == -1 :
            listOfIncompleteJobs.append(i)
            strListOfIncompletes += str(i)+", "
  # Print the incomplete jobs
    print dr[0]+": "+str(listOfIncompleteJobs.length())+" of "+dr[1]+" jobs processed. Incompletes: "+strListOfIncompletes


