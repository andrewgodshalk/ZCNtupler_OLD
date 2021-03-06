
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
    fileList = [ item for item in fileList if item.startswith("TestMuon_")]
    listOfIncompleteJobs = []
    strListOfIncompletes = ""
  # For each job number, check if there is a finished product, and add the mess to a list if no product is present.
    for i in range(numJobs) :
        try:
            fileList.index("TestMuon_"+str(i)+".root")
        except ValueError:
            listOfIncompleteJobs.append(i)
            strListOfIncompletes += str(i)+", "
  # Print the incomplete jobs
    print ds[0]+": "+str(len(listOfIncompleteJobs))+" of "+str(numJobs)+" jobs not processed."
    if strListOfIncompletes != "": print "    Incompletes: "+strListOfIncompletes

