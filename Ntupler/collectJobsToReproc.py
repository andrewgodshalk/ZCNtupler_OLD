
from os import system
from os import listdir

# Strings to find
THE_STRINGS = [
  "Aborted                 ./Ntupler ntuple_for_condor.py",
  "Segmentation fault      ./Ntupler ntuple_for_condor.py",
]


def findAndReplaceLineInFile(fileName, oldLine, newLine):
  # Open file, read in lines
    fileToMod = open(fileName, 'r')
    fileLines = fileToMod.readlines()
    fileToMod.close()
  # Find the line in question
    iToMod = fileLines.index(oldLine)
    #print iToMod
  # Change line to appropriate True/False
    fileLines[iToMod] = newLine
  # Read all lines out to a new file.
    fileToMod = open(fileName, 'w')
    for line in fileLines :
        fileToMod.write(line)
    fileToMod.close()



# Get list of datasets
countFile = open('formattedFileLists/file_counts_per_list.txt', 'r')
datasets = []
for line in countFile :
    datasets.append(line.split()[0])

redoDir = "condor_runDir_redo"

# In each direcotry...
for ds in datasets:
  # Check each stderr for THE STRING
    dirName = "condor_runDir_"+ds
  # Create a list of stderr files
    fileList = listdir(dirName)
    fileList = [ item for item in fileList if item.endswith(".stderr")]
  # If the file contains THE STRING...
    for errFile in fileList :
	for THE_STRING in THE_STRINGS :
            if THE_STRING in open(dirName+"/"+errFile).read() :
              # Save the corresponding config file to a new folder.
                jobNumber = errFile.split('_')[2].split('.')[0]
                system("cp "+dirName+"/ntuple_job_"+jobNumber+".py condor_runDir_do-over/ntuple_job_"+ds+"_"+jobNumber+".py")
                #system("cp "+dirName+"/"+errFile+" condor_runDir_do-over/"+ds+"_"+errFile)
                break

# Get list of config files to redo
fileList = [item for item in listdir(redoDir) if item.startswith("ntuple_job_")]

# Write list to file
configFileListFile = open(redoDir+"/configFileList.txt", 'w')
for fn in fileList :
    configFileListFile.write(fn)
configFileListFile.close()

# Compress all files into a tar.
system("tar czv --file="+redoDir+"/configFiles.tgz --files-from="+redoDir+"/configFileList.txt")

# Change final line of condor_config.script to show the correct number of jobs
findAndReplaceLineInFile(redoDir+"/condor_config.script","Queue 1\n", "Queue "+str(len(fileList))+"\n")
