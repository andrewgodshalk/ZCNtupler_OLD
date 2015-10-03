
from os import system
from os import listdir

# Function to open a file, find a line starting with something, and replace it with something else.
def FindLineInFileStartingWithAndReplaceWith(fileName,oldLineBeginning,newLine) :
    # Open file, read in lines
    fileToMod = open(fileName, 'r')
    fileLines = fileToMod.readlines()
    fileToMod.close()
    
    # Find the line in question
    iToMod = -1
    for line in fileLines:
        if line.startswith(oldLinePart) :
           iToMod = fileLines.index(line)
    
    # Change line to appropriate True/False
    fileLines[iToMod] = newLine
    
    # Read all lines out to a new file.
    fileToMod = open(fileName, 'w')
    for line in fileLines :
        fileToMod.write(line)
    fileToMod.close()







# Strings to find
THE_STRINGS = [
  "Aborted                 ./Ntupler ntuple_for_condor.py",
  "Segmentation fault      ./Ntupler ntuple_for_condor.py",
]

# Get list of datasets
countFile = open('../formattedFileLists/file_counts_per_list.txt', 'r')
datasets = []
for line in countFile :
    datasets.append(line.split()[0])

redoDir = "condor_runDir_redo"

# In each direcotry...
for ds in datasets:
  # Check each stderr for THE STRING
    dirName = "../condor_runDir_"+ds
  # Create a list of stderr files
    fileList = listdir(dirName)
    fileList = [ item for item in fileList if item.endswith(".stderr")]
  # If the file contains THE STRING...
    for errFile in fileList :
        for THE_STRING in THE_STRINGS :
            if THE_STRING in open(dirName+"/"+errFile).read() :
              # Save the corresponding config file to a new folder.
                jobNumber = errFile.split('_')[2].split('.')[0]
                newFileName = "ntuple_job_"+ds+"_"+jobNumber+".py"
                system("cp "+dirName+"/ntuple_job_"+jobNumber+".py "+newFileName)
                #system("cp "+dirName+"/"+errFile+" "+ds+"_"+errFile)
              # Change the output of the config file to something similar
                FindLineInFileStartingWithAndReplaceWith(newFileName, "fname = 'Test' + channel +", "fname = 'Test' + channel + '_"+ds+"_"+jobNumber+".root'\n")
                break

# Get list of config files to redo
fileList = [item for item in listdir('.') if item.startswith("ntuple_job_")]

# Write list to file
configFileListFile = open("configFileList.txt", 'w')
for fn in fileList :
    configFileListFile.write(fn+"\n")
configFileListFile.close()

# Compress all files into a tar.
tarCommand = "tar czv --file=configFiles.tgz --files-from=configFileList.txt"
print tarCommand
system(tarCommand)
system("rm ntuple_job_*")

# Change final line of condor_config.script to show the correct number of jobs
fileName = "condor_config.script"
oldLinePart = "Queue"
newLine = "Queue "+str(len(fileList))+"\n"
FindLineInFileStartingWithAndReplaceWith(fileName,oldLinePart,newLine)

