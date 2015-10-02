#
# 2015-10-02 - takes an input of file_counts_per_list and creates a
#    condor_runDir_XXX for each entry with the proper file modifications.
#
#

# Get list of inputs from file_counts_per_list.txt
# For each entry...
#   Make a new directory
#   Update ntupler_for_condor.py
#     Change isMC to true or false if muon/elec is contained in XXX
#   Update condor_config.script
#     Change number of entries to XXX/10 (rounded UP)
#   Update condor_executable.sh
#     Change input file to XXX.txt

from os import system
from math import ceil

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


# Get list of inputs from file_counts_per_list.txt
countFile = open('formattedFileLists/file_counts_per_list.txt', 'r')
datasets = []
for line in countFile :
    datasets.append(line.split())

#FOR TESTING PURPOSES:
#datasets = [datasets[2]]

# For each entry...
for ds in datasets :
  # Make a new directory from condor_runDir_template
    system("cp -r condor_runDir_template condor_runDir_"+ds[0]+"/")
  # Update ntupler_for_condor.py - Change isMC to true or false if muon/elec is contained in XXX
    if ds[0].find("muon") != -1 or ds[0].find("elec") != -1:
        isMCValue = "False"
    else :
        isMCValue = "True"
    print "For " + ds[0]+ ", changing ntupler_for_condor.py: isMC = "+isMCValue
    newLine = "    isMC 		= cms.bool("+isMCValue+"),\n"
    findAndReplaceLineInFile(
        "condor_runDir_"+ds[0]+"/ntuple_for_condor.py",
        "    isMC 		= cms.bool(True),\n",
        newLine
    )
  # Update condor_config.script - Change number of entries to XXX/10 (rounded UP)
    numJobs = int(ceil(float(ds[1])/10.0))
    findAndReplaceLineInFile(
        "condor_runDir_"+ds[0]+"/condor_config.script",
        "Queue 1\n",
        "Queue "+str(numJobs)+"\n"
    )
  # Update condor_executable.sh - Change input file to XXX.txt
    findAndReplaceLineInFile(
        "condor_runDir_"+ds[0]+"/condor_executable.sh",
        "fileName=\"formattedFileLists/dy.txt\"\n",
        "fileName=\"formattedFileLists/"+ds[0]+".txt\"\n"
    )
    



