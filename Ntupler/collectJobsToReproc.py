
# String to find:  "Aborted                 ./Ntupler ntuple_for_condor.py"


# Get list of datasets
countFile = open('formattedFileLists/file_counts_per_list.txt', 'r')
datasets = []
for line in countFile :
    datasets.append(line.split()[0])

# In each direcotry...
for ds in datasets:
  # Check each stderr for THE STRING
    
  # If they contain the string...
    # Save the corresponding config file to a new folder.


