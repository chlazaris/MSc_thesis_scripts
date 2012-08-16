#!/usr/bin/env python

"""
This is a script that is used for 
automatic BLAST of FASTA sequences in a
folder against formatted protein databases.
It accepts as input the directory with the 
databases and the protein files to be 
BLAST-searched against the databeses.
"""
import os

# Ask for the query file

# Ask for the directory containing 
# the files
dir_name = raw_input("Give the directory with databases:\n")
dir_name = dir_name.rstrip()

# List directory contents
dirList=os.listdir(dir_name)

#print dirList
print dirList

# Change path to the directory with files
os.chdir(dir_name)

# Create arrays to hold names
# for protein files and databases

protein_files = []
db_files = []

for file in dirList:
    if file.endswith(".fasta"):
        protein_files.append(file)
    elif file.endswith("blastp_db.psq"):
        db_files.append(file[0:19])

# print the lists to be sure
print protein_files
print db_files

# create a new directory to 
# save the BLAST results

new_dir = raw_input("BLAST result directory:\n")
new_dir = new_dir.rstrip()
os.makedirs(new_dir)


# perform BLAST for all protein files
# against all databases and export
# results
for protein_file in protein_files:
    for db_file in db_files:
        output = protein_file[0:8] + "_vs_" + db_file[0:10] + "blastp.txt"
        os.system("blastp -query %s -db %s -outfmt 7 | uniq > ./%s/%s" %(protein_file,db_file,new_dir,output))
        
# Print completion message
print("***The BLAST process is now complete***")


