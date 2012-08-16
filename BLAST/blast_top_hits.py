#!/usr/bin/env python

"""
This script is to be run outside a folder
containing .txt files with BLAST results.
The user provides the first part of the 
name describing the files: e.g Sc_Gal1p
to get the information from the files
where Sc Gal1 protein was BLASTED against
all the others. The user also provides
the name of the directory with the BLAST
files.
"""

import os
import csv

# show current directory
print os.getcwd()

# list the contents of the folder
dir_name = raw_input("What is the folder with the files?\n")
filenames = os.listdir(dir_name)
#print filenames

# move to this directory to
# be able to work with files
os.chdir(dir_name)
#print os.getcwd()

# read the right files
name_part = raw_input("What is the set of files?\n")
name_part = name_part.rstrip()

# check if the filename starts with this
# and if yes keep files in array
select_files = []
for filename in filenames:
    if filename.startswith(name_part) and filename.endswith("blastp.txt"):
        select_files.append(filename)

#an array to store everything needed
#for all selected files
csv_file = []

firstline = []
firstline.append("Species")
firstline.append("Bit_score_1")
firstline.append("Bit_score_2")
firstline = ",".join(firstline)

csv_file.append(firstline)

# open the selected files and 
# extract what you want

for file in select_files:
    
    # extract the species name
    split_name = file.split("_")
    species_name = split_name[3] + "_" + split_name[4]
    
    # and the name of query protein
    query_protein = split_name[0] + "_" + split_name[1]
    
    f = open(file, "r")
    file_content = f.readlines()
    f.close()

    #array to store lines 
    lines = []

    for line in file_content:
        if not line.startswith("#"):
            lines.append(line)
    
    #test the number of BLAST results and
    #report accordingly:
    
    if len(lines) == 0:
        print("No BLAST hits in file " + file)
    elif len(lines) == 1:
        # get the line
        line1 = lines[0].rstrip()
        
        #split it
        split_line1 = line1.split("\t")
        
        #create an array to keep all 
        #you want from the file
        report = []
        report.append(species_name)
        report.append(split_line1[11]) 
        # Give zero as second bit score
        report.append("0")
        
        # put all information in a single line
        final_line = ",".join(report)
    
        #append this line to total file
        csv_file.append(final_line)
    else:
        #get the lines corresponding to 2
        #highest bit scores
        line1 = lines[0].rstrip()
        line2 = lines[1].rstrip()
    
        #split them
        split_line1 = line1.split("\t")
        split_line2 = line2.split("\t")

    
        #create an array to keep all 
        #you want from the file
        report = []
        report.append(species_name)
        report.append(split_line1[11]) 
        report.append(split_line2[11])
    
        # put all information in a single line
        final_line = ",".join(report)
    
        #append this line to total file
        csv_file.append(final_line)
          
# print all what you got
# into a csv file
csv_filename = query_protein + "_" + "BLAST_top_hits" + ".csv"
csv_handle = open(csv_filename, 'wb')
i=0
for i in range(0,len(csv_file)):
    csv_handle.write(csv_file[i]+"\n")

