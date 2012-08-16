#!/usr/bin/env python

"""
This script retrieves the IDs of the
top hits of files containing BLAST
results and then uses these IDs to 
talk to Entrez and retrieve the 
corresponding protein sequences.

The top hit sequences are saved into
seperate files ending with:
"BLAST_top_hit_IDs.fasta" 

This script is to be run outside a folder
containing .txt files with BLAST results.
The user provides the name of the directory 
with the BLAST results and also the first part 
of the name describing the files: e.g Sc_Gal1p
to get the information from the files
where Sc Gal1 protein was BLAST-searched
against all the others. 
"""

import os
import csv
from Bio import Entrez, SeqIO

########## Functions ##############

# This function used as is from: 
#http://www.peterbe.com/plog/uniqifiers-benchmark/

def f2(seq): 
   # order preserving
   checked = []
   for e in seq:
       if e not in checked:
           checked.append(e)
   return checked
   
##################################

# show current directory
#print os.getcwd()

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

#create an array to keep all 
#you want from the file
ID_report = []

for file in select_files:
    
    # extract the species name
    split_name = file.split("_")
    species_name = split_name[3] + "_" + split_name[4]
    
    # and the name of query protein
    query_protein = split_name[0] + "_" + split_name[1]
    
    # open the file
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
        continue
    elif len(lines) == 1:
        continue
    else:
        #get the lines corresponding to 2
        #highest bit scores
        line1 = lines[0].rstrip()
        line2 = lines[1].rstrip()
    
        #split them
        split_line1 = line1.split("|")
        split_line2 = line2.split("|")
        
        ID_blast_hit_1 = split_line1[1]
        ID_blast_hit_2 = split_line2[1]

        #append to the ID_report
        ID_report.append(ID_blast_hit_1) 
        ID_report.append(ID_blast_hit_2)
          
# keep just the unique entries of ID_report
ID_report = f2(ID_report)

fasta_file = []

# Retrieve FASTA files
for id in ID_report:
    Entrez.email = "C.H.Lazaris@sms.ed.ac.uk" #Tell the NCBI who you are
    handle = Entrez.efetch(db="protein", id=id, rettype="fasta")
    record = handle.read()
    fasta_file.append(record)
    
# print all what you got
# into a .fasta file
fasta_filename = query_protein + "_" + "BLAST_top_hit_IDs" + ".fasta"
f = open("%s" %fasta_filename, "w")
f.write("\n".join(fasta_file))
f.close()


