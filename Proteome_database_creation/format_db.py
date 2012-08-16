#!/usr/bin/env python

"""
This is a script which is used to make
databases searchable by BLAST from FASTA
files. 
"""

import os

# Ask for the directory containing 
# the files
dir_name = raw_input("Give the directory name with FASTA files:\n")
dir_name = dir_name.rstrip()

# Get the path
path = os.path.realpath(dir_name)
os.chdir(path)
# list the contents of the directory 
dirList=os.listdir(os.getcwd())
#print dirList
for fname in dirList:
    if fname.endswith("refseq_prot.fasta"):
        #print fname
        dbname = fname[0:9]+"_blastp_db"
        #print dbname
        command = "makeblastdb -dbtype prot -in %s -out %s" %(fname, dbname)
        os.system(command)
    else:
        print(fname + " was not of the proper file type")
