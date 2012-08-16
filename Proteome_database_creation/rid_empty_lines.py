#!/usr/bin/env python

import os

"""
This is a script that is run outside a directory
containing not properly formatted (for rsd_search)
refseq_prot.fasta files. It formats them properly
(removes empty lines) and saves the new files into
a "converted_files" directory inside the parent 
directory.
"""

#Source used for "Get rid of empty lines" part: 
#http://ubuntuforums.org/showthread.php?t=302914

# Ask for folder with files

filefolder = raw_input("What is the folder with the files for conversion?\n")
filefolder = filefolder.rstrip()

# list the files in the folder
dirList = os.listdir(filefolder)

# change path to the folder
# with the files
os.chdir(filefolder)

# Create a new directory
os.mkdir("converted_files")

for file in dirList:
    
    if file.endswith("refseq_prot.fasta"):

        # open file
        filehandle = open(file)
        lines = filehandle.readlines()
        filehandle.close()

        # create array for not 
        # blank lines
        contents = []

        # Get rid of empty lines
        for line in lines:
            # Strip whitespace. If line 
            # was empty, nothing is left
            if not line.strip():
                continue
            # otherwise, save the line
            else:
                contents.append(line)
            
        # Print file with empty lines removed
        new_file = file.split(".")[-2] + "_new" + ".aa"
        filehandle = open(new_file, "w")
        filehandle.write("".join(contents))

        # write the converted files
        # into the converted_files dir
        os.system("mv %s ./converted_files" %new_file)
    
 
