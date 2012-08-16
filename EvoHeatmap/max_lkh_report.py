#!/usr/bin/env python

"""
Script to turn RoundUp report .txt file
into appropriate format (.csv) for
creation of the corresponding heatmap 
in R. It requires the file with the ids
for the genes and the directory with the
RoundUp results
"""

import sys
import os

from Bio import Entrez

########### Functions ############
# The code for the following two functions
# was derived from 
# http://biopython.org/wiki/Annotate_Entrez_Gene_IDs
# It has been used without modifications 
# *Always* tell NCBI who you are
Entrez.email = "C.H.Lazaris@sms.ed.ac.uk"
 
def retrieve_annotation(protein_id):
 
    """Annotates Entrez Gene IDs using Bio.Entrez, in particular epost (to
    submit the data to NCBI) and esummary to retrieve the information. 
    Returns a list of dictionaries with the annotations."""
 
    request = Entrez.epost("protein",id=protein_id)
    try:
        result = Entrez.read(request)
    except RuntimeError as e:
        #FIXME: How generate NAs instead of causing an error with invalid IDs?
        print "An error occurred while retrieving the annotations."
        print "The error returned was %s" % e
        sys.exit(-1)
 
    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    data = Entrez.esummary(db="protein", webenv=webEnv, query_key =
            queryKey)
    annotation = Entrez.read(data)
 
    return annotation
     
def print_data(annotation):
    for gene_data in annotation:
        gene_id = gene_data["Id"]
        gene_name = gene_data["Title"]
          
        return (gene_id, gene_name)

######### End of functions ###########

#open id file and save contents to array

id_file_name = sys.argv[1]
id_file_name = id_file_name.rstrip()

id_file = open(id_file_name,"r")
id_file_content = id_file.readlines()

# array to save ids
ids = []

#array to save protein names that 
#correspond to ids
protein_names = []

# save ids into the array
for line in id_file_content:
    line = line.rstrip()
    ids.append(line)

# retrieve protein names
# and save in array 
for protein_id in ids:
    annotation = retrieve_annotation(protein_id)
    prot_annotation = print_data(annotation)
    prot_name_parts = prot_annotation[1].split(" ")
    prot_name = prot_name_parts[0]
    protein_names.append(prot_name)

# directory containing the files
ortho_dir_name = sys.argv[2]
ortho_dir_name = ortho_dir_name.rstrip()

#Ask for score to replace missing values
na_value_score = raw_input("Score to replace missing values (X.XXXX):\n")
na_value_score = na_value_score.rstrip()

#list contents of the directory
ortho_dir_files = os.listdir(ortho_dir_name)

# change working dir to ortho_dir_name
os.chdir(ortho_dir_name)

# create array to save the report contents
max_lkh_report = []
# write ids in the first line
max_lkh_report.append("Species," + ",".join(protein_names))

#go through all the files
for ortho_dir_file in ortho_dir_files:
    if ortho_dir_file.endswith("orthologs.txt"):
        f = open(ortho_dir_file, "r")
        ortho_dir_file_content = f.readlines()
 
        # create dictionary 
        d = {}

        for line in ortho_dir_file_content:
            if line.startswith("PA"):
                line_parts = line.split("\t")
                query_prot = line_parts[2]
                query_species_parts = query_prot.split("_")
                query_species = query_species_parts[0] + "_" + query_species_parts[1]
            elif line.startswith("OR"):
                line_parts = line.split("\t")
                d_key = line_parts[1]
                d_value = line_parts[3].rstrip("\n")
                d[d_key] = d_value
        
        #check if d has key equal to the an element in ids

        final_line = []

        final_line.append(query_species)

        for i in range(0,len(ids)):
            key=ids[i]
            if d.has_key(key):
                final_line.append(d[key])
            else:
                final_line.append(na_value_score)
    
        #append the line with all the required 
        #information from the file to the final
        #report
        max_lkh_report.append(",".join(final_line))
    
#set the name for final report file
report_name_part = id_file_name.split("_")
report_name_part1 = report_name_part[0] + "_" + report_name_part[1] + "_" + report_name_part[2]
report_name_part2 = "_max_lkh_report.csv"

final_report_name = report_name_part1 + report_name_part2
     
#print the whole report content
#in the csv file
f = open(final_report_name, "w")
f.write("\n".join(max_lkh_report))





