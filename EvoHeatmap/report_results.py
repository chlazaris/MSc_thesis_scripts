#!/usr/bin/env python

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

#filename is the output of RoundUp
filename = sys.argv[1]
filename = filename.rstrip()

#show the path of the directory where filename is
filename_parent_dir = os.path.dirname(os.path.realpath(filename))

filehandle = open(filename, "r")
lines = filehandle.readlines()

# create an array to save 
# report content 
report_content = []
report_content.append("*****Query settings*****")

for line in lines:
    line = line.rstrip()
    
    #print query subject and settings
    
    if line.startswith("PA"):
        line_content = line.split("\t")
        query = line_content[1]
        subject = line_content[2]
        divergence = line_content[3]
        e_value = line_content[4]

        report_content.append("Query: %s" %query)
        report_content.append("Subject: %s" %subject) 
        report_content.append("Divergence: %s" %divergence)
        report_content.append("e-value: %s\n" %e_value)	

    # get query protein name and orthologs
    
    if line.startswith("OR"):
        line_content = line.split("\t")
        query_id = line_content[1]
        ortho_id = line_content[2]
        max_lkhd = line_content[3]

        # use query_id to retrieve the name 
        # of query protein
        query_annotation = retrieve_annotation(query_id)
        query_data = print_data(query_annotation)
        query_protein = query_data[1]

        
        # use ortho_id to retrieve the name
        # of orthologous protein
        ortholog_annotation = retrieve_annotation(ortho_id)
        ortholog_data = print_data(ortholog_annotation)
        ortholog_protein = ortholog_data[1]
        
        
        #print ("Orthology results")
        report_content.append("Query protein: %s" %query_protein)
        report_content.append("Orthologous protein: %s" %ortholog_protein)
        report_content.append("Maximum likelihood evolutionary distance: %s\n" %max_lkhd)
	
# print report on screen
print("\n".join(report_content))
        
# create final report file
report_name_parts = filename.split("_")
print report_name_parts
report_name = report_name_parts[0]  + "_" + report_name_parts[1] + "_" + report_name_parts[2]  + "_" + report_name_parts[3] + "_" + report_name_parts[4]  + "_" + report_name_parts[5] + "_" + report_name_parts[6]  + "_" + report_name_parts[7] + "_" + "orthologs" + "_" + "complete.txt"
 
#write the contents of report      
filehandle = open(report_name,"w")
filehandle.write("\n".join(report_content))

#Inform the user
print("The report results have been saved in %s" %(report_name))


        
        
