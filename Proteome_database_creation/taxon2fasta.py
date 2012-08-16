#!/usr/bin/env python

"""
This is a Python script where the Entrez taxon ID is used
as input and the output is a FASTA file containing
the whole proteome for the organism or organisms that
correspond to this taxon ID.
"""

# This Python code comes from
# http://www.biostars.org/post/show/17715/refseq-proteins-for-a-given-taxid/
# and Biopython manual just with slight modifications

from Bio import Entrez

taxon_id = raw_input("Please provide the NCBI taxon ID:\n")
taxon_id = taxon_id.rstrip()

entrezDbName = 'protein'
ncbiTaxId = taxon_id # Taxon ID provided
Entrez.email = 'history.C.H.Lazaris@sms.ed.ac.uk'

# Find entries matching the query
entrezQuery = 'txid%s[orgn] AND "refseq"[Filter]' %(ncbiTaxId)
searchResultHandle = Entrez.esearch(db=entrezDbName, 
                                    term=entrezQuery, 
                                    retmax=100000, 
                                    usehistory="y") # The value for returned 
                                                    # results set to max, 
                                                    # usehistory set to yes
search_result = Entrez.read(searchResultHandle)
searchResultHandle.close()

#print search_result

# Get the counts
gi_list = search_result["IdList"]
count = int(search_result["Count"])
assert count == len(gi_list)

print(count, "records will be downloaded")


# Add the necessary cookies to use the History and 
# retrieve data in steps

webenv = search_result["WebEnv"]
query_key = search_result["QueryKey"]

batch_size = 5
out_handle = open("txid%s_refseq_prot.fasta" %(ncbiTaxId), "w")
for start in range(0, count, batch_size):
    end = min(count, start+batch_size)
    print "Going to download record %i to %i" % (start+1, end)
    fetch_handle = Entrez.efetch(db="protein", rettype="fasta", retmode="text",
                                    retstart=start, retmax=batch_size,
                                    webenv=webenv, query_key=query_key)
    data = fetch_handle.read()
    fetch_handle.close()
    out_handle.write(data)
out_handle.close()

# print completion message
print("***Download is now complete-FASTA has been created***")