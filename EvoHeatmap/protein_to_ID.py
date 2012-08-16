#!/usr/bin/env python

"""
Script that retrieves gi IDs for corresponding
genes, from files (refseq) and prints them in
new file to be used by RoundUp
"""

import sys
import os
   

# import regular expressions
import re

# import module for autocompletion
import readline

# import module for adding arguments
import optparse

#######Dictionary of fungi species#######

fungi_species = {}

fungi_species["Ajellomyces capsulatus"]        = "Ajel_caps_txid5037_refseq_prot_new.aa"
fungi_species["Ajellomyces dermatitidis"]      = "Ajel_derm_txid5039_refseq_prot_new.aa"
fungi_species["Arthroderma benhamiae"]         = "Arth_benh_txid663331_refseq_prot_new.aa"
fungi_species["Arthroderma gypseum"]           = "Arth_gyps_txid535722_refseq_prot_new.aa"
fungi_species["Arthroderma otae"]              = "Arth_otae_txid554155_refseq_prot_new.aa"
fungi_species["Ashbya gossypii"]               = "Ashb_goss_txid284811_refseq_prot_new.aa"
fungi_species["Aspergillus clavatus"]          = "Aspe_clav_txid344612_refseq_prot_new.aa"
fungi_species["Aspergillus flavus"]            = "Aspe_flav_txid332952_refseq_prot_new.aa"
fungi_species["Aspergillus fumigatus"]         = "Aspe_fumi_txid746128_refseq_prot_new.aa"
fungi_species["Aspergillus nidulans"]          = "Aspe_nidu_txid227321_refseq_prot_new.aa"
fungi_species["Aspergillus niger"]             = "Aspe_nige_txid425011_refseq_prot_new.aa"
fungi_species["Aspergillus terreus"]           = "Aspe_terr_txid341663_refseq_prot_new.aa"
fungi_species["Botryotinia fuckeliana"]        = "Botr_fuck_txid332648_refseq_prot_new.aa"
fungi_species["Candida albicans"]              = "Cand_albi_txid5476_refseq_prot_new.aa"
fungi_species["Candida dubliniensis"]          = "Cand_dubl_txid573826_refseq_prot_new.aa"
fungi_species["Candida glabrata"]              = "Cand_glab_txid5478_refseq_prot_new.aa"
fungi_species["Candida tropicalis"]            = "Cand_trop_txid294747_refseq_prot_new.aa"
fungi_species["Chaetomium globosum"]           = "Chae_glob_txid306901_refseq_prot_new.aa"          
fungi_species["Chaetomium thermophilum"]       = "Chae_ther_txid759272_refseq_prot_new.aa"          
fungi_species["Clavispora lusitaniae"]         = "Clav_lusi_txid306902_refseq_prot_new.aa"         
fungi_species["Coccidioides immitis"]          = "Cocc_immi_txid246410_refseq_prot_new.aa"         
fungi_species["Coccidioides posadasii"]        = "Cocc_posa_txid199306_refseq_prot_new.aa"         
fungi_species["Coprinopsis cinerea"]           = "Copr_cine_txid240176_refseq_prot_new.aa"         
fungi_species["Cryptococcus gattii"]           = "Cryp_gatt_txid367775_refseq_prot_new.aa"         
fungi_species["Cryptococcus neoformans"]       = "Cryp_neof_txid5207_refseq_prot_new.aa"          
fungi_species["Debaryomyces hansenii"]         = "Deba_hans_txid284592_refseq_prot_new.aa"         
fungi_species["Encephalitozoon cuniculi"]      = "Ence_cuni_txid284813_refseq_prot_new.aa"         
fungi_species["Encephalitozoon intestinalis"]  = "Ence_inte_txid876142_refseq_prot_new.aa"         
fungi_species["Enterocytozoon bieneusi"]       = "Ente_bien_txid481877_refseq_prot_new.aa"         
fungi_species["Gibberella zeae"]               = "Gibb_zeae_txid229533_refseq_prot_new.aa"         
fungi_species["Kluyveromyces lactis"]          = "Kluy_lact_txid284590_refseq_prot_new.aa"         
fungi_species["Komagataella pastoris"]         = "Koma_past_txid644223_refseq_prot_new.aa"         
fungi_species["Laccaria bicolor"]              = "Lacc_bico_txid486041_refseq_prot_new.aa"         
fungi_species["Lachancea thermotolerans"]      = "Lach_ther_txid559295_refseq_prot_new.aa"         
fungi_species["Lodderomyces elongisporus"]     = "Lodd_elon_txid379508_refseq_prot_new.aa"         
fungi_species["Magnaporthe oryzae"]            = "Magn_oryz_txid242507_refseq_prot_new.aa"         
fungi_species["Malassezia globosa"]            = "Mala_glob_txid425265_refseq_prot_new.aa"         
fungi_species["Meyerozyma guilliermondii"]     = "Meye_guil_txid294746_refseq_prot_new.aa"         
fungi_species["Moniliophthora perniciosa"]     = "Moni_pern_txid554373_refseq_prot_new.aa"         
fungi_species["Myceliophthora thermophila"]    = "Myce_ther_txid573729_refseq_prot_new.aa"         
fungi_species["Nectria haematococca"]          = "Nect_haem_txid660122_refseq_prot_new.aa"         
fungi_species["Neosartorya fischeri"]          = "Neos_fisc_txid331117_refseq_prot_new.aa"         
fungi_species["Neurospora crassa" ]            = "Neur_cras_txid367110_refseq_prot_new.aa"         
fungi_species["Nosema ceranae"]                = "Nose_cera_txid578460_refseq_prot_new.aa"         
fungi_species["Paracoccidioides brasiliensis"] = "Para_bras_txid121759_refseq_prot_new.aa"         
fungi_species["Penicillium chrysogenum"]       = "Peni_chry_txid500485_refseq_prot_new.aa"         
fungi_species["Penicillium marneffei"]         = "Peni_marn_txid441960_refseq_prot_new.aa"         
fungi_species["Phaeosphaeria nodorum"]         = "Phae_nodo_txid321614_refseq_prot_new.aa"         
fungi_species["Postia placenta"]               = "Post_plac_txid561896_refseq_prot_new.aa"         
fungi_species["Puccinia graminis"]             = "Pucc_gram_txid418459_refseq_prot_new.aa"         
fungi_species["Pyrenophora teres"]             = "Pyre_tere_txid861557_refseq_prot_new.aa"         
fungi_species["Pyrenophora tritici"]           = "Pyre_trit_txid426418_refseq_prot_new.aa"         
fungi_species["Saccharomyces cerevisiae"]      = "Sacc_cere_txid4932_refseq_prot_new.aa"          
fungi_species["Scheffersomyces stipitis"]      = "Sche_stip_txid322104_refseq_prot_new.aa"         
fungi_species["Schizophyllum commune"]         = "Schi_comm_txid578458_refseq_prot_new.aa"         
fungi_species["Schizosaccharomyces japonicus"] = "Schi_japo_txid402676_refseq_prot_new.aa"       
fungi_species["Schizosaccharomyces pombe"]     = "Schi_pomb_txid4896_refseq_prot_new.aa"          
fungi_species["Sclerotinia sclerotiorum"]      = "Scle_scle_txid665079_refseq_prot_new.aa"         
fungi_species["Talaromyces stipitatus"]        = "Tala_stip_txid441959_refseq_prot_new.aa"         
fungi_species["Trichophyton rubrum"]           = "Tric_rubr_txid559305_refseq_prot_new.aa"         
fungi_species["Trichophyton verrucosum"]       = "Tric_verr_txid663202_refseq_prot_new.aa"         
fungi_species["Uncinocarpus reesii"]           = "Unci_rees_txid336963_refseq_prot_new.aa"         
fungi_species["Ustilago maydis"]               = "Usti_mayd_txid237631_refseq_prot_new.aa"         
fungi_species["Vanderwaltozyma polyspora"]     = "Vand_poly_txid436907_refseq_prot_new.aa"         
fungi_species["Verticillium alboatrum"]        = "Vert_albo_txid526221_refseq_prot_new.aa"         
fungi_species["Yarrowia lipolytica"]           = "Yarr_lipo_txid284591_refseq_prot_new.aa"         
fungi_species["Zygosaccharomyces rouxii"]      = "Zygo_roux_txid559307_refseq_prot_new.aa"
                                                                                                         
#########################################

############ Functions ##################

def openfile(filename): 
    """
    Function used to check if a file exists
    and open it.
    """
    try:
        openfile = open(filename, 'r') 
    except Exception:
        print "This file doesn't exist or cannot be read" 
        sys.exit(1)
        
def opendir(d):
    """
    Function used to check if a directory
    exists and list its contents.
    """
    try:
        dirList = os.listdir(d)
    except Exception:
        print("Directory not found")
        sys.exit(1)

#text autocompletion
#Function from:
#http://stackoverflow.com/questions/187621/
#how-to-make-a-python-command-line-program-
#autocomplete-arbitrary-things-not-int

def completer(text, state):
    options = [x for x in fungi_species.keys() if x.startswith(text)]
    try:
        return options[state]
    except IndexError:
        return None

readline.set_completer(completer)
readline.parse_and_bind("tab: complete")

#########################################   

# setup help menu for the program

desc="""
    
Help: Run protein_to_ID.py if you do not already have an ID file and you want to produce one by using the protein names. When asked for species, type the first
letters of the species you are interested in and press TAB. When asked for proteins insert them seperated by commas (eg,Gal1p,Gal2p,Gal3p...). A file with name of the form XX_protein_type_ID.txt
is the output of the program."""

parser = optparse.OptionParser(description=desc)
(args,opts) = parser.parse_args()

#check that the proteome directory exists
opendir("Fungi_refseq_proteomes")

#change to the directory with the proteomes
os.chdir("Fungi_refseq_proteomes")     

# give the name of the proteome folder
#proteome_folder_name = raw_input("Name of proteome folder:\n")
#proteome_folder_name = proteome_folder_name.rstrip()

#change path to proteome folder
#os.chdir(proteome_folder_name)

#provide species name
species_name = raw_input("Select fungi species:\n")
species_name = species_name.rstrip()

#check if this species name exists
#or return error
if fungi_species.has_key(species_name):
    print("Species was found")
    proteome_file_name = fungi_species[species_name]
else:
    print("Species not found")

#try if file exists
#open the file and read it
#f = openfile(proteome_file_name)
f = open(proteome_file_name, "r")
proteome_file_content = f.readlines()

#provide the protein names
protein_names_line = raw_input("Name the proteins (seperated by commas):\n")
protein_names_line = protein_names_line.rstrip()

#store in array
protein_name_array = []

#split the line with the protein names
protein_names = protein_names_line.split(",")

for protein_name in protein_names:
   protein_name_array.append(protein_name)

# create array to save protein IDs
protein_IDs = []

#get each line and split in words

for protein_name in protein_name_array:
    for line in proteome_file_content:
        if re.match("(.*)%s(.*)"%protein_name, line):
             line_parts = line.split("|")
             protein_IDs.append(line_parts[1])
         
#print protein IDs
#in a txt.file
print"""    
Name the file using the convention:
     
XX_protein_type_ID(s)
     
where XX the first two letters from species name 
(eg. Sc for S. cerevisiae) and protein_type, 
the type of proteins (histones etc).
"""

ID_filename = raw_input("ID filename:\n")
ID_filename = ID_filename.rstrip() + ".txt"

f = open("%s" %ID_filename, "w")
f.write("\n".join(protein_IDs))
f.close()

#move the file where the script is
os.system("mv %s ../%s" %(ID_filename,ID_filename))

#notify the user
print("The protein ID file has been created.")



