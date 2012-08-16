#!/usr/bin/env python

"""
This is a script to ask for settings, do the
reciprocal BLAST search against subject proteomes
in a directory and write results to files. The 
script is to be placed outside the directory
with the subject proteomes. The query file and
ids file are also in the same directory with
the script. A new directory with the results
is created in the directory with the subject
proteomes.
"""

import os
import sys
import select 
import termios 
import tty


# import module for autocompletion
import readline

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

########## Functions ##########

#Code from:
#http://stackoverflow.com/questions/1829/how-do-i-make-a-menu-in-python-
#that-does-not-require-the-user-to-press-enter-t
# This is a menu function

def getkey():
    old_settings = termios.tcgetattr(sys.stdin) 
    tty.setraw(sys.stdin.fileno())
    select.select([sys.stdin], [], [], 0)
    answer = sys.stdin.read(1)
    termios.tcsetattr(sys.stdin, termios.TCSADRAIN, old_settings) 
    return answer

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
    
########## End of functions ##########

#change to proteome directory
os.chdir("Fungi_refseq_proteomes")

# ask for query proteome
query_prot_species = raw_input("Please provide species name for query proteome:\n")
query_prot_species = query_prot_species.rstrip()

#notify the user about what is going to happen
print("\nThe proteome chosen will be compared to all proteomes in Fungi_refseq_proteomes\n")

# query_proteome
query_prot = fungi_species[query_prot_species]

#change path to where the
#script and ID file is
os.chdir(os.pardir)

#the subject directory
subject_prot_dir = "Fungi_refseq_proteomes"

# ask for IDs file:
ids = raw_input("Name of protein IDs file\n")
ids = ids.rstrip()

# ask for divergence
div = raw_input("Set divergence:\n")
div = div.rstrip()

# ask for e-value
e_val = raw_input("e-value:\n")
e_val = e_val.rstrip()

#directory to save the resulting files
output_dir_name_query_part = query_prot.split("_")
output_dir_name_query_part = output_dir_name_query_part[0] + "_" + output_dir_name_query_part[1]
output_dir_name_ID_part = ids.split("_")
output_dir_name_ID_part = output_dir_name_ID_part[1]
output_dir_name_settings_part = "_" + div + "_" + e_val

output_dir_name = output_dir_name_query_part + "_vs_all_" +  output_dir_name_ID_part + output_dir_name_settings_part + "_orthologs"

#create the directory
os.mkdir(output_dir_name)

# get the list of target proteomes
subject_prot_list = os.listdir(subject_prot_dir)

#change directory to proteomes
#in order to do the analysis
os.chdir(subject_prot_dir)

# for each target proteome
for subject_prot in subject_prot_list:
    
    if subject_prot.endswith("refseq_prot_new.aa"):
    
        # set name for output file
        query_name_part = query_prot.split("_")
        subject_name_part = subject_prot.split("_")

        output_query_part = query_name_part[0] + "_" + query_name_part[1] 
        output_subject_part = subject_name_part[0] + "_" + subject_name_part[1]
        settings_part = "_" + div + "_" + e_val
        output = output_query_part + "_vs_" + output_subject_part + settings_part + "_orthologs.txt"

        # run the search and produce resulting file
        os.system("rsd_search -q %s -s %s -o %s --ids ../%s --no-blast-cache --de %s %s " %(query_prot, subject_prot, output, ids, div, e_val))
        
        #move the  file to output directory
        os.system("mv %s ../%s" %(output,output_dir_name))
        
        #run the script to print report on screen
        #and save in file
        os.system("python ../report_results_folder.py ../%s/%s" %(output_dir_name,output))
        
#change path to where the
#script and ID file is
os.chdir(os.pardir)

#run the script for max_lkh_report
os.system("python max_lkh_report.py %s %s" %(ids,output_dir_name))

#Heatmap production to summarise results

#provide the .csv file
csv_file = raw_input("What is the .csv file for the heatmap?\n")
csv_file = csv_file.rstrip()

#Give the main title for the heatmap
heatmap_name = raw_input("Heatmap main title:\n")
heatmap_name = heatmap_name.rstrip()

#ask how the heatmap should be printed. With 
#clustering per row or not.

print """Menu 
         1) Heatmap - active value clustering 
         2) Heatmap - no active clustering
      """
         
option = getkey()

if "1" in option:
    os.system("Rscript heatmap_function.R %s/%s %s" %(output_dir_name,csv_file,heatmap_name))
    print("You have chosen to produce a heatmap with clustering (dendrograms)")
if "2" in option:
    os.system("Rscript ordered_heatmap_function.R %s/%s %s" %(output_dir_name,csv_file,heatmap_name))
    print("You have chosen to produce a heatmap with clustering (dendrograms)")
    
#Inform the user
print("The heatmap has been now generated")
print("EvoHeatmap run is complete.")



        
        