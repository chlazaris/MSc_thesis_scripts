#!/usr/bin/env python

"""
This is a script that is used to run all the
scripts required for average evolutionary 
distance (AED) heatmap production using 
RoundUp and R.
"""

import sys
import os
import select 
import termios 
import tty


############# Functions ###############

#Code from:
#http://stackoverflow.com/questions/1829/how-do-i-make-a-menu-in-python-
#that-does-not-require-the-user-to-press-enter-t

def getkey():
    old_settings = termios.tcgetattr(sys.stdin) 
    tty.setraw(sys.stdin.fileno())
    select.select([sys.stdin], [], [], 0)
    answer = sys.stdin.read(1)
    termios.tcsetattr(sys.stdin, termios.TCSADRAIN, old_settings) 
    return answer
    
#######################################

#Run script to produce the protein IDs file 
#required by RSD algorithm.

print """
      ##########  EvoHeatmap v0.1 ##########
      
      EvoHeatmap is a program which consists of Python and R scripts
      and performs putative ortholog detection of proteins used as 
      query. The search can be performed either by using the names
      of the proteins or a file containing the corresponding Entrez IDs. 
      The subject of the search can be either a single proteome or
      the whole set of proteomes covered in the local database (67
      species represented in August 2012). The ortholog detection is
      performed by RSD (Reciprocal Smallest Distance) algorithm and 
      PAML (Phylogenetic Analysis with Maximum Likelihood). The results
      are saved as text files and a comma-seperated file which contains
      the maximum likelihood evolutionary distance values. The final output
      is either a text file when a set of proteins is searched against 
      a single proteome, or heatmap (either clustered or unclustered)
      when a set of proteins is searched against the whole database.
      
      Copyright (c) 2012, Harris A. Lazaris
      All rights reserved.

      *** License information ***
     
      Redistribution and use in source and binary forms, with or without
      modification, are permitted provided that the following conditions are met: 

      1. Redistributions of source code must retain the above copyright notice, this
         list of conditions and the following disclaimer. 
      2. Redistributions in binary form must reproduce the above copyright notice,
         this list of conditions and the following disclaimer in the documentation
         and/or other materials provided with the distribution. 

      THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
      ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
      WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
      DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
      ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
      (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
      LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
      ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
      (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
      SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

      The views and conclusions contained in the software and documentation are those
      of the authors and should not be interpreted as representing official policies, 
      either expressed or implied, of the FreeBSD Project.
      """

print """
      Welcome to EvoHeatmap v.0.1
      
      To start, please select one of the options below: 
      Press '1' if you have the names of the proteins (eg. Gal1p) but not the IDs.
      Press '2' if you already have the protein ID file.
      Press '3' if you need help.
      Press '4' to quit.
      """
print """
      Menu
      You are about to run protein_to_ID.py:
      1) Run protein_to_ID.py
      2) Skip this step (if you have the ID file)
      3) Help (protein_to_ID.py usage)
      4) Quit"""
      
option = getkey()

if "1" in option:
    os.system("python protein_to_ID.py")
elif "2" in option:
    print("\n*****Continued to next step*****\n")
elif "3" in option:
    os.system("python protein_to_ID.py -h")
elif "4" in option:
    sys.exit() 
else:
    print("Please insert the right key")
    sys.exit() 
     
#Run RoundUp itself
#Provide two options. If you compare
#two proteomes (option 1) or a proteome
#with folder of proteomes (option 2).

print """

      Please select the type of search you would like to perform:
      Press '1' if you would like to search for orthologs against
      a certain proteome.
      Press '2' if you would like to search against the local 
      database.
      Press '3' to quit.
          
      Menu 
      1) Query protein(s) (Proteome A) vs. Proteome B. 
      2) Query protein(s) (Proteome A) vs. Local Proteome Database
      3) Quit program
      """
         
option = getkey()

if "1" in option:
    os.system("python roundup_rec_blast.py")
elif "2" in option:
    os.system("python roundup_rec_folder.py")
elif "3" in option:
    sys.exit()
else:
    print("Please insert the right key")
    sys.exit() 
    
    
         


