#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  8 15:03:36 2018

@author: damian
"""

import os as os
from Bio import SeqIO

#Specify working directory if necessary
os.chdir("/media/damian/Memorex USB/comp_mutualism")

file_list = os.listdir("./processed_data/longest_protein")
protein_list = []
for item in file_list:
    if "txtp" in item:
        protein_list.append(item)
    else:
        print "Not appending " + item
print protein_list
print len(protein_list)

#os.makedirs("./processed_data/refseq_longest/")
for item in protein_list:
    holder = item.split("_")
    FAA_holder = holder[0] + "_refseq_" + holder[1] + "_" + holder[2].replace(".txtp", ".FAA")
    max_protein = []
    max_file = open("./processed_data/longest_protein/" + item)
    for protein in max_file:
        if "\n" in protein:
            protein = protein.strip("\n")
            max_protein.append(protein)
        elif "\n" not in protein:
            max_protein.append(protein)
        else:
            print "Something went wrong making list of longest isoforms!"
    print item + " has this many lines: " + str(len(max_protein))
    records = SeqIO.parse("./input_data/refseq_protein/" + FAA_holder, "fasta")
    print records
    all_proteins = []
    for entry in records:
        all_proteins.append(entry)
    max_FAA_list = []
    FAA_check = []
    for record in all_proteins:
        holder = record.id
        holder2 = holder.split("|")
        accession = holder2[3]
        if accession in max_protein:
            max_FAA_list.append(record)
            FAA_check.append(accession)
        elif accession not in max_protein:
            pass
        else:
            print "Something went wrong selecting longest isoforms from FAA file."
    print "Longest isoforms in FAA file: " + str(len(max_FAA_list)) + " of " + str(len(max_protein)) + " (" + str((float(len(max_FAA_list)))/float(len(max_protein))) + ")" + " found."
    print ""
    missing_list = []
    for entry in max_protein:
        if entry not in FAA_check:
            missing_list.append(entry)
            print entry
        else:
            pass
    with open("./misc/missing_proteins/"+ item + "m", "w") as f:
        f.write("\n".join(missing_list))
    SeqIO.write(max_FAA_list, "./processed_data/refseq_longest/" + FAA_holder + "p", "fasta")