#!/usr/bin/env python

# Author: Ana Rodriguez Ronchel
# Date: 18-11-2020
# Description: This script generates 2 bed files with the coordinates of the
#              genome regions that could correspond with lost or retained loops
#              mediated by CTCF. Conditions to be fulfilled: 
#                 - Regions of less than 1 Mb
#                 - Regions flanked by two lost CBSs (for lost loops) or two
#                   retained CBSs (for retained loops).
#                 - The motif at the left CBS must be at the + strand and
#                   the motif at the right CBS at the - strand
#                 - A loop can not be formed by two motifs present in the same
#                   CBSs (too small distance).
# Input:
#       - Arg1: 4 column bed file with lost CBSs (chr, star, end, motif strand)
#       - Arg2: 4 column bed file with retained CBSs (")
# Output:
#       - 3 column bed file for potential lost loops (chr, start, end)
#       - 3 column bed file for potential retained loops (")
#
# Execution example:
# ./Potential_loop_selection_python.py ../Results/Annotated_files/Annotated_Lost_CBSs.bed ../Results/Annotated_files/Annotated_Retained_CBSs.bed

################################# MAIN #########################################

# Importing arguments and assigning variables:
import sys
lost_CBSs = sys.argv[1]
retained_CBSs = sys.argv[2]

# Make a dictionary of lost CBSs
chr_dict_lost = {}
with open(lost_CBSs, 'r') as lost_file:
    for line in lost_file:
        chr = line.strip().split()[0]
        if chr not in chr_dict_lost.keys():
            chr_dict_lost[chr]=[]
        chr_dict_lost[chr].append(line.strip().split()[1:4]) #Dict with a list
                                                             #of (start, end,
                                                             #strand) for each
                                                             #chr.

# Make a dictionary of retained CBSs
chr_dict_retained = {}
with open(retained_CBSs, 'r') as retained_file:
    for line in retained_file:
        chr = line.strip().split()[0]
        if chr not in chr_dict_retained.keys():
            chr_dict_retained[chr]=[]
        chr_dict_retained[chr].append(line.strip().split()[1:4])


# Open a new file to write potential lost loops:
lost_file_write = open("../Results/Potential_loop_file/potential_lost_loops.bed", 'w')

# LOOPS THAT START AND END WITH A LOST CBS:
for key in chr_dict_lost.keys():
    for cbs_1 in chr_dict_lost[key]: #First loop that select only + motifs
        if cbs_1[2] != "+":
            continue
        for cbs_2 in chr_dict_lost[key]: #Second loop that select - motifs
            if int(cbs_2[0]) - int(cbs_1[0]) <= 0: continue #To avoid selecting two motifs of the same CBS.
            if int(cbs_2[1]) - int(cbs_1[0]) >= 1000000: break
            if cbs_2[2] == "-":
                possible_loop=[key,cbs_1[0],cbs_2[1]]
                for column in possible_loop:
                    lost_file_write.write(column)
                    lost_file_write.write("\t")
                lost_file_write.write("\n")

lost_file_write.close()


# Open a new file to write potential retained loops:
retained_file_write = open("../Results/Potential_loop_file/potential_retained_loops.bed", 'w')

# LOOPS THAT START AND END WITH A RETAINED CBS:
for key in chr_dict_retained.keys():
    for cbs_1 in chr_dict_retained[key]: #First loop that select only + motifs
        if cbs_1[2] != "+":
            continue
        for cbs_2 in chr_dict_retained[key]: #Second loop that select - motifs
            if int(cbs_2[0]) - int(cbs_1[0]) <= 0: continue #To avoid selecting two motifs of the same CBS.
            if int(cbs_2[1]) - int(cbs_1[0]) >= 1000000: break
            if cbs_2[2] == "-":
                possible_loop=[key,cbs_1[0],cbs_2[1]]
                for column in possible_loop:
                    retained_file_write.write(column)
                    retained_file_write.write("\t")
                retained_file_write.write("\n")

retained_file_write.close()
