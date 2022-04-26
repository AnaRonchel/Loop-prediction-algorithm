#!/bin/bash

# Author: Ana Rodriguez Ronchel
#
# Date: 18-11-2020
#
###################################### DESCRIPTION ############################
#
# This script generates a bed file with the coordinates of the genome regions
# that could correspond with DNA loops that were present in the control condition
# but not in the CTCF-deficient condition. For more details about the processing
# of the files, go to the python script "Potential_loop_selection_python.py".
#
# INPUT:
#       - Arg1: 4 column bed file with lost peaks (chr, star, end, motif strand)
#       - Arg2: 4 column bed file with retained peaks (")
#
# OUTPUT:
#       - A 3 column bed file with potential loops coordinates (chr, start, end)
#
###################################### REQUIREMENTS ###########################
#
# The directory must be set to the file where this script is located.
#
# The parent directoy must contain the next directories:
#     - Results/Annotated_files: Annotated lost and retained peaks bed files.
#     - Scripts: this script and the script "Potential_loop_selection_pyhton.py"
#
# Installed programs: Python
#
##################################### RUN EXAMPLE #############################
#
# ./Potential_loop_selection.sh ../Results/Annotated_files/Annotated_Lost_CBSs.bed ../Results/Annotated_files/Annotated_Retained_CBSs.bed
#
######################################### MAIN ################################

# Arguments:

lost_file=$1
retained_file=$2

# Correct number of arguments assessment:

if ! [ $# -eq 2 ]; then
	echo -e "\nError. You must give 2 arguments to this program\n" >&2
	exit 1
fi

# Create the output directory

mkdir ../Results/Potential_loop_file

# 1. Run the python script

./Potential_loop_selection_python.py $lost_file $retained_file

# 2. Sort the resulting bed files:

sort -k1,1 -k2,2n ../Results/Potential_loop_file/potential_lost_loops.bed > ../Results/Potential_loop_file/Potential_lost_loops.bed
sort -k1,1 -k2,2n ../Results/Potential_loop_file/potential_retained_loops.bed > ../Results/Potential_loop_file/Potential_retained_loops.bed

# 3. Remove the unsorted files:

rm ../Results/Potential_loop_file/potential_lost_loops.bed
rm ../Results/Potential_loop_file/potential_retained_loops.bed
