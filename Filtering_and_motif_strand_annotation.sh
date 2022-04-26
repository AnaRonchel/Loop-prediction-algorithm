#!/bin/bash

# Author: Ana Rodriguez Ronchel
#
# Date: 18-11-2020
#
###################################### DESCRIPTION ############################
#
# This script generates a bed file with those CBS regions that have a CTCF 
# motif with a score > 8.7. It adds a 4th column with the strand in which the 
# motif has been found.
# If there is more than one motif within the same CBS and they are in opposite
# strands, the same CBS will be annotated two times to consider both motif
# orientations.
#
# INPUT:
#   - Arg 1: name of the lost CBSs bed file (3 columns: chromosome, start, end)
#   - Arg 2: name of the retained CBSs bed file (same 3 columns)
#   - Arg 3: name of the motif file (from homer, with .motif extension)
#
# OUTPUT:
#   - The same bed files but only with the CBSs with a high motif score and
#     with a fourth column with + and - indicating the strand in which the
#     motif has been found. These files will appear in the "Results" directoy,
#     in a new directory called "Annotated_files". They will have the same name
#     as the input but with the prefix "Annnotated_".
#
# This script is done to work with mm10 genome. To change the organism go to
# the point 3 and change "mm10" for other that must be previously configurated in
# your homer program.
# 
###################################### REQUIREMENTS ###########################
# 
# The directory must be set to the file where this script is located. 
#
# The parent directoy must contain the next directories: 
#     - Data: .bed files and .motif file 
#     - Results: empty
#     - Scripts: this script
#
# Installed programs: awk, Homer (it must be included in the path).
#
##################################### RUN EXAMPLE #############################
#
# ./Filtering_and_motif_strand_annotation.sh Lost_CBSs.bed Retained_CBSs.bed CTCF.motif
#
######################################### MAIN ################################

# Arguments:

lost_file=$1
retained_file=$2
motif_file=$3

# Correct number of arguments asseasment:

if ! [ $# -eq 3 ]; then
	echo -e "\nError. You must give 3 arguments to this program\n" >&2 
	exit 1                 
fi   

# Create the output directoy:

mkdir ../Results/Annotated_files

# 1. Sort the bed files by chromosome (lexically) and by region (nummerically). 

sort -k1,1 -k2,2n ../Data/$lost_file > ../Data/Sorted_$lost_file

sort -k1,1 -k2,2n ../Data/$retained_file > ../Data/Sorted_$retained_file

# 2. Add a 4th column with unique ID (line number), 5th column with a value
#    that will be ignored and a 6th column with the strand. As we don't know
#    the strand at this point it will have * values. This is done to to fulfil
#    the requirements for the input file of Homer program.

awk 'BEGIN{FS= "\t"; OFS = "\t"; line=0} {line+=1; print $0,line,"lost","*"}' ../Data/Sorted_$lost_file > ../Data/6_column_$lost_file

awk 'BEGIN{FS= "\t"; OFS = "\t"; line=0} {line+=1; print $0,line,"retained","*"}' ../Data/Sorted_$retained_file > ../Data/6_column_$retained_file


# 3. Generate a file only with the regions with a CTCF motif score > 8.7 and
#    with the strand info. 

findMotifsGenome.pl ../Data/6_column_$lost_file mm10 ../Results/Annotated_files -find ../Data/$motif_file -size given -mask > ../Results/Annotated_files/findMotifs_$lost_file

findMotifsGenome.pl ../Data/6_column_$retained_file mm10 ../Results/Annotated_files -find ../Data/$motif_file -size given -mask > ../Results/Annotated_files/findMotifs_$retained_file

# 4. Join the files from steps 2 and 3 to generate a new file containing the 
#    coordinates of the region and a 4th column with the motif strand.
# 
# FNR resets back to 1 for the first line of each file but NR keeps on 
# increasing, so NR==FNR is only true with the first file processed
# (findMotifs_file). With this command the first fiel, that contains strand info,
# is processed and, for each region the stand ($5) is stored in a variable "a"
# associated with this region. If the region is repeated in the file because 
# there are more than one motif within it, but the strand is the same,
# nothing happens. If the strand is different, this info is stored in
# another variable "b". Then, the second file is processed and for all the
# regions that were present in the first file (score > 8.7), their coordinates
# are printed ($1, $2, $3) and the strand of the motif too (a[$4]). If the 
# region had another motif orientation, the coordinates and this other strand
# are printed in a new line. 

awk -F'\t' -v OFS='\t' '{ if (NR==FNR && $1 in a){ if (a[$1]==$5){ next }
                                       if (a[$1]!=$5){ b[$1]=$5; next } }
    if (NR==FNR){ a[$1]=$5; next }
    if ($4 in a) { print $1, $2, $3, a[$4] }
    if ($4 in b) { print $1, $2, $3, b[$4] } }' ../Results/Annotated_files/findMotifs_$lost_file ../Data/6_column_$lost_file > ../Results/Annotated_files/Annotated_$lost_file

awk -F'\t' -v OFS='\t' '{ if (NR==FNR && $1 in a){ if (a[$1]==$5){ next }
                                       if (a[$1]!=$5){ b[$1]=$5; next } }
    if (NR==FNR){ a[$1]=$5; next }
    if ($4 in a) { print $1, $2, $3, a[$4] }
    if ($4 in b) { print $1, $2, $3, b[$4] } }' ../Results/Annotated_files/findMotifs_$retained_file ../Data/6_column_$retained_file > ../Results/Annotated_files/Annotated_$retained_file

# 5. Remove the intermediate files:

rm ../Data/Sorted_$lost_file
rm ../Data/Sorted_$retained_file
rm ../Data/6_column_$lost_file
rm ../Data/6_column_$retained_file
#rm ../Results/Annotated_files/findMotifs_$lost_file
#rm ../Results/Annotated_files/findMotifs_$retained_file
#rm ../Results/Annotated_files/motifFindingParameters.txt







