# Author: Ana Rodriguez Ronchel
#
# Date: 18-11-2020
#
################################## DESCRIPTION ################################
#
# This script process the BED files generated with the python script 
# "Loop_expression_annotation.py" and select the regions (potential loops) with
# a mean expression value in the top 10% (score <= 0.1) or in the bottom 10% 
# (score >= 0.9) of the corresponding Probability Density Function (PDF) based
# on the number of genes within the region.
#
################################## REQUIREMENTS ###############################
#
# The directory must be set to the file where this script is located.
#
# The parent directory must contain the next directories:
#     - Results/Potential_loop_file_annotated: Annotated lost and retained
#       potential loops (BED files).
#
#################################### MAIN #####################################

# LOAD DATA

df_retained <- read.csv("../Results/Potential_loop_file_annotated/Retained_loops_expression_annotated.bed", sep="\t", dec = ",", header = TRUE)
names(df_retained) <- c("Chr","Start","End","n_genes","Expression","Score")
df_lost <- read.csv("../Results/Potential_loop_file_annotated/Lost_loops_expression_annotated.bed", sep="\t", dec = ",", header = TRUE)
names(df_lost) <- c("Chr","Start","End","n_genes","Expression","Score")

# SELECT REGIONS

df_retained_select <- df_retained[which(df_retained$Score <= 0.1 | df_retained$Score >= 0.9), ]

df_retained_select_active <- df_retained_select[which(df_retained_select$Score <= 0.1), ]
df_retained_select_repressed <- df_retained_select[which(df_retained_select$Score >= 0.9), ]

df_lost_select <- df_lost[which(df_lost$Score <= 0.1 | df_lost$Score >= 0.9), ]

df_lost_select_active <- df_lost_select[which(df_lost_select$Score <= 0.1), ]
df_lost_select_repressed <- df_lost_select[which(df_lost_select$Score >= 0.9), ]

# WRITE FILES

write.table(df_retained_select_active[,0:3], file="../Results/Selected_loops/Retained_loops_active.bed", sep="\t", quote=F, row.names=F, col.names=F)
write.table(df_retained_select_repressed[,0:3], file="../Results/Selected_loops/Retained_loops_repressed.bed", sep="\t", quote=F, row.names=F, col.names=F)
write.table(df_lost_select_active[,0:3], file="../Results/Selected_loops/Lost_loops_active.bed", sep="\t", quote=F, row.names=F, col.names=F)
write.table(df_lost_select_repressed[,0:3], file="../Results/Selected_loops/Lost_loops_repressed.bed", sep="\t", quote=F, row.names=F, col.names=F)
  
