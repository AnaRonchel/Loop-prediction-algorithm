#!/usr/bin/env python

# Author: Ana Rodriguez Ronchel
# Date: 18-11-2020
# Description: This script annotates a BED file with the number of genes
#              included in a region (potential loop), their mean expression
#              and a score related to the relative position
#              of that mean expression in a theoretical expression distribution
#              adjusted for that number of genes. This score corresponds to the
#              integral between the mean expression value for that group of genes
#              and infinity. Thus, the most transcriptionally repressed regions 
#              (lower mean expression) will have a score closer to 1 while
#              the most active regions (higher mean expression) will have a
#              score closer to 0.
#
# Input:
#       - Arg1: Gene expression file with gene counts in column 14 and gene 
#               coordinates (chr, start and end) in columns 19-21.
#       - Arg2: 3 column bed file with the potential loops (chr, start, end)
#
# Output:
#       - 6 column bed file with the potential loops coordinates and
#         information about the expression of the genes in that region
#         (chr, start, end, n genes, mean expression, score)
#
# Execution example:
# For lost loops
# ./Loop_expression_annotation.py ../Data/gene_expression.csv ../Results/Potential_loop_file/Potential_lost_loops.bed
# For retained loops
# ./Loop_expression_annotation.py ../Data/gene_expression.csv ../Results/Potential_loop_file/Potential_retained_loops.bed

################################### MAIN ######################################

import sys
import math

# Data treatment
# ==============================================================================
import pandas as pd
import numpy as np

# Graphs
# ==============================================================================
import matplotlib.pyplot as plt
from matplotlib import style

# Distribution adjustment
# ==============================================================================
from scipy import stats
import scipy.integrate as integrate
import scipy.signal as signal
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import KFold

# Matplotlib configuration
# ==============================================================================
style.use('ggplot') or plt.style.use('ggplot')

# Warnings configuration
# ==============================================================================
import warnings
warnings.filterwarnings('ignore')

# Arguments assignment
# ==============================================================================
Genes_exp = sys.argv[1]
Potential_loops = sys.argv[2]




######################## 1. Processing expression data ########################

# Import data
Expression_data = []
with open(Genes_exp, 'r') as genes_file:
    for gene in genes_file:
        Gene=gene.strip().split(";")
        Gene_expression=float(Gene[13].replace(",",".")[0:10])
        Expression_data.append(Gene_expression)

# Convert the list into an array
Expression_data_array = np.array(Expression_data)

# Sort the data
data_sorted = np.sort(Expression_data)

# Histogram to visualize the real data:
#plt.style.use('ggplot')
#fig, ax = plt.subplots(figsize=(6,4))
#ax.hist(data_sorted, bins=100, density=True, color="#3182bd", alpha=0.5)
#ax.plot(data_sorted, np.full_like(data_sorted, -0.001), '|k', markeredgewidth=1)
#ax.set_title('Histograma (normalizado)')
#ax.set_xlabel('x')
#ax.set_ylabel('densidad');





################ 2. Kernel density estimation (KDE) ##################

# It aproximates the density function (sum of functions (kernel) of each 
# observation)

# Cross-validation for kernel and bandwidth identification
# ==============================================================================
param_grid = {'kernel': ['gaussian', 'epanechnikov', 'exponential', 'linear'],
              'bandwidth' : np.linspace(0.01, 3, 10)
             }

grid = GridSearchCV(
        estimator  = KernelDensity(),
        param_grid = param_grid,
        n_jobs     = -1,
        cv         = 10, 
        verbose    = 0
      )

# The result is assigned to _ so that it is not printed on the screen.
_ = grid.fit(X = data_sorted.reshape((-1,1)))

# Best hyperparameters by cross-validation
# ==============================================================================
#print("----------------------------------------")
#print("Best hyperparameters (cv)")
#print("----------------------------------------")
#print(grid.best_params_, ":", grid.best_score_, grid.scoring)

modelo_kde_final = grid.best_estimator_

# Density distribution graphs
# ==============================================================================
X_grid = np.linspace(0, 14000, 100000)
delta = (14000)/99999

# The previous numbers are subject to the corresponding input file

log_densidad_pred = modelo_kde_final.score_samples(X_grid.reshape((-1,1)))
densidad_pred_no_norm = np.exp(log_densidad_pred)

# Correction so that the value of the integral is 1 
result = integrate.simps(densidad_pred_no_norm, X_grid)
densidad_pred = densidad_pred_no_norm/result
result = integrate.simps(densidad_pred, X_grid)

# Illustrative graphs
# ==============================================================================
#fig, ax = plt.subplots(figsize=(7,4))
#ax.hist(data_sorted, bins=100, density=True, color="#3182bd", alpha=0.5)
#ax.plot(data_sorted, np.full_like(data_sorted, -0.001), '|k', markeredgewidth=1)
#ax.plot(X_grid, densidad_pred, color = 'red', label='Kernel: exponential \n bw:3.0')
#ax.set_title('Predicted density')
#ax.set_xlabel('x')
#ax.set_ylabel('Density')
#ax.legend();

# Illustrative example code
# ==============================================================================
#new_X = np.array([1])
#log_density_pred = modelo_kde_final.score_samples(X=new_X.reshape(-1, 1))
#density_pred = np.exp(log_density_pred)
#density_pred
#function = np.vstack([X_grid,densidad_pred])
#result = integrate.simps(densidad_pred, X_grid)
#result





########################### 3. Convolution operation ##########################

# Consecutive convolutions to find the PDF of the total expression of groups of
# genes containing from 1 to 70 genes.
# ==============================================================================

# The array "convoluciones" contains the probability values and the array 
# "grid" has the X-axis values (referring to the expression values) to which
# each of the probabilities contained in the array "convoluciones" refer.
# The size of this array increases with the number of convolutions because
# the convolution operation results in an expanded function with more X-axis 
# values.

convoluciones = np.zeros((70,100000*70))
convoluciones[0,0:len(densidad_pred)] = densidad_pred
grid = np.zeros((70,100000*70))
grid[0,0:len(densidad_pred)] = np.linspace(0,14000, len(densidad_pred))

for i in range(1,convoluciones.shape[0]):
    convoluciones[i,0:((i+1)*len(densidad_pred)-i)] = signal.fftconvolve(convoluciones[(i-1),0:(i*len(densidad_pred)-i+1)], densidad_pred) * delta
    grid[i,0:((i+1)*len(densidad_pred)-i)] = np.linspace(0,14000*(i+1),((i+1)*len(densidad_pred)-i))

 
# Modification of the resulting arrays to reflect the mean expression rather 
# than the total expression
# ==============================================================================
   
convoluciones_mean = np.copy(convoluciones)
grid_mean = np.copy(grid)

for i in range(convoluciones.shape[0]):
    convoluciones_mean[i,:] = convoluciones[i,:]*(i+1)
    grid_mean[i,:] = grid[i,:]/(i+1)


# Illustrative graphs
# ==============================================================================

# Plot the resulting functions to see how the mean expression distribution of
# groups with N genes looks like.
# It increasingly resembles a normal distribution centered on the mean
# distribution of all genes analysed (due to the central limit theorem).

#fig, ax = plt.subplots(figsize=(7,4))
#lista_plot = [0,4,9,19]
#coloreh = ["blue","red","green","black"]
#nombres = ["1 gen","5 genes","10 genes","20 genes"]
#for i in range(len(lista_plot)):
#    ax.plot(grid_mean[lista_plot[i],0:((lista_plot[i]+1)*len(densidad_pred)-lista_plot[i])], convoluciones_mean[lista_plot[i],0:((lista_plot[i]+1)*len(densidad_pred)-lista_plot[i])], color = coloreh[i], label=nombres[i])
#ax.set_title('Densidad predicha')
#ax.set_xlabel('x')
#ax.set_ylabel('densidad')
#plt.xlim([0,200]) 
#ax.legend();

# Illustrative example code
# ==============================================================================

# Code to find out the probability that 3 randomly chosen genes
# have a mean expression of 1000 or more.

#n_genes = 3
#mean_expression = 1000
#num = int((mean_expression)/(delta/n_genes))
#p_value = integrate.simps(convoluciones_mean[(n_genes-1),num:len(densidad_pred)*n_genes-(n_genes-1)],grid_mean[(n_genes-1),num:len(densidad_pred)*n_genes-(n_genes-1)])
#print(p_value)





############################ 4. Region annotation #############################

loops_file = open(Potential_loops, 'r')

output_file = open("../Results/Potential_loop_file_annotated/Loops_expression_annotated.bed", 'w')

for loop in loops_file:
    Loop=loop.strip().split()
    Loop_chr=Loop[0].replace("chr","")
    Loop_start=int(Loop[1])
    Loop_end=int(Loop[2])
    sum_C_expression=0
    n_genes=0
    with open(Genes_exp, 'r') as genes_file:
        for gene in genes_file:
            Gene=gene.strip().split(";")
            Gene_chr=Gene[18]
            Gene_start=int(Gene[19])
            Gene_end=int(Gene[20])
            C_expression=float(Gene[13].replace(",",".")[0:10])
            if Loop_chr == Gene_chr:
                if Loop_start<Gene_start and Loop_end>Gene_end:  #if gene inside loop
                    n_genes+=1
                    sum_C_expression+=C_expression #sum
    if n_genes==0:
        mean_C_expression="NA"
        integral="NA"
    else:
        mean_C_expression=sum_C_expression/n_genes
        num = int((mean_C_expression)/(delta/n_genes))
        integral = integrate.simps(convoluciones_mean[(n_genes-1),num:len(densidad_pred)*n_genes-(n_genes-1)],grid_mean[(n_genes-1),num:len(densidad_pred)*n_genes-(n_genes-1)])
    
    output_file.write(loop.strip())
    output_file.write("\t")
    output_file.write(str(n_genes))
    output_file.write("\t")
    output_file.write(str(mean_C_expression))
    output_file.write("\t")
    output_file.write(str(integral))
    output_file.write("\n")
    
loops_file.close()
output_file.close()
