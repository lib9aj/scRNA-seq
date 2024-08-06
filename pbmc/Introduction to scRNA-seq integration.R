
###Introduction to scRNA-seq integration

##Integration goals
#1.Identify cell subpopulations that are present in both datasets
#2.Obtain cell type markers that are conserved in both control and stimulated cells
#3.Compare the datasets to find cell-type specific responses to stimulation

#Setup the Seurat objects
library(Seurat)
library(SeuratData)
library(patchwork)

# install dataset

options(timeout = 300)  
InstallData("ifnb")
#The object contains data from human PBMC from two conditions, interferon-stimulated and control cells (stored in the stim column in the object metadata).
#We will aim to integrate the two conditions together, so that we can jointly identify cell subpopulations across datasets, and then explore how each group differs across conditions
#In Seurat v5, we keep all the data in one object, but simply split it into multiple ‘layers’. 