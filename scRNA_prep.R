#Preparing the environment with the necessary packages, colours and fontsizes
#Created by Markus Boesch

set.seed(1234567)
library(Seurat)
library(SeuratObject)
library(dplyr)
library(ggplot2)
library(cowplot)
library(clustree)
library(patchwork)
library(pheatmap)
library(Matrix)
library(readxl)
library(stringr)
library(SingleCellExperiment)
library(RColorBrewer)
library(biomaRt)
library(openxlsx)
library(ggpubr)
library("xlsx")
library(openxlsx)
library(gprofiler2)
library(ggpubr)
library(limma)
library(EnhancedVolcano)
library(miloR)
library(scater)

#Fontsize
fontsize<-14.5

##colours
colours_disease<-c("palegreen4", "steelblue3", "salmon3")
colours_ProgRec<-c("palegreen4", "steelblue3", "sienna3", "red4" )
colours_PBMC<-c("lightpink3" , "slategray3","mediumpurple2","navy","#00B6EB","turquoise4", "burlywood","#FB61D7")
colours_Mon<-c("lightpink3","rosybrown2","lightsteelblue3","burlywood")
colours_Tcells<-c("slategray3","navy","lightpink3","darkgreen","darkmagenta","mediumorchid2","slateblue2","steelblue4","#00B6EB")
colours_CD4<-c("slategray3","navy","lightpink3","darkgreen")
colours_CD8MAIT<-c("darkmagenta","mediumorchid2","slateblue2","steelblue4","#00B6EB")
colours_Bcell<-c("darkgreen","salmon1")
colours_NK<-c("dodgerblue3","orchid3")
colours_allsub<-c("lightpink3","rosybrown2","lightsteelblue3","burlywood","slategray3","navy","lightpink3","darkgreen","darkmagenta","mediumorchid2","slateblue2","steelblue4","#00B6EB","dodgerblue3","orchid3","darkgreen","salmon1","#FB61D7")