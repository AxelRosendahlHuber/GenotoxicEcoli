#=== UTILS ===# 
# Axel Rosendahl Huber
# Purpose: Load the commonly used packages and created functions for the Genotoxic E.coli study

# ==== Miscealaneous ===== # 
options(stringsAsFactors = F) # Prevent conversion of strings to factors in dataframe
library(tidyverse)
library(MutationalPatterns)
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)
library(reshape2)
library(ggdendro)
library(ggplot2)
library(tidyverse)
library(magrittr)
library(VariantAnnotation)

# ======== Retrieve and source all functions ===== # 
wd <- getwd()                       
setwd(paste(wd, "/Functions", sep = ""))

# Get all available functions in the Scripts/R_functions folder
functions <- dir(pattern = ".R")
for (i in functions) {
  source(i)
  cat("Loading function:", i, "\n" ,  sep = " ")
}
setwd(wd)

ref_genome <- "BSgenome.Hsapiens.UCSC.hg19"

