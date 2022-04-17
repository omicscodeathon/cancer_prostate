
install_github("TransBioInfoLab/pathwayMultiomics")
library(devtools)
library(tidyverse)
library(pathwayMultiomics)

alzheimersMultiOmics_df
adMiniMax_df <- 
  alzheimersMultiOmics_df %>% 
  dplyr::select(pathway, ends_with("Pval")) %>% 
  rename(
    SNP = snpPval, DNAm = dnamPval, RNAseq = rnaseqPval
  ) %>% 
  MiniMax()
adMiniMax_df
