#This Rscript was genereated on 07152020 for the microbiome tutorial series by Dr. Michael Jochum at Baylor College of Medicine
#and was originally designed to be used at the following binder repo:
#https://mybinder.org/v2/gh/MADscientist314/microbiome_binder_repo/master



#Import the libraries
library(phyloseq)
library(microbiome)
library(DirichletMultinomial)
library(reshape2)
library(magrittr)
library(dplyr)
library(microbiome)
library(DESeq2)
library(metacoder)
library(knitr)
library(tibble)

data(atlas1006) 
print(atlas1006)

write_phyloseq(atlas1006, type = "OTU", path = getwd())
write_phyloseq(atlas1006, type = "TAXONOMY", path = getwd())
write_phyloseq(atlas1006, type = "METADATA", path = getwd())


