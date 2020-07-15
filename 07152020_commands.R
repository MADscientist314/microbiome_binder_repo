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


