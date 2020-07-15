#This Rscript was genereated on 07152020 for the microbiome tutorial series by Dr. Michael Jochum at Baylor College of Medicine
#and was originally designed to be used at the following binder repo:
#https://mybinder.org/v2/gh/MADscientist314/microbiome_binder_repo/master

setwd("D:/github/microbiome_binder_repo")

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


# Note that this particular approach will be super slow.
# And take just as long everytime you edit your code
library(holepunch)
write_install() # Writes install.R with all your dependencies
write_runtime() # Writes the date your code was last modified. Can be overridden.
generate_badge() # Generates a badge you can add to your README. Clicking badge will launch the Binder.
# ----------------------------------------------
# At this time ???? push the code to GitHub ????
# ----------------------------------------------
# Then click the badge on your README or run
build_binder() # to kick off the build process
# ????????
data(atlas1006) 
print(atlas1006)

write_phyloseq(atlas1006, type = "OTU", path = getwd())
write_phyloseq(atlas1006, type = "TAXONOMY", path = getwd())
write_phyloseq(atlas1006, type = "METADATA", path = getwd())


