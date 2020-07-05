# load required libraries and sources
library("glmnet")
source("R/mutational.signatures_assignments.R")

# load discovered signatures to be used in the analysis
load("RData/SparseSignatures.RData")
load("processed_data/trinucleotides_counts_wgs.RData")
load("processed_data/trinucleotides_counts_wxs.RData")

# perform signatures assignments for WGS validation set
set.seed(12345)
SparseSignatures_wgs = signaturesAssignment(x=trinucleotides_counts_wgs,beta=SparseSignatures$beta)

# save the results
save(SparseSignatures_wgs,file="final_results/SparseSignatures_wgs.RData")

# perform signatures assignments for WXS validation set
set.seed(54321)
SparseSignatures_wxs = signaturesAssignment(x=trinucleotides_counts_wxs,beta=SparseSignatures$beta)

# save the results
save(SparseSignatures_wxs,file="final_results/SparseSignatures_wxs.RData")
