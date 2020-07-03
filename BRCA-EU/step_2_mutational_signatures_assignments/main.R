# load required libraries and sources
library("glmnet")
source("R/mutational.signatures_assignments.R")

# load discovered signatures to be used in the analysis
load("RData/SparseSignatures.RData")
load("RData/trinucleotides_counts.RData")

# perform signatures assignments
set.seed(12345)
SparseSignatures = signaturesAssignment(x=trinucleotides_counts,beta=SparseSignatures$beta)

# save the results
save(SparseSignatures,file="final_results/SparseSignatures.RData")
