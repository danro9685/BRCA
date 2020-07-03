# load required libraries and sources
library("NMF")
library("nnls")
library("nnlasso")
library("parallel")
source("R/signatures.discovery.lasso.R")

# load the data
load("final_data/trinucleotides_counts.RData")
trinucleotides_counts = ((trinucleotides_counts/rowSums(trinucleotides_counts))*2500) # normalize trinucleotides counts
load("final_data/background.RData")

# settings
K = 2:15
nmf_runs = 1000
my_seed_starting_beta = 12345

# fit the initial betas for each configuration
initial_betas = startingBetaEstimation(x=trinucleotides_counts,K=K,background_signature=background,normalize_counts=FALSE,nmf_runs=nmf_runs,seed=my_seed_starting_beta,verbose=TRUE)
save(initial_betas,file="RData/initial_betas.RData")
