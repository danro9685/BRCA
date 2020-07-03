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
load("RData/initial_betas.RData")

# settings
K = 2:15
lambda_values_alpha = c(0.000)
lambda_values_beta = c(0.000,0.025,0.050,0.100)
cross_validation_entries = 0.01
cross_validation_repetitions = 1000
num_processes = 100
my_seed_cv = 54321
log_file = ""

# performing cross-validation for a grid search to estimate the optimal number of signatures
cross_validation = nmfLassoCV(x=trinucleotides_counts,K=K,starting_beta=initial_betas,normalize_counts=FALSE,lambda_values_alpha=lambda_values_alpha,lambda_values_beta=lambda_values_beta,cross_validation_entries=cross_validation_entries,cross_validation_repetitions=cross_validation_repetitions,num_processes=num_processes,seed=my_seed_cv,verbose=TRUE,log_file=log_file)
save(cross_validation,file="RData/cross_validation.RData")
