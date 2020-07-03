# load required libraries and sources
library("nnls")
source("R/SignatureAnalyzer.PCAWG.function.R")

# load the data
load("final_data/trinucleotides_counts.RData")
trinucleotides_counts = ((trinucleotides_counts/rowSums(trinucleotides_counts))*2500) # normalize trinucleotides counts
trinucleotides_counts = as.matrix(t(trinucleotides_counts))

# settings
K_max = 15
repetitions = 1000
my_seed = 33333

# perform the inference
inference_results = list()
posteriors_runs = NULL
for(multiple_runs in 1:repetitions) {
    set.seed(multiple_runs)
    res = BayesNMF.L1W.L2H(trinucleotides_counts,200000,10,5,1.e-05,K_max,K_max,1)
    alpha = t(res[[2]])
    alpha = alpha[,which(colSums(res[[1]])>0)]
    colnames(alpha) = paste0("Signature_",1:ncol(alpha))
    beta = t(res[[1]])
    beta = beta[which(colSums(res[[1]])>0),]
    posterior = res[[4]]
    rownames(beta) = colnames(alpha)
    inference_results[[multiple_runs]] = list(alpha=alpha,beta=beta,posterior=posterior)
    posteriors_runs = c(posteriors_runs,posterior)
}

# save optimal set of signatures
best_solution = inference_results[[which(posteriors_runs==min(posteriors_runs))]]
beta = best_solution$beta/rowSums(best_solution$beta)
alpha = array(NA,c(ncol(trinucleotides_counts),nrow(beta)))
rownames(alpha) = colnames(trinucleotides_counts)
colnames(alpha) = rownames(beta)
for(j in 1:nrow(alpha)) {
    alpha[j,] = nnls(t(beta),as.vector(trinucleotides_counts[,j]))$x
}

# save the results
SignatureAnalyzer = list(inference_results=inference_results,best_solution=list(alpha=alpha,beta=beta))
save(SignatureAnalyzer,file="SignatureAnalyzer.RData")
