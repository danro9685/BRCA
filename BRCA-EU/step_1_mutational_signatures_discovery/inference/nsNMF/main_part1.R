# load required libraries and sources
library("NMF")
library("nnls")

# load the data
load("final_data/trinucleotides_counts.RData")
trinucleotides_counts = ((trinucleotides_counts/rowSums(trinucleotides_counts))*2500) # normalize trinucleotides counts

# settings
K = 1:15
nmf_runs = 1000
my_seed = 22222

# perform the inference
res_nmf = nmf(x=t(trinucleotides_counts),rank=K,method="nsNMF",nrun=nmf_runs,seed=my_seed)
beta = array(list(),c(length(K),1))
rownames(beta) = paste0("Rank ",K)
colnames(beta) = "Results"
alpha = beta
for(k in 1:length(res_nmf$fit)) {
    curr_beta = t(basis(res_nmf$fit[[k]]))
    rownames(curr_beta) = paste0("Signature ",1:K[k])
    beta[[k,1]] = curr_beta/rowSums(curr_beta)
    curr_alpha = array(NA,c(nrow(trinucleotides_counts),nrow(curr_beta)))
    rownames(curr_alpha) = rownames(trinucleotides_counts)
    colnames(curr_alpha) = paste0("Signature ",1:K[k])
    for(j in 1:nrow(curr_alpha)) {
        curr_alpha[j,] = nnls(t(beta[[k,1]]),as.vector(trinucleotides_counts[j,]))$x
    }
    alpha[[k,1]] = curr_alpha
}
goodness_fit = res_nmf$measures

# save the results
nmf_nsNMF = list(alpha=alpha,beta=beta,goodness_fit=goodness_fit)
save(nmf_nsNMF,file="nmf_nsNMF.RData")
