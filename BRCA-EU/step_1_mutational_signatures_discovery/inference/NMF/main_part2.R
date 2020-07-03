# load required libraries and sources
library("nnls")

# load data and results
load("final_data/trinucleotides_counts.RData")
load("nmf_brunet.RData")

# set the seed
set.seed(77777)

# fit final results
alpha = nmf_brunet$alpha
beta = nmf_brunet$beta
rank = as.matrix(nmf_brunet$goodness_fit[,c("evar","cophenetic","dispersion","silhouette.consensus")])
rownames(rank) = 1:nrow(rank)
colnames(rank) = c("Explained_Variance","Cophenetic_Coefficient","Dispersion_Coefficient","Silhouette_Consensus")
rank[1,"Cophenetic_Coefficient"] = 1
rank[1,"Silhouette_Consensus"] = 1
for(k in 1:15) {
    curr_beta = beta[[k]]
    curr_alpha = alpha[[k]]
    for(j in 1:nrow(curr_alpha)) {
        curr_alpha[j,] = nnls(t(curr_beta),as.vector(trinucleotides_counts[j,]))$x
    }
    alpha[[k]] = curr_alpha
}

# save results
NMF = list(alpha=alpha,beta=beta,rank=rank)
save(NMF,file="final_results/NMF.RData")
