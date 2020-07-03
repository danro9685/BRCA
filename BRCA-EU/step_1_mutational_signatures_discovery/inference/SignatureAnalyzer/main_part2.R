# load required libraries and sources
library("nnls")
library("lsa")

# load data and results
load("final_data/trinucleotides_counts.RData")
load("SignatureAnalyzer.RData")

# set the seed
set.seed(99999)

# fit final results
alpha_best = SignatureAnalyzer$best_solution$alpha
beta_best = SignatureAnalyzer$best_solution$beta
beta = array(list(),c(ncol(alpha_best),1))
rownames(beta) = paste0("Rank ",1:ncol(alpha_best))
colnames(beta) = "Results"
alpha = beta
goodness_fit_mean = NULL
goodness_fit_percent = NULL
valid = 1:ncol(alpha_best)
starting = NULL
for(i in 1:ncol(alpha_best)) {
    # consider each valid signature
    mean_fit = NULL
    for(j in valid) {
        # consider signature j
        curr_beta = rbind(starting,beta_best[j,,drop=FALSE])
        curr_alpha = array(NA,c(nrow(alpha_best),nrow(curr_beta)))
        rownames(curr_alpha) = rownames(trinucleotides_counts)
        # perform fit for alpha
        for(k in 1:nrow(curr_alpha)) {
            curr_alpha[k,] = nnls(t(curr_beta),as.vector(trinucleotides_counts[k,]))$x
        }
        colnames(curr_alpha) = paste0("Signature ",1:ncol(curr_alpha))
        rownames(curr_beta) = colnames(curr_alpha)
        # estimate goodness of fit
        predicted_counts = curr_alpha%*%curr_beta
        curr_goodness_fit = NULL
        for(k in 1:nrow(trinucleotides_counts)) {
            curr_goodness_fit = c(curr_goodness_fit,as.numeric(cosine(as.numeric(trinucleotides_counts[k,]),as.numeric(predicted_counts[k,]))))
        }
        mean_fit = c(mean_fit,mean(curr_goodness_fit,na.rm=TRUE))
    }
    # get signature with best predictive power (i.e., explaining most mutations)
    curr_best = valid[which(mean_fit==max(mean_fit))[1]]
    starting = rbind(starting,beta_best[curr_best,,drop=FALSE])
    valid = valid[which(valid!=curr_best)]
    # perform final fit
    curr_beta = starting
    curr_alpha = array(NA,c(nrow(alpha_best),nrow(curr_beta)))
    rownames(curr_alpha) = rownames(trinucleotides_counts)
    for(k in 1:nrow(curr_alpha)) {
        curr_alpha[k,] = nnls(t(curr_beta),as.vector(trinucleotides_counts[k,]))$x
    }
    colnames(curr_alpha) = paste0("Signature ",1:ncol(curr_alpha))
    rownames(curr_beta) = colnames(curr_alpha)
    alpha[[i,"Results"]] = curr_alpha
    beta[[i,"Results"]] = curr_beta
    # estimate goodness of fit
    predicted_counts = curr_alpha%*%curr_beta
    curr_goodness_fit = NULL
    for(k in 1:nrow(trinucleotides_counts)) {
        curr_goodness_fit = c(curr_goodness_fit,as.numeric(cosine(as.numeric(trinucleotides_counts[k,]),as.numeric(predicted_counts[k,]))))
    }
    goodness_fit_mean = c(goodness_fit_mean,mean(curr_goodness_fit))
    goodness_fit_percent = c(goodness_fit_percent,(length(which(curr_goodness_fit>0.95))/length(curr_goodness_fit)))
}
rank_estimation = cbind(goodness_fit_mean,goodness_fit_percent)
rownames(rank_estimation) = paste0("Rank ",1:nrow(rank_estimation))
colnames(rank_estimation) = c("Goodness of fit (mean cosine similarity)","Goodness of fit (percent >0.95 cosine similarity)")

# save results
SignatureAnalyzer = list(alpha=alpha,beta=beta,rank=rank_estimation)
save(SignatureAnalyzer,file="final_results/SignatureAnalyzer.RData")
