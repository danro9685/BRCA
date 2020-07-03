# load required libraries and sources
library("nnls")

# load the data
load("final_data/trinucleotides_counts.RData")
trinucleotides_counts = ((trinucleotides_counts/rowSums(trinucleotides_counts))*2500) # normalize trinucleotides counts
load("final_data/background.RData")
load("RData/cross_validation.RData")

# data structure to save final results
cross_validation_summary = array(NA,c(15,2))
rownames(cross_validation_summary) = paste0("Rank ",1:15)
colnames(cross_validation_summary) = c("Median_CV_Error","Best_Lambda")

# set the seed
set.seed(3764302)

# performing cross-validation for rank 1 solution (only background)
grid_search_mse = array(NA,c(1000,1))
rownames(grid_search_mse) = paste0("Iterarion_",1:1000)
colnames(grid_search_mse) = "LASSO_0.0"
for(i in 1:1000) {
    load(paste0("cv_entries/cv_entries",i,".RData"))
    # perform estimation
    x_cv = NULL
    for(j in 1:5) {
        if(j==1) {
            x_cv = trinucleotides_counts
            x_cv[cv_entries] = 0
        }
        else {
            predicted_counts = curr_iter_alpha %*% curr_iter_beta
            x_cv[cv_entries] = predicted_counts[cv_entries]
        }
        # perform the inference
        inference_beta = matrix(background,nrow=1)
        rownames(inference_beta) = "Background"
        colnames(inference_beta) = colnames(trinucleotides_counts)
        inference_alpha = array(NA,c(nrow(trinucleotides_counts),1))
        rownames(inference_alpha) = rownames(trinucleotides_counts)
        colnames(inference_alpha) = "Background"
        for(k in 1:nrow(trinucleotides_counts)) {
            inference_alpha[k,] <- nnls(t(inference_beta),as.vector(x_cv[k,]))$x
        }
        curr_iter_alpha = inference_alpha
        curr_iter_beta = inference_beta
    }
    # save cross validation error
    curr_predicted_counts = round(curr_iter_alpha%*%curr_iter_beta)
    curr_true_considered_counts = as.vector(trinucleotides_counts[cv_entries])
    curr_predicted_considered_counts = as.vector(curr_predicted_counts[cv_entries])
    error = mean((curr_true_considered_counts-curr_predicted_considered_counts)^2)
    grid_search_mse[paste0("Iterarion_",i),"LASSO_0.0"] = error
}

# save final results
cross_validation_summary[1,"Median_CV_Error"] = median(grid_search_mse[,"LASSO_0.0"])
cross_validation_summary[1,"Best_Lambda"] = 0.000
for(i in 1:14) {
    cv_medians = NULL
    for(j in names(cross_validation$grid_search_mse[,,i])) {
        cv_medians = c(cv_medians,median(cross_validation$grid_search_mse[,j,i][[1]]))
    }
    cross_validation_summary[(i+1),"Median_CV_Error"] = min(cv_medians)
    cross_validation_summary[(i+1),"Best_Lambda"] = as.numeric(gsub("_Lambda_Beta","",names(cross_validation$grid_search_mse[,,i])[which(cv_medians==min(cv_medians))[1]]))
}
save(cross_validation_summary,file="RData/cross_validation_summary.RData")
