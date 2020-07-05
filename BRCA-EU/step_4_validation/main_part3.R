# load required libraries
library("lsa")

# read SparseSignatures results
load("processed_data/trinucleotides_counts_wgs.RData")
load("processed_data/trinucleotides_counts_wxs.RData")
load("final_results/SparseSignatures_wgs.RData")
load("final_results/SparseSignatures_wxs.RData")

# estimate goodness of fit

# WGS
goodness_fit = NULL
observed_counts = trinucleotides_counts_wgs
predicted_counts = SparseSignatures_wgs$alpha%*%SparseSignatures_wgs$beta
for(i in rownames(observed_counts)) {
    goodness_fit = c(goodness_fit,cosine(as.numeric(observed_counts[i,]),as.numeric(predicted_counts[i,])))
}
names(goodness_fit) = rownames(observed_counts)
valid_samples_wgs = names(which(goodness_fit>0.85))
save(valid_samples_wgs,file="final_results/valid_samples_wgs.RData")

# WXS
goodness_fit = NULL
observed_counts = trinucleotides_counts_wxs
predicted_counts = SparseSignatures_wxs$alpha%*%SparseSignatures_wxs$beta
for(i in rownames(observed_counts)) {
    goodness_fit = c(goodness_fit,cosine(as.numeric(observed_counts[i,]),as.numeric(predicted_counts[i,])))
}
names(goodness_fit) = rownames(observed_counts)
valid_samples_wxs = names(which(goodness_fit>0.85))
save(valid_samples_wxs,file="final_results/valid_samples_wxs.RData")
