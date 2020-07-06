# load required libraries and sources
library("glmnet")
source("R/mutational.signatures_assignments.R")
library("lsa")

# load discovered signatures to be used in the analysis
load("RData/SparseSignatures.RData")
load("RData/trinucleotides_counts.RData")

# perform signatures assignments
set.seed(54321)
SparseSignatures = signaturesAssignment(x=trinucleotides_counts,beta=SparseSignatures$beta)

# estimate goodness of fit
predicted_counts = SparseSignatures$alpha%*%SparseSignatures$beta
goodness_fit = NULL
for(i in rownames(SparseSignatures$alpha)) {
    goodness_fit = c(goodness_fit,as.numeric(cosine(as.numeric(trinucleotides_counts[i,]),as.numeric(predicted_counts[i,]))))
}
names(goodness_fit) = rownames(SparseSignatures$alpha)
SparseSignatures$goodness_fit = goodness_fit

# save the results
save(SparseSignatures,file="final_results/SparseSignatures.RData")
