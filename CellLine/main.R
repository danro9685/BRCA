# load required libraries and sources
library("SparseSignatures")
library("data.table")
library("ggplot2")
library("gridExtra")
"signatures.plot.v2" = function (beta, useColNames = TRUE, mutation_categories = NULL, firstBackground = TRUE, xlabels = TRUE) {
    if (firstBackground) {
        rownames(beta) <- c("Background", paste0("Signature ", 
            1:(nrow(beta) - 1)))
    }
    else {
        rownames(beta) <- rownames(beta)
    }
    if (!useColNames) {
        colnames(beta) <- sort(mutation_categories$cat)
    }
    x <- as.data.table(melt(beta, varnames = c("signature", "cat")))
    x[, `:=`(Context, paste0(substr(cat,1,1), ".", substr(cat, 
        7, 7)))]
    x[, `:=`(alt, paste0(substr(cat,3,3), ">", substr(cat, 
        5, 5)))]
    glist <- list()
    for (i in 1:nrow(beta)) {
        plt <- ggplot(x[signature == rownames(beta)[i]]) + geom_bar(aes(x = Context, 
            y = value, fill = alt), stat = "identity", position = "identity") + 
            facet_wrap(~alt, nrow = 1, scales = "free_x") + theme(axis.text.x = element_text(angle = 90, 
            hjust = 1), panel.background = element_blank(), axis.line = element_line(colour = "black")) + 
            ggtitle(rownames(beta)[i]) + theme(legend.position = "none") + 
            ylab("Frequency of mutations")
        if (!xlabels) {
            plt <- plt + theme(axis.text.x = element_blank(), 
                axis.ticks.x = element_blank())
        }
        glist[[i]] <- plt
    }
    grid.arrange(grobs = glist, ncol = ceiling(nrow(beta)/3))
}
library("glmnet")
source("R/mutational.signatures_assignments.R")
library("lsa")

# load discovered signatures to be used in the analysis
load("raw_data/trinucleotides.RData")
load("RData/SparseSignatures.RData")

# read data
MCF7_T4 = read.table("raw_data/MCF7_T4_SNV_hg38_context.tsv",header=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
MCF7_A3BKO_T4 = read.table("raw_data/MCF7_A3BKO_T4_SNV_hg38_context.tsv",header=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
MDA_231_T4 = read.table("raw_data/MDA_231_T4_SNV_hg38_context.tsv",header=TRUE,check.names=FALSE,stringsAsFactors=FALSE)
trinucleotides_counts = rbind(MCF7_T4,MCF7_A3BKO_T4,MDA_231_T4)
rownames(trinucleotides_counts) = c("MCF7_T4","MCF7_A3BKO_T4","MDA_231_T4")
trinucleotides_counts = as.matrix(trinucleotides_counts[,trinucleotides])

# perform signatures assignments
set.seed(12345)
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

# make plots
data = rbind(trinucleotides_counts["MCF7_T4",,drop=FALSE],predicted_counts["MCF7_T4",,drop=FALSE],trinucleotides_counts["MCF7_A3BKO_T4",,drop=FALSE],predicted_counts["MCF7_A3BKO_T4",,drop=FALSE],trinucleotides_counts["MDA_231_T4",,drop=FALSE],predicted_counts["MDA_231_T4",,drop=FALSE])
rownames(data) = c("MCF7_T4 (Observations) - Cosine similarity 0.98","MCF7_T4 (Predictions) - Cosine similarity 0.98","MCF7_A3BKO_T4 (Observations) - Cosine similarity 0.98","MCF7_A3BKO_T4 (Predictions) - Cosine similarity 0.98","MDA_231_T4 (Observations) - Cosine similarity 0.97","MDA_231_T4 (Predictions) - Cosine similarity 0.97")
signatures.plot.v2(data,useColNames=TRUE,firstBackground=FALSE,xlabels=FALSE)

# process data
normalized_alpha = SparseSignatures$alpha/rowSums(SparseSignatures$alpha)
rownames(normalized_alpha) = c("cl1 - MCF7_T4","cl2 - MCF7_A3BKO_T4","cl3 - MDA_231_T4")
colnames(normalized_alpha) = c("SPS1","SPS2","SPS3","SPS4","SPS5","SPS6","SPS7","SPS8")

# make plot
EXPOSURE = NULL
SIGNATURE = NULL
SAMPLE = NULL
for(i in 1:nrow(normalized_alpha)) {
    for(j in 1:ncol(normalized_alpha)) {
        EXPOSURE = c(EXPOSURE,normalized_alpha[i,j])
        SIGNATURE = c(SIGNATURE,colnames(normalized_alpha)[j])
        SAMPLE = c(SAMPLE,rownames(normalized_alpha)[i])
    }
}
ggplot(data=data.frame(EXPOSURE=EXPOSURE,SIGNATURE=SIGNATURE,SAMPLE=SAMPLE),aes(x=SAMPLE,y=EXPOSURE,fill=SIGNATURE)) + geom_bar(stat="identity") + scale_fill_brewer(palette="Spectral") + ggtitle("Cell line dataset - Signatures Assignments")
