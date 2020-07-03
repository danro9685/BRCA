# load required libraries
library("SparseSignatures")
library("data.table")
library("ggplot2")
library("gridExtra")
library("lsa")
"signatures.plot.v2" = function (beta, useColNames = TRUE, mutation_categories = NULL, firstBackground = TRUE, xlabels = TRUE) {
    if (firstBackground) {
        rownames(beta) <- c("Background", paste0("Signature ", 
            1:(nrow(beta) - 1)))
    }
    else {
        ###rownames(beta) <- paste0("Signature ", 1:nrow(beta))
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

# read the results
load("final_data/trinucleotides_counts.RData")
load("results_inference/NMF.RData")
load("results_inference/nsNMF.RData")
load("results_inference/SignatureAnalyzer.RData")
load("results_inference/SigProfiler.RData")
load("results_inference/SparseSignatures.RData")

# read the signatures by COSMICv3.1
load("COSMIC/COSMIC_v3.1_June2020.RData")
cosmic_signatures = COSMIC_v3.1_June2020[,colnames(trinucleotides_counts)]

# make plots of optimal number of signatures
RANGE = 1:15
VALUE = as.numeric(NMF$rank[,"Silhouette_Consensus"])
ggplot(data=data.frame(VALUE=VALUE,RANGE=RANGE), aes(x=RANGE,y=VALUE,group=1)) + geom_line(linetype="dashed") + geom_point() + scale_x_continuous(breaks=RANGE,labels=as.character(RANGE),limits=c(min(RANGE),max(RANGE))) + ggtitle("NMF (Silhouette Consensus)") + xlab("Number of Signatures") + ylab("Mean Silhouette Consensus (1000 restarts)")

RANGE = 1:15
VALUE = as.numeric(NMF$rank[,"Explained_Variance"])
ggplot(data=data.frame(VALUE=VALUE,RANGE=RANGE), aes(x=RANGE,y=VALUE,group=1)) + geom_line(linetype="dashed") + geom_point() + scale_x_continuous(breaks=RANGE,labels=as.character(RANGE),limits=c(min(RANGE),max(RANGE))) + ggtitle("NMF (Explained Variance)") + xlab("Number of Signatures") + ylab("Mean Explained Variance (1000 restarts)")

RANGE = 1:15
VALUE = as.numeric(nsNMF$rank[,"Silhouette_Consensus"])
ggplot(data=data.frame(VALUE=VALUE,RANGE=RANGE), aes(x=RANGE,y=VALUE,group=1)) + geom_line(linetype="dashed") + geom_point() + scale_x_continuous(breaks=RANGE,labels=as.character(RANGE),limits=c(min(RANGE),max(RANGE))) + ggtitle("nsNMF (Silhouette Consensus)") + xlab("Number of Signatures") + ylab("Mean Silhouette Consensus (1000 restarts)")

RANGE = 1:15
VALUE = as.numeric(nsNMF$rank[,"Explained_Variance"])
ggplot(data=data.frame(VALUE=VALUE,RANGE=RANGE), aes(x=RANGE,y=VALUE,group=1)) + geom_line(linetype="dashed") + geom_point() + scale_x_continuous(breaks=RANGE,labels=as.character(RANGE),limits=c(min(RANGE),max(RANGE))) + ggtitle("nsNMF (Explained Variance)") + xlab("Number of Signatures") + ylab("Mean Explained Variance (1000 restarts)")

RANGE = 1:12
VALUE = as.numeric(SignatureAnalyzer$rank[,"Goodness of fit (mean cosine similarity)"])
ggplot(data=data.frame(VALUE=VALUE,RANGE=RANGE), aes(x=RANGE,y=VALUE,group=1)) + geom_line(linetype="dashed") + geom_point() + scale_x_continuous(breaks=RANGE,labels=as.character(RANGE),limits=c(min(RANGE),max(RANGE))) + ggtitle("SignatureAnalyzer (Mean Cosine Similarity)") + xlab("Number of Signatures") + ylab("Mean Cosine Similarity")

RANGE = 1:12
VALUE = as.numeric(SignatureAnalyzer$rank[,"Goodness of fit (percent >0.95 cosine similarity)"])
ggplot(data=data.frame(VALUE=VALUE,RANGE=RANGE), aes(x=RANGE,y=VALUE,group=1)) + geom_line(linetype="dashed") + geom_point() + scale_x_continuous(breaks=RANGE,labels=as.character(RANGE),limits=c(min(RANGE),max(RANGE))) + ggtitle("SignatureAnalyzer (Goodness of fit)") + xlab("Number of Signatures") + ylab("Goodness of fit (percentage of predictions with >0.95 cosine similarity)")

RANGE = 1:15
VALUE = as.numeric(SigProfiler$rank[,"Stability"])
ggplot(data=data.frame(VALUE=VALUE,RANGE=RANGE), aes(x=RANGE,y=VALUE,group=1)) + geom_line(linetype="dashed") + geom_point() + scale_x_continuous(breaks=RANGE,labels=as.character(RANGE),limits=c(min(RANGE),max(RANGE))) + ggtitle("SigProfiler (Stability)") + xlab("Number of Signatures") + ylab("Mean Stability (1000 restarts)")

RANGE = 1:15
VALUE = as.numeric(SigProfiler$rank[,"MeanL2"])
ggplot(data=data.frame(VALUE=VALUE,RANGE=RANGE), aes(x=RANGE,y=VALUE,group=1)) + geom_line(linetype="dashed") + geom_point() + scale_x_continuous(breaks=RANGE,labels=as.character(RANGE),limits=c(min(RANGE),max(RANGE))) + ggtitle("SigProfiler (L2 Norm)") + xlab("Number of Signatures") + ylab("Mean L2 Norm (1000 restarts)")

RANGE = 1:15
VALUE = as.numeric(SparseSignatures$rank[,"Median cross-validation error"])
ggplot(data=data.frame(VALUE=VALUE,RANGE=RANGE), aes(x=RANGE,y=VALUE,group=1)) + geom_line(linetype="dashed") + geom_point() + scale_x_continuous(breaks=RANGE,labels=as.character(RANGE),limits=c(min(RANGE),max(RANGE))) + ggtitle("SparseSignatures (Cross Validation Error)") + xlab("Number of Signatures") + ylab("Median Cross Validation Error (1000 restarts)")

RANGE = 1:15
VALUE = as.numeric(SparseSignatures$rank[,"Goodness of fit (percent >0.95 cosine similarity)"])
ggplot(data=data.frame(VALUE=VALUE,RANGE=RANGE), aes(x=RANGE,y=VALUE,group=1)) + geom_line(linetype="dashed") + geom_point() + scale_x_continuous(breaks=RANGE,labels=as.character(RANGE),limits=c(min(RANGE),max(RANGE))) + ggtitle("SparseSignatures (Goodness of fit)") + xlab("Number of Signatures") + ylab("Goodness of fit (percentage of predictions with >0.95 cosine similarity)")

# NMF
alpha = NMF$alpha["Rank 8",][[1]]
beta = NMF$beta["Rank 8",][[1]]
silhouette_consensus = as.numeric(NMF$rank[,"Silhouette_Consensus"])
names(silhouette_consensus) = paste0("Signature ",1:15)
explained_variance = as.numeric(NMF$rank[,"Explained_Variance"])
names(explained_variance) = paste0("Signature ",1:15)
mapping_val  = NULL
mapping_cor  = NULL
for(i in 1:nrow(beta)) {
    tmp = NULL
    for(j in 1:nrow(cosmic_signatures)) {
        tmp = c(tmp,as.numeric(cosine(as.numeric(beta[i,]),as.numeric(cosmic_signatures[j,]))))
    }
    mapping_val = c(mapping_val,rownames(cosmic_signatures)[which(tmp==max(tmp))][1])
    mapping_cor = c(mapping_cor,max(tmp))
}
indeces = c(7,5,6,8,1,3,4,2)
mapping_val = mapping_val[indeces]
mapping_cor = round(mapping_cor[indeces],digits=2)
print(mapping_val)
print(mapping_cor)
alpha = alpha[,indeces]
beta = beta[indeces,]
S17 = cosmic_signatures["SBS17a",]+cosmic_signatures["SBS17b",]
S17 = S17 / sum(S17)
print(round(as.numeric(cosine(as.numeric(beta[6,]),as.numeric(S17))),digits=2))
names = c("NMF1 (SBS1 - 0.96)","NMF2 (SBS2 - 0.99)","NMF3 (SBS3 - 0.91)","NMF4 (SBS5 - 0.91)","NMF5 (SBS13 - 0.96)","NMF6 (SBS17a+b - 0.86)","NMF7 (SBS18 - 0.94)","NMF8 (SBS26 - 0.85)")
colnames(alpha) = names
rownames(beta) = names
NMF = list()
NMF[["alpha"]] = alpha
NMF[["beta"]] = beta
NMF[["Silhouette_Consensus"]] = silhouette_consensus
NMF[["Explained_Variance"]] = explained_variance

# nsNMF
alpha = nsNMF$alpha["Rank 8",][[1]]
beta = nsNMF$beta["Rank 8",][[1]]
silhouette_consensus = as.numeric(nsNMF$rank[,"Silhouette_Consensus"])
names(silhouette_consensus) = paste0("Signature ",1:15)
explained_variance = as.numeric(nsNMF$rank[,"Explained_Variance"])
names(explained_variance) = paste0("Signature ",1:15)
mapping_val  = NULL
mapping_cor  = NULL
for(i in 1:nrow(beta)) {
    tmp = NULL
    for(j in 1:nrow(cosmic_signatures)) {
        tmp = c(tmp,as.numeric(cosine(as.numeric(beta[i,]),as.numeric(cosmic_signatures[j,]))))
    }
    mapping_val = c(mapping_val,rownames(cosmic_signatures)[which(tmp==max(tmp))][1])
    mapping_cor = c(mapping_cor,max(tmp))
}
indeces = c(3,7,5,8,4,2,1,6)
mapping_val = mapping_val[indeces]
mapping_cor = round(mapping_cor[indeces],digits=2)
print(mapping_val)
print(mapping_cor)
alpha = alpha[,indeces]
beta = beta[indeces,]
S17 = cosmic_signatures["SBS17a",]+cosmic_signatures["SBS17b",]
S17 = S17 / sum(S17)
print(round(as.numeric(cosine(as.numeric(beta[6,]),as.numeric(S17))),digits=2))
names = c("nsNMF1 (SBS1 - 0.98)","nsNMF2 (SBS2 - 0.87)","nsNMF3 (SBS3 - 0.87)","nsNMF4 (SBS5 - 0.78)","nsNMF5 (SBS13 - 0.76)","nsNMF6 (SBS17a+b - 0.65)","nsNMF7 (SBS18 - 0.96)","nsNMF8 (SBS26 - 0.90)")
colnames(alpha) = names
rownames(beta) = names
nsNMF = list()
nsNMF[["alpha"]] = alpha
nsNMF[["beta"]] = beta
nsNMF[["Silhouette_Consensus"]] = silhouette_consensus
nsNMF[["Explained_Variance"]] = explained_variance

# SignatureAnalyzer
alpha = SignatureAnalyzer$alpha["Rank 8",][[1]]
beta = SignatureAnalyzer$beta["Rank 8",][[1]]
mean_cosine_similarity = as.numeric(SignatureAnalyzer$rank[,"Goodness of fit (mean cosine similarity)"])
names(mean_cosine_similarity) = paste0("Signature ",1:12)
goodness_fit = as.numeric(SignatureAnalyzer$rank[,"Goodness of fit (percent >0.95 cosine similarity)"])
names(goodness_fit) = paste0("Signature ",1:12)
mapping_val  = NULL
mapping_cor  = NULL
for(i in 1:nrow(beta)) {
    tmp = NULL
    for(j in 1:nrow(cosmic_signatures)) {
        tmp = c(tmp,as.numeric(cosine(as.numeric(beta[i,]),as.numeric(cosmic_signatures[j,]))))
    }
    mapping_val = c(mapping_val,rownames(cosmic_signatures)[which(tmp==max(tmp))][1])
    mapping_cor = c(mapping_cor,max(tmp))
}
indeces = c(3,4,5,1,2,8,6,7)
mapping_val = mapping_val[indeces]
mapping_cor = round(mapping_cor[indeces],digits=2)
print(mapping_val)
print(mapping_cor)
alpha = alpha[,indeces]
beta = beta[indeces,]
print(round(as.numeric(cosine(as.numeric(beta[3,]),as.numeric(cosmic_signatures["SBS3",]))),digits=2))
print(round(as.numeric(cosine(as.numeric(beta[3,]),as.numeric(cosmic_signatures["SBS8",]))),digits=2))
S17 = cosmic_signatures["SBS17a",]+cosmic_signatures["SBS17b",]
S17 = S17 / sum(S17)
print(round(as.numeric(cosine(as.numeric(beta[6,]),as.numeric(S17))),digits=2))
names = c("SIA1 (SBS1 - 0.98)","SIA2 (SBS2 - 0.99)","SIA3 (SBS3 - 0.87)","SIA4 (SBS5 - 0.88)","SIA5 (SBS13 - 0.98)","SIA6 (SBS17a+b - 0.86)","SIA7 (SBS18 - 0.97)","SIA8 (SBS26 - 0.97)")
colnames(alpha) = names
rownames(beta) = names
SignatureAnalyzer = list()
SignatureAnalyzer[["alpha"]] = alpha
SignatureAnalyzer[["beta"]] = beta
SignatureAnalyzer[["Mean_Cosine_Similarity"]] = mean_cosine_similarity
SignatureAnalyzer[["Goodness_Fit"]] = goodness_fit

# SigProfiler
alpha = SigProfiler$alpha["Rank 8",][[1]]
beta = SigProfiler$beta["Rank 8",][[1]]
Stability = as.numeric(SigProfiler$rank[,"Stability"])
names(Stability) = paste0("Signature ",1:15)
MeanL2 = as.numeric(SigProfiler$rank[,"MeanL2"])
names(MeanL2) = paste0("Signature ",1:15)
mapping_val  = NULL
mapping_cor  = NULL
for(i in 1:nrow(beta)) {
    tmp = NULL
    for(j in 1:nrow(cosmic_signatures)) {
        tmp = c(tmp,as.numeric(cosine(as.numeric(beta[i,]),as.numeric(cosmic_signatures[j,]))))
    }
    mapping_val = c(mapping_val,rownames(cosmic_signatures)[which(tmp==max(tmp))][1])
    mapping_cor = c(mapping_cor,max(tmp))
}
indeces = c(3,6,2,1,5,8,4,7)
mapping_val = mapping_val[indeces]
mapping_cor = round(mapping_cor[indeces],digits=2)
print(mapping_val)
print(mapping_cor)
alpha = alpha[,indeces]
beta = beta[indeces,]
S17 = cosmic_signatures["SBS17a",]+cosmic_signatures["SBS17b",]
S17 = S17 / sum(S17)
print(round(as.numeric(cosine(as.numeric(beta[6,]),as.numeric(S17))),digits=2))
names = c("SIP1 (SBS1 - 0.97)","SIP2 (SBS2 - 0.99)","SIP3 (SBS3 - 0.91)","SIP4 (SBS5 - 0.91)","SIP5 (SBS13 - 0.97)","SIP6 (SBS17a+b - 0.86)","SIP7 (SBS18 - 0.95)","SIP8 (SBS26 - 0.85)")
colnames(alpha) = names
rownames(beta) = names
SigProfiler = list()
SigProfiler[["alpha"]] = alpha
SigProfiler[["beta"]] = beta
SigProfiler[["Stability"]] = Stability
SigProfiler[["Mean_L2"]] = MeanL2

# SparseSignatures
alpha = SparseSignatures$alpha["Rank 8",][[1]]
beta = SparseSignatures$beta["Rank 8",][[1]]
median_cv_error = as.numeric(SparseSignatures$rank[,"Median cross-validation error"])
names(median_cv_error) = paste0("Signature ",1:15)
goodness_fit = as.numeric(SparseSignatures$rank[,"Goodness of fit (percent >0.95 cosine similarity)"])
names(goodness_fit) = paste0("Signature ",1:15)
mapping_val  = NULL
mapping_cor  = NULL
for(i in 1:nrow(beta)) {
    tmp = NULL
    for(j in 1:nrow(cosmic_signatures)) {
        tmp = c(tmp,as.numeric(cosine(as.numeric(beta[i,]),as.numeric(cosmic_signatures[j,]))))
    }
    mapping_val = c(mapping_val,rownames(cosmic_signatures)[which(tmp==max(tmp))][1])
    mapping_cor = c(mapping_cor,max(tmp))
}
indeces = c(3,7,2,1,5,8,4,6)
mapping_val = mapping_val[indeces]
mapping_cor = round(mapping_cor[indeces],digits=2)
print(mapping_val)
print(mapping_cor)
alpha = alpha[,indeces]
beta = beta[indeces,]
S17 = cosmic_signatures["SBS17a",]+cosmic_signatures["SBS17b",]
S17 = S17 / sum(S17)
print(round(as.numeric(cosine(as.numeric(beta[6,]),as.numeric(S17))),digits=2))
names = c("SPS1 (SBS1 - 0.94)","SPS2 (SBS2 - 1.00)","SPS3 (SBS3 - 0.89)","SPS4 (SBS5 - 0.98)","SPS5 (SBS13 - 0.97)","SPS6 (SBS17a+b - 0.84)","SPS7 (SBS18 - 0.96)","SPS8 (SBS26 - 0.86)")
colnames(alpha) = names
rownames(beta) = names
SparseSignatures = list()
SparseSignatures[["alpha"]] = alpha
SparseSignatures[["beta"]] = beta
SparseSignatures[["Median_CV_Error"]] = median_cv_error
SparseSignatures[["Goodness_Fit"]] = goodness_fit

# save results
save(NMF,file="RData/NMF.RData")
save(nsNMF,file="RData/nsNMF.RData")
save(SignatureAnalyzer,file="RData/SignatureAnalyzer.RData")
save(SigProfiler,file="RData/SigProfiler.RData")
save(SparseSignatures,file="RData/SparseSignatures.RData")

# make plots of discovered signatures
signatures.plot.v2(NMF$beta,firstBackground=FALSE,xlabels=FALSE)
signatures.plot.v2(nsNMF$beta,firstBackground=FALSE,xlabels=FALSE)
signatures.plot.v2(SignatureAnalyzer$beta,firstBackground=FALSE,xlabels=FALSE)
signatures.plot.v2(SigProfiler$beta,firstBackground=FALSE,xlabels=FALSE)
signatures.plot.v2(SparseSignatures$beta,firstBackground=FALSE,xlabels=FALSE)
