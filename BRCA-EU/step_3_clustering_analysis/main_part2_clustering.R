# load required libraries
library("FactoMineR")
library("lsa")
library("factoextra")
library("tidyverse")
library("magrittr")
library("ggiraphExtra")
library("survival")
library("rms")
library("survminer")
library("ggplot2")

# set the seed
set.seed(12345)

# load data
load(file="data/SparseSignatures.RData")
load(file="processed_data/clinical_information.RData")
raw_alpha = SparseSignatures$alpha

# compute distance among samples based on signatures distributions
signatures_distance = cosine(t(raw_alpha))

# elbow method to assess number of clusters
fviz_nbclust(signatures_distance,kmeans,method="wss") + labs(subtitle="Elbow method")

# PCA analysis
res.pca = PCA(signatures_distance,graph=FALSE)

# visualize eigenvalues/variances
fviz_screeplot(res.pca,addlabels=TRUE,ylim=c(0,100))

# perform and visualize clustering
km = kmeans(res.pca$var$coord[,1:3],centers=3,nstart=1000)
clusters = t(t(km$cluster))
colnames(clusters) = "Clusters"
clusters[which(clusters[,"Clusters"]==1),"Clusters"] = "C1"
clusters[which(clusters[,"Clusters"]==2),"Clusters"] = "C2"
clusters[which(clusters[,"Clusters"]==3),"Clusters"] = "C3"
save(clusters,file="results/clusters.RData")
fviz_cluster(km,data=signatures_distance,ellipse.type ="convex") + theme_minimal()

# visualize features per cluster
data_df = as.data.frame(raw_alpha) %>% rownames_to_column()
cluster_pos = as.data.frame(km $cluster) %>% rownames_to_column()
colnames(cluster_pos) = c("rowname","cluster")
data_final = inner_join(cluster_pos,data_df)
colnames(data_final)[3:ncol(data_final)] = paste0("S",1:(ncol(data_final)-2))
ggRadar(data_final[-1],aes(group=cluster),rescale=FALSE,legend.position="none",size=1,interactive=FALSE,use.label=TRUE) + facet_wrap(~cluster) + scale_y_discrete(breaks=NULL) + theme(axis.text.x=element_text(size=10)) + scale_fill_manual(values=rep("#1c6193",nrow(data_final))) + scale_color_manual(values=rep("#1c6193",nrow(data_final))) + ggtitle("Signatures per Cluster")

# visualize normalized alpha per cluster
normalized_alpha = raw_alpha / rowSums(raw_alpha)
colnames(normalized_alpha) = c("SPS1","SPS2","SPS3","SPS4","SPS5","SPS6","SPS7","SPS8")
CLUSTER = NULL
SIGNATURE = NULL
EXPOSURE = NULL
for(i in 1:nrow(normalized_alpha)) {
    for(j in 1:ncol(normalized_alpha)) {
        CLUSTER = c(CLUSTER,clusters[rownames(normalized_alpha)[i],"Clusters"])
        SIGNATURE = c(SIGNATURE,colnames(normalized_alpha)[j])
        EXPOSURE = c(EXPOSURE, normalized_alpha[i,j])
    }
}
data = data.frame(CLUSTER=CLUSTER,SIGNATURE=SIGNATURE,EXPOSURE=EXPOSURE)
ggplot(data,aes(x=SIGNATURE,y=EXPOSURE,fill=CLUSTER)) + geom_boxplot()

# extract clinical features
AGE = NULL
GENDER = NULL
SIZE = NULL
CLUSTER = NULL
PAM50 = NULL
for(i in sort(unique(clusters[,"Clusters"]))) {
    curr_samples = rownames(raw_alpha[which(clusters[,"Clusters"]==i),])
    SIZE = c(SIZE,rep(length(curr_samples),length(curr_samples)))
    CLUSTER = c(CLUSTER,rep(i,length(curr_samples)))
    for(j in curr_samples) {
        if(length(which(clinical_information$PATIENT_ID==j))==1) {
            AGE = c(AGE,as.numeric(clinical_information$AGE[which(clinical_information$PATIENT_ID==j)]))
            GENDER = c(GENDER,as.character(clinical_information$GENDER[which(clinical_information$PATIENT_ID==j)]))
            PAM50 = c(PAM50,as.character(clinical_information$PAM50[which(clinical_information$PATIENT_ID==j)]))
        }
        else {
            AGE = c(AGE,NA)
            GENDER = c(GENDER,NA)
            PAM50 = c(PAM50,NA)
        }
    }
}

# plots for AGE and GENDER
data = data.frame(AGE=AGE,GENDER=GENDER,PAM50=PAM50,SIZE=SIZE,CLUSTER=CLUSTER)
ggplot(data,aes(x=CLUSTER,y=AGE,fill=CLUSTER)) + geom_boxplot()
ggplot(data,aes(x=CLUSTER,y=SIZE,fill=GENDER)) + geom_bar(stat="identity",position="fill")
ggplot(data,aes(x=CLUSTER,y=SIZE,fill=PAM50)) + geom_bar(stat="identity",position="fill")
