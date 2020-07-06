# load the required libraries
library("SparseSignatures")
data("mutation_categories")
library("BSgenome.Hsapiens.NCBI.GRCh38")

# make counts matrix for sample Golden
data = read.table("raw_data/data/Golden/deconstructSigs/Golden_deconstructSigs_input_split_vaf0.05.txt",header=TRUE,sep="\t",check.names=FALSE,stringsAsFactors=FALSE)
colnames(data) = c("sample","chromosome","chromosome_start","mutated_from_allele","mutated_to_allele","vaf")
data = data[which(as.matrix(data[,"mutated_from_allele"])%in%c("A","C","G","T")),]
data = data[which(as.matrix(data[,"mutated_to_allele"])%in%c("A","C","G","T")),]
data = data[,c("sample","chromosome","chromosome_start","mutated_from_allele","mutated_to_allele")]
data = unique(data)
data$sample = paste0("Golden_",data$sample)
data1 = data

# make counts matrix for sample MD
data = read.table("raw_data/data/MD/deconstructSigs/MD_deconstructSigs_input_split_vaf0.05.txt",header=TRUE,sep="\t",check.names=FALSE,stringsAsFactors=FALSE)
colnames(data) = c("sample","chromosome","chromosome_start","mutated_from_allele","mutated_to_allele","vaf")
data = data[which(as.matrix(data[,"mutated_from_allele"])%in%c("A","C","G","T")),]
data = data[which(as.matrix(data[,"mutated_to_allele"])%in%c("A","C","G","T")),]
data = data[,c("sample","chromosome","chromosome_start","mutated_from_allele","mutated_to_allele")]
data = unique(data)
data$sample = paste0("MD_",data$sample)
data2 = data

# make counts matrix for sample ML
data = read.table("raw_data/data/ML/deconstructSigs/ML_deconstructSigs_input_split_vaf0.05.txt",header=TRUE,sep="\t",check.names=FALSE,stringsAsFactors=FALSE)
colnames(data) = c("sample","chromosome","chromosome_start","mutated_from_allele","mutated_to_allele","vaf")
data = data[which(as.matrix(data[,"mutated_from_allele"])%in%c("A","C","G","T")),]
data = data[which(as.matrix(data[,"mutated_to_allele"])%in%c("A","C","G","T")),]
data = data[,c("sample","chromosome","chromosome_start","mutated_from_allele","mutated_to_allele")]
data = unique(data)
data$sample = paste0("ML_",data$sample)
data3 = data

# make counts matrix for sample C
data = read.table("raw_data/data/C/deconstructSigs/C_deconstructSigs_input_split_vaf0.05.txt",header=TRUE,sep="\t",check.names=FALSE,stringsAsFactors=FALSE)
colnames(data) = c("sample","chromosome","chromosome_start","mutated_from_allele","mutated_to_allele","vaf")
data = data[which(as.matrix(data[,"mutated_from_allele"])%in%c("A","C","G","T")),]
data = data[which(as.matrix(data[,"mutated_to_allele"])%in%c("A","C","G","T")),]
data = data[,c("sample","chromosome","chromosome_start","mutated_from_allele","mutated_to_allele")]
data = unique(data)
data$sample = paste0("C_",data$sample)
data4 = data

# make counts matrix for sample E
data = read.table("raw_data/data/E/deconstructSigs/E_deconstructSigs_input_split_vaf0.05.txt",header=TRUE,sep="\t",check.names=FALSE,stringsAsFactors=FALSE)
colnames(data) = c("sample","chromosome","chromosome_start","mutated_from_allele","mutated_to_allele","vaf")
data = data[which(as.matrix(data[,"mutated_from_allele"])%in%c("A","C","G","T")),]
data = data[which(as.matrix(data[,"mutated_to_allele"])%in%c("A","C","G","T")),]
data = data[,c("sample","chromosome","chromosome_start","mutated_from_allele","mutated_to_allele")]
data = unique(data)
data$sample = paste0("E_",data$sample)
data5 = data

# make counts matrix
data = rbind(data1,data2,data3,data4,data5)
data = unique(data)
data$chromosome = gsub("chr","",data$chromosome)
trinucleotides_counts = import.counts.data(input=data,bsg=BSgenome.Hsapiens.NCBI.GRCh38,mutation_categories=mutation_categories)

# save the results
save(trinucleotides_counts,file="processed_data/trinucleotides_counts.RData")
