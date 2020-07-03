# load required libraries and sources
library("SigProfilerExtractorR")

# save input data
load("final_data/trinucleotides_counts.RData")
trinucleotides_counts = ((trinucleotides_counts/rowSums(trinucleotides_counts))*2500) # normalize trinucleotides counts
data = cbind(t(t(colnames(trinucleotides_counts))),t(trinucleotides_counts))
colnames(data)[1] = "Mutation Types"
write.table(data,file="input_data/data.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=TRUE)

# set the seed
set.seed(44444)

# perform inference
sigprofilerextractor(input_type="table",output="results",inputdata="input_data/data.txt",minsigs=1,maxsigs=15,replicates=1000,mtype=c("96,DINUC,ID"),init="random",exome=FALSE,cpu=100)
