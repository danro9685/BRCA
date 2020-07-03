# load the required libraries
library("data.table")
library("SparseSignatures")
data("mutation_categories")
library("BSgenome.Hsapiens.1000genomes.hs37d5")

# process raw data
data = fread(input="raw_data/simple_somatic_mutation.open.tsv",sep="\t",header=TRUE)
data = data[,c("icgc_donor_id","project_code","chromosome","chromosome_start","chromosome_end","assembly_version","mutation_type","mutated_from_allele","mutated_to_allele","platform","sequencing_strategy")]
data = data[which(data[,"chromosome"]!="MT"),]
data = data[which(data[,"chromosome_start"]==data[,"chromosome_end"]),]
data = data[which(data[,"assembly_version"]=="GRCh37"),]
data = data[which(data[,"mutation_type"]=="single base substitution"),]
data = data[which(as.matrix(data[,"mutated_from_allele"])%in%c("A","C","G","T")),]
data = data[which(as.matrix(data[,"mutated_to_allele"])%in%c("A","C","G","T")),]
data = data[which(as.matrix(data[,"platform"])%in%c("Illumina HiSeq","Illumina HiSeq X Ten","Illumina GA sequencing")),]
data = data[which(data[,"sequencing_strategy"]=="WGS"),]

# make counts matrix
data = data[,c("icgc_donor_id","chromosome","chromosome_start","mutated_from_allele","mutated_to_allele")]
data = unique(data)
trinucleotides_counts = import.counts.data(input=data,bsg=BSgenome.Hsapiens.1000genomes.hs37d5,mutation_categories=mutation_categories)

# save the results
save(trinucleotides_counts,file="processed_data/trinucleotides_counts.RData")
