# load required libraries
library("ggplot2")

# load data
load(file="final_results/SparseSignatures.RData")
normalized_alpha = SparseSignatures$alpha/rowSums(SparseSignatures$alpha)
normalized_alpha = normalized_alpha[c("Golden_A_Truncal","Golden_A13_Primary","Golden_A15_Relapse","Golden_A18_Relapse","Golden_B_Truncal","Golden_B13_Primary","Golden_B15_Relapse","MD_Truncal","MD_Primary","MD_Relapse","ML_Truncal","ML_Primary","ML_Relapse","C_Truncal","C_Primary","C_Relapse","E_Truncal","E_Primary","E_Relapse"),]
rownames(normalized_alpha) = c("01A_T","02A13_P","03A15_R","04A18_R","05B_T","06B13_P","07B15_R","08MD_T","09MD_P","10MD_R","11ML_T","12ML_P","13ML_R","14C_T","15C_P","16C_R","17E_T","18E_P","19E_R")
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
ggplot(data=data.frame(EXPOSURE=EXPOSURE,SIGNATURE=SIGNATURE,SAMPLE=SAMPLE),aes(x=SAMPLE,y=EXPOSURE,fill=SIGNATURE)) + geom_bar(stat="identity") + scale_fill_brewer(palette="Spectral") + ggtitle("GOLDEN dataset - Signatures Assignments")
