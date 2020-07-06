# load required libraries
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

# load data
load("RData/trinucleotides_counts.RData")
load(file="final_results/SparseSignatures.RData")
predicted_counts = SparseSignatures$alpha%*%SparseSignatures$beta

# data structures to make plots
data1_obs = trinucleotides_counts[c("Golden_A_Truncal","Golden_A13_Primary","Golden_A15_Relapse","Golden_A18_Relapse"),]
rownames(data1_obs) = c("Golden_A_Truncal (0.82)","Golden_A13_Primary (0.99)","Golden_A15_Relapse (0.88)","Golden_A18_Relapse (0.95)")
data2_obs = trinucleotides_counts[c("Golden_B_Truncal","Golden_B13_Primary","Golden_B15_Relapse"),]
rownames(data2_obs) = c("Golden_B_Truncal (0.97)","Golden_B13_Primary (0.89)","Golden_B15_Relapse (0.96)")
data3_obs = trinucleotides_counts[c("MD_Truncal","MD_Primary","MD_Relapse"),]
rownames(data3_obs) = c("MD_Truncal (0.93)","MD_Primary (0.74)","MD_Relapse (0.59)")
data4_obs = trinucleotides_counts[c("ML_Truncal","ML_Primary","ML_Relapse"),]
rownames(data4_obs) = c("ML_Truncal (0.91)","ML_Primary (0.88)","ML_Relapse (0.92)")
data5_obs = trinucleotides_counts[c("C_Truncal","C_Primary","C_Relapse"),]
rownames(data5_obs) = c("C_Truncal (0.96)","C_Primary (0.92)","C_Relapse (0.99)")
data6_obs = trinucleotides_counts[c("E_Truncal","E_Primary","E_Relapse"),]
rownames(data6_obs) = c("E_Truncal (0.92)","E_Primary (0.94)","E_Relapse (0.93)")
data1_pre = predicted_counts[c("Golden_A_Truncal","Golden_A13_Primary","Golden_A15_Relapse","Golden_A18_Relapse"),]
rownames(data1_pre) = c("Golden_A_Truncal (0.82)","Golden_A13_Primary (0.99)","Golden_A15_Relapse (0.88)","Golden_A18_Relapse (0.95)")
data2_pre = predicted_counts[c("Golden_B_Truncal","Golden_B13_Primary","Golden_B15_Relapse"),]
rownames(data2_pre) = c("Golden_B_Truncal (0.97)","Golden_B13_Primary (0.89)","Golden_B15_Relapse (0.96)")
data3_pre = predicted_counts[c("MD_Truncal","MD_Primary","MD_Relapse"),]
rownames(data3_pre) = c("MD_Truncal (0.93)","MD_Primary (0.74)","MD_Relapse (0.59)")
data4_pre = predicted_counts[c("ML_Truncal","ML_Primary","ML_Relapse"),]
rownames(data4_pre) = c("ML_Truncal (0.91)","ML_Primary (0.88)","ML_Relapse (0.92)")
data5_pre = predicted_counts[c("C_Truncal","C_Primary","C_Relapse"),]
rownames(data5_pre) = c("C_Truncal (0.96)","C_Primary (0.92)","C_Relapse (0.99)")
data6_pre = predicted_counts[c("E_Truncal","E_Primary","E_Relapse"),]
rownames(data6_pre) = c("E_Truncal (0.92)","E_Primary (0.94)","E_Relapse (0.93)")

# make plots
signatures.plot.v2(data1_obs,useColNames=TRUE,firstBackground=FALSE,xlabels=FALSE)
signatures.plot.v2(data1_pre,useColNames=TRUE,firstBackground=FALSE,xlabels=FALSE)
signatures.plot.v2(data2_obs,useColNames=TRUE,firstBackground=FALSE,xlabels=FALSE)
signatures.plot.v2(data2_pre,useColNames=TRUE,firstBackground=FALSE,xlabels=FALSE)
signatures.plot.v2(data3_obs,useColNames=TRUE,firstBackground=FALSE,xlabels=FALSE)
signatures.plot.v2(data3_pre,useColNames=TRUE,firstBackground=FALSE,xlabels=FALSE)
signatures.plot.v2(data4_obs,useColNames=TRUE,firstBackground=FALSE,xlabels=FALSE)
signatures.plot.v2(data4_pre,useColNames=TRUE,firstBackground=FALSE,xlabels=FALSE)
signatures.plot.v2(data5_obs,useColNames=TRUE,firstBackground=FALSE,xlabels=FALSE)
signatures.plot.v2(data5_pre,useColNames=TRUE,firstBackground=FALSE,xlabels=FALSE)
signatures.plot.v2(data6_obs,useColNames=TRUE,firstBackground=FALSE,xlabels=FALSE)
signatures.plot.v2(data6_pre,useColNames=TRUE,firstBackground=FALSE,xlabels=FALSE)
