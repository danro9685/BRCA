# read the data
load("data/SparseSignatures.RData")
clinical_data = read.table(file="data/donor.tsv",header=TRUE,sep="\t",check.names=FALSE,stringsAsFactors=FALSE)
PAM50 = read.table(file="data/Supplementary Table 18.Expression.Subtyping.txt",header=TRUE,sep="\t",check.names=FALSE,stringsAsFactors=FALSE)
samples_data = read.table(file="data/sample.tsv",header=TRUE,sep="\t",check.names=FALSE,stringsAsFactors=FALSE)

# process clinical data
clinical_information = NULL
for(v in rownames(SparseSignatures$alpha)) {
    clinical_information = rbind(clinical_information,clinical_data[which(clinical_data$icgc_donor_id==v),])
}
clinical_information = clinical_information[,c("icgc_donor_id","donor_sex","donor_age_at_diagnosis","donor_age_at_diagnosis")]
rownames(clinical_information) = 1:nrow(clinical_information)
colnames(clinical_information) = c("PATIENT_ID","GENDER","AGE","PAM50")
clinical_information$PATIENT_ID = as.character(clinical_information$PATIENT_ID)
clinical_information$GENDER = toupper(as.character(clinical_information$GENDER))
clinical_information$AGE = as.numeric(clinical_information$AGE)
clinical_information$PAM50 = NA

# match PAM50 data to samples data
for(i in 1:nrow(PAM50)) {
    if(length(grep(PAM50[i,"sample_name"],samples_data$submitted_sample_id))>0) {
        match_position = unique(samples_data$icgc_donor_id[grep(PAM50[i,"sample_name"],samples_data$submitted_sample_id)])
        clinical_information[which(clinical_information$PATIENT_ID==match_position),"PAM50"] = PAM50[i,"PAM50.subtype"]
    }
}
clinical_information$PAM50[which(clinical_information$PAM50=="")] = NA

# save results
save(clinical_information,file="processed_data/clinical_information.RData")
