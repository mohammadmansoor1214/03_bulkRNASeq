
mymetadata <- read.csv("~/03_bulkRNASeq/bulk_RNAseq_data_qc_metrics.csv", header = TRUE, sep = ",")
load("~/03_bulkRNASeq/myQuantificationResults.Rdat")
geneCounts <- myCounts$counts
#checking if row names are gene names
row.names(geneCounts)
colnames(geneCounts)
#check if rownames are samples in the metadata
rownames(mymetadata)
rownames(mymetadata) <- mymetadata$short_sample_ID
