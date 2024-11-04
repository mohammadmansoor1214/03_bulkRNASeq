library(tidyverse)
library(DESeq2)


balbc_rep1 <- read.delim(file = "~/03_bulkRNASeq/Project2_BulkSeqAnalysis/balbc_rep1_selected_strandedness_counts.txt", header = TRUE, sep = "\t", row.names = 1)
balbc_rep2 <- read.delim(file = "~/03_bulkRNASeq/Project2_BulkSeqAnalysis/balbc_rep2_selected_strandedness_counts.txt", header = TRUE, sep = "\t", row.names = 1)
balbc_rep3 <- read.delim(file = "~/03_bulkRNASeq/Project2_BulkSeqAnalysis/balbc_rep3_selected_strandedness_counts.txt", header = TRUE, sep = "\t", row.names = 1)
SKG_rep1 <- read.delim(file = "~/03_bulkRNASeq/Project2_BulkSeqAnalysis/SKG_rep1_selected_strandedness_counts.txt", header = TRUE, sep = "\t", row.names = 1)
SKG_rep2 <- read.delim(file = "~/03_bulkRNASeq/Project2_BulkSeqAnalysis/SKG_rep2_selected_strandedness_counts.txt", header = TRUE, sep = "\t", row.names = 1)
SKG_rep3 <- read.delim(file = "~/03_bulkRNASeq/Project2_BulkSeqAnalysis/SKG_rep3_selected_strandedness_counts.txt", header = TRUE, sep = "\t", row.names = 1)

#reading in the the rna seq qc metrics and metadata
myMetadataFile <- read.csv(file = "~/03_bulkRNASeq/Project2_BulkSeqAnalysis/RNASeq_qc_metrics.csv", header = TRUE, sep = ",")


#merging the samples together to create a geneCount file
balbc_rep1 <- balbc_rep1 %>%
  rownames_to_column(var = "genes")
balbc_rep2 <- balbc_rep2 %>%
  rownames_to_column(var = "genes")
balbc_rep3 <- balbc_rep3 %>%
  rownames_to_column(var = "genes")

SKG_rep1 <- SKG_rep1 %>%
  rownames_to_column(var = "genes")
SKG_rep2 <- SKG_rep2 %>%
  rownames_to_column(var = "genes")
SKG_rep3 <- SKG_rep3 %>%
  rownames_to_column(var = "genes")

#merging all data frames together

alldfs <- list(balbc_rep1, balbc_rep2, balbc_rep3, SKG_rep1, SKG_rep2, SKG_rep3)
geneCount <- merge(balbc_rep1,balbc_rep2)
geneCount <- merge(geneCount,balbc_rep3) 
geneCount <- merge(geneCount, SKG_rep1)
geneCount <- merge(geneCount, SKG_rep2)
geneCount <- merge(geneCount, SKG_rep3)

#making the fist column the row names in the geneCount dataframe

rownames(geneCount) <- geneCount$genes
geneCount <- geneCount[, -1]

#changing the rownames to be samples names in the metadata file
rownames(myMetadataFile) <- myMetadataFile$short_sample_ID
myMetadataFile <- myMetadataFile[, -1]

rownames(myMetadataFile) <- c("BALBc_rep1", "BALBc_rep2", "BALBc_rep3", "SKG_rep1", "SKG_rep2", "SKG_rep3")

#STEP5: removing genes from the gene counts matrix that are have less than or 15 reads

keep <- rowSums(geneCount) >=15
filteredGeneCounts <- geneCount[keep, ]

#STEP5c data normalization via DESeq2


deObj <- DESeqDataSetFromMatrix(countData = filteredGeneCounts,
                                colData = myMetadataFile,
                                design = ~ biological_group)
deObj <- DESeq(deObj)

normCountsDf <- counts(deObj, normalized = TRUE)
rawCountsDF <- counts(deObj, normalized = FALSE)


#STEP5d PCA and covariate identification

pcs <- prcomp(t(normCountsDf), scale. = TRUE)
pc_variance_summary <- summary(pcs)
myMetadataFile <- cbind(pcs$x, myMetadataFile)

ggplot(myMetadataFile, aes(x = PC1, y = PC2, color = PCT_RIBOSOMAL_BASES)) +
  geom_point(size = 4) +
  xlab(paste("PC1 - ", pc_variance_summary$importance["Proportion of Variance", "PC1"]*100, "%", sep = " ")) +
  ylab(paste("PC2 - ", pc_variance_summary$importance["Proportion of Variance", "PC2"]*100, "%", sep = " "))


#testing how PC2 vs PC3 cluster together 
ggplot(myMetadataFile, aes(x = PC2, y = PC3, color = MEDIAN_5PRIME_BIAS)) +
  geom_point(size = 4) +
  xlab(paste("PC2 - ", pc_variance_summary$importance["Proportion of Variance", "PC2"]*100, "%", sep = " ")) +
  ylab(paste("PC3 - ", pc_variance_summary$importance["Proportion of Variance", "PC3"]*100, "%", sep = " "))


#STEP6: differerential gene analysis

myResults <- results(object = deObj, contrast = c("biological_group", "SKG", "BALBc"))
SKG_vs_BLABc <- as.data.frame(myResults)


#save all R data objects so i can go back to them later without having to reanalyze everything 

save(list = c("myMetadataFile", "deObj", "filteredGeneCounts", "geneCount", "pcs"),
     file = "Project2_BALBc_vs_SKG_bulk_RNAseq_analysys.Rdat")


#generating the volcano plots for which we need log2FC which we already have but also -log10FDR which we will have to calculate

SKG_vs_BLABc$negLog10FDR <- (-1)*log10(SKG_vs_BLABc$padj)
SKG_vs_BLABc$genes <- rownames(SKG_vs_BLABc)

#finally creating a volvano plot

ggplot(SKG_vs_BLABc, aes(x = log2FoldChange, y = negLog10FDR)) +
  geom_point() +
  geom_vline(xintercept = c(-1, 1), col = "red") +
  geom_hline (yintercept = (-1)*log10(0.10), col = "red") +
  geom_text_repel(data=subset(SKG_vs_BLABc, padj < 0.10), aes(label = genes))
  