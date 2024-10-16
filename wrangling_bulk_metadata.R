#read in metadata
mymetadata <- read.csv("~/03_bulkRNASeq/bulk_RNAseq_data_qc_metrics.csv", header = TRUE)
library(ggplot2)
ggplot2::ggplot(mymetadata, aes(x = replicate, y = PCT_RIBOSOMAL_BASES*100, fill = replicate)) +
                geom_boxplot()

