

# filter SampleTable based on column values:
SampleTableAnalysis <- SampleTable[ !(SampleTable$CHX == TRUE | SampleTable$treatment == "ATR"),]
SampleTableAnalysisAcute <- SampleTableAnalysis[ (SampleTableAnalysis$time != 336), ]
SampleTableAnalysisChronic <- SampleTableAnalysis[ (SampleTableAnalysis$time == 336), ]
SampleTableAnalysisAcute <- SampleTableAnalysisAcute[complete.cases(SampleTableAnalysisAcute),]
SampleTableAnalysisChronic <- SampleTableAnalysisChronic[complete.cases(SampleTableAnalysisChronic),]

# convert time to a factor:
SampleTableAnalysis$timeAsFactor <- factor(SampleTableAnalysis$time)
SampleTableAnalysisAcute$timeAsFactor <- factor(SampleTableAnalysisAcute$time)
SampleTableAnalysisChronic$timeAsFactor <- factor(SampleTableAnalysisChronic$time)

# R prefers that the control is the first level among the factors:
library("magrittr")
SampleTableAnalysis$treatment %<>% relevel("EtOH")
SampleTableAnalysisAcute$treatment %<>% relevel("EtOH")
SampleTableAnalysisChronic$treatment %<>% relevel("EtOH")

files <- file.path(  "/home/poxley/Documents/bioinformatics/Projects/Marcello RNASeq/Sample_data/quants" , samples, "quant.sf")
names(files) <- samples
all(file.exists(files))
files_acute <- files[c(as.numeric(rownames(SampleTableAnalysisAcute)))]
files_chronic <- files[c(as.numeric(rownames(SampleTableAnalysisChronic)))]

tx2gene_names <- read.table(tx2gene_names_file, header=TRUE, sep=",")
txi_named <- tximport(files, type = "salmon", tx2gene = tx2gene_names)
txi_acute <- tximport(files_acute, type = "salmon", tx2gene = tx2gene_names)
txi_chronic <- tximport(files_chronic, type = "salmon", tx2gene = tx2gene_names)

rownames(SampleTableAnalysisAcute) <- colnames(txi_acute$counts)
rownames(SampleTableAnalysisChronic) <- colnames(txi_chronic$counts)

#dds_taf <- DESeqDataSetFromTximport(txi_named, SampleTableAnalysis, ~ treatment + timeAsFactor)
#dds_taf_intx <- DESeqDataSetFromTximport(txi_named, SampleTableAnalysis, ~ treatment + timeAsFactor + treatment:timeAsFactor)
dds_taf_intx_acute <- DESeqDataSetFromTximport(txi_acute, SampleTableAnalysisAcute, ~ treatment + timeAsFactor + treatment:timeAsFactor)
dds_taf_chronic <- DESeqDataSetFromTximport(txi_chronic, SampleTableAnalysisChronic, ~ treatment)

#dds_lrt_taf <- DESeq(dds_taf, test="LRT", reduced = ~1)
#dds_lrt_taf_intx <- DESeq(dds_taf_intx, test="LRT", reduced=~ treatment + timeAsFactor)
dds_lrt_taf_intx_acute <- DESeq(dds_taf_intx_acute, test="LRT", reduced=~ treatment + timeAsFactor)
dds_lrt_taf_chronic <- DESeq(dds_taf_chronic, test="LRT", reduced=~ 1 )

#res_lrt_intx = results(dds_lrt_taf_intx)
res_lrt_intx_acute = results(dds_lrt_taf_intx_acute)
res_lrt_chronic = results(dds_lrt_taf_chronic)

#write.table(as.data.frame(res_lrt_intx), quote=FALSE, file="marcelo_results.lrt_intx.out.txt")
write.table(as.data.frame(res_lrt_intx_acute), quote=FALSE, file="marcelo_results.lrt_intx_acute.out.txt")
write.table(as.data.frame(res_lrt_chronic), quote=FALSE, file="marcelo_results.lrt_intx_chronic.out.txt")

#res_a2e_etoh <- results(dds_lrt_taf, alpha=0.05, contrast = c("treatment", "EtOH", "A2E"))
#res_atrd_etoh <- results(dds_lrt_taf, alpha=0.05, contrast = c("treatment", "EtOH", "ATRD"))
#res_atrd_a2e <- results(dds_lrt_taf, alpha=0.05, contrast = c("treatment", "A2E", "ATRD"))

#res_a2e_etoh_intx <- results(dds_lrt_taf_intx, alpha=0.05, contrast = c("treatment", "EtOH", "A2E"))
#res_atrd_etoh_intx <- results(dds_lrt_taf_intx, alpha=0.05, contrast = c("treatment", "EtOH", "ATRD"))
#res_atrd_a2e_intx <- results(dds_lrt_taf_intx, alpha=0.05, contrast = c("treatment", "A2E", "ATRD"))

vsd_acute <- vst(dds_lrt_taf_intx_acute, blind = TRUE)
vsd_chronic <- vst(dds_lrt_taf_chronic, blind = TRUE)
vstMat_acute <- assay(vsd_acute)
vstMat_chronic <- assay(vsd_chronic)
write.table(as.data.frame(vstMat_acute), file = "marcelo.vst_acute.txt", quote = FALSE)
write.table(as.data.frame(vstMat_chronic), file = "marcelo.vst_chronic.txt", quote = FALSE)

# plot the most significant gene:
library(ggplot2)
fiss <- plotCounts(dds_lrt_taf_intx_acute, which.min(res_lrt_intx_acute$padj), 
                    intgroup = c("treatment","timeAsFactor"), returnData = TRUE)
ggplot(fiss,
        aes(x = as.numeric(timeAsFactor), y = count, color = treatment, group = treatment))  
     geom_point() + geom_smooth(se = FALSE, method = "loess") + scale_y_log10()

# this step clusters the genes by looking at their log2foldchange values for all comparisons
betas <- coef(dds_lrt_taf_intx_acute)
colnames(betas)

topGenes <- head(order(res_lrt_intx_acute$padj),50)
mat <- betas[topGenes, -c(1,4,5,6)]
thr <- 3 
mat[mat < -thr] <- -thr
mat[mat > thr] <- thr
pheatmap(mat, breaks=seq(from=-thr, to=thr, length=101),
          cluster_col=FALSE)