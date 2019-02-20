# initialise environment
setwd("..")
library(tidyverse)
library(biomaRt)

######################################
### Prepare gene conversion tables ###
######################################

# translate the transcripts into genes
ensembl <- useMart("ensembl",dataset="mmusculus_gene_ensembl")
attributes <- listAttributes(ensembl)

# you can search for appropriate identifiers using the following:
hits <- grep(pattern="gene.name", x=attributes$description, ignore.case=T)
attributes[hits,]

# get list of all mRNA ids to convert:
eg_quant <- read_tsv("salmon-output/01-EC---3-7-2018_S25_R1_001.quant/quant.sf")
head(eg_quant)

# create data frame with name conversions
# transcript names must be stripped of their version ids for this db
tx2gene <- getBM(attributes=c('refseq_mrna', 'entrezgene', 'uniprot_gn', 'external_gene_name'),
                 filters = 'refseq_mrna',
                 values = sub("\\.\\d\\d?", "", eg_quant$Name),
                 mart = ensembl)
tx2gene <- as_tibble(tx2gene)
write.table(as.data.frame(tx2gene), file = "scripts/transcript-name-conversion.csv",
            quote = FALSE,
            sep = ',')
tx2gene

# create dataframe for tximport name conversion:
tx2gene_mart <- tx2gene %>% 
  dplyr::select(refseq_mrna, entrezgene) %>% 
  unique()

# old version - kept for archive 
#tx2gene_mart <- getBM(attributes=c('refseq_mrna', 'entrezgene'),
#                     filters = 'refseq_mrna',
#                     values = sub("\\.\\d\\d?", "", eg_quant$Name),
#                     mart = ensembl)

dim(tx2gene_mart)

# remove the rows with no corresponding geneid
tx2gene_mart_clean <- tx2gene_mart[complete.cases(tx2gene_mart),]
dim(tx2gene_mart_clean)


###################################
### Prepare Salmon quant files  ###
###################################

# set salmon quantification file locations and details
data_dir <- "salmon-output/"
list.files(data_dir)
samples <- read_tsv(file.path("docs", "sample_info.txt"), col_names = TRUE )
samples$sample_name <- paste(samples$condition, samples$cell, samples$replicate, sep="_")

# set condition and cell type as factors:
samples$condition <- factor(samples$condition, levels = c("WT","PBS", "Knull", "KO1", "KO2"))
samples$cell <- factor(samples$cell, levels = c("EC", "MG"))


samples_endothelial <-  samples %>%
  dplyr::filter(cell == "EC") %>% 
  dplyr::filter(condition != "Knull") %>% 
  dplyr::filter(sample_name != "WT_EC_C") %>% 

samples_microglial <-  samples %>%
  dplyr::filter(cell == "MG") %>% 
  dplyr::filter(condition != "Knull") %>% 
  dplyr::filter(sample_name != "PBS_MG_C")



files_ec <- paste(data_dir,samples_endothelial$filename,"/quant.sf.trim", sep="")
names(files_ec) <- paste(samples_endothelial$condition, 
                         samples_endothelial$cell, 
                         samples_endothelial$replicate, 
                         sep="_")
all(file.exists(files_ec))

files_mg <- paste(data_dir,samples_microglial$filename,"/quant.sf.trim", sep="")
names(files_mg) <- paste(samples_microglial$condition, 
                         samples_microglial$cell, 
                         samples_microglial$replicate, 
                         sep="_")
all(file.exists(files_mg))

# import into tximport data frame
library(tximport)
txi <- tximport(files_trimmed, type = "salmon", tx2gene = tx2gene_mart_clean)
names(txi)
head(txi$counts)

txi_ec <- tximport(files_ec, type = "salmon", tx2gene = tx2gene_mart_clean)
names(txi_ec)
head(txi_ec$counts)

txi_mg <- tximport(files_mg, type = "salmon", tx2gene = tx2gene_mart_clean)
names(txi_mg)
head(txi_mg$counts)

###################################
###         run DESeq2          ###
###################################

# for details, see http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html


library(DESeq2)


# create DESeq2 Data Set dds_ec
dds_ec <- DESeqDataSetFromTximport(txi_ec, samples_endothelial, ~ replicate + condition)
dds_mg <- DESeqDataSetFromTximport(txi_mg, samples_microglial, ~ replicate + condition)

dds_ec <- DESeq(dds_ec, test="LRT", reduced = ~ replicate)
dds_mg <- DESeq(dds_mg, test="LRT", reduced = ~ replicate)

res_ec <- results(dds_ec, alpha = 0.05)
res_mg <- results(dds_mg, alpha = 0.05)

summary(res_ec)
summary(res_mg)

write.csv(res_ec, file = "DESeq2_R-project/res_ec-results.csv")
write.csv(res_mg, file = "DESeq2_R-project/res_mg-results.csv")

# extract the variance stabilized transformed counts:
#vsd_ec <- vst(dds_ec, blind = TRUE)
vsd_ec <- varianceStabilizingTransformation(dds_ec, blind = FALSE)
vsd_ec_mat <- assay(vsd_ec)
write.table(as.data.frame(vsd_ec_mat), file = "DESeq2_R-project/counts-ec-vst.csv")

#vsd_mg <- vst(dds_mg, blind = FALSE)
vsd_mg <- varianceStabilizingTransformation(dds_mg, blind = FALSE)
vsd_mg_mat <- assay(vsd_mg)
write.table(as.data.frame(vsd_mg_mat), file = "DESeq2_R-project/counts-mg-vst.csv")

###########################################################################
### Wald tests for post-hoc analysis of DEGs identified by LRT analysis ###
###########################################################################
dds_ec_wald <- DESeqDataSetFromTximport(txi_ec, samples_endothelial, ~ replicate + condition)
dds_mg_wald <- DESeqDataSetFromTximport(txi_mg, samples_microglial, ~ replicate + condition)

dds_ec_wald <- DESeq(dds_ec_wald)
dds_mg_wald <- DESeq(dds_mg_wald)

res_ec_wald_PW <- results(dds_ec_wald, contrast=c("condition","PBS","WT"), alpha=0.05)
res_ec_wald_PK1 <- results(dds_ec_wald, contrast=c("condition","PBS","KO1"), alpha=0.05)
res_ec_wald_PK2 <- results(dds_ec_wald, contrast=c("condition","PBS","KO2"), alpha=0.05)
res_ec_wald_K1K2 <- results(dds_ec_wald, contrast=c("condition","KO1","KO2"), alpha=0.05)
res_ec_wald_WK1 <- results(dds_ec_wald, contrast=c("condition","WT","KO1"), alpha=0.05)
res_ec_wald_WK2 <- results(dds_ec_wald, contrast=c("condition","WT","KO2"), alpha=0.05)

write.csv(res_ec_wald_PW, file = "DESeq2_R-project/res_ec_wald_PW-results.csv")
write.csv(res_ec_wald_PK1, file = "DESeq2_R-project/res_ec_wald_PK1-results.csv")
write.csv(res_ec_wald_PK2, file = "DESeq2_R-project/res_ec_wald_PK2-results.csv")
write.csv(res_ec_wald_K1K2, file = "DESeq2_R-project/res_ec_wald_K1K2-results.csv")
write.csv(res_ec_wald_WK1, file = "DESeq2_R-project/res_ec_wald_WK1-results.csv")
write.csv(res_ec_wald_WK2, file = "DESeq2_R-project/res_ec_wald_WK2-results.csv")


res_mg_wald_PW <- results(dds_mg_wald, contrast=c("condition","PBS","WT"), alpha=0.05)
res_mg_wald_PK1 <- results(dds_mg_wald, contrast=c("condition","PBS","KO1"), alpha=0.05)
res_mg_wald_PK2 <- results(dds_mg_wald, contrast=c("condition","PBS","KO2"), alpha=0.05)
res_mg_wald_K1K2 <- results(dds_mg_wald, contrast=c("condition","KO1","KO2"), alpha=0.05)
res_mg_wald_WK1 <- results(dds_mg_wald, contrast=c("condition","WT","KO1"), alpha=0.05)
res_mg_wald_WK2 <- results(dds_mg_wald, contrast=c("condition","WT","KO2"), alpha=0.05)


write.csv(res_mg_wald_PW, file = "DESeq2_R-project/res_mg_wald_PW-results.csv")
write.csv(res_mg_wald_PK1, file = "DESeq2_R-project/res_mg_wald_PK1-results.csv")
write.csv(res_mg_wald_PK2, file = "DESeq2_R-project/res_mg_wald_PK2-results.csv")
write.csv(res_mg_wald_K1K2, file = "DESeq2_R-project/res_mg_wald_K1K2-results.csv")
write.csv(res_mg_wald_WK1, file = "DESeq2_R-project/res_mg_wald_WK1-results.csv")
write.csv(res_mg_wald_WK2, file = "DESeq2_R-project/res_mg_wald_WK2-results.csv")

# print the session info for reproducibility

plotMA(res_ec_wald_PK1, coef=4, ylim=c(-2,2))

sessionInfo()
