
## Methods
### RNASeq Quality check and alignment
RNASeq fastq files were analyzed with FastQC[1] v0.11.7 to assess sequence base quality, per-base sequence content, GC content, N content, and the sequence length distribution. Reads were subsequently trimmed using Trimmomatic[2] v0.38, to remove Illumina adapters, leading and trailing bases with score ≤ 3, all bases after the sliding window average ≤ 15, and all edited reads ≤ 36 bp.
Reads were aligned to _Mus musculus_ annotation GRCm38.p6 using Salmon[3] v0.10.0. Default parameters were used for building the index and for alignment.


### Differential expression analysis
DESeq2[4] v1.18.1 was used to assess differential gene expression using the likelihood ratio test, with the model `~ replicate + condition` analyzed against the reduced model `~ replicate`. Variance Stabilized Transformed gene counts were used to identify outliers, using both Principal Component Analysis and sample clustering (using the euclidean distance metric and the complete clustering method). Two samples - wildtype replicate C and PBS replicate C - were removed from further analysis.

To identify significant contrasts between treatments, a _post hoc_ analysis of genes differentially expressed according to LRT analysis was performed. DESeq2 was used to perform the nbinomial Wald test for contrasts between PBS, WT, KO1 and KO2. 

The focal gene set was identified as those genes in which:
* the likelihood ratio test was significant (*p* ≤ 0.05)
* WT and PBS were significantly differentially expressed (*p* ≤ 0.05)
* WT was significantly different from both KO1 and KO2 (*p* ≤ 0.05 in each contrast)
* KO1 and KO2 were concordantly up- or down- regulated with regard to WT

Log2(Fold change) values and _p_ values are reported according to the Wald tests.

Gene names were mapped to entrez gene identifiers using ensembl biomart[5], mouse version GRCm38.p6.


The R Script used for differential expression analysis can be found in supplementary file differential-expression-analysis.R The python script used for subsequent post hoc analysis can be found in the Jupyter notebooks endothelial_expression_analysis.ipynb and microglial_expression_analysis.ipynb.


### Bibliography
1. Andrews, S. FastQC: a quality control tool for high throughput sequence data. (2010). at <http://www.bioinformatics.babraham.ac.uk/projects/fastqc>
2. Bolger, A. M., Lohse, M. & Usadel, B. Trimmomatic: a flexible trimmer for Illumina sequence data. Bioinformatics 30, 2114–2120 (2014).
3. Patro, R., Duggal, G., Love, M. I., Irizarry, R. A. & Kingsford, C. Salmon provides fast and bias-aware quantification of transcript expression. Nat. Methods 14, 417–419 (2017).
4. Love, M. I., Huber, W. & Anders, S. Moderated estimation of fold change and  dispersion for RNA-seq data with DESeq2. Genome Biol. 15, 550 (2014).
5. Zerbino, D. R. et al. Ensembl 2018. Nucleic Acids Res. 46, D754–D761 (2018).


---
## Links to scripts, notebooks and results

### Description of contents of file [DESeq2_R-project](https://github.com/oxpeter/goncalo2019/tree/master/DESeq2_R-project)

|File|Description|
|---|---|
|[01_20180529_DESeq2-analysis.R](https://github.com/oxpeter/goncalo2019/blob/master/DESeq2_R-project/01_20180529_DESeq2-analysis.R)|R-Script for creating gene annotations, and executing DESeq2 initial and _post-hoc_ analyses|
|DESeq2_R-project.Rproj|R project file from analysis|
|counts-allcells-vst.csv|Variance Stabilized Transformed read counts, calculated with both microglial and endothelial samples combined |
|counts-ec-vst.csv|Variance Stabilized Transformed read counts, calculated with endothelial cells only|
|counts-mg-vst.csv|Variance Stabilized Transformed read counts, calculated with microglial cells only|
|[res_ec-results.csv](https://github.com/oxpeter/goncalo2019/blob/master/DESeq2_R-project/res_ec-results.csv)| DESeq2 Likelihood ratio test results for endothelial cells: see [notebooks](https://github.com/oxpeter/goncalo2019/tree/master/notebooks) folder for annotated result file |
|res_ec_wald_K1K2-results.csv| _post-hoc_ Wald test for differences between KO1 and KO2 treatments, in endothelial cells|
|res_ec_wald_PK1-results.csv|_post-hoc_ Wald test for differences between PBS and KO1 treatments, in endothelial cells|
|res_ec_wald_PK2-results.csv|_post-hoc_ Wald test for differences between PBS and KO2 treatments, in endothelial cells|
|res_ec_wald_PW-results.csv|_post-hoc_ Wald test for differences between PBS and WT treatments, in endothelial cells|
|res_ec_wald_WK1-results.csv|_post-hoc_ Wald test for differences between WT and KO1 treatments, in endothelial cells|
|res_ec_wald_WK2-results.csv|_post-hoc_ Wald test for differences between WT and KO2 treatments, in endothelial cells|
|[res_mg-results.csv](https://github.com/oxpeter/goncalo2019/blob/master/DESeq2_R-project/res_mg-results.csv)| DESeq2 Likelihood ratio test results for microglial cells: see [notebooks](https://github.com/oxpeter/goncalo2019/tree/master/notebooks) folder for annotated result file|
|res_mg_wald_K1K2-results.csv|_post-hoc_ Wald test for differences between KO1 and KO2 treatments, in microglial cells|
|res_mg_wald_PK1-results.csv|_post-hoc_ Wald test for differences between PBS and KO1 treatments, in microglial cells|
|res_mg_wald_PK2-results.csv|_post-hoc_ Wald test for differences between PBS and KO2 treatments, in microglial cells|
|res_mg_wald_PW-results.csv|_post-hoc_ Wald test for differences between PBS and WT treatments, in microglial cells|
|res_mg_wald_WK1-results.csv|_post-hoc_ Wald test for differences between WT and KO1 treatments, in microglial cells|
|res_mg_wald_WK2-results.csv|_post-hoc_ Wald test for differences between WT and KO2 treatments, in microglial cells|


[Jupyter notebook for _Post-hoc_ analysis of microglial cells](https://github.com/oxpeter/goncalo2019/blob/master/notebooks/microglial_expression_analysis_20181211.ipynb) (methods and results)

[Jupyter notebook for _Post-hoc_ analysis of endothelial cells](https://github.com/oxpeter/goncalo2019/blob/master/notebooks/endothelial_expression_analysis_20181204.ipynb) (methods and results)