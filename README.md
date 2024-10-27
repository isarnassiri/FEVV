# Functional Enrichment of Genomic Variants and Variations (FEVV)
FEVV is a tool for chromatin state and genomic feature eSNP enrichment. In our approach we fifteen-core chromatin states and 10 genomic features from the biomaRt (Ensembl) and UCSC were used to annotate a list of SNPs in eQTL profiles or genomic intervals.
Based on the type of input we follow two strategy for SNP enrichment:

(1) eQTL SNPs as input: for eSNP enrichment analysis, first, for each gene/transcript, we selected associated eSNPs as foreground (F) and SNPs within the 1Mb window around the TSS as background (B) SNP sets. We used the number of overlaps of foreground (f) and background (b) SNP sets in the genomic feature or chromatin state and calculate the enrichment score (z-score and odds ratio). 

(2) A query SNP as input: first, the 1MB window is defined for the query SNP, and the relevant variants captured from the 1000 genomes data in the European population background. Next, we split the SNPs based upon a LD (R² = 0.8 and MAF = 0.01) to foreground (query SNP and its LD proxies) and background SNPs sets. We used the number of overlaps of foreground (f) and background (b) SNP sets in the genomic feature or chromatin state and calculate the enrichment score (z-score and odds ratio).

<a name="installation"/>

### Installation
1. Install the R [(LINK)](https://cran.r-project.org/)
2. Install the free version of rStudio [(LINK)](https://www.rstudio.com/products/rstudio/download/)
3. Run the following command in rStudio to install scQCEA as an R package:

```{r,eval=FALSE}
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
    
library(devtools)
install_github("isarnassiri/FEVV")
```

#'@export
#'@import GenomicRanges
#'@import data.table
#'@export
#'@name eSNPsEnrichmentAnalysis
#'@title Functional Enrichment of eQTL SNPs
#'@description eQTL SNPs as input: for eSNP enrichment analysis, first, for each gene/transcript, we selected associated eSNPs as foreground (F) and SNPs within the 1Mb window around the TSS as background (B) SNP sets. We used the number of overlaps of foreground (f) and background (b) SNP sets in the genomic feature or chromatin state and calculate the enrichment score (z-score and odds ratio).
#'@author {Isar Nassiri, Benjamin Fairfax}
#'@param eQTL
#'eQTL profile including the following headers: seqnames, SNP_POS, SNP_POS, SNP_ID, gene_id
#'@param TranscriptName
#'Gene or isoform Ensembl ID
#'@param windowSize
#'window around the TSS (e.g. 1000000)
#'@param FDRthreshold
#'FDR threshold (e.g. 0.001)
#'@param BackendData_GenomicFeatures
#'fifteen-core chromatin states from https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/download/
#'@param BackendData_ChromatinStates
#'10 genomic features from the biomaRt (Ensembl) and UCSC
#'@param SNPs
#'List of all SNPs in human genome
#'@return You can find the results in R object under title of 'RESULTsGenomicFeatures' and 'RESULTsChromatinState'.
#'@examples
#'data(eQTL)
#'data(BackendData_GenomicFeatures)
#'data(BackendData_ChromatinStates)
#'data(SNPs)
#'eSNPsEnrichmentAnalysis(eQTL, TranscriptName = 'ENSG00000168310', windowSize=1000000, FDRthreshold = 0.001, BackendData_GenomicFeatures, BackendData_ChromatinStates, SNPs)
#'@export

eSNPsEnrichmentAnalysis <- NULL
eSNPsEnrichmentAnalysis <- function(eQTL, TranscriptName, windowSize, FDRthreshold, BackendData_GenomicFeatures, BackendData_ChromatinStates, SNPs)
{



Reference:
Isar Nassiri, James Gilchrist, Evelyn Lau, Sara Danielli, Hussein Al Mossawi, Jane Cheeseman, Matthew Neville, Julian C Knight, Benjamin P Fairfax. Genetic deter-minants of monocyte splicing are enriched for disease susceptibility loci includingfor COVID-19.

