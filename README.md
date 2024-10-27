Functional Enrichment of Genomic Variants and Variations (FEVV)
==========
* [Introduction](#introduction)
* [Installation](#Installation)
* [Usage](#Usage)
* [Citation](#citation)
<a name="introduction"/>

### Introduction

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
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("SNPlocs.Hsapiens.dbSNP144.GRCh38", quietly = TRUE)) BiocManager::install("SNPlocs.Hsapiens.dbSNP144.GRCh38")
    
library(devtools)
install_github("isarnassiri/FEVV")
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
```

#### Usage

##### A) Functional Enrichment of eQTL SNPs

It is easy to run the enrichment analysis for those who possess little or no programming language skills. To run the enrichment analysis on your computer install and call the FEVV using rStudio, select  `eSNPsEnrichmentAnalysis()` function, and click on the "Run" button at the top right of the Source tab. The enrichment analysis automatically will be performed, including .

```{r,eval=FALSE}

#########################################################################
# Please execute the code in the RStudio IDE (https://www.rstudio.com/) #
#########################################################################

library("FEVV")
InputDir=system.file("extdata", package = "FEVV")

start <- Sys.time()
eSNPsEnrichmentAnalysis(eQTL, TranscriptName = 'ENSG00000168310', windowSize=1000000, FDRthreshold = 0.001, BackendData_GenomicFeatures, BackendData_ChromatinStates, SNPs)
print( Sys.time() - start )

############################################################ 
#  Find the "Interactive QC Report" in the Outputs/ folder #
############################################################

```

# eQTL SNPs as input:

For eSNP enrichment analysis, first, for each gene/transcript, we selected associated eSNPs as foreground (F) and SNPs within the 1Mb window around the TSS as background (B) SNP sets. We used the number of overlaps of foreground (f) and background (b) SNP sets in the genomic feature or chromatin state and calculate the enrichment score (z-score and odds ratio).

Backend Data Genomic Features
Fifteen-core chromatin states has been downloaded from https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/download/ and are available with the package.


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
eSNPsEnrichmentAnalysis(eQTL, TranscriptName, windowSize, FDRthreshold, BackendData_GenomicFeatures, BackendData_ChromatinStates, SNPs)
{



Reference:
Isar Nassiri, James Gilchrist, Evelyn Lau, Sara Danielli, Hussein Al Mossawi, Jane Cheeseman, Matthew Neville, Julian C Knight, Benjamin P Fairfax. Genetic deter-minants of monocyte splicing are enriched for disease susceptibility loci includingfor COVID-19.

