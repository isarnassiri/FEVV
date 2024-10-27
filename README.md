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

library(devtools)
install_github("isarnassiri/FEVV")
```

<a name="usage"/>

### Usage

#### A) Functional Enrichment of eQTL SNPs

It is easy to run the enrichment analysis for those who possess little or no programming language skills. To run the enrichment analysis on your computer install and call the FEVV using rStudio, select  `eSNPsEnrichmentAnalysis()` function, and click on the "Run" button at the top right of the Source tab. The enrichment analysis will be carried out automatically according to the following.

```{r,eval=FALSE}

#########################################################################
# Please execute the code in the RStudio IDE (https://www.rstudio.com/) #
#########################################################################

library("FEVV")
InputDir=system.file("extdata", package = "FEVV")

start <- Sys.time()
eSNPsEnrichmentAnalysis(eQTL, TranscriptName = 'ENSG00000168310', windowSize=1000000, FDRthreshold = 0.001, BackendData_GenomicFeatures, BackendData_ChromatinStates, SNPs)
print( Sys.time() - start )

###################################################################################################### 
#  Find the "RESULTsChromatinState.txt" and "RESULTsGenomicFeatures.txt" files in the ~/FEVV/ folder #
######################################################################################################

```

The 'eQTL' object should have headers that include seqnames, SNP_POS, SNP_POS, SNP_ID, and gene_id for each column. Specifying the transcript name requires the Ensembl ID of the interest transcript. The search space around the transcription start site (TSS) is represented by the window size (e.g. 1000000). The threshold for FDR (FDR) can be described using the term FDR threshold (e.g. 0.001). "BackendData_GenomicFeatures" includes a fifteen-core chromatin states from https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/download/. 'BackendData_ChromatinStates' includes 10 genomic features from BiomaRt (Ensembl) and UCSC. A list of all SNPs in the human genome is represented by the 'SNPs' parameter.

#### B) Functional Enrichment of a query SNP

The functional enrichment of a query SNP will be carried out automatically according to the following.

```{r,eval=FALSE}

#########################################################################
# Please execute the code in the RStudio IDE (https://www.rstudio.com/) #
#########################################################################

library("FEVV")
InputDir=system.file("extdata", package = "FEVV")

start <- Sys.time()
e = "FEVV")
#'querySNPsEnrichmentAnalysis(SNP = 'rs13149699', mafThreshold = 0.039, windowSize = 1000000, BackendData_GenomicFeatures, BackendData_ChromatinStates, vcfMetaData = system.file("extdata", "Genotyping1000_samples_metatadata.txt", package="FEVV"), vcfPATH = 'http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz' )
print( Sys.time() - start )

###################################################################################################### 
#  Find the "RESULTsChromatinState.txt" and "RESULTsGenomicFeatures.txt" files in the ~/FEVV/ folder #
######################################################################################################

```
Specifying the SNP ID (SNP) requires the rs-ID of the interest SNP. The search space around the transcription start site (TSS) is represented by the window size (e.g. 1000000). The threshold for minor allele frequency (maf) can be described using the parameter "mafThreshold" (e.g. 0.039). "BackendData_GenomicFeatures" includes a fifteen-core chromatin states from https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/download/. 'BackendData_ChromatinStates' includes 10 genomic features from BiomaRt (Ensembl) and UCSC. A list of all SNPs in the human genome is represented by the 'SNPs' parameter. vcfMetaData is the metadata for phase 3 of the 1000 Genomes project, and it encompasses sample, pop, super_pop, and gender as a header. ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz contains genetic variation data for chromosome 4 from the 1000 Genomes Project Phase 3. This needs to be replaced based on the genomic location of the quety SNP. Other human chromosome profiles can also be found at the repository site, http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/.

<a name="citation"/>

### Citation:

Isar Nassiri, James Gilchrist, Evelyn Lau, Sara Danielli, Hussein Al Mossawi, Jane Cheeseman, Matthew Neville, Julian C Knight, Benjamin P Fairfax. Genetic deter-minants of monocyte splicing are enriched for disease susceptibility loci includingfor COVID-19.

