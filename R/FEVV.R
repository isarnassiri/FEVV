
#############################################################################################
###################################### eSNP enrichment ######################################
#############################################################################################

#'@export
#'@import GenomicRanges
#'@import data.table
#'@import motifbreakR
#'@import SNPlocs.Hsapiens.dbSNP144.GRCh38
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

  zScore <- function (cntA, totA, cntB, totB) {
    #calculate
    avgProportion <- (cntA + cntB) / (totA + totB)
    probA <- cntA/totA
    probB <- cntB/totB
    SE <- sqrt(avgProportion * (1-avgProportion)*(1/totA + 1/totB))
    zScore <- (probA-probB) / SE
    return (zScore)
  }

  #--------------------------- inputs
  states_bed <- BackendData_ChromatinStates

  chromatin_states <- GRanges(sample = states_bed[[4]],
                             seqnames = gsub('chr','',states_bed[[1]]),
                             ranges = IRanges(states_bed[[2]], states_bed[[3]]),
                             state = states_bed[[4]])

  genomic_regions <- BackendData_GenomicFeatures
  genomic_regions$seqnames <- gsub('chr','',genomic_regions$seqnames)
  genomic_regions <- genomic_regions[which(genomic_regions$seqnames %in% c(1:22, 'X')),]
  genomic_regions <- GRanges(sample = genomic_regions$sample,
                            seqnames = genomic_regions$seqnames,
                            ranges = IRanges(genomic_regions$start, genomic_regions$end),
                            state = genomic_regions$state)

  SNPs_Imputed <- SNPs
  SNPs_Imputed_subject <- GRanges(Rle(as.character(SNPs_Imputed$CHROM)), IRanges(SNPs_Imputed$POS, width=2), rsID = SNPs_Imputed$ID)
  SNPs_Imputed <- as.data.frame(SNPs_Imputed)

  eQTL <- eQTL[eQTL$FDR < FDRthreshold,]
  eQTL <- eQTL[!is.na(eQTL$seqnames),]
  eQTL <- as.data.frame(eQTL)

  #--- select eQTLs
  eQTL_sub <- eQTL[which(eQTL$gene_id == TranscriptName),]
  dim(eQTL_sub)

  #--- Foreground_rsIDs
  Foreground_rsIDs <- GRanges(paste0(eQTL_sub$seqnames,':',eQTL_sub$SNP_POS,'-',eQTL_sub$SNP_POS), rsID = eQTL_sub$SNP_ID)
  length(Foreground_rsIDs)

  #--- Background_rsIDs - select SNPs any a 1Mb windows
  posQuerySNP <- GRanges(paste0(eQTL_sub$seqnames[1],':',eQTL_sub$start[1],'-',eQTL_sub$start[1])) + windowSize
  fo <- findOverlaps(query=posQuerySNP, subject=SNPs_Imputed_subject, type="any")
  Background_rsIDs <- SNPs_Imputed_subject[subjectHits(fo),]

  #--- remove foreground from background
  Background_rsIDs <- Background_rsIDs[-which(Background_rsIDs$rsID %in% Foreground_rsIDs$rsID),]
  length(Background_rsIDs)

  #------------------- chromatin state

  #--- enrich foreground
  fo <- findOverlaps(query=Foreground_rsIDs, subject=chromatin_states, type="any")
  foreground_CS <- chromatin_states[subjectHits(fo),]

  #--- enrich background
  fo <- findOverlaps(query=Background_rsIDs, subject=chromatin_states, type="any")
  Background_CS <- chromatin_states[subjectHits(fo),]

  CSs <- unique(foreground_CS$sample)

  for(t in 1:length(CSs))
  {
    cntA<-length(which(foreground_CS$sample==CSs[t]))
    totA<-length(Foreground_rsIDs)
    cntB<-length(which(Background_CS$sample==CSs[t]))
    totB<-length(Background_rsIDs)

    table <- c(cntA, totA, cntB, totB)
    dim(table)<-c(2,2)
    fishert <- fisher.test(table)

    z_score <- zScore(cntA, totA, cntB, totB)

    temp <- data.frame(gene = TranscriptName, CS = CSs[t], p = fishert$p.value, oddsratio = fishert$estimate, CIl= fishert$conf.int[1], CIh = fishert$conf.int[2], zScore = z_score, FE = cntA, FA = totA, BE = cntB, BA = totB)

    if(t==1){RESULTsChromatinState = temp}else{RESULTsChromatinState = rbind(RESULTsChromatinState, temp)}

  }

  #------------------- genomic features

  #--- enrich foreground
  fo <- findOverlaps(query=Foreground_rsIDs, subject=genomic_regions, type="any")
  foreground_GR = genomic_regions[subjectHits(fo),]

  #--- enrich background
  fo <- findOverlaps(query=Background_rsIDs, subject=genomic_regions, type="any")
  Background_GR = genomic_regions[subjectHits(fo),]

  if(length(foreground_GR)>0)
  {
    GRs <- unique(foreground_GR$sample)
    for(t in 1:length(GRs))
    {
      cntA<-length(which(foreground_GR$sample==GRs[t]))
      totA<-length(foreground_GR)
      cntB<-length(which(Background_GR$sample==GRs[t]))
      totB<-length(Background_GR)

      table <- c(cntA, totA, cntB, totB)
      dim(table)<-c(2,2)
      fishert <- fisher.test(table)

      z_score <- zScore(cntA, totA, cntB, totB)

      temp <- data.frame(gene = TranscriptName, GR = GRs[t], p = fishert$p.value, oddsratio = fishert$estimate, CIl= fishert$conf.int[1], CIh = fishert$conf.int[2], zScore = z_score, FE = cntA, FA = totA, BE = cntB, BA = totB)

      if(t==1){RESULTsGenomicFeatures = temp}else{RESULTsGenomicFeatures = rbind(RESULTsGenomicFeatures, temp)}

    }
  }

  #------------------- save outputs

  Destiny_Folder <- system.file(package = "FEVV")
  Destiny_Folder <- paste(Destiny_Folder, "/RESULTsGenomicFeatures.txt", sep = "")

  write.table(
    RESULTsGenomicFeatures, Destiny_Folder, sep = "\t", row.names = FALSE, quote = FALSE
  )

  Destiny_Folder <- system.file(package = "FEVV")
  Destiny_Folder <- paste(Destiny_Folder, "/RESULTsChromatinState.txt", sep = "")

  write.table(
    RESULTsChromatinState, Destiny_Folder, sep = "\t", row.names = FALSE, quote = FALSE
  )

  print("You can find the results at: ")
  print(system.file(package="FEVV"))
}



#########################################################################################################
###################################### single query SNP enrichment ######################################
#########################################################################################################

#'@export
#'@import SNPlocs.Hsapiens.dbSNP144.GRCh38
#'@import BSgenome.Hsapiens.UCSC.hg38
#'@import motifbreakR
#'@import snpStats
#'@import GenomicRanges
#'@import data.table
#'@import biomaRt
#'@import GenomicRanges
#'@import VariantAnnotation
#'@name querySNPsEnrichmentAnalysis
#'@title Single query SNP enrichment
#'@description a query SNP as the input: a query SNP as input: first, the 1MB window is defined for the query SNP, and the relevant variants captured from the 1000 genomes data in the European population background. Next, we split the SNPs based upon a LD (R?? = 0.8 and MAF = 0.01) to foreground (query SNP and its LD proxies) and background SNPs sets. We used the number of overlaps of foreground (f) and background (b) SNP sets in the genomic feature or chromatin state and calculate the enrichment score (z-score and odds ratio).
#'@author {Isar Nassiri, Benjamin Fairfax}
#'@param windowSize
#'window around the TSS (e.g. 1000000)
#'@param BackendData_GenomicFeatures
#'fifteen-core chromatin states from https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/coreMarks/jointModel/final/download/
#'@param BackendData_ChromatinStates
#'10 genomic features from the biomaRt (Ensembl) and UCSC
#'@param SNP
#'a query SNP
#'@param mafThreshold
#'minor allele frequency (MAF) threshold
#'@param vcfMetaData
#'meta data of 1000 Genomes project (phase3) including 'sample', 'pop', 'super_pop', and 'gender' as a header.
#'@param vcfPATH
#'Pathe to the vcf file.
#'@return You can find the results in R object under title of 'RESULTsGenomicFeatures' and 'RESULTsChromatinState'.
#'@examples
#'data(BackendData_GenomicFeatures)
#'data(BackendData_ChromatinStates)
#'Destiny_Folder <- system.file(package = "FEVV")
#'querySNPsEnrichmentAnalysis(SNP = 'rs13149699', mafThreshold = 0.039, windowSize = 1000000, BackendData_GenomicFeatures, BackendData_ChromatinStates, vcfMetaData = system.file("extdata", "Genotyping1000_samples_metatadata.txt", package="FEVV"), vcfPATH = 'http://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz' )
#'@export

querySNPsEnrichmentAnalysis <- NULL
querySNPsEnrichmentAnalysis <- function(SNP, mafThreshold, windowSize, BackendData_GenomicFeatures, BackendData_ChromatinStates, vcfMetaData, vcfPATH)
{

  #########################
  ### Helpers functions ###
  #########################

  GetVariantsInWindow <- function(file, position, genome = "hg19", type = "vcf") {
    if (tolower(type) == "vcf") {
      if(!missing(position)) {
        if(is(position, "GRanges")) {
          params <- ScanVcfParam(which = position)
        } else if(position == "all") {
          params <- ScanVcfParam()
        } else if(is(position, "character")) {
          position <- as(position, "GRanges")
          params <- ScanVcfParam(which = position)
        } else {
          stop("I don't understand the position argument, must be unset, GRanges, or 'all'")
        }
      } else {
        stop("without a position argument, the full vcf will be imported into memory,\n",
             "this can be a very expensive operation. Set the position arg to 'all' to allow")
      }
      vcf <- TabixFile(file)
      variants <- getFILE(file, GetVariantsInWindowVCF, params, genome, N.TRIES = 3L)
      return(variants)
    } else if (tolower(type) == "bed") {
      if(!missing(position) & position != 'all') {
        variants <- import.bed(file, which = position, genome = genome)
      } else {
        variants <- import.bed(file, genome = genome)
      }
      return(variants)
    } else {
      stop("type ", type )
    }
  }

  SetPopulation <- function(vcf, sample_sheet) {
    population <- colData(vcf)
    if (!("Samples" %in% colnames(population))) {
      stop("vcf does not contain sample information")
    }
    population$Samples <- rownames(population)
    sample.col <- which(grepl(pattern = "sample",
                              x = colnames(sample_sheet),
                              ignore.case = TRUE))
    if (length(sample.col) > 0L) {
      sample_sheet <- DataFrame(sample_sheet)
      colnames(sample_sheet)[sample.col] <- "Samples"
      population <- S4Vectors::merge(population, sample_sheet, all.x = TRUE, all.y = FALSE, by = "Samples")
      rownames(population) <- population$Samples
      colData(vcf) <- population
      metadata(vcf)$source <- "funciVar"
      metadata(vcf)$ld <- "none"
      return(vcf)
    } else {
      stop("sample sheet does not contain sample column")
    }
  }

  CalcLD <- function(vcf, index, population, return = "valid", force = TRUE) {
    ## check input
    if(!is(vcf, "VCF")) {
      stop("vcf object must be of the class VCF")
    }
    if (index %in% rownames(vcf)) {
      if (!isSNV(vcf[rownames(vcf) %in% index, ])) {
        stop("funciVar can only calculate LD for SNVs at this time")
      }
    } else {
      stop("index snp ", index, " not found in current vcf")
    }
    if(missing(population)) {
      population <- "ALL"
      vcf.snv <- isSNV(vcf)
      if (return == "valid") {
        if (any(!vcf.snv)) {
          vcf <- vcf[vcf.snv, ]
        }
        vcf <- vcf[isSNV(vcf, singleAltOnly = TRUE), ]
      }
    } else {
      samples <- as.data.frame(colData(vcf))
      samples.col <- which(samples == population, arr.ind = TRUE)
      if (nrow(samples.col) > 0L) {
        ## metadata
        samples.col <- samples.col[1, "col"]
        pop.samples <- rownames(samples)[grepl(population, samples[, samples.col])]
        vcf <- vcf[, pop.samples]
        vcf.snv <- isSNV(vcf)
        if (return == "valid") {
          if (any(!vcf.snv)) {
            vcf <- vcf[vcf.snv, ]
          }
          vcf <- vcf[isSNV(vcf, singleAltOnly = TRUE), ]
        }
      } else {
        stop("your population was not found in your vcf.\nUse setPopulation() to add your sample descriptions, and double check that ",
             population, " is present")
      }
    }
    if(nrow(vcf) >= 1L) {
      vcf.ranges <- rowRanges(vcf)[, c("REF", "ALT")]
      mcols(vcf.ranges) <- c(mcols(vcf.ranges), snpSummary(vcf)[, c("a0Freq", "a1Freq", "HWEpvalue")])
      colnames(mcols(vcf.ranges)) <- c("ref", "alt", "refAlleleFreq", "altAlleleFreq", "HWEpvalue")
      mcols(vcf.ranges)[, "indexSNP"] <- index
      mcols(vcf.ranges)[, "population"] <- population
      mcols(vcf.ranges)[, "distanceToIndex"] <- abs(start(vcf.ranges[index, ]) - start(vcf.ranges))
      ## genotype
      vcf.geno <- genotypeToSnpMatrix(vcf)$genotypes
      vcf.ld <- ld(vcf.geno, vcf.geno[, index], stats = c('D.prime', 'R.squared'))
      mcols(vcf.ranges)$D.prime <- vcf.ld$D.prime[, 1]
      mcols(vcf.ranges)$R.squared <- vcf.ld$R.squared[, 1]
      rowRanges(vcf) <- vcf.ranges
      ## xXx rowData(vcf) <- cbind(rowData(vcf), mcols(vcf.ranges))
      if ("none" %in% metadata(vcf)$ld) {
        metadata(vcf)$ld <- index
      } else if(force) {
        metadata(vcf)$ld <- index
      } else if(!("ld" %in% names(metadata(vcf)))) {
        metadata(vcf)$source <- "funciVar"
        metadata(vcf)$ld <- index
      } else {
        stop("ld has already been calculated on this object for index snp: ", index)
      }
    }
    return(vcf)
  }

  GetBioFeatures <- function(files, genome) {
    if (length(files) < 1L) return(NULL)
    good.files <- sapply(files, function(x) file.exists(x))
    if (sum(good.files) < length(files)) {
      if (sum(!good.files) <= 5L) {
        stop(paste("cannot find the following", sum(!good.files), "files:\n"),
             paste0(files[!good.files], "\n"))
      } else {
        stop(paste("cannot find", sum(!good.files), "files. Some of which are:\n"),
             paste0(head(files[!good.files], n = 5L), "\n"))
      }
    } else {
      genome <- tryCatch(Seqinfo(genome = genome), error = NULL)
      bed.list <- lapply(files,
                         function(file, s.info) {
                           if (grepl(".narrowPeak", file, ignore.case = TRUE)) {
                             col.types <- cols_only(chr = col_character(),
                                                    start = col_integer(),
                                                    end = col_integer(),
                                                    name = col_character(),
                                                    score = col_integer(),
                                                    strand = col_character(),
                                                    signalValue = col_number())
                             col.names <- c("chr", "start", "end", "name", "score", "strand", "signalValue")
                             col.numbers <- c(1:7)
                           } else {
                             col.types <- cols_only(chr = col_character(),
                                                    start = col_integer(),
                                                    end = col_integer())
                             col.names <- c("chr", "start", "end")
                             col.numbers <- c(1:3)
                           }
                           if (any(grepl(".gz$", file))) {
                             xf <- suppressMessages(suppressWarnings(read_tsv(file = file, col_names = col.names,
                                                                              col_types = col.types,
                                                                              progress = FALSE)))
                           } else {
                             xf <- fread(input = file, sep = "\t", header = FALSE,
                                         select = col.numbers, skip = "chr",
                                         col.names = col.names,
                                         encoding = "UTF-8",
                                         stringsAsFactors = FALSE,
                                         data.table = FALSE, showProgress = FALSE)
                           }
                           if (ncol(xf) < 7) xf$signalValue <- NA
                           xf <- try(GRanges(seqnames = xf$chr,
                                             ranges = IRanges(start = xf$start + 1L,
                                                              end = xf$end),
                                             strand = "*",
                                             feature = base::rep.int(basename(file), nrow(xf)),
                                             signalValue = xf$signalValue,
                                             seqinfo = s.info))
                           return(xf)
                         }, s.info = genome)
    }
    if (is.null(bed.list)) {
      return(GRangesList())
    } else {
      bed.list <- GRangesList(bed.list)
      names(bed.list) <- basename(files)
      return(bed.list)
    }
  }

  GetSegmentations <- function(files) {
    bed.list <- sapply(files, function(file) {
      bed <- import.bed(file)
    })
    bed.names <- sapply(bed.list, function(bed) {
      bed.name <- tryCatch(bed@trackLine@name, error = NA)
    })
    bed.names[is.na(bed.names)] <- files[is.na(bed.names)]
    bed.list <- Map(format.bed, bed.list, bed.names)
    if (is.null(bed.list)) {
      return(GRangesList())
    } else {
      bed.list <- GRangesList(bed.list)
      names(bed.list) <- bed.names
      return(bed.list)
    }
  }

  format.bed <- function(bed, name) {
    bed.m <- mcols(bed)
    mcols(bed) <- NULL
    mcols(bed)$sample <- name
    mcols(bed)$state <- bed.m$name
    return(bed)
  }

  SplitVcfLd <- function(vcf, ld = list(metric = "R.squared", cutoff = 0.8, maf = 0.01), strict.subset = TRUE) {
    if (!is(vcf, "VCF")) {
      stop("parameter vcf must be a VCF object")
    }
    if (!all(names(ld) %in% c("metric", "cutoff", "maf"))) {
      stop("parameter ld must contain the fields 'metric', 'cutoff', and 'maf'")
    }
    if (!(ld[["metric"]] %in% colnames(mcols(rowRanges(vcf))))) {
      stop("in argument ld$metric: '", ld[['metric']],"' must be present in rowRanges(vcf)")
    }
    if (!(ld[["maf"]] >= 0 && ld[["maf"]] <= 0.5)) {
      stop("in argument ld$maf: '", ld[['maf']], "' must be between the values of 0 and 0.5")
    }
    orient <- mcols(rowRanges(vcf))[, "altAlleleFreq"] < mcols(rowRanges(vcf))[, "refAlleleFreq"]
    filter <- mcols(rowRanges(vcf))[, "altAlleleFreq"] >= ld[["maf"]]
    filter[!orient] <- mcols(rowRanges(vcf))[, "refAlleleFreq"][!orient] >= ld[["maf"]]
    fg.filter <- filter & mcols(rowRanges(vcf))[, ld[["metric"]]] >= ld[["cutoff"]]
    if (strict.subset) {
      bg.filter <- fg.filter | mcols(rowRanges(vcf))[, ld[["metric"]]] < ld[["cutoff"]]
    } else {
      bg.filter <- filter & !mcols(rowRanges(vcf))[, ld[["metric"]]] >= ld[["cutoff"]]
    }
    metadata(vcf)$strict.subset <- strict.subset
    bg.filter <- fg.filter | mcols(rowRanges(vcf))[, ld[["metric"]]] < ld[["cutoff"]]
    return(list(fg = vcf[fg.filter & !is.na(fg.filter), ], bg = vcf[bg.filter & !is.na(bg.filter), ]))
  }


  CalculateEnrichment <- function(variants, features, feature.type = "biofeatures",
                                  CI = 0.95, prior = c(a = 0.5, b = 0.5),
                                  strict.subset = "guess", return.overlaps = FALSE) {
    # set strict.subset
    if (strict.subset == "guess") {
      test.sub <- c(metadata(variants$fg)$strict.subset, metadata(variants$bg)$strict.subset)
      if (all(test.sub, !is.null(test.sub))) {
        strict.subset <- TRUE
      } else if (any(metadata(variants$fg)$strict.subset, metadata(variants$bg)$strict.subset)) {
        stop("foreground and background disagree about whether fg is a strict subset of bg\n",
             "create foreground and background with SplitVcfLD")
      } else if (all(rownames(variants$fg) %in% rownames(variants$bg))) {
        strict.subset <- TRUE
      } else {
        strict.subset <- FALSE
      }
    } else if (!is.logical(strict.subset)) {
      stop("strict.subset must be one of 'guess', TRUE, or FALSE")
    }
    # features
    if(!is(variants, "list")) stop("variants must be an object of class list")
    if (!missing(features)) {
      if (feature.type == "biofeatures") {
        if(inherits(features, "GRangesList")) features <- unlist(features, use.names = FALSE)
        variants$fg <- SetOverlaps(variants$fg, features)
        variants$bg <- SetOverlaps(variants$bg, features)
        fg.features <- colnames(mcols(ShowOverlaps(variants$fg)))
        bg.features <- colnames(mcols(ShowOverlaps(variants$bg)))
        all.features <- union(fg.features, bg.features)
        ## first fg
        if (any(!(is.element(all.features, fg.features)))) {
          if (is(variants$fg, "VCF")) {
            mcols(rowRanges(variants$fg))[, all.features[!(is.element(all.features, fg.features))]] <- 0L
          } else if (is(variants$fg, "GRanges")) {
            mcols(variants$fg)[, all.features[!(is.element(all.features, fg.features))]] <- 0L
          }
        }
        if (is(variants$bg, "VCF")) {
          mcols(rowRanges(variants$fg)) <- cbind(mcols(rowRanges(variants$fg))[, 1:metadata(variants$fg)$overlap.offset-1], mcols(rowRanges(variants$fg))[, all.features, drop = FALSE])
        } else if (is(variants$fg, "GRanges")) {
          mcols(variants$fg) <- cbind(mcols(variants$fg)[, 1:attributes(variants$fg)$metadata$overlap.offset-1], mcols(variants$fg)[, all.features, drop = FALSE])
        }
        ## then bg
        if (any(!(is.element(all.features, bg.features)))) {
          if (is(variants$bg, "VCF")) {
            mcols(rowRanges(variants$bg))[, all.features[!(is.element(all.features, bg.features))]] <- 0L
          } else if (is(variants$bg, "GRanges")) {
            mcols(variants$bg)[, all.features[!(is.element(all.features, bg.features))]] <- 0L
          }
        }
        if (is(variants$bg, "VCF")) {
          mcols(rowRanges(variants$bg)) <- cbind(mcols(rowRanges(variants$bg))[, 1:metadata(variants$bg)$overlap.offset-1], mcols(rowRanges(variants$bg))[, all.features, drop = FALSE])
        } else if (is(variants$bg, "GRanges")) {
          mcols(variants$bg) <- cbind(mcols(variants$bg)[, 1:attributes(variants$bg)$metadata$overlap.offset-1], mcols(variants$bg)[, all.features, drop = FALSE])
        }
      }
    }
    # feature.type = "biofeatures"
    if (feature.type == "biofeatures") {
      enrichment <- enrich.features(fg = variants$fg, bg = variants$bg, CI = CI, prior = prior, strict.subset = strict.subset)
      if(return.overlaps) {
        if(all(names(mcols(variants$fg)) == names(mcols(variants$bg)))) {
          list.fun <- GenomicRanges::GRangesList
        } else {
          list.fun <- S4Vectors::List
        }
        if (is(variants$fg, "GRanges") & is(variants$bg, "GRanges")) {
          overlaps <- list.fun(foreground.overlaps = variants$fg, background.overlaps = variants$bg)
        } else {
          overlaps <- list.fun(foreground.overlaps = rowRanges(variants$fg), background.overlaps = rowRanges(variants$bg))
        }
        return(list(overlaps = overlaps, enrichment = enrichment))
      } else {
        return(enrichment)
      }
    } else if (feature.type == "segmentations") {
      if (missing(features)) {
        stop("include segmentations as the 'features' argument")
      } else {
        if (is(variants$fg, "GRanges") & is(variants$bg, "GRanges")) {
          search.range <- GenomicRanges::union(range(variants$fg), range(variants$bg))
        } else {
          search.range <- GenomicRanges::union(range(rowRanges(variants$fg)), range(rowRanges(variants$bg)))
        }
        if (is(features, "GRangesList")) {
          features <- unlist(features, use.names = FALSE)
        }
        # nfeatures <- gaps(features)
        # nfeatures <- nfeatures[strand(nfeatures) == "*"]
        # mcols(nfeatures)$sample <- "none"
        # mcols(nfeatures)$state <- "unclassified"
        # features <- c(features, nfeatures)
        features <- keepSeqlevels(subsetByOverlaps(features, search.range), seqlevelsInUse(search.range))
        enrichment <- enrich.segments(fg = variants$fg,
                                      bg = variants$bg,
                                      features = features,
                                      CI = CI,
                                      prior = prior,
                                      strict.subset = strict.subset,
                                      return.overlaps = return.overlaps)
      }
    } else {
      stop("feature.type: ", feature.type, "is not availible; try 'biofeatures' or 'segmentations'")
    }
    return(enrichment)
  }


  SetOverlaps <- function(variants, features) {
    nfeatures <- gaps(features)
    nfeatures <- nfeatures[strand(nfeatures) == "*"]
    if(all(names(mcols(features)) == c("sample", "state"))) {
      mcols(nfeatures)$sample <- "none"
      mcols(nfeatures)$state <- "unclassified"
      id.name <- "sample"
    } else if (all(names(mcols(features)) == c("feature", "signalValue"))) {
      mcols(nfeatures)$feature <- "none"
      mcols(nfeatures)$signalValue <- 0L
      id.name <- "feature"
    }
    features <- c(features, nfeatures)
    if (is.null(names(variants))) {
      names(variants) <- make.unique(as.character(variants))
    }
    overlaps <- findOverlaps(variants, features, ignore.strand = TRUE)
    overlaps <- data.frame(from = names(variants)[from(overlaps)],
                           to = mcols(features)[to(overlaps), id.name],
                           stringsAsFactors = FALSE)
    resmatrix <- make.overlap.matrix(overlaps)
    resmatrix <- resmatrix[names(variants), !grepl("^none$",
                                                   colnames(resmatrix)), drop = FALSE]
    if(ncol(resmatrix) > 0L) {
      if(is(variants, "GRanges")) {
        overlap.offset <- ncol(mcols(variants)) + 1
        mcols(variants) <- cbind(mcols(variants), DataFrame(resmatrix))
        attributes(variants)$metadata$overlap.offset <- overlap.offset
      } else if (is(variants, "VCF")) {
        overlap.offset <- (ncol(mcols(rowRanges(variants))) + 1) - 4
        mcols(rowRanges(variants)) <- cbind(mcols(rowRanges(variants)), DataFrame(resmatrix))
        metadata(variants)$overlap.offset <- overlap.offset
      } else {
        stop("variants must be either a GRanges or VCF object, your object is: ", class(variants))
      }
    } else {
      if(is(variants, "GRanges")) {
        metadata(variants)$overlap.offset <- ncol(mcols(variants)) + 1
      } else if (is(variants, "VCF")) {
        metadata(variants)$overlap.offset <- ncol(mcols(rowRanges(variants))) + 1
      } else {
        stop("variants must be either a GRanges or VCF object, your object is: ", class(variants))
      }
    }
    return(variants)
  }

  ShowOverlaps <- function(variants, feature = NULL) {
    if (is(variants, "GRanges")) {
      offset <- attributes(variants)$metadata$overlap.offset
    } else if (is(variants, "VCF")) {
      offset <- metadata(variants)$overlap.offset
      variants <- rowRanges(variants)
      variants <- variants[, 1:(ncol(mcols(variants)) - 4)]
    } else {
      stop("variants must be either a VCF or GRanges object")
    }
    if (ncol(mcols(variants)) < offset) {
      return(DataFrame())
    }
    if (is.null(feature)) {
      if (!is.null(offset)) {
        return(variants[, offset:ncol(mcols(variants))])
      } else {
        stop("no overlap data is availible, run SetOverlaps() first")
      }
    } else {
      if (any(feature %in% colnames(mcols(rowRanges(variants))))) {
        return(variants[, feature])
      } else {
        stop("no overlap data is availible for feature search: ", feature)
      }
    }
  }

  PlotEnrichment <- function(variant.enrichment, value = "difference", block1 = NULL, block2 = NULL, color.by = NULL, colors = NULL, ncol = 1) {

    ## check args
    if (!is(variant.enrichment, "data.frame"))
      stop("variant.enrichment must be a data.frame")
    if (!is.null(block1)) {
      if (!(block1 %in% colnames(variant.enrichment)))
        stop("The block1 argument must either be NULL or a column of variant.enrichment.")
    }
    if (!is.null(block2)) {
      if (!(block2 %in% colnames(variant.enrichment)))
        stop("The block2 argument must either be NULL or a column of variant.enrichment.")
    }
    if (!is.null(color.by)) {
      if (!(color.by %in% colnames(variant.enrichment)))
        stop("The color.by argument must either be NULL or a column of variant.enrichment.")
    }
    if (!(value %in% c("log.odds.ratio", "difference"))) {
      stop("The y argument must be either 'log.odds.ratio' or 'difference'")
    }
    if(value == "log.odds.ratio") value.name <- "log1p(odds ratio)"
    if(value == "difference") value.name <- "difference of foreground \nadbackground distributions"
    if (is.null(color.by)) {
      variant.enrichment$color <- "#ececec"
      variant.enrichment[variant.enrichment$significant, "color"] <- "#ff0000"
      color.is <- "self"
    } else {
      color.vals <- variant.enrichment[, color.by]
      names(color.vals) <- rownames(variant.enrichment)
      color.name <- color.by
      if(!is.null(colors)) {
        if(length(color.vals) == length(colors)) {
          variant.enrichment$color <- colors
          color.is <- "self"
        } else {
          stop("arg 'colors' is not the same length as color.by")
        }
      }
      if (all(grepl("^#", color.vals) & (str_length(color.vals) == 7L))) {
        variant.enrichment$color <- color.vals
        variant.enrichment[!variant.enrichment$significant, "color"] <- alpha(variant.enrichment[!variant.enrichment$significant, "color"], 0.5)
        color.is <- "self"
      } else if (is.character(color.vals) | is.factor(color.vals) | is.logical(color.vals)) {
        color.vals <- as.factor(color.vals)
        if(nlevels(color.vals) < 10L) {
          levels(color.vals) <- brewer.pal(nlevels(color.vals), "Set1")
        } else if(nlevels(color.vals) < 20L) {
          kelly.colours <- c("gray95", "gray13", "gold2", "plum4",
                             "darkorange1", "lightskyblue2", "firebrick",
                             "burlywood3", "gray51", "springgreen4", "lightpink2",
                             "deepskyblue4", "lightsalmon2", "mediumpurple4",
                             "orange", "maroon", "yellow3", "brown4",
                             "yellow4", "sienna4", "chocolate", "gray19")
          levels(color.vals) <- kelly.colours[3:(nlevels(color.vals)+2)]
        } else {
          levels(color.vals) <- colorRampPalette(brewer.pal(9, "Spectral"))(nlevels(color.vals))
        }
        variant.enrichment$color <- as.character(color.vals)
        variant.enrichment[!variant.enrichment$significant, "color"] <- alpha(variant.enrichment[!variant.enrichment$significant, "color"], 0.5)
        color.is <- "self"
      } else {
        variant.enrichment$color <- variant.enrichment[, color.by]
        color.is <- "cont"
      }
    }

    ep <- ggplot(variant.enrichment, aes_string(x = "sample", y = value, group = "color")) + theme_minimal()

    if (n_distinct(variant.enrichment$sample) <= 7L) {
      if(value == "log.odds.ratio") {
        ep <- ep + geom_crossbar(aes(ymin = log.odds.lower, ymax = log.odds.upper, color = color), fatten = 2, width = 0.8)
      } else {
        ep <- ep + geom_crossbar(aes(ymin = lower, ymax = upper, color = color), fatten = 2, width = 0.8)
      }
    } else {
      if(value == "log.odds.ratio") {
        ep <- ep + geom_pointrange(aes(ymin = log.odds.lower, ymax = log.odds.upper, color = color), fatten = 0.5)
      } else {
        ep <- ep + geom_pointrange(aes(ymin = lower, ymax = upper, color = color), fatten = 0.5)
      }
    }
    ep <- ep + scale_x_discrete(name = "Sample")
    if(value == "log.odds.ratio") {
      my.max <- round(max(variant.enrichment$odds.upper), 2) + 0.01
      my.min <- round(min(variant.enrichment$odds.lower), 2) - 0.01
      ep <- ep + scale_y_continuous(name = value.name)
      ep <- ep + geom_hline(yintercept = 0, color = "#c4c4c4", alpha = 0.5) +
        annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=log1p(1), alpha=0.1, fill="black")
    } else {
      my.max <- round(max(variant.enrichment$upper), 2) + 0.01
      my.min <- round(min(variant.enrichment$lower), 2) - 0.01
      ep <- ep + scale_y_continuous(name = value.name, breaks = round(seq(from = my.min, to = my.max, length.out = 5), 2))
      ep <- ep + geom_hline(yintercept = 0, color = "#c4c4c4", alpha = 0.5) +
        annotate("rect", xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0, alpha=0.1, fill="black")
    }

    if (all(!is.null(block1), !is.null(block2))) {
      ep <- ep + theme(axis.text.x = element_blank(),
                       legend.position="bottom",
                       panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
                       strip.text.y = element_text(angle = 0, hjust = 0, vjust = 0.5, size = rel(1.5)),
                       strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = rel(1.5)))
    } else {
      ep <- ep + theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5, size = rel(1.5)),
                       legend.position="bottom",
                       panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(),
                       strip.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1, size = rel(1.5)))
    }
    if(!is.null(block1) & !is.null(block2)) {
      my.form <- as.formula(paste(block1, "~", block2))
      ep <- ep + facet_grid(my.form, scales = "free_x", space = "free_x", switch = "x")
    } else if (!is.null(block1) & is.null(block2)) {
      my.form <- as.formula(paste("~", block1))
      ep <- ep + facet_wrap(my.form, ncol = ncol, strip.position = ifelse(ncol > 1L, "top", "right"))
    } else if (is.null(block1) & !is.null(block2)) {
      my.form <- as.formula(paste("~", block2))
      ep <- ep + facet_wrap(my.form, ncol = ncol, strip.position = ifelse(ncol > 1L, "top", "right"))
    }

    if(color.is == "cont") {
      ep <- ep + scale_color_brewer(palette = "Greens")
    } else {
      ep <- ep + scale_color_identity()
    }
    plot(ep)
    return(invisible(ep))
  }

  GetVariantsInWindowVCF <- function(file, param, genome) {
    vcf <- readVcf(file = file, genome = genome, param = param)
    if (nrow(vcf) < 1L) stop("no variants found in interval; \n this is sometimes an error in fetching remote file, try with a local vcf")
    vcf <- keepSeqlevels(vcf, seqlevelsInUse(vcf))
    seqlevelsStyle(vcf) <- "UCSC"
    metadata(vcf)$overlap.offset <- NA
    return(vcf)
  }

  getFILE <- function(FILE, FUN, ..., N.TRIES=1L) {
    N.TRIES <- as.integer(N.TRIES)
    stopifnot(length(N.TRIES) == 1L, !is.na(N.TRIES))

    while (N.TRIES > 0L) {
      result <- tryCatch(FUN(FILE, ...), error = identity)
      if (!inherits(result, "error"))
        break
      N.TRIES <- N.TRIES - 1L
    }

    if (N.TRIES == 0L) {
      stop("'getFILE()' failed:",
           "\nFILE: ", FILE,
           "\nerror: ", conditionMessage(result))
    }

    result
  }

  make.overlap.matrix <- function(overlaps) {
    overlap.table <- table(overlaps)
    output.matrix <- matrix(overlap.table,
                            nrow = nrow(overlap.table),
                            dimnames = list(rownames(overlap.table),
                                            colnames(overlap.table)))

    if (!is(output.matrix, "matrix")) {
      dim(output.matrix) <- c(length(output.matrix), 1)
      dimnames(output.matrix) <- list(rownames(overlap.table),
                                      colnames(overlap.table))
    }
    output.matrix[output.matrix > 1L] <- 1L
    return(output.matrix)
  }

  enrich.segments <- function(fg, bg, features, CI, prior, strict.subset, return.overlaps) {
    if (is(features, "GRangesList")) {
      features <- unlist(features, use.names = FALSE)
    }
    if (all(c("sample", "state") %in% colnames(mcols(features)))) {
      states <- unique(mcols(features)$state)
      if (length(states) > 1L) {
        myapply <- pblapply
      } else {
        myapply <- lapply
      }
      enrichment <- myapply(states, function(state,
                                             local.features = features,
                                             local.fg = fg,
                                             local.bg = bg,
                                             local.CI = CI,
                                             local.prior = prior,
                                             local.strict.subset = strict.subset,
                                             local.return.overlaps = return.overlaps) {
        local.features <- local.features[mcols(local.features)$state %in% state, ]
        local.fg <- SetOverlaps(local.fg, local.features)
        local.bg <- SetOverlaps(local.bg, local.features)
        ## Equalize feature columns
        fg.features <- colnames(mcols(ShowOverlaps(local.fg)))
        bg.features <- colnames(mcols(ShowOverlaps(local.bg)))
        all.features <- union(fg.features, bg.features)
        ## first fg
        if (is(local.fg, "VCF")) {
          if(any(!(is.element(all.features, fg.features)))) {
            mcols(rowRanges(local.fg))[, all.features[!(is.element(all.features, fg.features))]] <- 0L
          }
          mcols(rowRanges(local.fg)) <- cbind(mcols(rowRanges(local.fg))[, 1:metadata(local.fg)$overlap.offset - 1],
                                              mcols(rowRanges(local.fg))[, all.features, drop = FALSE])
        } else if (is(local.fg, "GRanges")) {
          if (any(!(is.element(all.features, fg.features)))) {
            mcols(local.fg)[, all.features[!(is.element(all.features, fg.features))]] <- 0L
          }
          mcols(local.fg) <- cbind(mcols(local.fg)[, 1:attributes(local.fg)$metadata$overlap.offset - 1],
                                   mcols(local.fg)[, all.features, drop = FALSE])
        }
        ## then bg
        if (is(local.bg, "VCF")) {
          if(any(!(is.element(all.features, bg.features)))) {
            mcols(rowRanges(local.bg))[, all.features[!(is.element(all.features, bg.features))]] <- 0L
          }
          mcols(rowRanges(local.bg)) <- cbind(mcols(rowRanges(local.bg))[, 1:metadata(local.bg)$overlap.offset - 1],
                                              mcols(rowRanges(local.bg))[, all.features, drop = FALSE])
        } else if (is(local.bg, "GRanges")) {
          if (any(!(is.element(all.features, bg.features)))) {
            mcols(local.bg)[, all.features[!(is.element(all.features, bg.features))]] <- 0L
          }
          mcols(local.bg) <- cbind(mcols(local.bg)[, 1:attributes(local.bg)$metadata$overlap.offset - 1],
                                   mcols(local.bg)[, all.features, drop = FALSE])
        }
        enrich <- enrich.features(fg = local.fg,
                                  bg = local.bg,
                                  CI = local.CI,
                                  prior = local.prior,
                                  strict.subset = local.strict.subset)
        if (local.return.overlaps) {
          return(list(e = enrich, fg.o = ShowOverlaps(local.fg), bg.o = ShowOverlaps(local.bg)))
        } else {
          return(enrich)
        }
      })
      if(return.overlaps) {
        overlaps <- lapply(enrichment, function(x) {
          overlaps <- GRangesList(foreground.overlaps = x$fg.o, background.overlaps = x$bg.o)
          return(overlaps)
        })
        names(overlaps) <- states
        enrichment <- lapply(enrichment, function(x) {
          return(x$e)
        })
      }
      names(enrichment) <- states
      for (state in states) {
        enrichment[[state]]$state <- state
      }
      enrichment <- do.call("rbind", enrichment)
    } else {
      stop("Segmentations must have mcols with fields 'sample' and 'state'")
    }
    if(return.overlaps) {
      return(list(overlaps = overlaps, enrichment = enrichment))
    } else {
      return(enrichment)
    }
  }


  enrich.features <- function(fg, bg, CI, prior, strict.subset) {
    ## foreground overlaps
    fg.over.matrix <- as.matrix(mcols(ShowOverlaps(fg)))
    rownames(fg.over.matrix) <- rownames(fg)
    ## background overlaps
    bg.over.matrix <- as.matrix(mcols(ShowOverlaps(bg)))
    rownames(bg.over.matrix) <- rownames(bg)

    ## start enrichment
    sample.size <- dim(fg.over.matrix)[[1]]
    sample.stats <- colSums(fg.over.matrix)

    if (strict.subset) {
      total.size <- dim(bg.over.matrix)[[1]]
      total.stats <- colSums(bg.over.matrix)
    } else {
      total.size <- dim(bg.over.matrix)[[1]] + sample.size
      total.stats <- colSums(bg.over.matrix) + sample.stats
    }
    if (exists("prior") & !is.null(prior)) {
      a <- prior[["a"]]
      b <- prior[["b"]]
    } else {
      stop("no argument present for prior, needed for simulation")
    }
    enrichment <- data.frame(sample = character(),
                             fg.ratio = numeric(),
                             bg.ratio = numeric(),
                             probability = numeric(),
                             difference = numeric(),
                             lower = numeric(),
                             upper = numeric(),
                             fg.success = numeric(),
                             fg.total = numeric(),
                             bg.success = numeric(),
                             bg.total = numeric(),
                             stringsAsFactors = FALSE)
    CI <- (1 - CI)/2
    for (my.sample in seq_along(colnames(bg.over.matrix))) {
      total.success <- total.stats[my.sample]
      sample.success <- sample.stats[my.sample]

      n1 <- sample.size
      y1 <- sample.success
      n2 <- total.size - sample.size
      y2 <- total.success - sample.success
      # SIMULATION
      I = 100000 # simulations
      theta1 = rbeta(I, y1 + a, (n1 - y1) + b)
      theta2 = rbeta(I, y2 + a, (n2 - y2) + b)
      diff = theta1 - theta2

      # OUTPUT
      quantiles <- try(quantile(diff,c(CI, 0.5, 1 - CI)))
      #if(inherits(quantiles, "try-error")) browser()
      fisher.res <- fisher.test(matrix(c(sample.success,
                                         total.success - sample.success,
                                         sample.size - sample.success,
                                         (total.size - total.success) - sample.size + sample.success),
                                       2, 2),
                                conf.level = CI)

      result <- data.frame(sample = colnames(bg.over.matrix)[my.sample],
                           fg.ratio = sample.success/sample.size,
                           bg.ratio = total.success/total.size,
                           probability = mean(theta1 > theta2),
                           log.odds.ratio = log1p(fisher.res$estimate),
                           log.odds.lower = log1p(fisher.res$conf.int[[1]]),
                           log.odds.upper = log1p(fisher.res$conf.int[[2]]),
                           difference = quantiles[2],
                           lower = quantiles[1],
                           upper = quantiles[3],
                           fg.success = sample.success,
                           fg.total = sample.size,
                           bg.success = total.success,
                           bg.total = total.size,
                           stringsAsFactors = FALSE)
      enrichment <- rbind(enrichment, result)
    }
    enrichment$significant <- FALSE
    enrichment[enrichment$probability > (1 - CI) | enrichment$probability < CI, "significant"] <- TRUE
    return(enrichment)
  }

  clump_data <- function(dat, clump_kb=10000, clump_r2=0.001, clump_p1=1, clump_p2=1)
  {
    .Deprecated("ieugwasr::ld_clump()")

    if(missing(clump_r2))
    {
    }
    if(!is.data.frame(dat))
    {
      stop("Expecting data frame returned from format_data")
    }

    if(! "pval.exposure" %in% names(dat))
    {
      dat$pval.exposure <- 0.99
    }

    if(! "id.exposure" %in% names(dat))
    {
      dat$id.exposure <- random_string(1)
    }

    d <- data.frame(rsid=dat$SNP, pval=dat$pval.exposure, id=dat$id.exposure)
    out <- ieugwasr::ld_clump(d, clump_kb=clump_kb, clump_r2=clump_r2, clump_p=clump_p1)
    keep <- paste(dat$SNP, dat$id.exposure) %in% paste(out$rsid, out$id)
    return(dat[keep, ])
  }

  ld_matrix <- function(snps, with_alleles=TRUE)
  {
    .Deprecated("ieugwasr::ld_matrix()")
    ieugwasr::ld_matrix(variants=snps, with_alleles=with_alleles)
  }

  # setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  zScore <- function (cntA, totA, cntB, totB) {
    #calculate
    avgProportion <- (cntA + cntB) / (totA + totB)
    probA <- cntA/totA
    probB <- cntB/totB
    SE <- sqrt(avgProportion * (1-avgProportion)*(1/totA + 1/totB))
    zScore <- (probA-probB) / SE
    return (zScore)
  }

  ######################
  ### main functions ###
  ######################

  #--------------------------- read inputs

  #------- chromatin_states
  chromatin_states <- GRanges(sample = BackendData_ChromatinStates[[4]],
                             seqnames = gsub('chr','',BackendData_ChromatinStates[[1]]),
                             ranges = IRanges(BackendData_ChromatinStates[[2]], BackendData_ChromatinStates[[3]]),
                             state = BackendData_ChromatinStates[[4]])

  #------- BackendData_GenomicFeatures
  BackendData_GenomicFeatures$seqnames <- gsub('chr','',BackendData_GenomicFeatures$seqnames)
  BackendData_GenomicFeatures <- BackendData_GenomicFeatures[which(BackendData_GenomicFeatures$seqnames %in% c(1:22, 'X')),]
  BackendData_GenomicFeatures <- GRanges(sample = BackendData_GenomicFeatures$sample,
                                        seqnames = BackendData_GenomicFeatures$seqnames,
                                        ranges = IRanges(BackendData_GenomicFeatures$start, BackendData_GenomicFeatures$end),
                                        state = BackendData_GenomicFeatures$state)

  #--------------------------------- input
  #- query
  indexSNP <- SNP

  SNP_chr_POS <- snps.from.rsid(rsid = SNP, dbSNP = SNPlocs.Hsapiens.dbSNP144.GRCh38, search.genome = BSgenome.Hsapiens.UCSC.hg38)
  SNP_chr_POS <- data.frame(SNP_chr_POS, stringsAsFactors = FALSE)

  posQuerySNP <- GRanges(paste0(gsub('chr','',as.character(SNP_chr_POS$seqnames)),':',as.character(SNP_chr_POS$start),'-',as.character(SNP_chr_POS$end))) + windowSize
  #Using the 1MB window that we defined for the index variant, we capture the relevant variants from the VCF file for our LD calculations.
   
  if(!is.null(vcfPATH))
  { 
     
    #---- vcf is the address of the vcf file
    chr.remote.vcf <- vcfPATH
    vcfsubset <- NULL
    vcfsubset <- GetVariantsInWindow(file = chr.remote.vcf,position = posQuerySNP[1], type = "vcf")
    my.samples <- fread(vcfMetaData, stringsAsFactors = FALSE)
    vcfsubset <- SetPopulation(vcfsubset, sample_sheet = my.samples)
    row.names(vcfsubset) <- make.names(row.names(vcfsubset), unique = TRUE)
     
  }

  LD <- CalcLD(vcfsubset, index = indexSNP, population = "EUR")

  #----------------- foreground
  vcfsubsetsnps <- SplitVcfLd(vcf = LD, ld = c(metric = "R.squared", cutoff = 1, maf = mafThreshold), strict.subset = TRUE) #a strict subset cannot be the same set, that is, it cannot contain all of the elements that the other set does. Or in other words, a strict subset must be smaller, while a subset can be the same size.
  length( as.character(names(vcfsubsetsnps$fg)) )
  length( as.character(names(vcfsubsetsnps$bg)) )
  F1 <- as.character(names(vcfsubsetsnps$fg))

  fg_variants <- snps.from.rsid(rsid = F1, dbSNP = SNPlocs.Hsapiens.dbSNP144.GRCh38, search.genome = BSgenome.Hsapiens.UCSC.hg38)

  LD <- CalcLD(vcfsubset, index = indexSNP, population = "EUR")

  vcfsubsetsnps <- SplitVcfLd(vcf = LD, ld = c(metric = "R.squared", cutoff = 0.001, maf = mafThreshold), strict.subset = TRUE) #a strict subset cannot be the same set, that is, it cannot contain all of the elements that the other set does. Or in other words, a strict subset must be smaller, while a subset can be the same size.
  length( as.character(names(vcfsubsetsnps$fg)) )
  length( as.character(names(vcfsubsetsnps$bg)) )

  bg_variants <- snps.from.rsid(rsid = as.character(names(vcfsubsetsnps$bg)),
                                dbSNP = SNPlocs.Hsapiens.dbSNP144.GRCh38,
                                search.genome = BSgenome.Hsapiens.UCSC.hg38)

  length(names(vcfsubsetsnps$bg))

  foreground_sub <- as.data.frame(fg_variants)
  foreground_sub$seqnames <- gsub('chr', '', foreground_sub$seqnames)

  fGRhg38 <- GRanges(
    seqnames <- Rle(as.character(foreground_sub$seqnames)),
    ranges <- IRanges(start=as.numeric(foreground_sub$start), end=as.numeric(foreground_sub$end)+1, names=foreground_sub$SNP_ID),
    strand <- Rle(as.character(foreground_sub$strand))
  )

  background_sub <- as.data.frame(bg_variants)
  background_sub$seqnames <- gsub('chr', '', background_sub$seqnames)

  bGRhg38 <- GRanges(
    seqnames<-Rle(as.character(background_sub$seqnames)),
    ranges <- IRanges(start=as.numeric(background_sub$start), end=as.numeric(background_sub$end)+1, names=background_sub$SNP_ID),
    strand <- Rle(as.character(background_sub$strand))
  )

  #--------------------------------- CS enrichment

  #--- enrich foreground
  fo <- findOverlaps(query=fGRhg38, subject=chromatin_states, type="any")
  foreground_CS <- chromatin_states[subjectHits(fo),]

  #--- enrich background
  fo <- findOverlaps(query=bGRhg38, subject=chromatin_states, type="any")
  Background_CS <- chromatin_states[subjectHits(fo),]
  CSs <- unique(foreground_CS$sample)

  for(t in 1:length(CSs))
  {
    cntA<-length(which(foreground_CS$sample==CSs[t]))
    totA<-length(fGRhg38)
    cntB<-length(which(Background_CS$sample==CSs[t]))
    totB<-length(bGRhg38)

    table <- c(cntA, totA, cntB, totB)
    dim(table)<-c(2,2)
    fishert <- fisher.test(table)

    z_score <- zScore(cntA, totA, cntB, totB)

    temp <- data.frame(SNP = indexSNP, CS = CSs[t], p = fishert$p.value, oddsratio = fishert$estimate, CIl= fishert$conf.int[1], CIh = fishert$conf.int[2], zScore = z_score, FE = cntA, FA = totA, BE = cntB, BA = totB)

    if(t==1){RESULTsChromatinState = temp}else{RESULTsChromatinState = rbind(RESULTsChromatinState, temp)}

    print(z_score)
  }

  #------------------- genomic features

  #--- enrich foreground
  fo <- findOverlaps(query=fGRhg38, subject=BackendData_GenomicFeatures, type="any")
  foreground_GR <- BackendData_GenomicFeatures[subjectHits(fo),]

  #--- enrich background
  fo <- findOverlaps(query=bGRhg38, subject=BackendData_GenomicFeatures, type="any")
  Background_GR <- BackendData_GenomicFeatures[subjectHits(fo),]

  if(length(foreground_GR)>0)
  {
    GRs <- unique(foreground_GR$sample)
    for(t in 1:length(GRs))
    {
      cntA<-length(which(foreground_GR$sample==GRs[t]))
      totA<-length(foreground_GR)
      cntB<-length(which(Background_GR$sample==GRs[t]))
      totB<-length(Background_GR)

      table <- c(cntA, totA, cntB, totB)
      dim(table)<-c(2,2)
      fishert <- fisher.test(table)

      z_score <- zScore(cntA, totA, cntB, totB)

      temp <- data.frame(SNP = indexSNP, GR = GRs[t], p = fishert$p.value, oddsratio = fishert$estimate, CIl= fishert$conf.int[1], CIh = fishert$conf.int[2], zScore = z_score, FE = cntA, FA = totA, BE = cntB, BA = totB)

      if(t==1){RESULTsGenomicFeatures = temp}else{RESULTsGenomicFeatures = rbind(RESULTsGenomicFeatures, temp)}

      print(z_score)
    }
  }

  #------------------- save outputs

  Destiny_Folder <- system.file(package = "FEVV")
  Destiny_Folder <- paste(Destiny_Folder, "/RESULTsGenomicFeatures.txt", sep = "")

  write.table(
    RESULTsGenomicFeatures, Destiny_Folder, sep = "\t", row.names = FALSE, quote = FALSE
  )

  Destiny_Folder <- system.file(package = "FEVV")
  Destiny_Folder <- paste(Destiny_Folder, "/RESULTsChromatinState.txt", sep = "")

  write.table(
    RESULTsChromatinState, Destiny_Folder, sep = "\t", row.names = FALSE, quote = FALSE
  )

  print("You can find the results at: ")
  print(system.file(package="FEVV"))

}

