###
#   File name : RNASeqLociMapping.R
#   Author    : Hyunjin Kim
#   Date      : Feb 7, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Calculate how many rna-seq reads fall into the loci 
#
#   Instruction
#               1. Source("RNASeqLociMapping.R")
#               2. Run the function "rnaLociMapping()" - specify input BAM file, loci info, and output directory
#               3. A mapping info will be generated in the output directory
#
#   Example
#               > source("The_directory_of_RNASeqLociMapping.R/RNASeqLociMapping.R")
#               > rnaLociMapping(bamFilePath="./data/RNA_Seq/DJW4_mouse.bam", lociPath="./data/214_precursors.xlsx", outputDir="./results/loci_map_rna_seq/", overLapTh=1)
###

rnaLociMapping <- function(bamFilePath="./data/RNA_Seq/DJW4_mouse.bam", lociPath="./data/214_precursors.xlsx", outputDir="./results/loci_map_rna_seq/", overLapTh=1) {
  
  ### load library
  if(!require(xlsx)) {
    install.packages("xlsx")
    library(xlsx)
  }
  if(!require(GenomicAlignments)) {
    source("http://bioconductor.org/biocLite.R") 
    biocLite("GenomicAlignments") 
    library(GenomicAlignments)
  }
  if(!require(rtracklayer)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("rtracklayer")
    library(rtracklayer)
  }
  
  
  ### load datasets
  loci <- read.xlsx2(lociPath, sheetIndex=1, stringsAsFactors=FALSE)
  
  
  ### it turned out that the bam files are mm10 while the loci are mm9
  ###############################################################################
  ### change loci to granges form
  gr <- GRanges(Rle(loci$Chr), IRanges(as.numeric(loci$Start), as.numeric(loci$Stop)))
  
  
  ### url and file name for a chain file
  url <- "http://hgdownload.soe.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz"
  mm9ToMm10.chain <- "./data/mm9ToMm10over.chain"
  
  if (!file.exists(mm9ToMm10.chain)) {
    download.file(url, destfile = paste0(mm9ToMm10.chain, ".gz"))
    gunzip(paste0(mm9ToMm10.chain, ".gz"))
  }
  
  ### import chain file
  chain <- import.chain(mm9ToMm10.chain)
  
  ### change loci to mm10 form
  gr_mm10 <- liftOver(gr, chain)
  
  loci_mm10 <- data.frame(matrix(0, length(gr_mm10), 5))
  colnames(loci_mm10) <- c("Name", "Chr", "Start", "Stop", "Strand")
  loci_mm10$Name <- loci$Name
  loci_mm10$Chr <- loci$Chr
  loci_mm10$Strand <- loci$Strand
  
  ### A function to get lifted ranges
  getLiftedRanges <- function(chr_info, grangeObj) {
    nr <- c(0, 0)
    
    temp <- grangeObj[which(as.character(seqnames(grangeObj)) == chr_info)]
    
    nr[1] <- start(temp[1])
    nr[2] <- end(temp[length(temp)])
    
    return(nr)
  }
  
  ### get lifted loci in mm10
  for(i in 1:nrow(loci_mm10)) {
    temp <- getLiftedRanges(loci_mm10$Chr[i], gr_mm10[[i]])
    
    loci_mm10$Start[i] <- temp[1]
    loci_mm10$Stop[i] <- temp[2]
  }
  
  loci <- loci_mm10
  ###############################################################################
  
  
  ### loci mapping
  map <- cbind(loci, MapCnt=0, pMapCnt=0, mMapCnt=0)
  
  bf <- BamFile(bamFilePath)
  open(bf)
  chunk <- readGAlignments(bf)
  
  
  findStartIdx <- function(contList, fixedP) {
    startIdx <- 1
    sp <- round(length(contList)/2)

    repeat {
      if((contList[sp] < fixedP) && (contList[sp+1] >= fixedP)) {
        startIdx <- as.numeric(names(contList)[sp])
        break
      } else if(contList[sp] >= fixedP) {
        contList <- contList[1:sp]
        sp <- round(length(contList)/2)
      } else {
        contList <- contList[sp:length(contList)]
        sp <- round(length(contList)/2)
      }
    }

    return(startIdx)
  }
  
  
  getSharedLen <- function(r1_start, r1_end, r2_start, r2_end) {
    sharedLen <- 0
    
    if((r1_end <= r2_end) && (r1_start <= r2_start)) {
      sharedLen <- r1_end - r2_start + 1
    } else if((r1_end <= r2_end) && (r1_start > r2_start)) {
      sharedLen <- r1_end - r1_start + 1
    } else if(r1_start >= r2_start){
      sharedLen <- r2_end - r1_start + 1  
    } else {
      sharedLen <- r2_end - r2_start + 1
    }
    
    return(sharedLen)
  }
  
  
  for(i in 1:nrow(map)) {
    idx <- which(as.vector(seqnames(chunk)) == as.character(map$Chr[i]))
    
    starts <- start(chunk)[idx]
    ends <- end(chunk)[idx]
    strands <- as.character(strand(chunk)[idx])
    
    names(ends) <- 1:length(ends)
    
    startIdx <- 1
    startIdx <- findStartIdx(ends, as.numeric(map$Start[i]))
    
    for(j in startIdx:length(starts)) {
      if((ends[j] >= as.numeric(map$Start[i])) && (starts[j] <= as.numeric(map$Stop[i]))) {
        if(getSharedLen(starts[j], ends[j], as.numeric(map$Start[i]), as.numeric(map$Stop[i])) >= overLapTh) {
          map$MapCnt[i] <- map$MapCnt[i]+1
          
          if(strands[j] == "+") {
            map$pMapCnt[i] <- map$pMapCnt[i]+1 
          } else if(strands[j] == "-") {
            map$mMapCnt[i] <- map$mMapCnt[i]+1
          }
        }
      } else if(starts[j] > as.numeric(map$Stop[i])) {
        break
      }
    }
    
    writeLines(paste(i, "/", nrow(map)))
  }
  
  close(bf)
  
  if(!is.null(outputDir)) {
    fileName <- strsplit(bamFilePath, "/", fixed = TRUE)
    fileName <- fileName[[1]][length(fileName[[1]])]
    fileName <- paste0(substr(fileName, 1, nchar(fileName)-4), "_", overLapTh, "_mm10.csv")
    
    write.csv(map, paste0(outputDir, fileName))
  }
  
  return(map)
}
