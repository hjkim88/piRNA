###
#   File name : PiePlotWithBedGraph.R
#   Author    : Hyunjin Kim
#   Date      : Feb 17, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Load two bed graph files and loci info, remove noise (background), map bedgraph-loci, and make a pie plot
#
#   Instruction
#               1. Source("PiePlotWithBedGraph.R")
#               2. Run the function "piePlotBG()" - specify two input BED GRAPH file paths, loci path, and output file path
#               3. A pie plot of mapping info will be generated in the output file path
#
#   Example
#               > source("The_directory_of_PiePlotWithBedGraph.R/PiePlotWithBedGraph.R")
#               > piePlotBG(file1Path="./data/Chip_Seq/Pachytene_wt_Brdt_GA3305_norm10m-with-header.graph", file2Path="./data/Chip_Seq/Pachytene_wt_IgG_GA3306_norm10m-with-header.graph", lociPath="./data/214_precursors.xlsx", outputPath="./results/loci_map_chip_seq/GA3305_3306_PIE.pdf")
###


piePlotBG <- function(file1Path="./data/Chip_Seq/Pachytene_wt_Brdt_GA3305_norm10m-with-header.graph", file2Path="./data/Chip_Seq/Pachytene_wt_IgG_GA3306_norm10m-with-header.graph", lociPath="./data/214_precursors.xlsx", outputPath="./results/loci_map_chip_seq/GA3305_3306_PIE.pdf") {
  
  ### load library
  if(!require(xlsx)) {
    install.packages("xlsx")
    library(xlsx)
  }
  if(!require(rtracklayer)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("rtracklayer")
    library(rtracklayer)
  }
  if(!require(R.utils)) {
    install.packages("R.utils")
    library(R.utils)
  }
  
  
  ### load datasets
  ### bedgraph1 <- target, bedgraph2 <- background
  bedgraph1 <- as.data.frame(read.table(file1Path, sep="\t", skip=1, check.names=FALSE, stringsAsFactors=FALSE, quote=""))
  bedgraph2 <- as.data.frame(read.table(file2Path, sep="\t", skip=1, check.names=FALSE, stringsAsFactors=FALSE, quote=""))
  loci <- read.xlsx2(lociPath, sheetIndex=1, stringsAsFactors=FALSE)
  
  ### change loci to granges form
  gr <- GRanges(Rle(loci$Chr), IRanges(as.numeric(loci$Start), as.numeric(loci$Stop)))
  
  
  ### url and file name for a chain file
  url <- "http://hgdownload.soe.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm8.over.chain.gz"
  mm9ToMm8.chain <- "./data/mm9ToMm8over.chain"
  
  if (!file.exists(mm9ToMm8.chain)) {
    download.file(url, destfile = paste0(mm9ToMm8.chain, ".gz"))
    gunzip(paste0(mm9ToMm8.chain, ".gz"))
  }
  
  ### import chain file
  chain <- import.chain(mm9ToMm8.chain)
  
  ### change loci to mm8 form
  gr_mm8 <- liftOver(gr, chain)
  
  loci_mm8 <- data.frame(matrix(0, length(gr_mm8), 4))
  colnames(loci_mm8) <- c("Name", "Chr", "Start", "Stop")
  loci_mm8$Name <- loci$Name
  loci_mm8$Chr <- loci$Chr
  
  
  ### A function to get lifted ranges
  getLiftedRanges <- function(chr_info, grangeObj) {
    nr <- c(0, 0)
    
    temp <- grangeObj[which(as.character(seqnames(grangeObj)) == chr_info)]
    
    nr[1] <- start(temp[1])
    nr[2] <- end(temp[length(temp)])
    
    return(nr)
  }
  
  ### get lifted loci in mm8
  for(i in 1:nrow(loci_mm8)) {
    temp <- getLiftedRanges(loci_mm8$Chr[i], gr_mm8[[i]])
    
    loci_mm8$Start[i] <- temp[1]
    loci_mm8$Stop[i] <- temp[2]
  }
  
  
  ### bed graph1 - bed graph2: noise (background) removal
  list1 <- paste0(bedgraph1$V1, bedgraph1$V2, bedgraph1$V3)
  list2 <- paste0(bedgraph2$V1, bedgraph2$V2, bedgraph2$V3)
  
  targetIdx1 <- which(list1 %in% intersect(list1, list2))
  targetIdx2 <- which(list2 %in% intersect(list1, list2))
  
  bed <- bedgraph1
  bed$V4[targetIdx1] <- bedgraph1$V4[targetIdx1] - bedgraph2$V4[targetIdx2]
  bed <- bed[-which(bed$V4 < 0), ]
  
  
  ### loci mapping
  map <- cbind(loci_mm8, Map=0, MapPcnt=0)
  
  for(i in 1:nrow(map)) {
    temp <- bed[which(bed$V1 == map$Chr[i]),]
    
    startIdx <- 0
    endIdx <- 0
    
    for(j in 1:nrow(temp)) {
      if((temp$V2[j] <= as.numeric(map$Stop[i])) && (temp$V3[j] >= as.numeric(map$Start[i]))) {
        startIdx <- j
        break
      }
    }
    
    if(startIdx != 0) {
      map$Map[i] <- 1
      
      for(j in startIdx:nrow(temp)) {
        if((temp$V2[j] <= as.numeric(map$Stop[i])) && (temp$V3[j] >= as.numeric(map$Stop[i]))) {
          endIdx <- j
          break
        } else if((temp$V2[j] > as.numeric(map$Stop[i])) && (temp$V3[j] > as.numeric(map$Stop[i]))) {
          endIdx <- j-1
          break
        }
      }
      
      accum <- 0
      for(j in startIdx:endIdx) {
        if((temp$V3[j] <= as.numeric(map$Stop[i])) && (temp$V2[j] <= as.numeric(map$Start[i]))) {
          accum <- accum + (temp$V3[j] - as.numeric(map$Start[i]) + 1)
        } else if((temp$V3[j] <= as.numeric(map$Stop[i])) && (temp$V2[j] > as.numeric(map$Start[i]))) {
          accum <- accum + (temp$V3[j] - temp$V2[j] + 1)
        } else if(temp$V2[j] >= as.numeric(map$Start[i])){
          accum <- accum + (as.numeric(map$Stop[i]) - temp$V2[j] + 1)  
        } else {
          accum <- accum + (as.numeric(map$Stop[i]) - as.numeric(map$Start[i]) + 1)
        }
      }
      
      map$MapPcnt[i] <- accum / (as.numeric(map$Stop[i]) - as.numeric(map$Start[i]) + 1)
    }
  }
  
  
  ### pie plot
  slices <- c(sum(map$Map), nrow(map)-sum(map$Map))
  pct <- round(slices/sum(slices)*100)
  lbls <- c("Mapped", "Unmapped")
  lbls <- paste0(lbls, " ", slices, "/", sum(slices), " (", pct, "%)")
  
  fileName <- strsplit(outputPath, "/", fixed = TRUE)
  fileName <- fileName[[1]][length(fileName[[1]])]
  fileName <- paste0(substr(fileName, 1, nchar(fileName)-4))
  
  ### save plot
  pdf(outputPath)
  pie(slices, labels=lbls, col=rainbow(length(lbls)), main=fileName, cex.main=1, cex=0.7)
  dev.off()
  
}

