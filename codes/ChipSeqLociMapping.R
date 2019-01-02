###
#   File name : ChipSeqLociMapping.R
#   Author    : Hyunjin Kim
#   Date      : Feb 5, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Calculate how many chip-seq reads fall into the loci 
#
#   Instruction
#               1. Source("ChipSeqLociMapping.R")
#               2. Run the function "chipLociMapping()" - specify input BED file, loci info, and output directory
#               3. A mapping info will be generated in the output directory
#
#   Example
#               > source("The_directory_of_ChipSeqLociMapping.R/ChipSeqLociMapping.R")
#               > chipLociMapping(bedFilePath="./data/GSE39908/7_P_BD1_F3-F5-BC-Paired2-W200-G600-islands-summary-FDR001.bed", lociPath="./data/214_precursors.xlsx", outputDir="./results/loci_map_chip_seq/")
###

chipLociMapping <- function(bedFilePath="./data/GSE39908/7_P_BD1_F3-F5-BC-Paired2-W200-G600-islands-summary-FDR001.bed", lociPath="./data/214_precursors.xlsx", outputDir="./results/loci_map_chip_seq/") {
  
  ### load library
  if(!require(xlsx)) {
    install.packages("xlsx")
    library(xlsx)
  }
  
  
  ### load datasets
  bed <- as.data.frame(read.table(bedFilePath, sep="\t", stringsAsFactors=FALSE, quote=""))
  loci <- read.xlsx2(lociPath, sheetIndex=1, stringsAsFactors=FALSE)
  
  
  ### loci mapping
  map <- cbind(loci, Map=0, MapPcnt=0)
  
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
  
  if(!is.null(outputDir)) {
    fileName <- strsplit(bedFilePath, "/", fixed = TRUE)
    fileName <- fileName[[1]][length(fileName[[1]])]
    fileName <- paste0(substr(fileName, 1, nchar(fileName)-4), ".csv")
    
    write.csv(map, paste0(outputDir, fileName))
  }
  
  return(map)
}
