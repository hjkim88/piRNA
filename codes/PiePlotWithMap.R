###
#   File name : PiePlotWithMap.R
#   Author    : Hyunjin Kim
#   Date      : Feb 6, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Plot pie charts with multiple mapping info
#
#   Instruction
#               1. Source("PiePlotWithMap.R")
#               2. Run the function "piePlot()" - specify input BED file directory and output directory
#               3. A plot with many pie charts of the mapping info will be generated in the output directory
#
#   Example
#               > source("The_directory_of_PiePlotWithMap.R/PiePlotWithMap.R")
#               > piePlot(bedPath="./data/GSE39908/", lociPath="./data/214_precursors.xlsx", mappingCodePath="./codes/ChipSeqLociMapping.R", outputPath="./results/loci_map_chip_seq/GSE39908_PIE.pdf")
###

piePlot <- function(bedPath="./data/GSE39908/", lociPath="./data/214_precursors.xlsx", mappingCodePath="./codes/ChipSeqLociMapping.R", outputPath="./results/loci_map_chip_seq/GSE39908_PIE.pdf") {
  
  ### load codes
  source(mappingCodePath)
  
  
  ### collect files from the bedPath
  f <- list.files(bedPath)
  
  
  ### decide the number of rows in the pie plot
  plotRowNum <- round((length(f)+1)/2)
  
  
  ### save the results in pdf file
  pdf(outputPath)
  par(mfrow=c(plotRowNum, 2), mar=c(3.1, 5.1, 3.1, 5.1), oma=c(0,0,2,0))
  
  for(i in 1:length(f)) {
    map <- chipLociMapping(bedFilePath = paste0(bedPath, f[i]), lociPath = lociPath, outputDir = NULL)
    
    slices <- c(sum(map$Map), nrow(map)-sum(map$Map))
    pct <- round(slices/sum(slices)*100)
    lbls <- c("Mapped", "Unmapped")
    lbls <- paste0(lbls, " ", slices, "/", sum(slices), " (", pct, "%)")
    
    fileName <- paste0(substr(f[i], 1, nchar(f[i])-4))
    
    pie(slices, labels=lbls, col=rainbow(length(lbls)), main=fileName, cex.main=1, cex=0.7)
  }
  
  fileName <- strsplit(outputPath, "/", fixed = TRUE)
  fileName <- fileName[[1]][length(fileName[[1]])]
  fileName <- paste0(substr(fileName, 1, nchar(fileName)-4))
  title(fileName, outer = TRUE, cex = 1.5)
  dev.off()
  
}


