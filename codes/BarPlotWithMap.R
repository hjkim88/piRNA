###
#   File name : BarPlotWithMap.R
#   Author    : Hyunjin Kim
#   Date      : Feb 8, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Plot bar charts with multiple mapping info
#
#   Instruction
#               1. Source("BarPlotWithMap.R")
#               2. Run the function "barPlot()" - specify input mapping file directory and output directory
#               3. A plot with many bar charts of the mapping info will be generated in the output directory
#
#   Example
#               > source("The_directory_of_BarPlotWithMap.R/BarPlotWithMap.R")
#               > barPlot(mapFilePath="./results/loci_map_rna_seq/", outputPath="./results/loci_map_rna_seq/RNA-Seq_Map_BAR.pdf", plotNum=3)
###

barPlot <- function(mapFilePath="./results/loci_map_rna_seq/", outputPath="./results/loci_map_rna_seq/RNA-Seq_Map_BAR.pdf", plotNum=3) {
  
  ### load library
  if(!require(ggplot2)) {
    install.packages("ggplot2")
    library(ggplot2)
  }
  if(!require(reshape)) {
    install.packages("reshape")
    library(reshape)
  }
  
  source("./codes/multiplot.R")
  
  ### collect map files
  f <- list.files(mapFilePath)
  
  temp <- read.csv(paste0(mapFilePath, f[1]), row.names=1, check.names = FALSE, stringsAsFactors = FALSE)
  df <- data.frame(temp$MapCnt)
  rownames(df) <- temp$Name
  colnames(df)[1] <- strsplit(as.character(f), ".", fixed=TRUE)[[1]][1]
  
  if(length(f) > 1) {
    for(i in 2:length(f)) {
      temp <- read.csv(paste0(mapFilePath, f[i]), row.names=1, check.names = FALSE, stringsAsFactors = FALSE)
      df <- cbind(df, temp$MapCnt)
      colnames(df)[i] <- strsplit(as.character(f), ".", fixed=TRUE)[[i]][1]
    }
  }
  
  df <- cbind(df, Precursor=rownames(df))
  df <- df[order(-df[,1]),]
  df$Precursor <- factor(df$Precursor, levels = rownames(df))
  
  mdf <- list()
  if(plotNum > 1) {
    sp <- round(nrow(df)/plotNum)
    for(i in 1:(plotNum-1)) {
      mdf[[i]] <- df[(sp*(i-1)+1):(sp*i),]
    }
    mdf[[plotNum]] <- df[(sp*(plotNum-1)):nrow(df),]
  } else if(plotNum == 1) {
    mdf[[plotNum]] <- df
  } else {
    writeLines(paste("Wrong plotNum:", plotNum))
  }
  
  p <- list()
  for(i in 1:plotNum) {
    data.m <- melt(mdf[[i]], id.vars="Precursor")
    p[[i]] <- ggplot(data.m, aes(x=Precursor, y=value)) + labs(x="", y="Mapped Reads") + theme_classic(base_size = 16) + theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) + geom_bar(aes(fill = variable), position = "dodge", stat = "identity") + guides(fill=guide_legend(title=NULL))
  }
  
  fileName <- strsplit(outputPath, "/", fixed = TRUE)
  fileName <- fileName[[1]][length(fileName[[1]])]
  fileName <- paste0(substr(fileName, 1, nchar(fileName)-4))
  
  pdf(outputPath, width = 20, height = 10)
  multiplot(plotlist = p, cols = 1, title = fileName)
  dev.off()
  
}

