###
#   File name : OverlapTwoBedGraph.R
#   Author    : Hyunjin Kim
#   Date      : Feb 17, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Load two bed graph files, calculate peak overlap, and draw venn diagram with it
#
#   Instruction
#               1. Source("OverlapTwoBedGraph.R")
#               2. Run the function "overlap_venn()" - specify two input BED GRAPH file paths and output file path
#               3. A venn diagram of the overlap info will be generated in the output file path
#
#   Example
#               > source("The_directory_of_OverlapTwoBedGraph.R/OverlapTwoBedGraph.R")
#               > overlap_venn(file1Path="./data/Chip_Seq/Pachytene_wt_Brdt_GA3305_norm10m-with-header.graph", file2Path="./data/Chip_Seq/Pachytene_wt_IgG_GA3306_norm10m-with-header.graph", outputPath="./results/GA3305_3306_overlap.pdf")
###

overlap_venn <- function(file1Path="./data/Chip_Seq/Pachytene_wt_Brdt_GA3305_norm10m-with-header.graph", file2Path="./data/Chip_Seq/Pachytene_wt_IgG_GA3306_norm10m-with-header.graph", outputPath="./results/GA3305_3306_overlap.pdf") {
  
  ### load library
  if(!require(VennDiagram)) {
    install.packages("VennDiagram")
    library(VennDiagram)
  }
  if(!require(gridExtra)) {
    install.packages("gridExtra")
    library(gridExtra)
  }
  

  ### load datasets
  bedgraph1 <- as.data.frame(read.table(file1Path, sep="\t", skip=1, check.names=FALSE, stringsAsFactors=FALSE, quote=""))
  bedgraph2 <- as.data.frame(read.table(file2Path, sep="\t", skip=1, check.names=FALSE, stringsAsFactors=FALSE, quote=""))
  
  ### make lists with chr, start, end
  list1 <- paste0(bedgraph1$V1, bedgraph1$V2, bedgraph1$V3)
  list2 <- paste0(bedgraph2$V1, bedgraph2$V2, bedgraph2$V3)
  
  ### draw a venn diagram
  v <- venn.diagram(list(list1, list2), category.names = c("GA3305", "GA3306"), cat.cex = 1.5, cex = 1.5, filename = NULL)
  
  ### save the diaram as pdf
  pdf(outputPath)
  grid.arrange(gTree(children=v), top="Chip-Seq: GA3305 vs GA3306", bottom="")
  dev.off()
  
}
