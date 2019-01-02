###
#   File name : RNASeqLociMappingTSS.R
#   Author    : Hyunjin Kim
#   Date      : Feb 13, 2018
#   Email     : hk2990@cumc.columbia.edu
#   Purpose   : Calculate how many rna-seq reads fall into TSS region of the loci 
#
#   Instruction
#               1. Source("RNASeqLociMappingTSS.R")
#               2. Run the function "rnaLociMappingTSS()" - specify input BAM file, loci info, and output directory
#               3. A mapping info will be generated in the output directory
#
#   Example
#               > source("The_directory_of_RNASeqLociMappingTSS.R/RNASeqLociMappingTSS.R")
#               > rnaLociMappingTSS(bamFilePath="./data/RNA_Seq/DJW4_mouse.bam", lociPath="./data/214_precursors.xlsx", outputDir="./results/loci_map_rna_seq/")
###

rnaLociMappingTSS <- function(bamFilePath="./data/RNA_Seq/DJW4_mouse.bam", lociPath="./data/214_precursors.xlsx", outputDir="./results/loci_map_rna_seq/") {
  
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
  if(!require(BSgenome.Mmusculus.UCSC.mm9)) {
    source("http://bioconductor.org/biocLite.R") 
    biocLite("BSgenome.Mmusculus.UCSC.mm9") 
    library(BSgenome.Mmusculus.UCSC.mm9)
  }
  
  
  ### load datasets
  loci <- read.xlsx2(lociPath, sheetIndex=1, stringsAsFactors=FALSE)
  
  
  ### get TSS region of the loci
  genome <- BSgenome.Mmusculus.UCSC.mm9
  gr <- GRanges(Rle(loci$Chr), IRanges(as.numeric(loci$Start), as.numeric(loci$Stop)))
  
  loci_seq <- getSeq(genome, gr)
  
  vmatchPattern("ATG", loci_seq[1])
  
  
  ### loci mapping
  map <- cbind(loci, MapCnt=0)
  
  bf <- BamFile(bamFilePath)
  open(bf)
  chunk <- readGAlignments(bf)
  
  
  
}

