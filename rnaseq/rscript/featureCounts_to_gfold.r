cal_rpkm <- function(counts, gene_lens){
  r_size = sum(counts) # total read number
  rpm = counts*1e06/r_size
  rpkm = rpm*1e03/gene_lens
  return(rpkm)
}

count_table <- read.csv('featureCounts/featureCounts.txt', sep='\t', comment.char="#", check.names=F)
my.columns <- colnames(count_table) # read featureCounts/featureCounts.txt to make sure the order of samples are the sames. 
my.columns <- gsub("bowtie2/(.*)\\.bam", "\\1", my.columns, perl=T) 
print(my.columns)
colnames(count_table) <- my.columns
GeneSymbol <- count_table$Geneid
GeneName <- count_table$Geneid
for(i in colnames(count_table)[7:dim(count_table)[2]]){
  ReadCount <- count_table[, i]
  GeneLen <- count_table$Length
  RPKM <- cal_rpkm(ReadCount, GeneLen)
  # GeneSymbol, GeneName, ReadCount, GeneLen, RPKM
  df <- data.frame(GeneSymbol, GeneName, ReadCount, GeneLen, RPKM)
  write.table(df, paste0('gfold/', i, '.count'), sep='\t', row.names=F, col.names=F, quote=F)
}

