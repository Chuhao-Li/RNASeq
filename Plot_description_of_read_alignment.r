library('ggplot2')
library('reshape2')
data <- read.csv("featureCounts.txt.summary", sep='\t', check.names=F)
colnames(data) <- sub("bowtie2/(.*)\\.bam", "\\1", colnames(data))
data <- data[rowSums(data[-1]) != 0, ]

data <- melt(data, id.vars = "Status", variable.name = "sample", value.name = "reads_count")
data.sum <- aggregate(reads_count~sample, sum, data=data)
p <- ggplot(data) + 
  geom_bar(aes(sample, reads_count, fill=Status), stat="identity", position="stack") + 
  geom_label(data=data.sum, aes(as.factor(sample), reads_count, label=format(reads_count, big.mark=","), vjust=0, )) + 
  ggtitle("Description of read alignment")
  
ggsave("Description_of_read_alignment.png", p, width=6, height=5)
