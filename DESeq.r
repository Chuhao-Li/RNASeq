library("DESeq2")
library('sva')
# 绘图相关
library("ggplot2")
library("reshape2")
library("patchwork")
library("ggdendro")

my_data <- read.csv('featureCounts/count_table.xls', sep = '\t', row.names = 1)
coldata <- read.csv("DESeq/coldata.xls", row.names=1)

# 批次效应校正
if ("batch" %in% colnames(coldata) & length(unique(coldata$batch)) > 1){
    batch <- coldata$batch
    group <- as.numeric(as.factor(coldata$condition))
    my_data <- ComBat_seq(as.matrix(my_data), batch=batch, group=group)
}

# normalization
dds <- DESeqDataSetFromMatrix(countData = my_data,
                              colData = coldata, 
                              design = ~ condition)

dds <- DESeq(dds)
normalized <- counts(dds, normalized=TRUE) # 提取标准化后的read count
write.csv(normalized, 'DESeq/normalized_readCount.csv', quote=F)

#############################################
# 样品重复性相关图表：PCA、层次聚类、箱线图 #
#############################################

normalized.long <- melt(normalized)
vsd <- vst(dds)
p1 <- plotPCA(vsd, intgroup=c("condition")) + labs(title = "A")# PCA

# 层次聚类
normalized.diff <- normalized[colSums(normalized) != 0, ]
normalized.dist <- dist(t(normalized.diff), method = "euclidean")
normalized.hc <- hclust(d = normalized.dist, method = "ward.D2")
p1.1 <- ggdendrogram(normalized.hc) + labs(title = "B")
# 箱线图
p2 <- ggplot(normalized.long, aes(Var2, value)) + geom_boxplot() + labs(title = "C")
p3 <- ggplot(normalized.long, aes(Var2, value)) + geom_boxplot() + ylim(c(0, 5000)) + labs(title = "D")

##########################
# 组间比较。需要选择CK。 #
##########################
pvalue.cutoff <- 0.05
log2fc.cutoff <- 2
baseMean.cutoff <- 10

# 可能有多组。现在只考虑有1组的情况。
args = commandArgs(trailingOnly=TRUE)
contrast1 <- args[1]
contrast2 <- strsplit(contrast1, "_vs_")[[1]]
res <- results(dds, contrast=c("condition", contrast2[1], contrast2[2])) # contrast=c("condition", "", "")。野生型/CK应该在右边。
resOrdered <- as.data.frame(res[order(res$pvalue),])
anno <- read.csv('reference.idmap', sep='\t')
anno <- anno[, c('locusTag', 'product')]
resOrdered_anno <- merge(resOrdered, anno, by.x='row.names', by.y="locusTag", sort=F)
resOrdered_anno <- as.data.frame(resOrdered_anno)
resOrdered_anno.up <- resOrdered_anno[resOrdered_anno$baseMean > baseMean.cutoff & resOrdered_anno$log2FoldChange > log2fc.cutoff & resOrdered_anno$pvalue < pvalue.cutoff, ]
resOrdered_anno.down <- resOrdered_anno[resOrdered_anno$baseMean > baseMean.cutoff & resOrdered_anno$log2FoldChange < -log2fc.cutoff & resOrdered_anno$pvalue < pvalue.cutoff, ]

# 使用A_vs_B来对结果进行命名。
write.csv(resOrdered_anno, paste0('DESeq/', contrast1, '.csv'), quote=F, row.names=F) 
write.csv(resOrdered_anno.up, paste0('DESeq/', contrast1, '.up.csv'), quote=F, row.names=F)
write.csv(resOrdered_anno.down, paste0('DESeq/', contrast1, '.down.csv'), quote=F, row.names=F)

# MA plot
p4 <- ggplot(as.data.frame(res), aes(log2FoldChange, baseMean)) + 
    geom_point(aes(fill=ifelse(baseMean>10 & abs(log2FoldChange)>2 & pvalue <0.05, "red", "blue")), shape = 21, color="black") + 
    geom_vline(xintercept=2, linetype = "longdash") + 
    geom_vline(xintercept=-2, linetype = "longdash") + 
    scale_fill_identity() + 
    labs(title = "E")

p5 <- ggplot(as.data.frame(res), aes(log2FoldChange, baseMean)) + 
    geom_point(aes(fill=ifelse(baseMean>10 & abs(log2FoldChange)>2 & pvalue <0.05, "red", "blue")), shape = 21, color="black") + 
    ylim(c(0, 10000)) + 
    geom_vline(xintercept=2, linetype = "longdash") + 
    geom_vline(xintercept=-2, linetype = "longdash") + 
    scale_fill_identity() + labs(title = "F")

# 拼图
pall <- (p1+p1.1)/(p2+p3)/(p4+p5)
ggsave("DESeq/pall.png", pall)
