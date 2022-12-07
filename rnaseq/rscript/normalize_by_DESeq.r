library("DESeq2")

my_data <- read.csv('featureCounts/count_table.xls', sep = '\t', row.names = 1)
coldata.tmp <- read.csv("gfold/coldata.xls", row.names=1)
coldata <- coldata.tmp[colnames(my_data), , drop=F]

if(nrow(coldata.tmp) != nrow(coldata)){
    print("rownames of coldata should be the same as colnames of count table! ")
    quit()
}

# normalization
dds <- DESeqDataSetFromMatrix(countData = my_data,
                              colData = coldata,
                              design = ~ 1)

dds <- DESeq(dds)
normalized <- counts(dds, normalized=TRUE) # 提取标准化后的read count
write.csv(normalized, 'gfold/normalized_readCount.csv', quote=F)

