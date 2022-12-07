#!/usr/bin/env python3

import sys
import json
import os

config = json.load(open(sys.argv[1]))
outdir = "pollution"

if not os.path.exists(outdir):
    os.mkdir(outdir)

for sample, read1, read2 in config["reads_data"]:
    command = f"zcat {read1} |head -n 10000 |seqtk seq -A >{outdir}/{sample}.fa"
    # print(command)
    os.system(command)

'''
# shell: 
cd pollution

for i in $(ls *.fa); do
    blastn -query ${i} -db /mnt/sdb1/home/public_data/database/nt/nt -num_threads 30 -max_target_seqs 1 -max_hsps 1 -outfmt '6 std ssciname' >${i}.blast_nt.fmt6
done

for i in $(ls *.fmt6); do
    cut -f 13 ${i} |cut -d ' ' -f 1 |sort |uniq -c |awk '{print $2"\t"$1}' >${i}.stat
done

for i in $(ls *.fa |sed 's/.fa//'); do 
    awk -v a=${i} '{print a"\t"$1"\t"$2}' ${i}*.stat; 
done >merged.xls

# Rscript: 
data <- read.csv('merged.xls', sep='\t', header=F)
data2 <- data[data$V3>100, ]
colnames(data2) <- c("Sample", "Genus", "Read_count")
p <- ggplot(data2, aes(Sample, Read_count, fill=Genus)) + 
    geom_bar(stat="identity", position="stack") + 
    labs(title="各样品reads来源属分布") + 
    coord_flip()

ggsave("reads_source.png", p)
'''
