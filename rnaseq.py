#!/usr/bin/env python3

# todo :
# SNP calling
# 污染情况调查
# DESeq

# plot：火山图、PCA图
# batch effect
# KEGG、GO注释与富集分析
# report
# logging

import re
import os
import sys
import json
import argparse
from multiprocessing import Pool

script_base = os.path.dirname(os.path.abspath(__file__))

def check_ncbi_gff():
    '''检查是不是ncbi下载的gff文件。如果不是，则警告'''
    pass

def check_sequence_header():
    '''
    检查基因组序列的header与gff文件中的一致。
    '''
    pass

def check_dependency():
    '''
    检查依赖是否满足
    fastp, 
    bowtie2
    featureCounts
    gfold
    Rscript
    '''
    pass

def fastp(s):
    if os.path.exists("fastp/{s}_1.fq.gz".format(s=s)):
        return
    os.system("if [ ! -d fastp ]; then mkdir fastp; fi")
    read1, read2 = reads_data[s]
    command = '''fastp -i {read1} \
        -I {read2} \
        -o fastp/{s}_1.fq.gz \
        -O fastp/{s}_2.fq.gz \
        -j fastp/{s}.json \
        -h fastp/{s}.html
    '''.format(read1=read1, read2=read2, s=s)
    os.system(command)
    os.system("touch fastp/success")

def align(s):
    if os.path.exists("bowtie2/{s}.bam".format(s=s)):
        return
    os.system("if [ ! -d bowtie2 ]; then mkdir bowtie2; fi")
    index = "bowtie2/genome_index"
    os.system("bowtie2-build {fna} {index}".format(fna=fna, index=index))
    command = '''bowtie2 -x {index} \
        -1 fastp/{s}_1.fq.gz \
        -2 fastp/{s}_2.fq.gz \
        -p 16 |samtools sort -o bowtie2/{s}.bam
    '''.format(index=index, s=s)
    os.system(command)

def featureCount(samples):
    if os.path.exists("featureCounts/featureCounts.txt"):
        return
    os.system("if [ ! -d featureCounts ]; then mkdir featureCounts; fi")
    command = '''featureCounts \
      -a {gff} \
      -o featureCounts/featureCounts.txt \
      -s 2 -p -t CDS -g locus_tag -T 15 \
      {bams}
    '''.format(gff=gff, bams = " ".join(["bowtie2/{}.bam".format(i) for i in samples]))
    os.system(command)
    os.system("touch featureCounts/success")

######################
# Load configuration # 
######################
parser = argparse.ArgumentParser()
parser.description = "Bacterial RNASeq analysis pipeline"
parser.add_argument("-c", "--config", help="config file")
parser.add_argument("--no_dea", action="store_true", help="do not run differential expression analysis. ")
parser.add_argument("--no_gfold", action="store_true", help="do not run differential expression analysis by Gfold. ")
parser.add_argument("--no_deseq", action="store_true", help="do not run differential expression analysis by DESeq2. ")
parser.add_argument("-v", "--version", action="store_true", help="print pipeline version")
args = parser.parse_args()

if args.version:
    print("rnaseq v0.2.1")
    exit()

config_file = args.config
config = json.load(open(config_file))
samples = []
reads_data = {}
for i in config["reads_data"]:
    reads_data[i[0]] = (os.path.abspath(i[1]), os.path.abspath(i[2]))
    samples.append(i[0])

gff = os.path.abspath(config["genome_annotation"])
fna = os.path.abspath(config["genome_sequence"])

n_samples = len(samples)

workdir = os.path.abspath(config["output_dir"])
if not os.path.exists(workdir):
    os.mkdir(workdir)

os.chdir(workdir)
##################################
# Clean data, aligning and count #
##################################

print("running fastp...")
with Pool(processes=n_samples) as p:
    p.map(fastp, samples)

print("alinging...")
for s in samples:
    align(s)

print("featureCounts...")
featureCount(samples)

####################################
# Differential expression analysis #
####################################
if args.no_dea:
    exit()

os.system("python {}/gff2idmap.py {} >reference.idmap".format(script_base, gff))

if "design" in config and not args.no_deseq:
    # 有重复
    import pandas as pd
    if not os.path.exists("DESeq"):
        os.mkdir("DESeq")
    count_table = pd.read_csv("featureCounts/featureCounts.txt", comment="#", sep='\t')
    count_table = count_table.rename(columns=lambda x: re.sub(r"bowtie2/(.*)\.bam", r"\1", x))
    ncol = count_table.shape[1]
    target_col = [0] + list(range(6, ncol))
    count_table = count_table.iloc[:, target_col]
    count_table.to_csv('featureCounts/count_table.xls', sep='\t', index=False)
    with open("DESeq/coldata.xls", "w") as f:
        print("sample,condition,batch", file=f)
        for k,v in config["design"].items(): # k: condition name | v: list of sample
            for i in v: # i: sample name
                print("{i},{k},{b}".format(i=i, k=k, b=config["batch"][i]), file=f)
    contrasts = " ".join(['_vs_'.join(i) for i in config['contrast']])
    os.system("Rscript {}/DESeq.r {}".format(script_base, contrasts))

if "comparisons" in config and not args.no_gfold:
    # 无重复
    # prepare input files for gfold
    os.system("if [ ! -d gfold ]; then mkdir gfold; fi")
    os.system("Rscript {}/featureCounts_to_gfold.r".format(script_base)) # prepare .count file
    os.system("awk -F '\t' '{print $5\"\t\"$3}' reference.idmap |tail +2 >gfold/gene.des") # prepare .des file
    
    # run gfold
    comparisons = config["comparisons"] # a/b
    
    # run gfold
    os.chdir("gfold")
    for a,b in comparisons:
        os.system("gfold diff -s1 {b} -s2 {a} -suf .count -o {a}_vs_{b}.xls -d gene.des".format(a=a, b=b))
    
    # add gene product to end of each line. 
    mapping = {}
    for line in open("gene.des"):
        sline = line.rstrip().split('\t')
        mapping[sline[1]] = sline[0]
    
    for a,b in comparisons:
        with open(f"{a}_vs_{b}.anno.xls", "w") as f:
            print("GeneSymbol,GeneName,GFOLD(0.01),E-FDR,log2fdc,1stRPKM,2ndRPKM,product".replace(",", "\t"), file=f)
            for line in open(f"{a}_vs_{b}.xls"):
                if line.startswith('#'):
                    continue
                sline = line.rstrip().split('\t')
                sline.append(mapping[sline[0]])
                print('\t'.join(sline), file=f)
