#!/usr/bin/env python3

import os
import sys
import json
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
    if os.path.exists("fastp/success"):
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
    if os.path.exists("bowtie2/success"):
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
    os.system("touch bowtie2/success")

def featureCount(samples):
    if os.path.exists("featureCounts/success"):
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

if len(sys.argv) == 1 or '-h' in sys.argv:
    print("原核生物转录组测序数据分析。使用gfold进行差异表达分析。")
    print("使用方法：")
    print("rnaseq.py config.json")
    exit()


######################
# Load configuration # 
######################
config_file = sys.argv[1]
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
# prepare input files for gfold
os.system("if [ ! -d gfold ]; then mkdir gfold; fi")
os.system("Rscript {}/featureCounts_to_gfold.r".format(script_base)) # prepare .count file
os.system("python {}/gff2idmap.py {} >reference.idmap".format(script_base, gff))
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
        print("GeneSymbol,GeneName,GFOLD(0.01),E-FDR,log2fdc,1stRPKM,2ndRPKM,product", file=f)
        for line in open(f"{a}_vs_{b}.xls"):
            if line.startswith('#'):
                continue
            sline = line.rstrip().split('\t')
            sline.append(mapping[sline[0]])
            print('\t'.join(sline), file=f)
