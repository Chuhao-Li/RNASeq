#!/usr/bin/env python3
#
'''
0. 质控结果，multiqc页面一个。
1. 比对结果，一张比对情况图片就可以。 ok
2. 定量结果，文件，raw read count和定量之后的各一个文件。ok
3. 差异表达分析结果：
    无重复gfold：差异表达分析的结果文件，一对comparison对应一个文件。 ok
    有重复DESeq2：差异表达分析的结果文件，一对comparison对应一个文件。ok 
                  另外还有样品重复性图片一张。

差异基因筛选：火山图，每个比较组合，lfc1和lfc2各一张。
    差异基因列表：每个比较组合都有两个lfc阈值，然后每个lfc阈值分别得到上调和下调两个列表。

KEGG和GO富集分析：每个基因列表对应一个GO和一个KEGG的气泡图。

现在，gfold的结果没有输出火山图/差异基因。
DESeq2的结果只支持一对contrasts
DESeq2的差异表达分析和差异基因筛选分开会比较好。只输出差异表达分析的表格，让用户自己筛选。
DESeq2的重复性相关图片还没有输出到report目录中，以及index.html中。

输入、输出说明，以及工具流程说明都还没写。
筛选列表得到之后，可以用下游的工具进行图片绘制、富集分析等。
'''

import re
import os
import sys
import json
import argparse
import pandas as pd
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
    index = "bowtie2/genome_index"
    os.system("bowtie2-build {fna} {index}".format(fna=fna, index=index))
    command = '''bowtie2 -x {index} \
        -1 fastp/{s}_1.fq.gz \
        -2 fastp/{s}_2.fq.gz \
        -p 16 |samtools sort -o bowtie2/{s}.bam
    '''.format(index=index, s=s)
    os.system(command)

def featureCount(samples, attribute="locus_tag"):
    if os.path.exists("featureCounts/featureCounts.txt"):
        return
    os.system("if [ ! -d featureCounts ]; then mkdir featureCounts; fi")
    command = '''featureCounts \
      -a {gff} \
      -o featureCounts/featureCounts.txt \
      -s 2 -p -t CDS -g {attribute} -T 15 \
      {bams}
    '''.format(gff=gff, attribute=attribute, bams = " ".join(["bowtie2/{}.bam".format(i) for i in samples]))
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
    print("rnaseq v0.2.2")
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
os.system("if [ ! -d fastp ]; then mkdir fastp; fi")
with Pool(processes=n_samples) as p:
    p.map(fastp, samples)
    p.close()
    p.join()

print("alinging...")
os.system("if [ ! -d bowtie2 ]; then mkdir bowtie2; fi")
for s in samples:
    align(s)

print("featureCounts...")
featureCount(samples)

os.chdir(f'{workdir}/featureCounts')
os.system(f"Rscript {script_base}/rscript/Plot_description_of_read_alignment.r")
os.chdir(workdir)

####################################
# Differential expression analysis #
####################################
if args.no_dea:
    exit()

os.system("python3 {}/gff2idmap.py {} >reference.idmap".format(script_base, gff))
        

def create_count_table():
    count_table = pd.read_csv("featureCounts/featureCounts.txt", comment="#", sep='\t')
    count_table = count_table.rename(columns=lambda x: re.sub(r"bowtie2/(.*)\.bam", r"\1", x))
    ncol = count_table.shape[1]
    target_col = [0] + list(range(6, ncol))
    count_table = count_table.iloc[:, target_col]
    count_table.to_csv('featureCounts/count_table.xls', sep='\t', index=False)

if "design" in config and not args.no_deseq:
    # 有重复
    if not os.path.exists("DESeq/success"):
        if not os.path.exists("DESeq"):
            os.mkdir("DESeq")
        create_count_table()
        with open("DESeq/coldata.xls", "w") as f:
            if "batch" in config:
                print("sample,condition,batch", file=f)
                for k,v in config["design"].items(): # k: condition name | v: list of sample
                    for i in v: # i: sample name
                        print("{i},{k},{b}".format(i=i, k=k, b=config["batch"][i]), file=f)
            else:
                print("sample,condition", file=f)
                for k,v in config["design"].items(): # k: condition name | v: list of sample
                    for i in v: # i: sample name
                        print("{i},{k}".format(i=i, k=k), file=f)
                
        contrasts = ",".join(['_vs_'.join(i) for i in config['contrast']])
        command = (f"Rscript {script_base}/rscript/DESeq.r " + 
                   "--count_table featureCounts/count_table.xls " + 
                   f"--coldata DESeq/coldata.xls --contrast {contrasts} " + 
                   "-o DESeq --anno reference.idmap")
        r = os.system(command)
        if r == 0:
            os.system("touch DESeq/success")

if "comparisons" in config and not args.no_gfold:
    # 无重复
    if not os.path.exists("gfold/success"):
        # prepare input files for gfold
        os.system("if [ ! -d gfold ]; then mkdir gfold; fi")
        create_count_table()
        with open("gfold/coldata.xls", "w") as f:
            print("sample,condition", file=f)
            for i in samples:
                print(f"{i},{i}", file=f)
        os.system("Rscript {}/rscript/normalize_by_DESeq.r".format(script_base))
        os.system("Rscript {}/rscript/featureCounts_to_gfold.r".format(script_base)) # prepare .count file
        os.system("awk -F '\t' '{print $5\"\t\"$3}' reference.idmap |tail +2 >gfold/gene.des") # prepare .des file
        
        # run gfold
        comparisons = config["comparisons"] # a/b
        os.chdir("gfold")
        for a,b in comparisons:
            os.system("gfold diff -s1 {b} -s2 {a} -suf .count -o {a}_vs_{b}.xls -d gene.des".format(a=a, b=b))
        
        # add gene product to end of each line. 
        mapping = {}
        for line in open("gene.des"):
            sline = line.rstrip("\n").split('\t')
            mapping[sline[1]] = sline[0]
        
        for a,b in comparisons:
            with open(f"{a}_vs_{b}.anno.xls", "w") as f:
                print("GeneSymbol,GeneName,GFOLD(0.01),E-FDR,log2fdc,1stRPKM,2ndRPKM,product".replace(",", "\t"), file=f)
                for line in open(f"{a}_vs_{b}.xls"):
                    if line.startswith('#'):
                        continue
                    sline = line.rstrip("\n").split('\t')
                    sline.append(mapping[sline[0]])
                    print('\t'.join(sline), file=f)
        flag = 0
        for a,b in comparisons:
            if not os.path.exists(f"gfold/{a}_vs_{b}.anno.xls"):
                flag == 1
        if flag == 0:
            os.system("touch success")

# create report
report_dir = f"{workdir}/report"
if not os.path.exists(report_dir):
    os.mkdir(report_dir)
os.system(f"cp -r {script_base}/static {report_dir}")
fig_dir = f"{report_dir}/figure"
data_dir = f"{report_dir}/data"
if not os.path.exists(fig_dir):
    os.mkdir(fig_dir)
if not os.path.exists(data_dir):
    os.mkdir(data_dir)
os.system(f"cp {workdir}/featureCounts/Description_of_read_alignment.png {fig_dir}")

if "comparisons" in config and not args.no_gfold:
    os.system(f"cp {workdir}/gfold/normalized_readCount.csv {data_dir}")
    comparisons = config["comparisons"]
    for a,b in comparisons:
        os.system(f"cp {workdir}/gfold/{a}_vs_{b}.anno.xls {data_dir}")

if "design" in config and not args.no_deseq:
    os.system(f"cp {workdir}/DESeq/normalized_readCount.csv {data_dir}")
    os.system(f"cp {workdir}/DESeq/pall.png {fig_dir}")
    for a,b in config['contrast']:
        os.system(f"cp {workdir}/DESeq/{a}_vs_{b}.csv {data_dir}")

from jinja2 import FileSystemLoader, Environment
# 参考：https://www.jianshu.com/p/e1d36fe64e37

# 加载模板文件夹
loader = FileSystemLoader(searchpath=f'{script_base}/templates')
# 环境对象
enviroment = Environment(loader=loader)
# 指定模板文件
tpl = enviroment.get_template('index.html')

# 渲染模板
if "comparisons" in config:
    comparisons = [f"{a}_vs_{b}" for a,b in config["comparisons"]]
elif "design" in config:
    comparisons = [f"{a}_vs_{b}" for a,b in config["contrast"]]

repeat = "true" if "design" in config else "false"

output = tpl.render(comparisons=comparisons, repeat=repeat)
with open(f"{report_dir}/index.html", "w") as f:
    print(output, file=f)
