#!/usr/bin/env python3

import argparse
import os
import json

parser = argparse.ArgumentParser()
parser.description = "Call snp from RNASeq data"
parser.add_argument("-c", "--config", help="config file")
parser.add_argument("-o", "--outdir", default="snv", help="output directory")
parser.add_argument("-b", "--bamdir", default="bowtie2", help="bam file output directory")
parser.add_argument("-i", "--include", default="all", help="included samples, separated by comma. default all samples. ")
args = parser.parse_args()

config = json.load(open(args.config))
fa = config["genome_sequence"]
os.system(f"samtools faidx {fa}")

project_dir = os.path.abspath(config["output_dir"])

all_samples = [i[0] for i in config["reads_data"]]
if args.include == "all":
    samples = all_samples
else:
    samples = args.include.split(',')
    assert (len(samples) >=1) and (set(samples) & set(all_samples) == set(samples)), "some samples you specified not in config file! "

outdir = os.path.abspath(f"{project_dir}/{args.outdir}")
bamdir = os.path.abspath(f"{project_dir}/{args.bamdir}")
if not os.path.exists(outdir):
    os.mkdir(outdir)

for sample in samples:
    assert os.path.exists(f"{bamdir}/{sample}.bam"), f"bam file of sample {sample} not exist! "
    if not os.path.exists(f"{outdir}/{sample}.vcf"):
        os.system(f"samtools mpileup -g -f {fa} {bamdir}/{sample}.bam > {outdir}/{sample}.1.bcf")
        os.system(f"bcftools call -c --ploidy 1 {outdir}/{sample}.1.bcf> {outdir}/{sample}.2.bcf")
        os.system(f"bcftools view {outdir}/{sample}.2.bcf | vcfutils.pl varFilter - > {outdir}/{sample}.vcf")

# snv anno
gff = os.path.abspath(config['genome_annotation'])
annovar_path = "/home/lch/coding/pipeline/resequence/annovar"

for sample in samples:
    vcf_snp = f"{outdir}/{sample}.vcf"
    if not os.path.exists(f"{outdir}/{sample}/db"):
        os.makedirs(f"{outdir}/{sample}/db")
    genepred = f"{outdir}/{sample}/db/{sample}_refGene.txt"
    refGeneMrna = f"{outdir}/{sample}/db/{sample}_refGeneMrna.fa"
    avinput_snp = f"{outdir}/{sample}/{sample}_snp.avinput"
    func_anno_snp = f"{outdir}/{sample}/{sample}_snp.exonic_variant_function"
    os.system(f"{annovar_path}/gff3ToGenePred {gff} {genepred}")
    os.system(f"{annovar_path}/retrieve_seq_from_fasta.pl --format refGene --seqfile {fa} {genepred} --out {refGeneMrna}")
    os.system(f"{annovar_path}/convert2annovar.pl --format vcf4 {vcf_snp} >{avinput_snp}")
    os.system(f"{annovar_path}/annotate_variation.pl -buildver {sample} -out {outdir}/{sample}/{sample}_snp {avinput_snp} {outdir}/{sample}/db/")


