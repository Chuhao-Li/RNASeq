#!/usr/bin/env python3
# description: compare paired samples to get special snvs. 
import os
import re
import json
import argparse
import gffutils


def sample2vcf(sample):
    return f"{proj_dir}/snv/{sample}.vcf"

def sample2anno(sample):
    return f"{proj_dir}/snv/{sample}/{sample}_snp.exonic_variant_function"

def collect_snvs(samples):
    snvs_dict = {} # {'sample':('sample', chro, start, end, alt), ...}
    for sample in samples:
        snvs = open(sample2vcf(sample))
        snvs_dict[sample] = set()
        for snv in snvs:
            if snv.startswith('#'):
                continue
            snv = snv.split('\t')
            snv_tuple = (sample, snv[0], snv[1], snv[3], snv[4])
            snvs_dict[sample].add(snv_tuple)
    return snvs_dict

def create_gff_db(gfffile):
    dbfn=gfffile + '.db'
    if os.path.exists(dbfn):
        db = gffutils.FeatureDB(dbfn, keep_order=True)
    else:
        db = gffutils.create_db(gfffile, dbfn=dbfn, force=True, keep_order=True,
                merge_strategy='merge', sort_attribute_values=True)
    return db

def cluster_by_pos(snvs_dict, samples, gff, outfile):
    db = create_gff_db(gff)
    allsets = snvs_dict.values()
    all_snvs = [] # [('sample', chro, start, ref, alt), ...]
    anno = get_anno(samples)
    for i in allsets:
        for j in i:
            all_snvs.append(j)
    
    out = open(outfile, 'w')
    # print("# SNV(Single Nucleotide Variant) 结果文件。", 
    #         "# 每个样品都与对应的参考序列进行了比对，使用bcftools获得了SNVs。",
    #         "# 按照在参考基因组序列上的位置对SNVs进行了排列，间隔大于100bp的用空行隔开形成多个小组。同一个SNVs应该落在同一个小组中。", 
    #         "# 如果该小组中只包含了非野生型的SNVs，更有可能与表型相关。", 
    #         "# sample\tchromosome\tposition\tREF\tALT\tannotation", sep = '\n')
    print("sample\tchromosome\tposition\tREF\tALT\tannotation", file=out)
    last = 0
    last_chro = ''
    for i in sorted(all_snvs, key=lambda x: (x[1], int(x[2]))):
        this_chro = i[1]
        this = int(i[2])
        if this_chro != last_chro:
            print('-\n-{}\n-'.format(this_chro), file=out)
        elif abs(this - last) >= 100:
            print('-', file=out)
        if i in anno:
            gene = re.match(r'[^:#]+', anno[i][2]).group()
            try:
                product = next(db.children(gene, featuretype='CDS')).attributes['product'][0]
            except(StopIteration):
                product = 'Unknow'
            a = "SNVType={};ChangeGene={};SNVScore={};readDepth={};product={}".format(anno[i][1], 
                    gene, 
                    anno[i][-2],
                    anno[i][-1], 
                    product)
        else:
            a = 'SNVtype=intergenic'
        print('\t'.join(i) + '\t' + a, file=out)
        last_chro = this_chro
        last = this

def get_anno(samples):
    m = {} # {("chro", "pos", "ref", "alt"):"annotation", ...}
    for sample in samples:
        for line in open(sample2anno(sample)):
            line = line.rstrip().split('\t')
            m[(sample, line[3], line[4], line[6], line[7])] = line
    return m

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.description = "Assign SNVs in near position into the same group. "
    parser.add_argument("-o", "--outfile", default="snv/grouped_snv.txt", help="output file path")
    parser.add_argument("-c", "--config", help="config file")
    parser.add_argument("-i", "--included", default="all", help="included samples")
    args = parser.parse_args()
    config = json.load(open(args.config))
    gff = config["genome_annotation"]
    all_samples = [i[0] for i in config["reads_data"]]
    if args.included == "all":
        samples = all_samples
    else:
        samples = args.included.split(',')
    proj_dir = os.path.abspath(config["output_dir"])
    outfile = f"{proj_dir}/{args.outfile}"
    snvs_dict = collect_snvs(samples)
    cluster_by_pos(snvs_dict, samples, gff, outfile)

