#!/usr/bin/python3
import sys
import re
# usage: 
usage = 'python3 script/locus_to_refseqID.py in.gff >out.xls'

if len(sys.argv) == 1 or '-h' in sys.argv:
    print(usage)
    exit()

def locusID_to_refseqID(gff):
    result = {}
    this_gene = ''
    tag = 0
    for line in open(gff):
        if line.startswith('#'):
            continue
        line = line.rstrip().split('\t')
        if line[2] == 'region':
            continue
        if line[2] in ['gene', 'pseudogene']:
            geneid = re.search('(?<=ID=)[^;]*', line[8]).group()
            gene_name = re.search('(?<=Name=)[^;]+', line[8])
            if gene_name:
                gene_name = gene_name.group()
            else:
                gene_name = ''
            locustag = re.search('(?<=locus_tag=)[^;]*', line[8])
            if locustag: 
                locustag = locustag.group()
            else:
                locustag = ''
            result[geneid] = [gene_name, locustag]
            this_gene = geneid
            tag = 0 # 对于每个gene，只读取第一个CDS，获取其product等。
        elif tag == 0 and line[2] == "CDS":
            try:
                parent = re.search('(?<=Parent=)[^;]*', line[8]).group()
            except(AttributeError):
                print(line)
            try:
                Name = re.search('(?<=Name=)[^;]*', line[8]).group()
            except(AttributeError):
                Name = ''
            try:
                product = re.search('(?<=product=)[^;]*', line[8]).group()
            except(AttributeError):
                product = ''
            if parent.startswith(this_gene):
                result[geneid].append(Name)
                result[geneid].append(product)
                tag = 1
    return result

if __name__ == '__main__':
    gff = sys.argv[1]
    # write mapping file
    result = locusID_to_refseqID(gff)
    print('geneId\tgeneName\tlocusTag\trefseqId\tproduct')
    for key,value in result.items():
        print(key + '\t' + '\t'.join(value))
