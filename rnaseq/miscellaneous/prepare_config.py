#!/usr/bin/env python3

import glob
import json
import os
import argparse

script_dir = os.path.dirname(os.path.abspath(__file__))

parser = argparse.ArgumentParser()
parser.description = "准备config文件"

parser.add_argument("-e", "--email", default="a879942613@qq.com", help="email")
parser.add_argument("-i", "--indir", help="input directoty. a sample a folder, folder name is sample name, fastq file suffixed with _{1,2}.fq.gz")
parser.add_argument("--suffix", default="1.fq.gz", help="suffix of fastq file. ")
parser.add_argument("-o", "--outdir", help="output directoty")
parser.add_argument("--gff", help="gff file")
parser.add_argument("--seq", help="genome sequence file")
parser.add_argument("--nosnv", action="store_true", help="do not run snv calling pipeline")
parser.add_argument("--nosv", action="store_true", help="do not run sv calling pipeline")
parser.add_argument("--json", default="conf.json", help="output json file")

args = parser.parse_args()

out_json = {}

out_json["email"] = args.email
out_json["out_dir"] = os.path.abspath(args.outdir)
out_json["SNV"] = "true" if not args.nosnv else "false"
out_json["SV"] = "true" if not args.nosv else "false"

indir = os.path.abspath(args.indir)
samples = [i for i in os.listdir(indir) if os.path.isdir(f"{indir}/{i}")]
out_json["samples"] = []
for i in samples:
    fq1 = glob.glob(f"{indir}/{i}/*{args.suffix}")[0]
    fq2_suffix = args.suffix.replace("1", "2")
    fq2 = glob.glob(f"{indir}/{i}/*{fq2_suffix}")[0]
    out_json["samples"].append([i, fq1, fq2])


raw_gff = os.path.abspath(args.gff)
raw_seq = os.path.abspath(args.seq)

# check chromosome name of gff file and seq file are the same. 

def merge_file(outfile, f1, f2):
    with open(outfile, "w") as f:
        f1 = open(f1).read().rstrip("\n")
        f2 = open(f2).read().rstrip("\n")
        print(f1, file=f)
        print(f2, file=f)

new_gff = "added_plasmid.gff"
new_seq = "added_plasmid.fna"
merge_file(new_gff, raw_gff, f"{script_dir}/pBT20.gff")
merge_file(new_seq, raw_seq, f"{script_dir}/pBT20_true.fa")

out_json["reference_gff"] = os.path.abspath(new_gff)
out_json["reference_fasta"] = os.path.abspath(new_seq)

with open(args.json, "w") as f:
    json.dump(out_json, f, indent=4, sort_keys=True)
