#!/usr/bin/env python3

# 1. create output folder. 
# 2. copy files(use hard link to speed up) to output folder. 
# 3. create a new conf.json file, with all file paths changed into the paths inside docker container. 
# 4. pass output folder to docker entry. 

import os
import json
import sys

def change_path(s, p):
    '''
    change path, relative to output directory. 
    '''
    pass

if len(sys.argv) == 1 or '-h' in sys.argv:
    print(f"\nusage:\n\n{os.path.realpath(__file__)} conf.json\n")
    exit(1)

conf = json.load(open(sys.argv[1]))

outdir = os.path.abspath(conf["output_dir"])
inputdir = os.path.join(outdir, "input")
mount_dir = "/home/output"

# if os.path.exists(outdir):
#     sys.stderr.write(f"output directory [{outdir}] already exists! \n")
#     exit(1)
# else:
if not os.path.exists(outdir):
    os.mkdir(outdir)
    os.mkdir(inputdir)

# Structure of conf: 
# "output_dir": "/output/",
#   "genome_sequence": "rnaseq_example/GCF_001891105.1_ASM189110v1_genomic.fna",
#   "genome_annotation": "rnaseq_example/GCF_001891105.1_ASM189110v1_genomic.gff",
#   "reads_data": [

conf["output_dir"] = mount_dir
fna = conf["genome_sequence"]
gff = conf["genome_annotation"]

tmp = os.path.join(inputdir, "genome.fna")
os.system(f"cp -l {fna} {tmp}")
conf["genome_sequence"] = f"{mount_dir}/input/genome.fna"
tmp = os.path.join(inputdir, "genome.gff")
conf["genome_annotation"] = f"{mount_dir}/input/genome.gff"
os.system(f"cp -l {gff} {tmp}")

for i, r in enumerate(conf['reads_data']):
    n, r1, r2 = r
    tmp = os.path.join(inputdir, f"{n}.1.fq.gz")
    os.system(f"cp -l {r1} {tmp}")
    conf['reads_data'][i][1] = f"{mount_dir}/input/{n}.1.fq.gz"
    tmp = os.path.join(inputdir, f"{n}.2.fq.gz")
    os.system(f"cp -l {r2} {tmp}")
    conf['reads_data'][i][2] = f"{mount_dir}/input/{n}.2.fq.gz"

with open(f"{outdir}/conf.json", "w") as f:
    json.dump(conf, f, indent=4)


uid = os.getuid()
gid = os.getgid()

os.system(f"docker run -u {uid}:{gid} -v {outdir}:{mount_dir} lch/rnaseq_b:v0.1 /opt/rnaseq/rnaseq.py -c {mount_dir}/conf.json")
