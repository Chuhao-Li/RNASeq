import re
import sys

ingff = sys.argv[1]
outgff = sys.argv[2]

out = open(outgff, "w")

for line in open(ingff):
    if line.startswith("#"):
        print(line, end="", file=out)
        continue
    sline = line.rstrip('\n').split('\t')
    if sline[2] in {"gene", "pseudogene", "CDS"}:
        ID = re.search(r"(?<=ID=)[^;\n]+", sline[8]).group()
        print(line.rstrip(";\n") + ";locus_tag={}".format(ID.split(".")[0]), file=out)
    else:
        print(line, end="", file=out)

out.close()
