#!/usr/bin/env python3

# Description: 
# 1. 在每行末尾增加web_igv snapshot的链接。
# 2. 输出igv batch script
# 3. 转化为html表格。

# Input: 
# result file of special_SNV.py or special_SV.py

# Output: 
# 1. snapshot/snapshot.igv
# 2. index.html

# To do:
# 2. SNP报告
# 整合header
# 获取bam文件列表
# 4. multiqc报告

import os
import json
import sys
import glob
import argparse
from Bio import SeqIO

def row2tr(row, thead=False):
    '''
    row: a list
    '''
    html = ["<tr>"]
    if thead:
        cell = "th"
        for i in row:
            html.append('<{cell} class="{i}">{i}</{cell}>'.format(cell=cell, i=i))
        html.append("</tr>")
    else:
        cell = "td"
        for i in row:
            html.append('<{cell}>{i}</{cell}>'.format(cell=cell, i=i))
        html.append("</tr>")
    return html

def get_line_type(line):
    if line.startswith('-'):
        line = line.strip('-\n')
        if not line:
            return "blank"
        else:
            return "chrom"
    else:
        return "normal"

parser = argparse.ArgumentParser()
parser.description = "create SNV/SV report"
parser.add_argument("-c", "--config", help="config file")
parser.add_argument("-i", "--included", help="included samples")
parser.add_argument("--sv", help="grouped sv file. ")
parser.add_argument("--snv", default="snv/grouped_snv.txt", help="grouped snv file. ")
args = parser.parse_args()

# Initialize variables
config = json.load(open(args.config))
proj_dir = os.path.abspath(config["output_dir"])
outfile = f"{proj_dir}/snapshot/snapshot.igv"
outdir = os.path.dirname(outfile)
outdir_base = os.path.basename(outdir)
if not os.path.exists(outdir):
    os.mkdir(outdir)

seq = os.path.abspath(config["genome_sequence"])
gff = os.path.abspath(config["genome_annotation"])

all_samples = [i[0] for i in config["reads_data"]]
if args.included == "all":
    samples = all_samples
else:
    samples = args.included.split(',')
    assert (len(samples) >=1) and (set(samples) & set(all_samples) == set(samples)), "some samples you specified not in config file! "

bams = [f"{proj_dir}/bowtie2/{sample}.bam" for sample in samples]


# Prepare igv batch script header
chr_len_map = {}
for rec in SeqIO.parse(seq, 'fasta'):
    chr_len_map[rec.id] = len(rec.seq)

out_fh = open(outfile, 'w')

print("new", file=out_fh)
print("genome", seq, file=out_fh)
for bam in bams:
    base = os.path.basename(bam)
    print("load", bam, file=out_fh)
    print("viewaspairs", base, file=out_fh)
    print("collapse", base, file=out_fh)
print("load", gff, file=out_fh)
print("expand", os.path.basename(gff), file=out_fh)
print("snapshotDirectory", outdir, file=out_fh)

####################
# Output html file #
####################

# Prepare html header
document = []
script_dir = os.path.dirname(os.path.realpath(__file__))
print("script_dir", script_dir)
head = open('{}/report_header.html'.format(script_dir)).read()
document.append(head)

# Prepare html content
document.append('<body>')

# Prepare header
document.append('<h1 style="margin-left:150px">重测序变异检测分析报告</h1><hr>')

# Prapare TOC
toc = '''
<ul class="verticalnav">
    <li><a class="current" href="#sv_report">SV</a></li>
    <li><a href="#snv_report">SNP</a></li>
</ul>
'''
document.append(toc)

document.append('<div style="margin-left:150px; width:65%">')

#####################
# Prepare SV report #
#####################

# var = "SV"
# infile = 'report/{}/clustered.{}.txt'.format(var, var)
# if not os.path.exists(infile):
#     sys.stderr.write("input file not found! \n")
#     exit()
# data = open(infile)
# 
# # Prepare SV report title
# document.append('<div id="sv_report">')
# document.append('<h2>SV</h2>')
# document.append('''<p>每个样品都与对应的参考序列进行了比对，使用samtools得到低覆盖度区域(使用source=samtools标注的SV)，使用delly获得了SVs(使用source=delly标注的SV)。
# 按照在参考基因组序列上的位置对SVs进行了排列，间隔大于100bp的用空行隔开形成多个小组。同一个SV应该落在同一个小组中。
# 如果某个位置上，所有样品都有同样的SV，那可能是野生型本来就有的突变。</p>
# <p>各列的含义：</p>
# <ol>
# <li>chromosome: 变异所在的染色体ID</li>
# <li>start: 变异起始位点</li>
# <li>end: 变异终止位点</li>
# <li>length: 变异长度</li>
# <li>sv_type: 变异类型除了vcf文件定义的变异类型外，还添加了LOW和HIGH两种，分别代表低reads深度区域和高reads深度区域。</li>
# <li>others: 其他信息。
# <ul>
# <li>sample: 样品名</li>
# <li>source: 变异检测所用工具（delly/samtools）</li>
# <li>PE: 如果是delly，有PE属性，代表支持该变异的reads数</li>
# <li>depth: 如果是samtools，有depth参数，代表该区域的平均reads深度</li>
# </ul>
# </li>
# <li>detail: 包含了变异对应基因组区域的igv截图。长度大于8000的变异，只对左端(l)、中间(m)和右端(r)进行截图。</li>
# </ol>
# ''')
# 
# # Prepare table head
# document.append('<table class="pure-table pure-table-bordered" style="table-layout:fixed; width:100%">')
# document.append('<thead>')
# thead = next(data).rstrip() + "\tdetail"
# document.extend(row2tr(thead.split('\t'), thead=True))
# document.append('</thead>')
# 
# # Prepare table body and igv batch script actions. 
# # 若指定区域，igv最少显示11个碱基。如果指定的区域长度低于10，则不跳转。
# # 若指定单个位点，igv会自动往前后添加20个碱基。如果处于基因组前20或者后20个碱基，则不会跳转。
# # igv显示alignment的最大区间约为69000，大于这个区间，则要求放大才能显示。比较合适的区间是8000左右。
# 
# document.append('<tbody>')
# n = 1 # define variant ID
# for line in data:
#     if var == "SV":
#         lt = get_line_type(line)
#         if lt == "blank":
#             document.append('<tr><td colspan=7 style="background-color: #000000"></td></tr>')
#         elif lt == "chrom":
#             document.append('<tr><td colspan=7 style="background-color: #000000; color: white"><b>{}</b></td></tr>'.format(line.strip("\n-")))
#         else:
#             m = str(n).zfill(6)
#             # region: (chro, start, end)
#             sline = line.rstrip('\n').split('\t')
#             chro, start, end = sline[:3]
#             start = int(start)
#             end = int(end)
#             l = end - start + 1
#             chr_len = chr_len_map[chro]
#             # 两边都扩展500bp
#             start1 = start - 500
#             end1 = end + 500
#             # 若落在左端且小于1000，则取1000。若落在右端且
#             if end1 < 1000:
#                 end1 = 1000
#             if start1 > chr_len - 1000:
#                 start1 = chr_len - 1000
#             if start1 < 1:
#                 start1 = 1
#             if end1 > chr_len:
#                 end1 = chr_len
#             l2  = end1 - start1 + 1
#             # 若区间总长大于8000，则取8000，再加两端各1000
#             print("region {chro} {start} {end}".format(chro=chro, start=start, end=end), file=out_fh)
#             if l2 > 8000:
#                 left = (start1, start1+1000)
#                 mid = (start1+end1)//2
#                 middle = (mid-4000, mid+4000)
#                 right = (end1-1000, end1)
#                 cell = []
#                 for part, position in zip("lmr", [left, middle, right]):
#                     s, e = position
#                     print("goto {chro}:{start}-{end}".format(chro=chro, start=s, end=e), file=out_fh)
#                     fig_name = "{}_{}_{}.png".format(var, m, part)
#                     print("snapshot", fig_name, file=out_fh)
#                     cell.append('<a href="{outdir}/{fig_name}" style="padding:0 5px;" target="_blank">{part}</a>'.format(outdir=outdir_base, fig_name=fig_name, part=part))
#                 cell = "".join(cell)
#             else:
#                 print("goto {chro}:{start}-{end}".format(chro=chro, start=start1, end=end1), file=out_fh)
#                 fig_name = "{}_{}.png".format(var, m)
#                 print("snapshot", fig_name, file=out_fh)
#                 cell = '<a href="{outdir}/{fig_name}" style="padding:0 5px;" target="_blank">m</a>'.format(outdir=outdir_base, fig_name=fig_name)
#             sline[-1] = sline[-1].replace(';', '; ')
#             sline.append(cell)
#             document.extend(row2tr(sline))
#             n += 1
# 
# document.append('</tbody>')
# document.append('</table>')
# document.append('</div>')

######################
# Prepare SNV report #
######################
var = "SNV"
infile = f'{proj_dir}/{args.snv}'
if not os.path.exists(infile):
    sys.stderr.write("input file not found! \n")
    exit()
data = open(infile)

# Prepare SNV report title
document.append('<div id="snv_report">')
document.append('<h2>SNV</h2>')
document.append('''<p>结果说明： 每个样品都与对应的参考序列进行了比对，使用bcftools获得了SNVs。
按照在参考基因组序列上的位置对SNVs进行了排列，间隔大于100bp的用空行隔开形成多个小组。同一个SNV应该落在同一个小组中。
如果某个位置上，所有样品都有同样的SNV，那可能是野生型本来就有的突变。</p>
<p>各列的含义：</p>
<ol>
<li>sample: 样品名</li>
<li>chromosome: 变异所在染色体ID</li>
<li>position: 变异位置</li>
<li>REF: 参考基因组上对应位置的碱基</li>
<li>ALT: 样品基因组上对应位置的碱基信息</li>
<li>annotation：变异注释信息。
<ul>
<li>SNVType: SNV变异类型，分intergenic、synonymous、nonsynonymous几种。</li>
<li>ChangeGene: 变异所在基因</li>
<li>SNVScore: bcftools计算出来的变异分数，越高分代表SNV越可靠。</li>
<li>readDepth: SNV所在位点的reads深度</li>
<li>product: 变异所在基因的注释</li>
</ul></li>
<li>detail: 包含了变异对应基因组区域的igv截图。</li>
</ol>
''')

# Prepare table head
document.append('<table class="pure-table pure-table-bordered" style="table-layout:fixed; overflow: scroll; width:100%">')
document.append('<thead>')
thead = next(data).rstrip() + "\tdetail"
document.extend(row2tr(thead.split('\t'), thead=True))
document.append('</thead>')

# Prepare table body and igv batch script actions. 
document.append('<tbody>')
n = 1 # define variant ID
for line in data:
    if var == "SNV":
        lt = get_line_type(line)
        if lt == "blank":
            document.append('<tr><td colspan=7 style="background-color: #000000"></td></tr>')
        elif lt == "chrom":
            document.append('<tr><td colspan=7 style="background-color: #000000; color: white"><b>{}</b></td></tr>'.format(line.strip("\n-")))
        else:
            m = str(n).zfill(6)
            # region: (chro, start, end)
            sline = line.rstrip('\n').split('\t')
            chro, start = sline[1:3]
            start = int(start)
            chr_len = chr_len_map[chro]
            # 两边都扩展300bp
            start1 = start - 300
            end1 = start + 300
            # 若落在左端且小于600，则取600。若落在右端且
            if end1 < 600:
                end1 = 600
            if start1 > chr_len -600:
                start1 = chr_len - 600
            if start1 < 1:
                start1 = 1
            if end1 > chr_len:
                end1 = chr_len
            l2  = end1 - start1 + 1
            # 感兴趣的区域，前后各加2，避免遮挡。
            print("region {chro} {start} {end}".format(chro=chro, start=start-2, end=start+2), file=out_fh)
            print("goto {chro}:{start}-{end}".format(chro=chro, start=start1, end=end1), file=out_fh)
            fig_name = "{}_{}.png".format(var, m)
            print("snapshot", fig_name, file=out_fh)
            cell = '<a href="{outdir}/{fig_name}" style="padding:0 5px;" target="_blank">m</a>'.format(outdir=outdir_base, fig_name=fig_name)
            sline[-1] = sline[-1].replace(';', '; ')
            sline.append(cell)
            document.extend(row2tr(sline))
            n += 1

document.append('</tbody>')
document.append('</table>')
document.append('</div>') # end of SNV report
document.append('</div>') # end of tables
document.append('</body>')

script = '''
<script>
    toc = document.querySelector(".verticalnav"); 
    items = toc.querySelectorAll("a")
    for (i = 0; i < items.length; i++) {    
        items[i].onclick = function(){
            document.querySelector(".current").className = ""; 
            this.className = "current"; 
        }
    }

    //鼠标经过时，单元格里面的文本自动换行。只有遇到空格或者“-”的时候才会换行。
    tds = document.querySelectorAll("td"); 
    for (i=0; i<tds.length; i++){
        tds[i].onmouseover = function(){
            this.style["white-space"] = "normal";
        }
        tds[i].onmouseout = function(){
            this.style["white-space"] = "nowrap"; 
        }
    }
</script>
'''
document.append(script)
document.append('</html>')

with open(f'{proj_dir}/snv/resequence.html', 'w') as f:
    for i in document:
        f.write(i+'\n')
