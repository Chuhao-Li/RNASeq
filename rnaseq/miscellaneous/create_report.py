#!/usr/bin/env python3
'''
整理结果文件到report文件夹，并渲染report.html。
input: 
    src.json # 资源文件
    report.html # html内容
    header.html # html模板
output:
    report directory, which include: 
        data directory 
        figure directory
        table directory
        static directory
        index.html
'''

import os
import sys
import json
import argparse
import bs4

script_dir = os.path.dirname(os.path.abspath(__file__))

def render(info, index, tag_type):
    rendered = None
    if "rename" in info:
        path = tag_type + "/" + info["rename"]
    else:
        path = tag_type + "/" + os.path.basename(info["path"])
    if "title" in info:
        caption = info["title"]
    else:
        caption = ""
    if "desc" in info:
        desc = info["desc"]
    else:
        desc = ""
    if tag_type == "figure":
        if "max-width" in info:
            max_width = info["max-width"]
        else:
            max_width = "80%"
        Tag_type = tag_type.capitalize()
        caption_node = f'<div class="fig_caption"><b>{Tag_type}.&nbsp;{index}&nbsp;{caption}</b></div>'
        content_node = f'<div class="fig_image"><img src="{path}" style="max-width:{max_width};"></div>'
        description_node = f'<div class="fig_description"><p><b>NOTE:&nbsp;</b>{desc}</p></div>'
        rendered = caption_node + content_node + description_node
    elif tag_type == "data":
        if not caption:
            caption = os.path.basename(path)
        rendered = f'<a href="{path}">Supplementary file.&nbsp;{index}&nbsp;{caption}</a>'
    elif tag_type == "table":
        pass
    rendered = bs4.BeautifulSoup(rendered, "html.parser")
    return rendered

def copy_files(tag_info, target):
    origin = tag_info["path"]
    if not os.path.exists(origin):
        sys.stderr.write(f"file not found: {origin}, skipped. ")
    else:
        if "rename" in tag_info:
            new_name = tag_info["rename"]
            target = f"{target}/{new_name}"
        os.system(f"cp -l {origin} {target}")

# parse arguments
parser = argparse.ArgumentParser()
parser.description = "整理结果文件到report文件夹，并渲染report.html"
parser.add_argument("--header", help="报告标题")
parser.add_argument("-o", "--outdir", default="report", help="输出目录")
parser.add_argument("-r", "--report", default="report.html", help="报告内容.html")
parser.add_argument("-s", "--src", default="src.json", help="标签信息.json")
parser.add_argument("-t", "--template", default=f"{script_dir}/template.html", help="报告模板.html")
args = parser.parse_args()

# load input files
template = bs4.BeautifulSoup(open(args.template), "html.parser")
body = bs4.BeautifulSoup(open(args.report), "html.parser") # 加上html.parser可以避免自动补全<html>和<head>标签。
src = json.load(open(args.src))
outdir = os.path.abspath(args.outdir)

if not os.path.exists(outdir):
    os.mkdir(outdir)

# render tags
for tag_type in ["figure", "data", "table"]:
    tags = body.find_all("div", tag_type) # 找出所有标签
    tag_outdir = f"{outdir}/{tag_type}"
    if len(tags) > 0 and not os.path.exists(tag_outdir):
        os.mkdir(tag_outdir)
    for index, tag in enumerate(tags):
        index += 1
        tag_id = tag.attrs['id']
        tag_info = src[tag_type][tag_id]
        copy_files(tag_info, tag_outdir)
        content = render(tag_info, index, tag_type)
        tag.append(content)

template.find_all("div", "article-body")[0].append(body)
with open(f'{outdir}/index.html', 'w') as f:
    print(template.prettify(), file=f)
