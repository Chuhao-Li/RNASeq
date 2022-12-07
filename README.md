# Bacterial RNASeq analysis pipeline

---

version v0.2.3

update: fix some bugs. 

# Install: 

This pipeline is developed in ubuntu. 

First, you need to clone the code: 

```
git clone https://github.com/Chuhao-Li/RNASeq.git
```

You are recommended to install with Dockerfile now. I did not push the image into the dockerhub yet, 
so you need to build it from Dockerfile: 

```bash
docker build -f Dockerfile -t lch/rnaseq_b:v0.1 .
```

Once it is finished, you can check it by `docker image ls`: 

```
REPOSITORY                 TAG       IMAGE ID       CREATED             SIZE
lch/rnaseq_b               v0.1      952ccaccc2e7   54 minutes ago      1.59GB
```

If you want to install all the dependency by yourself, you may need to ensure the 
following softwares are ready in your environment:

- fastp
- samtools
- bowtie2
- featureCounts
- gfold
- python3
- python packages: jinja2 pandas
- Rscript
- R packages: ggplot2 patchwork reshape2 ggdendro argparse DESeq2 sva

# How to use?

If you install by docker, you can run this pipeline like this: 

```bash
python rnaseq_docker_entry.py config.json
```

If you install it manually, you can run by: 

``` bash
./rnaseq/rnaseq.py -c config.json
```

Here is an example conf.json file: 

```
{
  "output_dir": "/path/to/output_directory",
  "genome_sequence": "path/to/genome_sequence.fna",
  "genome_annotation": "path/to/genome_annotation.gff",
  "reads_data": [
      ["X-1", "/path/to/read1.fq.gz", "/path/to/read2.fq.gz"],
      ["CK-1", "/path/to/read1.fq.gz", "/path/to/read2.fq.gz"],
      ["CK-2", "/path/to/read1.fq.gz", "/path/to/read2.fq.gz"],
      ["X-2", "/path/to/read1.fq.gz", "/path/to/read2.fq.gz"]
  ],
  "batch": {"X-2": 2, "CK-1": 1, "X-1": 1, "CK-2": 2},
  "design": {
      "X": ["X-1", "X-2"],
      "CK": ["CK-1", "CK-2"]
  },
  "contrast": [
      ["X", "CK"]
  ]
}

```

**Mandatory parameters**: 

- output_dir: output directoty. Result files will output to this directory. 
- genome_sequence: path to genome sequence file in fasta format. 
- genome_annotation: path to genome annotation file in gff format. Please ensure that the sequence ids in genome sequence file are the same as that in genome annotation file. Each CDS feature should contatin a `locus_tag` attribute. 
- reads_data: One sample a line, the three columns represents [sample name, path to read1 NGS data, path to read2 NGS data]. Here, the sample name should be unique. 

**Other parameters**: 

Other parameters depends on your experiment design. In this example, we have 2 condition. Two repeat in each conditions. 
For expriments that has repeat, you should provide the `design` and `contrast` parameters: 

- **design**: `"X": ["X-1", "X-2"]` means sample `X-1` and `X-2` is belong to `X` group. 
- **contrast**: `["X", "CK"]` means I want to contrast `X` group and `CK` group. As the `comparisons` parameter, order is meaningful when calculating the log2foldchange. 

If you want to correct batch effect, you need to specify the `batch` parameter: 
- **batch**: `{"X-2": 2, "CK-1": 1, "X-1": 1, "CK-2": 2}` means `CK-1` and `X-1` are from first batch, `X-2` and `CK-2` are from second batch. 


If your experiment has no repeat, you should provide the `comparisons` parameter: 

In this case, your conf.json file should looked like this: 
```
{
  "output_dir": "/path/to/output_directory",
  "genome_sequence": "path/to/genome_sequence.fna",
  "genome_annotation": "path/to/genome_annotation.gff",
  "reads_data": [
      ["X-1", "/path/to/read1.fq.gz", "/path/to/read2.fq.gz"],
      ["CK-1", "/path/to/read1.fq.gz", "/path/to/read2.fq.gz"]
  ],
  "comparisons": [
      ["X-1", "CK-1"]
  ]
}

```

- **comparisons**: One contrast one line. For example, `["X-1", "CK-1"]` means I want to compare X-1 and CK-1 to get differential expression gene. The order is meaningful when calculating the log2foldchange, the first sample is devided by the second sample. 


# Output

Output files are in `report` directory, which includes alignment quality, sample PCA, and differential expression analysis result. You can double click the index.html and read it in your web browser. 

