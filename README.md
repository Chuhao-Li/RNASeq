# Bacterial RNASeq analysis pipeline

version v0.2.1

### Prequisite: 

Please ensure the following softwares are ready in your environment. 
- python3
- fastp
- bowtie2
- featureCounts
- gfold
- Rscript

### Install: 
```
git clone https://github.com/Chuhao-Li/RNASeq.git
```

### How to use?

The first step to use this pipeline is to edit the config file: 

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
  "comparisons": [
      ["X-1", "CK-1"],
      ["X-2", "CK-2"]
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
The config.json is a template of config file. 

Mandatory parameters: 
- output_dir: output directoty. Result files will output to this directory. 
- genome_sequence: path to genome sequence file in fasta format. 
- genome_annotation: path to genome annotation file in gff format. Please ensure that the sequence ids in genome sequence file are the same as that in genome annotation file. Each CDS feature should contatin a `locus_tag` attribute. 
- reads_data: One sample a line, the three columns represents sample name, path to read1 NGS data, path to read2 NGS data. 

If your experiment has no repeat, you should provide the `comparisons` parameter: 
- comparisons: One contrast one line. For example, `["X-1", "CK-1"]` means I want to compare X-1 and CK-1 to get differential expression gene. The order is meaningful when calculating the log2foldchange, the first sample is devided by the second sample. 

If your experiment has repeat, you should provide the `design` and `contrast` parameters: 
- design: `"X": ["X-1", "X-2"]` means sample `X-1` and `X-2` is belong to `X` group. 
- contrast: `["X", "CK"]` means I want to contrast `X` group and `CK` group. As the `comparisons` parameter, order is meaningful when calculating the log2foldchange. 

If you want to correct batch effect, you need to specify the `batch` parameter: 
- batch: `{"X-2": 2, "CK-1": 1, "X-1": 1, "CK-2": 2}` means `CK-1` and `X-1` are from first batch, `X-2` and `CK-2` are from second batch. 

As the config.json is ready, you can run the following command to run the pipeline: 
``` bash
./rnaseq.py -c config.json
```

### Other toolkits

`Plot_description_of_read_alignment.r` can be used to plot aligned ratio: 
``` bash
cd featureCounts
Rscript /path/to/Plot_description_of_read_alignment.r
```

`pollution.py` can be used to identified pollution source, but you need to modify this script to adapt to your need. 

`call_snv.py`, `group_snv.py` and `create_snv_report.py` can be used to identify snv. 
``` bash 
python call_snv.py -c config.json # call snv and annotation snv
python group_snv.py -c config.json # group snv 
python create_snv_report.py -c config.json # create html report 
```
