# gtf_features.py

A Python script for extracting feature information from a GTF file.

### input file

* GTF file containing `exon` features with `gene_id` and `transcript_id` attributes

### usage

```
usage: gtf_features.py [-h] {bed,count,length} ...

Extract feature information from GTF file

positional arguments:
  {bed,count,length}
    bed               Extract a BED3 file for the selected feature
    count             Count features per gene (i.e. exon, intron, transcript)
    length            Extract feature lengths (i.e. exon, intron, transcript, gene)

optional arguments:
  -h, --help          show this help message and exit
```

### `bed` module

```
usage: gtf_features.py bed [-h] [--feature {exon,intron,transcript,gene}] gtf

Extract a BED3 file for the selected feature

positional arguments:
  gtf                   path of input GTF file

optional arguments:
  -h, --help            show this help message and exit
  --feature {exon,intron,transcript,gene}
                        feature of interest
```

### `count` module

```
$ python gtf_features.py count --help
usage: gtf_features.py count [-h] [--summary] gtf

Count features per gene (i.e. exon, intron, transcript)

positional arguments:
  gtf         path of input GTF file

optional arguments:
  -h, --help  show this help message and exit
  --summary   print summary statistics
```

### `length` module

```
usage: gtf_features.py length [-h] [--summary] gtf

Extract feature lengths (i.e. exon, intron, transcript, gene)

positional arguments:
  gtf         path of input GTF file

optional arguments:
  -h, --help  show this help message and exit
  --summary   print summary statistics
```
