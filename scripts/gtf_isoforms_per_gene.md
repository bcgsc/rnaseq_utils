# gtf_isoforms_per_gene.py

A Python script for counting the number of isoforms for each gene in a GTF file.

### input file

* GTF file containing `gene_id` and `transcript_id` attributes

### usage

```
usage: gtf_isoforms_per_gene.py [-h] [--summary] gtf

Count the number isoforms for each gene in a GTF file

positional arguments:
  gtf         path of input GTF file

optional arguments:
  -h, --help  show this help message and exit
  --summary   print summary statistics
```

### example usage

```
$ python gtf_isoforms_per_gene.py Homo_sapiens.GRCh38.103.gtf
ENSG00000179818	239
ENSG00000109339	192
ENSG00000226674	189
ENSG00000249859	176
ENSG00000241469	164
ENSG00000227195	146
ENSG00000215386	143
ENSG00000242086	142
ENSG00000229140	142
ENSG00000224078	142
...
```

```
$ python gtf_isoforms_per_gene.py Homo_sapiens.GRCh38.103.gtf --summary
genes:   	60666
isoforms:	234393
  mean:  	3.8636633369597466
  stdev: 	6.733126550886392
  min:   	1
  Q1:    	1
  median:	1.0
  Q3:    	4
  max:   	239
```
