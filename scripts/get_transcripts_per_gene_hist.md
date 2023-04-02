# get_transcripts_per_gene_hist.py

A Python script for extracting the histogram of transcripts per gene for single-transcript genes (stg) and multi-transcript genes (mtg).

### input files

* a text file containing a list of transcript IDs to be evaluated
* a GTF file
* a text file containing a list of grouth truth transcript IDs for filtering the list of transcripts IDs 

### usage

```
usage: get_transcripts_per_gene_hist.py [-h] [--truth TRUTH] tids gtf

Extract the histogram of transcripts per gene

positional arguments:
  tids           path of transcript IDs
  gtf            path of GTF for matching transcript to gene

optional arguments:
  -h, --help     show this help message and exit
  --truth TRUTH  path of ground truth transcript IDs for filtering
```

### example usage

```
$ python get_transcripts_per_gene_hist.py rnabloom.txt mouse_mm10_ensembl.gtf --truth mouse_sim_2M.txt
category	transcripts_per_gene	frequency
stg	1	980
stg	2	52
stg	3	13
...
mtg	1	5681
mtg	2	1686
mtg	3	430
...
```
