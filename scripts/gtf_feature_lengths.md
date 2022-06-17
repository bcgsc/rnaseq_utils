# gtf_feature_lengths

A Python script for extracting the lengths of transcript, gene, exon, and introns from a GTF file.

### input file

* GTF file containing `gene_id` and `transcript_id` attributes

### usage

```
usage: gtf_feature_lengths.py [-h] [--summary] gtf

Extract the lengths of transcript, gene, exon, and intron.

positional arguments:
  gtf         path of input GTF file

optional arguments:
  -h, --help  show this help message and exit
  --summary   print summary statistics
```

