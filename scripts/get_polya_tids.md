# get_polya_tids.py

A Python script to extract the IDs of polyadenylated reference transcripts based on [PolyASite](https://polyasite.unibas.ch/atlas) BED file.

## input files

* Ensembl annotation GTF file
* PolyASite BED file

## usage

```
$ python get_polya_tids.py --help
usage: get_polya_tids.py [-h] GTF BED

Extract polyA transcript IDs.

positional arguments:
  GTF         annotation GTF file
  BED         polyA site atlas BED file

optional arguments:
  -h, --help  show this help message and exit
```

