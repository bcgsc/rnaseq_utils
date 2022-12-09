# gtf_filter.py

A Python script to filter a GTF file based on a list of transcript IDs.

## input files

* GTF annotation file
* file containing a list of transcript IDs (one ID per line)

## usage

```
$ python gtf_filter.py --help
usage: gtf_filter.py [-h] [--fix] gtf tids

Filter a GTF file based on a list of transcript IDs

positional arguments:
  gtf         path of input GTF file
  tids        path of transcript ID list file

optional arguments:
  -h, --help  show this help message and exit
  --fix       fix transcript IDs by removing `.` and trailing characters
```

