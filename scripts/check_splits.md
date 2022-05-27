# check_split.py

A Python script to check split-alignment of contigs for paired read support.

## input files

* contigs-to-genome alignment BAM file
* paired-reads-to-genome alignment BAM file

## example usage

```
python check_splits.py contigs.bam paired_reads.bam > check_splits.txt
```

## interpreting results

```
$ grep qq check_splits.txt | awk '$9>0' | wc -l
30000
$ grep qq check_splits.txt | awk '$9>0' | cut -f2 -d' ' | uniq | wc -l
10000
$ grep qq check_splits.txt | awk '$9>0' | cut -f2 -d' ' | uniq -c | sed 's/^[ ]*//' | awk '$1==1' | wc -l
6000
```

30,000 splits within 10,000 contigs have >= 1 read pair support each.

6,000 contigs with only 1 split each, where each split has >= 1 read pair support.

