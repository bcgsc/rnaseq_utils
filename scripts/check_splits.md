# check_split.py

A Python script to check split-alignment of contigs for paired read support.

## input

* contigs-to-genome alignment BAM file
* paired-reads-to-genome alignment BAM file

## usage

```
python check_splits.py contigs.bam paired_reads.bam > check_splits.txt

```

## output

```
$ grep qq check_splits.txt | awk '$9>0' | wc -l
30000
$ grep qq check_splits.txt | awk '$9>0' | cut -f2 -d' ' | uniq | wc -l
10000
$ grep qq check_splits.txt | awk '$9>0' | cut -f2 -d' ' | uniq -c | sed 's/^[ ]*//' | awk '$1==1' | wc -l
6000
```

10,000 transcripts with 30,000 splits supported by >= 1 read pair support.
6,000 transcripts with 1 split pair that has >= 1 read pair support.

