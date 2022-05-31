# tns_eval.py

A Python script for evaluating transcriptome assembly quality.

### input files

* PAF file of assembly alignments against the reference transcriptome (**not genome**)
* text file of grouth truth transcript IDs (one ID per line)
  * use `tns_get_tids.sh` to create this file
* GTF file of reference annotation
* path prefix of output files
* TSV file of Trans-NanoSim transcript expression levels (optional)
  * 3 columns: `target_id`, `est_counts`, `tpm`

### usage

```
$ python tns_eval.py --help
usage: tns_eval.py [-h] [--full_prop FLOAT] [--aln_pid FLOAT] [--aln_len INT] [--aln_indel INT] [--tpm TSV] paf truth gtf outprefix

Evaluate transcriptome assembly quality

positional arguments:
  paf                path of input PAF file
  truth              path of truth transcript IDs
  gtf                path of gtf
  outprefix          path of output prefix

optional arguments:
  -h, --help         show this help message and exit
  --full_prop FLOAT  minimum length proportion for full-length transcripts (default: 0.95)
  --aln_pid FLOAT    minimum alignment percent identity (default: 0.9)
  --aln_len INT      minimum alignment length (default: 150)
  --aln_indel INT    maximum alignment indel (default: 70)
  --tpm TSV          path of transcript expression TSV
```

### example usage

```
# align the assembly
minimap2 -x map-ont -c -t 12 reference_transcripts.fasta assembly.fasta | gzip -c > aln.paf.gz

# evaluate the assembly
python tns_eval.py aln.paf.gz truth.txt annotation.gtf ./results_ --tpm transnanosim_quant.tsv > ./results_summary.txt
```
