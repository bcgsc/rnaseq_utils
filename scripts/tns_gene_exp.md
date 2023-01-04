# tns_gene_exp.py

A Python script for extracting gene expression from Trans-NanoSim quantification file.

### input arguments

* TSV file of Trans-NanoSim transcript expression levels
  * 3 columns: `target_id`, `est_counts`, `tpm`
* GTF file of reference annotation
  * Must contain attributes for `gene_id` and `transcript_id`

### usage

```
usage: tns_gene_exp.py [-h] tpm gtf

Extract gene expression from Trans-NanoSim quantification file

positional arguments:
  tpm         path of input TPM file
  gtf         path of input GTF file

optional arguments:
  -h, --help  show this help message and exit
```

