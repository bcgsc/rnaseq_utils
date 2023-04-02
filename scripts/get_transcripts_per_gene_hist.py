import argparse

parser = argparse.ArgumentParser(description='Extract the histogram of transcripts per gene')
parser.add_argument('tids', help='path of transcript IDs')
parser.add_argument('gtf', help='path of GTF for matching transcript to gene')
parser.add_argument('--truth', help='path of ground truth transcript IDs for filtering')
args = parser.parse_args()

# function to open both gzip'd and regular files
def gzopen(file_path, mode='rt', compresslevel=6):
    if file_path.lower().endswith('.gz'):
        return gzip.open(file_path, mode=mode, compresslevel=compresslevel)
    return open(file_path, mode)

# fix ENSEMBL gene/transcript names    
def fix_name(name):
    if name.startswith('ENS'):
        return name.split('.')[0]
    return name

def get_gene_map(gtf):
    gene_map = dict()
    with gzopen(gtf) as fh:
        for line in fh:
            line = line.strip()
            if line[0] != '#':
                cols = line.split('\t')
                if cols[2] == 'transcript' or cols[2] == 'exon':
                    gene = None
                    transcript = None
                    for info in cols[8].split(';'):
                        key, val = info.strip().split()
                        if key == 'gene_id':
                            gene = fix_name(val.strip('"'))
                        elif key == 'transcript_id':
                            transcript = fix_name(val.strip('"'))
                        if gene and transcript:
                            gene_map[transcript] = gene
                            break
    return gene_map


gene_map = get_gene_map(args.gtf)

gene_transcript_count_map = dict()
for tid, gid in gene_map.items():
    if gid in gene_transcript_count_map:
        gene_transcript_count_map[gid] += 1
    else:
        gene_transcript_count_map[gid] = 1

# check whether gene ID is a multi-transcript gene
def is_gid_mtg(gid):
    global gene_transcript_count_map
    num_txpt = gene_transcript_count_map[gid]
    assert num_txpt > 0
    return num_txpt > 1

# check whether transcript ID is from a multi-transcript gene
def is_tid_mtg(tid):
    global gene_map, gene_transcript_count_map
    gid = gene_map[tid]
    num_txpt = gene_transcript_count_map[gid]
    assert num_txpt > 0
    return num_txpt > 1

# extract transcript IDs
has_truth = args.truth is not None
truth_tids = set()
tids = list()

if has_truth:
    with gzopen(args.truth) as fh:
        for line in fh:
            truth_tids.add(line.strip())

    with gzopen(args.tids) as fh:
        for line in fh:
            t = line.strip()
            if t in truth_tids:
                tids.append(t)
else:
    with gzopen(args.tids) as fh:
        for line in fh:
            tids.append(line.strip())

# count number of transcripts per gene
# split into two categories
# 1. mtg : multi-transcript genes
# 2. stg : single-transcript genes
stg_counts = dict()
mtg_counts = dict()
for tid in tids:
    gene = gene_map[tid]
    if is_tid_mtg(tid):
        if gene in mtg_counts:
            mtg_counts[gene] += 1
        else:
            mtg_counts[gene] = 1
    else:
        if gene in stg_counts:
            stg_counts[gene] += 1
        else:
            stg_counts[gene] = 1

# extract histograms
stg_hist = dict()
for val in stg_counts.values():
    if val in stg_hist:
       stg_hist[val] += 1
    else:
       stg_hist[val] = 1

mtg_hist = dict()
for val in mtg_counts.values():
    if val in mtg_hist:
       mtg_hist[val] += 1
    else:
       mtg_hist[val] = 1

# print histograms
print('category', 'transcripts_per_gene', 'frequency', sep='\t')
for key,value in sorted(stg_hist.items()):
    print('stg', key, value, sep='\t')

for key,value in sorted(mtg_hist.items()):
    print('mtg', key, value, sep='\t')

