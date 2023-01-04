import argparse

parser = argparse.ArgumentParser(description='Extract gene expression from Trans-NanoSim quantification file')
parser.add_argument('tpm', help='path of input TPM file')
parser.add_argument('gtf', help='path of input GTF file')
args = parser.parse_args()

# function to extract gene ID and transcript ID
def get_gid_tid_from_attribute_col(col):
    gid = None
    tid = None
    for info in col.split(';'):
        key, val = info.strip().split()
        if key == 'transcript_id':
            tid = val.strip('"')
        elif key == 'gene_id':
            gid = val.strip('"')
        if gid and tid:
            break
    return (gid, tid)

tid_gid_dict = dict()

# parse GTF and extract transcript-gene ID pairings
with open(args.gtf) as fh:
    for line in fh:
        line = line.strip()
        if len(line.strip()) > 0 and line[0] != '#':
            cols = line.split('\t')
            if cols[2] == 'exon' or cols[2] == 'transcript':                
                gid, tid = get_gid_tid_from_attribute_col(cols[8])
                if gid and tid:
                    tid_gid_dict[tid] = gid

gene_counts = dict()
gene_tpms = dict()

# parse transcript expression file and tally gene expression
with open(args.tpm) as fh:
    fh.readline() # read header line
    for line in fh:
        tid, count, tpm = line.split('\t')
        gid = tid_gid_dict[tid]
        if gid in gene_tpms:
            gene_tpms[gid] += float(tpm)
            gene_counts[gid] += float(count)
        else:
            gene_tpms[gid] = float(tpm)
            gene_counts[gid] = float(count)

# print gene expression values
print('ID', 'count', 'TPM', sep='\t')
for gid in gene_tpms:
    print(gid, gene_counts[gid], gene_tpms[gid], sep='\t')

