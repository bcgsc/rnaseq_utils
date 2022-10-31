import argparse
import gzip
import logging
import re

# writtern by Ka Ming Nip @kmnip

num_redundant = 0
assigned_txpts = dict()

# function to open both gzip'd and regular files
def gzopen(file_path, mode='rt', compresslevel=6):
    if file_path.lower().endswith('.gz'):
        return gzip.open(file_path, mode=mode, compresslevel=compresslevel)
    return open(file_path, mode)

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

def get_tpm_bin_map(tsv, truth_set):
    tpm_map = dict()
    tpms = list()
    with gzopen(tsv) as fh:
        # read header
        headers = fh.readline().strip().split('\t')
        tid_col = headers.index('ID')
        tpm_col = headers.index('TPM')
        
        # read rest of file
        for line in fh:
            cols = line.strip().split('\t')
            tid = fix_name(cols[tid_col])
            if tid in truth_set:
                tpm = float(cols[tpm_col])
                tpm_map[tid] = tpm
                tpms.append(tpm)
    
    tpms.sort()
    length = len(tpms)
    q1 = tpms[int(length / 4)]
    m = tpms[int(length / 2)]
    q3 = tpms[int(length * 3 / 4)]
    
    tpm_bin_map = dict()
    for tid, tpm in tpm_map.items():
        if tpm <= q1:
            tpm_bin_map[tid] = 0
        elif tpm <= m:
            tpm_bin_map[tid] = 1
        elif tpm <= q3:
            tpm_bin_map[tid] = 2
        else:
            tpm_bin_map[tid] = 3
    
    return tpm_bin_map, (round(tpms[0], 2), round(q1, 2), round(m, 2), round(q3, 2), round(tpms[-1], 2))

def get_tpm_quartile_size(tid_recon_list, tpm_bin_map):
    bin_sizes = [0, 0, 0, 0]
    for t in tid_recon_list:
        b = tpm_bin_map[t[0]]
        assert b is not None
        bin_sizes[b] += 1
    return bin_sizes

# fix ENSEMBL gene/transcript names    
def fix_name(name):
    if name.startswith('ENS'):
        return name.split('.')[0]
    return name

def get_combined_length(start1, end1, start2, end2):
    if start2 <= end1 and end2 >= end1:
        # dovetail: 1 -> 2
        return end2 - start1
    elif start1 <= end2 and end1 >= end2:
        # dovetail: 2 -> 1
        return end1 - start2
    elif start2 > end1 or end2 < start1:
        # no overlap: 1,2 or 2,1
        return (end1 - start1) + (end2 - start2)
    elif start2 >= start1 and end2 <= end1:
        # contained: 2 in 1
        return end1 - start1
    elif start1 >= start2 and end1 <= start1:
        # contained: 1 in 2
        return end2 - start2
    return None

def get_paf_cigar(cols):
    for i in range(12, len(cols)):
        tag, tmp, val = cols[i].split(':')
        if tag == 'cg':
            return val
    return None

def get_max_indel(cigar):
    max_indel = 0
    if cigar:
        for op in re.findall(r'(\d+)[DI]', cigar):
            max_indel = max(max_indel, int(op))
    return max_indel

def process_batch(batch, txpt_recon_props, min_aln_len, truth_ids, gene_map, full_prop):
    # find the best record
    best_record = None
    best_nmatch = 0
    
    for cols in batch:
        tname = cols[5]
        cols[5] = tname
        tlen = int(cols[6])
        nmatch = int(cols[9])
        blen = int(cols[10])
        
        if nmatch > best_nmatch:
            best_record = cols
            best_nmatch = nmatch
        elif nmatch == best_nmatch:
            if tname in truth_ids and best_record[5] not in truth_ids:
                best_record = cols
    
    if best_record:
        qname = best_record[0]
        best_tname = best_record[5]
        best_tlen = int(best_record[6])
        best_tstart = int(best_record[7])
        best_tend = int(best_record[8])
        best_gene = gene_map[best_tname]
        
        if len(batch) > 1:
            # attempt to find misassembly
            best_qlen = int(best_record[1])
            best_qstart = int(best_record[2])
            best_qend = int(best_record[3])
            alt_best_record = None
            # the alt record must have at least `min_aln_len` nucleotides that are not already covered by the best record
            min_len = best_qend - best_qstart + min_aln_len
            merged_length = best_qend - best_qstart
            for r in batch:
                if r is not best_record:
                    qstart = int(r[2])
                    qend = int(r[3])
                    m = get_combined_length(best_qstart, best_qend, qstart, qend)
                    if m >= min_len and m > merged_length:
                        alt_best_record = r
                        merged_length = m
            if alt_best_record:
                alt_best_tname = alt_best_record[5]
                return ('MISASSEMBLY', qname, best_tname, alt_best_tname, gene_map[alt_best_tname] == best_gene)
        
        # not a misassembly; calculate reconstruction
        trp = float(best_tend - best_tstart)/best_tlen
        assert trp <= 1.0
        if best_tname not in txpt_recon_props:
            txpt_recon_props[best_tname] = trp
        else:
            best_trp = txpt_recon_props[best_tname]
            if trp > best_trp:
                txpt_recon_props[best_tname] = trp
            if trp >= full_prop and best_trp >= full_prop:
                global num_redundant
                num_redundant += 1
        
        if trp >= full_prop:
            if best_tname in assigned_txpts:
                cids = assigned_txpts[best_tname]
            else:
                cids = []
                assigned_txpts[best_tname] = cids
            cids.append(qname)
        
        return ('RECONSTRUCTION', qname, best_tname, trp)
        
    return None
        
        

parser = argparse.ArgumentParser(description='Evaluate transcriptome assembly quality')
parser.add_argument('assembly',
                    help='path of assembly FASTA file')
parser.add_argument('paf',
                    help='path of input PAF file')
parser.add_argument('truth',
                    help='path of truth transcript IDs')
parser.add_argument('gtf',
                    help='path of GTF')
parser.add_argument('outprefix',
                    help='path of output prefix')
parser.add_argument('--full_prop', dest='full_prop', default='0.95', metavar='FLOAT', type=float,
                    help='minimum length proportion for full-length transcripts (default: %(default)s)')
parser.add_argument('--aln_pid', dest='aln_pid', default='0.9', metavar='FLOAT', type=float,
                    help='minimum alignment percent identity (default: %(default)s)')
parser.add_argument('--aln_len', dest='aln_len', default='150', metavar='INT', type=int,
                    help='minimum alignment length (default: %(default)s)')
parser.add_argument('--aln_indel', dest='aln_indel', default='70', metavar='INT', type=int,
                    help='maximum alignment indel (default: %(default)s)')
parser.add_argument('--tpm', dest='tpm', metavar='TSV', type=str,
                    help='path of transcript expression TSV')
args = parser.parse_args()

logging.basicConfig(
    format='%(asctime)s %(levelname)-8s %(message)s',
    level=logging.INFO,
    datefmt='%Y-%m-%d %H:%M:%S')

min_aln_pid = args.aln_pid
min_aln_len = args.aln_len
min_full_prop = args.full_prop
max_aln_indel = args.aln_indel

truth_ids = set()
logging.info('parsing truth file...')
with gzopen(args.truth) as fh:
    for line in fh:
        truth_ids.add(line.strip())

tpm_bin_map = None
tpm_quantiles = None
if args.tpm:
    logging.info('parsing abundance file...')
    tpm_bin_map, tpm_quantiles = get_tpm_bin_map(args.tpm, truth_ids)
    logging.info('TPM quantiles:')
    logging.info('min\tq1\tM\tq3\tmax')
    logging.info(str(tpm_quantiles[0]) +
        '\t' + str(tpm_quantiles[1]) +
        '\t' + str(tpm_quantiles[2]) +
        '\t' + str(tpm_quantiles[3]) +
        '\t' + str(tpm_quantiles[4]))

logging.info('parsing GTF file...')
gene_map = get_gene_map(args.gtf)

txpt_lengths = dict()
txpt_recon_props = dict()
logging.info('parsing PAF file...')
"""
1	string	Query sequence name
2	int	Query sequence length
3	int	Query start (0-based; BED-like; closed)
4	int	Query end (0-based; BED-like; open)
5	char	Relative strand: "+" or "-"
6	string	Target sequence name
7	int	Target sequence length
8	int	Target start on original strand (0-based)
9	int	Target end on original strand (0-based)
10	int	Number of residue matches
11	int	Alignment block length
12	int	Mapping quality (0-255; 255 for missing)
"""

prev_qname = None
intragene_misassemblies = list()
intergene_misassemblies = list()
num_complete_contigs = 0
num_partial_contigs = 0
num_misassembled_contigs = 0
num_false_pos_contigs = 0

with open(args.outprefix + 'reconstruction.tsv', 'wt') as fw:
    with gzopen(args.paf) as fh:
        batch = list()
        for line in fh:
            cols = line.strip().split('\t')
            
            qname = cols[0]
            tname = fix_name(cols[5])
            cols[5] = tname
            tlen = int(cols[6])
            nmatch = int(cols[9])
            blen = int(cols[10])
            
            txpt_lengths[tname] = tlen

            
            if prev_qname and prev_qname != qname and len(batch) > 0:
                result = process_batch(batch, txpt_recon_props, min_aln_len, truth_ids, gene_map, min_full_prop)
                if result:
                    result_type = result[0]
                    if result_type == 'MISASSEMBLY':
                        if result[-1]:
                           intragene_misassemblies.append(result[1:])
                        else:
                           intergene_misassemblies.append(result[1:])
                        num_misassembled_contigs += 1
                    elif result_type == 'RECONSTRUCTION':
                        rtype, cid, tid, reconstruction = result
                        fw.write(cid + '\t' + tid + '\t' + str(reconstruction) + '\n')
                        if tid in truth_ids:
                            if reconstruction >= min_full_prop:
                                num_complete_contigs += 1
                            else:
                                num_partial_contigs += 1
                        else:
                            num_false_pos_contigs += 1
                batch = list()
                
            if blen >= min_aln_len and \
                float(nmatch)/blen >= min_aln_pid and \
                get_max_indel(get_paf_cigar(cols)) <= max_aln_indel:
                batch.append(cols)
                
            prev_qname = qname

        # process the last read's alignments
        result = process_batch(batch, txpt_recon_props, min_aln_len, truth_ids, gene_map, min_full_prop)
        if result:
            result_type = result[0]
            if result_type == 'MISASSEMBLY':
                if result[-1]:
                   intragene_misassemblies.append(result[1:])
                else:
                   intergene_misassemblies.append(result[1:])
                num_misassembled_contigs += 1
            elif result_type == 'RECONSTRUCTION':
                rtype, cid, tid, reconstruction = result
                fw.write(cid + '\t' + tid + '\t' + str(reconstruction) + '\n')
                if tid in truth_ids:
                    if reconstruction >= min_full_prop:
                        num_complete_contigs += 1
                    else:
                        num_partial_contigs += 1
                else:
                    num_false_pos_contigs += 1

# parse assembly FASTA
logging.info('parsing assembly file...')
assembly_cid_seq_dict = dict()
with gzopen(args.assembly) as fh:
    cid = None
    seq = ''
    for line in fh:
        if line[0] == '>':
            if cid:
                # store previous seq
                assembly_cid_seq_dict[cid] = seq
            cid = line[1:].strip().split(' ', 1)[0]
            seq = ''
        else:
            seq += line.strip()
    if cid:
        # store final seq
        assembly_cid_seq_dict[cid] = seq
num_contigs = len(assembly_cid_seq_dict)

print("contigs", num_contigs, sep='\t')
print("complete contigs", num_complete_contigs, sep='\t')
print("partial contigs", num_partial_contigs, sep='\t')
print("misassembled contigs", num_misassembled_contigs, sep='\t')
print("false pos. contigs", num_false_pos_contigs, sep='\t')
print("unclassified contigs", num_contigs - num_complete_contigs - num_partial_contigs - num_misassembled_contigs - num_false_pos_contigs, sep='\t')

# tally results
complete = list()
partial = list()
missing = list()
for t in truth_ids:
    if t in txpt_recon_props:
        p = txpt_recon_props[t]
        assert p <= 1.0
        if p >= min_full_prop:
            complete.append((t, p))
        else:
            partial.append((t, p))
    else:
        missing.append((t, 0))

false_pos = list()
for fp in txpt_recon_props.keys() - truth_ids:
    false_pos.append((fp, txpt_recon_props[fp]))

# print results
num_complete = len(complete)
print("complete transcripts", num_complete, sep='\t')
if tpm_bin_map:
    q = 1
    for val in get_tpm_quartile_size(complete, tpm_bin_map):
        print("complete transcripts (Q" + str(q) +")", val, sep='\t')
        q += 1

print("partial transcripts", len(partial), sep='\t')
if tpm_bin_map:
    q = 1
    for val in get_tpm_quartile_size(partial, tpm_bin_map):
        print("partial transcripts (Q" + str(q) +")", val, sep='\t')
        q += 1
    
print("missing transcripts", len(missing), sep='\t')
if tpm_bin_map:
    q = 1
    for val in get_tpm_quartile_size(missing, tpm_bin_map):
        print("missing transcripts (Q" + str(q) +")", val, sep='\t')
        q += 1

print("false pos. transcripts", len(false_pos), sep='\t')

num_intragene_mis = len(intragene_misassemblies)
num_intergene_mis = len(intergene_misassemblies)
num_misassemblies = num_intragene_mis + num_intergene_mis
print("intra-gene misassemblies", num_intragene_mis, sep='\t')
print("inter-gene misassemblies", num_intergene_mis, sep='\t')
print("misassemblies", num_misassemblies, sep='\t')
print("misassemblies / complete", float(num_misassemblies)/num_complete, sep='\t')

# write results
names = ['complete', 'partial', 'missing', 'false_pos']
lists = [complete, partial, missing, false_pos]

for i in range(0, len(names)):
    n = names[i]
    l = lists[i]
    l.sort(key=lambda tup: tup[1], reverse=True)
    with open(args.outprefix + n + '.tsv', 'wt') as fh:
        for p in l:
            fh.write(p[0] + '\t' + str(p[1]) + '\n')

with open(args.outprefix + 'intergene_misassemblies.tsv', 'wt') as fh:
    for m in intergene_misassemblies:
        fh.write('\t'.join(m[:-1]) + '\n')

with open(args.outprefix + 'intragene_misassemblies.tsv', 'wt') as fh:
    for m in intragene_misassemblies:
        fh.write('\t'.join(m[:-1]) + '\n')

with open(args.outprefix + 'redundant.tsv', 'wt') as fh:
    for ref, names in assigned_txpts.items():
        num_names = len(names)
        if num_names > 1:
            fh.write(ref + '\t' + str(num_names) + '\t' + ' '.join(names) + '\n')

