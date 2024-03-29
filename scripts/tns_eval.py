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

def evaluate_batch(batch, txpt_recon_props, min_aln_len, min_aln_pid, max_aln_indel,
                  truth_ids, gene_map, full_prop):
    # find the best record
    best_record = None
    best_nmatch = 0
    has_skipped_record = False
    
    for cols in batch:
        nmatch = int(cols[9])
        blen = int(cols[10])
        
        if get_max_indel(get_paf_cigar(cols)) <= max_aln_indel:
            tname = cols[5]
            
            if nmatch > best_nmatch:
                best_record = cols
                best_nmatch = nmatch
            elif nmatch == best_nmatch:
                if tname in truth_ids and best_record[5] not in truth_ids:
                    best_record = cols
        else:
            has_skipped_record = True
        
    if has_skipped_record and not best_record:
        for cols in batch:
            tname = cols[5]
            nmatch = int(cols[9])
            
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
        best_pid = float(best_record[9])/float(best_record[10])
        
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
        
        max_indel = get_max_indel(get_paf_cigar(best_record))
        if max_indel > max_aln_indel:
            # indel too large
            #return ('MISASSEMBLY', qname, best_tname, best_tname, True)
            return ('LARGEINDEL', qname, best_tname, max_indel)
        
        if best_pid < min_aln_pid:
            # percent identity too low
            return ('LOWQUALITY', qname, best_tname, best_pid)
        
        # not a misassembly and not low quality; calculate reconstruction
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
        
        return ('RECONSTRUCTION', qname, best_tname, trp, best_pid)
        
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
parser.add_argument('--aln_pid', dest='aln_pid', default='0.95', metavar='FLOAT', type=float,
                    help='minimum alignment percent identity (default: %(default)s)')
parser.add_argument('--aln_len', dest='aln_len', default='100', metavar='INT', type=int,
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

num_complete_contigs = 0
num_partial_contigs = 0
num_misassembled_contigs = 0
num_false_pos_contigs = 0
num_low_qual_contigs = 0
num_large_indel_contigs = 0
classified_contigs = set()
unclassified_contigs = set()
intragene_misassemblies = list()
intergene_misassemblies = list()

prev_qname = None
with gzopen(args.paf) as fh, \
    open(args.outprefix + 'reconstruction.tsv', 'wt') as fw, \
    open(args.outprefix + 'lowquality.tsv', 'wt') as fw2, \
    open(args.outprefix + 'largeindel.tsv', 'wt') as fw3:
    
    fw.write('contig_id\ttranscript_id\treconstruction\tpercent_identity\n')
    fw2.write('contig_id\ttranscript_id\tpercent_identity\n')
    fw3.write('contig_id\ttranscript_id\tmax_indel\n')
    
    batch = list()
    
    def process_batch():
        result = evaluate_batch(batch, txpt_recon_props, min_aln_len, min_aln_pid,
                     max_aln_indel, truth_ids, gene_map, min_full_prop)
        if result:
            result_type = result[0]
            if result_type == 'MISASSEMBLY':
                classified_contigs.add(prev_qname)
                if result[-1]:
                   intragene_misassemblies.append(result[1:])
                else:
                   intergene_misassemblies.append(result[1:])
                global num_misassembled_contigs
                num_misassembled_contigs += 1
            elif result_type == 'LARGEINDEL':
                classified_contigs.add(prev_qname)
                cid, tid, maxindel = result[1:]
                fw3.write(cid + '\t' + tid + '\t' + str(maxindel) + '\n')
                global num_large_indel_contigs
                num_large_indel_contigs += 1
            elif result_type == 'RECONSTRUCTION':
                classified_contigs.add(prev_qname)
                cid, tid, reconstruction, pid = result[1:]
                fw.write(cid + '\t' + tid + '\t' + str(reconstruction) + '\t' + str(pid) + '\n')
                if tid in truth_ids:
                    # not a false positive
                    if reconstruction >= min_full_prop:
                        # a "complete" reconstruction
                        global num_complete_contigs
                        num_complete_contigs += 1
                    else:
                        # a "partial" reconstruction
                        global num_partial_contigs
                        num_partial_contigs += 1
                else:
                    # a false positive
                    global num_false_pos_contigs
                    num_false_pos_contigs += 1
            elif result_type == 'LOWQUALITY':
                classified_contigs.add(prev_qname)
                cid, tid, pid = result[1:]
                fw2.write(cid + '\t' + tid + '\t' + str(pid) + '\n')
                global num_low_qual_contigs
                num_low_qual_contigs += 1
    
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
            process_batch()
            batch = list()
            
        if blen >= min_aln_len:
            batch.append(cols)
            
        prev_qname = qname

    # process the last read's alignments
    process_batch()

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
            if cid not in classified_contigs:
                unclassified_contigs.add(cid)
        else:
            seq += line.strip()
    if cid:
        # store final seq
        assembly_cid_seq_dict[cid] = seq
num_contigs = len(assembly_cid_seq_dict)

print("total contigs", num_contigs, sep='\t')
print("complete contigs", num_complete_contigs, sep='\t')
print("partial contigs", num_partial_contigs, sep='\t')
print("misassembled contigs", num_misassembled_contigs, sep='\t')
print("false-positive contigs", num_false_pos_contigs, sep='\t')
print("large-indel contigs", num_large_indel_contigs, sep='\t')
print("low-quality contigs", num_low_qual_contigs, sep='\t')

num_unclassified_contigs = num_contigs - num_complete_contigs \
                           - num_partial_contigs - num_misassembled_contigs \
                           - num_false_pos_contigs - num_low_qual_contigs \
                           - num_large_indel_contigs
assert num_unclassified_contigs == len(unclassified_contigs)
print("unclassified contigs", num_unclassified_contigs, sep='\t')

with open(args.outprefix + 'unclassified_contigs.fa', 'wt') as fh:
    for cid in sorted(unclassified_contigs):
        fh.write('>' + cid + '\n' + assembly_cid_seq_dict[cid] + '\n')

# tally all results
complete = list()
partial = list()
missing = list()

# transcripts of single transcript genes
complete_stg = list()
partial_stg = list()
missing_stg = list()

# transcripts of multi-transcript genes
complete_mtg = list()
partial_mtg = list()
missing_mtg = list()

for t in truth_ids:
    is_mtg = is_tid_mtg(t)
    if t in txpt_recon_props:
        p = txpt_recon_props[t]
        assert p <= 1.0
        if p >= min_full_prop:
            complete.append((t, p))
            if is_mtg:
                complete_mtg.append((t, p))
            else:
                complete_stg.append((t, p))
        else:
            partial.append((t, p))
            if is_mtg:
                partial_mtg.append((t, p))
            else:
                partial_stg.append((t, p))
    else:
        missing.append((t, 0))
        if is_mtg:
            missing_mtg.append((t, 0))
        else:
            missing_stg.append((t, 0))

false_pos = list()
false_pos_stg = list()
false_pos_mtg = list()
for fp in txpt_recon_props.keys() - truth_ids:
    false_pos.append((fp, txpt_recon_props[fp]))
    if is_tid_mtg(fp):
        false_pos_mtg.append((fp, txpt_recon_props[fp]))
    else:
        false_pos_stg.append((fp, txpt_recon_props[fp]))

# check results
assert len(complete) == len(complete_stg) + len(complete_mtg)
assert len(partial) == len(partial_stg) + len(partial_mtg)
assert len(missing) == len(missing_stg) + len(missing_mtg)
assert len(false_pos) == len(false_pos_stg) + len(false_pos_mtg)

# print results

print("complete transcripts", len(complete), sep='\t')
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

print("false-positive transcripts", len(false_pos), sep='\t')

num_intragene_mis = len(intragene_misassemblies)
num_intergene_mis = len(intergene_misassemblies)
num_misassemblies = num_intragene_mis + num_intergene_mis
print("intra-gene misassemblies", num_intragene_mis, sep='\t')
print("inter-gene misassemblies", num_intergene_mis, sep='\t')
print("total misassemblies", num_misassemblies, sep='\t')

num_intergene_mis_stg = 0
num_intergene_mis_mtg = 0
for m in intergene_misassemblies:
    tid1 = m[1]
    tid2 = m[2]
    if is_tid_mtg(tid1) or is_tid_mtg(tid2):
        num_intergene_mis_mtg += 1
    else:
        num_intergene_mis_stg += 1
        
num_intragene_mis_stg = 0
num_intragene_mis_mtg = 0
for m in intragene_misassemblies:
    tid1 = m[1]
    tid2 = m[2]
    if is_tid_mtg(tid1) or is_tid_mtg(tid2):
        num_intragene_mis_mtg += 1
    else:
        num_intragene_mis_stg += 1

# check results
assert num_intragene_mis == num_intragene_mis_stg + num_intragene_mis_mtg
assert num_intergene_mis == num_intergene_mis_stg + num_intergene_mis_mtg

# transcripts from single-transcript genes
print("STG complete transcripts", len(complete_stg), sep='\t')
if tpm_bin_map:
    q = 1
    for val in get_tpm_quartile_size(complete_stg, tpm_bin_map):
        print("STG complete transcripts (Q" + str(q) +")", val, sep='\t')
        q += 1

print("STG partial transcripts", len(partial_stg), sep='\t')
if tpm_bin_map:
    q = 1
    for val in get_tpm_quartile_size(partial_stg, tpm_bin_map):
        print("STG partial transcripts (Q" + str(q) +")", val, sep='\t')
        q += 1
    
print("STG missing transcripts", len(missing_stg), sep='\t')
if tpm_bin_map:
    q = 1
    for val in get_tpm_quartile_size(missing_stg, tpm_bin_map):
        print("STG missing transcripts (Q" + str(q) +")", val, sep='\t')
        q += 1

print("STG false-positive transcripts", len(false_pos_stg), sep='\t')
print("STG intra-gene misassemblies", num_intragene_mis_stg, sep='\t')
print("STG inter-gene misassemblies", num_intergene_mis_stg, sep='\t')
print("STG total misassemblies", num_intragene_mis_stg + num_intergene_mis_stg, sep='\t')

# transcripts from multi-transcript genes
print("MTG complete transcripts", len(complete_mtg), sep='\t')
if tpm_bin_map:
    q = 1
    for val in get_tpm_quartile_size(complete_mtg, tpm_bin_map):
        print("MTG complete transcripts (Q" + str(q) +")", val, sep='\t')
        q += 1

print("MTG partial transcripts", len(partial_mtg), sep='\t')
if tpm_bin_map:
    q = 1
    for val in get_tpm_quartile_size(partial_mtg, tpm_bin_map):
        print("MTG partial transcripts (Q" + str(q) +")", val, sep='\t')
        q += 1
    
print("MTG missing transcripts", len(missing_mtg), sep='\t')
if tpm_bin_map:
    q = 1
    for val in get_tpm_quartile_size(missing_mtg, tpm_bin_map):
        print("MTG missing transcripts (Q" + str(q) +")", val, sep='\t')
        q += 1

print("MTG false-positive transcripts", len(false_pos_mtg), sep='\t')
print("MTG intra-gene misassemblies", num_intragene_mis_mtg, sep='\t')
print("MTG inter-gene misassemblies", num_intergene_mis_mtg, sep='\t')
print("MTG total misassemblies", num_intragene_mis_mtg + num_intergene_mis_mtg, sep='\t')

# write results
names = ['complete', 'partial', 'missing', 'false_pos']
lists = [complete, partial, missing, false_pos]

for i in range(0, len(names)):
    n = names[i]
    l = lists[i]
    l.sort(key=lambda tup: tup[1], reverse=True)
    with open(args.outprefix + n + '.tsv', 'wt') as fh:
        fh.write('transcript_id\tmax_reconstruction\n')
        for p in l:
            fh.write(p[0] + '\t' + str(p[1]) + '\n')

with open(args.outprefix + 'intergene_misassemblies.tsv', 'wt') as fh:
    fh.write('contig_id\ttranscript_id1\ttranscript_id2\n')
    for m in intergene_misassemblies:
        fh.write('\t'.join(m[:-1]) + '\n')

with open(args.outprefix + 'intragene_misassemblies.tsv', 'wt') as fh:
    fh.write('contig_id\ttranscript_id1\ttranscript_id2\n')
    for m in intragene_misassemblies:
        fh.write('\t'.join(m[:-1]) + '\n')

with open(args.outprefix + 'redundant.tsv', 'wt') as fh:
    fh.write('transcript_id\tnum_contigs\tcontig_ids\n')
    for ref, names in assigned_txpts.items():
        num_names = len(names)
        if num_names > 1:
            fh.write(ref + '\t' + str(num_names) + '\t' + ' '.join(names) + '\n')

