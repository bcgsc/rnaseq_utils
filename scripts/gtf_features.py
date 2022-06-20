import argparse
from statistics import mean, median, stdev

# Written by Ka Ming Nip @kmnip

class FEATURE(str):
    EXON = 'exon'
    INTRON = 'intron'
    TRANSCRIPT = 'transcript'
    GENE = 'gene'

class MODE(str):
    BED = 'bed'
    COUNT = 'count'
    LENGTH = 'length'

def get_summary(lst):
    lst.sort()
    num = len(lst)
    min_val = lst[0]
    q1_val = lst[int(num/4)]
    median_val = median(lst)
    q3_val = lst[int(num*3/4)]
    max_val = lst[-1]
    mean_val = mean(lst)
    stdev_val = stdev(lst)
    return [min_val, q1_val, median_val, q3_val, max_val, mean_val, stdev_val]

def get_tid_gid_from_attribute_col(col):
    tid = None
    gid = None
    for info in col.rstrip(';').split(';'):
        key, val = info.strip().split(' ', 1)
        if key == 'transcript_id':
            tid = val.strip('"')
        elif key == 'gene_id':
            gid = val.strip('"')
    return tid, gid

def exon_generator(gtf):
    with open(gtf) as fh:
        for line in fh:
            line = line.strip()
            if len(line) > 0 and line[0] != '#':
                cols = line.split('\t')
                if cols[2] == FEATURE.EXON:
                    #print(line)
                    chrom = cols[0]
                    start = int(cols[3])
                    end = int(cols[4])
                    strand = cols[6]
                    tid, gid = get_tid_gid_from_attribute_col(cols[8].rstrip(';'))
                    yield (chrom, start, end, strand, tid, gid)

def extract_introns(exon_chain):
    introns = list()
    num_exons = len(exon_chain)
    if num_exons >= 2:
        exon_chain.sort()
        #print(exon_chain)
        for i in range(1, num_exons):
            chrom1, start1, end1, strand1, gid1 = exon_chain[i-1]
            chrom2, start2, end2, strand2, gid2 = exon_chain[i]
            assert chrom1 == chrom2
            assert strand1 == strand2
            assert gid1 == gid2
            introns.append((chrom1, end1+1, start2-1, strand1, gid1))
    return introns

def extract_feature_lengths(gtf, summarize=False):
    transcript_lengths = []
    gene_lengths = []
    exon_lengths = []
    intron_lengths = []

    if not summarize:
        print('name', 'feature', 'length', sep='\t')
    
    prev_gid = None
    prev_tid = None
    txpt_len = 0
    min_gid_start = float('inf')
    max_gid_end = -float('inf')
    exons = set()
    introns = set()
    exon_chain = list()
                        
    for chrom, start, end, strand, tid, gid in exon_generator(gtf):
        if tid == prev_tid:
            txpt_len += end - start + 1
        else:
            if prev_tid:
                if summarize:
                    transcript_lengths.append(txpt_len)
                else:
                    print(prev_tid, FEATURE.TRANSCRIPT, str(txpt_len), sep='\t')
            introns.update(extract_introns(exon_chain))
            exon_chain = list()
            txpt_len = end - start + 1
        
        if gid == prev_gid:
            min_gid_start = min(min_gid_start, start)
            max_gid_end = max(max_gid_end, end)
        else:
            if prev_gid:
                if summarize:
                    gene_lengths.append(max_gid_end - min_gid_start + 1)
                    for _chrom, _start, _end, _strand, _gid in exons:
                        exon_lengths.append(_end - _start + 1)
                    for _chrom, _start, _end, _strand, _gid in introns:
                        intron_lengths.append(_end - _start + 1)
                else:
                    print(prev_gid, FEATURE.GENE, str(max_gid_end - min_gid_start + 1), sep='\t')
                    for _chrom, _start, _end, _strand, _gid in sorted(exons):
                        print(_chrom + ':' + str(_start) + '-' + str(_end) + ':' + _strand, FEATURE.EXON, str(_end - _start + 1), sep='\t')
                    for _chrom, _start, _end, _strand, _gid in sorted(introns):
                        print(_chrom + ':' + str(_start) + '-' + str(_end)+ ':' + _strand, FEATURE.INTRON, str(_end - _start + 1), sep='\t')
                
            min_gid_start = start
            max_gid_end = end
            prev_gid = gid
            exons = set()
            introns = set()
        
        exons.add((chrom, start, end, strand, gid))
        exon_chain.append((chrom, start, end, strand, gid))
        prev_end = end
        prev_tid = tid
    
    
    introns.update(extract_introns(exon_chain))
    if summarize:
        # process final record
        transcript_lengths.append(txpt_len)
        gene_lengths.append(max_gid_end - min_gid_start + 1)
        for _chrom, _start, _end, _strand, _gid in exons:
            exon_lengths.append(_end - _start + 1)
        for _chrom, _start, _end, _strand, _gid in introns:
            intron_lengths.append(_end - _start + 1)
        
        # print summary
        print('feature', 'n', 'min', 'q1', 'median', 'q3', 'max', 'mean', 'stdev', sep='\t')
        
        row = [FEATURE.TRANSCRIPT]
        row.append(str(len(transcript_lengths)))
        for val in get_summary(transcript_lengths):
            row.append(str(val))
        print('\t'.join(row))
        
        row = [FEATURE.GENE]
        row.append(str(len(gene_lengths)))
        for val in get_summary(gene_lengths):
            row.append(str(val))
        print('\t'.join(row))
        
        row = [FEATURE.EXON]
        row.append(str(len(exon_lengths)))
        for val in get_summary(exon_lengths):
            row.append(str(val))
        print('\t'.join(row))
        
        row = [FEATURE.INTRON]
        row.append(str(len(intron_lengths)))
        for val in get_summary(intron_lengths):
            row.append(str(val))
        print('\t'.join(row))
    else:
        # print final record
        if prev_tid:
            print(prev_tid, FEATURE.TRANSCRIPT, str(txpt_len), sep='\t')
        if prev_gid:
            print(prev_gid, FEATURE.GENE, str(max_gid_end - min_gid_start + 1), sep='\t')
            for _chrom, _start, _end, _strand, _gid in sorted(exons):
                print(_chrom + ':' + str(_start) + '-' + str(_end) + ':' + _strand, FEATURE.EXON, str(_end - _start + 1), sep='\t')
            for _chrom, _start, _end, _strand, _gid in sorted(introns):
                print(_chrom + ':' + str(_start) + '-' + str(_end) + ':' + _strand, FEATURE.INTRON, str(_end - _start + 1), sep='\t')

def count_features_per_gene(gtf, summarize=False):
    transcript_counts = []
    exon_counts = []
    intron_counts = []
    
    prev_gid = None
    prev_tid = None
    exons = set()
    introns = set()
    transcripts = set()
    exon_chain = list()
    
    if not summarize:
        print('gene_name', 'feature', 'count_per_gene', sep='\t')
    
    for chrom, start, end, strand, tid, gid in exon_generator(gtf):
        if tid != prev_tid:
            introns.update(extract_introns(exon_chain))
            exon_chain = list()
        if gid != prev_gid:
            if prev_gid:
                if summarize:
                    transcript_counts.append(len(transcripts))
                    exon_counts.append(len(exons))
                    intron_counts.append(len(introns))
                else:
                    print(gid, FEATURE.TRANSCRIPT, str(len(transcripts)))
                    print(gid, FEATURE.EXON, str(len(exons)))
                    print(gid, FEATURE.INTRON, str(len(introns)))
            
            transcripts = set()
            exons = set()
            introns = set()
        
        exons.add((chrom, start, end, strand, gid))
        exon_chain.append((chrom, start, end, strand, gid))
        transcripts.add(tid)
        prev_tid = tid
        prev_gid = gid  
    
    introns.update(extract_introns(exon_chain))
                    
    if prev_gid:
        # process final record
        if summarize:
            transcript_counts.append(len(transcripts))
            exon_counts.append(len(exons))
            intron_counts.append(len(introns))
        else:
            print(gid, FEATURE.TRANSCRIPT, str(len(transcripts)))
            print(gid, FEATURE.EXON, str(len(exons)))
            print(gid, FEATURE.INTRON, str(len(introns)))

    if summarize:
        # print summary
        print('feature', 'n', 'min', 'q1', 'median', 'q3', 'max', 'mean', 'stdev', sep='\t')
        row = [FEATURE.TRANSCRIPT]
        row.append(str(sum(transcript_counts)))
        for val in get_summary(transcript_counts):
            row.append(str(val))
        print('\t'.join(row))
        
        row = [FEATURE.EXON]
        row.append(str(sum(exon_counts)))
        for val in get_summary(exon_counts):
            row.append(str(val))
        print('\t'.join(row))
        
        row = [FEATURE.INTRON]
        row.append(str(sum(intron_counts)))
        for val in get_summary(intron_counts):
            row.append(str(val))
        print('\t'.join(row))

def extract_bed3(gtf, feature):
    if feature == FEATURE.EXON:
        exons = set()
        for chrom, start, end, strand, tid, gid in exon_generator(gtf):
            exons.add((chrom, str(start), str(end), strand, gid))
        for chrom, start, end, strand, gid in sorted(exons):
            print(chrom, str(start), str(end), sep='\t')
            
    elif feature == FEATURE.INTRON:
        prev_tid = None
        introns = set()
        exon_chain = list()
        for chrom, start, end, strand, tid, gid in exon_generator(gtf):
            if tid != prev_tid:
                introns.update(extract_introns(exon_chain))
                exon_chain = list()
            prev_tid = tid
            exon_chain.append((chrom, start, end, strand, gid))
        introns.update(extract_introns(exon_chain))
        for chrom, start, end, strand, gid in sorted(introns):
            print(chrom, str(start), str(end), sep='\t')
            
    elif feature == FEATURE.TRANSCRIPT:
        transcript_itv = dict()
        for chrom, start, end, strand, tid, gid in exon_generator(gtf):
            if tid in transcript_itv:
                t_chrom, t_start, t_end = transcript_itv[tid]
                assert t_chrom == chrom
                transcript_itv[tid] = (chrom, min(start, t_start), max(end, t_end))
            else:
                transcript_itv[tid] = (chrom, start, end)
        for chrom, start, end in transcript_itv.values():
            print(chrom, str(start), str(end), sep='\t')
            
    elif feature == FEATURE.GENE:
        gene_itv = dict()
        for chrom, start, end, strand, tid, gid in exon_generator(gtf):
            if gid in gene_itv:
                g_chrom, g_start, g_end = gene_itv[gid]
                assert g_chrom == chrom
                gene_itv[gid] = (chrom, min(start, g_start), max(end, g_end))
            else:
                gene_itv[gid] = (chrom, start, end)
        for chrom, start, end in gene_itv.values():
            print(chrom, str(start), str(end), sep='\t')
            
    else:
        raise ValueError('Unknown feature value ' + feature)

parser = argparse.ArgumentParser(description='Extract feature information from GTF file')
subparsers = parser.add_subparsers(dest='mode')

parser_bed_help = "Extract a BED3 file for the selected feature"
parser_bed = subparsers.add_parser(MODE.BED,
    description=parser_bed_help, help=parser_bed_help)
parser_bed.add_argument('gtf', help='path of input GTF file')
parser_bed.add_argument('--feature', help='feature of interest',
    required=True,
    choices=[FEATURE.EXON, FEATURE.INTRON, FEATURE.TRANSCRIPT, FEATURE.GENE])

parser_count_help = "Count features per gene (i.e. exon, intron, transcript)"
parser_count = subparsers.add_parser(MODE.COUNT,
    description=parser_count_help, help=parser_count_help)
parser_count.add_argument('gtf', help='path of input GTF file')
parser_count.add_argument('--summary', action='store_true', help='print summary statistics')

parser_length_help = "Extract feature lengths (i.e. exon, intron, transcript, gene)"
parser_length = subparsers.add_parser(MODE.LENGTH,
    description=parser_length_help, help=parser_length_help)
parser_length.add_argument('gtf', help='path of input GTF file')
parser_length.add_argument('--summary', action='store_true', help='print summary statistics')

args = parser.parse_args()

if args.mode == MODE.BED:
    extract_bed3(args.gtf, args.feature)
elif args.mode == MODE.LENGTH:
    extract_feature_lengths(args.gtf, args.summary)
elif args.mode == MODE.COUNT:
    count_features_per_gene(args.gtf, args.summary)

