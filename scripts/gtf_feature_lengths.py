import argparse
from statistics import mean, median, stdev

# Written by Ka Ming Nip @kmnip

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
    return [num, min_val, q1_val, median_val, q3_val, max_val, mean_val, stdev_val]

def get_tid_gid_from_attribute_col(col):
    tid = None
    gid = None
    for info in col.split(';'):
        key, val = info.strip().split()
        if key == 'transcript_id':
            tid = val.strip('"')
        elif key == 'gene_id':
            gid = val.strip('"')
    return tid, gid

parser = argparse.ArgumentParser(description='Extract the lengths of transcript, gene, exon, and intron.')
parser.add_argument('gtf', help='path of input GTF file')
parser.add_argument('--summary', action='store_true', help='print summary statistics')
args = parser.parse_args()

summarize = args.summary

transcript_lengths = []
gene_lengths = []
exon_lengths = []
intron_lengths = []

with open(args.gtf) as fh:
    if not summarize:
        print('name', 'feature', 'length', sep='\t')

    prev_gid = None
    prev_tid = None
    prev_end = None
    txpt_len = 0
    min_gid_start = float('inf')
    max_gid_end = -float('inf')
    exons = set()
    introns = set()
    
    for line in fh:
        line = line.strip()
        if len(line) > 0 and line[0] != '#':
            cols = line.split('\t')
            if cols[2] == 'exon':
                chrom = cols[0]
                start = int(cols[3])
                end = int(cols[4])
                
                tid, gid = get_tid_gid_from_attribute_col(cols[8].rstrip(';'))
                                
                if tid == prev_tid:
                    txpt_len += end - start
                    introns.add((chrom, prev_end+1, start-1))
                else:
                    if prev_tid:
                        if summarize:
                            transcript_lengths.append(txpt_len)
                        else:
                            print(prev_tid, 'transcript', str(txpt_len), sep='\t')
                    txpt_len = end - start
                    prev_tid = tid      
                
                if gid == prev_gid:
                    min_gid_start = min(min_gid_start, start)
                    max_gid_end = max(max_gid_end, end)
                else:
                    if prev_gid:
                        if summarize:
                            gene_lengths.append(max_gid_end - min_gid_start)
                            for c,s,e in exons:
                                exon_lengths.append(e-s)
                            for c,s,e in introns:
                                intron_lengths.append(e-s)
                        else:
                            print(prev_gid, 'gene', str(max_gid_end - min_gid_start), sep='\t')
                            for c,s,e in sorted(exons):
                                print(c + ':' + str(s) + '-' + str(e), 'exon', str(e-s), sep='\t')
                            for c,s,e in sorted(introns):
                                print(c + ':' + str(s) + '-' + str(e), 'intron', str(e-s), sep='\t')
                        
                    min_gid_start = start
                    max_gid_end = end
                    prev_gid = gid
                    exons = set()
                    introns = set()
                
                exons.add((chrom, start, end))
                prev_end = end
    
    if summarize:
        transcript_lengths.append(txpt_len)
        gene_lengths.append(max_gid_end - min_gid_start)
        for c,s,e in exons:
            exon_lengths.append(e-s)
        for c,s,e in introns:
            intron_lengths.append(e-s)
        
        print('feature', 'n', 'min', 'q1', 'median', 'q3', 'max', 'mean', 'stdev', sep='\t')
        row = ['transcript']
        for val in get_summary(transcript_lengths):
            row.append(str(val))
        print('\t'.join(row))
        row = ['gene']
        for val in get_summary(gene_lengths):
            row.append(str(val))
        print('\t'.join(row))
        row = ['exon']
        for val in get_summary(exon_lengths):
            row.append(str(val))
        print('\t'.join(row))
        row = ['intron']
        for val in get_summary(intron_lengths):
            row.append(str(val))
        print('\t'.join(row))
    else:
        print(prev_tid, 'transcript', str(txpt_len), sep='\t')
        print(prev_gid, 'gene', str(max_gid_end - min_gid_start), sep='\t')
        for c,s,e in sorted(exons):
            print(c + ':' + str(s) + '-' + str(e), 'exon', str(e-s), sep='\t')
        for c,s,e in sorted(introns):
            print(c + ':' + str(s) + '-' + str(e), 'intron', str(e-s), sep='\t')

