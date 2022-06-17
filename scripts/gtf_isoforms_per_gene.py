import argparse
from statistics import mean, median, stdev
from operator import itemgetter

# Written by Ka Ming Nip @kmnip

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

parser = argparse.ArgumentParser(description='Count the number isoforms for each gene in a GTF file')
parser.add_argument('gtf', help='path of input GTF file')
parser.add_argument('--summary', action='store_true', help='print summary statistics')
args = parser.parse_args()

# store a set of `transcript_id` for each `gene_id`
gid_tids_dict = dict()
with open(args.gtf) as fh:
    for line in fh:
        line = line.strip()
        if len(line.strip()) > 0 and line[0] != '#':
            cols = line.split('\t')
            if cols[2] == 'exon' or cols[2] == 'transcript':                
                gid, tid = get_gid_tid_from_attribute_col(cols[8])
                if gid and tid:
                    if gid in gid_tids_dict:
                        gid_tids_dict[gid].add(tid)
                    else:
                        tids = set()
                        tids.add(tid)
                        gid_tids_dict[gid] = tids

# count the number of isoforms for each gene
gid_counts = list((gid, len(tids)) for gid, tids in gid_tids_dict.items())

if args.summary:
    counts = list(c for gid, c in gid_counts)
    counts.sort()
    
    num_genes = len(counts)
    num_isoforms = sum(counts)
    
    min_count = counts[0]
    q1_count = counts[int(num_genes/4)]
    median_count = median(counts)
    q3_count = counts[int(num_genes*3/4)]
    max_count = counts[-1]
    
    mean_count = mean(counts)
    stdev_count = stdev(counts) 
    
    print('genes:   ', str(num_genes), sep='\t')
    print('isoforms:', str(num_isoforms), sep='\t')
    print('  mean:  ', str(mean_count), sep='\t')
    print('  stdev: ', str(stdev_count), sep='\t')
    print('  min:   ', str(min_count), sep='\t')
    print('  Q1:    ', str(q1_count), sep='\t')
    print('  median:', str(median_count), sep='\t')
    print('  Q3:    ', str(q3_count), sep='\t')
    print('  max:   ', str(max_count), sep='\t')
else:
    # reverse-sort by the number of isoforms
    gid_counts.sort(key=itemgetter(1),reverse=True)
    
    # print `gene_id` and its isoform count     
    for d in gid_counts:
        print(d[0], d[1], sep='\t')

