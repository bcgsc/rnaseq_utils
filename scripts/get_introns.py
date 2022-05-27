import argparse

# Written by Ka Ming Nip @kmnip

def get_tid_from_attribute_col(col):
    for info in col.split(';'):
        key, val = info.strip().split()
        if key == 'transcript_id':
            return val.strip('"')
    return None

parser = argparse.ArgumentParser(description='Extract a 3-column BED file for introns from an input GTF')
parser.add_argument('gtf', help='path of input GTF file')
args = parser.parse_args()

intron_set = set()
with open(args.gtf) as fh:
    prev_tid = None
    prev_end = None
    for line in fh:
        line = line.strip()
        if len(line.strip()) > 0 and line[0] != '#':
            cols = line.split('\t')
            if cols[2] == 'exon':
                chrom = cols[0]
                start = cols[3]
                end = cols[4]
                
                tid = get_tid_from_attribute_col(cols[8])
                if prev_tid and tid and tid == prev_tid:
                    intron_set.add((chrom, int(prev_end)+1, int(start)-1))
                
                prev_tid = tid
                prev_end = end
                
for i in sorted(intron_set):
    assert i[1] < i[2]
    print(i[0], i[1], i[2], sep='\t')

