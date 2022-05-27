import pysam
import sys
from collections import defaultdict
import re

def get_secondary_alignments(aln):
    mappings = defaultdict(list)
    if aln.has_tag('SA'):
        for mapping in aln.get_tag('SA').split(';'):
            if mapping:
                chrom, start, strand, cigar, mapq, nm = mapping.split(',')
                mappings[chrom].append([int(start), strand, cigar, int(mapq), int(nm)])

    return mappings

def find_splits(aln):
    ends = []
    if aln.cigartuples[-1][0] >= 4:
        ends.append('end')
    if aln.cigartuples[0][0] >= 4:
        ends.append('start')

    return ends

def find_splits_old(aln, min_split_size=100, closeness_in_size=200):
    """ aln must have been checked if it has SA """
    '''
    clipped = None
    partner_start = None
    clipped_size1 = None
    clipped_size2 = None
    '''
    clipped_start_regex = re.compile('^(\d+)[S|H]')
    clipped_end_regex = re.compile('(\d+)[S|H]$')

    strand = '-' if aln.is_reverse else '+'
    splits = []
    sas = get_secondary_alignments(aln)
    for chrom in sorted(sas.keys()):
        if chrom == aln.reference_name:
            continue
        
        if aln.cigartuples[-1][0] >= 4 and aln.cigartuples[-1][1] >= min_split_size:
            clipped_size1 = aln.cigartuples[-1][1]
            for sa in sas[chrom]:
                if sa[1] != strand:
                    continue
                    
                match_start = clipped_start_regex.search(sa[2])
                if match_start:
                    clipped_size2 = int(match_start.group(1))

                    if abs(clipped_size1 - clipped_size2) <= closeness_in_size:
                        clipped = 'end'
                        partner_start = sa[0]
                        print('gg', aln.query_name, strand, clipped, partner_start, aln.cigarstring, sa)
                        splits.append((clipped, partner_start))

        if aln.cigartuples[0][0] >= 4 and aln.cigartuples[0][1] >= min_split_size:
            clipped_size1 = aln.cigartuples[0][1]
            for sa in sas[chrom]:
                if sa[1] != strand:
                    continue

                match_end = clipped_end_regex.search(sa[2])
                if match_end:
                    clipped_size2 = int(match_end.group(1))

                    if abs(clipped_size1 - clipped_size2) <= closeness_in_size:
                        clipped = 'start'
                        partner_start = sa[0]
                        print('gg', aln.query_name, strand, clipped, partner_start, aln.cigarstring, sa)
                        splits.append((clipped, partner_start))

    return splits

def get_reads(bam, ref, start, end, strand):
    reads = {}
    for aln in bam.fetch(ref, start, end):
        #if aln.is_proper_pair:
        #    continue
        if (strand == '+' and not aln.is_reverse) or (strand == '-' and aln.is_reverse):
            mate = 1 if aln.is_read1 else 2
            reads[aln.query_name] = mate

    return reads

c2g_file = sys.argv[1]
c2g = pysam.AlignmentFile(c2g_file)
r2c_file = sys.argv[2]
r2c = pysam.AlignmentFile(r2c_file)

splits = {}
for aln in c2g.fetch(until_eof=True):
    if aln.is_unmapped:
        continue
    if not aln.query_name in splits:
        splits[aln.query_name] = defaultdict(list)
    ends = find_splits(aln)
    for end in ends:
        splits[aln.query_name][end].append(aln)

w = 200
check_reads = []
for ctg in sorted(splits.keys()):
    if 'start' in splits[ctg] and 'end' in splits[ctg]:
        for aln1 in splits[ctg]['end']:
            for aln2 in splits[ctg]['start']:
                ref1_rlen = c2g.get_reference_length(aln1.reference_name)
                ref2_rlen = c2g.get_reference_length(aln2.reference_name)
                if aln1 == aln2 or aln1.is_reverse != aln2.is_reverse or aln1.reference_name == aln2.reference_name:
                    continue
                print('ww', ctg, aln1.reference_name, aln1.reference_start, aln1.reference_end, ref1_rlen, aln1.is_reverse, aln1.cigarstring, aln2.reference_name, aln2.reference_start, aln2.reference_end, ref2_rlen, aln2.is_reverse, aln2.cigarstring)

                check_reads.append((ctg,\
                                    (aln1.reference_name, max(0, aln1.reference_end - w), aln1.reference_end),\
                                    (aln2.reference_name, aln2.reference_start, min(ref2_rlen, aln2.reference_start + w)),\
                                    (aln1.cigarstring, aln2.cigarstring)))

for ctg, span1, span2, cigars in check_reads:
    reads1 = get_reads(r2c, span1[0], span1[1], span1[2], '+')
    reads2 = get_reads(r2c, span2[0], span2[1], span2[2], '-')

    paired = []
    for read in reads1.keys():
        if read in reads2 and reads1[read] != reads2[read]:
            paired.append(read)
            #print('rr', ctg, read, reads1[read], reads2[read])
    
    print('qq', ctg, '{}:{}-{}'.format(span1[0], span1[1], span1[2]), cigars[0], '{}:{}-{}'.format(span2[0], span2[1], span2[2]), cigars[1], len(reads1), len(reads2), len(paired))

