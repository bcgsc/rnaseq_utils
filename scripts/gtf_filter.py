import argparse

# Written by Ka Ming Nip @kmnip

def get_tid_from_attribute_col(col, fix_id):
    for info in col.split(';'):
        key, val = info.strip().split()
        if key == 'transcript_id':
            tid = val.strip('"')
            if fix_id:
                tid = tid.split('.')[0]
            return tid
    return None
    
parser = argparse.ArgumentParser(description='Filter a GTF file based on a list of transcript IDs')
parser.add_argument('gtf', help='path of input GTF file')
parser.add_argument('tids', help='path of transcript ID list file')
parser.add_argument('--fix', action='store_true', help='fix transcript IDs by removing `.` and trailing characters')
args = parser.parse_args()

fix_id = args.fix
tids_set = set()
with open(args.tids) as fh:
    if fix_id:
        for line in fh:
            tids_set.add(line.strip().split('.')[0])
    else:
        for line in fh:
            tids_set.add(line.strip())

with open(args.gtf) as fh:
    features = ['exon', 'transcript', 'start_codon', 'stop_codon', 'CDS', 'UTR']

    for line in fh:
        line = line.strip()
        if len(line.strip()) > 0 and line[0] != '#':
            cols = line.split('\t')
            if cols[2] in features :
                tid = get_tid_from_attribute_col(cols[8], fix_id)
                if tid and tid in tids_set:
                    print(line)
            else:
                print(line)
        else:
            print(line)

