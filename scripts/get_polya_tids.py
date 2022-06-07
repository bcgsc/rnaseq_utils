from pybedtools import BedTool
import argparse

parser = argparse.ArgumentParser(description='Extract polyA transcript IDs.')
parser.add_argument('gtf', metavar='GTF',
                    help='annotation GTF file')
parser.add_argument('atlas', metavar='BED',
                    help='polyA site atlas BED file')
args = parser.parse_args()

gtf = args.gtf # e.g. `Mus_musculus.GRCm39.105.gtf.gz` (from Ensembl)
atlas = args.atlas # e.g. `atlas.clusters.2.0.GRCm38.96.Mm39_liftover.bed` (liftover from PolyASite)

# extract number of exons for each transcript
num_exons_dict = dict()
for t in BedTool(gtf).filter(lambda x : x[2] == 'exon'):
    tid = t.attrs['transcript_id']
    if tid in num_exons_dict:
        num_exons_dict[tid] += 1
    else:
        num_exons_dict[tid] = 1

terminus_exons = BedTool(gtf).filter(lambda x : x[2] == 'exon' and int(x.attrs['exon_number']) == num_exons_dict[x.attrs['transcript_id']])

polya_sites = BedTool(atlas)

polya_exons = terminus_exons.intersect(polya_sites, s=True)

polya_transcript_ids = set()

for r in polya_exons:
    polya_transcript_ids.add(r.attrs['transcript_id'])

for t in sorted(polya_transcript_ids):
    print(t)
