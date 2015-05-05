import sys
# try to import from project first
from os.path import dirname
sys.path.insert(1, dirname(dirname(__file__)))
from gff3 import Gff3

gff = Gff3('annotations.gff')
type_map = {'exon': 'pseudogenic_exon', 'transcript': 'pseudogenic_transcript'}
pseudogenes = [line for line in gff.lines if line['type'] == 'pseudogene']
for pseudogene in pseudogenes:
    # convert types
    for line in gff.descendants(pseudogene):
        if line['type'] in type_map:
            line['type'] = type_map[line['type']]
    # find overlapping gene
    overlapping_genes = [line for line in gff.lines if line['type'] == 'gene' and gff.overlap(line, pseudogene)]
    if overlapping_genes:
        # move pseudogene children to overlapping gene
        gff.adopt(pseudogene, overlapping_genes[0])
        # remove pseudogene
        gff.remove(pseudogene)
gff.write('annotations_fixed.gff')