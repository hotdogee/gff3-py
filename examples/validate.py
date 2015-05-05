import sys
# try to import from project first
from os.path import dirname
sys.path.insert(1, dirname(dirname(__file__)))
from gff3 import Gff3

# initialize a Gff3 object
gff = Gff3()
# parse GFF3 file and do syntax checking, this populates gff.lines and gff.features
# if an embedded ##FASTA directive is found, parse the sequences into gff.fasta_embedded
gff.parse('annotations.gff')
# parse the external FASTA file into gff.fasta_external
#gff.parse_fasta_external('annotations.fa')
# Check seqid, bounds and the number of Ns in each feature using one or more reference sources
gff.check_reference(allowed_num_of_n=0, feature_types=['CDS'])
# Checks whether child features are within the coordinate boundaries of parent features
gff.check_parent_boundary()
# Calculates the correct phase and checks if it matches the given phase for CDS features
gff.check_phase()