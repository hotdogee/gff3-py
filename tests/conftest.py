# -*- coding: utf-8 -*-
import logging
from pathlib import Path

import pytest

from gff3 import Gff3

EXAMPLES_DIR = Path(__file__).resolve().parent.parent / 'examples'

# 60 bp reference sequence for ctg1 with a run of Ns at positions 41-44
CTG1_SEQ = 'ATGGCGTACGTTAGCATCGATCGTAGCTAGCTAGGATCCANNNNACGTAGCTAGCATGCA'
assert len(CTG1_SEQ) == 60

# A small, valid GFF3 file: one gene with a two-exon mRNA (with correct CDS phases),
# a multi-parent exon, a multi-line CDS, a non-reserved attribute, and an embedded FASTA
SAMPLE_GFF3 = """\
##gff-version 3
##sequence-region ctg1 1 60
ctg1\texample\tgene\t1\t60\t.\t+\t.\tID=gene1;Name=EDEN;note_lowercase=free text
ctg1\texample\tmRNA\t1\t60\t.\t+\t.\tID=mRNA1;Parent=gene1;Name=EDEN.1
ctg1\texample\tmRNA\t1\t60\t.\t+\t.\tID=mRNA2;Parent=gene1;Name=EDEN.2
ctg1\texample\texon\t1\t20\t.\t+\t.\tID=exon1;Parent=mRNA1,mRNA2
ctg1\texample\texon\t31\t60\t.\t+\t.\tID=exon2;Parent=mRNA1
ctg1\texample\tCDS\t5\t20\t.\t+\t0\tID=cds1;Parent=mRNA1
ctg1\texample\tCDS\t31\t39\t.\t+\t2\tID=cds1;Parent=mRNA1
ctg1\texample\tmatch\t2\t30\t42.5\t-\t.\tID=match1
###
##FASTA
>ctg1 test contig
{0}
{1}
""".format(CTG1_SEQ[:30], CTG1_SEQ[30:])


@pytest.fixture
def quiet_logger():
    """A logger that swallows the validation messages Gff3 logs while parsing."""
    log = logging.getLogger('gff3.tests.quiet')
    log.propagate = False
    if not log.handlers:
        log.addHandler(logging.NullHandler())
    return log


@pytest.fixture
def sample_path(tmp_path):
    path = tmp_path / 'sample.gff3'
    path.write_text(SAMPLE_GFF3, encoding='utf-8', newline='\n')
    return path


@pytest.fixture
def sample(sample_path, quiet_logger):
    gff = Gff3(logger=quiet_logger)
    gff.parse(str(sample_path))
    return gff


@pytest.fixture
def parse_text(tmp_path, quiet_logger):
    """Parse GFF3 text from a file and return the Gff3 object."""
    def _parse_text(text, name='input.gff3'):
        path = tmp_path / name
        path.write_text(text, encoding='utf-8', newline='\n')
        gff = Gff3(logger=quiet_logger)
        gff.parse(str(path))
        return gff
    return _parse_text
