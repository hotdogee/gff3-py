# -*- coding: utf-8 -*-
"""Tests for the FASTA helpers and external FASTA handling."""
import io

from gff3 import Gff3
from gff3.gff3 import fasta_dict_to_file, fasta_file_to_dict

FASTA = '>seq1 first sequence\nACGTacgt\nNNNN\n>seq2\nTTTTGGGG\n'


def test_fasta_file_to_dict_from_path(tmp_path):
    path = tmp_path / 'seqs.fa'
    path.write_text(FASTA, encoding='utf-8')
    fasta_dict, count = fasta_file_to_dict(str(path))
    assert count == 2
    assert list(fasta_dict) == ['seq1', 'seq2']
    assert fasta_dict['seq1']['header'] == '>seq1 first sequence'
    assert fasta_dict['seq1']['id'] == 'seq1'
    # sequences are joined and upper-cased
    assert fasta_dict['seq1']['seq'] == 'ACGTACGTNNNN'
    assert fasta_dict['seq2']['seq'] == 'TTTTGGGG'


def test_fasta_file_to_dict_from_file_object_and_pathlib(tmp_path):
    path = tmp_path / 'seqs.fa'
    path.write_text(FASTA, encoding='utf-8')
    from_path, _ = fasta_file_to_dict(path)
    from_stream, _ = fasta_file_to_dict(io.StringIO(FASTA))
    assert from_path == from_stream


def test_fasta_dict_to_file_round_trip(tmp_path):
    fasta_dict, _ = fasta_file_to_dict(io.StringIO(FASTA))
    out_path = tmp_path / 'out.fa'
    fasta_dict_to_file(fasta_dict, str(out_path))
    assert out_path.read_bytes() == b'>seq1 first sequence\nACGTACGTNNNN\n>seq2\nTTTTGGGG\n'
    reread, count = fasta_file_to_dict(str(out_path))
    assert count == 2
    assert reread == fasta_dict


def test_fasta_dict_to_file_line_char_limit():
    fasta_dict, _ = fasta_file_to_dict(io.StringIO(FASTA))
    out = io.StringIO()
    fasta_dict_to_file(fasta_dict, out, line_char_limit=5)
    assert out.getvalue() == '>seq1 first sequence\nACGTA\nCGTNN\nNN\n>seq2\nTTTTG\nGGG\n'


def test_parse_fasta_external_from_path(tmp_path, quiet_logger):
    path = tmp_path / 'ref.fa'
    path.write_text(FASTA, encoding='utf-8')
    gff = Gff3(logger=quiet_logger)
    gff.parse_fasta_external(str(path))
    assert gff.fasta_external['seq2']['seq'] == 'TTTTGGGG'


def test_external_fasta_is_written_as_embedded_fasta(tmp_path, quiet_logger):
    fasta_path = tmp_path / 'ref.fa'
    fasta_path.write_text(FASTA, encoding='utf-8')
    gff_path = tmp_path / 'in.gff3'
    gff_path.write_text('##gff-version 3\nseq2\t.\tgene\t1\t8\t.\t+\t.\tID=g1\n', encoding='utf-8')
    gff = Gff3(str(gff_path), fasta_external=str(fasta_path), logger=quiet_logger)
    out_path = tmp_path / 'out.gff3'
    gff.write(str(out_path))
    reparsed = Gff3(str(out_path), logger=quiet_logger)
    assert reparsed.fasta_embedded['seq2']['seq'] == 'TTTTGGGG'
    assert reparsed.fasta_embedded['seq1']['seq'] == 'ACGTACGTNNNN'
    # the ##sequence-region directive is generated from the FASTA sequence length
    assert reparsed.lines[1]['directive'] == '##sequence-region'
    assert (reparsed.lines[1]['seqid'], reparsed.lines[1]['start'], reparsed.lines[1]['end']) == ('seq2', 1, 8)
