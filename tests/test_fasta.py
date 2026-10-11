# -*- coding: utf-8 -*-
"""Tests for the FASTA helpers and external FASTA handling."""
import io
import logging

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


def test_url_escaped_id_is_indexed_under_both_ids_and_written_once():
    fasta_dict, count = fasta_file_to_dict(io.StringIO('>chr%201 escaped\nACGT\n>chr2\nGGGG\n'))
    assert count == 2
    assert list(fasta_dict) == ['chr%201', 'chr 1', 'chr2']
    assert fasta_dict['chr 1'] is fasta_dict['chr%201']
    assert fasta_dict['chr%201']['unescaped_id'] == 'chr 1'
    # ids without escapes are not given a redundant unescaped_id
    assert 'unescaped_id' not in fasta_dict['chr2']
    out = io.StringIO()
    fasta_dict_to_file(fasta_dict, out)
    assert out.getvalue() == '>chr%201 escaped\nACGT\n>chr2\nGGGG\n'


def test_blank_lines_are_ignored():
    fasta_dict, count = fasta_file_to_dict(io.StringIO('\n>seq1\nACGT\n\nTTTT\n\n>seq2\nGG\n\n'))
    assert count == 2
    assert fasta_dict['seq1']['seq'] == 'ACGTTTTT'
    assert fasta_dict['seq2']['seq'] == 'GG'


def test_embedded_fasta_with_blank_line_after_directive(parse_text):
    gff = parse_text('##gff-version 3\nctg1\t.\tgene\t1\t4\t.\t+\t.\tID=g1\n##FASTA\n\n>ctg1\nACGT\n')
    assert gff.fasta_embedded['ctg1']['seq'] == 'ACGT'


def test_duplicate_ids_in_a_stream_are_logged(caplog):
    with caplog.at_level(logging.WARNING, logger='gff3.gff3'):
        fasta_dict, count = fasta_file_to_dict(io.StringIO('>seq1\nAAAA\n>seq1\nCCCC\n'))
    assert count == 2
    assert fasta_dict['seq1']['seq'] == 'CCCC'
    assert '<stream> : Line 4 : Duplicate id [seq1]' in caplog.text


def test_round_trip_embedded_fasta_with_escaped_seqid(parse_text, tmp_path, quiet_logger):
    gff = parse_text('##gff-version 3\nctg%201\t.\tgene\t1\t4\t.\t+\t.\tID=g1\n##FASTA\n>ctg%201\nACGT\n')
    out_path = tmp_path / 'out.gff3'
    gff.write(str(out_path))
    text = out_path.read_text(encoding='utf-8')
    assert text.count('>ctg%201') == 1
    assert text.endswith('##FASTA\n>ctg%201\nACGT\n')
