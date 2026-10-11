.. :changelog:

History
-------

1.1.0 (unreleased)
---------------------

* Requires Python 3.10 or newer. Python 2 and Python 3.3 to 3.9 are no longer supported; gff3 1.0.1 remains available for them.
* Fixed ``Gff3.write()`` raising ``TypeError: a bytes-like object is required, not 'str'`` (#9).
* Fixed ``Gff3.parse_fasta_external()`` and ``fasta_file_to_dict()`` failing on file paths, and ``fasta_dict_to_file()`` not closing the file it opened.
* GFF3 and FASTA files given as paths are read and written as UTF-8 on every platform, bytes that are not valid UTF-8 are kept unchanged, written files always use ``\n`` line endings, and paths may be ``pathlib.Path`` objects. Open file objects are used as they are, so open them with ``encoding='utf-8', errors='surrogateescape'`` (and ``newline='\n'`` for writing) to get the same behaviour.
* ``Gff3.write()`` no longer drops children listed before their parent or the ``##sequence-region`` directives of files without FASTA, writes the ``Target`` attribute correctly, writes ``.`` instead of an empty column 9 for features without attributes, writes each sequence of the embedded FASTA once, and copies lines it could not parse instead of raising.
* Forward ``Parent`` references are linked to their parents, and ``Parent`` IDs that are never defined are reported as errors.
* Fixed crashes in ``parse()`` on a duplicate ID after a directive or a line without an ID, in ``check_phase()`` on a ``CDS`` phase of ``.`` (#6) or inconsistent ``CDS`` strands, in ``check_reference()`` on a second feature with a seqid missing from the reference, and in the validation checks on features whose start or end is not an integer.
* Blank lines are kept in ``Gff3.lines`` with the ``line_type`` ``'blank'``, so ``line_index`` and the line numbers in error messages match the file after a blank line. A blank first line is now reported as ``"##gff-version" missing from the first line``.
* ``Gff3.remove(line_data, root_type=...)`` removes the whole tree of the ancestor whose type is ``root_type``, as documented. Before, ``root_type`` never matched and only ``line_data`` and its descendants were removed.
* Fixed parsing of blank lines in FASTA sequences.
* ``Gff3.sequence()`` returns the sequence of a feature, or of its children of ``child_type``, instead of ``None``.
* The example scripts run on Python 3.
* Packaging moved to ``pyproject.toml``. Tests run with pytest on GitHub Actions, and releases are published to PyPI with trusted publishing.
* The package metadata now declares the ISC license, which is the license in ``LICENSE``. Earlier releases were labelled BSD by mistake.

1.0.1 (2021-11-12)
---------------------

* new source distribution uploaded

1.0.0 (2018-12-01)
---------------------

* Fix Python3 issues
* Added sequence functions: complement(seq) and translate(seq)
* Added fasta write function: fasta_dict_to_file(fasta_dict, fasta_file, line_char_limit=None)
* Added Gff method to return the sequence of line_data: sequence(self, line_data, child_type=None, reference=None)
* Gff.write no longer prints redundent '###' when the whole gene is marked as removed


0.3.0 (2015-03-10)
---------------------

* Fixed phase checking.

0.2.0 (2015-01-28)
---------------------

* Supports python 2.6, 2.7, 3.3, 3.4, pypy.
* Don't report empty attributes as errors.
* Improved documentation.

0.1.0 (2014-12-11)
---------------------

* First release on PyPI.
