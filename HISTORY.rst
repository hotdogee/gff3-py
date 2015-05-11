.. :changelog:

History
-------

0.4.0 (2015-05-05)
---------------------

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
