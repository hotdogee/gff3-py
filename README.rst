===============================
gff3-py
===============================

.. image:: https://badge.fury.io/py/gff3.png
   :target: http://badge.fury.io/py/gff3

.. image:: https://travis-ci.org/hotdogee/gff3-py.png?branch=master
   :target: https://travis-ci.org/hotdogee/gff3-py

.. image:: https://pypip.in/d/gff3/badge.png
   :target: https://pypi.python.org/pypi/gff3


Manipulate genomic features and validate the syntax and reference sequence of your |GFF3|_ files.

* Free software: BSD license
* Documentation: https://gff3-py.readthedocs.org.

Features
--------

* Parses a |GFF3|_ file into a structure composed of simple python |dict| and |list|.
* Validates the |GFF3|_ syntax on parse, and saves the error messages in the parsed structure.
* Uses the python |logging| library to log error messages with support for custom loggers.
* Parses embeded or external |FASTA| sequences to check bounds and number of ``N``s.
* Check and correct the phase for ``CDS`` features.
* Tree traversal methods ``ancestors`` and ``descendants`` returns a simple ``list`` in Breadth-first search order.
* Transfer children and parents using the ``adopt`` and ``adopted`` methods.
* Test for overlapping features using the ``overlap`` method.
* Remove a feature and its associated features using the ``remove`` method.
* Write the modified structure to a GFF3 file using the ``write`` mthod.

.. |GFF3| replace:: ``GFF3``
.. |dict| replace:: ``dict``
.. |list| replace:: ``list``
.. |logging| replace:: ``logging``
.. |FASTA| replace:: ``FASTA``

.. _GFF3: http://www.sequenceontology.org/gff3.shtml
.. _dict: https://docs.python.org/2/tutorial/datastructures.html#dictionaries
.. _list: https://docs.python.org/2/tutorial/datastructures.html#more-on-lists
.. _logging: https://docs.python.org/2/library/logging.html
.. _FASTA: http://en.wikipedia.org/wiki/FASTA_format