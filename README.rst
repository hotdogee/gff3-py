===============================
gff3-py
===============================

.. image:: https://badge.fury.io/py/gff3.png
   :target: http://badge.fury.io/py/gff3

.. image:: https://travis-ci.org/hotdogee/gff3-py.png?branch=master
   :target: https://travis-ci.org/hotdogee/gff3-py

.. image:: https://pypip.in/d/gff3/badge.png
   :target: https://pypi.python.org/pypi/gff3


Manipulate genomic features and validate the syntax and reference sequence of your ```GFF3`_`` files.

* Free software: BSD license
* Documentation: https://gff3-py.readthedocs.org.

Features
--------

* Parses a ```GFF3```_ file into a structure composed of simple python ``dict`` and ``list``.
* Validates the ```GFF3 <http://www.sequenceontology.org/gff3.shtml>`_`` syntax on parse, and saves the error messages in the parsed structure.
* Uses the python ```logging``<http://www.sequenceontology.org/gff3.shtml>`_ library to log error messages with support for custom loggers.
* Parses embeded or external ``FASTA`` sequences to check bounds and number of ``N``s.
* Check and correct the phase for ``CDS`` features.
* Tree traversal methods ``ancestors`` and ``descendants`` returns a simple ``list`` in Breadth-first search order.
* Transfer children and parents using the ``adopt`` and ``adopted`` methods.
* Test for overlapping features using the ``overlap`` method.
* Remove a feature and its associated features using the ``remove`` method.
* Write the modified structure to a GFF3 file using the ``write`` mthod.


.. _GFF3: http://www.sequenceontology.org/gff3.shtml