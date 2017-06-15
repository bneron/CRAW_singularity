.. Counter RNAseq Window documentation master file, created by
   sphinx-quickstart on Tue Sep 27 13:03:09 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Counter RNAseq Window's documentation!
=================================================

User Guide
==========

.. toctree::
   :maxdepth: 2

   overview
   installation
   quickstart
   inputs-outputs

Developer Guide
===============

Overview
--------

Scripts are located in *bin* directory, and use some modules located in *craw* directory.

* *craw_coverage* use module :mod:`craw.annotation` to handle annotation file and module :mod:`craw.coverage`
  to compute coverage this module rely on pysam.
* *craw_htmp* read coverage file generate by *craw_coverage* and produce graphical representation of data.
  This script use functions in module :mod:`craw.heatmap` in the form of heatmap.
  The module :mod:`craw.heatmap` have some capabilities to sort, crop, normalize data before represent them.
  this module rely on numpy, pandas to manipulate data (:mod:`craw.heatmap.sort`, :mod:`craw.heatmap.crop_matrix`,
  :mod:`craw.heatmap.lin_norm`, ...) and matplotlib and/or pillow to generate images (:mod:`craw.heatmap.draw_heatmap` ,
  :mod:`craw.heatmap.draw_raw_image`)

reference API
-------------

.. toctree::
   :maxdepth: 5

   craw
   annotation
   wig
   coverage
   heatmap
   argparse_util

.. only:: html

   unit tests coverage
   -------------------

   http://bneron.pages.pasteur.fr/craw/htmlcov


Release notes
-------------

.. toctree::
   :maxdepth: 2

   release_notes


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

