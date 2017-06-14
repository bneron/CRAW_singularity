.. _release_notes:


=============
Release notes
=============

Release 0.9.0
=============

BUG FIX
-------

* **craw_coverage** When feature is located on reverse strand the coverage must be reverse to be read in 5'->3'

 (`issue 27 <https://gitlab.pasteur.fr/bneron/craw/issues/27>`_)


New Feature
-----------

* **craw_coverage** add support of 2 wig files one for the forward strand the other for the reverse.

  (`issue 26 <https://gitlab.pasteur.fr/bneron/craw/issues/26>`_)

* **craw_htmp** allow to draw a vertical line in a given color at a given position.

  see \-\-mark option (`issue 25 <https://gitlab.pasteur.fr/bneron/craw/issues/25>`_)


Release 0.8.0
=============

This release implements the support of wig files.

* Support of wig file as input data.
  The user can not only provide bam file but also wig file with --wig option.
* The annotation file support all columns file, not only tabular file (tsv).
  The tabular format is still the default, but by specifying the columns separator
  with the --sep option the annotation file can be in any format
  for csv file use --sep ','
* the control of verbosity has been extend with -v to increase it and -q to decrease it.
  These options are cumulative so it possible to use -vv or -qq
* improve traceability by displaying craw version and all python dependencies version
  pysam, pandas, numpy, matplotlib, ...