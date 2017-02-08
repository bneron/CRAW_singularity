Counter RnAseq Window
=====================


Installation
============

requirements
------------

- python > 3
- pysam >= 0.9.1.4
- pandas >= 0.17.1
- numpy >= 1.11.2
- matplotlib >= 1.5.3

from package
------------

using pip

`pip install craw-x.x.x.tar.gz`

if you use virtualenv do not forget to configure the matplotlib backend


from repository
---------------

clone the project and install with the setup.py

`git clone https://gitlab.pasteur.fr/bneron/craw.git`

`cd craw`

`python3 setup.py install`

You can also use the package without install it.
You have to export the **CRAW_HOME** environment variable.
Then you can use it directly

Testing my installation
-----------------------

The package come from with some functional tests.
to test if everything is correctly installed.

`python3 tests/run_tests.py -vvv`


Quickstart
==========

.. _inputs-outputs:

================
Inputs / Outputs
================


craw_coverage
=============

Inputs
------

bam file

*craw_coverage* need a file of alignment reads called bam file.
a bam file is a short DNA sequence read alignments in the Binary Alignment/Map format (.bam).
*craw_coverage* needs also the corresponding index file (bai). The index file must be located beside the bam file
with the same name instead to have the *.bam* extension it end by *.bai* extension.
If you have not  the index file you have to create it.

To index a bam file you need samtools. The command line is ::

    samtools index file.bam

For more explanation see http://www.htslib.org/doc/ .


annotation file

The annotation file is a `tsv` file. It's mean that it is a text file with value separated by tabulation (not spaces).
The first line of the file must be the name of the columns
the other lines the values. Each line represent a row. ::

    name    gene    chromosome      strand  Position
    YEL072W RMD6    chrV    +       14415
    YEL071W DLD3    chrV    +       17845
    YEL070W DSF1    chrV    +       21097


All lines starting with '#' character will be ignored. ::

    # This is the annotation file for Wild type
    # bla bla ...
    name    gene    chromosome      strand  Position
    YEL072W RMD6    chrV    +       14415
    YEL071W DLD3    chrV    +       17845
    YEL070W DSF1    chrV    +       21097


mandatory columns


There is 3 mandatory columns in the annotation file.

columns with fixed name

two with a fixed name:

* **strand** indicate on which strand is located the region of interest. The authorized values for this columns are +/- , 1/-1 or for/rev.
* **chromosome** the chromosome name where is located the region of interest.

columns with variable name


In addition of these two columns the column to define the position of reference is mandatory too, but the name of this
column can be specified by the user. If it's not craw_coverage will use a column name 'position'.

If we want to compute coverage on variable window size, 2 extra columns whose name must be specified by the user by the following option:

* \-\-start-col to define the beginning of the window (this position is included in the window)
* \-\-stop-col to define the end of the window (this position is included in the window)

::

    name    gene    type    chromosome      strand  annotation_start        annotation_end  has_transcript  transcription_end       transcription_start
    YEL072W RMD6    gene    chrV    1       13720   14415   1       14745   13569
    YEL071W DLD3    gene    chrV    1       16355   17845   1       17881   16177
    YEL070W DSF1    gene    chrV    1       19589   21097   1       21197   19539



::

    craw_coverage --bam file.bam --annot annot.txt --ref-col annotation_start --start-col annotation_start --stop-col annotation_end


The position of reference must be between start and end.
The authorized values are positive integers.

.. note::
    the position of reference can be used to define the reference and the start ot the end of the window. ::

        craw_coverage --bam file.bam --annot annot.txt --ref-col annotation_start --start-col annotation_start --stop-col annotation_end

All other columns are not necessary but will be reported as is in the coverage file.



Outputs
-------

coverage_file
^^^^^^^^^^^^^

It's a `tsv` file with all columns found in annotation file plus the result of coverage position by position centered
on the reference position define for each line. for instance ::

    craw_coverage -bam=../data/craw_data_test/WTE1.bam --annot=../data/craw_data_test/annotations.txt
    --ref-col=annotation_start --before=0  --after=2000

In the command line above, the column '0' correspond to the annotation_start position the column '1' to annotation_start + 1
on so on until '2000' (here we display only the first 3 columns of the coverage). ::

    # Running Counter RnAseq Window
    # Version: craw NOT packaged, it should be a development version | Python 3.4
    # With the following arguments:
    # --after=2000
    # --annot=../data/craw_data_test/annotations.txt
    # --bam=../data/craw_data_test/WTE1.bam
    # --before=0
    # --output=WTE1_0+2000.new.cov
    # --qual-thr=15
    # --ref-col=annotation_start
    # --suffix=cov
    sense   name    gene    type    chromosome      strand  annotation_start        annotation_end  has_transcript  transcription_end       transcription_start     0       1       2
    S       YEL072W RMD6    gene    chrV    +       13720   14415   1       14745   13569   7       7       7
    AS      YEL072W RMD6    gene    chrV    +       13720   14415   1       14745   13569   0       0       0
    S       YEL071W DLD3    gene    chrV    +       16355   17845   1       17881   16177   31      33      33


The line starting with '#' are comments and will be ignored for further processing.
But in traceability/reproducibility concern, in the comments `craw_coverage` indicate
the version of the program and the arguments used for this experiment.



craw_htmp
=========

Inputs
------

see craw_covergae output

Outputs
-------

The default output of *craw_htmp* (if --out is omitted) is grapical window on the screen.
The figure display on the screen can be saved using the window menu.
It is also possible to generate directly a image file in various format by specifying the --out option.
The output format will be deduced form the filename extension provide to --out option. ::

  --out foo.jpeg  for jpeg image or --out foo.png  for png image

The supported format vary in function of the matloblib backend used (see ).


command line options
--------------------

TODO