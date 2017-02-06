.. _inputs-outputs:

================
Inputs / Outputs
================


craw_coverage
=============

Inputs
------


bam file
^^^^^^^^

*craw_coverage* need a file of alignment reads called bam file.
a bam file is a short DNA sequence read alignments in the Binary Alignment/Map format (.bam).
*craw_coverage* needs also the corresponding index file (bai). The index file must be located beside the bam file
with the same name instead to have the *.bam* extension it end by *.bai* extension.
If you have not  the index file you have to create it.

To index a bam file you need samtools. The command line is ::

    samtools index file.bam

For more explanation see http://www.htslib.org/doc/ .


annotation file
^^^^^^^^^^^^^^^

The annotation file is a `tsv` file. It's mean that it is a text file with value separated by tabulation (not spaces).
The first line of the file must be the name of the columns
the other lines the values. Each line represent a row. ::

    name    gene    chromosome      strand  Position
    YEL072W RMD6    chrV    +       14415
    YEL071W DLD3    chrV    +       17845
    YEL070W DSF1    chrV    +       21097
    YEL066W HPA3    chrV    +       27206
    YEL065W SIT1    chrV    +       29543
    YEL062W NPR2    chrV    +       36254
    YEL058W PCM1    chrV    +       44925
    YEL056W HAT2    chrV    +       48373


All lines starting with '#' character will be ignored. ::

    # This is the annotation file for Wild type
    # bla bla ...
    name    gene    chromosome      strand  Position
    YEL072W RMD6    chrV    +       14415
    YEL071W DLD3    chrV    +       17845
    YEL070W DSF1    chrV    +       21097
    YEL066W HPA3    chrV    +       27206
    YEL065W SIT1    chrV    +       29543
    YEL062W NPR2    chrV    +       36254
    YEL058W PCM1    chrV    +       44925
    YEL056W HAT2    chrV    +       48373


mandatory columns
"""""""""""""""""

columns with fixed name
'''''''''''''''''''''''


columns with variable name
''''''''''''''''''''''''''

::

    name    gene    type    chromosome      strand  annotation_start        annotation_end  has_transcript  transcription_end       transcription_start
    YEL072W RMD6    gene    chrV    1       13720   14415   1       14745   13569
    YEL071W DLD3    gene    chrV    1       16355   17845   1       17881   16177
    YEL070W DSF1    gene    chrV    1       19589   21097   1       21197   19539
    YEL066W HPA3    gene    chrV    1       26721   27206   1       27625   26137
    YEL065W SIT1    gene    chrV    1       27657   29543   1       29601   27625
    YEL062W NPR2    gene    chrV    1       34407   36254   1       36401   34321
    YEL058W PCM1    gene    chrV    1       43252   44925   1       44993   43217
    YEL056W HAT2    gene    chrV    1       47168   48373   1       48457   47105
    YEL052W AFG1    gene    chrV    1       56571   58100   1       58105   56537


All other columns are not necessary but will report as is in coverage file.


.. _cov_out:

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

see `cov_out`_

Outputs
-------

The default output of *craw_htmp* (if --out is omitted) is grapical window on the screen.
The figure display on the screen can be saved using the window menu.
It is also possible to generate directly a image file in various format by specifying the --out option.
The output format will be deduced form the filename extension provide to --out option. ::

  --out foo.jpeg  for jpeg image or --out foo.png  for png image

The supported format vary in function of the matloblib backend used (see :ref:`matplotlibrc`).
